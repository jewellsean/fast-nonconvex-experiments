# optimize_parameters
# Optimize over l0 and l1 solution paths

optimize_parameters <- function(train_dat, data_params, optim_params, cost_params) {
    l0_optimal <- optimize_l0(train_dat, data_params, optim_params, cost_params)
    l1_optimal <- optimize_l1(train_dat, data_params, optim_params, cost_params)
    return(list(l0 = l0_optimal, l1 = l1_optimal))
}

optimize_l0 <- function(train_dat, data_params, optim_params, cost_params) {
    out <- grid_min(min_spike_l0, train_dat = train_dat, data_params = data_params, optim_params = optim_params, cost_params = cost_params)
    return(out)
}

optimize_l1 <- function(train_dat, data_params, optim_params, cost_params) {
    out <- l1_paths(train_dat, data_params, optim_params, cost_params)
    return(out)
}

l1_paths <- function(train_dat, data_params, optim_params, cost_params) {
    df <- NULL
    lams <- optim_params$l1_lams
    nLams <- length(lams)

    thresholds <- optim_params$l1_thresholds
    nThresholds <- length(thresholds)

    prefactor <- 1
    if (data_params$measure == "Cor") {
        prefactor <- -1
    }

    for (i in 1:nLams) {
        fit <- oasis(train_dat$trace, data_params$gamma, lams[i], "penalized", "test")
        n_thresholds <- length(thresholds)
        for (j in 1:nThresholds) {
            out_fit <- l1_make_into_spikes(train_dat, data_params, fit, thresholds[j])
            df <- rbind(df,
            data.frame(dist =
            prefactor * calculate_spike_distance(out_fit$train_spikes,
            out_fit$fit_spikes,
            data_params$measure, cost_params),
            threshold = thresholds[j],
            lam = lams[i]))
        }

    }

    min_val <- min(df$dist)
    min_inds <- (df$dist == min_val)
    min_arr <- df[min_inds, ]
    out <- list(minimum = mean(min_arr$lam), objective = mean(min_arr$dist),
    threshold = mean(min_arr$threshold))

    return(out)
}

grid_min <- function(f, train_dat, data_params, optim_params, cost_params) {
    lams <- optim_params$l0_lams
    nLams <- length(lams)
    eval_f <- numeric(nLams) + Inf

    for (lam_i in 1:nLams) {
        eval_f[lam_i] <- f(lams[lam_i], train_dat = train_dat, data_params = data_params, optim_params = optim_params, cost_params = cost_params)
    }

    min_val <- min(eval_f)
    min_inds <- (eval_f == min_val)
    out <- list(minimum = mean(lams[min_inds]), objective = mean(eval_f[min_inds]))
    return(out)
}


min_spike_l1 <- function(lam, train_dat, params, thresholds) {
    df <- l1_paths(lam, train_dat, params, thresholds)
    return(min(df$dist))
}


## L0 problem
min_spike_l0 <- function(lam, train_dat, data_params, optim_params, cost_params) {
    fit <- estimate_spikes(train_dat$trace, data_params$gamma, lam,
                            constraint = optim_params$l0_constraint,
                            estimate_calcium = TRUE,
                            EPS = optim_params$l0_min_cal_conc)
    out_fit <- l0_make_into_spikes(train_dat, data_params, fit)
    prefactor <- 1
    if (data_params$measure == "Cor") {prefactor <- -1}

    return(prefactor * calculate_spike_distance(out_fit$train_spikes, out_fit$fit_spikes, data_params$measure, cost_params))
}

l0_make_into_spikes <- function(train_dat, data_params, fit) {
    fps <- data_params$fps
    if (data_params$measure == "Cor") {
        train_spikes <- train_dat$spikes
        fit_spikes <- 0 * numeric(length(fit$estimated_calcium))
        st <- get_spike_magnitudes(fit)
        fit_spikes[fit$spikes] <- st
    } else {
        train_spikes <- which(train_dat$spikes > 0) * (1 / fps)
        fit_spikes <- fit$spikes * (1 / fps)
    }
    return(list(train_spikes = train_spikes, fit_spikes = fit_spikes))
}



l1_make_into_spikes <- function(train_dat, params, fit, threshold) {
    spikes <-  which(fit$st > threshold)
    fps <- params$fps

    if (params$measure == "Cor") {
        train_spikes <- train_dat$spikes
        fit_spikes <- fit$st * (fit$st > threshold)
    } else {
        train_spikes <- which(train_dat$spikes > 0) * (1 / fps)
        fit_spikes <- spikes * (1 / fps)
    }

    return(list(train_spikes = train_spikes,
    fit_spikes = fit_spikes))
}



evaluate_performance <- function(test_dat, optimal_params, optimal_params_test, data_params, optim_params, cost_params) {
    ## Calculate fits
    fit_l0 <- estimate_spikes(test_dat$trace, data_params$gamma, optimal_params$l0$minimum,
                                constraint = optim_params$l0_constraint,
                                estimate_calcium = TRUE,
                                EPS = optim_params$l0_min_cal_conc)
    fit_l0_out <- l0_make_into_spikes(test_dat, data_params, fit_l0)
    eval_l0 <- calculate_spike_distance(fit_l0_out$train_spikes, fit_l0_out$fit_spikes, data_params$measure, cost_params)
    l0_summary <- list(fit = fit_l0, eval = eval_l0, max_n_intervals = max(fit_l0$n_intervals))

    fit_l1 <- oasis(test_dat$trace, data_params$gamma, optimal_params$l1$minimum, "penalized", "test")
    fit_l1_out <- l1_make_into_spikes(test_dat, data_params, fit_l1, optimal_params$l1$threshold)
    eval_l1 <- calculate_spike_distance(fit_l1_out$train_spikes, fit_l1_out$fit_spikes, data_params$measure, cost_params)

    ## write out spike times for l1
    fit_l1$spikes <- which(fit_l1_out$st > optimal_params$l1$threshold)
    l1_summary <- list(fit = fit_l1, eval = eval_l1)

    write_results <- data.frame(l0_optim = eval_l0, l0_lam = optimal_params$l0$minimum, l0_max_n_intervals_test = l0_summary$max_n_intervals,
    l0_optim_test = optimal_params_test$l0$objective, l0_lam_test = optimal_params_test$l0$minimum,
    l1_optim = eval_l1, l1_lam = optimal_params$l1$minimum, l1_threshold = optimal_params$l1$threshold,
    l1_optim_test = optimal_params_test$l1$objective, l1_lam_test = optimal_params_test$l1$minimum,
    l1_threshold_test = optimal_params_test$l1$threshold)

    return(list(l0 = l0_summary, l1 = l1_summary, out = write_results))
}
