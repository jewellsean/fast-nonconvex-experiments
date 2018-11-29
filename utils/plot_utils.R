required_packages <- c("latex2exp", "magrittr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(latex2exp)
library(magrittr)

plot_estimates <- function(times, calcium_concentration, spike_times, estimated_spikes, ...) {
    ## some default settings, to add the spike lines nicely
    ylim <- c(-2, 10)
    ys <- -0.2
    lineSize <- 1
    ye <- ys - lineSize
    buff <- 0.1
    incr <- lineSize + buff

    ## base plot of fl
    plot(times, calcium_concentration, pch = 20, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, ylab = "", col = "darkgrey", yaxt='n', ...)
    axis(side = 2, labels = c(0, 5, 10), at = c(0, 5, 10))

    ## add estimated spikes
    for (spike_i in estimated_spikes) {
        segments(x0 = times[spike_i], x1 = times[spike_i], y0 = ye, y1 = ys, col = "blue", ...)
    }

    ## reset start and end vertical bounds for spike segments
    ys <- ys - incr
    ye <- ye - incr

    ## add true spikes
    for (spike_i in spike_times) {
        segments(x0 = spike_i, y0 = ye, x1 = spike_i, y1 = ys , col = "black", ...)
    }
}


plot_spikefinder <- function(test_set_results, dat, optimal_params, path_file_name) {
    l0_fit <- test_set_results$l0$fit
    l1_fit <- test_set_results$l1$fit

    pdf(path_file_name, height = 10, width = 10)
    par(mfrow = c(2, 1))
    plot(l0_fit)
    abline(v = which(dat$spikes > 0))
    abline(v = l0_fit$spikes, col = "red")

    plot(l1_fit$fittedValues)
    abline(v = which(dat$spikes > 0))
    abline(v = which(l1_fit$st > optimal_params$l1$threshold), col = "red")
    dev.off()

}

plot_whole_data <- function(dat_cell, path_file_name) {
    pdf(path_file_name, height = 10, width = 10)
    n_train <- length(dat_cell$train$trace)
    n_test <- length(dat_cell$test$trace)
    trace <- c(dat_cell$train$trace, dat_cell$test$trace)
    x <- 1:(n_train + n_test)
    color_transparent <- adjustcolor("grey", alpha.f = 0.3)

    plot(x, trace, cex = 0.2, pch = 20, col = ifelse(x > n_train, "darkblue", "darkred"))
    abline(v = which(c(dat_cell$train$spikes, dat_cell$test$spikes) > 0), col = color_transparent)
    dev.off()


}


plot.chen <- function(times, aibs, spikeTimes, l0Spikes = NULL, fitSpikes = NULL, xlim, pch=3, cex=0.5, l = 1, fitCol, xlab,
llbls) {
    ylim <- c(-8, 10)
    ## base plot of fl
    plot(times, aibs, cex = cex, pch = 20, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, xlab = xlab, ylab = "", xlim = xlim, col = "darkgrey", ylim = ylim, yaxt='n')
    axis(side = 2, labels = c(0, 5, 10), at = c(0, 5, 10))

    ## add oasis points
    i <- 1
    ys <- -0.2
    lineSize <- 1
    ye <- ys - lineSize
    buff <- 0.1
    incr <- lineSize + buff
    scaleOver <- 0.3 * (xlim[2] - xlim[1]) / 2
    for (thresh in fitSpikes) {
        for (spikey in thresh) {
            segments(x0 = times[spikey], x1 = times[spikey], y1 = ys, y0 = ye, col = fitCol[i],
            lwd = l)
        }
        i <- i + 1
        ys <- ys - incr
        ye <- ye - incr
    }

    ## add l0 points
    for (spikey in l0Spikes) {
        segments(x0 = times[spikey], x1 = times[spikey], y0 = ye, y1 = ys, col = "blue",
        lwd = l)
    }

    ys <- ys - incr
    ye <- ye - incr

    for (spikey in spikeTimes) {
        segments(x0 = spikey, y0 = ye, x1 = spikey, y1 = ys , col = "black", lwd = l)
    }
}

plot_error_vs_tuning_val <- function(df, thresholds, log_y = TRUE, ylab, colSpec, legend_txt) {
    pos = "topleft"
    if (log_y) {
        logscale = "y"
    } else {
        logscale = ""
        df = df %>% mutate(dist = 1 - dist)
    }

    n_thresholds <- length(thresholds)

    rgn <-range(df$dist)
    colSpec <- c('#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026')

    df %>% filter(method == "ar-fpop") %>%
        arrange(nSpikes) %$%
    plot(nSpikes, dist, type = 'b', cex = 1, pch = 20, log = logscale,
    ylim = rgn,
    col = "blue",
    xlab = "Number of Estimated Spikes",
    ylab = ylab, cex.axis = 1.5, cex.lab = 1.5,
    xlim = c(0, 300))

    dat <- df %>% filter(method == "oasis")
    for (k in 1:n_thresholds) {
        dat %>% filter(method == "oasis", threshold == thresholds[k]) %>%
            arrange(nSpikes) %$%
        lines(nSpikes, dist, type = 'b', cex = 0.5, col = colSpec[k], pch = 20)
    }
    legend(pos, lwd = 2, legend_txt, col = c(colSpec, "blue"))
}

plot_fig5 <- function(df, fit_l0, lam_l0, fitList, measure_i, thresholds, threshLabels, l1_type, output_dir) {

    # extra measure_i dataframe from stored results
    df <- df %>% filter(measure == measure_i)

    # fix formatting
    measure_names <- data.frame(measure = c("vanRossum", "VictorPurpura" ), name = c("van Rossum", "Victor-Purpura"))
    colSpec <- c('#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026')

    # update labels
    threshLabels <- c(threshLabels, TeX(sprintf("Soln (1.2): ($\\lambda = %.2f, L = %.2f)", lam_l0, 0)))
    threshLabels <- c(threshLabels, "True Spike Times")


    pdf(paste0(output_dir, "/example-", l1_type, "-", measure_i, ".pdf"), height = 6, width = 20)
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1,3), heights=c(1,1))

    # plot (tuning parameter value, error in terms of measure_i)
    if (measure_i != "downSampleCor") {
        plot_error_vs_tuning_val(df, thresholds, TRUE,
        ylab = paste0("Error in Spike Detection (", measure_names$name[measure_names$measure == measure_i], ")"),
        colSpec, threshLabels)
    } else {
        plot_error_vs_tuning_val(df, thresholds, FALSE, ylab = "Error in Spike Detection (1 - Correlation)",
        colSpec, threshLabels)
    }

    # plot estimated spikes based on optimal tuning values selected on training set
    timeFrame <- 1:length(test_data$trace) * (1 / fps)
    spikeTimes <- which(test_data$spikes != 0) * (1 / fps)
    plot.chen(timeFrame, test_data$trace, spikeTimes, l0Spikes = fit_l0$spikes,
    fitSpikes = fitList, xlim = c(0, 100), fitCol = colSpec, xlab = "Time (s)", llbls = threshLabels)

    dev.off()
}