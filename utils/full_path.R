required_packages <- c("dplyr", "magrittr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(magrittr)
library(dplyr)
library(FastLZeroSpikeInference)

source("utils/spike_distance_utils.R")

l1_solution_path <- function(cal_conc, true_spikes, gam, l1_path, params) {
    n_thresholds <- length(thresholds)
    n_l1_lam <- length(l1_path$rng)
    df <- NULL
    for (i in 1:n_l1_lam) {
        fit <- oasis(cal_conc, gam, l1_path$rng[i], l1_path$type, "test")
        for (j in 1:n_thresholds) {
                estimated_spikes <- list(spike_times = which(fit$st > thresholds[j]) * (1 / fps),
                spike_locations = fit$st * (fit$st > thresholds[j]))
                df <- rbind(df,
                data.frame(calcDist(true_spikes, estimated_spikes, params),
                method = "oasis", lam = l1_path$rng[i], threshold = thresholds[j], constrained = NA, min_cal_conc = NA))
        }
    }
    return(df)
}

l0_solution_path <- function(cal_conc, true_spikes, gam, l0_path, params) {
    constrained <- (l0_path$type == "constrained")
    min_cal_conc <- l0_path$min_cal_conc

    df <- NULL

    # lams
    lams <- seq(l0_path$min_lam, l0_path$max_lam, length.out = l0_path$n_lam)
    n_l0_lam <- length(lams)

    for (i in 1:n_l0_lam) {
        fit <- estimate_spikes(cal_conc, gam, lams[i], constrained, TRUE, min_cal_conc)
        fitted_spike_locations <- 0 * numeric(length(cal_conc))
        st <- get_spike_magnitudes(fit)
        fitted_spike_locations[fit$spikes] <- st

        estimated_spikes <- list(spike_times = fit$spikes * (1 / fps),
        spike_locations = fitted_spike_locations)

        df <- rbind(df,
        data.frame(calcDist(true_spikes, estimated_spikes, params), method = "ar-fpop",
        constrained = constrained, lam = lams[i], threshold = NA, min_cal_conc = min_cal_conc))
    }
    return(df)
}

## calcium_d preprocessed calcium data
## spike_locations array of spikes where spike_locations[i] denote the # of spikes at timepoint i
## eg. 0 0 1 2 0 0 2
compare_l0_l1_solutions_over_path <- function(cal_conc, spike_locations, gam, l1_path, l0_path, params) {

    stopifnot(length(spike_locations) == length(cal_conc))

    # setup ground truth for distance calculations
    spike_times <- which(spike_locations > 0) * (1 / fps)
    true_spikes <- list(spike_times = spike_times, spike_locations = spike_locations)


    # l1_path: contains information about the l1 path
    # rng - range of values for tuning values
    # thresholds - number of thresholds (only for penalized problem)
    df_l1 <- l1_solution_path(cal_conc, true_spikes, gam, l1_path, params) # penalized version


    # l0_path: contains information about the l0 path
    # type - constrained or unconstrained
    # rng - range of tuning parameters lambda
    # min_cal_conc - minimum calcium concentration in optimization problem
    df_l0 <- l0_solution_path(cal_conc, true_spikes, gam, l0_path, params)

   return(rbind(df_l1, df_l0))

}
