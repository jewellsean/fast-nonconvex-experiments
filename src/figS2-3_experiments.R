## This script reproduces Fig 3 of Jewell et al. (2018)
library(FastLZeroSpikeInference)
library(LZeroSpikeInference)
library(yaml)


# CLI
args <- commandArgs(trailingOnly = TRUE)

# varying arguments from the CL
n = as.numeric(args[1])
poisMean = as.numeric(args[2])

# global configs
config_file <- args[3]
configs <- yaml.load_file(config_file)
gam <- configs$gam
sigma <- configs$sigma
min_cal_conc <- configs$min_cal_conc
number_simulations <- configs$number_simulations
lam <- configs$lam
op_threshold <- configs$op_threshold


## setup write directory and file
write_directory <- configs$local_exp_directory
dir.create(write_directory, showWarnings = FALSE)
file_name <- paste0("n_", n, "-poisMean_", poisMean, ".csv")


## store version info
info <- devtools::session_info()
fast_version <- info$packages$source[info$packages$package == "FastLZeroSpikeInference"]
slow_version <- info$packages$source[info$packages$package == "LZeroSpikeInference"]
fpop_version <- fast_version
pelt_version <- slow_version

times_df <- NULL


for (gam_i in gam) {
    for (lam_i in lam) {
        for (min_cal_conc_i in min_cal_conc) {
            for (seed_i in 1:number_simulations) {
                set.seed(seed_i)
                sim <- simulate_ar1(n = n, gam = gam_i, poisMean = poisMean, sd = sigma, seed = seed_i)

                if (n < op_threshold) {
                    tim1 <- system.time(fit1 <- estimateSpikes(sim$fl, gam = gam_i, lambda = lam_i, type = "ar1", hardThreshold = TRUE, calcFittedValues = FALSE))
                    tim2 <- system.time(fit2 <- estimateSpikes(sim$fl, gam = gam_i, lambda = lam_i, type = "ar1", hardThreshold = TRUE, pelt = FALSE, calcFittedValues = FALSE))
                    tim3 <- system.time(fit3 <- estimate_spikes(sim$fl, gam_i, lam_i, constraint = FALSE, estimate_calcium = FALSE, EPS = min_cal_conc_i))
                    tim4 <- system.time(fit4 <- estimate_spikes(sim$fl, gam_i, lam_i, constraint = TRUE, estimate_calcium = FALSE, EPS = min_cal_conc_i))

                    times <- c(tim1[[3]], tim2[[3]], tim3[[3]], tim4[[3]])
                    n_intervals <- c(max(fit1$nIntervals), max(fit2$nIntervals), max(fit3$n_intervals), max(fit4$n_intervals))

                    df <- data.frame(alg = c("pelt", "op", "fast-uncon", "fast-poscon"), timed = times,
                    n_intervals = n_intervals)
                } else {
                    tim1 <- system.time(fit1 <- estimateSpikes(sim$fl, gam = gam_i, lambda = lam_i, type = "ar1", hardThreshold = TRUE, calcFittedValues = FALSE))
                    tim3 <- system.time(fit3 <- estimate_spikes(sim$fl, gam_i, lam_i, constraint = FALSE, estimate_calcium = FALSE, EPS = min_cal_conc_i))
                    tim4 <- system.time(fit4 <- estimate_spikes(sim$fl, gam_i, lam_i, constraint = TRUE, estimate_calcium = FALSE, EPS = min_cal_conc_i))

                    times <- c(tim1[[3]], tim3[[3]], tim4[[3]])
                    n_intervals <- c(max(fit1$nIntervals), max(fit3$n_intervals), max(fit4$n_intervals))

                    df <- data.frame(alg = c("pelt", "fast-uncon", "fast-poscon"), timed = times, n_intervals = n_intervals)
                }

                df$seed <- seed_i
                df$gam <- gam_i
                df$poisMean <- poisMean
                df$lam <- lam_i
                df$sigma <- sigma
                df$n <- n
                df$min_cal_conc <- min_cal_conc_i

                times_df <- rbind(times_df, df)
                write.csv(times_df, file = paste0(write_directory, file_name), row.names = FALSE, quote = FALSE)
            }
        }
    }
}




