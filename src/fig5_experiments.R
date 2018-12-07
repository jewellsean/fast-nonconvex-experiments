## Arguments
## 1. temp directory for results
## 2. code directory to properly source everything
## 3. measure (VR, VP, Correlation)
## 4. Dataset
## 5. Data directory

## Dependencies
required_packages <- c("tidyverse", "magrittr", "yaml")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

require(tidyverse)
require(magrittr)
require(FastLZeroSpikeInference)
require(yaml)

source("utils/data_utils.R")
source("utils/oasis_utils.R")
source("utils/plot_utils.R")
source("utils/optimization_utils.R")
source("utils/spike_distance_utils.R")

args <- commandArgs(trailingOnly = TRUE)
measure <- args[1]
dataset <- as.numeric(args[2])
cell_i <-  as.numeric(args[3])

# global configs
config_file <-  args[5]
configs <- yaml.load_file(config_file)

data_dir <- configs$data$dat_dir
if (args[4] == "constrained") {
    l0_constraint <- TRUE
} else {
    l0_constraint <- FALSE
}

# data parameters
data_params <- data.frame(measure = measure,
                    dataset = dataset,
                    cell_i = cell_i, 
                    fps = configs$data$fps)

# optimization parameters
n_l1_lam <- configs$optimization$grid_size
min_lam_log <- configs$optimization$l1$min_lam_log
max_lam_log <- configs$optimization$l1$max_lam_log
l1_lams <- 10 ^ seq(min_lam_log, max_lam_log, length.out = n_l1_lam)

n_thresholds <- configs$optimization$l1$n_threshold
min_threshold <- as.numeric(configs$optimization$l1$min_threshold)
max_threshold <- configs$optimization$l1$max_threshold
stopifnot(min_threshold > 0)
thresholds <- seq(min_threshold, max_threshold, length.out = n_thresholds)

n_l0_lam <- configs$optimization$grid_size
min_lam_l0 <- configs$optimization$l0$min_lam
max_lam_l0 <- configs$optimization$l0$max_lam
l0_lams <- seq(min_lam_l0, max_lam_l0, length.out = n_l0_lam)
min_cal_conc <- as.numeric(configs$optimization$l0$min_cal_conc)

optim_params <- list( l0_constraint = l0_constraint,
                            l0_min_cal_conc = min_cal_conc,
                            l0_lams = l0_lams,
                            l1_lams = l1_lams,
                            l1_thresholds = thresholds)

# evaluation params
cost_params <- data.frame(fps = configs$data$fps, time_scale = configs$evaluation$vr_timescale,
vp_cost = configs$evaluation$vp_cost, downsample = configs$evaluation$cor_downsample)

## Set gamma parameter
names <- readr::read_csv(paste0(data_dir, "/data_names.txt"))
data_params$gamma <- names$gamma[dataset]

## 1. Read data
## Returns a list of datasets for each cell in selected dataset
## Each element of the list is itself a list of training and test data
## within each of training a test there is the FL trace and true spikes
## Structure
## Cells
## ---> Training
##  ---> Fl trace
##  ---> Spikes
## ---> Test
##  ---> Fl trace
##  ---> Spikes
dat <- read_spikefinder_data(data_dir, dataset)

## For each cell in data, preform the following operations, based on input measure,
## 1. Optimize parameters on training set
## 2. Report the value of the measure, at the optimal parameters, on the test set
## 3. Plot fits on training and test data
## 4. Save spike magnitudes as raw (easier for post transformations)

n_cells <- length(dat)
if (cell_i <= n_cells) {
    dat_cell <- dat[[cell_i]]
    train_dat <- dat_cell$train
    test_dat <- dat_cell$test

    ## 1. Optimize parameters on training set
    optimal_params <- optimize_parameters(train_dat, data_params, optim_params, cost_params)

    ## For comarison, optimize on test set, too
    optimal_params_test <- optimize_parameters(test_dat, data_params, optim_params, cost_params)

    file_name_base <- paste0(configs$output$file_directory, "dataset-", dataset, "-cell_i-", cell_i, "-measure-", measure, "-constraint-", l0_constraint)

    ## 2. Report the value of the measure, at the optimal parameters, on the test set
    test_set_results <- evaluate_performance(test_dat, optimal_params, optimal_params_test, data_params, optim_params, cost_params)
    df_out <- cbind(test_set_results$out, data_params, cost_params, data.frame(constrained = optim_params$l0_constraint, min_cal_conc = optim_params$l0_min_cal_conc))

    write_directory <- paste0(configs$output$file_directory, "results/")
    dir.create(write_directory, showWarnings = FALSE)
    write_out_filename <- paste0(write_directory, "dataset-", dataset, "-cell_i-", cell_i, "-measure-", measure, "-constraint-", l0_constraint)
    write_csv(df_out, path = paste0(write_out_filename, "-opt_test_results_", ".csv"))

    ## 3. Plot fits test data
    plot_sub_dir <- paste0(configs$output$file_directory, "plots/")
    dir.create(plot_sub_dir, showWarnings = FALSE)
    file_plot_base <- paste0(plot_sub_dir, "dataset-", dataset, "-cell_i-", cell_i, "-measure-", measure, "-constraint-", l0_constraint)
    plot_spikefinder(test_set_results, test_dat, optimal_params, path = paste0(file_plot_base, "-test_results_", ".pdf"))

    ## 4. Save spike magnitudes
    magnitude_sub_dir <- paste0(configs$output$file_directory, "magnitudes/")
    dir.create(magnitude_sub_dir, showWarnings = FALSE)
    file_magnitude_base <- paste0(magnitude_sub_dir, "dataset-", dataset, "-cell_i-", cell_i, "-measure-", measure, "-constraint-", l0_constraint)
    save_spike_magnitudes(test_set_results$l0$fit, test_dat, path = paste0(file_magnitude_base, "-spike_mag_", ".csv"))


    ## 5. Plot full data series
    plot_whole_data(dat_cell, path = paste0(file_plot_base, "-whole_data_", ".pdf"))
}


