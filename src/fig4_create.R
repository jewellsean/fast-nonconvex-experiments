 ## This is the code used to reproduce the two panel figure
## left panel is the van rossum curve for one dataset in the chen data
## right panel is the estimated spikes from each different method
required_packages <- c("tidyverse", "magrittr", "yaml")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)


library(tidyverse)
library(magrittr)
library(FastLZeroSpikeInference)
library(yaml)

source("utils/data_utils.R")
source("utils/plot_utils.R")
source("utils/full_path.R")
source("utils/oasis_utils.R")

# configs
config_file <- "../configs/configs_fig4.yml"
configs <- read_yaml(config_file)

# data parameters
dat_dir <- configs$data$dat_dir
fps <- configs$data$fps
dataset <- configs$data$data_set
cell_i <- configs$data$cell_i
train_subset <- 1:configs$data$train_subset_size
test_subset <- (configs$data$train_subset_size + 1):configs$data$test_subset_size

# model parameters
names <- readr::read_csv(paste0(dat_dir, "/data_names.txt"))
gam <- names$gamma[dataset]

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

min_cal_conc <- as.numeric(configs$optimization$l0$min_cal_conc)

# evaluation params
cost_params <- data.frame(fps = fps, time_scale = configs$evaluation$vr_timescale,
vp_cost = configs$evaluation$vp_cost, downsample = configs$evaluation$cor_downsample)

# output parameters
output_dir <- paste0(configs$output$figure_directory, "dataset_", dataset, "-cell_", cell_i)
dir.create(output_dir, showWarnings = FALSE)


# read data and split into train and test sets
dat <- read_spikefinder_data(dat_dir, dataset, split = FALSE)
train_test_dat <- dat_split_test_train(dat[[cell_i]], train_subset, test_subset)

train_data <- train_test_dat$train
test_data <- train_test_dat$test


# calculate best fits
l1_path <- list(type = "penalized", rng = l1_lams, thresholds = thresholds)
l0_path <- list(type = "constrained", min_lam = min_lam_l0, max_lam = max_lam_l0, n_lam = n_l0_lam, min_cal_conc = min_cal_conc)

df <- compare_l0_l1_solutions_over_path(train_data$trace, train_data$spikes, gam, l1_path, l0_path, cost_params)


unique_measures <- df %>% select(measure) %>% unique() %$% as.vector(measure)
for (measure_i in unique_measures) {
    ## Find optimal finds and plot
    if (measure_i != "downSampleCor") {
        d <- df %>% filter(measure == measure_i)
        optimal_params <- d %>% group_by(method, threshold) %>% summarize(dist = min(dist)) %>% left_join(d)
    } else {
        d <- df %>% filter(measure == measure_i)
        optimal_params <- df %>% filter(measure == measure_i) %>% group_by(method, threshold) %>% summarize(dist = max(dist)) %>% left_join(d)
    }

    opt_l0s <- optimal_params %>% ungroup() %>% filter(method == "ar-fpop") %>% select(lam)
    
    if (dim(opt_l0s)[[1]] > 1) {
      lam_l0 <- (opt_l0s %>% arrange(-lam) %>% as.matrix())[1]
    } else {
      lam_l0 <- opt_l0s %>% as.numeric()  
    }

    fit_l0 <- estimate_spikes(test_data$trace, gam, lam_l0, TRUE)

    fitList <- list()
    threshLabels <- c()
    for (i in 1:n_thresholds) {
        lam <- optimal_params %>% ungroup() %>% filter(method == "oasis", threshold == thresholds[i]) %>% select(lam)
        lam <- lam$lam[1]
        fit_oasis <- oasis(test_data$trace, gam, lam, "penalized", "test")
        fit_oasis$spikes <-  which(fit_oasis$st > thresholds[i])
        fitList[[i]] <- fit_oasis$spikes
        threshLabels <- c(threshLabels,
        TeX(sprintf("Soln (3.18): ($\\lambda = %.2f, L = %.2f)", lam, thresholds[i])) )

    }

    plot_fig5(df, fit_l0, lam_l0, fitList, measure_i, thresholds, threshLabels, "penalized", output_dir)



}
