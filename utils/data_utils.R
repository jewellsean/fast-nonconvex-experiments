## Utilities to create Figures and Tables from
## Jewell et al. (2018) arXiv:1802.07380
required_packages <- c("dplyr", "magrittr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

require(magrittr)
require(dplyr)

col_renamer <- function(x) {return(paste0("x_", x))}

subset_data <- function(calcium_dat, spike_dat, subset_is, cell_i) {
  calcium_d <- calcium_dat[subset_is , cell_i] %>% as.matrix() %>% as.numeric()
  calcium_d <- calcium_d[complete.cases(calcium_d)]
  spikes <- spike_dat[subset_is , cell_i] %>% as.matrix() %>% as.numeric()
  spikes <- spikes[complete.cases(spikes)]  
  return(list(calcium_d = calcium_d, spikes = spikes))
}

read_csvs_in_directory <- function(dir) {
  df <- NULL
  experiment_files <- list.files(dir, pattern = ".csv")

  if (is.null(experiment_files)) {
    stop(paste0('This merges csv files directory: ', dir,
    ", but there are no files currently in the directory."))
  }

  for (experiment_file in experiment_files) {
    file_i <- paste0(dir, "/", experiment_file)
    df <- rbind(df, read.csv(file_i))
  }
  return(df)
}


read_spikefinder_data <- function(data_dir, dataset, split = TRUE) {
  data_str <- paste0(data_dir, "/spikefinder.train/", dataset, ".train.")

  calcium_dat <- read.csv(paste0(data_str, "calcium.csv"))
  spike_dat <- read.csv(paste0(data_str, "spikes.csv"))

  calcium_dat %<>% rename_all(col_renamer)
  spike_dat %<>% rename_all(col_renamer)

  dat <- list()
  n_cells <- dim(calcium_dat)[[2]]
  for (cell_i in 1:n_cells) {
    dat_slice <- list(trace = calcium_dat[ ,cell_i], spikes = spike_dat[ ,cell_i])
    dat_slice <- dat_clean(dat_slice)
    if (split) {
        dat[[cell_i]] <- dat_split_test_train(dat_slice)
    } else {
        dat[[cell_i]] <- dat_slice
    }
  }
  return(dat)
}

## clean data by taking only complete cases
dat_clean <- function(dat) {
  return(list(trace =  dat$trace[complete.cases(dat$trace)],
  spikes =  dat$spikes[complete.cases(dat$spikes)]))
}

## dat is an array
dat_split_test_train <- function(dat, train_ind = NULL, test_ind = NULL) {
  stopifnot(length(dat$trace) == length(dat$spikes))
  if (is.null(train_ind)) {
      n <- length(dat$trace)
      ind <- 1:n
      train_ind <- ind <= floor(n / 2)
      test_ind <- ind >= floor(n / 2) + 1
  }
  return(list(train = list(trace =  dat$trace[train_ind], spikes =  dat$spikes[train_ind]),
  test = list(trace =  dat$trace[test_ind], spikes =  dat$spikes[test_ind])))
}

save_spike_magnitudes <- function(l0_fit, test_dat, path) {

  len <- length(l0_fit$estimated_calcium)

  spike_mag <- l0_fit$estimated_calcium[2:len] - l0_fit$gam * l0_fit$estimated_calcium[1:(len - 1)]

  spike_loc <- numeric(len)
  spike_loc[l0_fit$spikes] <- rep(1, length(l0_fit$spikes))

  out <- data.frame(spike_mag = c(0, spike_mag),
  spike_loc = spike_loc,
  true_spikes = test_dat$spikes)

  write_csv(out, path = path)
}
