# Python wrapper for OASIS calls
oasis <- function(y, gam, theta, type, tmp_dir) {
    prefix <- "/"
    if (file.exists(paste0(tmp_dir, "/y_dat.csv"))) {
        prefix <- paste0("/random_prefix-", sample.int(10^8, size = 1), "-")
    }

    # write y data to tmp file
    y_file <- paste0(tmp_dir, prefix, "y_dat.csv")
    out_file <- paste0(tmp_dir, prefix, "out_dat.csv")

    readr::write_csv(data.frame(y = y), y_file, col_names = FALSE)


    cmd_with_args <- paste("utils/run_oasis.sh", y_file, gam, theta, type, out_file ,sep = " ")
    oasis_out <- system(cmd_with_args, intern = TRUE, wait = TRUE, input = NULL)


    # read in tmp out file from oasis
    df_out <- readr::read_csv(out_file)
    spikes <- df_out$spikes
    calConc <- df_out$calcium

    if (type == "l1") {
        out <- list(spikes = spikes, fittedValues = calConc,
        dat = y, type = "oasis", changePts = NA,
        gam = gam,
        lambda = theta,
        cost = NA,
        nIntervals = NA,
        keepChgPts = NA,
        table = NA,
        calConc = calConc,
        type = type,
        st = spikes)
    } else {
        out <- list(spikes = spikes, fittedValues = calConc,
        dat = y, type = "oasis", changePts = NA,
        gam = gam,
        lambda = NA,
        s_min = theta,
        type = type,
        cost = NA,
        nIntervals = NA,
        keepChgPts = NA,
        table = NA,
        calConc = calConc,
        st = spikes)
    }

    class(out) <- "estimatedSpikes"


    # delete tmp files
    file.remove(y_file)
    file.remove(out_file)

    return(out)

}