## spike_distances
## code for the spike distances, slightly refactored for server runs
## Spike distance metrics
calculate_spike_distance <- function(x, y, measure, cost_params) {
    if (measure == "VP") {
        return(victorPurpuraDist(x, y, cost_params$vp_cost))
    } else if (measure == "VR") {
        return(vanRossumDist(x, y, cost_params$time_scale))
    } else if (measure == "Cor") {
        return(corr_metric(x, y, cost_params$downsample))
    } else {
        stop("Could not find specified metric")
    }
}

calcDist <- function(true_spikes, estimated_spikes, params) {
    df <- NULL
    nSpikes <- length(estimated_spikes$spike_times)
    # vanRossum
    # requires spike times
    df <- rbind(df,
    data.frame(measure = "vanRossum",
    dist = vanRossumDist(true_spikes$spike_times,
    estimated_spikes$spike_times,
    params$time_scale),
    nSpikes = nSpikes))

    # victorPurpuraDist
    df <- rbind(df,
    data.frame(measure = "VictorPurpura",
    dist = victorPurpuraDist(true_spikes$spike_times,
    estimated_spikes$spike_times,
    params$vp_cost),
    nSpikes = nSpikes))

    # downsampled correlation measure
    df <- rbind(df,
    data.frame(measure = "downSampleCor",
    dist = corr_metric(true_spikes$spike_locations,
    estimated_spikes$spike_locations,
    params$downsample),
    nSpikes = nSpikes))

    return(df)
}

## http://www-users.med.cornell.edu/~jdvicto/spkdm.html
victorPurpuraDist <- function(tli,tlj,cost) {

    # d=spkd(tli,tlj,cost) calculates the "spike time" distance
    # (Victor & Purpura 1996) for a single cost
    #
    # tli: vector of spike times for first spike train
    # tlj: vector of spike times for second spike train
    # cost: cost per unit time to move a spike
    #
    # Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
    # Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
    # Translated to R by Sean Jewell from Matlab code by Daniel Reich

    nspi=length(tli)
    nspj=length(tlj)

    if (cost==0) {
        return(abs(nspi-nspj))
    }
    if (cost==Inf) {
        return(nspi+nspj)
    }

    scr <- matrix(nrow = nspi+1, ncol = nspj+1)
    #
    #     INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
    #
    scr[ ,1] <- c(0:nspi)
    scr[1, ] <- c(0:nspj)

    if (nspi & nspj) {
        for (i in 2:(nspi + 1)) {
            for (j in 2:(nspj+1)) {
                scr[i,j] =
                min(c(scr[i - 1, j] + 1,
                scr[i, j - 1] + 1,
                scr[i - 1, j - 1] + cost * abs(tli[i - 1] - tlj[j - 1])))
            }
        }

    }
    d <- scr[nspi + 1, nspj + 1]
    return(d)
}

vanRossumDist <- function(u, v, tau) {
    lenU <- length(u)
    lenV <- length(v)
    d2 <- 0
    s1 <- 0
    for (i in 1:lenU) {
        for (j in 1:lenU) {
            s1 <- s1 + exp(-abs(u[i] - u[j]) / tau)
        }
    }
    s3 <- 0
    for (i in 1:lenU) {
        for (j in 1:lenV) {
            s3 <- s3 - 2 * exp(-abs(u[i] - v[j]) / tau)
        }
    }
    s2 <- 0
    for (i in 1:lenV) {
        for (j in 1:lenV) {
            s2 <- s2 + exp(-abs(v[i] - v[j]) / tau)
        }
    }
    return(sum(c(s1, s2, s3)))
}

downsample <- function(y, factor) {
    x <- rep(1, factor)
    out <- convolve(y, rev(x), type = 'filter')
    return(out[seq(1, by = factor, to = length(out))])
}

corr_metric <- function(times_x, times_y, factor) {
    times_x_down <- downsample(times_x, factor)
    times_y_down <- downsample(times_y, factor)
    if (sum(times_x_down) == 0 || sum(times_y_down) == 0) {
        return(0)
    }
    return(cor(times_x_down, times_y_down))
}

get_spike_magnitudes <- function(fit) {
    gam <- fit$gam
    st <- fit$estimated_calcium[fit$spikes] - gam * fit$estimated_calcium[fit$spikes - 1]
    return(st)
}