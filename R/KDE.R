KDE <- function(spike_times) {
    sort(spike_time)    # sort
    max_value <- max(spike_times)
    min_value <- min(spike_times)
    T <- max_value - min_value

    diff_spike <- sort(diff(spike_times))
    dt_samp <- diff_spike[which(diff_spike != 0)][1]

    tin <- seq(min_value, max_value, length=min(ceiling(T / dt_samp), 1e3))
    spike_ab <- spike_times[which((spike_times >= min(tin)) & (spike_times <= max(tin)))]

    dt <- min(diff(tin))

    y_hist <- hist(x=spike_ab, plot=FALSE, breaks=(c(tin, max_value)))$counts
    L <- length(y_hist)
    N <- sum(y_hist)
    y_hist <- y_hist / (N * dt)

    Wmin <- 2 * dt
    Wmax <- 1 * (max_value - min_value)

    tol <- 1e-5
    phi <- (sqrt(5) + 1) / 2

    a <- iLogExp(Wmin)
    b <- iLogExp(Wmax)

    c1 <- (phi - 1) * a + (2 - phi) * b
    c2 <- (2 - phi) * a + (phi - 1) * b

    f1 <- costFunction(y_hist, N, logExp(c1), dt)[[1]]
    f2 <- costFunction(y_hist, N, logExp(c2), dt)[[1]]

    k <- 1
    W <- numeric(20)
    C <- numeric(20)

    # repeat 20 times if c1+c2 < (difference between a and b)

    while((abs(b - a) > (tol * (abs(c1) + abs(c2)))) & (k <= 20)) {
        if (f1 < f2) {
            b <- c2
            c2 <- c1

            c1 <- (phi - 1) * a + (2 - phi) * b

            f2 <- f1
            res <- costFunction(y_hist, N, logExp(c1), dt)
            f1 <- res[[1]]
            yh1 <- res[[2]]

            W[k] <- logExp(c1)
            C[k] <- f1
            optw <- logExp(c1)
            y <- yh1 / sum(yh1 * dt)
        } else {
            a <- c1
            c1 <- c2

            c2 <- (2 - phi) * a + (phi - 1) * b

            f1 <- f2
            res <- costFunction(y_hist, N, logExp(c2), dt)
            f2 <- res[[1]]
            yh2 <- res[[2]]

            W[k] <- logExp(c2)
            C[k] <- f2
            optw <- logExp(c2)
            y <- yh2 / sum(yh2 * dt)
        }

        k <- k + 1
    }

    y <- y * length(spike_times)

    drawKDE(y, tin)

    return (list(y, tin, optw))
}

logExp <- function(x) {
    if (x < 1e2) {
        return (log(1 + exp(x)))
    } else {
        return (x)
    }
}

iLogExp <- function(x) {
    if (x < 1e2) {
        return (log(exp(x) - 1))
    } else {
        return (x)
    }
}

costFunction <- function(y_hist, N, w, dt) {
    yh <- fftKernel(y_hist, w / dt) # density

    # formula for density
    C <- sum(yh * yh) * dt - 2 * sum(yh * y_hist) * dt + 2 * 1 / (sqrt(2 * pi) * w * N)
    C <- C * N * N

    return (list(C, yh))
}

fftKernel <- function(x, w) {
    # y = fftkernel(x, w)
    # 
    # Function `fftkernel' applies the Gauss kernel smoother to an input signal using FFT algorithm.
    #
    # Input argument
    # x : Sample signal vector
    # w : Kernel bandwidth (the standard deviation) in unit of the sampling resolution of x.
    # Output argument
    # y : Smoothed signal.
    #
    # MAY 5 / 23, 2012 Author Hideaki Shimazaki
    # RIKEN Brain Science Institute
    # http://2000.jukuin.keio.ac.jp/shimazaki
    # 
    # (New correction in version 1)
    # y-axis was multiplied by the number of data, so that
    # y is a time histogram representing the density of spikes.

    L <- length(x)
    Lmax <- max(1.0, floor(L + 3.0 * w))
    n <- floor(2^(nextPow2(Lmax)))

    X <- fft(c(x, numeric(n - L)))

    f <- (0:(n-1)) / n
    f <- c(-f[1:(floor(n / 2) + 1)], f[(floor(n / 2) + 1):3])

    K <- exp(-0.5 * ((w * 2 * pi * f)^2))

    y <- fft(X * K, inverse=TRUE) / n

    return (Re(y[1:L]))
}

nextPow2 <- function(n) {
    if(n < 0) {
        return(0)
    } else {
        m <- ceiling(log2(n))
        return(m)
    }
}

drawKDE <- function(y, t) {
    plot(t, y, type="l", ylim=range(0, y))
}

    
    
