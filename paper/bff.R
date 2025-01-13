## ----"main-setup", include = FALSE--------------------------------------------
## knitr options
library(knitr)
opts_chunk$set(fig.height = 4,
               echo = FALSE,
               warning = FALSE,
               message = FALSE,
               cache = TRUE,
               eval = TRUE)

## should sessionInfo be printed at the end?
Reproducibility <- TRUE

## packages
library(metadat) # data set from Bartos (2023)
library(metabf) # BFFs for meta-analysis (not on CRAN, install from GitHub repo)
library(brms) # Bayesian regression models with MCMC
library(INLA) # Bayesian regression models with INLA
library(ggplot2) # plotting
library(haven) # handling SPSS imported data set


## ----"pFun-and-BFFun", echo = FALSE, fig.height = 3.25, fig.align = "center"----
## Create cartoon that shows p-value function and support curve
par("mar" = c(2.75, 3.5, 1.5, 0.5), mfrow = c(1, 2))
null <- seq(-2.5, 2.5, length.out = 500)
yi <- 0
sei <- 0.8
t0 <- -0.5

## p-value function
pval <- function(yi, sei, null) {2*pnorm(q = abs(yi - null)/sei, lower.tail = FALSE)}
p <- pval(yi, sei, null)
p0 <- pval(yi, sei, t0)
alpha <- 0.2
ci <- yi + c(-1, 1)*sei*qnorm(p = 1 - alpha/2)
plot(null, p, xaxt = "n", yaxt = "n", las = 1, type = "n", bty = "n",
xlab = "", ylab = "",
main = bquote(italic(P) * "-value function"))
axis(side = 1, at = c(-3, 3, t0), labels = c("", "", ""))
mtext(text = expression(theta[0]), side = 1, at = t0, line = 0.5, cex = 0.8)
axis(side = 2, at = c(alpha, 0, 1), labels = c(expression(alpha), 0, 1), las = 1)
mtext(text = "Tested parameter value", side = 1, line = 1.5)
mtext(text = bquote(italic(P) * "-value"), side = 2, line = 1.5)
abline(h = alpha, lty = 2, col = "darkgrey")
arrows(x0 = ci[1], x1 = ci[2], y0 = alpha, lwd = 2.5,
code = 3, length = 0.05, angle = 90, col = "darkgrey")
arrows(x0 = yi + 0.6, x1 = yi + 0.1, y1 = 1, y0 = 1, length = 0.05)
points(x = yi, y = 1, col = "darkgrey", pch = 20)
text(x = yi + 0.45, y = 1, label = "Point estimate", cex = 0.7, pos = 4)
arrows(x0 = yi,  y1 = alpha - 0.02, y0 = alpha - 0.13, length = 0.05)
text(x = yi, y = alpha - 0.17, label = bquote((1 - alpha) * 100 * "% confidence interval"),
     cex = 0.7)
segments(x0 = t0, y0 = -1, y1 = p0, col = "darkgrey", lty = 3, lwd = 1.5)
segments(x0 = -100, x1 = t0, y0 = p0, col = "darkgrey", lty = 3, lwd = 1.5)
arrows(x0 = -2.6 + 0.4, x1 = -2.6 + 0.05, y1 = p0 + 0.05, y0 = p0 + 0.15, length = 0.05)
text(x = -2.8, y = p0 + 0.18, label = bquote(italic(P) * "-value for" ~ theta[0]),
     cex = 0.7, pos = 4)
lines(null, p, lwd = 2)

## support curve
bffun <- function(yi, sei, muPrior, sdPrior, null) {
    sqrt(1 +sdPrior^2/sei^2)*
        exp(-0.5*((yi - null)^2/sei^2 - (yi - muPrior)^2/(sdPrior^2 + sei^2)))
}
sdPrior <- 2
muPrior <- 0.5
bf <- bffun(yi, sei, muPrior, sdPrior, null)
bflogbks <- c(1/100, 1/30, 1/10, 1/3, 1, 3, 10)
bfloglabs <- c("1/100", "1/30", "1/10", "1/3", "1", "3", "10")
bf0 <- bffun(yi, sei, muPrior, sdPrior, t0)
k <- 0.5
kMEE <- bffun(yi, sei, muPrior, sdPrior, yi)
si <- yi + c(-1, 1)*sei*sqrt(log(1 + sdPrior^2/sei^2) +
                             (yi - muPrior)^2/(sei^2 + sdPrior^2) - 2*log(k))
plot(null, bf, xaxt = "n", yaxt = "n", las = 1, type = "n", bty = "n",
xlab = "", ylab = "",
main = bquote("Support curve" * " "))
axis(side = 1, at = c(-3, 3, t0), labels = c("", "", ""))
axis(side = 2, at = c(0, 100, k), labels = c(0, "", expression(italic(k))), las = 1)
mtext(text = expression(atop(atop(infinity, " "  %up%  " "), atop(" ", " "))),
      las = 1, side = 2, at = kMEE*0.95, line = 0, cex = 1.5)
mtext(text = expression(theta[0]), side = 1, at = t0, line = 0.5, cex = 0.8)
mtext(text = "Tested parameter value", side = 1, line = 1.5)
mtext(text = "Bayes factor", side = 2, line = 1.5)
abline(h = k, lty = 2, col = "darkgrey")
arrows(x0 = si[1], x1 = si[2], y0 = k, lwd = 2.5,
code = 3, length = 0.05, angle = 90, col = "darkgrey")
arrows(x0 = yi + 0.6, x1 = yi + 0.1, y1 = kMEE + 0.02, y0 = kMEE + 0.02, length = 0.05)
points(x = yi, y = kMEE, col = "darkgrey", pch = 20)
text(x = yi + 0.45, y = kMEE, label = "Point estimate", cex = 0.7, pos = 4)
arrows(x0 = yi,  y1 = k - 0.05, y0 = k - 0.3, length = 0.05)
text(x = yi, y = k - 0.4, label = bquote(italic(k) ~ "support interval"), cex = 0.7)
segments(x0 = t0, y0 = -1, y1 = bf0, col = "darkgrey", lty = 3, lwd = 1.5)
segments(x0 = -100, x1 = t0, y0 = bf0, col = "darkgrey", lty = 3, lwd = 1.5)
arrows(x0 = -2.6 + 0.4, x1 = -2.6 + 0.05, y1 = bf0 - 0.1, y0 = bf0 - 0.3, length = 0.05)
text(x = -2.8, y = bf0 - 0.4, label = bquote("Bayes factor for" ~ theta[0]),
     cex = 0.7, pos = 4)
lines(null, bf, lwd = 2)


## ----"RECOVERY-example"-------------------------------------------------------
## Baricitinib in patients admitted to hospital with COVID-19 (RECOVERY): a
## randomised, controlled, open-label, platform trial and updated meta-analysis
## https://doi.org/10.1016/S0140-6736(22)01109-6

## Overall, 514 (12%) of 4148 patients allocated to baricitinib versus 546 (14%)
## of 4008 patients allocated to usual care died within 28 days (age-adjusted
## rate ratio 0.87; 95% CI 0.77–0.99; p=0.028).
HR <- 0.87
ciHR <- c(0.77, 0.99)
p <- 0.028
y <- log(HR)
se <- diff(log(ciHR))/2/qnorm(p = 0.975)
se2 <- se^2
## This 13% proportional reduction in mortality was somewhat smaller than that
## seen in a meta-analysis of eight previous trials of a JAK inhibitor
## (involving 3732 patients and 425 deaths), in which allocation to a JAK
## inhibitor was associated with a 43% proportional reduction in mortality (rate
## ratio 0.57; 95% CI 0.45–0.72)
HRprior <- 0.57
ciHRprior <- c(0.45, 0.72)
m <- log(HRprior)
s <- diff(log(ciHRprior))/2/qnorm(p = 0.975)
v <- s^2
t0seq <- seq(-0.5, 0.1, 0.01)
delta <- 0.1

k <- 1
## normal prior
bf1 <- function(t0) sqrt(1 + v/se2)*exp(-0.5*((y - t0)^2/se2 - (y - m)^2/(se2 + v)))
bf1b <- function(t0) sqrt(1 + v/se2)*exp(-0.5*((y - t0)^2/se2))
## local normal prior
bf2 <- function(t0) sqrt(1 + v/se2)*exp(-0.5*((y - t0)^2/(se2*(1 + se2/v))))
## pathological prior
bf3 <- function(t0) exp((-delta*(y - t0) + delta^2*0.5)/se2)
bfs <- sapply(X = list(bf1, bf1b, bf2, bf3), FUN = function(f) f(t0seq))

## SI for normal prior
si1 <- function(k) y + c(-1, 1)*sqrt(se2)*sqrt(log(1 + v/se2) + (y - m)^2/(se2 + v) - log(k^2))
si1b <- function(k) y + c(-1, 1)*sqrt(se2)*sqrt(log(1 + v/se2) - log(k^2))
## SI for local normal prior
si2 <- function(k) y + c(-1, 1)*sqrt(se2)*sqrt((log(1 + v/se2) - log(k^2))*(1 + se2/v))
## SI for pathological prior
si3 <- function(k) c(y + se2*log(k)/delta - delta/2, 1e+100) #Inf)
sis <- t(sapply(X = list(si1, si1b, si2, si3), FUN = function(f) f(k)))


## ----"normal-mean-example", fig.height = 4------------------------------------
colors <- palette.colors(n = 5, palette = "Okabe-Ito", alpha = 0.9)[2:5]
bfbks <- c(1/1000, 1/100, 1/10, 1, 10, 100, 1000)
bflabs <- c("1/1000", "1/100", "1/10", "1", "10", "100", "1000")
par(mar = c(4, 5, 1, 1.5))
matplot(t0seq, bfs, type = "l", lty = 1, col = colors, lwd = 1.5, las = 1,
        xlab = bquote("Log hazard ratio" ~ theta), ylab = "Bayes factor",
        log = "y", yaxt = "n", ylim = c(1/1000, 1000),
        panel.first = graphics::grid(lty = 3, equilogs = FALSE))
axis(side = 2, at = bfbks, labels = bflabs, las = 1)
abline(h = 1, lty = 2, col = adjustcolor(col = "black", alpha.f = 0.3))
arrows(x0 = sis[,1], x1 = sis[,2], y0 = k + c(-0.0025, 0.0025, 0, 0), col = colors,
       code = 3, length = 0.075, angle = 90)
points(x = c(y, y, y), y = c(bf1(y), bf1b(y), bf2(y)), col = colors[1:3], pch = 20)
legend("topleft", lty = 1, col = colors, lwd = 1.5, cex = 0.7,
       title = expression("Prior under"~ italic(H)[1]),
       legend = c(bquote(theta ~ "~ N(" * .(round(m, 2)) * "," ~ .(round(sqrt(v), 2))^2 * ")"),
                  bquote(theta ~ "~ N(" * .(round(y, 2)) * "," ~ .(round(sqrt(v), 2))^2 * ")"),
                  bquote(theta ~ "~ N(" * theta[0] * "," ~ .(round(sqrt(v), 2))^2 * ")"),
                  bquote(theta == theta[0] + .(delta))),
       bg = "white")
mtext(text = bquote("" %->% "Harm"), side = 1, line = 0, at = 0.045, cex = 0.9)
mtext(text = bquote("Benefit" %<-% ""), side = 1, line = 0, at = -0.05, cex = 0.9)
arrows(x0 = min(t0seq), y0 = c(1.5, 1/1.5), y1 = c(5, 1/5),
       col = adjustcolor("black", alpha.f = 0.8), length = 0.05)
text(x = min(t0seq), y = c(2, 1/2),
     labels = c(expression("Support for" ~ theta),
                expression("Support for" ~ italic(H)[1])),
     pos = 4, cex = 0.75, col = adjustcolor("black", alpha.f = 0.8))


## ----"distribution-k"---------------------------------------------------------
## ## generate data under H1, what is the distribution of k?
## ## the theoretical density
## fk <- function(k, v, se) sqrt(2/pi*(1 + v/se^2)/k^4/log(k^2/(1 + v/se^2)))

## set.seed(42)
## ysim <- rnorm(n = 100000, mean = m, sd = sqrt(v + se^2))
## ksim <- bffun(yi = ysim, sei = se, muPrior = m, sdPrior = sqrt(v), null = ysim)

## hist(ksim, probability = TRUE, xlim = c(1, 20),
##      breaks = seq(0, 100000, 0.1))

## kseq <- seq(0, 100, 0.01)
## lines(kseq, fk(kseq), lty = 2, lwd = 2)


## ----"BF-asymptotic", fig.height = 4.5, fig.width = 8-------------------------
## BFF for normal mean
bff <- function(t0, y, kappa, n, m, v) {
    sqrt(1 + v*n/kappa^2)*
        exp(-0.5*((y - t0)^2*n/kappa^2 - (y - m)^2/(kappa^2/n + v)))
}

## CDF of the BF
pbf. <- function(q, t0, kappa, n, m, v, tstar) {
    lambda <- (tstar - (t0 - m)*kappa^2/n/v - t0)^2*n/kappa^2
    A <- log(1 + n*v/kappa^2) + (t0 - m)^2/v - 2*log(q)
    if (A <= 0) p <- 1
    else {
        p <- pchisq(q = A*(1 + kappa^2/v/n), df = 1,
                    ncp = lambda, lower.tail = FALSE)
    }
    return(p)
}
pbf <- Vectorize(FUN = pbf.)

## quantile function of the BF
qbf. <- function(p, t0, kappa, n, m, v, tstar) {
    upper <- bff(t0 = t0, y = t0, kappa = kappa, n = n, m = m, v = v)
    if (p == 0) {
        q <- 0
    } else if (p == 1) {
        q <- upper
    } else {
        rootFun <- function(q) {
            pbf(q = q, t0 = t0, kappa = kappa, n = n, m = m, v = v,
                tstar = tstar) - p
        }
        res <- try(uniroot(f = rootFun,
                           interval = c(.Machine$double.eps, upper),
                           tol = .Machine$double.eps^0.75)$root,
                   silent = TRUE)
        if (inherits(res, "try-error")) q <- NaN
        else q <- res
    }
    return(q)
}
qbf <- Vectorize(FUN = qbf.)

## different sample sizes
ns <- c(50, 500, 5000)
cols <- hcl.colors(n = length(ns), palette = "viridis", alpha = 0.9,
                   rev = TRUE)
cols2 <- adjustcolor(col = cols, alpha.f = 0.4)
cols3 <- adjustcolor(col = cols, alpha.f = 0.1)
t0seq <- seq(-1, 1, length.out = 500)
tstar <- 0
kappa <- 2
v <- 4
p <- 0.5
plower <- 0.05
pupper <- 0.95
bfmedians <- sapply(X = ns, FUN = function(n) {
    qbf(p = p, t0 = t0seq, kappa = kappa, n = n, m = t0seq, v = v,
        tstar = tstar)
})
bfupper <- sapply(X = ns, FUN = function(n) {
    qbf(p = pupper, t0 = t0seq, kappa = kappa, n = n, m = t0seq, v = v,
        tstar = tstar)
})
bflower <- sapply(X = ns, FUN = function(n) {
    qbf(p = plower, t0 = t0seq, kappa = kappa, n = n, m = t0seq, v = v,
        tstar = tstar)
})

par(mar = c(4, 5, 1, 1.5))
bfbks <- c(1/1000, 1/100, 1/10, 1, 10, 100)
bflabs <- c("1/1000", "1/100", "1/10", "1", "10", "100")
matplot(t0seq, y = bfupper, type = "l", lty = 2,
        col = cols2, log = "y", ylim = c(1/1000, 120),
        las = 1, xlab = bquote("Mean" ~ theta),
        ylab = bquote("SC distribution (" *
                      .(round(plower*100, 1)) * "/" * .(round(p*100, 1)) * "/" *
                          .(round(pupper*100, 1)) * "% quantiles)"),
        yaxt = "n",
        panel.first = graphics::grid(lty = 3, equilogs = FALSE))
axis(side = 2, at = bfbks, labels = bflabs, las = 1)
for (i in seq(1, length(ns))) {
    polygon(x = c(t0seq, rev(t0seq)), y = c(bfupper[,i], rev(bflower[,i])),
            col = cols3[i], border = FALSE)
}
matlines(t0seq, y = bflower, type = "l", lty = 2,
         col = cols2)
matlines(t0seq, y = bfmedians, type = "l", lwd = 1.5, lty = 1,
         col = cols)
abline(h = 1, lty = 3, col = adjustcolor(col = "black", alpha.f = 0.1))
abline(v = tstar, lty = 3, col = adjustcolor(col = "black", alpha.f = 0.1))
legend("topright", legend = rev(ns), title = expression("Sample size" ~ italic(n)),
       lty = 1, lwd = 1.5, col = rev(cols), bg = "white")
arrows(x0 = min(t0seq), y0 = c(1.5, 1/1.5), y1 = c(5, 1/5),
       col = adjustcolor("black", alpha.f = 0.8), length = 0.05)
text(x = min(t0seq), y = c(2.5, 1/2.5),
     labels = c(expression("Support for" ~ theta),
                expression("Support for" ~ italic(H)[1])),
     pos = 4, cex = 0.75, col = adjustcolor("black", alpha.f = 0.8))


## ----"Bartos-data"------------------------------------------------------------
## data from Bartos
y <- 178079
n <- 350757
a <- 5100
b <- 4900
u <- 1
l <- 0.5


## ----"Bartos-analysis", fig.height = 4----------------------------------------
## BFF for binomial proportion based on truncated beta prior under the alternative
BFFbinomial <- function(p, y, n, a, b, l = 0.5, u = 1, log = FALSE) {
    logbf <- y*log(p) + (n - y)*log(1 - p) - lbeta(a + y, b + n - y) + lbeta(a, b) +
        log(pbeta(q = u, shape1 = a, shape2 = b) -
            pbeta(q = l, shape1 = a, shape2 = b)) -
        log(pbeta(q = u, shape1 = a + y, shape2 = b + n - y) -
            pbeta(q = l, shape1 = a + y, shape2 = b + n - y))
    if (log == TRUE) return(logbf)
    else return(exp(logbf))
}

## support interval for binomial proportion
SIbinomial <- function(k, y, n, a, b, l = 0.5, u = 1) {
    mee <- y/n
    rootFun <- function(p) {
        BFFbinomial(p = p, y = y, n = n, a = a, b = b, l = l, u = u) - k
    }
    lower <- try({uniroot(f = rootFun, interval = c(0, mee))$root})
    upper <- try({uniroot(f = rootFun, interval = c(mee, 1))$root})
    if (inherits(lower, "try-error")) res <- c(NaN, NaN, NaN)
    else res <- c(lower, mee, upper)
    names(res) <- c("lower", "mee", "uppper")
    return(res)
}

## compute MEE
mee <- y/n
kME <- BFFbinomial(p = mee, y = y, n = n, a = a, b = b)

## compute k = 1 support interval
si <- SIbinomial(k = 1, y = y, n = n, a = a, b = b)

## compute BF for p = 0.5
bf05 <- BFFbinomial(p = 0.5, y = y, n = n, a = a, b = b)

## compute BFF
p0seq <- seq(from = 0.5, to = 0.515, length.out = 500)
bff <- BFFbinomial(p = p0seq, y = y, n = n, a = a, b = b)

## plot results
bks <- c(10^seq(-18, -3, 3), 1)
labs <- c(expression(10^-18), expression(10^-15), expression(10^-12),
          expression(10^-9), expression(10^-6), expression(10^-3), "1")
transpblack <- adjustcolor(col = "black", alpha.f = 0.5)
par(mar = c(4, 5, 2, 1.5))
plot(x = p0seq, y = bff,
     type = "l", log = "y", lwd = 1.5,
     ylab = bquote("Bayes factor"),
     xlab = bquote("Probability of coin landing on same side" ~ theta),
     ylim = c(1/10^18, 10),
     yaxt = "n",
     panel.first = graphics::grid(lty = 3, equilogs = FALSE))
axis(side = 2, at = bks, labels = labs, las = 1, cex = 0.6)
abline(h = 1, lty = 2, col = "lightgrey")
arrows(x0 = si[1], x1 = si[3], y0 = kME, code = 3, angle = 90, length = 0,
       col = 1)
points(x = mee, y = kME, pch = 20, col = 1)
arrows(x0 = min(p0seq), y0 = c(2, 1/2), y1 = c(40, 1/40),
       col = adjustcolor("black", alpha.f = 0.8), length = 0.05)
text(x = min(p0seq), y = c(7, 1/7),
     labels = c(expression("Support for" ~ theta),
                expression("Support for" ~ italic(H)[1])),
     pos = 4, cex = 0.75, col = adjustcolor("black", alpha.f = 0.8))


## ----"bartos-meta-analysis1", cache = TRUE------------------------------------
## flipper-wise data from Bartos (2023)
dat <- dat.bartos2023

## compute proportions and their variances
dat$yi <- with(dat, same/flips)
dat$vi <- with(dat, yi*(1 - yi)/flips)
dat$lower <- dat$yi - qnorm(p = 0.975)*sqrt(dat$vi)
dat$upper <- dat$yi + qnorm(p = 0.975)*sqrt(dat$vi)

## perform random effects meta-analysis
prior <- function(x) {
    dbeta(x = x, shape1 = a, shape2 = b) /
        (pbeta(q = u, shape1 = a, shape2 = b) - pbeta(q = l, shape1 = a, shape2 = b)) *
        as.numeric(l <= x & x <= u)
}
scale <- 0.02
res <- metabf(yi = dat$yi, sei = sqrt(dat$vi), labels = dat$person,
              theta1 = prior, tau1 = function(x) dnorm(x, sd = scale)*2,
              control = list(thetaLim = c(l, u), thetaSEmultSearch = 5,
                             tauUpperSearch = 5, subdivisions = 100L,
                             rel.tol = .Machine$double.eps^0.5,
                             abs.tol = .Machine$double.eps^0.5,
                             tol = .Machine$double.eps^0.5, maxiter = 1000))
plotdat <- plot(res, thetaRange = c(0.5, 0.52), tauRange =  c(0, 0.04),
                plot = FALSE, ngrid = 200)

## perform same analysis for different scale parameters for the heterogeneity prior
scales <- c(0.005, 0.01, scale, 0.03, 0.04)
sensitivityList <- lapply(X = scales, FUN = function(scale) {
    ma <- metabf(yi = dat$yi, sei = sqrt(dat$vi), labels = dat$person,
                 theta1 = prior, tau1 = function(x) dnorm(x, sd = scale)*2,
                 control = list(thetaLim = c(l, u), thetaSEmultSearch = 5,
                                tauUpperSearch = 5, subdivisions = 100L,
                                rel.tol = .Machine$double.eps^0.5,
                                abs.tol = .Machine$double.eps^0.5,
                                tol = .Machine$double.eps^0.5, maxiter = 1000))
    plotdat <- plot(ma, thetaRange = c(0.5, 0.52), tauRange =  c(0, 0.04),
                    plot = FALSE)
    list("tauDF" = data.frame(plotdat$tauDF, scale = scale),
         "meetau" = data.frame(unname(as.data.frame(t(plotdat$MEEtau))),
                               scale = scale, k = ma$ktau),
         "thetaDF" = data.frame(plotdat$thetaDF, scale = scale),
         "meetheta" = data.frame(unname(as.data.frame(t(plotdat$MEEtheta))),
                                 scale = scale, k = ma$ktheta)
         )
})
tausens <- do.call("cbind", (lapply(X = sensitivityList,
                                    FUN = function(x) x$tauDF$bf)))
taumee <- do.call("rbind", (lapply(X = sensitivityList,
                                   FUN = function(x) x$meetau)))
thetasens <- do.call("cbind", (lapply(X = sensitivityList,
                                      FUN = function(x) x$thetaDF$bf)))
thetamee <- do.call("rbind", (lapply(X = sensitivityList,
                                     FUN = function(x) x$meetheta)))


## ----"bartos-meta-analysis2", fig.height = 10, fig.width = 8------------------
par(mfrow = c(2, 2), mar = c(4.9, 5, 4.1, 1.5))

## forest plot
ybreaks <- seq(length(dat$yi), 1, -1)
plot(x = dat$yi, y = ybreaks, type = "n",
     xlim = c(min(dat$lower), max(dat$upper)),
     yaxt = "n", ylim = c(0.5, length(dat$yi) + 0.5),
     xlab = "Probability estimate with 95% CI",
     ylab = "",
     panel.first = graphics::grid(ny = NA, lty = 3),
     main = bquote("Data from Bartos et al. (2023)" * ""))
axis(side = 2, at = ybreaks, labels = dat$person, las = 2, cex.axis = 0.5)
mtext(text = "Participant", side = 2, line = 4)
arrows(x = dat$lower, x1 = dat$upper, y0 = ybreaks, angle = 90, code = 3,
       length = 0)
points(x = dat$yi, y = ybreaks, pch = 20, cex = 1)

## 2D BFF plot
bfMatrix <- matrix(plotdat$contourDF$bf, ncol = 200, byrow = TRUE)
image(x = plotdat$tauDF$tau, y = plotdat$thetaDF$theta, z = bfMatrix,
      col = grDevices::hcl.colors(n = 100, palette = "Blues 3", rev = TRUE),
      las = 1,
      xlab = bquote("Heterogeneity" ~ tau),
      ylab = bquote("Probability of coin landing on same side" ~ theta),
      main = bquote("Bayes factor surface" ~ ""))
contour(x = plotdat$tauDF$tau, y = plotdat$thetaDF$theta, z = bfMatrix,
        levels = c(1/1000, 1/100, 1/10, 3, 10, 30, 100), add = TRUE,
        col = "#00000080", lty = 3)
contour(x = plotdat$tauDF$tau, y = plotdat$thetaDF$theta, z = bfMatrix,
        levels = 1, add = TRUE, col = "#00000080", lty = 2, lwd = 1.5)
contour(x = plotdat$tauDF$tau, y = plotdat$thetaDF$theta, z = bfMatrix,
        levels = res$k, add = TRUE, col = transpblack, lty = 2, lwd = 1.5)
points(x = res$MEEjoint[2], y = res$MEEjoint[1], pch = 20, col = 1)
legend("topright", legend = "", title = paste("HN prior scale", scale),
       bty = "n", cex = 0.85)

## BFF plot for probability
colors <- hcl.colors(n = length(scales), alpha = 0.3)
colors[4:5] <- colors[3:4]
colors[3] <- "black"
lty <- 1
bfbks <- c(1/1000, 1/100, 1/10, 1, 10, 100)
bflabs <- c("1/1000", "1/100", "1/10", "1", "10", "100")
matplot(sensitivityList[[1]]$thetaDF$theta, thetasens, type = "l", lty = lty,
        ylim = c(1/1000, 100), lwd = 1.5, col = colors, las = 1, log = "y",
        main = bquote(tau ~ "is the nuisance parameter"),
        xlab = bquote("Probability of coin landing on same side" ~ theta),
        ylab = "Bayes factor",
        panel.first = graphics::grid(lty = 3, ny = NA, equilogs = FALSE),
        yaxt = "n")
axis(side = 2, at = bfbks, labels = bflabs, las = 1)
abline(h = bfbks, lty = 3, col = adjustcolor(col = 1, alpha = 0.1))
arrows(x0 = thetamee$X1, x1 = thetamee$X3, y0 = thetamee$k, col = colors, code = 3,
       angle = 90, length = 0, lty = lty)
points(x = thetamee$X2, y = thetamee$k, col = colors, pch = 20)
abline(h = 1, lty = 2, col = transpblack)
legend("topright", title = "HN prior scale", bg = "white", pch = 20,
       legend = rev(scales), col = rev(colors), lty = rev(lty), lwd = 1.5, cex = 0.85)
arrows(x0 = 0.5, y0 = c(1.5, 1/1.5), y1 = c(3, 1/3),
       col = adjustcolor("black", alpha.f = 0.8), length = 0.05)
text(x = 0.5, y = c(2, 1/2),
     labels = c(expression("Support for" ~ theta),
                expression("Support for" ~ italic(H)[1])),
     pos = 4, cex = 0.6, col = adjustcolor("black", alpha.f = 0.8))

## BFF plot for heterogeneity
matplot(sensitivityList[[1]]$tauDF$tau, tausens, type = "l", lty = lty,
        ylim = c(1/1000, 100), lwd = 1.5, col = colors, las = 1, log = "y",
        main = bquote(theta ~ "is the nuisance parameter"),
        xlab = bquote("Heterogeneity" ~ tau), ylab = "Bayes factor",
        panel.first = graphics::grid(lty = 3, ny = NA, equilogs = FALSE),
        yaxt = "n")
axis(side = 2, at = bfbks, labels = bflabs, las = 1)
abline(h = bfbks, lty = 3, col = adjustcolor(col = 1, alpha = 0.1))
arrows(x0 = taumee$X1, x1 = taumee$X3, y0 = taumee$k, col = colors, code = 3,
       angle = 90, length = 0, lty = lty)
points(x = taumee$X2, y = taumee$k, col = colors, pch = 20)
abline(h = 1, lty = 2, col = transpblack)
arrows(x0 = 0, y0 = c(1.5, 1/1.5), y1 = c(3, 1/3),
       col = adjustcolor("black", alpha.f = 0.8), length = 0.05)
text(x = 0, y = c(2, 1/2),
     labels = c(expression("Support for" ~ tau),
                expression("Support for" ~ italic(H)[1])),
     pos = 4, cex = 0.6, col = adjustcolor("black", alpha.f = 0.8))


## ----"prior-for-tau"----------------------------------------------------------
tauseq <- seq(0, 0.1, 0.001)
## CDF of the halfnormal distribution
phn. <- function(x, s) integrate(f = function(x) 2*dnorm(x, sd = s),
                                 lower = 0, upper = x)$value
phn <- Vectorize(FUN = phn.)
## plot(tauseq, phn(x = tauseq, s = scale), type = "l")
## phn(x = 0.05, s = scale)


## ----"data-protzko", eval = FALSE---------------------------------------------
## ## load data from Protzko et al. (2023)
## protzko <- protzko2020
## labels <- subset(protzko, experiment == "Labels")
## yo <- subset(labels, type == "original")$smd
## so <- subset(labels, type == "original")$se
## yr <- subset(labels, type == "external-replication")$smd
## sr <- subset(labels, type == "external-replication")$se
## lab <- as.numeric(as.character(subset(labels, type == "external-replication")$lab))


## ----"data-wagenmakers"-------------------------------------------------------
## Data from 10.1177/1745691616674458
lower <- c(-0.05, -0.4, -0.45, -0.52, -0.43, -0.69, -0.65, -0.35, -0.73, -0.65,
           -0.19, -0.76, -1.41, -0.76, -0.51, -0.27, -0.74, -0.18)
upper <- c(1.69, 0.67, 0.75, 0.49, 0.67, 0.46, 0.54, 0.75, 0.35, 0.69, 0.92,
           0.28, 0.26, 0.49, 0.56, 0.57, 0.34, 0.88)
mdRecalc <- rowMeans(x = cbind(lower, upper))
se <- apply(X = cbind(lower, upper), MARGIN = 1,
            FUN = function(x) diff(x)/(2*qnorm(p = 0.975)))
md <- c(0.82, 0.14, 0.16, -0.02, 0.12, -0.11, -0.05, 0.2, -0.19, 0.02, 0.37,
        -0.24, -0.58, -0.13, 0.02, 0.15, -0.2, 0.35)
site <- c("Original", "Albohn", "Allard", "Benning", "Bulnes", "Capaldi",
          "Chasten", "Holmes", "Koch", "Korb", "Lynott", "Oosterwijk",
          "Ozdogru", "Pacheco-Unguetti", "Talarico", "Wagenmakers", "Wayand",
          "Zeelenberg")
wagenmakers2016 <- data.frame(site, md, se)
yo <- md[1]
so <- se[1]
yr <- md[-1]
sr <- se[-1]
srpool <- sqrt(1/sum(1/sr^2))
yrpool <- sum(yr/sr^2)*srpool^2
cipool <- yrpool + c(-1, 1)*srpool*qnorm(p = 0.975)
smallest <- which.min(yr)


## ----"replication-analysis", fig.height = 6.5---------------------------------
## replication BF
repBFF <- function(null, yr, sr, yo, so, log = FALSE) {
    logbf <- dnorm(x = yr, mean = null, sd = sr, log = TRUE) -
        dnorm(x = yr, mean = yo, sd = sqrt(sr^2 + so^2), log = TRUE)
    if (log == TRUE) return(logbf)
    else return(exp(logbf))
}

## support interval based on replication BF
SIrep <- function(k, yr, sr, yo, so, log = FALSE) {
    si <- yr + c(-1, 1)*sr*
        sqrt(log(1 + so^2/sr^2) + (yr - yo)^2/(sr^2 + so^2) - 2*log(k))
    res <- c(si[1], yr, si[2])
    names(res) <- c("lower", "mee", "upper")
    return(res)
}


## compute BFF, k=1 support interval, and MEE
mdseq <- c(seq(-2, 0, length.out = 500), seq(0.001, 2, length.out = 500))
bff <- sapply(X = seq(1, length(yr)), FUN = function(i) {
    repBFF(null = mdseq, yr = yr[i], sr = sr[i], yo = yo, so = so)
})
bffpool <- repBFF(null = mdseq, yr = yrpool, sr = srpool, yo = yo, so = so)
si <- t(sapply(X = seq(1, length(yr)), FUN = function(i) {
    SIrep(k = 1, yr = yr[i], sr = sr[i], yo = yo, so = so)
}))
sipool <- SIrep(k = 1, yr = yrpool, sr = srpool, yo = yo, so = so)
kme <- sapply(X = seq(1, length(yr)), FUN = function(i) {
    repBFF(null = yr[i], yr = yr[i], sr = sr[i], yo = yo, so = so)
})
kmepool <- repBFF(null = yrpool, yr = yrpool, sr = srpool, yo = yo, so = so)

## plot inferences
par(mar = c(4, 4, 1, 1), mfrow = c(2, 1))
cols <- adjustcolor(col = 1, alpha.f = 0.2)
bks <- 10^seq(-6, 3, 1)
labs <- c(expression(10^-6), expression(10^-5), expression(10^-4),
          expression(10^-3), expression(10^-2), expression(10^-1), "1",
          expression(10^1), expression(10^2), expression(10^3))
matplot(mdseq, bff, type = "l", lty = 1, ylim = c(1/10^5, 10^3),
        lwd = 1.5, col = cols, las = 1, log = "y",
        xlab = bquote("Mean difference" ~ theta),
        ylab = "Bayes factor",
        panel.first = graphics::grid(lty = 3, equilogs = FALSE),
        yaxt = "n")
axis(side = 2, at = bks, labels = labs, las = 1)
abline(h = 1, lty = 2, col = adjustcolor(col = 1, alpha = 0.3))
## arrows(x0 = si[,1], x1 = si[,3], y0 = kme, col = cols, code = 3,
##        angle = 90, length = 0, lty = 1)
points(x = si[,2], y = kme, col = cols, pch = 20)
lines(mdseq, bffpool, col = 2, lwd = 1.5)
arrows(x0 = sipool[1], x1 = sipool[3], y0 = kmepool, col = 2, code = 3,
       angle = 90, length = 0, lty = 1)
points(x = sipool[2], y = kmepool, col = 2, pch = 20)
arrows(x0 = min(mdseq), y0 = c(1.5, 1/1.5), y1 = c(5, 1/5),
       col = adjustcolor("black", alpha.f = 0.8), length = 0.05)
text(x = min(mdseq), y = c(3, 1/3),
     labels = c(expression("Support for" ~ theta),
                expression("Support for" ~ italic(H)[1])),
     pos = 4, cex = 0.7, col = adjustcolor("black", alpha.f = 0.8))
legend("topright", lty = c(1, 1), col = c(cols, 2),
       legend = c("Support curves based on individual replications",
                  "Support curve based on pooled replications"),
       ## legend = expression(theta ~ "|" ~ italic(H)[1] ~ "~ N(" * italic(y)["o"] * "," ~ sigma["o"]^2 *")"),
       bg = "white", cex = 0.6, lwd = c(1, 1.5))

## plot posterior
prior <- function(md) dnorm(x = md, mean = yo, sd = so)
posterior <- sapply(X = seq(1, length(yr)), FUN = function(i) {
    repBFF(null = mdseq, yr = yr[i], sr = sr[i], yo = yo, so = so)*prior(mdseq)
})
posteriorpool <- repBFF(null = mdseq, yr = yrpool, sr = srpool, yo = yo, so = so)*prior(mdseq)
cri <- t(sapply(X = seq(1, length(yr)), FUN = function(i) {
    vpost <- 1/(1/so^2 + 1/sr[i]^2)
    mpost <- (yr[i]/sr[i]^2 + yo/so^2)*vpost
    cri <- mpost + c(-1, 1)*sqrt(vpost)*qnorm(p = 0.975)
    height <- dnorm(mpost, mpost, sqrt(vpost))
    res <- c(cri[1], mpost, cri[2], height)
    names(res) <- c("lower", "mee", "upper", "height")
    return(res)
}))
vpostpool <- 1/(1/so^2 + 1/srpool^2)
mpostpool <- (yrpool/srpool^2 + yo/so^2)*vpostpool
cripool <- mpostpool + c(-1, 1)*sqrt(vpostpool)*qnorm(p = 0.975)
height <- dnorm(mpostpool, mpostpool, sqrt(vpostpool))
res <- c(cripool[1], mpostpool, cripool[2], height)
names(res) <- c("lower", "mee", "upper", "height")
matplot(mdseq, posterior, type = "l", lty = 1,
        lwd = 1.5, col = cols, las = 1,
        xlab = bquote("Mean difference" ~ theta),
        ylab = "Posterior density",
        panel.first = graphics::grid(lty = 3, equilogs = FALSE), ylim = c(0, height))
lines(mdseq, prior(mdseq), lty = 2, col = adjustcolor("black", alpha.f = 0.8),
      lwd = 1.5)
lines(mdseq, posteriorpool, lwd = 1.5, col = 2)
## arrows(x0 = cri[,1], x1 = cri[,3], y0 = cri[,4], col = cols, code = 3,
##        angle = 90, length = 0, lty = 1)
points(x = cri[,2], y = cri[,4], col = cols, pch = 20)
arrows(x0 = res[1], x1 = res[3], y0 = res[4], col = 2, code = 3, lwd = 1.5,
       angle = 90, length = 0, lty = 1)
points(x = res[2], y = res[4], col = 2, pch = 20)
legend("topright", lty = c(2, 1, 1), col = c(1, cols, 2),
       legend = c("Prior based on original study",
                  "Posteriors based on individual replications",
                  "Posterior based on pooled replications"),
       ## legend = expression(theta ~ "|" ~ italic(H)[1] ~ "~ N(" * italic(y)["o"] * "," ~ sigma["o"]^2 *")"),
       bg = "white", cex = 0.6, lwd = c(1.5, 1, 1.5))


## ----"logistic-regression", fig.height = 6------------------------------------
## load data from Ejbye-Ernst et al. (2023)
## <https://doi.org/10.1007/s11292-023-09602-9>
## <https://osf.io/nb56d>
load(file = "../data/data.RData")
dat <- data_frame
class(dat$dow) <- "factor"
dat$dealers_presence <- haven::as_factor(dat$dealers_presence)
dat$dow <- factor(dat$dow, levels = c("wo", "do", "vr", "za"),
                  labels = c("Wednesday", "Thursday", "Friday", "Saturday"))
dat$intervention <- factor(dat$intervention_factor,
                           levels = c("Before intervention", "After intervention"),
                           labels = c("Before", "After"))
dat$camera <- factor(dat$camera)

## fit logistic regression with ML
glm1 <- glm(dealers_presence ~ crowding + dow + intervention + camera +
                hours_after_eight, family = binomial(link = "logit"), data = dat)
## summary(glm1)

## create prior for brms
pm <- 0 # prior mean
psd <- sqrt(1/2) # prior standard deviation
prior <- paste0("normal(", pm, ", ", psd, ")")

## Bayesian logistic regression model
set.seed(4242)
nmcmc <- 1000000
nchain <- 10
if ("stanmcmc.RData" %in% list.files(path = "../data/")) {
    load(file = "../data/stanmcmc.RData")
} else {
    bglm1 <- brm(dealers_presence ~ crowding + dow + intervention + camera +
                     hours_after_eight, data = dat,
                 prior = c(set_prior("", class = "Intercept"), # flat prior
                           set_prior(prior, class = "b")), iter = nmcmc/nchain,
                 chains = nchain, warmup = 500, cores = 12, seed = 4242,
                 family = "bernoulli")
    save(bglm1, file = "../data/stanmcmc.RData")
}
## ## MCMC diagnostics
## summary(bglm1)
## exp(fixef(bglm1))[-1,]
## plot(bglm1, N = 8)

## priors for INLA
pm <- 0 # prior mean
prec <- 2 # prior precision
coefs <- names(coef(glm1))[-1]
pmeans <- as.list(rep(pm, length(coefs)))
names(pmeans) <- coefs
pprec <- as.list(rep(prec, length(coefs)))
names(pprec) <- coefs

## fit model with INLA
dat$dealers_presence_bin <- as.numeric(dat$dealers_presence == "yes")
model <- dealers_presence_bin ~ crowding + dow + intervention + camera +
    hours_after_eight
bglm2 <- inla(formula = model, data = dat, family = "binomial",
              control.fixed = list(mean = pmeans, prec = pprec))
## summary(bglm2)

coefs <- c("Crowding" = "crowding",
           "Day: Thursday" = "dowThursday",
           "Day: Friday" = "dowFriday",
           "Day: Saturday" = "dowSaturday",
           "After intervention" = "interventionAfter",
           "Camera 2" = "camera1",
           "Hours after 8PM" = "hours_after_eight")
logORseq <- seq(log(1/30), log(100), length.out = 2^10)
ORseq <- exp(logORseq)
bfbks <- c(1/1000, 1/100, 1/10, 1, 10, 100, 1000)
bflabs <- c("1/1000", "1/100", "1/10", "1", "10", "100", "1000")
orbks <- c(bfbks, 1/30, 1/3, 3, 30)
orlabs <- c(bflabs, "1/30", "1/3", "3", "30")
glmdat <- summary(glm1)$coefficients

bffres <- lapply(X = seq(1, length(coefs)), FUN = function(i) {
    coef <- coefs[i]
    ## Univariate normal analysis
    glmest <- glmdat[rownames(glmdat) == coef, 1]
    glmse <- glmdat[rownames(glmdat) == coef, 2]
    ## MCMC analysis
    mcmcdraws <- brms::as_draws_matrix(bglm1, variable = paste0("b_", coef))
    priordraws <- rnorm(n = length(mcmcdraws), mean = pm, sd = psd)
    postdens <- density(mcmcdraws,  from = min(logORseq), to = max(logORseq),
                        n = length(logORseq))

    ## compute SCs
    priordens <- dnorm(x = logORseq, mean = pm, sd = psd)
    bffmcmc <- postdens$y/priordens
    bffglm <- dnorm(x = glmest, mean = logORseq, sd = glmse) /
        dnorm(x = glmest, mean = pm, sd = sqrt(psd^2 + glmse^2))
    BFFinla <- function(logOR) {
        inla.dmarginal(x = logOR, marginal = bglm2$marginals.fixed[[coef]])/
            dnorm(x = logOR, mean = pm, sd = psd)
    }
    bffinla <- BFFinla(logORseq)
    bffDF <- rbind(data.frame(logOR = logORseq, bff = bffmcmc, method = "MCMC"),
                   data.frame(logOR = logORseq, bff = bffinla, method = "INLA"),
                   data.frame(logOR = logORseq, bff = bffglm,
                              method = "Univariate normal"))
    bffDF$coef <- names(coefs)[i]

    ## compute MEEs
    meeinla <- optim(par = 0, fn = function(logOR) -log(BFFinla(logOR)),
                     method = "Brent", lower = log(1/100), upper = log(100))$par
    kinla <- BFFinla(meeinla)
    meemcmc <- logORseq[which.max(bffmcmc)]
    kmcmc <- max(bffmcmc)
    meeglm <- glmest
    kglm <- max(bffglm)

    k <- 1
    SInorm <- glmest + glmse*c(-1, 1)*
        sqrt(log(1 + psd^2/glmse^2) + (glmest - pm)^2/(psd^2 + glmse^2 - log(k)))
    lowerSInorm <- SInorm[1]
    upperSInorm <- SInorm[2]
    ## HACKY way to compute the MCMC SIs
    lowerSImcmc <- logORseq[which.min(abs(log(bffmcmc[
        logORseq < meemcmc]) - log(k)))]
    upperSImcmc <- logORseq[logORseq > meemcmc][
        which.min(abs(log(bffmcmc[logORseq > meemcmc]) - log(k)))]

    out <- list(bffDF,
                rbind(data.frame(coef = names(coefs)[i], lowerSI = lowerSImcmc,
                                 upperSI = upperSImcmc,
                                 mee = meemcmc, k = kmcmc, method = "MCMC"),
                      data.frame(coef = names(coefs)[i], lowerSI = NaN,
                                 upperSI = NaN, mee = meeinla, k = kinla,
                                 method = "INLA"),
                      data.frame(coef = names(coefs)[i], lowerSI = lowerSInorm,
                                 upperSI = upperSInorm,
                                 mee = meeglm, k = kglm, method = "Univariate normal"))
                )
    return(out)
})

## ----"logistic-regression-plot", fig.height = 6, dependson = "logistic-regression"----
meeDF <- do.call("rbind", lapply(X = bffres, FUN = function(x) x[[2]]))
bffDF <- do.call("rbind", lapply(X = bffres, FUN = function(x) x[[1]]))
cols <- palette.colors(n = 4, palette = "Okabe-Ito")[2:4]
names(cols) <- c("MCMC", "INLA", "Univariate normal")
bffDF$method <- factor(bffDF$method, levels = c("MCMC", "INLA", "Univariate normal"))
bfbks <- c(1/1000, 1/100, 1/10, 1, 10, 100, 1000)
bflabs <- c("1/1000", "1/100", "1/10", "1", "10", "100", "1000")
orbks <- c(32, 16, 8, 4, 2, 1, 1/2, 1/4, 1/8, 1/16, 1/32)
orlabs <- c("32", "16", "8", "4", "2", "1", "1/2", "1/4", "1/8", "1/16", "1/32")
ggplot(data = subset(bffDF, bff != 0),
       aes(x = exp(logOR), y = bff, color = method)) +
    facet_wrap(~ coef, ncol = 2) +
    geom_hline(yintercept = 1, col = 1, linetype = "dashed", alpha = 0.2) +
    ## geom_vline(xintercept = 1, col = 1, linetype = "dashed", alpha = 0.2) +
    geom_line(aes(linetype = method), linewidth = 0.8, alpha = 0.8) +
    geom_point(data = meeDF, aes(x = exp(mee), y = k), show.legend = FALSE, alpha = 0.8) +
    coord_cartesian(xlim = c(1/10, 30), ylim = c(1/300, 300)) +
    scale_x_log10(breaks = orbks, labels = orlabs) +
    scale_y_log10(breaks = bfbks, labels = bflabs) +
    scale_color_manual(values = cols) +
    labs(x = "Odds ratio", y = "Bayes factor",
         color = "Computational method", linetype = "Computational method") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.position = c(0.85, 0.075),
          legend.key.size = unit(1.5, "lines"))



## ----"sessionInfo2", echo = Reproducibility, results = Reproducibility--------
cat(paste(Sys.time(), Sys.timezone(), "\n"))
sessionInfo()

