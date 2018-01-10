library(tidyverse)

prob_dist <- function(surv, lcl, ucl, ci = 0.95, arm) {

  lci <- log(-log(lcl))
  uci <- log(-log(ucl))
  lsurv <- log(-log(surv))

  z <- qnorm((1 + ci)/2)
  lsd <- abs(lsurv - lci)/z
  usd <- abs(uci - lsurv)/z

  sd <- (lsd + usd)/2
  xs <- log(-log(seq(from = 0, to = 1, length.out = 1000)))

  samp <- dnorm(x = xs, mean = lsurv, sd = sd)

  tibble::tibble(arm = arm, xs = exp(-exp(xs)), samp = samp)
}


sim_event <- function(surv, lcl, ucl, ci = 0.95, N_sim = 1000) {
  # ci = log(-log(surv)) + z * (sqrt(Var))
  # lcl = exp(-exp(lci))
  # lci = log(-log(lcl))
  lci <- log(-log(lcl))
  uci <- log(-log(ucl))
  lsurv <- log(-log(surv))

  z <- qnorm((1 + ci)/2)
  lsd <- abs(lsurv - lci)/z
  usd <- abs(uci - lsurv)/z

  sd <- (lsd + usd)/2

  pop_surv <- exp(-exp(rnorm(N_sim, mean = lsurv, sd = sd)))

  events <- vector("numeric", N_sim)

  for (i in seq_len(N_sim)) {
    events[i] <- sample(c(0, 1), size = 1, replace = TRUE,
                        prob = c(pop_surv[i], 1 - pop_surv[i]))
  }

  events
}

simulation <- function(s1, lcl1, ucl1, s2, lcl2, ucl2, n_sim) {

  foo <- function() {
    ev1 <- sim_event(s1, lcl1, ucl1, 0.95, 1000)
    ev2 <- sim_event(s2, lcl2, ucl2, 0.95, 1000)

    c(prop_01 = mean(1*((ev1 == 0) & (ev2 == 1))),
      prop_10 = mean(1*((ev1 == 1) & (ev2 == 0))),
      prop_11 = mean(1*((ev1 == 1) & (ev2 == 1))),
      prop_00 = mean(1*((ev1 == 0) & (ev2 == 0))))
  }

  rerun(n_sim, foo())
}

