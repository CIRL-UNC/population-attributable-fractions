library(tidyverse)
library(multiroot)
# Data Generating Process ----
dgp = function(seed, n=1200, m=12000,
               study_exposure=0.4, study_w=0.25,
               target_exposure=0.4, target_w=0.25){
  study = data.frame(
    a = rbinom(n, 1, study_exposure),
    w = rbinom(n, 1, study_w)
  ) |>
    mutate(
      ya0 = rbinom(n, 1, 0.1),
      ya1 = rbinom(n, 1, 0.3 + 0.1 * w),
      y = a * ya1 + (1-a) * ya0,
      pop = 0
    )

  target = data.frame(
    a = rbinom(m, 1, target_exposure),
    w = rbinom(m, 1, target_w)
  ) |>
    mutate(
      ya0 = rbinom(m, 1, 0.1),
      ya1 = rbinom(m, 1, 0.3 + 0.1 * w),
      y = a * ya1 + (1-a) * ya0,
      pop = 1
    )

  bind_rows(study, target)
}

dgp(seed=1) |>
  group_by(pop, a, w) |>
  summarise(mean(ya0), mean(ya1), mean(y))

dgp(seed=1, target_exposure = 0.2, target_w = 0.5) |>
  group_by(pop, a, w) |>
  summarise(
    n=n(),
    y0 = mean(ya0),
    y1 = mean(ya1),
    y = mean(y)
  )


# M-Estimation for Population Attributable Fraction

paf_eefun = function(data, eval=T){
  # Return curried function
  function(theta){
    pr_ahat_l = with(data, (pop == 0) * plogis(theta[2] + theta[3] * w))
    weight_a0 = (1-theta[1]) / (1 - pr_ahat_l)

    pr_shat_aw = with(data, plogis(theta[4] + theta[5] * a + theta[6] * w + theta[7] * a * w))
    weight_iosw_aw = with(data, (pop == 0) * (pr_shat_aw) / (1-pr_shat_aw))

    pr_shat_w = with(data, plogis(theta[8] + theta[9] * w))
    weight_iosw_w = with(data, (pop == 0) * (pr_shat_w) / (1-pr_shat_w))


    ee_eval = data |>
      transmute(
        theta1 = (pop == 0) * (a - theta[1]),
        theta2 = (pop == 0) * (a - pr_ahat_l),
        theta3 = (pop == 0) * (a - pr_ahat_l) * w,

        theta4 = (pop - pr_shat_aw),
        theta5 = (pop - pr_shat_aw) * a,
        theta6 = (pop - pr_shat_aw) * w,
        theta7 = (pop - pr_shat_aw) * a * w,

        theta8 = (pop - pr_shat_w),
        theta9 = (pop - pr_shat_w) * w,

        # For y
        theta10 = (pop == 0) * (y * weight_iosw_aw - theta[10] * weight_iosw_aw),

        # For y0
        theta11 = (pop == 0) * (y * (1-a) * weight_a0 * weight_iosw_w - theta[11] * (1-a) * weight_a0 * weight_iosw_w),

        # For PAF
        theta12 = theta[10] - theta[11] - theta[10] * theta[12],

        # For incorrect y0
        theta13 = (pop == 0) * (y * (1-a) * weight_a0 - theta[13] * (1-a) * weight_a0),

        # For incorrect y
        theta14 = (pop == 0) * (y - theta[14]),
        # For incorrect paf
        theta15 = theta[14] - theta[13] - theta[14] * theta[15]
      )

    if(eval){
      return(ee_eval |> colSums())
    } else {
      return(ee_eval)
    }
  }
}

library(rootSolve)
library(numDeriv)

true_paf = function(data, mypop=1){
  data |>
    filter(pop == mypop) |>
    summarise(
      pr_y = mean(y),
      pr_y0 = mean(ya0),
      paf = (pr_y - pr_y0) / pr_y
    )
}

# Simulation 1: Random Sample ----
single_pass_sim1 = function(iteration){
  d = dgp(seed=iteration)

  theta1 = d |> select(a) |> colMeans()
  theta23 = glm(pop ~ w, data=d |> filter(pop == 0)) |> coef()
  theta4567 = glm(pop ~ a * w, data = d) |> coef()
  theta89 = glm(pop ~ w, data=d) |> coef()

  roots = multiroot(paf_eefun(d), start=c(theta1, theta23, theta4567, theta89, 0.2, 0.1, 0.5, 0.2, 0.1, 0.5))

  ddot = jacobian(paf_eefun(d), roots$root)
  n = nrow(d)
  ainv = solve(-ddot / n)
  b = paf_eefun(d, eval=F)(roots$root)
  bread = t(b) %*% as.matrix(b) / n
  sandwich = ainv %*% bread %*% t(ainv) / n

  paf = true_paf(d)
  se_right = sandwich |> diag() |> nth(n=12) |> sqrt()
  bounds_right = roots$root[12] + c(-1.96, 1.96) * se_right
  se_wrong = sandwich |> diag() |> nth(n=15) |> sqrt()
  bounds_wrong = roots$root[15] + c(-1.96, 1.96) * se_wrong
  contained_right = between(paf$paf, bounds_right[1], bounds_right[2])
  contained_wrong = between(paf$paf, bounds_wrong[1], bounds_wrong[2])
  return(c(paf,
           calc_paf=roots$root[12],
           bounds_right=bounds_right,
           contained_right=contained_right,
           calc_wrong=roots$root[15],
           bounds_wrong=bounds_wrong,
           contained_wrong=contained_wrong))
}


single_pass_sim1(10)
many_paf_passes = as.data.frame(do.call(rbind, lapply(1:5000, single_pass_sim1)))


many_paf_passes |>
  unnest(cols=everything()) |>
  mutate(bias = calc_paf - overall_true_paf,
         bias_wrong = calc_wrong - overall_true_paf) |>
  summarise(
    avg_paf = mean(paf),
    avg_estimate = mean(calc_paf),
    avg_estimate_wrong = mean(calc_wrong),
    avg_bias = mean(bias),
            avg_bias_wrong = mean(bias_wrong),
            ese = sd(bias),
            ese_wrong = sd(bias_wrong),
            coverage_right = mean(contained_right),
            coverage_wrong = mean(contained_wrong)) |> view()


# Simulation 2: Not Random ----
single_pass_sim2 = function(iteration){
  d = dgp(seed=iteration, study_exposure = 0.2, study_w = 0.5)

  theta1 = d |> select(a) |> colMeans()
  theta23 = glm(pop ~ w, data=d |> filter(pop == 0)) |> coef()
  theta4567 = glm(pop ~ a * w, data = d) |> coef()
  theta89 = glm(pop ~ w, data=d) |> coef()

  roots = multiroot(paf_eefun(d), start=c(theta1, theta23, theta4567, theta89, 0.2, 0.1, 0.5, 0.2, 0.1, 0.5))

  ddot = jacobian(paf_eefun(d), roots$root)
  n = nrow(d)
  ainv = solve(-ddot / n)
  b = paf_eefun(d, eval=F)(roots$root)
  bread = t(b) %*% as.matrix(b) / n
  sandwich = ainv %*% bread %*% t(ainv) / n

  paf = true_paf(d)
  se_right = sandwich |> diag() |> nth(n=12) |> sqrt()
  bounds_right = roots$root[12] + c(-1.96, 1.96) * se_right
  se_wrong = sandwich |> diag() |> nth(n=15) |> sqrt()
  bounds_wrong = roots$root[15] + c(-1.96, 1.96) * se_wrong
  contained_right = between(paf$paf, bounds_right[1], bounds_right[2])
  contained_wrong = between(paf$paf, bounds_wrong[1], bounds_wrong[2])
  return(c(paf,
           calc_paf=roots$root[12],
           bounds_right=bounds_right,
           contained_right=contained_right,
           calc_wrong=roots$root[15],
           bounds_wrong=bounds_wrong,
           contained_wrong=contained_wrong))
}

single_pass_sim2(1)
many_paf_passes_sim2 = as.data.frame(do.call(rbind, lapply(1:5000, single_pass_sim2)))
#ESE, ASE, Bias
# 2 simulations: where W is equal across, and where W differs across

many_paf_passes_sim2 |>
  unnest(cols=everything()) |>
  mutate(bias = calc_paf - overall_true_paf,
         bias_wrong = calc_wrong - overall_true_paf,
         covered_right = (overall_true_paf > bounds_right1) * (overall_true_paf < bounds_right2),
         covered_wrong = (overall_true_paf > bounds_wrong1) * (overall_true_paf < bounds_wrong2)) |>
  summarise(
    avg_paf = mean(paf),
    avg_estimate = mean(calc_paf),
    avg_estimate_wrong = mean(calc_wrong),
    avg_bias = mean(bias),
            avg_bias_wrong = mean(bias_wrong),
            ese = sd(calc_paf),
            ese_wrong = sd(calc_wrong),
            coverage = mean(covered_right),
            coverage_wrong = mean(covered_wrong)) |>
  view()

