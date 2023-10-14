# 0. Libraries ----
library(tidyverse)
library(rms)
library(geex)
library(broom)
library(rootSolve)
library(ggthemes)
library(numDeriv)
library(MASS)
library(table1)

# 1. Data Cleaning - WIHS dataset ----

wihs = read_csv("data/lau.csv") |>
  mutate(
    # Categorize age to match CDC estimates
    age = case_when(13 <= ageatfda & ageatfda < 25 ~ "13-24",
                    25 <= ageatfda & ageatfda < 35 ~ "25-34",
                    35 <= ageatfda & ageatfda < 45 ~ "35-44",
                    45 <= ageatfda & ageatfda < 55 ~ "45-54",
                    55 <= ageatfda  ~ "55+"),
    # Right censor events after 10 years
    d10 = ifelse(t > 10, 0, dthev)
  ) |>
  mutate(
    ageknots = rcs(ageatfda) %>% as.data.frame() %>% as.matrix(),
    cd4knots = rcs(cd4nadir) %>% as.data.frame() %>% as.matrix(),
    Cases = 1
  ) |>
  rename(aa_indicator=black,
         idu=BASEIDU)


set.seed(17)
wihs_subset = bind_rows(
  wihs |> filter(aa_indicator == 1),
  wihs |> filter(aa_indicator == 0) |> slice_sample(prop = 0.5)
)



# 2. Data Cleaning - CDC dataset ----

cdc_2021 = read_csv("data/hiv_diagnoses_2021_2023.csv") |>
  filter(Year == 2021) |>
  mutate(
    idu = ifelse(`Transmission Category` == "Injection drug use", 1, 0),
    aa_indicator = ifelse(`Race/Ethnicity` == "Black/African American", 1, 0)
  ) |>
  rename(
    race_ethnicity = `Race/Ethnicity`,
    age = `Age Group`) |>
  select(-c(Indicator, Year)) |>
  group_by(age, idu, aa_indicator) |>
  summarise(
    cases = sum(Cases)
  )

cdc_sample_n = function(n=500, seed=17, dataset=cdc_2021) {
  set.seed(seed)

  cdc = dataset |>
    uncount(weights=cases) |>
    ungroup() |>
    slice_sample(n=n, replace=F) |>
    mutate(
      population = 1
    ) |>
    group_by(age, idu, aa_indicator) |>
    count(name="cases")
  return(cdc)
}

cdc_sample_200 = cdc_sample_n(n=200)
cdc_sample_500 = cdc_sample_n(n=500)
cdc_sample_1000 = cdc_sample_n(n=1000)
cdc_sample_5000 = cdc_sample_n(n=5000)

# 3. Table 1 Generation ----
table1_data = bind_rows(
  target = cdc_2021 |> uncount(cases),
  study = wihs |> select(age, aa_indicator, idu),
  .id = "pop"
)

table1(~factor(age) + factor(aa_indicator) + factor(idu) | pop, data = table1_data)


# 4. Both IOSW model generation ----
generate_iosw_model = function(target_data, study_data, iosw_formula) {

  td = target_data |>
    uncount(cases) |>
    ungroup()

  sd = study_data |>
    dplyr::select(age, aa_indicator, idu, d10)

  iosw_model_data = bind_rows(
    td,
    sd,
    .id = "pop"
  ) |>
    mutate(s = ifelse(pop == "1", 1, 0))

  iosw_model = iosw_model_data |>
    glm(formula=iosw_formula,  family=binomial())


  lau_data = sd |>
    mutate(
      # Not inverted because predicting odds of being in target population
      iosw = predict(iosw_model, newdata = sd) |> exp()
    )

  return(lau_data)
}

iosw_w = generate_iosw_model(cdc_2021, wihs, s ~ age * aa_indicator)
iosw_aw = generate_iosw_model(cdc_2021, wihs, s ~ age * aa_indicator * idu)

iosw_w.subset = generate_iosw_model(cdc_2021, wihs_subset, s ~ age * aa_indicator)
iosw_aw.subset = generate_iosw_model(cdc_2021, wihs_subset, s ~ age * aa_indicator * idu)

# Subsets
{
  iosw_w.200 = generate_iosw_model(cdc_sample_200, wihs, s ~ age * aa_indicator)
  iosw_aw.200 = generate_iosw_model(cdc_sample_200, wihs, s ~ age * aa_indicator * idu)

  iosw_w.500 = generate_iosw_model(cdc_sample_500, wihs, s ~ age * aa_indicator)
  iosw_aw.500 = generate_iosw_model(cdc_sample_500, wihs, s ~ age * aa_indicator * idu)

  iosw_w.1000 = generate_iosw_model(cdc_sample_1000, wihs, s ~ age * aa_indicator)
  iosw_aw.1000 = generate_iosw_model(cdc_sample_1000, wihs, s ~ age * aa_indicator * idu)

  iosw_w.5000 = generate_iosw_model(cdc_sample_5000, wihs, s ~ age * aa_indicator)
  iosw_aw.5000 = generate_iosw_model(cdc_sample_5000, wihs, s ~ age * aa_indicator * idu)
}


# 5. IPTW model generation ----
generate_iptw_model = function(study_data, iptw_formula, stabilize_formula){
  iptw_model = glm(iptw_formula, data=study_data, family=binomial())
  iptw_stabilizer = glm(stabilize_formula, data=study_data, family=binomial())


  iptw_data = study_data |>
    mutate(
      ps = predict(iptw_model, type="response"),
      numerator = predict(iptw_stabilizer, type="response"),
      iptw = (1 - idu) * (1 - numerator) / (idu * ps + (1 - idu) * (1-ps))
    )

  return(iptw_data)

}

iptw_0 = generate_iptw_model(wihs,
                    idu ~ aa_indicator + cd4knots + ageknots,
                    idu ~ 1)

iptw_0.subset = generate_iptw_model(wihs_subset,
                                    idu ~ aa_indicator + cd4knots + ageknots,
                                    idu ~ 1)

# 6. PAF calculation ----
generate_paf_data = function(study_data, iosw_w, iosw_aw, iptw){
  data = study_data |>
    mutate(
      iosw_w = iosw_w$iosw,
      iosw_aw = iosw_aw$iosw,
      iptw = iptw$iptw
    )
}

paf_data = generate_paf_data(wihs, iosw_w, iosw_aw, iptw_0)
paf_data.subset = generate_paf_data(wihs_subset, iosw_w.subset, iosw_aw.subset, iptw_0.subset)

# Subsets
{
  paf_data.200 = generate_paf_data(wihs, iosw_w.200, iosw_aw.200, iptw_0)
  paf_data.500 = generate_paf_data(wihs, iosw_w.500, iosw_aw.500, iptw_0)
  paf_data.1000 = generate_paf_data(wihs, iosw_w.1000, iosw_aw.1000, iptw_0)
  paf_data.5000 = generate_paf_data(wihs, iosw_w.5000, iosw_aw.5000, iptw_0)
}

paf_calculation = function(data, outcome, iosw_w, iosw_aw, iptw){
  data |>
    summarise(
      y0 = weighted.mean({{outcome}}, {{iosw_w}} * {{iptw}}),
      y = weighted.mean({{outcome}}, {{iosw_aw}}),
      paf = 1 - y0 / y
    )
}

# Intermediate Results ----
paf_calculation(paf_data, d10, iosw_w, iosw_aw, iptw)

paf_data |>
  mutate(
    iosw_w = 1,
    iosw_aw = 1
  ) |>
  paf_calculation(d10, iosw_w, iosw_aw, iptw)


paf_calculation(paf_data.subset, d10, iosw_w, iosw_aw, iptw)
paf_data.subset |>
  mutate(
    iosw_w = 1,
    iosw_aw = 1
  ) |>
  paf_calculation(d10, iosw_w, iosw_aw, iptw)

# Subsets
bind_rows(
  p200=paf_calculation(paf_data.200, d10, iosw_w, iosw_aw, iptw),
  p500=paf_calculation(paf_data.500, d10, iosw_w, iosw_aw, iptw),
  p1000=paf_calculation(paf_data.1000, d10, iosw_w, iosw_aw, iptw),
  p5000=paf_calculation(paf_data.5000, d10, iosw_w, iosw_aw, iptw)
)

# 7. Nonparametric Bootstrap ----
bootstrap_step = function(iteration,
                          target_pop = cdc_2021,
                          study_pop = wihs) {

  target_pop_size = sum(target_pop$cases)

  # Get random sample of CDC, as frequency table.
  target = target_pop |>
    ungroup() |>
    slice_sample(n = target_pop_size, weight_by = cases, replace=T) |>
    count(age, aa_indicator, idu, name = "cases")

  # Get random sample of study population
  study = study_pop |>
    slice_sample(n = nrow(study_pop), replace = T)

  # Make the models for IOSW_W, IOSW_AW, IPTW
  iosw_w = generate_iosw_model(target_data = target, study_data = study, s ~ age * aa_indicator)
  iosw_aw = generate_iosw_model(target_data = target, study_data = study, s ~ age * aa_indicator * idu)
  iptw_0 = generate_iptw_model(study,
                               idu ~ aa_indicator + cd4knots + ageknots,
                               idu ~ 1)

  # Calculate the PAF
  bs_data = study |>
    mutate(
      iosw_w = iosw_w$iosw,
      iosw_w_naive = 1,
      iosw_aw = iosw_aw$iosw,
      iosw_aw_naive = 1,
      iptw = iptw_0$iptw
    )

  return(
    c(i=iteration,
      true=paf_calculation(bs_data, d10, iosw_w, iosw_aw, iptw),
      naive=paf_calculation(bs_data, d10, iosw_w_naive, iosw_aw_naive, iptw)
    ))
}

bootstrap_step(1) |> view()

bootstrap_runner = function(B=1000, bs_step, ...){
  sapply(1:B, function(x) bs_step(x, ...))
}



bootstrap = bootstrap_runner(B=1000, bootstrap_step, target_pop=cdc_2021, study_pop=wihs)
bootstrap.subset = bootstrap_runner(B=1000, bootstrap_step, target_pop=cdc_2021, study_pop=wihs_subset)

bootstrap.200 = bootstrap_runner(B=1000, bootstrap_step, target_pop=cdc_sample_200, study_pop=wihs)
bootstrap.500 = bootstrap_runner(B=1000, bootstrap_step, target_pop=cdc_sample_500, study_pop=wihs)
bootstrap.1000 = bootstrap_runner(B=1000, bootstrap_step, target_pop=cdc_sample_1000, study_pop=wihs)
bootstrap.5000 = bootstrap_runner(B=1000, bootstrap_step, target_pop=cdc_sample_5000, study_pop=wihs)


bs.se = sd(bootstrap["true.paf", ] |> unlist())
bs.se.subset = sd(bootstrap.subset["true.paf", ] |> unlist())
bs.se.naive = sd(bootstrap["naive.paf", ] |> unlist())

bs.se.200 = bootstrap.200["true.paf", ] |> unlist() |> sd(na.rm=T)
bs.se.500 = bootstrap.500["true.paf", ] |> unlist() |> sd(na.rm=T)
bs.se.1000 = bootstrap.1000["true.paf", ] |> unlist() |> sd(na.rm=T)
bs.se.5000 = bootstrap.5000["true.paf", ] |> unlist() |> sd(na.rm=T)

bootstrap_cis = function(paf, bs.se) {
  return(paf + c(-1.96, 1.96) * bs.se)
}

quantile_cis = function(bootstraps, var="true.paf"){
  return(bootstraps[var, ] |> unlist() |> quantile(c(0.025, 0.975), na.rm=T))
}

# 8. M-Estimation - root finding ----
paf_mestimation = function(theta,
                           target_pop = cdc_2021,
                           study_pop = wihs,
                           solve = T){

  target = target_pop |>
    uncount(cases)

  stacked_data = bind_rows(
    target,
    study_pop,
    .id="pop"
  ) |>
    replace_na()

  # IOSW calculation
  X.pop = stacked_data |>
    ungroup() |>
    transmute(
      population = ifelse(pop == "1", 1, 0)
    )|>
    as.matrix()

  X.iosw.w = stacked_data |>
    model.matrix(object = ~age * aa_indicator)

  X.iosw.aw = stacked_data |>
    model.matrix(object = ~age * aa_indicator * idu)

  X.iosw.w.study = X.iosw.w |>
    tail(n = nrow(study_pop))

  X.iosw.aw.study = X.iosw.aw |>
    tail(n = nrow(study_pop))

  # IPTW calculation
  # Note; because target pop has NA for values, not included in model matrix.
  X.exposure = model.matrix(~idu - 1, data=study_pop)
  X.iptw = stacked_data |>
    model.matrix(object = ~aa_indicator + cd4knots + ageknots)

  # Y0, Y, PAF calculation
  X.outcome = model.matrix(~d10 - 1, data=study_pop)

  # Estimating Equation construction
  ee_iptw = rbind(matrix(0, nrow=nrow(target), ncol=10),
                  X.iptw * (((1 - X.exposure) - plogis(X.iptw %*% theta[1:10])) %*% matrix(1, ncol=10))
                  )


  ee_iptw_numerator = rbind(matrix(0, nrow=nrow(target), ncol=1),
                            (1 - X.exposure) - theta[11])

  ee_iosw_w = X.iosw.w * ((X.pop - plogis(X.iosw.w %*% theta[12:21])) %*% matrix(1, ncol=10))

  ee_iosw_aw = X.iosw.aw * ((X.pop - plogis(X.iosw.aw %*% theta[22:41])) %*% matrix(1, ncol=20))

  weight = replace_na(
    (((1 - X.exposure) *  theta[11]) / (plogis(X.iptw %*% theta[1:10]))) * # IPTW
      exp(X.iosw.w.study %*% theta[12:21]),
    0
  )

  ee_y0 = rbind(matrix(0, nrow=nrow(target), ncol=1),
                (X.outcome - theta[42]) * weight)

  ee_y = rbind(matrix(0, nrow=nrow(target), ncol=1),
               (X.outcome - theta[43]) * exp(X.iosw.aw.study %*% theta[22:41]))


  ee_paf = matrix((theta[42] - theta[43] * (1 - theta[44])), nrow = nrow(stacked_data))

  ees = cbind(
    iptw=ee_iptw,
    iptw_num=ee_iptw_numerator,
    iosw_w=ee_iosw_w,
    iosw_aw=ee_iosw_aw,
    y0=ee_y0,
    y=ee_y,
    paf=ee_paf
  )

  if(solve) {
    return(colSums(ees))
  } else {
    return(ees)
  }
}

paf_mestimation_naive = function(theta,
                                 study_pop = wihs,
                                 solve = T){
  X.exposure = model.matrix(~idu - 1, data=study_pop)
  X.iptw = model.matrix(object = ~aa_indicator + cd4knots + ageknots, data=study_pop)

  # Y0, Y, PAF calculation
  X.outcome = model.matrix(~d10 - 1, data=study_pop)

  # Estimating Equation construction
  ee_iptw = X.iptw * (((1 - X.exposure) - plogis(X.iptw %*% theta[1:10])) %*% matrix(1, ncol=10))
  ee_iptw_numerator = (1 - X.exposure) - theta[11]

  weight = replace_na(
    (((1 - X.exposure) *  theta[11]) / (plogis(X.iptw %*% theta[1:10]))),
    0
  )

  ee_y0 = (X.outcome - theta[12]) * weight


  ee_y = (X.outcome - theta[13])


  ee_paf = matrix((theta[12] - theta[13] * (1 - theta[14])), nrow = nrow(study_pop))

  ees = cbind(
    iptw=ee_iptw,
    iptw_num=ee_iptw_numerator,
    y0=ee_y0,
    y=ee_y,
    paf=ee_paf
  )

  if(solve) {
    return(colSums(ees))
  } else {
    return(ees)
  }
}


vec_start = jitter(rep(0, times=44))
vec_start_naive = jitter(rep(0, times=14))
mr = multiroot(paf_mestimation,
               vec_start,
               target_pop = cdc_2021,
               study_pop = wihs,
               solve = T)

mr.subset = multiroot(paf_mestimation,
                      vec_start,
                      target_pop = cdc_2021,
                      study_pop = wihs_subset,
                      solve = T)

mr.naive = multiroot(paf_mestimation_naive,
                     vec_start_naive,
                     study_pop = wihs,
                     solve=T)

mr.naive2 = multiroot(paf_mestimation,
                      vec_start,
                      target_pop = wihs |> group_by(age, idu, aa_indicator) |> count(name="cases"),
                      study_pop = wihs,
                      solve = T)

mr.200 = multiroot(paf_mestimation,
               vec_start,
               target_pop = cdc_sample_200,
               study_pop = wihs,
               solve = T)

mr.500 = multiroot(paf_mestimation,
                   vec_start,
                   target_pop = cdc_sample_500,
                   study_pop = wihs,
                   solve = T)

mr.1000 = multiroot(paf_mestimation,
                   vec_start,
                   target_pop = cdc_sample_1000,
                   study_pop = wihs,
                   solve = T)

mr.5000 = multiroot(paf_mestimation,
                   vec_start,
                   target_pop = cdc_sample_5000,
                   study_pop = wihs,
                   solve = T)



# 9. M-Estimation - variance ----
paf_sandwich_variance = function(f, xval, ...){
  a0 = -1 * jacobian(func=f, x=xval, method="simple")
  evaluated = f(xval, ...)
  b0 = t(evaluated) %*% evaluated

  v0 = solve(a0) %*% b0 %*% t(solve(a0))
  return(v0)
}

paf.vcov = paf_sandwich_variance(paf_mestimation,
                      mr$root,
                      target_pop = cdc_2021,
                      study_pop = wihs,
                      solve=F)

paf.vcov.subset = paf_sandwich_variance(paf_mestimation,
                                        mr.subset$root,
                                        target_pop = cdc_2021,
                                        study_pop = wihs_subset,
                                        solve=F)

paf.vcov.200 = paf_sandwich_variance(paf_mestimation,
                                     mr.200$root,
                                     target_pop = cdc_sample_200,
                                     study_pop = wihs,
                                     solve=F)
paf.vcov.500 = paf_sandwich_variance(paf_mestimation,
                                     mr.500$root,
                                     target_pop = cdc_sample_500,
                                     study_pop = wihs,
                                     solve=F)
paf.vcov.1000 = paf_sandwich_variance(paf_mestimation,
                                     mr.1000$root,
                                     target_pop = cdc_sample_1000,
                                     study_pop = wihs,
                                     solve=F)
paf.vcov.5000 = paf_sandwich_variance(paf_mestimation,
                                     mr.5000$root,
                                     target_pop = cdc_sample_5000,
                                     study_pop = wihs,
                                     solve=F)
# Note - no difference between estimating naive population with both IOSW models.
paf.vcov.naive = paf_sandwich_variance(paf_mestimation_naive,
                                       mr.naive$root,
                                       study_pop = wihs,
                                       solve = F)
paf.vcov.naive2 = paf_sandwich_variance(paf_mestimation,
                                        mr.subset$root,
                                        target_pop = wihs |> group_by(age, idu, aa_indicator) |> count(name="cases"),
                                        study_pop = wihs,
                                        solve=F)

get_paf_variance = function(variance_matrix, rf, rootnum=44){
  paf_mest = rf$root[rootnum]
  paf_mest.variance = diag(variance_matrix)[rootnum]
  paf_mest.cis = paf_mest + c(-1.96, 1.96) * sqrt(paf_mest.variance)

  return(c(se=paf_mest.variance |> sqrt(), cis=paf_mest.cis))
}

# 10. Results ----
paf = c(
  paf = mr$root[44],
  paf.cis = get_paf_variance(paf.vcov,
                             mr),
  paf.bootstrap.cis = bootstrap_cis(mr$root[44], bs.se),
  paf.bootstrap.quantile = quantile_cis(bootstrap)
)

paf.subset = c(
  paf = mr.subset$root[44],
  paf.cis = get_paf_variance(paf.vcov.subset,
                             mr.subset),
  paf.bootstrap.cis = bootstrap_cis(mr.subset$root[44], bs.se.subset),
  paf.bootstrap.quantile = quantile_cis(bootstrap.subset)
)

paf.naive = c(
  paf = mr.naive$root[14],
  paf.cis = get_paf_variance(paf.vcov.naive,
                             mr.naive, rootnum = 14),
  paf.bootstrap.cis = bootstrap_cis(mr.naive$root[14], bs.se.naive),
  paf.bootstrap.quantile = quantile_cis(bootstrap, var="naive.paf")
)

# Subset results
{
  paf.200 = c(
    paf = mr.200$root[44],
    paf.cis = get_paf_variance(paf.vcov.200, mr.200),
    paf.bootstrap.cis = bootstrap_cis(mr.200$root[44], bs.se.200),
    paf.bootstrap.quantile = quantile_cis(bootstrap.200)
  )

  paf.500 = c(
    paf = mr.500$root[44],
    paf.cis = get_paf_variance(paf.vcov.500, mr.500),
    paf.bootstrap.cis = bootstrap_cis(mr.500$root[44], bs.se.500),
    paf.bootstrap.quantile = quantile_cis(bootstrap.500)
  )

  paf.1000 = c(
    paf = mr.1000$root[44],
    paf.cis = get_paf_variance(paf.vcov.1000, mr.1000),
    paf.bootstrap.cis = bootstrap_cis(mr.1000$root[44], bs.se.1000),
    paf.bootstrap.quantile = quantile_cis(bootstrap.1000)
  )

  paf.5000 = c(
    paf = mr.5000$root[44],
    paf.cis = get_paf_variance(paf.vcov.5000, mr.5000),
    paf.bootstrap.cis = bootstrap_cis(mr.5000$root[44], bs.se.5000),
    paf.bootstrap.quantile = quantile_cis(bootstrap.5000)
  )
}

# 11. Additional Deatils ----
# 10 year risk of AIDS/death in WIHS cohort
wihs.risk = c(
  naive.risk = mr.naive$root[43],
  naive.cis = get_paf_variance(paf.vcov.naive, mr.naive, rootnum = 43),
  y.risk = mr$root[43],
  y.cis = get_paf_variance(paf.vcov, mr, rootnum = 43),
  y0.risk = mr$root[42],
  y0.cis = get_paf_variance(paf.vcov, mr, rootnum = 42)
)

wihs.risk
