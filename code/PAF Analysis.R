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


# Right censored within 10 years:
wihs |> mutate(admin = t >= 10) |> group_by(eventtype, admin) |> summarise(n=n(), percent = n/1164)

# Assess Subsets (not included in paper)
set.seed(17)
wihs_subset = bind_rows(
  wihs |> filter(aa_indicator == 1),
  wihs |> filter(aa_indicator == 0) |> slice_sample(prop = 0.5)
)



# 2. Data Cleaning - CDC dataset ----

cdc_2008 = read_csv("data/cdc_2008.csv") |>
  rename(
    age = `AgeGroup`) |>
  dplyr::select(age, idu, aa_indicator, cases)

cdc_sample_n = function(n=500, seed=17, dataset=cdc_2008) {
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

cdc_sample_200 = cdc_sample_n(n=200, seed=18)
cdc_sample_500 = cdc_sample_n(n=500, seed=17)
cdc_sample_1000 = cdc_sample_n(n=1000, seed=17)
cdc_sample_5000 = cdc_sample_n(n=5000, seed=17)

# 3. Table 1 Generation ----
table1_data = bind_rows(
  target = cdc_2008 |> uncount(cases),
  study = wihs |> dplyr::select(age, aa_indicator, idu),
  .id = "pop"
)

table1(~factor(age) + factor(aa_indicator) + factor(idu) | pop, data = table1_data)
table1(~factor(age) + factor(aa_indicator) + factor(idu), data=cdc_sample_200 |> uncount(cases))

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

iosw_w = generate_iosw_model(cdc_2008, wihs, s ~ age * aa_indicator)
iosw_aw = generate_iosw_model(cdc_2008, wihs, s ~ age * aa_indicator * idu)

iosw_w.subset = generate_iosw_model(cdc_2008, wihs_subset, s ~ age * aa_indicator)
iosw_aw.subset = generate_iosw_model(cdc_2008, wihs_subset, s ~ age * aa_indicator * idu)

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
                          target_pop = cdc_2008,
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

bootstrap_runner = function(B=2000, bs_step, ...){
  sapply(1:B, function(x) bs_step(x, ...))
}



bootstrap = bootstrap_runner(B=2000, bootstrap_step, target_pop=cdc_2008, study_pop=wihs)
bootstrap.naive = bootstrap_runner(B=2000, bootstrap_step, target_pop=wihs, study_pop=wihs)
bootstrap.subset = bootstrap_runner(B=2000, bootstrap_step, target_pop=cdc_2008, study_pop=wihs_subset)

bootstrap.200 = bootstrap_runner(B=2000, bootstrap_step, target_pop=cdc_sample_200, study_pop=wihs)
bootstrap.500 = bootstrap_runner(B=2000, bootstrap_step, target_pop=cdc_sample_500, study_pop=wihs)
bootstrap.1000 = bootstrap_runner(B=2000, bootstrap_step, target_pop=cdc_sample_1000, study_pop=wihs)
bootstrap.5000 = bootstrap_runner(B=2000, bootstrap_step, target_pop=cdc_sample_5000, study_pop=wihs)


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

make_stacked_data = function(target, study) {
  bind_rows(
    target |> uncount(cases),
    study,
    .id="pop"
  ) |> mutate(
    population = ifelse(pop == "1", 1, 0)
  )
}

stacked_data = make_stacked_data(cdc_2008, wihs)


paf_mestimation = function(stacked_data=stacked_data,
                           solve = T){
  function(theta){
    population = stacked_data |>
      ungroup() |>
      transmute(
        population = ifelse(pop == "1", 1, 0)
      )|>
      pull(population)

    iosw.w = stacked_data |>
      model.matrix(object = ~age * aa_indicator)

    iosw.aw = stacked_data |>
      model.matrix(object = ~age * aa_indicator * idu)

    iptw.lc = with(stacked_data, plogis(theta[1] + theta[2] * aa_indicator + cd4knots %*% theta[3:6] + ageknots %*% theta[7:10]))

    iptw.weights = with(stacked_data, (1-idu) * theta[11] / (1-iptw.lc))

    iosw.w.weights = plogis(iosw.w %*% theta[12:21])
    iosw.aw.weights = plogis(iosw.aw %*% theta[22:41])

    iosw.w.odds = exp(iosw.w %*% theta[12:21])
    iosw.aw.odds = exp(iosw.aw %*% theta[22:41])

    ee_eval = stacked_data |>
      ungroup() |>
      mutate(
        across(c(ageknots, cd4knots), ~replace_na(.x, 0))
      ) |>
      transmute(
        theta1 = (population == 0) * (idu - iptw.lc),
        theta2 = (population == 0) * (idu - iptw.lc) * aa_indicator,
        theta3 = (population == 0) * (idu - iptw.lc) * cd4knots[,1],
        theta4 = (population == 0) * (idu - iptw.lc) * cd4knots[,2],
        theta5 = (population == 0) * (idu - iptw.lc) * cd4knots[,3],
        theta6 = (population == 0) * (idu - iptw.lc) * cd4knots[,4],
        theta7 = (population == 0) * (idu - iptw.lc) * ageknots[,1],
        theta8 = (population == 0) * (idu - iptw.lc) * ageknots[,2],
        theta9 = (population == 0) * (idu - iptw.lc) * ageknots[,3],
        theta10 = (population == 0) * (idu - iptw.lc) * ageknots[,4],

        theta11 = (population == 0) * (idu - theta[11]),
        theta12 = (population - iosw.w.weights),
        theta13 = (population - iosw.w.weights) * (age=="25-34"),
        theta14 = (population - iosw.w.weights) * (age=="35-44"),
        theta15 = (population - iosw.w.weights) * (age=="45-54"),
        theta16 = (population - iosw.w.weights) * (age=="55+"),
        theta17 = (population - iosw.w.weights) * aa_indicator,
        theta18 = (population - iosw.w.weights) * (age=="25-34") * aa_indicator,
        theta19 = (population - iosw.w.weights) * (age=="35-44") * aa_indicator,
        theta20 = (population - iosw.w.weights) * (age=="45-54") * aa_indicator,
        theta21 = (population - iosw.w.weights) * (age=="55+") * aa_indicator,

        theta22 = (population - iosw.aw.weights),
        theta23 = (population - iosw.aw.weights) * (age=="25-34"),
        theta24 = (population - iosw.aw.weights) * (age=="35-44"),
        theta25 = (population - iosw.aw.weights) * (age=="45-54"),
        theta26 = (population - iosw.aw.weights) * (age=="55+"),
        theta27 = (population - iosw.aw.weights) * aa_indicator,
        theta28 = (population - iosw.aw.weights) * idu,
        theta29 = (population - iosw.aw.weights) * (age=="25-34") * aa_indicator,
        theta30 = (population - iosw.aw.weights) * (age=="35-44") * aa_indicator,
        theta31 = (population - iosw.aw.weights) * (age=="45-54") * aa_indicator,
        theta32 = (population - iosw.aw.weights) * (age=="55+") * aa_indicator,
        theta33 = (population - iosw.aw.weights) * (age=="25-34") * idu,
        theta34 = (population - iosw.aw.weights) * (age=="35-44") * idu,
        theta35 = (population - iosw.aw.weights) * (age=="45-54") * idu,
        theta36 = (population - iosw.aw.weights) * (age=="55+") * idu,
        theta37 = (population - iosw.aw.weights) * aa_indicator * idu,
        theta38 = (population - iosw.aw.weights) * (age=="25-34") * aa_indicator * idu,
        theta39 = (population - iosw.aw.weights) * (age=="35-44") * aa_indicator * idu,
        theta40 = (population - iosw.aw.weights) * (age=="45-54") * aa_indicator * idu,
        theta41 = (population - iosw.aw.weights) * (age=="55+") * aa_indicator * idu,
        theta42 = (population == 0) * (d10 - theta[42])  * iptw.weights * iosw.w.odds,
        theta43 = (population == 0) * (d10 - theta[43]) * iosw.aw.odds,
        theta44 = theta[42] - theta[43] * (1 - theta[44])
      )

    if(solve){
      return(ee_eval |> colSums(na.rm=T))
    } else {
      return(ee_eval |> (\(.) {replace(., is.na(.), 0)})() |> as.matrix())
    }
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


theta1..10 = glm(idu ~ aa_indicator + cd4knots + ageknots, data=wihs,family = binomial()) |> coef()
theta11 = 0.4
theta12..21 = glm(population ~ age * aa_indicator, data=stacked_data) |> coef()
theta22..41 = glm(population ~ age * aa_indicator * idu, data=stacked_data) |> coef()
theta42 = 0.2
theta43 = 0.24
theta44 = 0.14

vec_start = c(theta1..10, theta11, theta12..21, theta22..41, theta42, theta43, theta44) |> unname()


mr = multiroot(paf_mestimation(make_stacked_data(cdc_2008, wihs), solve=T),
              vec_start)


## 8a. Test against other ----

#
{
mr.subset = multiroot(paf_mestimation(make_stacked_data(cdc_2008, wihs_subset), solve=T), vec_start)

mr.200 = multiroot(paf_mestimation(make_stacked_data(cdc_sample_200, wihs), solve=T), vec_start)

mr.500 =multiroot(paf_mestimation(make_stacked_data(cdc_sample_500, wihs), solve=T), vec_start)

mr.1000 = multiroot(paf_mestimation(make_stacked_data(cdc_sample_1000, wihs), solve=T), vec_start)

mr.5000 = multiroot(paf_mestimation(make_stacked_data(cdc_sample_5000, wihs), solve=T), vec_start)
}


# 9. M-Estimation - variance ----
paf_sandwich_variance = function(f, data, xval){
  n = nrow(data)
  ee_fun = f(stacked_data = data, solve=T)

  a0 = - gradient(f=ee_fun, x=xval) / n
  evaluated = f(stacked_data=data, solve=F)(xval)
  b0 = t(evaluated) %*% evaluated / n

  v0 = solve(a0) %*% b0 %*% t(solve(a0)) / n
  return(v0)
}

paf_sandwich_variance(paf_mestimation, stacked_data, mr$root) |> diag() |> sqrt()


get_paf_variance = function(variance_matrix, rf, rootnum=44){
  paf_mest = rf$root[rootnum]
  paf_mest.variance = diag(variance_matrix)[rootnum]
  paf_mest.cis = paf_mest + c(-1.96, 1.96) * sqrt(paf_mest.variance)

  return(c(paf=rf$root[rootnum], se=paf_mest.variance |> sqrt(), cis=paf_mest.cis))
}


get_paf_variance(
  paf_sandwich_variance(paf_mestimation,
                        make_stacked_data(cdc_2008, wihs),
                        mr$root),
  rf=mr
)

targets = read_csv(file.choose())

### 1996 ----
stacked_1996 = make_stacked_data(
  target=targets |> mutate(cases = cases1996),
  study=wihs
)

mr_1996 = multiroot(paf_mestimation(stacked_1996), vec_start)

get_paf_variance(
  paf_sandwich_variance(paf_mestimation,
                        stacked_1996, mr_1996$root),
  rf=mr_1996
)

### 2006 ----
stacked_2006 = make_stacked_data(
  target=targets |> mutate(cases = cases2006),
  study=wihs
)

mr_2006 = multiroot(paf_mestimation(stacked_2006), vec_start)

get_paf_variance(
  paf_sandwich_variance(paf_mestimation,
                        stacked_2006, mr_2006$root),
  rf=mr_2006
)

### 2021 ----
stacked_2021 = make_stacked_data(
  target=targets |> mutate(cases = cases2021),
  study=wihs
)

mr_2021 = multiroot(paf_mestimation(make_stacked_data(cdc_2021, wihs)), vec_start)

get_paf_variance(
  paf_sandwich_variance(paf_mestimation,
                        stacked_data, mr_2021$root),
  rf=mr_2021
)

### 2008 - HIV Incident population ----
stacked_2008 = make_stacked_data(
  target = target_2008 |>
    select(age, idu, aa_indicator, cases) |>
    drop_na(),
  study=wihs
)

mr_2008 = multiroot(paf_mestimation(stacked_2008), vec_start)

get_paf_variance(
  paf_sandwich_variance(paf_mestimation,
                        stacked_2008, mr_2008$root),
  rf=mr_2008,
  rootnum=42
)

### WIHS ----
stacked_wihs = make_stacked_data(
  target=targets |> mutate(cases = n_wihs),
  study=wihs
)

mr_wihs = multiroot(paf_mestimation(stacked_wihs), vec_start)

get_paf_variance(
  paf_sandwich_variance(paf_mestimation,
                        stacked_wihs, mr_wihs$root),
  rf=mr_wihs,
  rootnum=42
)


# 10. Results ----
paf = c(
  paf = mr$root[44],
  paf.cis = get_paf_variance(
    paf_sandwich_variance(paf_mestimation,
                          make_stacked_data(cdc_2008, wihs),
                          mr$root),
    rf=mr
  ),
  paf.bootstrap.cis = bootstrap_cis(mr$root[44], bs.se),
  paf.bootstrap.quantile = quantile_cis(bootstrap)
)


# Subset results
{
  paf.200 = c(
    paf = mr.200$root[44],
    paf.cis = get_paf_variance(
      paf_sandwich_variance(paf_mestimation,
                            make_stacked_data(cdc_sample_200, wihs),
                            mr.200$root),
      rf=mr.200
    ),
    paf.bootstrap.cis = bootstrap_cis(mr.200$root[44], bs.se.200),
    paf.bootstrap.quantile = quantile_cis(bootstrap.200)
  )

  paf.500 = c(
    paf = mr.500$root[44],
    paf.cis = get_paf_variance(
      paf_sandwich_variance(paf_mestimation,
                            make_stacked_data(cdc_sample_500, wihs),
                            mr.500$root),
      rf=mr.500
    ),
    paf.bootstrap.cis = bootstrap_cis(mr.500$root[44], bs.se.500),
    paf.bootstrap.quantile = quantile_cis(bootstrap.500)
  )

  paf.1000 = c(
    paf = mr.1000$root[44],
    paf.cis = get_paf_variance(
      paf_sandwich_variance(paf_mestimation,
                            make_stacked_data(cdc_sample_1000, wihs),
                            mr.1000$root),
      rf=mr.1000
    ),
    paf.bootstrap.cis = bootstrap_cis(mr.1000$root[44], bs.se.1000),
    paf.bootstrap.quantile = quantile_cis(bootstrap.1000)
  )

  paf.5000 = c(
    paf = mr.5000$root[44],
    paf.cis = get_paf_variance(
      paf_sandwich_variance(paf_mestimation,
                            make_stacked_data(cdc_sample_5000, wihs),
                            mr.5000$root),
      rf=mr.5000
    ),
    paf.bootstrap.cis = bootstrap_cis(mr.5000$root[44], bs.se.5000),
    paf.bootstrap.quantile = quantile_cis(bootstrap.5000)
  )
}

bind_rows(paf, paf.5000, paf.1000, paf.500, paf.200) |> view()
