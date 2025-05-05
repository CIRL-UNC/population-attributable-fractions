library(tidyverse)

atlas = read_csv(file.choose())

table1 = atlas |>
  mutate(
    aa_indicator = Race == "Black/African American",
    idu = `Transmission Category` == "Injection drug use"
  ) |>
  rename(age = `Age Group`) |>
  group_by(aa_indicator, age, idu) |>
  summarise(
    cases1996 = sum(`Cases 1996`),
    cases2006 = sum(`Cases 2006`),
    cases2021 = sum(`Cases 2021`)
  ) |>
  ungroup() |>
  mutate(
    p1996 = cases1996 / sum(cases1996),
    p2006 = cases2006 / sum(cases2006),
    p2021 = cases2021 / sum(cases2021)
  )

atlas |>
  mutate(
    aa_indicator = Race == "Black/African American",
    idu = `Transmission Category` == "Injection drug use"
  ) |>
  rename(age = `Age Group`) |>
  group_by(aa_indicator, age) |>
  summarise(
    cases1996 = sum(`Cases 1996`),
    cases2006 = sum(`Cases 2006`),
    cases2021 = sum(`Cases 2021`)
  ) |>
  ungroup() |>
  mutate(
    p1996 = cases1996 / sum(cases1996),
    p2006 = cases2006 / sum(cases2006),
    p2021 = cases2021 / sum(cases2021)
  ) |> view()

table1_wihs = wihs |>
  group_by(aa_indicator, age, idu) |>
  count() |>
  ungroup() |>
  mutate(
    p = n / sum(n)
  )

table1_joined = left_join(table1, table1_wihs,
          by=c("aa_indicator", "age", "idu"), suffix=c("", "_wihs"))
