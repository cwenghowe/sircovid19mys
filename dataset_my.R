# alternative local datasets
# already corrected the data points <= 29/2/2020 to official MOH announcement
# also state data from earliest offiical MOH announcement

library(tidyverse)
covid_my_full = read_csv("https://wnarifin.github.io/covid-19-malaysia/covid-19_my_full.csv")
covid_my_state = read_csv("https://wnarifin.github.io/covid-19-malaysia/covid-19_my_state.csv")