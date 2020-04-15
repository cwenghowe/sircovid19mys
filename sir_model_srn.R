# ====================================================
# title: "Covid-19 SIR Modeling Malaysia"
# author: "Chan Weng Howe (UTM), Wan Nor Arifin (USM)"
# Fit both S and R
# Find right % of susceptible to match the data
# ====================================================


## Libraries
library(tidyverse)
library(deSolve)
library(lubridate)
library(ggpubr)
library(plotly)


## Data
# Sourced from local datasets.
mys_data <- read_csv("https://wnarifin.github.io/covid-19-malaysia/covid-19_my_full.csv") %>% select(1:10) %>% rename(Date = "date")  # uncomment
# mys_data <- read_csv("covid-19_my_full.csv") %>% select(1:10)  %>% rename(Date = "date")  # comment out later
# rename date to Date for compatability with original codes

# add active variable
mys_data$active <- with(mys_data, total_cases - total_deaths - total_recover)
names(mys_data)


## Functions

# SIR model
# S - Susceptible, I - Infected, R - Removed (recovered + death)
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- (-beta * I * S) / n
    dI <- ((beta * I * S) / n) - (gamma * I)
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# RSS 
RSS <- function(parameters) {
  names(parameters) <- c("beta","gamma")
  out <- ode(y=init, times=Day, func=SIR, parms=parameters)
  fit1 <- out[,3]
  fit2 <- out[,4]
  # sum((Active - fit1)^2) + sum((Removed - fit2)^2) # find min RSS on original scale 
  sum((log(Active) - log(fit1))^2) + sum((log(Removed) - log(fit2))^2)  # find min RSS on log scale
}


## Parameters

# Uncomment under each headers to perform SIR for the period

# Pre MCO #
# name = "Pre MCO"
# start_date  <- "2020-03-01"
# end_date    <- "2020-03-17"

# MCO #
# name = "MCO"
# start_date  <- "2020-03-18"
# end_date    <- "2020-04-28"

# MCO 1 #
# name = "MCO 1"
# start_date  <- "2020-03-18"
# end_date    <- "2020-03-31"

# MCO 2 #
# name = "MCO 2"
# start_date  <- "2020-04-01"
# end_date    <- "2020-04-14"

# MCO 2a #
name = "MCO 2a"
start_date  <- "2020-04-01"
end_date    <- "2020-04-07"

# MCO 2b #
# name = "MCO 2b"
# start_date  <- "2020-04-08"
# end_date    <- "2020-04-14"

# Basic info from data
Infected   <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_cases)
Recovered  <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_recover)
Death      <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_deaths)
Active     <- Infected - Recovered - Death
Removed    <- Recovered + Death
Date       <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(Date)
Day        <- 1:(length(Date))

# Initialization

# Total population
N <- 32.68E6  # 4th quarter, 2019
# initial susceptible population of Malaysia
p_start = 0.0001
p_end = 0.1  # max p of population
p_step = 0.00005
susceptible = seq(p_start, p_end, p_step)
N <- N * susceptible
inits <- data.frame(S=N-Infected[1]-Recovered[1]-Death[1], I=Infected[1]-Death[1]-Recovered[1], R=Recovered[1]+Death[1])

## Run SIR model

# SIR parameters
parameters_values <- c(1/2, 1/14)  # set reasonable start, R0 = 7
parameters_values_lower <- c(1/100, 1/42)  # days recover 6 weeks
parameters_values_upper <- c(1, 1/11)  # updated to min 11 days

t = 1:max(Day)
Opts = vector("list", 10)
for (i in 1:length(N)) {
  n = N[i]
  init = setNames(as.numeric(inits[i,]), c("S", "I", "R"))
  Opt_ = optim(parameters_values, RSS, method = "L-BFGS-B", lower = parameters_values_lower, upper = parameters_values_upper)
  Opts[i]  = list(Opt_)
}
# Opt <- optim(parameters_values, RSS, method = "L-BFGS-B", lower = parameters_values_lower, upper = parameters_values_upper)
Opt_value = sapply(Opts, function(x) x$value)
Opt_min_loc = which(Opt_value == min(Opt_value))
diff(Opt_value)

susceptible[Opt_min_loc]/denom
N[Opt_min_loc]  # Susceptible population
Opt = Opts[[Opt_min_loc]]
Opt$message  # make sure converge
Opt$value
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
R0 <- (Opt_par['beta']/Opt_par['gamma']); names(R0) <- "R0"
parameters_values_lower; parameters_values; parameters_values_upper  # just to check whether values stuck
cat("beta = ", Opt_par[['beta']], ", infectious contact rate (/person/day)\n",
    "gamma = ", Opt_par[['gamma']], ", recovery rate (/day)\n",
    "R_0 = ", R0, " number infected/person\n",
    "Recovery days = ", 1/Opt_par[['gamma']], " days",
    sep = "")


## Fit Data

# time in days for predictions
t = 1:max(Day)
n = N[Opt_min_loc]
init = setNames(as.numeric(inits[Opt_min_loc,]), c("S", "I", "R"))
# get the fitted values from our SIR model
fitted_projected <- data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add a Date & Active
fitted_projected$Date = Date
fitted_projected$A = Active
fitted_projected$Rm = Removed
fitted_projected

# Fit measure, log scales
tss1 = sum((log(fitted_projected$A) - mean(log(fitted_projected$A)))^2); tss1
rss1 = sum((log(fitted_projected$A) - log(fitted_projected$I))^2); rss1
R2_1 = 1 - (rss1 / tss1); R2_1
tss2 = sum((log(fitted_projected$Rm) - mean(log(fitted_projected$Rm)))^2); tss2
rss2 = sum((log(fitted_projected$Rm) - log(fitted_projected$R))^2); rss2
R2_2 = 1 - (rss2 / tss2); R2_2

# Plots

# color settings
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", "Observed Active" = "orange")

# plot fit the data
fitplot1 = ggplot(fitted_projected, aes(x = Date)) + geom_line(aes(y = I, color = "Infectious")) + 
  geom_point(aes(y = A, color = "Observed Active")) + geom_line(aes(y = R, color = "Recovered")) + 
  geom_point(aes(y = Rm, color = "Observed Active")) +
  labs(y = "Infectious", title = paste("COVID-19 fitted vs Observed Active Cases in Malaysia,", name),
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", n, "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n")) +
  scale_colour_manual(values = colors)
fitplot1
# ggsave(paste0("plots_srn/fit", name, ".png"), width = 16, height = 10)

# plot fit the data, in log10
fitplot1_log = ggplot(fitted_projected, aes(x = Date)) + geom_line(aes(y = I, color = "Infectious")) + 
  geom_point(aes(y = A, color = "Observed Active")) + geom_line(aes(y = R, color = "Recovered")) + 
  geom_point(aes(y = Rm, color = "Observed Active")) + scale_y_log10() + 
  labs(y = "Infectious", title = paste("COVID-19 fitted vs Observed Active Cases in Malaysia,", name, "log10"),
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", n, "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n")) +
  scale_colour_manual(values = colors)  
fitplot1_log
# ggsave(paste0("plots_srn/fit", name, "_log.png"), width = 16, height = 10)


## Projected Data

# last date to project
last_date = "2020-09-30"
# time in days for predictions
t = 1:as.integer(ymd(last_date) + 1 - ymd(start_date))
# get the fitted values from our SIR model
fitted_projected <- data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add a Date & Active
fitted_projected$Date = ymd(start_date) + days(t - 1)
fitted_projected$A = c(Active, rep(NA, length(t) - length(Active)))
fitted_projected$Rm = c(Removed, rep(NA, length(t) - length(Active)))
head(fitted_projected, 10); tail(fitted_projected, 10)
# date peak
max_I = which(round(fitted_projected$I) == round(max(fitted_projected$I)))
max_date = fitted_projected$Date[max_I]

total_infected = fitted_projected$I + fitted_projected$R
total_infected[14]
max(total_infected)

# color settings
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", "Observed Active" = "orange")

# plot projection data
sirplot1 = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_point(aes(y = A, color = "Observed Active")) +
  geom_point(aes(y = Rm, color = "Observed Active")) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Cumulative incidence", title = paste("COVID-19 SIR model Malaysia,", name), 
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", n, "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n")) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = 30E6, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot1
# ggsave(paste0("plots_srn/sir", name, ".png"), width = 16, height = 10)

# plot projection data, in log10
sirplot1_log = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_point(aes(y = A, color = "Observed Active")) +
  geom_point(aes(y = Rm, color = "Observed Active")) +
  scale_y_log10(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Cumulative incidence", title = paste("COVID-19 SIR model Malaysia,", name, "log10"), 
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", n, "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n")) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = 1E5, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot1_log
# ggsave(paste0("plots_srn/sir", name, "_log.png"), width = 16, height = 10)
