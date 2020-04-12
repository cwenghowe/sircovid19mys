library(tidyverse)
library(deSolve)
library(lubridate)
library(ggpubr)
library(plotly)

mys_data <- read_csv("https://wnarifin.github.io/covid-19-malaysia/covid-19_my_full.csv") %>% select(1:10) %>% rename(Date = "date")  # uncomment
# mys_data <- read_csv("covid-19_my_full.csv") %>% select(1:10)  %>% rename(Date = "date")  # comment out later
# rename date to Date for compatability with original codes

# add active variable
mys_data$active <- with(mys_data, total_cases - total_deaths - total_recover)
names(mys_data)

# SIR model
# S - Susceptible, I - Infected, R - Removed (recovered + death)
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- (-beta * I * S) / N
    dI <- ((beta * I * S) / N) - (gamma * I)
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# RSS 
RSS <- function(parameters) {
  names(parameters) <- c("beta","gamma")
  out <- ode(y=init, times=Day, func=SIR, parms=parameters)
  fit <- out[,3]
  sum((Active - fit)^2)  # find min RSS on original scale
  # sum((log(Active) - log(fit))^2)  # find min RSS on log scale
}

# Total population
N <- 32.68E6  # 4th quarter, 2019

# range of period to be analyzed
# Pre MCO
name = "Pre MCO"
start_date  <- "2020-03-01"
end_date    <- "2020-03-17"
# percentage of population implement social distancing (assumption)
social_distance <- 0  # no social distance

# MCO 1
# name = "MCO 1"
# start_date  <- "2020-03-18" 
# end_date    <- "2020-03-31"
# social_distance <- 0.7

# MCO 2
# name = "MCO 2"
# start_date  <- "2020-04-01" 
# end_date    <- "2020-04-14"
# social_distance <- 0.9

# initial population of Malaysia considering exclude the population who implemented social distancing
N <- N * (1 - social_distance)

Infected   <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_cases)
Recovered  <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_recover)
Death      <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_deaths)
Active     <- Infected - Recovered - Death
Date       <- mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(Date)
Day        <- 1:(length(Date))

# Initialization
init <- c(S=N-Infected[1]-Recovered[1]-Death[1], I=Infected[1]-Death[1]-Recovered[1], R=Recovered[1]+Death[1])

# SIR parameters
parameters_values <- c(1/7, 1/21)  # seem no longer stuck at local minima if take middle values?
parameters_values_lower <- c(1/14, 1/42)
parameters_values_upper <- c(1, 1/7)

# SIR model
Opt <- optim(parameters_values, RSS, method = "L-BFGS-B", lower = parameters_values_lower, upper = parameters_values_upper)
Opt$message  # make sure converge
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
R0 <- (Opt_par['beta']/Opt_par['gamma']); names(R0) <- "R0"
parameters_values_lower; parameters_values; parameters_values_upper  # just to check whether values stuck
cat("beta = ", Opt_par[[1]], ", infectious contact rate (/person/day)\n",
    "gamma = ", Opt_par[[2]], ", recovery rate (/day)\n",
    "R_0 = ", R0, " number infected/person",
    sep = "")

# Projected Data & Model fit

# fit data
# time in days for predictions
t = 1:max(Day)
# get the fitted values from our SIR model
fitted_projected <- data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add a Date & Active
fitted_projected$Date = Date
fitted_projected$A = Active
fitted_projected

# fit measures
tss = sum((fitted_projected$A - mean(fitted_projected$A))^2); tss
rss = sum((fitted_projected$A - fitted_projected$I)^2); rss
R2 = 1 - (rss / tss); R2
# or in log
tss = sum((log(fitted_projected$A) - mean(log(fitted_projected$A)))^2); tss
rss = sum((log(fitted_projected$A) - log(fitted_projected$I))^2); rss
R2 = 1 - (rss / tss); R2
# or better this way? bcs VAR(X) = mean(X) if X is distributed as Poisson
mean(fitted_projected$I) / mean(fitted_projected$A)

# Plots
# color settings
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", "Observed Active" = "orange")

# plot fit the data
fitplot1 = ggplot(fitted_projected, aes(x = Date)) + geom_line(aes(y = I, color = "Infectious")) + 
  geom_point(aes(y = A, color = "Observed Active")) +
  labs(y = "Infectious", title = paste("COVID-19 fitted vs Observed Active Cases in Malaysia,", name, "log10")) +
  scale_colour_manual(values = colors, name = "")
fitplot1
ggsave("sir_mco2.png", width = 16, height = 10)

# plot fit the data, in log10
fitplot1_log = ggplot(fitted_projected, aes(x = Date)) + geom_line(aes(y = I, color = "Infectious")) + 
  geom_point(aes(y = A, color = "Observed Active")) + scale_y_log10() + 
  labs(y = "Infectious", title = paste("COVID-19 fitted vs Observed Active Cases in Malaysia,", name, "log10")) +
  scale_colour_manual(values = colors, name = "")  
fitplot1_log
ggsave("sir_mco2_log.png", width = 16, height = 10)

# projection
# last date to project
last_date <- "2020-09-30"
# time in days for predictions
t = 1:as.integer(ymd(last_date) + 1 - ymd(start_date))
# get the fitted values from our SIR model
fitted_projected <- data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add a Date & Active
fitted_projected$Date = ymd(start_date) + days(t - 1)
fitted_projected$A = c(Active, rep(NA, length(t) - length(Active)))
head(fitted_projected, 10); tail(fitted_projected, 10)
# date peak
max_I = which(round(fitted_projected$I) == round(max(fitted_projected$I)))
max_date = fitted_projected$Date[max_I]

# color settings
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", "Observed Active" = "orange")

# plot projection data
sirplot1 = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_point(aes(y = A, color = "Observed Active")) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Cumulative incidence", title = paste("COVID-19 SIR model Malaysia,", name), 
       color = paste0("R square = ", round(R2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Max Active = ", round(max(fitted_projected$I)), "\n")) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = 30E6, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot1
# ggsave("sir_mco2.png", width = 16, height = 10)

# plot projection data, in log10
sirplot1_log = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_point(aes(y = A, color = "Observed Active")) +
  scale_y_log10(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Cumulative incidence", title = paste("COVID-19 SIR model Malaysia,", name, "log10"), 
       color = paste0("R square = ", round(R2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Max Active = ", round(max(fitted_projected$I)), "\n")) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = 1E5, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot1_log
# ggsave("sir_mco2_log.png", width = 16, height = 10)
