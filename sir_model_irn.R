# ===================================================#
# title: "Covid-19 SIR Modeling Malaysia"
# author: "Chan Weng Howe (UTM), Wan Nor Arifin (USM)"
# Fit both I and R
# Find right n of susceptible to match the data
# ===================================================#


## Libraries
library(tidyverse)
library(deSolve)
library(lubridate)
library(ggpubr)
library(plotly)


## Data
# Sourced from local datasets.
mys_data = read_csv("https://wnarifin.github.io/covid-19-malaysia/covid-19_my_full.csv") %>% select(1:10) %>% rename(Date = "date")
# mys_data = read_csv("covid-19_my_full.csv") %>% select(1:10)  %>% rename(Date = "date")  # comment out to use local data
# rename date to Date for compatability with original codes

# add active variable
mys_data$active = with(mys_data, total_cases - total_deaths - total_recover)
names(mys_data)


## Functions

# SIR model
# S - Susceptible, I - Infected, R - Removed (recovered + death), n - initial susceptible population
SIR = function(time, state, parameters) {
  par = as.list(c(state, parameters))
  with(par, {
    dS = (-beta * I * S) / n
    dI = ((beta * I * S) / n) - (gamma * I)
    dR = gamma * I
    list(c(dS, dI, dR))
  })
}

# RSS 
RSS = function(parameters) {
  names(parameters) = c("beta","gamma")
  out = ode(y=init, times=Day, func=SIR, parms=parameters)
  fit1 = out[,3]
  fit2 = out[,4]
  w1 = 4
  w2 = 1
  # w1*sum((Active - fit1)^2) + w2*sum((Removed - fit2)^2) # find min RSS on original scale, weighted
  w1*sum((log(Active) - log(fit1))^2) + w2*sum((log(Removed) - log(fit2))^2)  # find min RSS on log scale, weighted
  # give weight to less well predicted curve
}


## Parameters

# Uncomment under each header to perform SIR for the period

# Pre MCO #
# name = "Pre MCO"
# start_date  = "2020-03-01"
# end_date    = "2020-03-17"

# MCO All -> Today #
# name = "MCO All"
# start_date  = "2020-03-18"
# end_date    = today()  # will analyze up to max available date

# MCO All minus 1 week #
# name = "MCO All minus 1 week"
# start_date  = "2020-03-18"
# end_date    = "2020-04-07"

# MCO week 2 -> Today # to take into account MCO effect after 1 week
name = "MCO week 2 -> Today"
start_date  = "2020-03-25"
end_date    = today()
# So far, this provides good balance between MCOs

# MCO 1 #
# name = "MCO 1"
# start_date  = "2020-03-18"
# end_date    = "2020-03-31"

# MCO 2 #
# name = "MCO 2"
# start_date  = "2020-04-01"
# end_date    = "2020-04-14"

# MCO 2 -> Today #
# name = "MCO 2 -> Today"
# start_date  = "2020-04-01"
# end_date    = today()

# MCO 2a #
# name = "MCO 2a"
# start_date  = "2020-04-01"
# end_date    = "2020-04-07"

# MCO 2b #
# name = "MCO 2b"
# start_date  = "2020-04-08"
# end_date    = "2020-04-14"

# Basic info from data
Infected   = mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_cases)
Recovered  = mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_recover)
Death      = mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(total_deaths)
Active     = Infected - Recovered - Death
Removed    = Recovered + Death
Date       = mys_data %>% filter(Date>=ymd(start_date), Date<=ymd(end_date)) %>% pull(Date)
Day        = 1:(length(Date))


## Run Optimization for SIR model

# SIR parameters
parameters_values       = c(1/2, 1/14)  # set reasonable start, R0 = 7
parameters_values_lower = c(1/100, 1/19)  # days recover 6 weeks, reduce bound bcs not many severe, observed show quick recovery
parameters_values_upper = c(1, 1/11)  # updated to min 11 days

# Placeholder to ind optimal susceptible population
N_in_steps = data.frame(step=1:7, N=rep(NA,7), Loc=rep(NA,7))

# Initial values
N           = 32.68E6/2  # 4th quarter, 2019
Opt_min_loc = 1  # Optimum minimum n location in output vector
step_i      = 1  # initialize step

# Step 1-7
for (step_i in 1:7) {
  cat("=== Step ", step_i, ": Finding optimal values of beta, gamma and n ===\n", sep ="")
  N = N[Opt_min_loc]  # Susceptible population from previous Step
  # initial susceptible population of Malaysia
  p_start = 0.1
  p_end = 2  # max p of population, also include up direction
  p_step = 0.1
  susceptible = seq(p_start, p_end, p_step)
  N = N * susceptible
  inits = data.frame(S=N-Infected[1]-Recovered[1]-Death[1], I=Infected[1]-Death[1]-Recovered[1], R=Recovered[1]+Death[1])
  Opts = vector("list", length(N))
  for (i in 1:length(N)) {
    n = N[i]
    init = setNames(as.numeric(inits[i,]), c("S", "I", "R"))
    Opt_ = optim(parameters_values, RSS, method = "L-BFGS-B", lower = parameters_values_lower, upper = parameters_values_upper)
    Opts[i]  = list(Opt_)
  }
  Opt_value = sapply(Opts, function(x) x$value)
  Opt_min_loc = which(Opt_value == min(Opt_value))
  N_in_steps[step_i, "N"] = N[Opt_min_loc]
  N_in_steps[step_i, "Loc"] = Opt_min_loc
  cat("=== Found n = ", N[Opt_min_loc], " at location ", Opt_min_loc, " in vector N ===\n\n", sep = "")
  step_i = step_i + 1
  if (step_i == 8) {
    cat("=== Finalizing results =====\n")
    cat("============================\n")
    print(N_in_steps)
  }
}

# Saving optimized parameters
Opt = Opts[[Opt_min_loc]]
Opt$message  # make sure converge
Opt_par = setNames(Opt$par, c("beta", "gamma"))
R0 = (Opt_par['beta']/Opt_par['gamma']); names(R0) = "R0"
parameters_values_lower; parameters_values; parameters_values_upper  # just to check whether values in range
# Print final parameters
cat("beta = ", Opt_par[['beta']], ", infectious contact rate (/person/day)\n",
    "gamma = ", Opt_par[['gamma']], ", recovery rate (/day)\n",
    "R_0 = ", R0, " number infected/person\n",
    "Recovery days = ", 1/Opt_par[['gamma']], " days",
    sep = "")


## Fit Data

# time in days for fitting
t = 1:max(Day)
n = N[Opt_min_loc]
init = setNames(as.numeric(inits[Opt_min_loc,]), c("S", "I", "R"))
# get the fitted values from our SIR model
fitted_projected = data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add Date, Active, Removed
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
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", 
           "Observed Active" = "orange", "Observed Recovered" = "blue")

# plot fit the data
fitplot1 = ggplot(fitted_projected, aes(x = Date)) + geom_line(aes(y = I, color = "Infectious")) + 
  geom_point(aes(y = A, color = "Observed Active")) + geom_line(aes(y = R, color = "Recovered")) + 
  geom_point(aes(y = Rm, color = "Observed Recovered")) +
  labs(y = "Number of cases", title = paste("COVID-19 SIR model Malaysia, fitted and observed", name),
       subtitle = paste("Projection from data:", start_date, "to", fitted_projected$Date[max(Day)]),
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", round(n), "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n")) +
  scale_colour_manual(values = colors)
fitplot1
# uncomment ggsave line to save plot in png format, make sure folder "plots_srn" is created
# ggsave(paste0("plots_srn/fit", name, ".png"), width = 12, height = 9)

# plot fit the data, in log10
fitplot1_log = ggplot(fitted_projected, aes(x = Date)) + geom_line(aes(y = I, color = "Infectious")) + 
  geom_point(aes(y = A, color = "Observed Active")) + geom_line(aes(y = R, color = "Recovered")) + 
  geom_point(aes(y = Rm, color = "Observed Recovered")) + scale_y_log10() + 
  labs(y = "Number of cases", title = paste("COVID-19 SIR model Malaysia, fitted and observed,", name, "log10"),
       subtitle = paste("Projection from data:", start_date, "to", fitted_projected$Date[max(Day)]),
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", round(n), "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n")) +
  scale_colour_manual(values = colors)  
fitplot1_log
# uncomment ggsave line to save plot in png format, make sure folder "plots_srn" is created
# ggsave(paste0("plots_srn/fit", name, "_log.png"), width = 12, height = 9)


## Projected Data

# last date to project
last_date = "2020-09-30"
# time in days for predictions
t = 1:as.integer(ymd(last_date) + 1 - ymd(start_date))
# get the fitted values from our SIR model
fitted_projected = data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add add Date, Active, Removed
fitted_projected$Date = ymd(start_date) + days(t - 1)
fitted_projected$A = c(Active, rep(NA, length(t) - length(Active)))
fitted_projected$Rm = c(Removed, rep(NA, length(t) - length(Active)))
head(fitted_projected, 10); tail(fitted_projected, 10)
# date of peak active cases
# max_I = which(round(fitted_projected$I) == round(max(fitted_projected$I)))  # at times this works better
max_I = which(fitted_projected$I == max(fitted_projected$I))
max_date = fitted_projected$Date[max_I]
# add cumulative infected cases
fitted_projected$total_infected = fitted_projected$I + fitted_projected$R
# predicted new cases today
new_today = (fitted_projected[fitted_projected$Date == today(), ] - fitted_projected[fitted_projected$Date == today()-1, ])$total_infected
# maximum cumulative cases, date. May add to plot.
fitted_projected$Date[min(which(round(fitted_projected$total_infected) == max(round(fitted_projected$total_infected))))]
fitted_projected[min(which(round(fitted_projected$total_infected) == max(round(fitted_projected$total_infected)))),]

# color settings
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", "Observed Active" = "orange", "Observed Recovered" = "blue")

# plot projection data
sirplot1 = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_point(aes(y = A, color = "Observed Active")) +
  geom_point(aes(y = Rm, color = "Observed Recovered")) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Number of cases", title = paste("COVID-19 SIR model Malaysia,", name), 
       subtitle = paste("Projection from data:", start_date, "to", fitted_projected$Date[max(Day)]),
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", round(n), "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n",
                      "Maximum Total Infected = ", round(max(fitted_projected$total_infected)))) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = n*1.3, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  geom_vline(xintercept = as.numeric(as.Date(today())), linetype = "dotted", color = "red") +
  annotate(geom = "text", x = as.Date(today())+25, y = n*1.2, 
           label = paste0("Today's Prediction (", format(today(), "%d/%m/%y"), ")\n",
                          "Total Cases = ", round(fitted_projected[fitted_projected$Date == today(), "total_infected"]),
                          "\nNew Cases = ", round(new_today)), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot1
# uncomment ggsave line to save plot in png format, make sure folder "plots_srn" is created
# ggsave(paste0("plots_srn/sir", name, ".png"), width = 12, height = 9)

# plot projection data, in log10
sirplot1_log = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_point(aes(y = A, color = "Observed Active")) +
  geom_point(aes(y = Rm, color = "Observed Recovered")) +
  scale_y_log10(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Number of cases", title = paste("COVID-19 SIR model Malaysia,", name, "log10"), 
       subtitle = paste("Projection from data:", start_date, "to", fitted_projected$Date[max(Day)]),
       color = paste0("R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n",
                      "Susceptible = ", round(n), "\n",
                      "Peak Active = ", round(max(fitted_projected$I)), "\n",
                      "Maximum Total Infected = ", round(max(fitted_projected$total_infected)))) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = n*1.3, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  geom_vline(xintercept = as.numeric(as.Date(today())), linetype = "dotted", color = "red") +
  annotate(geom = "text", x = as.Date(today())+25, y = n*0.7, 
           label = paste0("Today's Prediction (", format(today(), "%d/%m/%y"), ")\n",
                          "Total Cases = ", round(fitted_projected[fitted_projected$Date == today(), "total_infected"]),
                          "\nNew Cases = ", round(new_today)), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot1_log
# uncomment ggsave line to save plot in png format, make sure folder "plots_srn" is created
# ggsave(paste0("plots_srn/sir", name, "_log.png"), width = 12, height = 9)

## Projection + Asymptomatic

# if 5%, 50%  asymptomatic
# https://www.cebm.net/covid-19/covid-19-what-proportion-are-asymptomatic/

# 1st situation --- #
# get the fitted values from our SIR model
p_asym = .05  # China figure, 5%
n1 = n/(1-p_asym)  # inflate by %
init = c(S=n1-Infected[1]-Recovered[1]-Death[1], I=Infected[1]-Death[1]-Recovered[1], R=Recovered[1]+Death[1])
fitted_projected1 = data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add a Date & Active
fitted_projected1$Date = ymd(start_date) + days(t - 1)
fitted_projected1$A = c(Active, rep(NA, length(t) - length(Active)))
fitted_projected1$Rm = c(Removed, rep(NA, length(t) - length(Active)))
head(fitted_projected1, 10); tail(fitted_projected1, 10)
# add cumulative cases
fitted_projected1$total_infected = fitted_projected1$I + fitted_projected1$R
# maximum cumulative cases, date. May add to plot.
fitted_projected1$Date[min(which(round(fitted_projected1$total_infected) == max(round(fitted_projected1$total_infected))))]
fitted_projected1[min(which(round(fitted_projected1$total_infected) == max(round(fitted_projected1$total_infected)))),]
# --- 2nd situation --- #
# get the fitted values from our SIR model
p_asym = .5  # Iceland figure, 50%
n2 = n/(1-p_asym)  # inflate by %
init = c(S=n2-Infected[1]-Recovered[1]-Death[1], I=Infected[1]-Death[1]-Recovered[1], R=Recovered[1]+Death[1])
fitted_projected2 = data.frame(ode(y=init, times=t, func=SIR, parms=Opt_par))
# add a Date & Active
fitted_projected2$Date = ymd(start_date) + days(t - 1)
fitted_projected2$A = c(Active, rep(NA, length(t) - length(Active)))
fitted_projected2$Rm = c(Removed, rep(NA, length(t) - length(Active)))
head(fitted_projected2, 10); tail(fitted_projected2, 10)
# add cumulative cases
fitted_projected2$total_infected = fitted_projected2$I + fitted_projected2$R
# maximum cumulative cases, date. May add to plot.
fitted_projected2$Date[min(which(round(fitted_projected2$total_infected) == max(round(fitted_projected2$total_infected))))]
fitted_projected2[min(which(round(fitted_projected2$total_infected) == max(round(fitted_projected2$total_infected)))),]

# color settings
colors = c("Susceptible" = "black", "Recovered" = "green", "Infectious" = "red", 
           "Susceptible + 5% (dashed)" = "black", "Recovered + 5% (dashed)" = "green", "Infectious + 5% (dashed)" = "red", 
           "Susceptible + 50% (dot-dash)" = "black", "Recovered + 50% (dot-dash)" = "green", "Infectious + 50% (dot-dash)" = "red", 
           "Observed Active" = "orange", "Observed Recovered" = "blue")

# plot projection data
sirplot2 = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_line(data=fitted_projected1, aes(y = I, color = "Infectious + 5% (dashed)"), linetype = "dashed") + 
  geom_line(data=fitted_projected1, aes(y = S, color = "Susceptible + 5% (dashed)"), linetype = "dashed") + 
  geom_line(data=fitted_projected1, aes(y = R, color = "Recovered + 5% (dashed)"), linetype = "dashed") +
  geom_line(data=fitted_projected2, aes(y = I, color = "Infectious + 50% (dot-dash)"), linetype = "dotdash") + 
  geom_line(data=fitted_projected2, aes(y = S, color = "Susceptible + 50% (dot-dash)"), linetype = "dotdash") + 
  geom_line(data=fitted_projected2, aes(y = R, color = "Recovered + 50% (dot-dash)"), linetype = "dotdash") +
  geom_point(aes(y = A, color = "Observed Active")) +
  geom_point(aes(y = Rm, color = "Observed Recovered")) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Cumulative incidence", title = paste("COVID-19 SIR model Malaysia,", name), 
       subtitle = paste("Projection from data:", start_date, "to", fitted_projected$Date[max(Day)]),
       color = paste0("Model fit:\n",
                      "R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n\n",
                      "SIR parameters:\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n\n",
                      "Susceptible Individuals:\n",
                      "No asymptomatic = ", round(n), "\n",
                      "Asymptomatic 5% = ", round(n1), "\n",
                      "Asymptomatic 50% = ", round(n2), "\n\n",
                      "Peak Active Cases:\n",
                      "No asymptomatic = ", round(max(fitted_projected$I)), "\n",
                      "Asymptomatic 5% = ", round(max(fitted_projected1$I)), "\n",
                      "Asymptomatic 50% = ", round(max(fitted_projected2$I)), "\n\n",
                      "Maximum Total Infected:\n",
                      "No asymptomatic = ", round(max(fitted_projected$total_infected)), "\n",
                      "Asymptomatic 5% = ", round(max(fitted_projected1$total_infected)), "\n",
                      "Asymptomatic 50% = ", round(max(fitted_projected2$total_infected)), "\n")) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = n2*1.3, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  geom_vline(xintercept = as.numeric(as.Date(today())), linetype = "dotted", color = "red") +
  annotate(geom = "text", x = as.Date(today())+25, y = n2*1.2, 
           label = paste0("Today's Prediction (", format(today(), "%d/%m/%y"), ")\n",
                          "Total Cases = ", round(fitted_projected[fitted_projected$Date == today(), "total_infected"]),
                          "\nNew Cases = ", round(new_today)), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot2
# uncomment ggsave line to save plot in png format, make sure folder "plots_srn" is created
# ggsave(paste0("plots_srn/sir", name, "_asymp.png"), width = 12, height = 9)

# plot projection data, in log10
sirplot2_log = ggplot(fitted_projected, aes(x = Date)) + 
  geom_line(aes(y = I, color = "Infectious")) + 
  geom_line(aes(y = S, color = "Susceptible")) + 
  geom_line(aes(y = R, color = "Recovered")) +
  geom_line(data=fitted_projected1, aes(y = I, color = "Infectious + 5% (dashed)"), linetype = "dashed") + 
  geom_line(data=fitted_projected1, aes(y = S, color = "Susceptible + 5% (dashed)"), linetype = "dashed") + 
  geom_line(data=fitted_projected1, aes(y = R, color = "Recovered + 5% (dashed)"), linetype = "dashed") +
  geom_line(data=fitted_projected2, aes(y = I, color = "Infectious + 50% (dot-dash)"), linetype = "dotdash") + 
  geom_line(data=fitted_projected2, aes(y = S, color = "Susceptible + 50% (dot-dash)"), linetype = "dotdash") + 
  geom_line(data=fitted_projected2, aes(y = R, color = "Recovered + 50% (dot-dash)"), linetype = "dotdash") +
  geom_point(aes(y = A, color = "Observed Active")) +
  geom_point(aes(y = Rm, color = "Observed Recovered")) +
  scale_y_log10(labels = scales::comma) +
  scale_x_date(date_breaks = "14 day", date_labels = "%d/%m/%y") + 
  scale_colour_manual(values = colors) +
  labs(y = "Cumulative incidence", title = paste("COVID-19 SIR model Malaysia,", name), 
       subtitle = paste("Projection from data:", start_date, "to", fitted_projected$Date[max(Day)]),
       color = paste0("Model fit:\n",
                      "R square 1 = ", round(R2_1,3), "\n",
                      "R square 2 = ", round(R2_2,3), "\n\n",
                      "SIR parameters:\n",
                      "R0 = ", round(R0, 3), "\n",
                      "beta = ", round(Opt_par[1], 3), "\n",
                      "gamma = ", round(Opt_par[2], 3), "\n\n",
                      "Susceptible Individuals:\n",
                      "No asymptomatic = ", round(n), "\n",
                      "Asymptomatic 5% = ", round(n1), "\n",
                      "Asymptomatic 50% = ", round(n2), "\n\n",
                      "Peak Active Cases:\n",
                      "No asymptomatic = ", round(max(fitted_projected$I)), "\n",
                      "Asymptomatic 5% = ", round(max(fitted_projected1$I)), "\n",
                      "Asymptomatic 50% = ", round(max(fitted_projected2$I)), "\n\n",
                      "Maximum Total Infected:\n",
                      "No asymptomatic = ", round(max(fitted_projected$total_infected)), "\n",
                      "Asymptomatic 5% = ", round(max(fitted_projected1$total_infected)), "\n",
                      "Asymptomatic 50% = ", round(max(fitted_projected2$total_infected)), "\n")) +
  geom_vline(xintercept = as.numeric(as.Date(max_date)), linetype = "dotted") +
  annotate(geom = "text", x = as.Date(max_date)+20, y = n2*1.3, 
           label = paste0("Peak on ", format(max_date, "%d/%m/%y")), angle = 0) +
  geom_vline(xintercept = as.numeric(as.Date(today())), linetype = "dotted", color = "red") +
  annotate(geom = "text", x = as.Date(today())+25, y = n2*0.7, 
           label = paste0("Today's Prediction (", format(today(), "%d/%m/%y"), ")\n",
                          "Total Cases = ", round(fitted_projected[fitted_projected$Date == today(), "total_infected"]),
                          "\nNew Cases = ", round(new_today)), angle = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sirplot2_log
# uncomment ggsave line to save plot in png format, make sure folder "plots_srn" is created
# ggsave(paste0("plots_srn/sir", name, "_log_asymp.png"), width = 12, height = 9)
