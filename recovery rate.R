# data
my_data = read.csv("https://wnarifin.github.io/covid-19-malaysia/covid-19_my_full.csv")
my_data$date = as.Date(my_data$date)

# before MCO
my_data1 = subset(my_data, date >= "2020-02-28" & date < "2020-03-18")
my_data1$active = my_data1$total_cases - my_data1$total_recover - my_data1$total_deaths
# shorter dataset
data_my1 = as.data.frame(with(my_data1, cbind(total_cases, recover, total_recover, new_deaths, total_deaths, active)))
str(data_my1)
data_my1$recover_death = data_my1$recover + data_my1$new_deaths
data_my1$days = 1:dim(data_my1)[1]  # generate days
data_my1
# by glm
model1 = glm(recover_death ~ days, data = data_my1, family = "poisson", offset = log(active))
summary(model1)
gamma = exp(coef(model1)[[2]]) - 1; gamma
recover_days = 1/gamma; recover_days

# MCO 1
my_data2 = subset(my_data, date >= "2020-03-18" & date < "2020-04-01")
my_data2$active = my_data2$total_cases - my_data2$total_recover - my_data2$total_deaths
# shorter dataset
data_my2 = as.data.frame(with(my_data2, cbind(total_cases, recover, total_recover, new_deaths, total_deaths, active)))
str(data_my2)
data_my2$recover_death = data_my2$recover + data_my2$new_deaths
data_my2$days = 1:dim(data_my2)[1]  # generate days
data_my2
# by glm
model2 = glm(recover_death ~ days, data = data_my2, family = "poisson", offset = log(active))
summary(model2)
gamma = exp(coef(model2)[[2]]) - 1; gamma
recover_days = 1/gamma; recover_days

# MCO 2
my_data3 = subset(my_data, date >= "2020-04-01")
my_data3$active = my_data3$total_cases - my_data3$total_recover - my_data3$total_deaths
# shorter dataset
data_my3 = as.data.frame(with(my_data3, cbind(total_cases, recover, total_recover, new_deaths, total_deaths, active)))
str(data_my3)
data_my3$recover_death = data_my3$recover + data_my3$new_deaths
data_my3$days = 1:dim(data_my3)[1]  # generate days
data_my3
# by glm
model3 = glm(recover_death ~ days, data = data_my3, family = "poisson", offset = log(active))
summary(model3)
gamma = exp(coef(model3)[[2]]) - 1; gamma
recover_days = 1/gamma; recover_days