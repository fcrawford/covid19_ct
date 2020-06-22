##
## Script to smooth mobility and testing data
##
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(forecast)

file.mobi <- "../data/ct_mobility.csv"
data.mobi <- read.csv(file.mobi)
data.mobi$county <- data.mobi$polygon_name
pop <- read.csv("../data/ct_population.csv")
pop$population <- pop$population / sum(pop$population)
data.mobi <- data.mobi %>% 
			left_join(pop[, -1]) %>%  
			group_by(ds) %>% 
			summarize(stay_put = sum(all_day_ratio_single_tile_users * population))
data.mobi$date <- as.Date(data.mobi$ds)
data.mobi <- data.mobi[order(data.mobi$date), ]
bl_mobile_prop = 1 - mean(data.mobi$stay_put[1:7])
data.mobi$mobile_prop = 1 - data.mobi$stay_put
data.mobi$relative_mobility = 1 - (bl_mobile_prop - data.mobi$mobile_prop)/bl_mobile_prop
data.mobi$t <- as.numeric(data.mobi$date - data.mobi$date[1]) + 1
# sp <- smooth.spline(data.mobi$t, data.mobi$relative_mobility, nknots=round(max(data.mobi$t)/15))
# data.mobi$smooth <- sp$y
sp <- forecast::ma(data.mobi$relative_mobility, order = 7)
data.mobi$smooth <- sp
g0 <- ggplot(data.mobi, aes(x = date, y = relative_mobility)) + geom_line(alpha = 0.5) + geom_line(aes(y = smooth), color='red')+ theme_bw()



file.test <- "../data/ct_cum_pcr_tests.csv"
data.test <- read.csv(file.test)
data.test$date <- as.Date(data.test$date, format = "%Y-%m-%d")
data.test$daily_tests <- c(NA, diff(data.test$cum_tests))
data.test$t <- as.numeric(data.test$date - data.test$date[1]) + 1
# sp <- smooth.spline(data.test$t, (data.test$cum_tests), nknots=round(max(data.test$t)/15))
# data.test$smooth.cum <- (sp$y)
sp <- forecast::ma(data.test$cum_tests, order = 7)
data.test$smooth.cum <- sp
data.test$smooth <- c(diff(c(0, data.test$smooth.cum)))
data.test$smooth[1] <- data.test$smooth[2]

g1 <- ggplot(data.test, aes(x = date, y = daily_tests)) + geom_line(alpha=0.5)+ geom_line(aes(y = smooth), color='red') + theme_bw()
g2 <- ggplot(data.test, aes(x = date, y = cum_tests)) + geom_line(alpha=0.5)+ geom_line(aes(y = smooth.cum), color='red') + theme_bw()
g0 + g1 + g2
 

