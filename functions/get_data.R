get_ct_data <- function(day0 = ymd("2020-03-01")){

nyt <- "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"
data <- try(read.csv(nyt), TRUE)
if(is(data, "try-error")){
	# resolve to local version if NYT site is down
	data <- read.csv("../data/us-counties.csv")
}

dat_ct_county <- subset(data, state == "Connecticut")
dat_ct_county$date <- ymd(dat_ct_county$date)

## caclulate state-wide observations ##
dat2 <- subset(dat_ct_county, select=c(date, cases, deaths))
dat_ct_state <- aggregate(. ~ date, data=dat2, FUN=function(x){sum(x)})
dat_ct_state <- dat_ct_state[order(dat_ct_state$date),]
dat_ct_state$time <- round(as.numeric(difftime(dat_ct_state$date, day0, units="days")),0) # add time variable that indexes time 

# make corrections for May 2 and May 3 data to match CT DPH reports
dat_ct_state$cases[dat_ct_state$date == "2020-05-02"] <- 29287
dat_ct_state$cases[dat_ct_state$date == "2020-05-03"] <- 29287

dat_ct_state$deaths[dat_ct_state$date == "2020-05-02"] <- 2436
dat_ct_state$deaths[dat_ct_state$date == "2020-05-03"] <- 2495


# calculate daily new reported cases and deaths 
#nyt_ct$new_cases <- nyt_ct$new_deaths <- 0
#nyt_ct$new_cases[1] <- nyt_ct$cases[1] 
#for (i in 2:nrow(nyt_ct)){
#   nyt_ct$new_cases[i] <- nyt_ct$cases[i] - nyt_ct$cases[i-1]
#   nyt_ct$new_deaths[i] <- nyt_ct$deaths[i] - nyt_ct$deaths[i-1]
# }

## add state-level hospitalizations (current counts): this csv file needs to be updated manually ##
ct.hosp <- read.csv('../data/ct_hosp.csv')
ct.hosp$date <- mdy(ct.hosp$date)
ct.hosp$time <- round(as.numeric(difftime(ct.hosp$date, day0, units="days")),0)

# merge cases and death counts from NYT with hospitalization data
dat_ct_state <- merge(dat_ct_state, ct.hosp, by='time', all=T)
dat_ct_state$date.y <- NULL
names(dat_ct_state)[names(dat_ct_state) == "date.x"] <- "date"

## add county-level hospitalization (current counts): this csv file needs to be updated using get_data_capacity.R file ##
ct.hosp.county <- read.csv('../data/ct_current_hosp.csv')
ct.hosp.county$time <- round(as.numeric(difftime(ymd(ct.hosp.county$Date), day0, units="days")),0)
ct.hosp.county$cur_hosp <- ct.hosp.county$value
ct.hosp.county$county <- ct.hosp.county$County
dat_ct_county$time <- round(as.numeric(difftime(ymd(dat_ct_county$date), day0, units="days")),0)
dat_ct_county <- merge(dat_ct_county, ct.hosp.county[, c("time", "county", "cur_hosp")], by=c("time","county"), all=T)
dat_ct_county$state <-  dat_ct_county$state[1]
dat_ct_county$date <- day0 +  dat_ct_county$time
dat_ct_county$cur_hosp[is.na(dat_ct_county$cur_hosp)] <- 0
# make sure when hosp data is lagged behind, we use NA
hosp_last_day <- max(ct.hosp.county$time)
dat_ct_county$cur_hosp[dat_ct_county$time > hosp_last_day] <- NA

dat_ct_county$cases[is.na(dat_ct_county$cases)] <- dat_ct_county$cur_hosp[is.na(dat_ct_county$cases)]

# make corrections for May 2 and May 3 data to match CT DPH reports
dat_ct_county$cases[dat_ct_county$date == "2020-05-02" & dat_ct_county$county == "Unknown" ] <- 312
dat_ct_county$deaths[dat_ct_county$date == "2020-05-02" & dat_ct_county$county == "Unknown" ] <- 1

dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Fairfield" ] <- 886
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Hartford" ] <- 756
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Litchfield" ] <- 92
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Middlesex" ] <- 93
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "New Haven" ] <- 580
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "New London" ] <- 43
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Tolland" ] <- 40
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Windham" ] <- 3
dat_ct_county$deaths[dat_ct_county$date == "2020-05-03" & dat_ct_county$county == "Unknown" ] <- 2


CTmap <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp", verbose=FALSE)
adj = read.csv("../map/CT_adj_matrix.csv", stringsAsFactors=FALSE)
colnames(adj)[-1] <- rownames(adj) <- adj$X
adj = adj[,-1]
adj = as.matrix(adj)


dat_ct_capacity <- read.csv("../data/ct_hosp_cap.csv")
# set county level capacities 
county_capacities = list()
for(nm in colnames(adj)) {
  county_cap = filter(dat_ct_capacity, County==nm)
  county_days = as.numeric(ymd(county_cap$Date) - day0)
  county_cap_fun = approxfun(county_days, county_cap$Capacity, rule=2)
  county_capacities[[nm]] = county_cap_fun
}



## get smooth mobility data ## 
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
sp <- smooth.spline(data.mobi$t, data.mobi$relative_mobility, nknots=round(max(data.mobi$t)/15))
data.mobi$smooth <- sp$y
mob_state = data.mobi[, c("date", "smooth")]






## get smooth testing data ## 
# assume linear increase in daily number of tests between day0 (March 1) and the first day of test numbers reporting
file.test <- "../data/ct_cum_pcr_tests.csv"
data.test <- read.csv(file.test)
data.test$date = ymd(data.test$date)
if(sum(diff(data.test$date) != 1) > 0) stop("Missing dates")

add.dates = seq(day0, data.test$date[1], by="day")
b = data.test$cum_tests[1]/(0.5*length(add.dates)*(length(add.dates)+1))
add_daily_tests = c(1:length(add.dates))*b
add_cum_tests = cumsum(add_daily_tests)

add.data.tests = tibble(add.dates, add_cum_tests)
colnames(add.data.tests) = colnames(data.test)

data.test = rbind(add.data.tests, data.test[-1,])
#data.test$daily_tests = c(diff(c(0, data.test$cum_tests)))

data.test$t <- as.numeric(data.test$date - data.test$date[1]) + 1
sp <- smooth.spline(data.test$t, log(data.test$cum_tests), nknots=round(max(data.test$t)/15))
data.test$smooth.cum <- exp(sp$y)
data.test$smooth <- c(diff(c(0, data.test$smooth.cum)))
data.test$smooth[1] <- data.test$smooth[2]

testing_state = data.test[, c("date", "smooth")]

# county populations
pop = read.csv('../data/ct_population.csv', stringsAsFactors=FALSE)
populations = pop$population[match(colnames(adj), pop$county)]

return(list(dat_ct_state=dat_ct_state, dat_ct_county=dat_ct_county, CTmap=CTmap, adj=adj, dat_ct_capacity=dat_ct_capacity, county_capacities=county_capacities, mob=mob_state, testing = testing_state, populations=populations))
}

