get_ct_data <- function(day0 = ymd("2020-03-01")){

nyt <- "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"

data0 <- read.csv("../data/us-counties.csv")
data0$date <- as.Date(data0$date)
out_date <- as.numeric(Sys.Date()  - max(data0$date))
if(out_date > 7){
  warning("Local data file out dated by 7 days or more: download from NYT github repository again...", immediate. = TRUE)
  data <- try(read_csv(nyt), TRUE)
  if(is(data, "try-error")){
    # resolve to local version if NYT site is down
    data <- data0
    warning("NYT server not available, using local data.", immediate. = TRUE)
  }else{
    write.csv(data, file = "../data/us-counties.csv")
  }
}else{
  data <- data0
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



## add state-level hospitalizations and deaths: this csv file needs to be updated from CHA reports ##
ct.hosp <- read.csv('../data/ct_hosp.csv')
ct.hosp$date <- ymd(ct.hosp$date)
ct.hosp$time <- round(as.numeric(difftime(ct.hosp$date, day0, units="days")),0)


# merge cases and death counts from NYT with hospitalization data
# if dataset contain different number of observations, use all hospitalizations data
dat_ct_state$date = NULL
   if ( max(dat_ct_state$time) > max(ct.hosp$time) ) 
      {dat_ct_state <- merge(ct.hosp, dat_ct_state, by='time', all=F)} else 
      {dat_ct_state <- merge(ct.hosp, dat_ct_state, by='time', all=T)}
dat_ct_state$total_deaths <- dat_ct_state$deaths
dat_ct_state$deaths <- dat_ct_state$hosp_death
dat_ct_state = subset(dat_ct_state, select = c('time', 'date', 'cases', 'total_deaths', 'deaths', 'hosp_death', 'cur_hosp', 'cum_hosp'))





# estimated hospitalizations and deaths coming from congregate settings
hosp_cong = read.csv("../data/estimate_congregate_hospitalizations.csv", stringsAsFactors = FALSE)
hosp_cong$date = ymd(hosp_cong$date)


# last data point in hosp_cong should be the same as dat_ct_state
# time_range <- seq(min(dat_ct_state$time), max(dat_ct_state$time), by = 1) 
# missing_times = time_range[!time_range %in% hosp_cong$time] 
if (max(dat_ct_state$time) != max(hosp_cong$time)) {stop("estimate_congregate_hospitalizations.csv should have the same last data point as DAT_CT_STATE")}










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



# adjacency matrix for counties
# CTmap <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp", verbose=FALSE)
adj = read.csv("../map/CT_adj_matrix.csv", stringsAsFactors=FALSE)
colnames(adj)[-1] <- rownames(adj) <- adj$X
adj = adj[,-1]
adj = as.matrix(adj)



# hospital capacity at the county level
dat_ct_capacity <- read.csv("../data/ct_hosp_cap.csv")
# set county level capacities 
county_capacities = list()
for(nm in colnames(adj)) {
  county_cap = filter(dat_ct_capacity, County==nm)
  county_days = as.numeric(ymd(county_cap$Date) - day0)
  county_cap_fun = approxfun(county_days, county_cap$Capacity, rule=2)
  county_capacities[[nm]] = county_cap_fun
}










## get smooth hospital death hazard ##

# death function is based on estimated community hCFR, which is computed separately and saved in a CSV file
file.dhaz <- "../data/ct_hdeath_haz.csv"
smooth_hdeath_haz <- read.csv(file.dhaz)
smooth_hdeath_haz$date = ymd(smooth_hdeath_haz$date)


# code for estimation of instantaneous hospital death hazard: deaths[t]/hosp[t-1]
if (FALSE){
### ALL deaths and hospitalizations (community and congregate) ###

if (length(which(is.na(dat_ct_state$cur_hosp) ) ) > 0){
   d = dat_ct_state[1:( min(which(is.na(dat_ct_state$cur_hosp))) -1 ), ]
} else {d = dat_ct_state}

### Community deaths and hospitalizations ###
if (FALSE){
dd = hosp_cong
dd$date = NULL
d = merge(d, dd, by='time')
d$hosp_death = d$deaths - d$hosp_death_cong
d$cur_hosp_all = d$cur_hosp
d$cur_hosp = d$cur_hosp_all - d$cur_hosp_cong
}
###



# daily deaths without smoothing
d$daily_hdeath = c(diff(c(0, d$hosp_death)))

# smooth cumulative deaths
sp <- smooth.spline(d$time, d$hosp_death, nknots=round(nrow(d)/7))
d$smooth.cum_hdeath <- sp$y
# ggplot(d, aes(x = time, y = hosp_death)) + geom_line(alpha = 0.5) + geom_line(aes(y = smooth.cum_hdeath), color='red')+ theme_bw()

# smooth daily deaths: first order difference of smooth cumulative deaths
d$smooth.daily_hdeath <- c(diff(c(0, d$smooth.cum_hdeath)))

# remove early observations with small counts
d = subset(d, time > 26) # total deaths: 26; community deaths 30

# compute hospital death hazard: smooth daily deaths on day k / current hospitalizations on day (k-1)
d$haz = NA
for (k in 2:nrow(d)){
   d$haz[k] = d$smooth.daily_hdeath[k]/d$cur_hosp[k-1]
}
d$haz[1] = d$haz[2]

# smooth hospital death hazard: spline on death hazard
sp <- smooth.spline(d$time, d$haz, nknots=round(nrow(d)/60)) # 45 for community
d$smooth.haz = sp$y
#ggplot(d, aes(x = time, y = haz)) + geom_line(alpha = 0.5) + geom_line(aes(y = smooth.haz), color='red')+ theme_bw()

# compute relative hospital death hazard: relative to the average of first 14 days
d$rel_haz = d$haz / mean(d$haz[1:14]) # 1 for community; 7, 14
d$smooth.rel_haz = d$smooth.haz / mean(d$smooth.haz[1:14]) # 1 for community; 7, 14

#ggplot(d, aes(x = time, y = rel_haz)) + geom_line(alpha = 0.5) + geom_line(aes(y = smooth.rel_haz), color='red')+ theme_bw()


# if smoothing at the end is not reasonable, carry forward the last reasonable value (> 0.15)
if(FALSE){
    for (i in nrow(d):(nrow(d)-15) ){
       if (d$smooth.rel_haz[i] < 0.15) {d$smooth.rel_haz[i]=NA}
    }
    
    ii = which(is.na(d$smooth.rel_haz))
    
    if (length(ii)>0){
       d$smooth.rel_haz[is.na(d$smooth.rel_haz)] = d$smooth.rel_haz[(min(ii)-1)]
    }
}

# subset and save
smooth_hdeath_haz = subset(d, select=c(time, date, daily_hdeath, smooth.daily_hdeath, haz, smooth.haz, rel_haz, smooth.rel_haz))
#ggplot(smooth_hdeath_haz, aes(x = time, y = rel_haz)) + geom_line(alpha = 0.5) + geom_line(aes(y = smooth.rel_haz), color='red')+ theme_bw()
}













## get smooth mobility (contact) data ## 
########################################

## updated with close contact data ##

file.mobi <- "../data/ct_mobility.csv"
data.mobi <- read.csv(file.mobi)

data.mobi$date <- ymd(data.mobi$date)
data.mobi <- data.mobi[order(data.mobi$date), ]

date_range <- seq(min(data.mobi$date), max(data.mobi$date), by = 1) 
missing_dates = date_range[!date_range %in% data.mobi$date] 

if (length(missing_dates)>0){
   for (k in 1:length(missing_dates)){
      data.mobi[nrow(data.mobi) + 1,] = list(ymd(missing_dates[k]), NA)
   }
   
   data.mobi <- data.mobi[order(data.mobi$date), ]
   data.mobi = na_interpolation(data.mobi)
}

data.mobi$time <- round(as.numeric(difftime(ymd(data.mobi$date), day0, units="days")),0)

bl_mobility = mean(data.mobi$mcont[1:41]) # major drop happened in 03-13 (from 1 to 0.5), all time points before this date contribute to 'baseline'
data.mobi$rel_mob = data.mobi$mcont/bl_mobility

# 7-day moving average
xx <- forecast::ma(data.mobi$rel_mob, order = 7)
xx[1:3] = xx[4]
xx[(length(xx)-2):length(xx)] = xx[(length(xx)-3)]
data.mobi$movavg = xx
# ggplot(mb, aes(x = time, y = rel_mob)) + geom_line(alpha = 0.5) + geom_line(aes(y = movavg) , color='red') + theme_bw()

# smooth spline of moving average
#sp <- smooth.spline(data.mobi$time, data.mobi$rel_mob, nknots=round(nrow(data.mobi)/14))
sp <- smooth.spline(data.mobi$time, data.mobi$movavg, nknots=round(nrow(data.mobi)/21)) #14
data.mobi$smooth = sp$y
#data.mobi$smooth[1:41] = 1

# ggplot(data.mobi, aes(x = time, y = rel_mob)) + geom_line(alpha = 0.5) + geom_line(aes(y = smooth) , color='red') + theme_bw()

mob_state = data.mobi[, c("time", "date", "rel_mob", "movavg", "smooth")]



















## get smooth testing volume data ## 
####################################
file.test <- "../data/ct_cum_pcr_tests.csv"
data.test <- read.csv(file.test)
data.test$date = ymd(data.test$date)

# assume linear increase in daily number of tests between March 8 and the first day of test numbers reporting
#add.dates = seq(ymd("2020-03-08"), data.test$date[1], by="day")
#b = data.test$cum_tests[1]/(0.5*length(add.dates)*(length(add.dates)+1))
#add_daily_tests = c(1:length(add.dates))*b
#add_cum_tests = cumsum(add_daily_tests)

#add.data.tests = tibble(add.dates, add_cum_tests)
#colnames(add.data.tests) = colnames(data.test)

#data.test = rbind(add.data.tests, data.test[-1,])

# add time
data.test$t <- as.numeric(data.test$date - data.test$date[1]) + 1

# compute daily number of tests
data.test$daily_tests = c(NA, diff(data.test$cum_tests))
data.test$daily_tests[1] = data.test$daily_tests[2]


# smooth cumulative number of tests
sp <- smooth.spline(data.test$t, data.test$cum_tests, nknots=round(nrow(data.test)/7))
data.test$cum_tests.smooth = sp$y

# ggplot(data.test, aes(x = t, y = cum_tests)) + geom_line(alpha = 0.5) + geom_line(aes(y = cum_tests.smooth), color='red')+ theme_bw()

# smooth daily tests
data.test$daily_tests.cumsmooth = c(NA, diff(data.test$cum_tests.smooth))
data.test$daily_tests.cumsmooth[1] = data.test$daily_tests.cumsmooth[2]

# ggplot(data.test, aes(x = t, y = daily_tests)) + geom_line(alpha = 0.5) + geom_line(aes(y = daily_tests.cumsmooth), color='red')+ theme_bw()

sp <- smooth.spline(data.test$t, data.test$daily_tests.cumsmooth, nknots=round(nrow(data.test)/45))
data.test$daily_tests.smooth <- sp$y
data.test$daily_tests.smooth[data.test$daily_tests.smooth<0]=1

#ggplot(data.test, aes(x = t, y = daily_tests)) + geom_line(alpha = 0.5) + geom_line(aes(y = daily_tests.smooth), color='red')+ theme_bw()
#ggplot(data.test, aes(x = t, y = log(daily_tests))) + geom_line(alpha = 0.5) + geom_line(aes(y = log(daily_tests.smooth)), color='red')+ theme_bw()

data.test$smooth = data.test$daily_tests.smooth
data.test$time <- round(as.numeric(difftime(ymd(data.test$date), day0, units="days")),0)

testing_state = data.test[, c("time", "date", "daily_tests", "smooth")]

# spline
#sp <- smooth.spline(data.test$t, log(data.test$cum_tests), nknots=round(max(data.test$t)/15))
#data.test$smooth.cum <- exp(sp$y)

# 7-day moving average smoothed by spline with biweekly knots
# cum_tests = c(c(0,0,0), data.test$cum_tests)
# sp <- forecast::ma(cum_tests, order = 7)
# 
# sp = sp[4:length(sp)]
# last_daily = sp[length(sp)-3] - sp[length(sp)-4]
# sp[(length(sp)-2):length(sp)] = sp[(length(sp)-3)] + last_daily*c(1,2,3)
# 
# data.test$movavg <- sp
# data.test$daily_movavg = c(diff(c(0, data.test$movavg)))
# data.test$daily_movavg[1] = data.test$daily_tests[1]
# 
# # smooth moving average with spline
# sp <- smooth.spline(data.test$t, log(data.test$movavg), nknots=round(nrow(data.test)/21))
# data.test$smooth.cum <- exp(sp$y)
# 
# data.test$smooth <- c(diff(c(0, data.test$smooth.cum)))
# data.test$smooth[1:3] <- data.test$daily_tests[1:3]
# 
# 
# data.test$time <- round(as.numeric(difftime(ymd(data.test$date), day0, units="days")),0)
# 
# testing_state = data.test[, c("time", "date", "daily_tests", "daily_movavg", "smooth")]
# 







## get smooth severity data ##
##############################
# changes in proportion of cases > 60 y.o after 05/01 #
file.severity = "../data/ct_severity.csv"
data.severity = read.csv(file.severity)
data.severity$date = ymd(data.severity$date)








## get combined testing positive proportion and CLI ED visits data
##############################
# this is a combined measure of incidence
file.incidence = "../data/ct_incidence.csv"
data.incidence = read.csv(file.incidence)
data.incidence$date = ymd(data.incidence$date)









## get smooth hospital LOS ##
#############################
file.hlos = "../data/ct_hlos.csv"
data.hlos = read.csv(file.hlos)
data.hlos$date = ymd(data.hlos$date)









## get vaccinated counts  ##
############################
file.vac = "../data/ct_vaccination.csv"
data.vac = read.csv(file.vac)
data.vac$date = mdy(data.vac$date)










# county populations
pop = read.csv('../data/ct_population.csv', stringsAsFactors=FALSE)
populations = pop$population[match(colnames(adj), pop$county)]
populations = round(populations * (1-0.018)) # 1.8% is an estimated proportion residing in congregate settings (about 65,000)



# state hospitalizations and deaths as proportion of the population
# dat_ct_state$total_deaths.prop = dat_ct_state$total_deaths/sum(populations)
# dat_ct_state$hosp_deaths.prop = dat_ct_state$hosp_death/sum(populations)
# dat_ct_state$cur_hosp.prop = dat_ct_state$cur_hosp/sum(populations)
# dat_ct_state$cum_hosp.prop = dat_ct_state$cum_hosp/sum(populations)


return(list(dat_ct_state=dat_ct_state, dat_ct_county=dat_ct_county, hosp_cong=hosp_cong, incidence = data.incidence,
            adj=adj, dat_ct_capacity=dat_ct_capacity, county_capacities=county_capacities, 
            mob=mob_state, testing = testing_state, severity = data.severity, hlos = data.hlos,
            smooth_hdeath_haz = smooth_hdeath_haz, vac = data.vac, populations=populations))
}

