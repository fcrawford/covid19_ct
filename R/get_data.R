library(lubridate)



#####################################
# CT Map 

CTmap <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp")


######### Get data on reported counts: updated daily #########


# Can we replace this with CT hospital association data? 


# Get COVID-19 Data from NYT #

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

first.data.time <- round(as.numeric(difftime(dat_ct_state$date[1], day0, units="days")),0)

dat_ct_state$time <- c(first.data.time:(nrow(dat_ct_state)+first.data.time-1)) # add time variable that indexes time 

#time.date <- subset(nyt_ct, select=c(date,time)) # for future matching use
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
first.hosp.data.time <- round(as.numeric(difftime(ct.hosp$date[1], day0, units="days")),0)

ct.hosp$time <- c(first.hosp.data.time:(nrow(ct.hosp)+first.hosp.data.time-1)) 

dat_ct_state <- merge(dat_ct_state, ct.hosp, by='time', all=T)
dat_ct_state$date.y <- NULL
names(dat_ct_state)[names(dat_ct_state) == "date.x"] <- "date"

return(list(dat_ct_state=dat_ct_state, dat_ct_county=dat_ct_county))
}


######### Get data on hospital capacity: updated daily #########
dat_ct_capacity <- read.csv("../data/ct_hosp_cap.csv")
