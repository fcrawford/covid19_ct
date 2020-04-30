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

pop = read.csv('../data/ct_population.csv', stringsAsFactors=FALSE)
populations = pop$population[match(colnames(adj), pop$county)]


return(list(dat_ct_state=dat_ct_state, dat_ct_county=dat_ct_county, CTmap=CTmap, adj=adj, dat_ct_capacity=dat_ct_capacity, county_capacities=county_capacities, populations=populations))
}

