# Get COVID-19 Data from NYT
nyt <- read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
nyt <- subset(nyt, state == "Connecticut")
date <- as.character(nyt[dim(nyt)[1],  1])
write.csv(nyt, file = paste0("../data/CT_cases_deaths_NYT_", date, ".csv"), row.names=FALSE)


# Get ACS 2018 CT data: population and age distribution
library(tidycensus)
library(tidyr)
dat <- get_acs(geography = "county", table = "S0101", state = "CT")
dat <- data.frame(dat)
var <- c(paste0("S0101_C01_00", 1:9), paste0("S0101_C01_0", 10:19))
label <- c("population", paste(seq(0, 85, by=5), seq(0, 85, by=5)+4, sep="-"))
label[19] <-  "85 and over"
dat$label <- NA
for(i in 1:length(var)){
	dat$label[dat$variable == var[i]] <- label[i]
}
demog <- subset(dat, !is.na(label))
demog$county <- gsub(" County, Connecticut", "", demog$NAME)
demog <- demog[, c("county", "label",  "estimate")]
demog$label <-  factor(demog$label, levels = label)
demog.wide <- spread(demog, label, estimate)
write.csv(demog.wide, file = "../data/CT_demog_ACS2018_5year_estimates.csv", row.names=FALSE)
