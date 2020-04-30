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


if(FALSE){
	library(ggplot2)
	library(reshape2)
	demog <- read.csv("../data/CT_demog_ACS2018_5year_estimates.csv", check.names=FALSE)
	# demog[, -1] <- demog[, -1] / demog[,]
	demoglong <- melt(demog, id.vars=c("county"))
	ggplot(subset(demoglong, variable != "population"), aes(x=variable, y=value, color=county, group=county)) + geom_point() + geom_line() + scale_color_brewer(palette="Set1") + xlab("Age") + ylab("Population")
}