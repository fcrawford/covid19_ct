##
##	Read and organize data from COVID Data Extract (4-17-20)
##	Write ct_hosp_cap.csv  (long format data frame)
##  Write ct_current_hosp.csv 
##  Write ct_cum_hosp.csv 

library(data.table)
library(ggplot2)
library(reshape2)

# data_stream <- "4-17-20"
data_stream <- "5-01-20"
raw <- fread(paste0("../data/COVID Data Extract (", data_stream, ").csv"))
dat <- transpose(raw[,-1])
colnames(dat) <- as.character(as.matrix(raw)[, 1])
dat <- cbind(Date = as.Date(gsub("X", "",  colnames(raw)[-1]),  format = "%m/%d/%y"), dat)
message(which(diff(dat$Date) < 0))

# fix bug in spreadsheet received
if(data_stream %in% c("4-21-20")){
	#if the message above is 40
	print(dat$Date[31:51])
	dat$Date[40:47] <- dat$Date[40]
}
tab <- matrix(c(0, dat$Date-18338), ncol=8, byrow=TRUE)
check1 <- apply(tab, 1, function(x){length(unique(x))}) 
check2 <-apply(tab, 2, function(x){length(unique(x))})

if(sum(!(check1==1)) != 0 || length(unique(check2)) > 1){
	stop("Manual inspection of the spreadsheet is needed")
}else{
dat.long <-  melt(dat, id.vars = c("Date", "County"))
dat.long$value <- as.numeric(dat.long$value)


##
##	Capacity
##
subset <- c("Negative pressure total available beds",
		"Negative pressure beds - Occupied",
		"Negative pressure additional surge beds",
		"ICU/critical care total available beds",
		"ICU/critical care beds - Occupied",
		"ICU/critical care additional surge beds",
		"Isolation total available beds",
		"Isolation beds - Occupied",
		"Isolation additional surge beds",
		"Total available beds",
		"Addition surge beds")
# ggplot(subset(dat.long, variable %in% subset[c(2,5,8,10)]), aes(x = Date, y = value, color =  variable)) + geom_line() + facet_wrap(~County, ncol =  4, scale='free') 
dat.long$value[is.na(dat.long$value)] <- 0

# subset1 <- c("Total available beds",
# 			"Addition surge beds",
# 			"Negative pressure beds - Occupied",
# 			"ICU/critical care beds - Occupied",
# 			"Isolation beds - Occupied",
# 			"Inpatient COVID-positive census"
# 			)
# dat.capacity <- subset(dat.long, variable %in% subset1)
# dat.capacity <- data.frame(dcast(dat.capacity, Date  + County ~  variable,  value.var="value"), check.names=FALSE)
# dat.capacity$Capacity <- dat.capacity[, subset1[1]] + dat.capacity[, subset1[2]] - 
# 						 dat.capacity[, subset1[3]] - dat.capacity[, subset1[4]] - 
# 						 dat.capacity[, subset1[5]] + dat.capacity[, subset1[6]]   

dat.capacity <- subset(dat.long, variable %in% "Total available beds")
dat.capacity$Capacity <- dat.capacity$value
write.csv(dat.capacity[, c("Date", "County", "Capacity")], 
		  "../data/ct_hosp_cap.csv", row.names=FALSE, quote=FALSE)
# ggplot(dat.capacity) + geom_line(aes(x = Date, y = `Total available beds`), color =  "blue")+ geom_line(aes(x = Date, y = `Addition surge beds`), color =  "darkgreen") + geom_line(aes(x = Date, y = Capacity), color =  "red") + facet_wrap(~County, ncol =  4, scale='free') 



##
##	Current hospitalization
##
current <- subset(dat.long, variable == "Inpatient COVID-positive census")
alldates <-  min(current$Date) + c(0:(max(current$Date)-min(current$Date)))
for(i in 1:length(alldates)){
	if(alldates[i] %in% current$Date==FALSE){
		counties <- unique(current$County)
		previous <- subset(current, Date==alldates[i] - 1)
		values <- previous$value[match(counties, previous$County)]
		current <- rbind(current, data.frame(Date=alldates[i], 
											 County = counties,
											 variable="Inpatient COVID-positive census",
											 value = values))
		message("\n", alldates[i], " imputed from previous day")
	}
}
write.csv(current, 
		  "../data/ct_current_hosp.csv", row.names=FALSE, quote=FALSE)



##
##	New hospitalization
##
new <- subset(dat.long, variable == "COVID-associated admissions in last 24 hours")
# 04/08 data is likely mislabeled 04/07
newdate <- as.Date("2020-04-08")
counties <- unique(new$County)
previous <- subset(new, Date==newdate - 1)
values <- previous$value[match(counties, previous$County)]
new <- rbind(new, data.frame(Date=newdate, 
									 County = counties,
									 variable="COVID-associated admissions in last 24 hours",
									 value = values))
new$value[is.na(new$value)] <- 0
# This is from https://portal.ct.gov/-/media/Coronavirus/CTDPHCOVID19summary4162020.pdf?la=en
# Windham is manually altered to be sum of cumulative new admissions instead of 16 in the report
dat_4162020 <- data.frame(County = c("Fairfield","Hartford","Litchfield","Middlesex","New Haven","New London","Tolland","Windham"))
dat_4162020$cum_hosp <- c(1852, 932, 123, 251, 1194, 41, 24, 26)
dat_4162020$Date <- as.Date("2020-04-16")
new <- merge(new, dat_4162020, all=TRUE)
pre <- as.numeric(difftime(dat_4162020$Date[1], as.Date("2020-04-07"), unit='days'))
for(c in dat_4162020$County){	
	for(i in 1:pre){
		row <- which(new$County==c &  new$Date==dat_4162020$Date[1]-i)	
		row.next <- which(new$County==c &  new$Date==dat_4162020$Date[1]-i+1)	
		new$cum_hosp[row] <- new$cum_hosp[row.next] -  new$value[row.next]
	}
}
post <- as.numeric(difftime(max(new$Date), dat_4162020$Date[1], unit='days'))
for(c in dat_4162020$County){	
	for(i in 1:post){
		row <- which(new$County==c &  new$Date==dat_4162020$Date[1]+i)	
		row.prev <- which(new$County==c &  new$Date==dat_4162020$Date[1]+i-1)	
		new$cum_hosp[row] <- new$cum_hosp[row.prev] +  new$value[row]
	}
}
message("\n==== Sample cum hospitalization data ====")
print(tail(subset(data.frame(new)[,-3], County==dat_4162020$County[8])))
write.csv(new[,c("Date", "County", "cum_hosp")], 
		  "../data/ct_cum_hosp.csv", row.names=FALSE, quote=FALSE)

message("\n**** Hospital data updated ****")
}


