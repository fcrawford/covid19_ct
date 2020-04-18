##
##	Read and organize data from COVID Data Extract (4-17-20)
##	Write ct_hosp_cap.csv  (long format data frame)

library(data.table)
library(ggplot2)
library(reshape2)

raw <- fread("../data/COVID Data Extract (4-17-20).csv")
dat <- transpose(raw[,-1])
colnames(dat) <- as.character(as.matrix(raw)[, 1])
dat <- cbind(Date = as.Date(gsub("X", "",  colnames(raw)[-1]),  format = "%m/%d/%y"), dat)
dat.long <-  melt(dat, id.vars = c("Date", "County"))
dat.long$value <- as.numeric(dat.long$value)
head(dat.long)
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
ggplot(subset(dat.long, variable %in% subset[1:3]), aes(x = Date, y = value, color =  variable)) + geom_line() + facet_wrap(~County, ncol =  4, scale='free') 
ggplot(subset(dat.long, variable =="Total available beds"), aes(x = Date, y = value, color =  variable)) + geom_line() + facet_wrap(~County, ncol =  4, scale='free') 
dat.long$value[is.na(dat.long$value)] <- 0
write.csv(subset(dat.long, variable == "Total available beds"), 
		  "../data/ct_hosp_cap.csv", row.names=FALSE, quote=FALSE)