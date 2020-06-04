library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(patchwork)
w <- read.csv("../data/dat_weights.csv")
w$date <- as.POSIXct(w$date, format="%m/%d/%y")
color_dat <- data.frame(time = seq(min(w$date), max(w$date), by = "1 hour"))
color_dat$date <- as.character(format(color_dat$time, "%Y-%m-%d"))
color_dat$dobs_wt <- w$dobs_wt[match(color_dat$date, as.character(w$date))]
color_dat$hobs_wt <- w$hobs_wt[match(color_dat$date, as.character(w$date))]

mypalette <- rev(colorRampPalette(brewer.pal(8, "Greys"))(100)[1:70]) 
w <- rbind(cbind(w, obs = w$deaths, type="Cumulative deaths"), 
		   cbind(w, obs = w$cur_hosp, type="Hospitalizations"))
g1 <- ggplot(data=w, aes(x=date, y=obs)) + 
			    geom_segment(data = color_dat, aes(x = time, xend = time,
                             y = -Inf, yend = Inf, color = hobs_wt)) + 
				scale_color_gradientn("Weights", colors=mypalette, lim = c(0, 1)) + 
				geom_point(aes(fill=type), shape=21, stroke=0, size=2) + 
				scale_fill_manual("", values = c("#e41a1c", "#377eb8")) + 
				scale_x_datetime(expand=c(0,0), limits = range(w$date)) + 
				xlab("") + ylab("People") +
			    theme_bw() + 
			    theme(panel.border = element_blank(), legend.title=element_text(size=9), 
			    	 axis.line = element_line(colour = "gray30"))
g2 <- ggplot(data=w, aes(x = date, y = hobs_wt)) + geom_line() + 
 				xlab("") + ylab("Observation weight") + theme_bw() + 
			    theme(panel.border = element_blank(), legend.title=element_text(size=9), 
			    	 panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
			    	 axis.line = element_line(colour = "gray30"))

ggsave(g2 + g1, file = "weights.pdf", width=11*.9, height = 4*.9)