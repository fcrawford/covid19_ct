library(rgdal)
library(ggplot2)
library(sp)
library(SUMMER)
library(reshape2)
library(patchwork)

# Load shapefile
CTmap <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp")

# create and save adj matrix
adj <- getAmat(CTmap, CTmap$NAME10)
write.csv(adj, file = "../map/CT_adj_matrix.csv", quote=FALSE, row.names=TRUE)

# Quick plot
plot(CTmap)

# Plot with names
mapPlot(data=NULL, geo=CTmap, by.geo="NAME10")

# Make map and adj side by side
g1 <- mapPlot(data=NULL, geo=CTmap, by.geo="NAME10") + theme(panel.border = element_blank()) +  theme(plot.margin = margin(-0.5,-0.5,0,-0.5, "cm"))

adj.long <- cbind(name=rownames(adj), data.frame(adj, check.names=FALSE))
adj.long <- melt(adj.long, id.vars = "name")
adj.long$value <- factor(adj.long$value)
names.plot <- sort(colnames(adj))
adj.long$name <- factor(adj.long$name, levels=names.plot)
adj.long$variable <- factor(adj.long$variable, levels=rev(names.plot))
g2 <- ggplot(adj.long, aes(x=name, y=variable, fill= value)) + geom_tile(color="black") + scale_fill_manual(values=c("1"="gray30", "0"="white")) + scale_x_discrete(position = "top") + theme_minimal()  + theme(axis.ticks = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 30, hjust=0.1), legend.position = "none") + xlab("") + ylab("") +  theme(plot.margin = margin(-0.5,0.5,0,0, "cm"))

 pdf("../../covid19_intvx_report/figures/map_adj.pdf", width=6, height=3)
 print(g1 + g2+  plot_layout(widths = c(1.5,1)))
 dev.off()



