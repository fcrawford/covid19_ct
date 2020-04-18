library(rgdal)
library(ggplot2)
library(sp)
library(SUMMER)

# Load shapefile
CTmap <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp")

# create and save adj matrix
adj <- getAmat(CTmap, CTmap$NAME10)
write.csv(adj, file = "../map/CT_adj_matrix.csv", quote=FALSE, row.names=TRUE)

# Quick plot
plot(CTmap)

# Plot with names
mapPlot(data=NULL, geo=CTmap, by.geo="NAME10")


