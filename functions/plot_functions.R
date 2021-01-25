####################################
# Four generic plotting functions
# 
# 1. plot_ct_region: plot trajectory
# 2. plot_ct_region_list: plot a list of trajectories (sharing same y axis)
# 3. mapplot_ct_region: plot estimates on maps
# 4. mapplot_ct_region: plot a list of estimates on maps (sharing same y axis)
#
####################################


####################
# plotting 
# @param data column include [variable, time, mean, lower, upper]
# @param which.plot the prefix of the compartment to plot, e.g., S, E, I_s, I_m. If a vector of more than one specified, it plots multiple lines
# @param title the tile of the figure
# @param xlab the xlab of the figure
# @param ylab the ylab of the figure
# @param color the color of the line


plot_ct_region = function(data=NULL, 
                          region_name = "Connecticut", 
                          which.plot = "D", 
                          color =  NULL,
                          title=NULL, xlab=NULL, ylab=NULL,
                          #tmax.plot = tmax,
                          #start_day = day0,
                          end_day=NULL, # pass in daymax
                          capacity_func = COUNTY_CAPACITIES,
                          obs_state = DAT_CT_STATE,
                          obs_county = DAT_CT_COUNTY, 
                          ymax = NULL,
                          avg = "mean",
                          sentence=TRUE,
                          show.data=TRUE, 
                          title.override=NULL,
                          goodness = FALSE,
                          add.cong = FALSE,
                          ref_day=Sys.Date()) {

lab.table <- data.frame(compartment=c("D","rD", "H","rH",
                                       "Hbar", "rHbar", "rHsum","rcum_modH",
                                       "S","E","I_s","I_m",
                                       "A", "dailyI", "cum_modI", "alive_cum_incid_prop", 
                                       "alive_cum_incid_num", "R_eff", "currentI", 
                                       "contact_pattern"),
                          color=c('#e41a1c','#e41a1c','#377eb8','#377eb8', 
                                  '#4daf4a','#4daf4a','#377eb8', '#984ea3',
                                  '#ff7f00','#ffff33', '#a65628','#f781bf',
                                  '#999999', '#a65628', '#388a72', "#006d2c", 
                                  "#09543e", "#6a3d9a", "#ff7f00", 
                                  "#a34e7e"),
                          labels=c("Cumulative deaths", "Cumulative deaths",  "Hospitalizations", "Hospitalizations",
                                    "Hospital overflow", "Hospital overflow",  "Required hospitalizations",  "Cumulative hospitalizations",
                                    "Susceptible population", "Exposed population", "Severe infections", "Mild infections",
                                    "Asymptomatic infections",   "Daily new infections", "Cumulative infections",  "Cumulative incidence proportion", 
                                    "Cumulative incidence among living", "R_eff", "Current Infections", 
                                    "Overall contact pattern"))

  if(which.plot %in% lab.table$compartment == FALSE){
    lab.table <- rbind(lab.table, data.frame(compartment=which.plot, color = "#006d2c", labels=which.plot))
  }
  if(!is.null(region_name)){
    variable_name = paste0(which.plot,".",region_name)
    toplot <- paste(rep(which.plot,each=length(region_name)),
                rep(region_name, length(which.plot)),sep=".")
  }else{
    variable_name = which.plot
    region_name = ""
    toplot <- which.plot
  }

  if(goodness){
    if(which.plot %in% c("rD", "rH", "rHsum", "rcum_modH", "Inc")){
      end_day <- ref_day
      # ymax <- NULL
    }else{
      goodness <- FALSE
    }
  }
  
 #dayseq = seq(day0, daymax, by="day")
  start_day = day0
  tmax.plot = as.numeric(difftime(end_day, day0, units="days"))

  if(goodness) {
    monthseq = seq(start_day, ymd(end_day+1), by="week")
    lab_show = format(monthseq, "%b %d")
    lab_where = difftime(monthseq, start_day, units="days")
    lab_cex <- 0.95
  } else {
    monthseq = seq(start_day, ymd(end_day+1), by="month")
    lab_show = format(monthseq, "%b %Y")
    lab_where = difftime(monthseq, start_day, units="days")
    lab_cex <- 0.95
  }
  

  which.plot.ci <- which.plot
  add <- FALSE
  if("rH" %in% which.plot) add <- TRUE
  if("rHsum" %in% which.plot) add <- TRUE

  par(mar=c(2.5,4,2.5,0), bty="n")
  sir_result_region= filter(data, variable%in%toplot)

  comm <- ifelse(add.cong, "", "Community ")
  if(!add.cong){
    title0 <- paste0("Community ", tolower(lab.table$labels[lab.table$compartment==which.plot[1]]))
  }else{
    title0 <- paste0(lab.table$labels[lab.table$compartment==which.plot[1]])
  }

  if(region_name != ""){
    title <- paste0(title0, " in ", region_name, " ", title)
  }else{
    title <- title0
  }

  if(!is.null(title.override)) {
    title = title.override
  }


  tmin.plot <- 0
  if("Inc" %in% which.plot && goodness) tmin.plot <- min(obs_state$time[obs_state$rel_inc > 0])
    

  if(add && (!goodness)){
    if(region_name %in% names(capacity_func)){
        sub.add <- data.frame(time = 0:tmax.plot, 
                              Capacity = capacity_func[[region_name]](0:tmax.plot))
    }else if(region_name == "Connecticut"){
      cap.state <- rep(0, tmax.plot+1)
      for(nm in 1:length(capacity_func)){
          cap.state <- cap.state + capacity_func[[nm]](0:tmax.plot)
      }
      sub.add <- data.frame(time = 0:tmax.plot, Capacity = cap.state)
    }else{
      stop(paste0(region_name, " not recognized in capacity function"))
    }
    # ymax <- max(ymax, sub.add$Capacity)
  }

  # get observations
  if(region_name == "Connecticut") {
    points.add <- obs_state
    # points(obs_state$time, obs_state$deaths, pch=16, cex=0.6, col=col.line) 
  }else if(region_name != "") {
    obs.region <- subset(obs_county, county == region_name)
    obs.region$date <- ymd(obs.region$date)
    first.region.time <- round(as.numeric(difftime(obs.region$date[1], start_day, units="days")),0)
    obs.region$time <- c(first.region.time:(nrow(obs.region)+first.region.time-1)) # add time variable that indexes time 
    points.add <- obs.region
    # points(obs.region$time, obs.region$deaths, pch=16, cex=0.6, col=col.line)
  }



  lab.table$color <- as.character(lab.table$color)
  col.line <- lab.table$color[which(lab.table$compartment==which.plot)]
  col.polygon <- adjustcolor(col.line, alpha.f = 0.5 - 0.1*goodness)
  sir_result_region_sub <- filter(sir_result_region, variable==variable_name)
  sir_result_region_sub <- subset(sir_result_region_sub, time <= tmax.plot)
  
  if(add.cong && which.plot == "rHsum"){
    for(i in 1:dim(sir_result_region_sub)[1]){
      cong <- HOSP_CONG$cur_hosp_cong[which(HOSP_CONG$time == sir_result_region_sub$time[i])]
      if(length(cong) == 0) cong <- 0

      sir_result_region_sub[i, "mean"] <- sir_result_region_sub[i, "mean"] + cong
      sir_result_region_sub[i, "median"] <- sir_result_region_sub[i, "median"] + cong
      sir_result_region_sub[i, "lower"] <- sir_result_region_sub[i, "lower"] + cong
      sir_result_region_sub[i, "upper"] <- sir_result_region_sub[i, "upper"] + cong
    }
    
    ## use projected estimated proportion of community / congregare hospitalizations to apply to model-projected dynamics 
    # predict hospitalizations for projections
      d = merge(DAT_CT_STATE, HOSP_CONG, by = 'time', all=F)
      d$prop_hosp_cong = d$cur_hosp_cong/d$cur_hosp

      if (FALSE){
      # fit regression line to estimate changes in prop_cong_hosp during the last 90 days
      d2 = subset(d, time > (d$time[nrow(d)] - 90) , select = c('time', 'prop_hosp_cong'))
      rg = lm(prop_hosp_cong ~ time, data=d2)
      
      # calculate predicted values
      d3 = merge(sir_result_region_sub, d2, by='time', all=T)
      d3$pred.prop_hosp_cong = predict(rg, newdata=d3)
      d3$cong_mean = d3$mean * d3$pred.prop_hosp_cong/(1-d3$pred.prop_hosp_cong)
      d3$cong_median = d3$median * d3$pred.prop_hosp_cong/(1-d3$pred.prop_hosp_cong)
      # adjust predicted values by the difference beween estimated and predicted values at last observation to make a smooth transition
      adjh_mean = HOSP_CONG$cur_hosp_cong[nrow(HOSP_CONG)] - d3$cong_mean[d3$time == max(HOSP_CONG$time)]
      adjh_median = HOSP_CONG$cur_hosp_cong[nrow(HOSP_CONG)] - d3$cong_median[d3$time == max(HOSP_CONG$time)]
      
      # remove values for observed period
      d3$cong[d3$time <= max(HOSP_CONG$time)] = NA
      }
      
      # use mean / median proportion of congregate hospitalizations during the last 30 days for projections into the future
      pred.prop_hosp_cong = mean(d$prop_hosp_cong[(nrow(d)-30):nrow(d)])
      d3 = sir_result_region_sub
      d3$cong_mean = d3$mean * pred.prop_hosp_cong/(1-pred.prop_hosp_cong)
      d3$cong_median = d3$median * pred.prop_hosp_cong/(1-pred.prop_hosp_cong)
      
      # remove values for observed period
      d3$cong_mean[d3$time <= max(HOSP_CONG$time)] = NA
      d3$cong_median[d3$time <= max(HOSP_CONG$time)] = NA
      
      # set adjustment to 0
      adjh_mean = adjh_median = 0
      
    # add predicted values to projections
    for (i in (max(HOSP_CONG$time)+2):dim(sir_result_region_sub)[1]){
      cong_mean <- d3$cong_mean[i] + adjh_mean
      cong_median <- d3$cong_median[i] + adjh_median
      sir_result_region_sub[i, "mean"] <- sir_result_region_sub[i, "mean"] + cong_mean
      sir_result_region_sub[i, "median"] <- sir_result_region_sub[i, "median"] + cong_median
      sir_result_region_sub[i, "lower"] <- sir_result_region_sub[i, "lower"] + cong_mean
      sir_result_region_sub[i, "upper"] <- sir_result_region_sub[i, "upper"] + cong_mean
    }
  }


  if(add.cong && which.plot == "rD"){
    cong.data <- data.frame(HOSP_CONG[, c("hosp_death_cong", "time")]) 
    tmp <- data.frame(time = DAT_CT_STATE$time, 
                      diff = DAT_CT_STATE$total_deaths - DAT_CT_STATE$hosp_death)
    cong.data <- merge(cong.data, tmp, by = "time", all.x = TRUE)
    cong.data$diff[is.na(cong.data$diff) & cong.data$time<100] <- 0
    cong.data$diff[is.na(cong.data$diff) & cong.data$time>=100] <- max(cong.data$diff, na.rm = TRUE)
    cong.data$cong <- cong.data$hosp_death_cong + cong.data$diff


    for(i in 1:dim(sir_result_region_sub)[1]){
      cong <- cong.data$cong[which(cong.data$time == sir_result_region_sub$time[i])]
      if(length(cong) == 0) cong <- cong.data$cong[which.max(cong.data$time)]

      sir_result_region_sub[i, "mean"] <- sir_result_region_sub[i, "mean"] + cong
      sir_result_region_sub[i, "median"] <- sir_result_region_sub[i, "median"] + cong    
      sir_result_region_sub[i, "lower"] <- sir_result_region_sub[i, "lower"] + cong
      sir_result_region_sub[i, "upper"] <- sir_result_region_sub[i, "upper"] + cong
    }
  

  ## add predicted deaths from congregate hospitalizations to projections ###
  # estimate additional hospitalizations and project into the future  
  hsp = filter(data, variable=="rHsum.Connecticut")
      for(i in 1:dim(hsp)[1]){
      cong <- HOSP_CONG$cur_hosp_cong[which(HOSP_CONG$time == hsp$time[i])]
      if(length(cong) == 0) cong <- 0

      hsp[i, "mean"] <- hsp[i, "mean"] + cong
      hsp[i, "median"] <- hsp[i, "median"] + cong      
      hsp[i, "lower"] <- hsp[i, "lower"] + cong
      hsp[i, "upper"] <- hsp[i, "upper"] + cong
    }
    
    ## use projected estimated proportion of community / congregare hospitalizations to apply to model-projected dynamics 
    # predict congregate hospitalizations for projections
      d = merge(DAT_CT_STATE, HOSP_CONG, by = 'time', all=F)
      d$prop_hosp_cong = d$cur_hosp_cong/d$cur_hosp

      if (FALSE){
      # fit regression line to estimate changes in prop_cong_hosp during the last 90 days
      d2 = subset(d, time > (d$time[nrow(d)] - 90) , select = c('time', 'prop_hosp_cong'))
      rg = lm(prop_hosp_cong ~ time, data=d2)
      
      # calculate predicted values
      d3 = merge(hsp, d2, by='time', all=T)
      d3$pred.prop_hosp_cong = predict(rg, newdata=d3)
      d3$cong = d3$mean * d3$pred.prop_hosp_cong/(1-d3$pred.prop_hosp_cong)
      d3$cong[d3$time <= max(DAT_CT_STATE$time)] = NA
      d3 = subset(d3, !is.na(d3$cong))
      d3$deaths_cong = (1/HLOS$smooth.alos[nrow(HLOS)])*d3$cong*0.38
      d3$cum_deaths_cong = cumsum(d3$deaths_cong)
      }
      
      # use mean proportion of congregate hospitalizations during the last 30 days for projections 
      pred.prop_hosp_cong = mean(d$prop_hosp_cong[(nrow(d)-30):nrow(d)])
      d3 = hsp
      d3$cong_mean = d3$mean * pred.prop_hosp_cong/(1-pred.prop_hosp_cong)
      d3$cong_median = d3$median * pred.prop_hosp_cong/(1-pred.prop_hosp_cong)
      
      # remove values for observed period
      d3$cong_mean[d3$time <= max(HOSP_CONG$time)] = NA
      d3$cong_median[d3$time <= max(HOSP_CONG$time)] = NA
      
      d3 = subset(d3, !is.na(d3$cong_mean))
      
      # daily deaths = estimated daily discharges * current cong hospitalizations * hCFR among congregate residents
      d3$deaths_cong_mean = (1/HLOS$smooth.alos[nrow(HLOS)])*d3$cong_mean*0.38 
      d3$cum_deaths_cong_mean = cumsum(d3$deaths_cong_mean)
      
      d3$deaths_cong_median = (1/HLOS$smooth.alos[nrow(HLOS)])*d3$cong_median*0.38
      d3$cum_deaths_cong_median = cumsum(d3$deaths_cong_median)

  cng = merge(sir_result_region_sub, d3, by='time', all=T)
   
  # add estimated numbers to projections
   for (i in (max(HOSP_CONG$time)+2):dim(sir_result_region_sub)[1]){
     cong_mean = cng$cum_deaths_cong_mean[i]
     cong_median = cng$cum_deaths_cong_median[i]
       sir_result_region_sub[i, "mean"] <- sir_result_region_sub[i, "mean"] + cong_mean
       sir_result_region_sub[i, "median"] <- sir_result_region_sub[i, "median"] + cong_median      
       sir_result_region_sub[i, "lower"] <- sir_result_region_sub[i, "lower"] + cong_mean
       sir_result_region_sub[i, "upper"] <- sir_result_region_sub[i, "upper"] + cong_mean
     }
  }

  
  
  
  if(is.null(ymax)) ymax <- max(sir_result_region_sub$upper[sir_result_region$time <= tmax.plot & sir_result_region$time >= tmin.plot], na.rm=TRUE)

  if("rD" %in% which.plot){
    points.add$toplot <- points.add$total_deaths
    ymax <- max(c(ymax, points.add$toplot), na.rm=TRUE)
  }
  if("rH" %in% which.plot || "rHsum" %in% which.plot){
    points.add$toplot <- points.add$cur_hosp
    ymax <- max(c(ymax, points.add$toplot), na.rm=TRUE)
  }
   if("rcum_modH" %in% which.plot){
    points.add$toplot <- points.add$cum_hosp
    ymax <- max(c(ymax, points.add$toplot), na.rm=TRUE)
  }
   if("Inc" %in% which.plot){
    points.add$toplot <- points.add$rel_inc
    ymax <- max(c(ymax, points.add$toplot[points.add$rel_inc > 0]), na.rm=TRUE)
  }


   if(is.null(xlab)) xlab <- ""
   if(is.null(ylab)) ylab <- "People"
    plot(0, type="n", xlab=xlab, ylab=ylab, main=title, col="black", 
         ylim=c(0,1.05*ymax), xlim=c(tmin.plot,1.1*tmax.plot), axes=FALSE)
    axis(1,at=lab_where, lab=lab_show, cex.axis=lab_cex)
    axis(2)

    if(!goodness) abline(v=ref_day-start_day, col="gray", lty=2)

    time.print <- tmax.plot + 0.5
    polygon(c(sir_result_region_sub$time, rev(sir_result_region_sub$time)), c(sir_result_region_sub$lower, rev(sir_result_region_sub$upper)), col=col.polygon, border=NA)
    if(!goodness){
      if (avg=='median') {text(sir_result_region_sub$time[tmax.plot+1], sir_result_region_sub$median[time.print], format(sir_result_region_sub$median[time.print],digits=2, big.mark=","), pos=4, col=col.line)} else
      {text(sir_result_region_sub$time[tmax.plot+1], sir_result_region_sub$mean[time.print], format(sir_result_region_sub$mean[time.print],digits=2, big.mark=","), pos=4, col=col.line)}
      text(sir_result_region_sub$time[tmax.plot+1], sir_result_region_sub$lower[time.print], format(sir_result_region_sub$lower[time.print],digits=2, big.mark=","), pos=4, col=col.polygon)
      text(sir_result_region_sub$time[tmax.plot+1], sir_result_region_sub$upper[time.print], format(sir_result_region_sub$upper[time.print],digits=2, big.mark=","), pos=4, col=col.polygon)
    }
  
    if (avg=='median') {lines(sir_result_region_sub$time, sir_result_region_sub$median, col=col.line)} else
    {lines(sir_result_region_sub$time, sir_result_region_sub$mean, col=col.line)}
  

  if(add && (!goodness)){
    lines(sub.add$time, sub.add$Capacity, col='gray30', lty  = 2,  lwd=1.2)
  }

  # Add observed deaths
  if(which.plot %in% "rD"){
    # col.line <- lab.table$color[which(lab.table$compartment=="D")]
    points.add <- subset(points.add, time <= ref_day-start_day)
    points(points.add$time, points.add$toplot, pch=16, cex=0.6, col=col.line)
  }
  if(which.plot %in% c("rH", "rHsum")){
    # col.line <- lab.table$color[which(lab.table$compartment=="H")]
    points.add <- subset(points.add, time <= ref_day-start_day)
    points(points.add$time, points.add$toplot, pch=16, cex=0.6, col=col.line)
  }
  if(which.plot %in% "rcum_modH"){
    # col.line <- lab.table$color[which(lab.table$compartment=="rcum_modH")]
    points.add <- subset(points.add, time <= ref_day-start_day)
    points(points.add$time, points.add$toplot, pch=16, cex=0.6, col=col.line)
  }
 if(which.plot %in% "Inc"){
    # col.line <- lab.table$color[which(lab.table$compartment=="currentI")]
    points.add <- subset(points.add, time <= ref_day-start_day)
    points(points.add$time, points.add$toplot, pch=16, cex=0.6, col=col.line)
  }


  if("D" %in% which.plot || "rD" %in% which.plot || "rHsum" %in% which.plot || "dailyI" %in% which.plot || "Inc" %in% which.plot){
      sir_result_region_sub <- filter(sir_result_region, variable==variable_name)
  }

  region_summary <- NULL
  count.min <- count.max <- NA
  if("D" %in% which.plot || "rD" %in% which.plot){
      count <- sir_result_region_sub$mean[sir_result_region_sub$time==tmax.plot]
      count.min <- sir_result_region_sub$lower[sir_result_region_sub$time==tmax.plot]
      count.max <- sir_result_region_sub$upper[sir_result_region_sub$time==tmax.plot]
      region_summary = paste(region_summary, "On ", format(end_day, "%B %d"),
                       " projections show ", format(count, digits=2, big.mark=","),
                       " cumulative deaths in ", region_name,
                       " with 90% uncertainty interval between ", format(count.min, digits=2, big.mark=","),
                       " and ", format(count.max, digits=2, big.mark=","),
                       ".",
                       sep="")    
  }
  if("rHsum" %in% which.plot || "dailyI" %in% which.plot){
      count <- max(sir_result_region_sub$mean, na.rm=TRUE)
      peak <- which.max(sir_result_region_sub$mean)
      count.min <- sir_result_region_sub$lower[peak]
      count.max <- sir_result_region_sub$upper[peak]
      name <- ifelse(which.plot ==  "rHsum",  "required hospitalizations", "daily infections")
      region_summary = paste(region_summary, "On ", format(end_day, "%B %d"),
                       " projections show a peak of ", format(count, digits=2, big.mark=","),
                       " ", name,  " reported in ", region_name,
                       " with 90% uncertainty interval between ", format(count.min, digits=2, big.mark=","), 
                       " and ", format(count.max, digits=2, big.mark=","),
                       ". The dashed line shows historical and projected hospital capacity. ",
                       sep="")    
  }
  if(!sentence){
    region_summary <- list(date=format(end_day, "%B %d"), 
                           count=format(count, digits=2, big.mark=","),
                           region_name=region_name, 
                           lower=format(count.min, digits=2, big.mark=","),
                           upper=format(count.max, digits=2, big.mark=","))
  }
  return(region_summary)
}




















#####################################

plot_ct_region_list = function(data=NULL, 
                          region_name = "Connecticut", 
                          which.plot = "D", 
                          color =  NULL,
                          title=NULL, xlab=NULL, ylab=NULL,
                          #tmax.plot = tmax,
                          #start_day = day0,
                          end_day=NULL, # pass in daymax
                          capacity_func = COUNTY_CAPACITIES,
                          obs_state = DAT_CT_STATE,
                          obs_county = DAT_CT_COUNTY, 
                          sentence=FALSE,
                          description=NULL){
  #  Get common y axis
  start_day = day0
  tmax.plot = as.numeric(difftime(end_day, day0, units="days"))

  toplot <- paste(rep(which.plot,each=length(region_name)),
                  rep(region_name, length(which.plot)),sep=".")
  ymax <- 0
  for(i in 1:length(data)){
    sir_result_region= filter(data[[i]], variable%in%toplot)
    ymax <- max(c(ymax, sir_result_region$mean[sir_result_region$time <= tmax.plot]), na.rm=TRUE)
  }
  # run multiple models
  out <- NULL
  for(i in 1:length(data)){
    out[[i]] <- plot_ct_region(data=data[[i]], region_name=region_name, which.plot = which.plot, color=color, title=paste0("\n",title[[i]]), xlab=xlab, ylab=ylab, end_day=end_day, capacity_fun=capacity_func, obs_state=obs_state, obs_county=obs_county, ymax=ymax, sentence=sentence)
  }
  
  # form a sentence
  unit <- "reported deaths"
  if(which.plot=="rHsum") unit <- "reported hospitalization required"
  if(which.plot=="dailyI") unit <- "reported daily Infections"

  out.print=NULL

  for(i in 1:length(data)){
    out.print[[i]] <- paste0(description[[i]], ", ")
    out.print[[i]] <- paste0(out.print[[i]], "projections show ", out[[i]]$count, " ", unit, 
                            " in ", out[[i]]$region_name, " on ", out[[i]]$date, 
                            ", with 90% uncertainty interval between ",
                            out[[i]]$lower, " and ", out[[i]]$upper, ".")
  }  
  return(out.print)
}



#####################################
#
# @param which.plot the prefix of the compartment to plot, e.g., S, E, I_s, I_m. If a vector of more than one specified, it takes the sum of the compartment (only the mean!)
mapplot_ct_region = function(data, which.plot = "D", label = "Cumulative Deaths", palette="Reds", subtitle=NULL, ...) {
  map = CT_MAP
  region_names <- as.character(map$NAME10)
  toplot <- paste(rep(which.plot,each=length(region_names)),
                  rep(region_names, length(which.plot)),sep=".")
  sir_result_internal = data.frame(filter(data, variable%in%toplot))
  sir_result_internal$County <- sir_result_internal$variable
  for(i in which.plot) sir_result_internal$County <- gsub(paste0(i,"\\."), "", sir_result_internal$County)
  
  # take the sum if necessary, remove lower and upper though
  if(length(which.plot) >  1){
    sir_result_internal <- aggregate(mean ~ County+time, data=sir_result_internal, FUN=function(x){mean(x)})
  }
  sir_result_internal$Date <- format(sir_result_internal$time + day0, "%B")
  sir_result_internal$Date <- factor(sir_result_internal$Date, unique(sir_result_internal$Date))
  # Plot last day cumulative at each month? Or take diff?
  t.index <- unique(sir_result_internal$time[format(sir_result_internal$time + day0, "%d")=="01"]) 
  t.index <- (t.index - 1)[-1]

  g <- mapPlot(data=subset(sir_result_internal, time %in% t.index), geo=map,
    by.data="County", by.geo = "NAME10",
    variables = "Date", values = "mean", is.long=TRUE, 
    legend.label = label, ...) #, direction=-1,  ...) 
  if(is.null(subtitle)) subtitle<-waiver()
  suppressMessages(
    g <- g + scale_fill_distiller(label, palette=palette, direction=1) + ggtitle(paste("Projected", tolower(gsub("\n"," ",label)), "in Connecticut counties by month"), subtitle=subtitle)
  )
  return(g)
}


mapplot_ct_region_list =  function(data, which.plot = "D", label = "Cumulative Deaths", palette="Reds", subtitle=NULL, ...) {
  ylim <- NULL
  map = CT_MAP
  region_names <- as.character(map$NAME10)
  for(i in 1:length(data)){
      toplot <- paste(rep(which.plot,each=length(region_names)),
                  rep(region_names, length(which.plot)),sep=".")
      sir_result_internal = data.frame(filter(data[[i]], variable%in%toplot))
      t.index <- unique(sir_result_internal$time[format(sir_result_internal$time + day0, "%d")=="01"]) 
      t.index <- (t.index - 1)[-1]
      ylim <- range(c(ylim, subset(sir_result_internal, time %in% t.index)$mean))
  }
  g <- NULL
  for(i in 1:length(data)){
    g[[i]] <- mapplot_ct_region(data=data[[i]], which.plot=which.plot, label=label, palette=palette, subtitle = subtitle[[i]], ...)
    suppressMessages(
      g[[i]] <- g[[i]] + scale_fill_distiller(label, palette=palette, direction=1, limits=ylim)
     )
  }
  return(g)
}



plot_interventions = function(sir_results, daymax, subtitle=NULL, ref_day=Sys.Date()) {


  dayseq = seq(day0, daymax, by="day")
  tmax = as.numeric(difftime(daymax, day0, units="days"))

  monthseq = seq(day0, daymax, by="month")
  monthseq_lab = format(monthseq, "%b %Y")
  daymonthseq = difftime(monthseq, day0, units="days")

  par(mar=c(3,3,3,0), bty="n")
  nn <- 4
  #layout(matrix(c(1:nn),nrow=,nn), heights=rep(2, nn))

  title.print <- "Overall contact intervention"
  if(!is.null(subtitle)) title.print <- paste0(title.print, "\n", subtitle)

  plot(sir_results[[1]]$intervention_pattern, ylim=c(0,1), xlim=c(0,1.05*tmax), type="n", ylab="", xlab="", main=title.print, axes=FALSE)
  axis(1, at=daymonthseq, lab=monthseq_lab)
  axis(2)
  polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_pattern, 0), col="orange", border=NA)
  abline(v=ref_day-day0, col="gray", lty=2)
}
