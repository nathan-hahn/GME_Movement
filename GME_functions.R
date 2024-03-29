##### GME_Movement Functions #####
# Nathan Hahn


##### rollstat #####
# function to return stats from roll apply. 
# Used inside window_stats
rollstat <- function(x, nas = TRUE) {
  # x     = numeric vector
  # na.rm = boolean, whether or not to remove NA's
  
  m  <- mean(x, na.rm = nas)
  s  <- sd(x, na.rm = nas)
  hi <- m + 2*s
  lo <- m - 2*s
  ag.window <- max(x) # all points within the window will be labeled 1 - signals a "raiding period" for activity budgets/plotting
  
  ret <- c(mean = m, stdev = s, hi.95 = hi, lo.95 = lo, ag.window = ag.window) 
  return(ret)
}

##### window_stats #####
# Get rollstats for a list of individuals. Supply list of individuals and windows(s) to apply
window_stats <- function(df.list, window, align = "center"){
  require(zoo)
  require(xts)
  system.time({
    # Calculate moving window stats
    ts <- lapply(df.list, as.ts, start = df.list[[i]]$date, frequency = 1) # frequency set with 1hr fixes
    tx <- lapply(ts, as.xts)
    res <- lapply(seq_along(tx), function(i) rollapply(tx[[i]]$ag.used, width = window, align = align,
                                                       by.column = FALSE, FUN = rollstat))
    res <- lapply(res, as.data.frame) # converting inline with Map() doesn't work
    
    # merge outputs into data frame for each individual
    merge <- Map(cbind, df.list, res)
    output <- do.call("rbind", merge)
    
    return(output)
  })
}

##### staging stats #####

r_stat <- function(x, nas = TRUE) {
  # x     = numeric vector
  # na.rm = boolean, whether or not to remove NA's
  
  m  <- mean(x, na.rm = nas)
  s  <- sd(x, na.rm = nas)
  hi <- m + 2*s
  lo <- m - 2*s
  window <- max(x) # all points within the window will be labeled 1 - signals a "raiding period" for activity budgets/plotting
  
  ret <- c(mean = m, stdev = s, hi.95 = hi, lo.95 = lo, window = window) 
  return(ret)
}

# Get rollstats for a list of individuals. Supply list of individuals and windows(s) to apply
stage_stats <- function(df.list, window, col, align = "center"){
  require(zoo)
  require(xts)
  system.time({
    # Calculate moving window stats
    ts <- lapply(df.list, as.ts, start = df.list[[i]]$date, frequency = 1) # frequency set with 1hr fixes
    tx <- lapply(ts, as.xts)
    res <- lapply(seq_along(tx), function(i) rollapply(tx[[i]][,col], width = window, align = align,
                                                       by.column = FALSE, FUN = r_stat))
    res <- lapply(res, as.data.frame) # converting inline with Map() doesn't work
    
    # merge outputs into data frame for each individual
    merge <- Map(cbind, df.list, res)
    output <- do.call("rbind", merge)
    
    return(output)
  })
}



##### plot_budget #####
# function for plotting activity budgets
plot_budget <- function(t=t, facet, title){
  require(tidyverse)
  require(lubridate)
  ggplot(t, aes(hour(date))) +
    facet_grid(facet) +
    geom_density(aes(x = hour(date),
                     fill = factor(viterbi), colour = factor(viterbi)), alpha = 0.3, adjust = 1.5) +
    # add day/night lines - code for shading below
    geom_vline(aes(xintercept=6),
               color="dark grey", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=18),
               color="dark grey", linetype="dashed", size=1) +
    
    # add colors
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
    xlab("hour (0-23)") + ggtitle(title)
}

##### velocity #####
velocity <- function(df, dateTime) {
  require(rlist)
  # calculate time differences
  tdiff <- diff(dateTime)
  units(tdiff) <- "hours"
  tdiff <- as.numeric(tdiff)
  tdiff <- list.append(tdiff, 0) # add speed zero for last fix
  
  # step lengths
  step <- df$step 
  
  # calculate velocity at each time step
  speed <- as.vector(step/tdiff)
  
  # attach to the data.frame
  df$step <- speed
  return(df)
}


##### log.velocity #####
log_velocity <- function(df) {
  require(rlist)
  # calculate time differences
  tdiff <- diff(df$date)
  units(tdiff) <- "hours"
  tdiff <- as.numeric(tdiff)
  tdiff <- list.append(tdiff, 0) # add speed zero for last fix
  
  # step lengths
  step <- df$step 
  
  # calculate velocity at each time step
  speed <- as.vector(step/tdiff)
  
  # attach to the data.frame
  df$step <- log(speed + 0.001) # add constant for zero steps
  return(df)
}


##### mask.poly.raster #####
mask_poly_raster <- function(poly, raster) {
  # crop raster to polygon extent
  cr <- crop(raster, extent(poly), snap = "out")
  # mask the raster by mcp 
  lr <- mask(x = cr, mask = poly)
  return(lr)
}


##### withold #####
# withold data
withold <- function(x, cut, type) {
  n <- nrow(x)*cut
  
  if (type == "train") {
    y <- head(x, (nrow(x)-n))
  }
  
  if (type == "test") {
    y <- tail(x, n)
  }
  
  return(y)
}


##### cluster cutpoints #####
#To extract the cutoff between any two clusters from an mclust object
cluster_cutpoint<-function(object, comp, grph=T) {
  
  mu1<-object$parameters$mean[comp[1]]
  mu2<-object$parameters$mean[comp[2]]
  sd1<-sqrt(object$parameters$variance$sigmasq[comp[1]])
  sd2<-sqrt(object$parameters$variance$sigmasq[comp[2]])
  
  if(length(object$parameters$variance$sigmasq)>1) {sd2<-sqrt(object$parameters$variance$sigmasq[comp[2]])}
  
  ss<-seq(0,.5, 0.001)
  cut<-ss[max(which(pnorm(ss, mu1, sd1, lower.tail=F)>pnorm(ss, mu2, sd2)))]
  
  if (grph) {plot(ss, pnorm(ss, mu1, sd1, lower.tail=F),col="red", type="l", xlab="Mean Ag Use", ylab="Probability", main="Prob of being in cluster")
    lines(ss,  pnorm(ss, mu2, sd2), col="blue")
    abline(v=cut, lwd=2, lty=2)}
  
  return(cut)
  
}


