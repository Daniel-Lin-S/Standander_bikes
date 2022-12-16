# Title:        Some preprocessing on the Santander Bikes data
# Author:       Francesco Sanna Passino
# Affiliation:  Department of Mathematics, Imperial College London
# Email:        f.sannapassino@imperial.ac.uk

## required library
## for qqexp
if(!require('SMPracticals')) { 
  install.packages('SMPracticals')
  library('SMPracticals')
}


####################
### Basic import ###
####################
{
  ## List all the .csv file names (and sort them in case they are not ordered by number)
  file_names = sort(list.files('/Users/linyuhang/Downloads/UROP/Data/santander_summaries/', pattern='.csv'))
  ## file_names = file_names[which(grepl('256',file_names)):which(grepl('261',file_names))]
  file_names = file_names[which(grepl('221', file_names)):which(grepl('236', file_names))]
  start_day = as.Date(strsplit(strsplit(file_names[1], '_')[[1]][2], '-')[[1]][1], '%d%b%Y')
  ## Total number of files
  n_weeks = length(file_names)
  
  ## Import
  weekly_data = list()
  for(week in 1:n_weeks){
    weekly_data[[week]] = read.table(paste('/Users/linyuhang/Downloads/UROP/Data/santander_summaries/',
                                           file_names[week],sep=''), sep=',', header=FALSE, 
                                     col.names=c('start_id','end_id','start_time','duration'))
  }
}

####################
### Basic Handle(included in the report) ###
####################
{
  library(dplyr)
  df = dplyr::bind_rows(weekly_data)
  df = transform(df, end_time = start_time + duration)
  ## change unit of duration from s to min
  df = transform(df, duration = duration / 60)
  
  ## Import stations
  stations = read.table('/Users/linyuhang/Downloads/UROP/Data/santander_locations.csv', sep=',', header=TRUE)
  
  ## Calculate geodesic distance
  library(geodist)
  # distances in km
  distances = geodist(stations[,c('longitude','latitude')], measure='vincenty') / 1000
  ## Station IDs  
  ## if entry of xx is 0, station missing
  xx = numeric(max(stations$Station.Id))
  for(i in 1:length(xx)){
    kk = which(stations$Station.Id == i)
    if(length(kk) > 0){
      xx[i] = kk
    }
  }
  ## V represents number of stations
  V <- length(xx) 
  ## Travel Distances in km
  vv = apply(df, MARGIN=1, FUN=function(x) distances[xx[x[1]], xx[x[2]]])
  df = transform(df, dist=vv)
  
  ## travel speed
  df <- transform(df, speed = dist / duration * 60)  # unit of speed: km/h
  
  ## filter out abnormal travel speed
  orig_length <- length(df$speed)
  speed_threshold <- median(df$speed) + 3*IQR(df$speed)
  # 25.7km/h, about 16mph
  df <- df[df$speed < speed_threshold, ]
  ## check not too large portion of data eliminated 
  ratio_preserved <- length(df$speed) / orig_length
}


#################################################
### analyse travels for leisure (included in the report) ###
#################################################
{
  ## package for fitting distribution
  if(!require('fitdistrplus')) { 
    install.packages('fitdistrplus')
    library('fitdistrplus')
  }
  leis_df <- df[df$dist == 0, ]
  ## drop useless columns
  leis_df <- subset(leis_df, select=c('start_id', 'start_time', 'duration', 'end_time'))
  leis_duration <- leis_df$duration
  (leis_duration_quantile <- quantile(leis_duration, prob=c(.01, .05, .1, .25, .5, .75, .9, .95, .99, .995)))
  ## more than 10% of data are below 3min, that could possibly mean taking the bike by mistake.
  ## 0.5% of data are above 440min (more than 7.3 hours) and it is possible that the users forget to return the bikes 
  ## so treat these as outliers

  leis_duration_trimmed <- leis_duration[(leis_duration > 10) & (leis_duration < quantile(leis_duration, prob=.95))]

  sink(file = "quantiles_of_duration_leisure.txt")
  quantile(leis_duration_trimmed, prob=c(.01, .05, .1, .25, .5, .75, .9, .95, .99, .995))
  sink(file = NULL)
  
  hist(leis_duration_trimmed, main="histogram of duration of journeys for leisure",  xlab="duration/min", breaks=100, xaxp = c(0, 150, 5), xlim = c(0, max(leis_duration_trimmed)))
  
  
  ## histogram looks like gamma
  fit_gamma_leis <- fitdist(leis_duration_trimmed, distr="gamma", method="mme")
  sink(file = "gamma_fit_leisure.txt")
  summary(fit_gamma_leis)
  sink(file = NULL)
  
  h <- hist(leis_duration_trimmed, breaks = 100, plot=FALSE) #generate hist
  plot(h, col="grey", main="Test of gamma fit",  xlab="duration/min", xaxp = c(0, 150, 5), xlim = c(0, max(leis_duration_trimmed))) #plot hist
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines, shape = fit_gamma_leis$estimate[1], rate = fit_gamma_leis$estimate[2]) *length(leis_duration_trimmed)*diff(h$breaks)[1], col="red")
  
  par(mfrow=c(1, 1))
  plot(fit_gamma_leis)
  ## fails, not gamma
}
 
{
  ## try discrete gamma
  f_disc <- function(x) {
    n <- 1000
    a <- x[1]
    b <- x[2]
    # for mean
    d1 <- mean(leis_duration_trimmed)
    # for variance
    d2 <- var(leis_duration_trimmed)

    for (k in 1:n) {
      I <- pgamma(k, a, rate = b) - pgamma(k-1, a, rate = b)
      d1 <- d1 - k * I
      d2 <- d2 - k^2 * I
    }
    
    d2 <- d2 + (mean(leis_duration_trimmed) - d1)^2
    return (d1^2 + d2^2)
  }
  
  ## decide numbers of estimates for infinite sum using difference 
  f_disc_diff <- function(x) {
    tol <- 1e-2
    k <- 1
    a <- x[1]
    b <- x[2]
    # for mean
    d1 <- mean(leis_duration_trimmed)
    # for variance
    d2 <- var(leis_duration_trimmed)
    
    # record previous and current value
    record <- c(d1^2 + d2^2, 0)
    diff <- abs(record[2] - record[1])
    
    prev_I <- 0
    
    while(diff > tol) {
      if (k == 5000) {
        print("maximum iteration: 5000 reached")
        break
      }
      
      I <- pgamma(k, a, rate = b) - pgamma(k-1, a, rate = b)
      d1 <- d1 - k * I
      d2 <- d2 - k^2 * I
      record[2] <- d1^2 + d2^2
      diff <- abs(record[2] - record[1])
      if (prev_I <= I) {
        diff <- 1
      }
      prev_I <- I
      record[1] <- record[2]
      k <- k + 1
    }
    
    # add the E(d)^2 to second difference
    d2 <- d2 + (mean(leis_duration_trimmed) - d1)^2
    return (d1^2 + d2^2)
  }
  
  gamma_par <- optim(c(1, 2), f_disc, method="L-BFGS-B", lower=c(1e-10, 1e-10), upper=c(Inf, Inf))
  
  f_disc_diff(gamma_par$par)
  f_disc(gamma_par$par)
  
  gamma_par_2 <- optim(c(1, 2), f_disc_diff, method="L-BFGS-B", lower=c(1e-10, 1e-10), upper=c(Inf, Inf))
  
  
  h <- hist(leis_duration_trimmed, plot=FALSE, breaks=100) #generate hist
  par(mfrow=c(1, 1))
  plot(h, col="grey") #plot hist
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines,gamma_par$par[1],rate = gamma_par$par[2]) *length(leis_duration_trimmed)*diff(h$breaks)[1], col="red")
  
  gamma_par_2[1:2]
  
  par(mfrow=c(1, 1))
  plot(h, col="grey") #plot hist
  lines(x = xlines,y=dgamma(xlines,gamma_par_2$par[1],rate = gamma_par_2$par[2]) *length(leis_duration_trimmed)*diff(h$breaks)[1], col="red")
  
  Q_disc <- function(x) {
    n <- 1000
    a <- exp(x[1])
    b <- exp(x[2])
    # for mean
    d1 <- mean(leis_duration_trimmed)
    # for variance
    d2 <- var(leis_duration_trimmed)
    
    for (k in 1:n) {
      I <- pgamma(k, a, rate = b) - pgamma(k-1, a, rate = b)
      d1 <- d1 - k * I
      d2 <- d2 - k^2 * I
    }
    
    d2 <- d2 + (mean(leis_duration_trimmed) - d1)^2
    return (d1^2 + d2^2)
  }
  
  gamma_par_Nelder <- optim(c(1, 1), Q_disc, method="Nelder-Mead")
  plot(h, col="grey") #plot hist
  lines(x = xlines,y=dgamma(xlines,exp(gamma_par_Nelder$par[1]),rate = exp(gamma_par_Nelder$par[2])) *length(leis_duration_trimmed)*diff(h$breaks)[1], col="red")
  
  Q_disc_diff <- function(x) {
    tol <- 1e-2
    k <- 1
    a <- exp(x[1])
    b <- exp(x[2])
    # for mean
    d1 <- mean(leis_duration_trimmed)
    # for variance
    d2 <- var(leis_duration_trimmed)
    
    # record previous and current value
    record <- c(d1^2 + d2^2, 0)
    diff <- abs(record[2] - record[1])
    
    prev_I <- 0
    
    while(diff > tol) {
      if (k == 5000) {
        print("maximum iteration: 5000 reached")
        break
      }
      
      I <- pgamma(k, a, rate = b) - pgamma(k-1, a, rate = b)
      d1 <- d1 - k * I
      d2 <- d2 - k^2 * I
      record[2] <- d1^2 + d2^2
      diff <- abs(record[2] - record[1])
      if (prev_I <= I) {
        diff <- 1
      }
      prev_I <- I
      record[1] <- record[2]
      k <- k + 1
    }
    
    # add the E(d)^2 to second difference
    d2 <- d2 + (mean(leis_duration_trimmed) - d1)^2
    return (d1^2 + d2^2)
  }
  
  gamma_par_Nelder_2 <- optim(c(1, 1), Q_disc_diff, method="Nelder-Mead")
  plot(h, col="grey") #plot hist
  lines(x = xlines,y=dgamma(xlines,exp(gamma_par_Nelder_2$par[1]),rate = exp(gamma_par_Nelder_2$par[2])) *length(leis_duration_trimmed)*diff(h$breaks)[1], col="red")
  
  ## return the parameters
  exp(gamma_par_Nelder_2$par)
}


## try to treat durations in every 30-minute interval as different distributions
{
  leis_duration_1 <- leis_duration_trimmed[leis_duration_trimmed < 30]
  leis_duration_2 <- leis_duration_trimmed[(leis_duration_trimmed >= 30) & (leis_duration_trimmed < 60)]
  leis_duration_3 <- leis_duration_trimmed[(leis_duration_trimmed >= 30) & (leis_duration_trimmed < 60)]
  leis_duration_4 <- leis_duration_trimmed[(leis_duration_trimmed >= 60) & (leis_duration_trimmed < 90)]
  
  fit_gamma_leis_1 <- fitdist(leis_duration_1, distr="gamma", method="mme")
  summary(fit_gamma_leis_1)
  
  fit_gamma_leis_1_mle <- fitdist(leis_duration_1, distr="gamma", method="mle")
  summary(fit_gamma_leis_1_mle)
  
  h <- hist(leis_duration_1, breaks = 100, plot=FALSE) #generate hist
  plot(h, col="grey", main="Test of gamma fit to durations below 30 minutes(MME)",  xlab="duration/min", xlim = c(0, 30))
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines, shape = fit_gamma_leis_1$estimate[1], rate = fit_gamma_leis_1$estimate[2]) *length(leis_duration_1)*diff(h$breaks)[1], col="red")
  plot(h, col="grey", main="Test of gamma fit to durations below 30 minutes(MLE)",  xlab="duration/min", xlim = c(0, 30)) 
  lines(x = xlines,y=dgamma(xlines, shape = fit_gamma_leis_1_mle$estimate[1], rate = fit_gamma_leis_1_mle$estimate[2]) *length(leis_duration_1)*diff(h$breaks)[1], col="red")
  
  fit_gamma_leis_2 <- fitdist(leis_duration_2, distr="gamma", method="mme")
  summary(fit_gamma_leis_2)
  
  fit_gamma_leis_2_mle <- fitdist(leis_duration_2, distr="gamma", method="mle")
  summary(fit_gamma_leis_2_mle)
  
  h <- hist(leis_duration_2, breaks = 100, plot=FALSE) #generate hist
  plot(h, col="grey", main="Test of gamma fit to durations 30-60 mins(MME)",  xlab="duration/min", xlim = c(30, 60)) #plot hist
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines, shape = fit_gamma_leis_2$estimate[1], rate = fit_gamma_leis_2$estimate[2]) *length(leis_duration_2)*diff(h$breaks)[1], col="red")
  plot(h, col="grey", main="Test of gamma fit to durations 30-60 mins(MLE)",  xlab="duration/min", xlim = c(30, 60)) #plot hist
  lines(x = xlines,y=dgamma(xlines, shape = fit_gamma_leis_2_mle$estimate[1], rate = fit_gamma_leis_2_mle$estimate[2]) *length(leis_duration_2)*diff(h$breaks)[1], col="red")
}



#################################################
### distribution of duration against distance ###
#################################################

## build distance bins and trim data (included in the report)
{ 
  duration_by_distance <- list()
  ## div represents number of partitions for 1km distance grid
  ## 1/div should have finite decimal places for precision e.g. div = 2, 4, 5, 10
  div <- 5
  
  ## 2h will be length of each grid
  h <- 1 / (2 * div)
  ## first entry: vector of durations for distance < 1/2div km 
  duration_by_distance[[1]] <- df[(df$dist < h),]$duration
  ## number of grids 
  num <-  ceiling(div * max(df$dist))
  for (i in 1:num) {
    duration_by_distance[[i + 1]] <- df[(df$dist < i / div + h) & (df$dist >= i / div - h),]$duration
  }
  
  # vector for x-axis for plotting
  dist_values <- seq(from = 0, to = num / div, by = 1 / div)
  
  ## find mean of durations by distance
  duration_by_distance_mean <- list()
  for (i in 1:length(duration_by_distance)) {
    duration_by_distance_mean[[i]] <- mean(duration_by_distance[[i]])
  }
  
  plot(dist_values, duration_by_distance_mean, xlab = "distance/km", ylab="mean of duration/min")
  
  ## based on this plot, I will decide a threshold for trimming data, by aborting extremely long duration
  ## first three points have different trend so they may be affected by extremely long duration. 
  duration_by_distance_mean_trimmed <- duration_by_distance_mean[4:63]
  y_pred = unlist(duration_by_distance_mean_trimmed)
  x_pred = dist_values[4:63]
  lm_dist_dur <- lm(y_pred ~ x_pred)
  sink(file = "lm_dist_dur_output.txt")
  summary(lm_dist_dur)
  sink(file = NULL)
  duration_threshold <- lm_dist_dur$coefficients[1] + lm_dist_dur$coefficients[2] * (num / div)
  ## multiply some coefficient for tolerance, may change the number 1.5
  tolerance <- 1.5
  duration_threshold <- duration_threshold * tolerance
  
  ## trim the data
  ## for travel less than 2 min, that can be a bike taken by mistake. 
  for (i in 1:num) {
    duration_by_distance[[i]] <- duration_by_distance[[i]][(duration_by_distance[[i]] < duration_threshold) & (duration_by_distance[[i]] > 2)]
  }
  ## also trim the whole df
  trimmed_df <- df[(df$duration < duration_threshold) & (df$duration > 2), ]
  ratio_preserved_2 <- nrow(trimmed_df) / nrow(df)
  
  ## now remove any distance that have less than 10 data points for accuracy
  count <- 0
  for (i in 1:(num+1)) {
    if (length(duration_by_distance[[i - count]]) <= 10) {
      duration_by_distance[[i - count]] <- NULL
      dist_values <- dist_values[-(i - count)]
      count <- count + 1
    }
  }
  print(count)
}


## summary of duration and dist
sink(file = "duration_and_dist_summary.txt")
summary(df_one_side_trimmed[c("duration", "dist")])
sink(file = NULL)
 

## minor adjustment ## (included in the report)
## WARNING: you must run last section before running this section
{
  ## plot histogram of n'th entry of duration_by_distance, with each bin of width h. 
  # prob should be a boolean indicating plot density or not.  
  plot_dur_hist <- function(n, width, prob, xlim = "default") {
    if (xlim == "default") {
      xlim <- c(min(duration_by_distance[[n]]), max(duration_by_distance[[n]]))
    }
    h <- 0.1
    if (n == 1) {
      title = paste("distribution of durations for 0-", paste(toString(h), "km journeys",sep=""), sep="")
    } else{    
      title = paste("distribution of durations for", toString(dist_values[n] - h), sep = " ")
      title = paste(title, paste(toString(dist_values[n] + h), "km journeys", sep=""), sep="-")
    }
    if (!prob) { hist(duration_by_distance[[n]], main=title, xlab="duration/min", ylab="frequency", xlim = xlim, breaks=seq(0, max(duration_by_distance[[n]]) + width, by=width)) } 
    else { hist(duration_by_distance[[n]], main=title, xlab="duration/min", ylab="density", xlim = xlim, prob=TRUE) }
  }
  
  plot_dur_hist(1, 5, FALSE)
  plot_dur_hist(2, 3, FALSE)
  plot_dur_hist(3, 2, FALSE)
  plot_dur_hist(4, 2, FALSE)
  plot_dur_hist(10, 2, FALSE)
  ## for div = 5, there is a mixture for journeys between 0.1 to 0.5km.  try to cut finer
  ## no need to trim again as the whole data frame is trimmed
  duration_by_distance_copy <- duration_by_distance
  duration_by_distance_copy[[2]] <- trimmed_df[(trimmed_df$dist < 0.15) & (trimmed_df$dist >= 0.1), ]$duration
  duration_by_distance_copy[[3]] <- trimmed_df[(trimmed_df$dist < 0.2) & (trimmed_df$dist >= 0.15), ]$duration
  duration_by_distance_copy[[4]] <- trimmed_df[(trimmed_df$dist < 0.25) & (trimmed_df$dist >= 0.2), ]$duration
  duration_by_distance_copy[[5]] <- trimmed_df[(trimmed_df$dist < 0.3) & (trimmed_df$dist >= 0.25), ]$duration
  duration_by_distance_copy[[6]] <- trimmed_df[(trimmed_df$dist < 0.35) & (trimmed_df$dist >= 0.3), ]$duration
  duration_by_distance_copy[[7]] <- trimmed_df[(trimmed_df$dist < 0.4) & (trimmed_df$dist >= 0.35), ]$duration
  duration_by_distance_copy[[8]] <- trimmed_df[(trimmed_df$dist < 0.45) & (trimmed_df$dist >= 0.4), ]$duration
  duration_by_distance_copy[[9]] <- trimmed_df[(trimmed_df$dist < 0.5) & (trimmed_df$dist >= 0.45), ]$duration
  for (i in 4:length(duration_by_distance)) {
    duration_by_distance_copy[[i+6]] <- duration_by_distance[[i]]
  }
  hist(duration_by_distance_copy[[2]], main="distribution of durations for 0.1-0.15km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[3]], main="distribution of durations for 0.15-0.2km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[4]], main="distribution of durations for 0.2-0.25km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[5]], main="distribution of durations for 0.25-0.3km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[6]], main="distribution of durations for 0.3-0.35km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[7]], main="distribution of durations for 0.35-0.4km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[8]], main="distribution of durations for 0.4-0.45km journeys", xlab="duration/min", breaks=100)
  hist(duration_by_distance_copy[[9]], main="distribution of durations for 0.45-0.5km journeys", xlab="duration/min", breaks=100)
  ## the trend for mixture persists, clearly for 0.1-0.3km distance, some are using for commuting, someone else are using for leisure
  
  ## investigate commuting purpose now, journeys >= 0.5km
  ## we should trim the extreme low speed, as we are investigating commuting purpose. 
  ## lower bound: 2km/h, considering that distance between stations is not actual distance
  ## and there can be traffic signal waiting times
  ratio_preserved_vec <- c()
  ## duration_by_distance_trimmed will only include those with travel dist > 0.3 !!!! 
  duration_by_distance_trimmed <- vector(mode="list", length=length(duration_by_distance)-1)
  ## first index in duration by distance that represents journeys for commuting
  comm_index <- 4
  for (i in comm_index:length(duration_by_distance)) {
    ## upper bound for duration
    duration_bd <- dist_values[i] * 30
    duration_by_distance_trimmed[[i-comm_index+1]] <- duration_by_distance[[i]][duration_by_distance[[i]] < duration_bd]
    ratio_preserved_vec <- c(ratio_preserved_vec, length(duration_by_distance_trimmed[[i-comm_index+1]]) / length(duration_by_distance[[i]]))
  }
  ratio_preserved_vec
  ## except 0.5-0.7, we preserved over 96% of data. That is not a great loss considering the size of data. 
  ## use this trimmed data for journeys >= 0.7km
  
  for (i in 5:length(duration_by_distance)) {
    duration_by_distance[[i]] <- duration_by_distance_trimmed[[i-4]]
  }
  plot_dur_hist(4, 2, FALSE)
  plot_dur_hist(4, 1, FALSE, xlim = c(0, 40))
  plot_dur_hist(5, 1, FALSE)
  plot_dur_hist(6, 2, FALSE)
  plot_dur_hist(7, 2, FALSE)
  plot_dur_hist(8, 2, FALSE)
  plot_dur_hist(30, 2, FALSE)
  plot_dur_hist(50, 2, FALSE)
}


#### investigate mixtures in small distances ####
{
  if(!require('mixtools')) { 
    install.packages('mixtools')
    library('mixtools')
  }
  
  # test on 0.1-0.15km journeys
  gammamix2 <- gammamixEM(duration_by_distance_copy[[2]], lambda=c(1/2, 1/2), alpha=NULL, beta=NULL, k=2, maxit=1000, maxrestarts = 20)
  gammamix2[2:3]
  
  h <- hist(duration_by_distance_copy[[2]], breaks=100)
  par(mfrow=c(1, 1))
  plot(h, col="grey", main="distribution of durations for 0.1-0.15km journeys", xlab="duration/min")
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines,gammamix2$gamma.pars[1, 1],rate = gammamix2$gamma.pars[2, 1]) *length(duration_by_distance_copy[[2]])* diff(h$breaks)[1] * gammamix2$lambda[1], col="red")
  lines(x = xlines,y=dgamma(xlines,gammamix2$gamma.pars[1, 2],rate = gammamix2$gamma.pars[2, 2]) *length(duration_by_distance_copy[[2]])* diff(h$breaks)[1] * gammamix2$lambda[2], col="red")
  
  # 0.1 - 0.3 km
  gammamix <- gammamixEM(duration_by_distance[[2]], lambda=c(1/2, 1/2), alpha=NULL, beta=NULL, k=2, maxit=1000, maxrestarts = 20)
  gammamix[2:3]
  
  h <- hist(duration_by_distance[[2]], breaks=100)
  par(mfrow=c(1, 1))
  plot(h, col="grey", main="distribution of durations for 0.1-0.3km journeys", xlab="duration/min")
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines,gammamix$gamma.pars[1, 1],rate = gammamix$gamma.pars[2, 1]) *length(duration_by_distance[[2]])* diff(h$breaks)[1] * gammamix$lambda[1], col="red")
  lines(x = xlines,y=dgamma(xlines,gammamix$gamma.pars[1, 2],rate = gammamix$gamma.pars[2, 2]) *length(duration_by_distance[[2]])* diff(h$breaks)[1] * gammamix$lambda[2], col="red")
  
  h <- hist(duration_by_distance[[2]], breaks=100, prob=TRUE, main="distribution of durations for 0.1-0.3km journeys", xlab="duration/min")
  par(mfrow=c(1, 1))
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100) #seq of x for pdf
  lines(x = xlines,y=dgamma(xlines,gammamix$gamma.pars[1, 1],rate = gammamix$gamma.pars[2, 1]) * gammamix$lambda[1], col="red")
  lines(x = xlines,y=dgamma(xlines,gammamix$gamma.pars[1, 2],rate = gammamix$gamma.pars[2, 2]) * gammamix$lambda[2], col="red")
  
  plot(0, 0, xlim = c(0, 10), ylim = c(0, 1), type = "n")
  curve(dgamma(x, shape = gammamix$gamma.pars[1, 1], rate = gammamix$gamma.pars[2, 1]), from = 0, to = 10, col = i, add = TRUE)
  curve(dgamma(x, shape = gammamix$gamma.pars[1, 2], rate = gammamix$gamma.pars[2, 2]), from = 0, to = 10, col = i, add = TRUE)
}
  
#### investigate frequency against distance(included in the report) ####
{
  
  frequency_by_distance <- lapply(duration_by_distance, length)
  # spike at 0 + a distribution(looks like gamma)
  plot(dist_values, frequency_by_distance, main="Frequency of journeys by distance", xlab="distance/km", ylab="frequency", pch = 20)
  pi_hat <- length(trimmed_df[trimmed_df$dist == 0, ]$dist) / length(trimmed_df$dist)
  ## find mean distance for usage for commuting 
  mean_dist_com <- mean(trimmed_df[trimmed_df$dist > 0, ]$dist)
  var_dist_com <- var(trimmed_df[trimmed_df$dist > 0, ]$dist)
  alpha_hat <- mean_dist_com^2 / var_dist_com
  beta_hat <- var_dist_com / mean_dist_com
  frequency_by_distance <- unlist(frequency_by_distance)
  density_by_distance <- frequency_by_distance / sum(frequency_by_distance)
  plot(dist_values[-1], density_by_distance[-1], main="Density of journeys by distance", xlab="distance/km", ylab="Density", pch = 20)
  xlines <- seq(0, max(dist_values), length.out=100)
  lines(x = xlines, y = dgamma(xlines, alpha_hat, rate = beta_hat) * 0.17, col="red")
}


#### mode of duration against distance(included in the report) ####
{
  ## for mode function
  if(!require('modeest')) {
    install.packages('modeest')
    library('modeest')
  }
  ## find mode of durations by distance
  duration_by_distance_mode <- list()
  for (i in 1:length(duration_by_distance)) { 
    duration_by_distance_mode[[i]] <- mlv(duration_by_distance[[i]], method="mfv")
    
    ## for small dataset, multiple modes may exist, use the mean of modes
    if (length(duration_by_distance_mode[[i]]) > 1) {
      duration_by_distance_mode[[i]] <- mean(duration_by_distance_mode[[i]])
    }
  }
  
  ## plot mode of duration against distance
  plot(dist_values, duration_by_distance_mode, xlab = "distance/km", ylab="mode of duration/min")
  
  num <- length(dist_values)
  
  ## investigate relation between mode of duration and distance
  duration_by_distance_mode_vec <- unlist(duration_by_distance_mode[2:num])
  dist_values_nozero <- dist_values[2:num]
  lmod_dur_mode_dist <- lm(duration_by_distance_mode_vec ~ dist_values_nozero)
  summary(lmod_dur_mode_dist)
  
  ## plot fitted model
  par(mfrow=c(1, 1))
  plot(dist_values_nozero, duration_by_distance_mode_vec, main="fitted model", xlab = "distance/km", ylab="mode of duration/min", pch = 20)
  abline(lmod_dur_mode_dist$coef, col="red")
  
  ## diagnostic plot
  par(mfrow=c(2, 2))
  plot(lmod_dur_mode_dist)
  par(mfrow=c(1, 1))
  
  ## remove outliers, strong evidence shows that 84 is outlier.
  dist_values_nozero_trimmed <- dist_values_nozero[-c(83, 84)]
  duration_by_distance_mode_vec_trimmed <- duration_by_distance_mode_vec[-c(83, 84)]
  lmod_dur_mode_dist <- lm(duration_by_distance_mode_vec_trimmed ~ dist_values_nozero_trimmed)
  par(mfrow=c(2, 2))
  plot(lmod_dur_mode_dist)
  par(mfrow=c(1, 1))
  
  sink(file="lmod_dur_mode_dist.txt")
  summary(lmod_dur_mode_dist)
  sink(file=NULL)
  ## better than before
  ## but issue of heteroscedasticity and non-normality of residuals persist so errors could be Cauchy distributed.
  
  ## try robust regression
  require(MASS)
  ltsmod_dur_mode_dist <- ltsreg(duration_by_distance_mode_vec ~ dist_values_nozero)
  coef(ltsmod_dur_mode_dist)
  
  ## plot two models together
  par(mfrow=c(1, 1))
  plot(dist_values_nozero, duration_by_distance_mode_vec, main="comparison of robust regression and LS model", xlab = "distance/km", ylab="mode of duration/min", pch = 20)
  abline(lmod_dur_mode_dist$coef, col="red", lty=5)
  abline(ltsmod_dur_mode_dist$coefficients, col="red")
  
  ## there is no much difference between two models, so we can just use the conclusion of LS model. (with the three outliers removed)
}

#### IQR, sd of duration against distance(included in the report) ####
{
  ## find IQR of durationsby distance
  duration_by_distance_IQR <- list()
  for (i in 1:length(duration_by_distance)) {
    duration_by_distance_IQR[[i]] <- IQR(duration_by_distance[[i]])
  }
  
  plot(dist_values, duration_by_distance_IQR, xlab = "distance/km", ylab="IQR/min")
  
  ## there are clearly three segments 
  duration_by_distance_IQR_vector_1 <- unlist(duration_by_distance_IQR[1:4])
  duration_by_distance_IQR_vector_2 <- unlist(duration_by_distance_IQR[5:26])
  duration_by_distance_IQR_vector_3 <- unlist(duration_by_distance_IQR[27:num])
  lmod_dur_IQR_dist_1 <- lm(duration_by_distance_IQR_vector_1 ~ dist_values[1:4])
  # 5 -- 0.8km
  lmod_dur_IQR_dist_2 <- lm(duration_by_distance_IQR_vector_2 ~ dist_values[5:26])
  # 27 -- 5km
  lmod_dur_IQR_dist_3 <- lm(duration_by_distance_IQR_vector_3 ~ dist_values[27:num])
  plot(dist_values, duration_by_distance_IQR, main = "fitted model", xlab = "distance/km", ylab="IQR/min")
  abline(v=0.8, lty=5, col="grey")
  abline(v=5, lty=5, col="grey")
  segments(0, lmod_dur_IQR_dist_1$coef[1] + lmod_dur_IQR_dist_1$coef[2] * 0, 0.8, lmod_dur_IQR_dist_1$coef[1] + lmod_dur_IQR_dist_1$coef[2] * 0.8, col="red")
  segments(5, lmod_dur_IQR_dist_2$coef[1] + lmod_dur_IQR_dist_2$coef[2] * 5, 0.8, lmod_dur_IQR_dist_2$coef[1] + lmod_dur_IQR_dist_2$coef[2] * 0.8, col="red")
  segments(5, lmod_dur_IQR_dist_3$coef[1] + lmod_dur_IQR_dist_3$coef[2] * 5, max(dist_values), lmod_dur_IQR_dist_3$coef[1] + lmod_dur_IQR_dist_3$coef[2] * max(dist_values), col="red")
  
  ## the last segment seems to be distorted by some outliers 
  par(mfrow=c(2, 2))
  plot(lmod_dur_IQR_dist_3)
  par(mfrow=c(1, 1))
  ## 49, 50, 51(index within 27-77=num of original IQR vector) turns out to be outliers
  duration_by_distance_IQR_vector_3 <- duration_by_distance_IQR_vector_3[-c(58, 55)]
  dist_values_3 <- dist_values[27:num]
  dist_values_3 <- dist_values_3[-c(58, 55)]
  lmod_dur_IQR_dist_3_trimmed <- lm(duration_by_distance_IQR_vector_3 ~ dist_values_3)
  
  ## diagnostic of trimmed model
  par(mfrow=c(2, 2))
  plot(lmod_dur_IQR_dist_3_trimmed)
  par(mfrow=c(1, 1))
  # there are still three outliers. Try robust regression(Huber's method)
  require(MASS)
  rlmod_dur_IQR_dist_3 <- rlm(duration_by_distance_IQR_vector_3 ~ dist_values_3)
  summary(rlmod_dur_IQR_dist_3)
  summary(lmod_dur_IQR_dist_3_trimmed)
  ## significance confirmed, magnitude almost the same 
  
  plot(dist_values, duration_by_distance_IQR, main = "fitted model", xlab = "distance/km", ylab="IQR/min")
  abline(v=0.8, lty=5, col="grey")
  abline(v=5, lty=5, col="grey")
  segments(0, lmod_dur_IQR_dist_1$coef[1] + lmod_dur_IQR_dist_1$coef[2] * 0, 0.8, lmod_dur_IQR_dist_1$coef[1] + lmod_dur_IQR_dist_1$coef[2] * 0.8)
  segments(5, lmod_dur_IQR_dist_2$coef[1] + lmod_dur_IQR_dist_2$coef[2] * 5, 0.8, lmod_dur_IQR_dist_2$coef[1] + lmod_dur_IQR_dist_2$coef[2] * 0.8)
  segments(5, rlmod_dur_IQR_dist_3$coef[1] + rlmod_dur_IQR_dist_3$coef[2] * 5, max(dist_values), lmod_dur_IQR_dist_3$coef[1] + lmod_dur_IQR_dist_3$coef[2] * max(dist_values))
  
  ## find sd of duration by distance
  duration_by_distance_sd <- list()
  for (i in 1:length(duration_by_distance)) {
    duration_by_distance_sd[[i]] <- sd(duration_by_distance[[i]])
  }

  plot(dist_values, duration_by_distance_sd, xlab = "distance/km", ylab="sd of duration/min")
  
  ## again, investigate relationship of distance and sd of duration
  require(splines)
  knots <- c(0, 0, 0, 0, 2.5, 10, 16.8, 16.8, 16.8, 16.8)
  bx <- splineDesign(knots, dist_values)
  duration_by_distance_sd <- unlist(duration_by_distance_sd)
  lmodb_dur_sd_dist <- lm(duration_by_distance_sd ~ bx -1)
  matplot(dist_values, cbind(duration_by_distance_sd, lmodb_dur_sd_dist$fit), type="pl", xlab = "distance/km", ylab="sd of duration/min", pch=20, lty=1, col=1)
  sink(file="spline_dur_sd_dist.txt")
  summary(lmodb_dur_sd_dist)
  sink(file=NULL)
  matplot(dist_values, bx, type="l", col=1, xlab="distance/km", ylab="base function value")
}


#### duration against distance(included in the report) ####
{
  ## package for 2D histogram
  if(!require('hexbin')) {
    install.packages('hexbin')
    library('hexbin')
  }
  ## draw 2D histogram of duration against distance (finer bins)
  ## for all the stations
  plot(hexbin(log(df_one_side_trimmed$duration) ~ df_one_side_trimmed$dist), mincnt = 100, main="all stations", xlab="distance/km", ylab="log(duration/min)")
  hist(trimmed_df[trimmed_df$dist == 0, ]$duration, breaks=100)
  hist(trimmed_df[(trimmed_df$dist >= 1) & (trimmed_df$dist < 1.1), ]$duration, breaks=1000, xlim=c(0, 50))
  
  ## study the cluster close to y-axis
  plot(hexbin(log(df_one_side_trimmed$duration[(df_one_side_trimmed$dist <= 0.5) & (df_one_side_trimmed$dist > 0)]) ~ df_one_side_trimmed$dist[(df_one_side_trimmed$dist <= 0.5) & (df_one_side_trimmed$dist > 0)]), mincnt = 100, main="0 < distance < 0.5", xlab="distance/km", ylab="log(duration/min)")
  
  ## check this trend for particular stations(randomly chosen)
  set.seed(2022)
  rnd_stations <- sample(stations$Station.Id, size = 5, replace=FALSE)
  sample1_df <- trimmed_df[trimmed_df$start_id == rnd_stations[1], ]
  sample2_df <- trimmed_df[trimmed_df$start_id == rnd_stations[2], ]
  sample3_df <- trimmed_df[trimmed_df$start_id == rnd_stations[3], ]
  sample4_df <- trimmed_df[trimmed_df$start_id == rnd_stations[4], ]
  sample5_df <- trimmed_df[trimmed_df$start_id == rnd_stations[5], ]
  
  plot(hexbin(log(sample1_df$duration) ~ sample1_df$dist), mincnt = 5, main="station 233", xlab="distance/km", ylab="log(duration/min)")
  plot(hexbin(log(sample2_df$duration) ~ sample2_df$dist), mincnt = 5, main="station 459", xlab="distance/km", ylab="log(duration/min)")
  plot(hexbin(log(sample3_df$duration) ~ sample3_df$dist), mincnt = 5, main="station 748", xlab="distance/km", ylab="log(duration/min)")
  plot(hexbin(log(sample4_df$duration) ~ sample4_df$dist), mincnt = 5, main="station 738", xlab="distance/km", ylab="log(duration/min)")
  plot(hexbin(log(sample5_df$duration) ~ sample5_df$dist), mincnt = 5, main="station 805", xlab="distance/km", ylab="log(duration/min)")
  ## indeed same trend for all start stations
  
  ## the relation looks like quadratic, try to plot quadratic
  plot(hexbin(log(trimmed_df$duration)^2 ~ trimmed_df$dist), mincnt = 7000, main="all stations", xlab="distance/km", ylab="log(duration/min)^2")
  ## not really quadratic, require polynomial model here
  
  ## we only consider usage for commuting here
  ## threshold is just 0 as in the frequency plot before. 
  lmod_dur_dist_squared <- lm(dist ~ log(duration) + I(log(duration)^2), trimmed_df[trimmed_df$dist > 0, ])
  summary(lmod_dur_dist_squared)
  lmod_dur_dist_cube <- lm(dist ~ log(duration) + I(log(duration)^2) + I(log(duration)^3), trimmed_df[trimmed_df$dist > 0, ])
  summary(lmod_dur_dist_cube)
  lmod_dur_dist_4 <- lm(dist ~ log(duration) + I(log(duration)^2) + I(log(duration)^3) + I(log(duration)^4), trimmed_df[trimmed_df$dist > 0, ])
  sink(file = "lmod_dur_dist_4.txt")
  summary(lmod_dur_dist_4)
  sink(file=NULL)
  ## although fourth power term is significant, the magnitude is small compared to the other coefficients. 
  ## reassure by orthogonal polynomials
  lmod_dur_dist_ortho <- lm(dist ~ poly(duration, 10), trimmed_df[trimmed_df$dist > 0, ])
  ## up to 9'th power, could be overfitting 
  summary(lmod_dur_dist_ortho)
}



#### duration conditional on distance(included in the report) ####
{
  require(fitdistrplus)
  ## fit gamma, Weibull for entry n of duration by distance
  fit_gamma_dur_dist <- function(n) {
    fit_dur_dist_gamma <- fitdist(duration_by_distance[[n]], distr="gamma")
    print(paste("gamma likelihood: ", toString(fit_dur_dist_gamma$loglik)))
    fit_dur_dist_weib <- fitdist(duration_by_distance[[n]], distr="weibull")
    print(paste("weibull likelihood: ", toString(fit_dur_dist_weib$loglik)))
    par(mfrow=c(1, 1))
    plot_dur_hist(n, 2, FALSE)
    plot(fit_dur_dist_gamma, pch=20)
    plot(fit_dur_dist_weib, pch=20)
  }
  fit_gamma_dur_dist(36)
  fit_dur_dist_weib_36 <- fitdist(duration_by_distance[[36]], distr="weibull")
  summary(fit_dur_dist_weib_36)
  fit_gamma_dur_dist(50)
  fit_gamma_dur_dist(10)
  ## fit inverse gamma, weibull for entry n of duration by distance and compare
  fit_icomp_dur_dist <- function(n, sep = 2, plot=TRUE) {
    inv_dur <- 1 / duration_by_distance[[n]]
    fitG <- fitdist(inv_dur, distr="gamma")
    fitW <- fitdist(inv_dur, distr="weibull")
    if (plot) {
      par(mfrow=c(1, 1))
      plot_dur_hist(n, sep, FALSE)
      cdfcomp(list(fitW, fitG), legendtext=c("inverse Weibull","inverse gamma"))
      denscomp(list(fitW, fitG), legendtext=c("inverse Weibull", "inverse gamma"))
      qqcomp(list(fitW, fitG), legendtext=c("inverse Weibull", "inverse gamma"), pch=20)
    }
    
    list("gamma_fit" = fitG, "Weibull_fit" = fitW)
  }
  fit_icomp_dur_dist(36)
  fit_icomp_dur_dist(50)
  result_15 <- fit_icomp_dur_dist(15, sep=1, plot=FALSE)
  result_15$gamma_fit$loglik
  result_15$Weibull_fit$loglik
  
  ## positive loglike difference if in favour of Weibull
  n <- length(duration_by_distance)
  choice_GW_loglik <- c()
  choice_GW_sd <- c()
  parameters_Weib <- vector(mode="list", length=n)
  for (i in 1:n) {
    print(paste("testing on i = ", toString(i)))
    GW_fitness <- fit_icomp_dur_dist(i, plot = FALSE)
    loglik_diff <- GW_fitness$Weibull_fit$loglik - GW_fitness$gamma_fit$loglik
    G_sd <- GW_fitness$gamma_fit$sd[[1]] + GW_fitness$gamma_fit$sd[[2]]
    W_sd <- GW_fitness$Weibull_fit$sd[[1]] + GW_fitness$Weibull_fit$sd[[2]]
    Weibull_better_sd <- (W_sd < G_sd)
    choice_GW_loglik <- c(choice_GW_loglik, loglik_diff)
    choice_GW_sd <- c(choice_GW_sd, Weibull_better_sd)
    parameters_Weib[[i]] <- GW_fitness$Weibull_fit$estimate
  }

  plot(dist_values, choice_GW_loglik, main="loglikelihood difference", xlab="distance/km")
  length(which(choice_GW_sd))
  dur_dist_Weib_shape <- c()
  dur_dist_Weib_scale <- c()
  for (i in 1:n) {
    dur_dist_Weib_shape <- c(dur_dist_Weib_shape, parameters_Weib[[i]][1])
    dur_dist_Weib_scale <- c(dur_dist_Weib_scale, parameters_Weib[[i]][2])
  }
  plot(dist_values, dur_dist_Weib_shape, main="shape parameter of Weibull component", xlab="distance/km", pch = 20)
  plot(dist_values, dur_dist_Weib_scale, main="scale parameter of Weibull component", xlab="distance/km", pch = 20)
  hat_Weib_mean <- dur_dist_Weib_scale * gamma(1 + 1 / dur_dist_Weib_shape)
  hat_Weib_var <- dur_dist_Weib_scale ^ 2 * (gamma(1 + 2 / dur_dist_Weib_shape) - (gamma(1 + 1 / dur_dist_Weib_shape))^2)
  plot(dist_values, hat_Weib_mean, main="mean of 1 / duration", xlab="distance/km", pch = 20)
  plot(dist_values, hat_Weib_var, main="variance of 1 / duration", xlab="distance/km", pch = 20)
}


###################################################################
### distribution of duration according to start and end station(included in the report) ###
###################################################################


## aborted method: study distribution of individual edge. Model fit fails

{
  library(dplyr)
  duration_summary_by_station <- trimmed_df %>% group_by(start_id, end_id) %>% summarise(
    mean_duration = mean(duration), 
    sd_duration = sd(duration), 
    number = n(),
    dist = mean(dist)
  )
  
  popularity_order_duration_summary_by_station <- duration_summary_by_station[order(duration_summary_by_station$number, decreasing = TRUE), ]
  
  ## return the ratio of edges with start_id = end_id for the n most popular edges
  leis_ratio_by_station <- function(n) {
    pop_duration_summary <- head(popularity_order_duration_summary_by_station, n)
    length(pop_duration_summary[(pop_duration_summary$start_id == pop_duration_summary$end_id),]$start_id) / length(pop_duration_summary$start_id)
  }
  
  leis_ratio_by_station(100)
  leis_ratio_by_station(500)
  popularity <- seq(from = 50, to = 3000, by = 50)
  plot(popularity, sapply(popularity, leis_ratio_by_station), ylab="ratio", main = "ratio of self-connecting edges(start_id = end_id)  \n among the most popular edges", pch=20)
  
  nrow(duration_summary_by_station[duration_summary_by_station$number<100, ]) / nrow(duration_summary_by_station)
  ## more than 99% percent of the edges have less than 10 journeys. But we can analyse the popular edges
  
  plot_dur_hist_by_station <- function(start_id, end_id, breaks=100) {
    dur_data <- df_one_side_trimmed[(df_one_side_trimmed$start_id == start_id) & (df_one_side_trimmed$end_id == end_id), ]$duration
    k = ceiling(max(dur_data) / 30)
    x_upper <-k  * 30
    hist(dur_data, breaks = breaks, xlab = "duration/min", main= paste("Start station: ", toString(start_id), ", end station: ", toString(end_id), sep=""), xlim = c(0, x_upper), xaxp = c(0, x_upper, k))
  }
  
  plot_inverse_dur_hist_by_station <- function(start_id, end_id, breaks=100) {
    dur_data <- df_one_side_trimmed[(df_one_side_trimmed$start_id == start_id) & (df_one_side_trimmed$end_id == end_id), ]$duration
    hist(1/dur_data, breaks = breaks, xlab = "1/duration", main= paste("Start station: ", toString(start_id), ", end station: ", toString(end_id), sep=""))
  }
  
  plot_dur_hist_by_station(785, 785)
  plot_dur_hist_by_station(191, 191)
  plot_dur_hist_by_station(111, 111)
  plot_inverse_dur_hist_by_station(111, 111)
  plot_dur_hist_by_station(248, 248)
  plot_dur_hist_by_station(789, 789)
  plot_dur_hist_by_station(350, 350)
  plot_dur_hist_by_station(303, 303)
  plot_dur_hist_by_station(787, 787)
  plot_dur_hist_by_station(786, 786)
  plot_dur_hist_by_station(114, 114)
  plot_dur_hist_by_station(99, 99)
  plot_dur_hist_by_station(532, 532)
  plot_dur_hist_by_station(811, 811)
  plot_dur_hist_by_station(547, 547)
  

  plot_dur_hist_by_station(785, 789)
  plot_dur_hist_by_station(789, 785)
  
  
  plot_dur_hist_by_station(191, 248)
  require(fitdistrplus)
  durations_191_248 <- df_one_side_trimmed[(df_one_side_trimmed$start_id == 191) & (df_one_side_trimmed$end_id == 248), ]$duration
  fit_gamma_191_248 <- fitdist(durations_191_248, distr="gamma", method="mle")
  fit_gamma_191_248$loglik
  plot(fit_gamma_191_248, pch=20)
  h <- hist(durations_191_248, plot=FALSE, breaks=100)
  plot(h, col="grey", xaxp = c(0, 180, 6), main="Start_id: 191, end_id:248")
  xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
  lines(x = xlines, y=dgamma(xlines, shape = fit_gamma_191_248$estimate[1], rate = fit_gamma_191_248$estimate[2]) *length(durations_191_248)*diff(h$breaks)[1], col="red", lty=2)
  # try inverse gamma
  fit_igamma_191_248 <- fitdist(1 / durations_191_248, distr="gamma", method="mle")
  plot(fit_igamma_191_248, pch=20)
  fit_igamma_191_248$loglik
  ## better fit than gamma 
  
  durations_768_833 <- trimmed_df[(trimmed_df$start_id == 768) & (trimmed_df$end_id == 833), ]$duration
  hist_dur_768_833 <- hist(durations_768_833, main = "histogram from station 768 to station 833", xlab="duration/min", breaks=100)
  ## cut off data above 60min
  durations_768_833 = durations_768_833[durations_768_833 < 60]
  hist_dur_768_833 <- hist(durations_768_833, main = "histogram from station 768 to station 833", xlab="duration/min", breaks=50)
  fit_gamma_768_833 <- fitdist(durations_768_833, distr="gamma", method="mle")
  summary(fit_gamma_768_833)
  plot(fit_gamma_768_833, pch=20)
  ## not good fit, try moment match
  fit_gamma_768_833 <- fitdist(durations_768_833, distr="gamma", method="mme")
  summary(fit_gamma_768_833)
  plot(fit_gamma_768_833, pch=20)
  fit_gamma_768_833$loglik
  
  ## compare with other distributions with similar shape
  fit_lnorm_768_833 <- fitdist(durations_768_833, distr="lnorm")
  fit_weib_768_833 <- fitdist(durations_768_833, distr="weibull")
  cdfcomp(list(fit_weib_768_833, fit_gamma_768_833, fit_lnorm_768_833), legendtext=c("Weibull", "gamma", "lognormal"))
  denscomp(list(fit_weib_768_833, fit_gamma_768_833, fit_lnorm_768_833), legendtext=c("Weibull", "gamma", "lognormal"))
  qqcomp(list(fit_weib_768_833, fit_gamma_768_833, fit_lnorm_768_833), legendtext=c("Weibull", "gamma", "lognormal"))
  ppcomp(list(fit_weib_768_833, fit_gamma_768_833, fit_lnorm_768_833), legendtext=c("Weibull", "gamma", "lognormal"))
  
  ## try inverse gamma
  fit_igamma_768_833 <- fitdist(1/durations_768_833, distr="gamma")
  plot(fit_igamma_768_833, pch=20)
  
  ## fit gamma and inverse gamma for journeys between station 
  fit_igamma_gamma <- function(start_id, end_id, plot = FALSE) {
    durations_ij <- trimmed_df[(trimmed_df$start_id == start_id) & (trimmed_df$end_id == end_id), ]$duration
    fit_igamma_ij <- fitdist(1/durations_ij, distr="gamma")
    fit_gamma_ij <- fitdist(durations_ij, distr="gamma")
    if (plot) {
      plot(fit_igamma_ij, pch=20)
      plot(fit_gamma_ij, pch=20)
    }
    
    list("inverse_gamma_fit" = fit_igamma_ij, "gamma_fit" = fit_gamma_ij)
  }
  
  ## test other stations
  fitness_igamma_191_248 <- fit_igamma_gamma(191, 248)
  fitness_igamma_191_248$inverse_gamma_fit$
  fit_igamma_gamma(785, 789)
  fit_igamma_gamma(248, 191)
  fit_igamma_gamma(786, 785)
  fit_igamma_gamma(789, 787)
  
  pop500_commuting_edges <- head(popularity_order_duration_summary_by_station, 500)
  pop500_commuting_edges_start <- pop500_commuting_edges[pop500_commuting_edges$start_id != pop500_commuting_edges$end_id, ]$start_id
  pop500_commuting_edges_end <- pop500_commuting_edges[pop500_commuting_edges$start_id != pop500_commuting_edges$end_id, ]$end_id
  
  loglik_igamma_stations <- c()
  loglik_gamma_stations <- c()
  sd_igamma_stations <- c()
  sd_gamma_stations <- c()
  for (i in 1:length(pop500_commuting_edges_end)) {
    fitness_result <- fit_igamma_gamma(pop500_commuting_edges_start[i], pop500_commuting_edges_end[i])
    inv_gamma_fit <- fitness_result$inverse_gamma_fit
    gamma_fit <- fitness_result$gamma_fit
    
    loglik_igamma_stations <- c(loglik_igamma_stations, inv_gamma_fit$loglik)
    loglik_gamma_stations <- c(loglik_gamma_stations, gamma_fit$loglik)
    sd_igamma_stations <- c(sd_igamma_stations, inv_gamma_fit$sd[1]^2 + inv_gamma_fit$sd[2]^2)
    sd_gamma_stations <- c(sd_gamma_stations, gamma_fit$sd[1]^2 + gamma_fit$sd[2]^2)
  }
  
  plot(loglik_igamma_stations, pch=20, xlab="index", ylab="log likelihood", main="log likelihood for fitting inverse gamma")
  plot(loglik_gamma_stations, pch=20, xlab="index", ylab="log likelihood", main="log likelihood for fitting gamma")
  # remove extreme points
  sd_igamma_stations <- sd_igamma_stations[-c(118, 193, 228, 256, 266, 276, 286)]
  plot(sd_igamma_stations, pch=20, xlab="index", ylab="sum squared of sds", main="sum of standard deviations squared \n for fitting inverse gamma")
  sd_gamma_stations <- sd_gamma_stations[-c(118, 125, 150, 197, 254, 293)]
  plot(sd_gamma_stations, pch=20, xlab="index", ylab="sum squared of sds", main="sum of standard deviations squared for fitting gamma")
  
  
}





#########################################
### distance conditioned on duration(included in the report) ####
#########################################
{
  ## some examples 
  dist_5min <- trimmed_df[trimmed_df$duration == 5, ]$dist
  hist(dist_5min, breaks=50, main="travel distances for travels with duration 5min", xlab="distance/km")
  
  dist_10min <- trimmed_df[trimmed_df$duration == 10, ]$dist
  hist(dist_10min, breaks=50, main="travel distances for travels with duration 10min", xlab="distance/km")
  
  dist_60min <- trimmed_df[trimmed_df$duration == 60, ]$dist
  hist(dist_60min, breaks=100, main="travel distances for travels with duration 60min", xlab="distance/km")
  
  dist_120min <- trimmed_df[trimmed_df$duration == 120, ]$dist
  hist(dist_120min, breaks=100, main="travel distances for travels with duration 120min", xlab="distance/km")
  
  ## restore values for duration = 1/2 min
  df_one_side_trimmed <- df[(df$duration < duration_threshold), ]
  
  dist_1min <- df_one_side_trimmed[df_one_side_trimmed$duration == 1, ]$dist
  hist(dist_1min, breaks=50, main="travel distances for travels with duration 1min", xlab="distance/km")
  
  dist_1min_non_zero <- df_one_side_trimmed[(df_one_side_trimmed$duration == 1) & (df_one_side_trimmed$dist > 0), ]$dist
  hist(dist_1min_non_zero, breaks=50, main="travel distances(non-zero) for travels with duration 1min", xlab="distance/km")
  z_dist_1min <- (dist_1min_non_zero - mean(dist_1min_non_zero)) / sd(dist_1min_non_zero)
  qqnorm(z_dist_1min)
  abline(0, 1, col="red")
  
  dist_2min <- df_one_side_trimmed[df_one_side_trimmed$duration == 2, ]$dist
  hist(dist_2min, breaks=50, main="travel distances for travels with duration 2min", xlab="distance/km")
  dist_2min_non_zero <- dist_2min[dist_2min > 0]
  hist(dist_2min_non_zero, breaks=100, main="travel distances(non-zero) for travels with duration 2min", xlab="distance/km")
  z_dist_2min <- (dist_2min_non_zero - mean(dist_2min_non_zero)) / sd(dist_2min_non_zero)
  qqnorm(z_dist_2min)
  abline(0, 1, col="red")
  
  ## investigate relation of pi_0 (possibility of dist = 0) with duration
  pi_0 <- c()
  t_values <- min(df_one_side_trimmed$duration):max(df_one_side_trimmed$duration)
  for (t in t_values) {
    dist_con_dur <- df_one_side_trimmed[df_one_side_trimmed$duration == t, ]$dist
    dist_zero_count <- length(dist_con_dur[dist_con_dur == 0])
    pi_0 <- c(pi_0, dist_zero_count / length(dist_con_dur))
  }
  
  plot(t_values, pi_0, pch=20, main="proportion of zero distance against duration", xlab="duration/s")
  
  ## pi_0 vs duration
  require(splines)
  max_dur <- max(df_one_side_trimmed$duration)
  min_dur <- min(df_one_side_trimmed$duration)
  knots <- c(min_dur, min_dur, min_dur, min_dur, 5, 60, max_dur, max_dur, max_dur, max_dur)
  bx <- splineDesign(knots, t_values)
  lmod_spline_pi0_vs_duration <- lm(pi_0 ~ bx -1)
  matplot(t_values, cbind(pi_0, lmod_spline_pi0_vs_duration$fit), type="pl", main = "fit test of pi_0", ylab=expression(pi_0), pch=20, lty=1, col=c(1,"red"))
  ## check base functions
  matplot(t_values, bx, type="l", col=1)
  sink(file="lmod_spline_pi0_vs_duration.txt")
  summary(lmod_spline_pi0_vs_duration)
  sink(file=NULL)
  
  ## investigate relation of pi_0 (possibility of dist = 0), mean and variance (of non-zero distance travels) with duration
  mu <- c()
  var <- c()
  for (t in t_values) {
    dist_con_dur <- df_one_side_trimmed[df_one_side_trimmed$duration == t, ]$dist
    dist_con_dur_non_zero <- dist_con_dur[dist_con_dur > 0]
    mu <- c(mu, mean(dist_con_dur_non_zero))
    var <- c(var, var(dist_con_dur_non_zero))
  }
  
  plot(t_values, mu, pch=20, main="mean distance against duration", xlab="duration/s", ylab="distance/km")
  ## mu vs duration
  knots <- c(min_dur, min_dur, min_dur, min_dur, 25, 60,  max_dur, max_dur, max_dur, max_dur)
  bx <- splineDesign(knots, t_values)
  lmod_spline_mean_vs_duration <- lm(mu ~ bx -1)
  matplot(t_values, cbind(mu, lmod_spline_mean_vs_duration$fit), type="pl", main = "fit test", ylab="mean distance", xlab="duration/s", pch=20, lty=1, col=c(1, "red"))
  sink(file="lmod_spline_mean_vs_duration.txt")
  summary(lmod_spline_mean_vs_duration)
  sink(file=NULL)
  
  
  plot(t_values, var, pch=20, main="variance of distance against duration", xlab="duration/s", ylab="distance/km")
  
  ## var vs duration
  knots <- c(min_dur, min_dur, min_dur, min_dur, 60, 100, max_dur, max_dur, max_dur, max_dur)
  bx <- splineDesign(knots, t_values)
  lmod_spline_var_vs_duration <- lm(var ~ bx -1)
  matplot(t_values, cbind(var, lmod_spline_var_vs_duration$fit), type="pl", main = "fit test", ylab="variance of distance", xlab="duration/s", pch=20, lty=1, col=c(1, "red"))
  sink(file="lmod_spline_var_vs_duration.txt")
  summary(lmod_spline_var_vs_duration)
  sink(file=NULL)
}





### fitting normal distribution to non-zero distances(included in the report)  ###
{
  if(!require('hexbin')) {
    install.packages('hexbin')
    library('hexbin')
  }
  
  fit_normal_dist <- function(duration) {
    t <- duration
    dist_con_dur <- df_one_side_trimmed[(df_one_side_trimmed$duration == t) & (df_one_side_trimmed$dist > 0),]$dist
    hist(dist_con_dur, breaks=100, main=cbind("travel distances(non-zero) for travels with duration", toString(t), "min"), xlab="distance/km")
    z_dist <- (dist_con_dur - mean(dist_con_dur)) / sd(dist_con_dur)
    qqnorm(z_dist, pch=20)
    abline(0, 1, col="red")
  }
  
  fit_normal_dist(10)
  fit_normal_dist(12)
  fit_normal_dist(15)
  fit_normal_dist(20)
  fit_normal_dist(30)
  fit_normal_dist(40)
  fit_normal_dist(60)
  ## for 60min, there seems to be an exponential distribution involved 
  
  ## ks_score against duration
  ks_scores <- c()
  ks_p_vals <- c()
  for (t in 1:max_dur) {
    dist_con_dur <- df_one_side_trimmed[(df_one_side_trimmed$duration == t) & (df_one_side_trimmed$dist > 0),]$dist
    ksresult <- ks.test(dist_con_dur, "pnorm", mean=mean(dist_con_dur), sd = sd(dist_con_dur))
    ks_scores <- c(ks_scores, ksresult$statistic)
    ks_p_vals <- c(ks_p_vals, ksresult$p.value)
  }
  dur_values = min_dur:max_dur
  plot(dur_values, ks_scores, pch=20, main="ks scores against duration", xlab="duration/min", ylab="ks score")
  plot(dur_values, log(ks_p_vals), pch=20, main="ks test p values against duration", xlab="duration/min", ylab="log(p value)")
  abline(0.05, 0, col="red")
  
}

##############################################
### try mixture models, using EM algorithm(included in the report) ###
##############################################
{
  ## calculate MLE for mixture of two distributions where one is gamma (the other one to be decided)
  ## given the estimate in teh last round alpha, the data x 
  ## and estimates of conditional expectations z where z_i = E(xi belongs to this gamma distribution|xi)
  ## Warning: the new_beta returned is the scale, not rate
  MLE_gamma <- function(alpha, z, x) {
    if (length(z) != length(x)) {
      stop("Length of z, x should be the same")
    }
    
    if (alpha <= 0) {
      stop("alpha should be positive")
    }
    
    new_alpha <- alpha - (log(alpha) - digamma(alpha) - log(sum(z * x) / sum(z)) + sum(z * log(x)) / sum(z)) / (1/alpha - trigamma(alpha))
    new_beta <- sum(z * x) / (new_alpha * sum(z))
    
    c(new_alpha, new_beta)
  }
  
  ## calculate MLE for mixture of two distributions where one is Weibull (the other one to be decided)
  ## given the estimate in the last round beta, the data x 
  ## and estimates of conditional expectations z where z_i = E(xi belongs to this gamma distribution|xi)
  MLE_Weibull <- function(beta, z, x) {
    if (length(z) != length(x)) {
      stop("Length of z, x should be the same")
    }
    
    if (beta <= 0) {
      stop("alpha should be positive")
    }
    
    if (length(beta) > 1) {
      stop("The first argument is point estimation of beta")
    }
    
    A <- sum(z * log(x)) / sum(z)
    B <- sum(z * (x^beta))
    C <- sum(z * (x^beta) * log(x))
    D <- sum(z * (x^beta) * (log(x))^2)
    
    new_beta <- beta + (A + (1/beta) - (C / B)) / ((1/beta)^2 + (B * D - C^2) / B^2)
    new_alpha <- (sum(z * (x^new_beta)) / sum(z))^(1 / new_beta)
    
    c(new_alpha, new_beta)
  }
  
  MLE_exponential <- function(z, x) {
    if (length(z) != length(x)) {
      stop("Length of z, x should be the same")
    }
    ## estimator for lambda 
    sum(z * x) / sum(z)
  }
  
  MLE_normal <- function(z, x) {
    if (length(z) != length(x)) {
      stop("Length of z, x should be the same")
    }
    
    ## effective number of points that follows normal distribution
    N <- sum(z)
    new_mu <- sum(z * x) / N
    new_var <- sum(z * (x-new_mu)^2) / N
    
    c(new_mu, new_var)
  }
  
  ## enter the data x, initial values of three parameters lambda, alpha, beta(scale) 
  ## and initial weight pi is the proportion of data that belongs to the first distribution
  ## tol sets the threshold to stop iterations (convergence is L_2 norm of the parameters)
  ## stop_method must be one of "parameter" and "loglik"
  EM_exponential_gamma <- function(x, lambda, alpha, beta, pi, tol = 1e-5, stop_method = "parameter") {
    ## diff measures difference between estimates of two iterations (using Euler distance)
    diff <- 1 
    n <- length(x)
    ## r is number of iterations 
    r <- 0
    
    ## store previous estimates 
    prev <- c(lambda, alpha, beta, pi)
    
    if(!(is.element(stop_method, c("parameter", "loglik")))) {
      stop("stop_method must be one of 'parameter' and 'loglik'")
    }
    
    if ((pi > 1) | (pi < 0)) {
      stop("initial weight pi should be a probability! (between 0, 1)")
    }
    
    ## used to compute log likelihood in EM iterations
    ## par is a vector c(lambda, alpha, beta, pi)
    log_lik <- function(x, par) {
      sum(log(par[4] * dexp(x, rate = par[1]) + (1 - par[4]) * dgamma(x, par[2], scale=par[3])))
    }
    
    ll <- c()
    
    while (diff > tol) {
      r = r + 1
      ## E-step
      # vector of probabilities of each data point belonging to the first distribution(exponential)
      z_1 <- prev[4] * dexp(x, rate=prev[1]) / (prev[4] * dexp(x, rate=prev[1]) + (1-prev[4]) * dgamma(x, prev[2], scale=prev[3]))
      # conditional probability on the second distribution
      z_2 <- 1 - z_1
      
      ## M-step
      new_pi <- mean(z_1)
      new_lambda <- MLE_exponential(z_1, x)
      gamma_par <- MLE_gamma(prev[2], z_2, x)
      new_alpha <- gamma_par[1]
      new_beta <- gamma_par[2]
      
      current <- c(new_lambda, new_alpha, new_beta, new_pi)
      
      ll <- c(ll, log_lik(x, current))
      
      ## update difference and prev
      if (stop_method == "parameter") {
        diff <- sum((current - prev)^2)
      } else {
        diff <- abs(diff(tail(ll, 2)))
      }
      
      prev <- current
    }
    
    print(paste("number of iterations: ", toString(r), sep=""))
    list("par" = current, "log_lik" = ll, "it_nums" = r)
  }
  
  ## draw histogram to test fitness of exponential-gamma mixture models 
  ## par must be the result from the function EM_exponential_gamma
  fitness_test_exp_gamma <- function(t, x, par, breaks=50, plot=TRUE, test_set = c()) {
    
    if (plot) {
      h <- hist(x, plot=FALSE, breaks=breaks)
      plot(h, col="grey", xlab="distance/km", main=paste("histogram of travel distances(duration = ", toString(t), "min)", sep=""), ylim=c(0, ceiling(max(h$counts) * 1.15)))
      xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
      lines(x = xlines, y=dexp(xlines, rate=par[1]) *length(x)*diff(h$breaks)[1] * par[4],col="red", lty=2)
      lines(x = xlines, y=dgamma(xlines, par[2], scale=par[3]) *length(x)*diff(h$breaks)[1] * (1-par[4]), col="blue", lty=2)
      lines(x = xlines, y=dgamma(xlines, par[2], scale=par[3]) *length(x)*diff(h$breaks)[1] * (1-par[4]) + dexp(xlines, rate=par[1]) *length(x)*diff(h$breaks)[1] * par[4])
      legend(x = "topleft", legend = c("exponential", "gamma"), lty = c(2, 2), col = c("red", "blue"), lwd=2, cex=0.8)
    }
    
    ## perform ks test
    n <- length(x)
    n_first <- floor(par[4] * n)
    sample <- c(rexp(n_first, par[1]), rgamma(n - n_first, par[2], scale=par[3]))
    ksresult <- ks.test(x, sample)

    
    if (length(test_set) > 0) {
      n <- length(test_set)
      n_first <- floor(par[5] * n)
      sample <- c(rexp(n_first, par[1]), rgamma(n - n_first, par[2], scale=par[3]))
      test_set_ksresult <- ks.test(test_set, sample)
      return(list("ks" = ksresult$statistic, "p_val" = ksresult$p.value, "test_ks" = test_set_ksresult$statistic))
    }
    
    return(list("ks" = ksresult$statistic, "p_val" = ksresult$p.value))
  }
  
  ## enter the data x, initial values of three parameters lambda, alpha, beta
  ## and initial weight pi is the proportion of data that belongs to the first distribution
  ## tol sets the threshold to stop iterations (convergence is L_2 norm of the parameters)
  EM_exponential_weibull <- function(x, lambda, alpha, beta, pi, tol = 1e-5, max_it = 1000, stop_method = "parameter") {
    if(!(is.element(stop_method, c("parameter", "loglik")))) {
      stop("stop_method must be one of 'parameter' and 'loglik'")
    }
    
    if ((pi > 1) | (pi < 0)) {
      stop("initial weight pi should be a probability! (between 0, 1)")
    }
    
    ## diff measures difference between estimates of two iterations (using Euler distance)
    diff <- 1 
    n <- length(x)
    ## r is number of iterations 
    r <- 0
    
    ## store previous estimates 
    prev <- c(lambda, alpha, beta, pi)
    
    log_lik <- function(x, par) {
      sum(log( par[4] * dexp(x, rate = par[1]) + (1-par[4]) * dweibull(x, par[2], par[3]) ))
    }
    
    
    ll <- c()
    
    while (diff > tol) {
      if (r > max_it) {
        print(paste("maximum iterations:", toString(max_it), "is reached", sep=" "))
        return(list("par" = current, "log_lik" = ll, "it_nums" = r))
      }
      r = r + 1
      ## E-step
      # vector of probabilities of each data point belonging to the first distribution(exponential)
      z_1 <- prev[4] * dexp(x, rate=prev[1]) / (prev[4] * dexp(x, rate=prev[1]) + (1-prev[4]) * dweibull(x, prev[2], prev[3]))
      # conditional probability on the second distribution
      z_2 <- 1 - z_1
      
      ## M-step
      new_pi <- mean(z_1)
      new_lambda <- MLE_exponential(z_1, x)
      weibull_par <- MLE_Weibull(prev[3], z_2, x)
      new_alpha <- weibull_par[1]
      new_beta <- weibull_par[2]
      
      current <- c(new_lambda, new_alpha, new_beta, new_pi)
      
      ll <- c(ll, log_lik(x, current))
      
      ## update difference and prev
      if (stop_method == "parameter") {
        diff <- sum((current - prev)^2)
      } else {
        diff <- abs(diff(tail(ll, 2)))
      }
      
      prev <- current
    }
    
    print(paste("number of iterations: ", toString(r), sep=""))
    list("par" = current, "log_lik" = ll, "it_nums" = r)
  }
  
  ## draw histogram to test fitness of exponential-weibull mixture models 
  ## par must be the result from the function EM_exponential_weibull
  fitness_test_exp_weib <- function(x, par, breaks=50) {
    h <- hist(x, plot=FALSE, breaks=breaks)
    plot(h, col="grey")
    xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
    lines(x = xlines, y=dexp(xlines, rate=par[1]) *length(x)*diff(h$breaks)[1] * par[4], col="red", lty=2)
    lines(x = xlines, y=dweibull(xlines, par[2], par[3]) *length(x)*diff(h$breaks)[1] * (1-par[4]), col="blue", lty=2)
    lines(x = xlines, y=dweibull(xlines, par[2], par[3]) *length(x)*diff(h$breaks)[1] * (1-par[4]) + dexp(xlines, rate=par[1]) *length(x)*diff(h$breaks)[1] * par[4], col="blue", lty=2)
    
    ## perform ks test
    n <- length(x)
    n_first <- floor(par[4] * n)
    sample <- c(rexp(n_first, par[1]), rweibull(n - n_first, par[2], par[3]))
    ksresult <- ks.test(x, sample)
    list("ks" = ksresult$statistic, "p_val" = ksresult$p.value)
  }
  
  ## EM on a mixture of a gamma and a Weibull distribution
  EM_gamma_weibull <- function(x, alpha1, beta1, alpha2, beta2, pi, tol = 1e-5, max_it = 1000, stop_method = "parameter") {
    if(!(is.element(stop_method, c("parameter", "loglik")))) {
      stop("stop_method must be one of 'parameter' and 'loglik'")
    }
    
    if ((pi > 1) | (pi < 0)) {
      stop("initial weight pi should be a probability! (between 0, 1)")
    }
    
    ## diff measures difference between estimates of two iterations (using Euler distance)
    diff <- 1 
    n <- length(x)
    ## r is number of iterations 
    r <- 0
    
    ## store previous estimates 
    prev <- c(alpha1, beta1, alpha2, beta2, pi)
    
    log_lik <- function(x, par) {
      sum(log( par[5] * dgamma(x, par[1], scale=par[2]) + (1-par[5]) * dweibull(x, par[3], par[4]) ))
    }
    
    
    ll <- c()
    
    while (diff > tol) {
      if (r > max_it) {
        print(paste("maximum iterations:", toString(max_it), "is reached", sep=" "))
        return(list("par" = current, "log_lik" = ll, "it_nums" = r))
      }
      r = r + 1
      ## E-step
      # vector of probabilities of each data point belonging to the first distribution(exponential)
      z_1 <- prev[5] * dgamma(x, prev[1], scale=prev[2]) / (prev[5] * dgamma(x, prev[1], scale=prev[2]) + (1-prev[5]) * dweibull(x, prev[3], prev[4]))
      # conditional probability on the second distribution
      z_2 <- 1 - z_1
      
      ## M-step
      new_pi <- mean(z_1)
      gamma_par <- MLE_gamma(prev[1], z_1, x)
      weibull_par <- MLE_Weibull(prev[4], z_2, x)
      new_alpha1 <- gamma_par[1]
      new_beta1 <- gamma_par[2]  
      new_alpha2 <- weibull_par[1]
      new_beta2 <- weibull_par[2]
      
      current <- c(new_alpha1, new_beta1, new_alpha2, new_beta2, new_pi)
      
      ll <- c(ll, log_lik(x, current))
      
      ## update difference and prev
      if (stop_method == "parameter") {
        diff <- sum((current - prev)^2)
      } else {
        diff <- abs(diff(tail(ll, 2)))
      }
      
      prev <- current
    }
    
    print(paste("number of iterations: ", toString(r), sep=""))
    list("par" = current, "log_lik" = ll, "it_nums" = r)
  }
  
  ## draw histogram to test fitness of gamma-weibull mixture models 
  ## par must be the result from the function EM_gamma_weibull
  fitness_test_gamma_weibull <- function(t, x, par, breaks=50, plot=TRUE, test_set = c()) {
    
    if (plot) {
      h <- hist(x, plot=FALSE, breaks=breaks)
      plot(h, col="grey", xlab="distance/km", main=paste("histogram of travel distances(duration = ", toString(t), "min)", sep=""), ylim=c(0, ceiling(max(h$counts) * 1.15)))
      xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
      lines(x = xlines, y=dgamma(xlines, par[1], scale=par[2]) *length(x)*diff(h$breaks)[1] * par[5], col="red", lty=2)
      lines(x = xlines, y=dweibull(xlines, par[3], par[4]) *length(x)*diff(h$breaks)[1] * (1-par[5]), col="blue", lty=2)
      lines(x = xlines, y=dweibull(xlines, par[3], par[4]) *length(x)*diff(h$breaks)[1] * (1-par[5]) + 
              dgamma(xlines, par[1], scale=par[2]) *length(x)*diff(h$breaks)[1] * par[5])
      legend(x = "topleft", legend = c("gamma", "weibull"), lty = c(2, 2), col = c("red", "blue"), lwd=2, cex=0.8)
    }
    
    ## perform ks test
    n <- length(x)
    n_first <- floor(par[5] * n)
    sample <- c(rgamma(n_first, par[1], scale=par[2]), rweibull(n - n_first, par[3], par[4]))
    ksresult <- ks.test(x, sample)
    list("ks" = ksresult$statistic, "p_val" = ksresult$p.value)
  }
  
  ## EM on a mixture of a gamma and a normal distribution
  EM_gamma_normal <- function(x, alpha, beta, mu, var, pi, tol = 1e-5, max_it = 1000, stop_method = "parameter") {
    if(!(is.element(stop_method, c("parameter", "loglik")))) {
      stop("stop_method must be one of 'parameter' and 'loglik'")
    }
    
    if ((pi > 1) | (pi < 0)) {
      stop("initial weight pi should be a probability! (between 0, 1)")
    }
    
    ## diff measures difference between estimates of two iterations (using Euler distance)
    diff <- 1 
    n <- length(x)
    ## r is number of iterations 
    r <- 0
    
    ## store previous estimates 
    prev <- c(alpha, beta, mu, var, pi)
    
    log_lik <- function(x, par) {
      sum(log(par[5] * dgamma(x, par[1], scale=par[2]) + (1-par[5]) * dnorm(x, mean=par[3], sd=sqrt(par[4]))))
    }
    
    
    ll <- c()
    
    while (diff > tol) {
      if (r > max_it) {
        print(paste("maximum iterations:", toString(max_it), "is reached", sep=" "))
        return(list("par" = current, "log_lik" = ll, "it_nums" = r, "max_it_reached" = TRUE))
      }
      r = r + 1
      ## E-step
      # vector of probabilities of each data point belonging to the first distribution(exponential)
      z_1 <- prev[5] * dgamma(x, prev[1], scale=prev[2]) / (prev[5] * dgamma(x, prev[1], scale=prev[2]) + (1-prev[5]) * dnorm(x, mean=prev[3], sd = sqrt(prev[4])))
      # conditional probability on the second distribution
      z_2 <- 1 - z_1
      
      ## M-step
      new_pi <- mean(z_1)
      gamma_par <- MLE_gamma(prev[1], z_1, x)
      normal_par <- MLE_normal(z_2, x)
      new_alpha <- gamma_par[1]
      new_beta <- gamma_par[2]  
      new_mu <- normal_par[1]
      new_var <- normal_par[2]
      
      current <- c(new_alpha, new_beta, new_mu, new_var, new_pi)
      
      ll <- c(ll, log_lik(x, current))
      
      ## update difference and prev
      if (stop_method == "parameter") {
        diff <- sum((current - prev)^2)
      } else {
        diff <- abs(diff(tail(ll, 2)))
      }
      
      prev <- current
    }
  
    list("par" = current, "log_lik" = ll, "it_nums" = r, "max_it_reached" = FALSE)
  }
  
  ## draw histogram to test fitness of gamma-normal mixture models 
  ## par must be the result from the function EM_gamma_normal
  ## test set is used to detect over-fitting, a test ks will be returned if test_set is not empty
  fitness_test_gamma_norm <- function(t, x, par, breaks=50, plot=TRUE, test_set = c()) {
    if (plot) {
      h <- hist(x, plot=FALSE, breaks=breaks)
      plot(h, col="grey", xlab="distance/km", main=paste("histogram of travel distances(duration = ", toString(t), "min)", sep=""), ylim=c(0, ceiling(max(h$counts) * 1.15)))
      xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
      lines(x = xlines, y= dgamma(xlines, par[1], scale=par[2]) *length(x)*diff(h$breaks)[1] * par[5], col="red", lty=2)
      lines(x = xlines, y= dnorm(xlines, mean = par[3], sd = sqrt(par[4])) *length(x)*diff(h$breaks)[1] * (1-par[5]), col="blue", lty=2)
      lines(x = xlines, y= dnorm(xlines, mean = par[3], sd = sqrt(par[4])) *length(x)*diff(h$breaks)[1] * (1-par[5]) + 
              dgamma(xlines, par[1], scale=par[2]) *length(x)*diff(h$breaks)[1] * par[5])
      legend(x = "topleft", legend = c("gamma", "normal"), lty = c(2, 2), col = c("red", "blue"), lwd=2, cex=0.8)
    }
    
    ## perform ks test
    n <- length(x)
    n_first <- floor(par[5] * n)
    sample <- c(rgamma(n_first, par[1], scale=par[2]), rnorm(n - n_first, mean=par[3], sd=sqrt(par[4])))
    ksresult <- ks.test(x, sample)
    
    
    if (length(test_set) > 0) {
      n <- length(test_set)
      n_first <- floor(par[5] * n)
      sample <- c(rgamma(n_first, par[1], scale=par[2]), rnorm(n - n_first, mean=par[3], sd=sqrt(par[4])))
      test_set_ksresult <- ks.test(test_set, sample)
      return(list("ks" = ksresult$statistic, "p_val" = ksresult$p.value, "test_ks" = test_set_ksresult$statistic))
    }
    
    return(list("ks" = ksresult$statistic, "p_val" = ksresult$p.value))
  }
  
  
  ## test on duration = 40 min distances
  test_set <- df_one_side_trimmed[(df_one_side_trimmed$duration == 40) & (df_one_side_trimmed$dist > 0),]$dist
  
  ## exponential-gamma
  test_EM_eg <- EM_exponential_gamma(test_set, 1, 1, 5, 0.5)
  test_EM_eg$par
  plot(1:test_EM_eg$it_nums, test_EM_eg$log_lik, xlab="iteraion", ylab="log likelihood")
  test_fitness_eg <- fitness_test_exp_gamma(40, test_set, test_EM_eg$par, breaks=200)
  test_fitness_eg
  
  ## exponential-Weibull
  test_EM_ew <- EM_exponential_weibull(test_set, 1, 1, 1, 0.5)
  ## fails
  
  ## gamma-Weibull
  test_EM_gw <- EM_gamma_weibull(test_set, 1, 1, 5, 5, 0.5)
  test_EM_gw$par
  plot(1:test_EM_gw$it_nums, test_EM_gw$log_lik)
  test_fitness_gw <- fitness_test_gamma_weibull(40, test_set, test_EM_gw$par, breaks=200)
  test_fitness_gw
  
  ## try mixture of gamma (2-components)
  ## initial values chosen based on moment matching and Finch's method
  {
    pi_0 <- runif(1, min=0, max=1)
    n <- length(test_set)
    n_first <- floor(n * pi_0)
    mu_1_0 <- mean(test_set) - sd(test_set) / 2
    mu_2_0 <- mean(test_set) + sd(test_set) / 2
    var_0 <- (sum((test_set[1:n_first] - mu_1_0)^2) + sum((test_set[(n_first+1):n] - mu_2_0)^2)) / (n - 2)
    alpha_0 <- c(mu_1_0, mu_2_0)^2 / var_0
    beta_0 <- c(mu_1_0, mu_2_0) / var_0
  }
  
  ## gamma-gamma
  library(mixtools)
  test_EM_gg <- gammamixEM(test_set, lambda = c(pi_0, 1-pi_0), alpha=alpha_0, beta=beta_0, k=2, maxit=1000, maxrestarts = 20)
  ## fails, not mixture of gamma
  
  ## normal-normal
  test_EM_nn <- normalmixEM(test_set, lambda = c(pi_0, 1-pi_0), mu = c(mu_1_0, mu_2_0), sigma = c(sqrt(var_0), sqrt(var_0)))
  plot(test_EM_nn, loglik = TRUE)
  length(test_EM_nn$all.loglik)
  
  ## lambda is the mixing proportion
  fitness_test_norm_norm <- function(t, x, lambda, mu, sigma, breaks=50, plot=TRUE, test_set = c()) {
    if (plot) {
      h <- hist(x, plot=FALSE, breaks=breaks)
      plot(h, col="grey", xlab="distance/km", main=paste("histogram of travel distances(duration = ", toString(t), "min)", sep=""), ylim=c(0, ceiling(max(h$counts) * 1.3)))
      xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
      lines(x = xlines, y= dnorm(xlines, mu[1], sigma[1]) *length(x)*diff(h$breaks)[1] * lambda[1], col="red", lty=2)
      lines(x = xlines, y= dnorm(xlines, mean = mu[2], sd = sigma[2]) *length(x)*diff(h$breaks)[1] * lambda[2], col="blue", lty=2)
      lines(x = xlines, y= dnorm(xlines, mu[1], sigma[1]) *length(x)*diff(h$breaks)[1] * lambda[1] + 
              dnorm(xlines, mean = mu[2], sd = sigma[2]) *length(x)*diff(h$breaks)[1] * lambda[2])
      legend(x = "topleft", legend = c("first normal", "second normal"), lty = c(2, 2), col = c("red", "blue"), lwd=2, cex=0.8)
    }
    
    ## perform ks test
    n <- length(x)
    n_first <- floor(lambda[1] * n)
    sample <- c(rnorm(n_first, mu[1], sigma[1]), rnorm(n - n_first, mu[2], sigma[2]))
    ksresult <- ks.test(x, sample)
    
    
    if (length(test_set) > 0) {
      n <- length(test_set)
      n_first <- floor(lambda[1] * n)
      sample <- c(rnorm(n_first, mu[1], sigma[1]), rnorm(n - n_first, mu[2], sigma[2]))
      test_set_ksresult <- ks.test(test_set, sample)
      return(list("ks" = ksresult$statistic, "p_val" = ksresult$p.value, "test_ks" = test_set_ksresult$statistic))
    }
    
    return(list("ks" = ksresult$statistic, "p_val" = ksresult$p.value))
  }
  
  fitness_test_norm_norm(40, test_set, test_EM_nn$lambda, test_EM_nn$mu, test_EM_nn$sigma, breaks=200)
  
  ## gamma-normal
  test_EM_gn <- EM_gamma_normal(test_set, mu_1_0^2 / var_0, mu_1_0 / var_0, mu_2_0, var_0, 0.5, max_it=150)
  test_EM_gn$par
  test_EM_gn$it_nums
  plot(1:test_EM_gn$it_nums, test_EM_gn$log_lik, ylab="log likelihood", xlab="iteration")
  fitness_test_gamma_norm(40, test_set, test_EM_5$par, breaks=200)
}



## Based on the above test, use gamma-normal to study (included in the report) 
{
  
  ## sample_size returned is the size of training set!! 
  ## initial proportion is the way of choosing initial mixing proportion pi_0
  mixture_gamma_normal <- function(t, initial_proportion = "even", max_it = 200, max_restart = 10, plot=FALSE, breaks = 100, tol = 1e-3, train_proportion = 0.8) {
    data <- df_one_side_trimmed[(df_one_side_trimmed$duration == t) & (df_one_side_trimmed$dist > 0),]$dist
    n <- length(data)
    train_length <- floor(n * train_proportion)
    train_set <- data[1:train_length]
    test_set <- data[(train_length+1):n]
    
    restart <- FALSE
    restart_num <- -1
    
    repeat {
      if (restart_num >= max_restart) {
        print(paste("maximum restart number:", toString(max_restart), "reached", sep = " "))
        break
      }
      
      restart_num <- restart_num + 1
      #choose initial parameter values 
      if (initial_proportion == "even") {
        pi_0 <- 0.5
      } else if (initial_proportion == "random") {
        pi_0 <- runif(1)
      } else if (initial_proportion == "GMM") {
        ## use mixing proportion from gaussian mixture model 
        n <- length(train_set)
        n_first <- floor(n * 0.5)
        mu_1_0 <- mean(train_set) - sd(train_set) / 2
        mu_2_0 <- mean(train_set) + sd(train_set) / 2
        pi_0 <- normalmixEM(test_set, lambda = c(pi_0, 1-pi_0), mu = c(mu_1_0, mu_2_0), sigma = c(sqrt(var_0), sqrt(var_0)))$lambda[1]
      } else {
        stop("initial_proportion must be one of: even, random, GMM")
      }
      n <- length(train_set)
      n_first <- floor(n * pi_0)
      mu_1_0 <- mean(train_set) - sd(train_set) / 2
      mu_2_0 <- mean(train_set) + sd(train_set) / 2
      var_0 <- (sum((train_set[1:n_first] - mu_1_0)^2) + sum((train_set[(n_first+1):n] - mu_2_0)^2)) / (n - 2)
      
      EM_result <- EM_gamma_normal(train_set, mu_1_0^2 / var_0, mu_1_0 / var_0, mu_2_0, var_0, pi_0, max_it= max_it, tol=tol)
      if ((!EM_result$max_it_reached) || (initial_proportion == "even")) {
        break
      }
    }
    
    fitness_result <- fitness_test_gamma_norm(t, train_set, EM_result$par, plot = plot, breaks = breaks, test_set = test_set)
    list("par" = EM_result$par, "log_lik" = tail(EM_result$log_lik, 1), "ks" = fitness_result$ks, "test_ks" = fitness_result$test_ks, "sample_size" = train_length, "it_nums" = EM_result$it_nums, "restarts" = restart_num)
  }
  
  mixture_gamma_normal(15, initial_proportion = "even", plot=TRUE, max_it=150)
  mixture_gamma_normal(15, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(15, initial_proportion = "random", plot=TRUE, max_it=150)
  mixture_gamma_normal(18, plot=TRUE, max_it=150)
  mixture_gamma_normal(18, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(18, initial_proportion = "random", plot=TRUE, max_it=150)
  mixture_gamma_normal(24, plot=TRUE, max_it=150)
  mixture_gamma_normal(24, initial_proportion = "random", plot=TRUE, max_it=150)
  mixture_gamma_normal(24, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(25, initial_proportion = "random", plot=TRUE, max_it=150)
  mixture_gamma_normal(25, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(30, plot=TRUE, max_it=150)
  mixture_gamma_normal(30, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(30, initial_proportion = "random", plot=TRUE, max_it=150)
  mixture_gamma_normal(40, plot=TRUE, max_it=150)
  mixture_gamma_normal(40, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(40, initial_proportion = "random", plot=TRUE, max_it=150)
  mixture_gamma_normal(120, plot=TRUE, max_it=150)
  mixture_gamma_normal(120, initial_proportion = "GMM", plot=TRUE, max_it=150)
  mixture_gamma_normal(120, initial_proportion = "random", plot=TRUE, max_it=150)
  
  
  
  
  ### running the programme on the whole dataset ### 
  ## initialise
  {
    ks_gamma_normal <- c()
    test_ks_gamma_normal <- c()
    par_gamma_normal <- matrix(0, nrow = max_dur, ncol = 5, byrow=TRUE)
    colnames(par_gamma_normal) <- c("shape(gamma)", "scale(gamma)", "mean(normal)", "variance(normal)", "proportion of gamma")
    rownames(par_gamma_normal) <- 1:max_dur
    it_nums_gamma_normal <- c()
    restarts_gamma_normal <- c()
    sample_sizes <- c()
  }
  
  ## must run the previous part before this part
  ## main program: takes about 4-5min
  for (t in 1:max_dur) { 
    ## there is no sign of mixture for duration < 15min, so let initial mixing proportion be random
    if (t < 15){
      ## sample size too large, reduce max_restarts
      ## actually accurate results will still be produced 
      mix_result <- mixture_gamma_normal(t, initial_proportion = "random", plot=FALSE, max_it=150, max_restart = 4)
    } else {
      mix_result <- mixture_gamma_normal(t, initial_proportion = "random", plot=FALSE, max_it=150, max_restart = 15)
    }
    
    ks_gamma_normal <- c(ks_gamma_normal, mix_result$ks)
    test_ks_gamma_normal <- c(test_ks_gamma_normal, mix_result$test_ks)
    par_gamma_normal[t, ] <- mix_result$par
    it_nums_gamma_normal <- c(it_nums_gamma_normal, mix_result$it_nums)
    sample_sizes <- c(sample_sizes, mix_result$sample_size)
    restarts_gamma_normal <- c(restarts_gamma_normal, mix_result$restarts)
    print(paste("duration = ", toString(t), "finished", sep=" "))
  }
  
  
  
  
  ## show those cases where max_restarts reached
  check_fitness_gamma_normal <- function(t) {
    data <- df_one_side_trimmed[(df_one_side_trimmed$duration == t) & (df_one_side_trimmed$dist > 0),]$dist
    fitness_result <- fitness_test_gamma_norm(t, data, par_gamma_normal[t, ], plot=TRUE, breaks=100)
    print(paste("ks : ", fitness_result$ks, sep=""))
  }
  
  for (t in which(it_nums_gamma_normal %in% c(151, 201))) {
    check_fitness_gamma_normal(t)
  }
  
  
  plot(sample_sizes, ks_gamma_normal, pch=20, xlab="sample size", ylab="ks-score")
  plot(1:max_dur, sample_sizes, pch=20, xlab="duration/min", ylab="sample size", xaxp = c(0, 180, 6))
  plot(1:max_dur, ks_gamma_normal, pch=19, xlab="duration/min", ylab="ks-scores", col="black", ylim=c(0, max(test_ks_gamma_normal) + 0.05))
  points(1:max_dur, test_ks_gamma_normal, pch=20, col="red")
  legend(x = "topleft", legend = c("train set", "test set"), col = c("black", "red"), pch=c(19, 20))
  
  plot(1:max_dur, par_gamma_normal[,5], xlab="duration/min", ylab="mixing proportion of gamma", xaxp = c(0, 180, 6), pch=20)
  sink(file="mixing_proportions_EM_1.txt")
  head(par_gamma_normal[,5], 25)
  sink(file=NULL)
  ## it seems first 13 (except duration = 1min) data are almost purely normal. So we ignore the gamma part of them
  check_fitness_gamma_normal(14)
  check_fitness_gamma_normal(16)
  check_fitness_gamma_normal(18)
  check_fitness_gamma_normal(19)
  check_fitness_gamma_normal(21)
  check_fitness_gamma_normal(26)
  check_fitness_gamma_normal(29)
  check_fitness_gamma_normal(30)
  check_fitness_gamma_normal(31)
  check_fitness_gamma_normal(32)
  
  plot(1:max_dur, par_gamma_normal[,3], xlab="duration/min", ylab="mean of normal", xaxp = c(0, 180, 6), pch=20)
  plot(1:max_dur, par_gamma_normal[,4], xlab="duration/min", ylab="variance of normal", xaxp = c(0, 180, 6), pch=20)
  plot(19:max_dur, par_gamma_normal[19:max_dur,1], xlab="duration/min", ylab="shape of gamma", xaxp = c(0, 180, 6), pch=20)
  plot(15:max_dur, par_gamma_normal[15:max_dur,2], xlab="duration/min", ylab="scale of gamma", xaxp = c(0, 180, 6), pch=20)
  plot(1:max_dur, par_gamma_normal[,5], xlab="duration/min", ylab="mixing proportion of gamma", xaxp = c(0, 180, 6), pch=20)
  
  
  ### same program, but start with equal mixing proportion
  ## initialise
  {
    ks_gamma_normal_2 <- c()
    test_ks_gamma_normal_2 <- c()
    par_gamma_normal_2 <- matrix(0, nrow = max_dur, ncol = 5, byrow=TRUE)
    colnames(par_gamma_normal_2) <- c("shape(gamma)", "scale(gamma)", "mean(normal)", "variance(normal)", "proportion of gamma")
    rownames(par_gamma_normal_2) <- 1:max_dur
    it_nums_gamma_normal_2 <- c()
    restarts_gamma_normal_2 <- c()
    sample_sizes_2 <- c()
  }
  
  ## must run the previous part before this part
  ## main program: takes about 4-5min
  for (t in 1:max_dur) { 
    ## there is no sign of mixture for duration < 15min, so let initial mixing proportion be random
    if (t < 15) {
      mix_result <- mixture_gamma_normal(t, plot=FALSE, max_it=150)
    } else {
      mix_result <- mixture_gamma_normal(t, plot=FALSE, max_it=200)
    }
    ks_gamma_normal_2 <- c(ks_gamma_normal_2, mix_result$ks)
    test_ks_gamma_normal_2 <- c(test_ks_gamma_normal_2, mix_result$test_ks)
    par_gamma_normal_2[t, ] <- mix_result$par
    it_nums_gamma_normal_2 <- c(it_nums_gamma_normal_2, mix_result$it_nums)
    sample_sizes_2 <- c(sample_sizes_2, mix_result$sample_size)
    restarts_gamma_normal_2 <- c(restarts_gamma_normal_2, mix_result$restarts)
    print(paste("duration = ", toString(t), "finished", sep=" "))
  }
  
  df_gamma_normal_2 <- data.frame("duration" = 1:max_dur,"size" = sample_sizes_2, "ks" = ks_gamma_normal_2, "test_ks" = test_ks_gamma_normal_2, "mean_gamma" = par_gamma_normal_2[, 1] * par_gamma_normal_2[, 2], "var_gamma" = par_gamma_normal_2[, 1] * par_gamma_normal_2[, 2]^2, 
                                  "mean_normal" = par_gamma_normal_2[, 3], "var_normal" = par_gamma_normal_2[, 4], "proportion" = par_gamma_normal_2[, 5])
  
  plot(df_gamma_normal_2$duration, df_gamma_normal_2$proportion, xlab="duration/min", ylab="mixing proportion of gamma", xaxp = c(0, 180, 6), pch=20)
  sink(file="mixing_proportions_EM_2.txt")
  head(par_gamma_normal_2[, 5], 25)
  sink(file=NULL)
  
  check_fitness_gamma_normal_2 <- function(t) {
    data <- df_one_side_trimmed[(df_one_side_trimmed$duration == t) & (df_one_side_trimmed$dist > 0),]$dist
    fitness_result <- fitness_test_gamma_norm(t, data, par_gamma_normal_2[t, ], plot=TRUE, breaks=100)
    print(paste("ks : ", fitness_result$ks, sep=""))
  }
  
  plot(df_gamma_normal_2$size, df_gamma_normal_2$ks, pch=20, xlab="sample size", ylab="ks-score")
  
  plot(df_gamma_normal_2$duration, df_gamma_normal_2$ks, pch=19, xlab="duration/min", ylab="ks-scores", col="black", ylim=c(0, max(test_ks_gamma_normal_2) + 0.05), xaxp = c(0, 180, 6))
  points(df_gamma_normal_2$duration, df_gamma_normal_2$test_ks, pch=20, col="red")
  legend(x = "topleft", legend = c("train set", "test set"), col = c("black", "red"), pch=c(19, 20))
  
  library(ggplot2)
  ggplot(df_gamma_normal_2, aes(x = duration, y = mean_normal, size = log(size))) + geom_point(shape=1) + 
    ggtitle("mean of normal component against duration") + labs(x = "duration/min", y = "mean/km")
  
  ## linear regression on the means
  plot(df_gamma_normal_2$duration, df_gamma_normal_2$mean_normal^2, xlab="duration/min", ylab="squared mean of normal", xaxp = c(0, 180, 6), pch=20)
  ## model 1: squared mean
  lmod_mean_normal <- lm(mean_normal^2 ~ duration, df_gamma_normal_2, weights = size)
  abline(lmod_mean_normal$coefficients[1], lmod_mean_normal$coefficients[2], col="red")
  # diagnostic
  plot(lmod_mean_normal$residuals ~ duration, df_gamma_normal_2, ylab="residual", xlab="duration/min", pch=20)
  abline(0, 0, col="red", lty=3)
  
  ## model 2: spline
  require(splines)
  lmodb_mean_normal = lm(mean_normal ~ bs(duration, knots = c(18, 37, 60)), data = df_gamma_normal_2, subset = (duration <= 80))
  matplot(1:80, cbind(df_gamma_normal_2$mean_normal[min_dur:80], lmodb_mean_normal$fitted.values), type="pl", xlab="duration/min", ylab="mean of normal", pch=20, lty=1, col=c("black", "red"), xaxp = c(0, 180, 6), main="B-spline")
  abline(v = 18, lty = 5)
  abline(v = 37, lty = 5)
  abline(v = 60, lty = 5)
  # diagnostic
  plot(1:80, lmodb_mean_normal$residuals, xlab="duration/min", ylab="residuals", pch=20, xlim=c(0, 90), xaxp = c(0, 90, 3))
  abline(0, 0, col="red", lty=3)
  abline(v = 18, lty = 5)
  abline(v = 37, lty = 5)
  abline(v = 60, lty = 5)
  qqnorm(lmodb_mean_normal$residuals, ylab="residuals", main="qq plot of spline(for mean of normal)")
  qqline(lmodb_mean_normal$residuals)


  
  ggplot(df_gamma_normal_2, aes(x = duration, y = var_normal, size = log(size))) + geom_point(shape=1) + 
    ggtitle("variance of normal component against duration") + labs(x = "duration/min", y = "var")
  plot(df_gamma_normal_2$duration^2, df_gamma_normal_2$var_normal, xlab="duration^2/min^2", ylab="variance of normal", pch=20)
  
  ## due to stochastic behaviour of variance for durations > 60min, we only model the first 60 data points
  lmod_var_normal <- lm(var_normal ~ I(duration^2), df_gamma_normal_2, subset = (duration <= 60), weights = size)
  plot(df_gamma_normal_2$duration^2, df_gamma_normal_2$var_normal, xlab="duration^2/min^2", ylab="variance of normal", pch=20)
  abline(lmod_var_normal$coefficients[1], lmod_var_normal$coefficients[2], col="red")
  abline(v = 3600, lty = 5)
  # diagnostic
  plot((1:60)^2, lmod_var_normal$residuals, xlab="duration squared", ylab="residuals")
  abline(0, 0, col="red")
  # non-linearity present, use spline model instead
  # this time model up to duration = 80min
  lmod_cs_var_normal <- lm(df_gamma_normal_2$var_normal[min_dur:80] ~ d + I(d^2) + I(d^3) + I((d - 18)^3*(d >= 18)) + I((d - 37)^3*(d >= 37)) + I((d - 60)^3*(d >= 60)))
  matplot(d, cbind(df_gamma_normal_2$var_normal[min_dur:80], lmod_cs_var_normal$fitted.values), type="pl", xlab="duration/min", ylab="variance of normal", pch=20, lty=1, col=c("black", "red"), xaxp = c(0, 180, 6), main="cubic spline")
  abline(v = 18, lty = 5)
  abline(v = 37, lty = 5)
  abline(v = 60, lty = 5)
  plot(df_gamma_normal_2$duration[min_dur:80], lmod_cs_var_normal$residuals, xlab="duration/min", ylab="residuals")
  abline(v = 18, lty = 5)
  abline(v = 37, lty = 5)
  abline(v = 60, lty = 5)
  qqnorm(lmod_cs_var_normal$residuals, ylab="residuals", main="qq plot of cubic spline(for variance of normal)")
  qqline(lmod_cs_var_normal$residuals)
  
  ## model 3: discontinuous spline
  lmodb1 <- lm(var_normal ~ duration + I(duration^2) + I(duration^3), df_gamma_normal_2, subset = (duration <= 18))
  lmodb2 <- lm(var_normal ~ bs(duration, knots = c(60)), df_gamma_normal_2, subset = (duration > 18) & (duration <= 80))
  fitted <- c(lmodb1$fitted.values, lmodb2$fitted.values)
  matplot(1:80, cbind(df_gamma_normal_2$var_normal[1:80], fitted), type="pl", xlab="duration/min", ylab="variance of normal", pch=20, lty=1, col=c("black", "red"), xaxp = c(0, 180, 6), main="discontinuous spline")
  plot(1:18, lmodb1$residuals, xlab="duration/min", ylab="residuals", main="residuals of first part")
  abline(0, 0, col="red")
  plot(19:80, lmodb2$residuals, xlab="duration/min", ylab="residuals", main="residuals of second part")
  abline(0, 0, col="red")
  
  
  plot(df_gamma_normal_2$duration[15:max_dur], par_gamma_normal_2[15:max_dur, 1] , xlab="duration/min", ylab="shape of gamma", xlim = c(0, 180), xaxp = c(0, 180, 6), pch=20)
  # first few points very extreme, probably due to low mixing proportions 
  plot(df_gamma_normal_2$duration[19:max_dur], df_gamma_normal_2$mean_gamma[19:max_dur], xlab="duration/min", ylab="mean of gamma", xlim = c(0, 180), xaxp = c(0, 180, 6), pch=20)
  library(dplyr)
  ggplot(df_gamma_normal_2 %>% slice(19:max_dur), aes(x = duration, y = mean_gamma, size = log(size))) + geom_point(shape=1) + 
    ggtitle("mean of gammma component against duration") + labs(x = "duration/min", y = "mean/km")
  
  ## model: broken stick with 3 segments
  ## turning points: 24， 38
  lmod1 <- lm(mean_gamma ~ duration, df_gamma_normal_2, subset = ((duration > 18) & (duration <= 24)), weights = size)
  lmod2 <- lm(mean_gamma ~ duration, df_gamma_normal_2, subset = ((duration >= 24) & (duration <= 38)), weights = size)
  lmod3 <- lm(mean_gamma ~ duration, df_gamma_normal_2, subset = ((duration >= 38) & (duration <= 60)), weights = size)
  # fitness check
  plot(df_gamma_normal_2$duration[19:max_dur], df_gamma_normal_2$mean_gamma[19:max_dur], xlab="duration/min", ylab="mean of gamma", xlim = c(0, 180), xaxp = c(0, 180, 6), pch=20)
  abline(v = 24, lty = 5)
  abline(v = 38, lty = 5)
  abline(v = 60, lty = 5)
  segments(19, lmod1$coefficients[1] + lmod1$coefficients[2] * 19, 24, lmod1$coefficients[1] + lmod1$coefficients[2] * 24, col="red")
  segments(38, lmod2$coefficients[1] + lmod2$coefficients[2] * 38, 24, lmod2$coefficients[1] + lmod2$coefficients[2] * 24, col="red")
  segments(38, lmod3$coefficients[1] + lmod3$coefficients[2] * 38, 60, lmod3$coefficients[1] + lmod3$coefficients[2] * 60, col="red")
  # diagnostic
  plot(19:24, lmod1$residuals, main="residuals of first segment")
  abline(0, 0, col="red")
  plot(24:38, lmod2$residuals, main="residuals of second segment")
  abline(0, 0, col="red")
  plot(38:60, lmod3$residuals, main="residuals of third segment")
  abline(0, 0, col="red")
  ## due to stochastic behaviour of mean for duration > 60, we do not perform regression on them. 
  
  
  plot(df_gamma_normal_2$duration[19:max_dur], df_gamma_normal_2$var_gamma[19:max_dur], xlab="duration/min", ylab="variance of gamma", xlim = c(0, 180), xaxp = c(0, 180, 6), pch=20)
  ## model: broken stick with 3 segments
  ## turning points: 25， 38
  lmod1 <- lm(var_gamma ~ duration, df_gamma_normal_2, subset = ((duration > 18) & (duration <= 25)), weights = size)
  lmod2 <- lm(var_gamma ~ duration, df_gamma_normal_2, subset = ((duration >= 25) & (duration <= 38)), weights = size)
  lmod3 <- lm(var_gamma ~ duration, df_gamma_normal_2, subset = ((duration >= 38) & (duration <= 60)), weights = size)
  # fitness check
  plot(df_gamma_normal_2$duration[19:max_dur], df_gamma_normal_2$var_gamma[19:max_dur], xlab="duration/min", ylab="variance of gamma", xlim = c(0, 180), xaxp = c(0, 180, 6), pch=20)
  abline(v = 25, lty = 5)
  abline(v = 38, lty = 5)
  abline(v = 60, lty = 5)
  segments(19, lmod1$coefficients[1] + lmod1$coefficients[2] * 19, 25, lmod1$coefficients[1] + lmod1$coefficients[2] * 25, col="red")
  segments(38, lmod2$coefficients[1] + lmod2$coefficients[2] * 38, 25, lmod2$coefficients[1] + lmod2$coefficients[2] * 25, col="red")
  segments(38, lmod3$coefficients[1] + lmod3$coefficients[2] * 38, 60, lmod3$coefficients[1] + lmod3$coefficients[2] * 60, col="red")
  # diagnostic
  plot(19:25, lmod1$residuals, main="residuals of first segment")
  abline(0, 0, col="red")
  plot(25:38, lmod2$residuals, main="residuals of second segment")
  abline(0, 0, col="red")
  plot(38:60, lmod3$residuals, main="residuals of third segment")
  abline(0, 0, col="red")
  ## again we ignore the stochastic part 
  
  
  plot(df_gamma_normal_2$duration[19:max_dur], df_gamma_normal_2$proportion[19:max_dur], xlab="duration/min", ylab="proportion of gamma", xlim = c(0, 180), xaxp = c(0, 180, 6), pch=20)
  ## model: borken stick with 3 segments
  ## turning points: 23, 38, 80
  lmod1 <- lm(proportion ~ duration, df_gamma_normal_2, subset = ((duration > 18) & (duration <= 23)), weights = size)
  lmod2 <- lm(proportion ~ duration, df_gamma_normal_2, subset = ((duration >= 23) & (duration <= 38)), weights = size)
  lmod3 <- lm(proportion ~ duration, df_gamma_normal_2, subset = ((duration >= 38) & (duration <= 80)), weights = size)
  plot(df_gamma_normal_2$duration[19:80], df_gamma_normal_2$proportion[19:80], xlab="duration/min", ylab="proportion of gamma", xlim = c(0, 90), xaxp = c(0, 90, 3), pch=20)
  abline(v = 23, lty = 5)
  abline(v = 38, lty = 5)
  segments(19, lmod1$coefficients[1] + lmod1$coefficients[2] * 19, 23, lmod1$coefficients[1] + lmod1$coefficients[2] * 23, col="red")
  segments(38, lmod2$coefficients[1] + lmod2$coefficients[2] * 38, 23, lmod2$coefficients[1] + lmod2$coefficients[2] * 23, col="red")
  segments(38, lmod3$coefficients[1] + lmod3$coefficients[2] * 38, 80, lmod3$coefficients[1] + lmod3$coefficients[2] * 80, col="red")
  # diagnostic
  plot(19:23, lmod1$residuals, main="residuals of the first segment")
  abline(0, 0, col="red")
  plot(23:38, lmod2$residuals, main="residuals of the second segment")
  abline(0, 0, col="red")
  plot(38:80, lmod3$residuals, main="residuals of the third segment")
  abline(0, 0, col="red")
  qqnorm(lmod1$residuals, ylab="residuals", main="first segment")
  qqline(lmod1$residuals)
  qqnorm(lmod2$residuals, ylab="residuals", main="second segment")
  qqline(lmod2$residuals)
  qqnorm(lmod3$residuals, ylab="residuals", main="third segment")
  qqline(lmod3$residuals)
}

#####################################
## To study the relationship of the coefficients and durations more closely, perform mixture of regression
#####################################
## according to the paper, CEM converges faster. Though SEM, EM are more accurate, this can be compensated by the large sample size here. 
{
  require(MASS)
  ## mean parametrisation of gamma distribution as a member of exponential family
  # phi means dispersion
  dgamma_exp_family <- function(x, mu, phi) {
    theta <- 0 - 1 / mu
    pdf_val <- exp((theta * x - log(mu)) / phi + (1 - phi) * log(x) / phi - log(gamma(1 / phi)))
    if (any(is.nan(pdf_val))) {
      stop("dgamma pdf is NaN")
    }
    return(pdf_val)
  }
  
  ## auxiliary fuction 1: calculate the means of two components
  ## MR means mixture of regressions
  ## coef: is a list of length 2 containing the parameters. 
  ## Each entry with length d_j + p_j + 1
  ## mean function modeled by each component be in the order: gamma, normal(do not give names)
  ## knots: same structure as coef, each entry with length p_j 
    ## knots will not be updated during CEM, no repetition allowed!
  ## d is a VECTOR of length 2 with degrees of truncated power base(for each parameter)
  ## p is a VECTOR of length 2 with number of knots of truncated power base(for each parameter) 
    ## if p_j = 0, polynomial base used. 
  ## return n x 2 matrix. (first row: gamma, second row: normal)
  par_MR <- function(x, y, coef, knots, d, p) {
    ## build the means(first column: gamma, second column: normal) using coefficients
    n <- length(x)
    mu <- matrix(0, n, 2)
    for (j in 1:2) {
      for (i in 0:d[j]) {
        mu[, j] <- mu[, j] + coef[[j]][i+1] * (x^i)
      }
      if (p[j] > 0) {
        for (i in 1:p[j]) {
          mu[, j] <- mu[, j] + coef[[j]][d[j] + i + 1] * (ifelse(x > knots[[j]][i], x - knots[[j]][i], 0))^d[j]
        }
      }
    }
    
    for (j in 1:2) {
      for (i in 1:n) {
        if ((mu[i, j] < 0) || is.nan(mu[i, j])) {
          stop("some mu is negative or is NaN")
        }
      }
    }
    return(mu)
  }
  
  ## auxiliary function 2: calculate required log likelihood
  # phi is dispersion parameter (for gamma as a member of exponential family)
  # var is variance of normal component
  loglik_MR <- function(x, y, coef, knots, pi, d, p, var, phi) {
    n <- length(x)
      
    for (j in 1:2) {
      if (length(coef[[j]]) != (d[j] + 1 + p[j])) {
        stop("make sure d + p + 1 = length of coef")
      }
      
      if (p[j] != length(unlist(knots[j]))) {
        stop("length of knots should be p")
      }
    }
    
    if (length(x) != length(y)) {
      stop("length of x, y should be the same")
    }
    
    mu <- par_MR(x, y, coef, knots, d, p)
    
    
    
    loglik <- 0
    for (i in 1:n) {
      loglik <- loglik + log(pi * dgamma_exp_family(y[i], mu=mu[i, 1], phi=phi) + (1-pi) * dnorm(y[i], mean = mu[i, 2], sd = sqrt(var)))
    }
    
    return(loglik)
  }
  
  
  ## Model of first component: gamma regression (assume initial alpha = 1)
  ## Model of second component: Ordinary Linear Regression
  ## assumption: constant variance, constant mixing proportion
  ## parameters: modeled by truncated power base, with x being duration of journey
  ## knots, d, p same as above. 
  ## pi: initial mixing proportion
  ## init has the same structure as coef in the above function, but with initial parameters
  MR_gamma_normal <- function(x, y, init, knots, pi = 0.5, d, p, tol = 1e-5, max_it = 250) {
    ### load necessary package ##
    
    #### basic checks ####
    if (length(x) != length(y)) {
      stop("length of x, y should be the same")
    }
    
    if ((pi > 1) | (pi < 0)) {
      stop("initial weight pi should be a probability! (between 0, 1)")
    }
    
    coef <- init
    
    for (j in 1:2) {
      if (length(coef[[j]]) != (d[j] + 1 + p[j])) {
        stop("make sure d + p + 1 = length of coef")
      }
      
      if (p[j] != length(knots[[j]])) {
        stop("length of knots should be p")
      }
    }
    
    #### basic set up ####
    n <- length(x)
    ## predictors
    X_1 <- matrix(0, nrow = n, ncol = (d[1] + p[1] + 1))
    for (i in 0:d[1]) {
      X_1[, i + 1] <- x^i
    }
    if (p[1] > 0) {
      for (i in 1:p[1]) {
        X_1[, i + d[1] + 1] <- (ifelse(x > knots[[1]][i], x - knots[[1]][i], 0))^d[1]
      }
    }
    
    X_2 <- matrix(0, nrow = n, ncol = (d[2] + p[2] + 1))
    for (i in 0:d[2]) {
      X_2[, i + 1] <- x^i
    }
    if (p[2] > 0) {
      for (i in 1:p[2]) {
        X_2[, i + d[2] + 1] <- (ifelse(x > knots[[2]][i], x - knots[[2]][i], 0))^d[2]
      }
    }
    
    var_normal <- 0
    for (i in 1:n) {
      var_normal <- var_normal + ( y[i] - sum(X_2[i, ]*coef[[2]]) ) ^ 2
    }
    var_normal <- var_normal / n
    
    ## diff measures difference between log likelihoods of two iterations
    diff <- 1 
    
    ## r is number of iterations 
    r <- 0
    
    ## store log-likelihoods
    ## assume alpha = 1
    phi <- 1
    ll <- c(loglik_MR(x, y, coef, knots, pi, d, p, var_normal, phi))
    
    ## formulas for LM and GLM 
    formula_GLM <- "y ~ X_1[,2]"
    for (i in 2:(d[1] + p[1])) {
      formula_GLM <- paste(formula_GLM, " + X_1[,", toString(i+1), "]", sep="")
    }
    formula_GLM <- formula(formula_GLM)
    
    formula_LM <- "y ~ X_2[,2]"
    for (i in 2:(d[2] + p[2])) {
      formula_LM <- paste(formula_LM, " + X_2[,", toString(i+1), "]", sep="")
    }
    formula_LM <- formula(formula_LM)
    
    #### iterations ####
    while (diff > tol) {
      if (r >= max_it) {
        print(paste("maximum iterations:", toString(max_it), "is reached", sep=" "))
        return(list("coef" = coef, "degrees" = d, "log_lik" = ll, "max_it_reached" = TRUE))
      }
      
      r = r + 1
      

      #### E-step ####
      
      mu <- par_MR(x, y, coef, knots, d, p)
      
      # w_1 is the vector of probabilities of each data point y belonging to gamma component conditioned on x
      w_1 <- c()
      w_1 <- (pi * dgamma_exp_family(y, mu=mu[, 1], phi=phi)) / 
        (pi * dgamma_exp_family(y, mu=mu[, 1], phi=phi) + (1-pi) * dnorm(y, mean=mu[, 2], sd = sqrt(var_normal)))
      # probabilities of belonging to normal component
      w_2 <- 1 - w_1
      
      if ((length(w_1[w_1 == 0]) == length(w_1)) || (length(w_1[w_1 == 1]) == length(w_1))) {
        stop("fails, one component eliminated")
      }
      
      #### M-step #### 
      # update coefficients of the first component (X_1[,1] is intercept, not included in the formula)
      ###################### may have to change the weights ######################
      gamma_GLM <- glm(formula_GLM, weights = w_1, family = Gamma(link = "inverse"), control = glm.control(maxit = 150))
      coef[[1]] <- unname(gamma_GLM$coefficients)
      est_shape <- gamma.shape(gamma_GLM)
      phi <- 1 / est_shape$alpha
      
      # WLS for second component
      lm_result <- lm(formula_LM, weights = w_2)
      coef[[2]] <- lm_result$coefficients
      
      # update mixing proportion
      pi <- sum(w_1) / n
      
      # update variance (of normal component)
      var_normal <- 0 
      for (i in 1:n) {
        var_normal <- var_normal + w_2[i] * ( y[i] - sum(X_2[i, ]*coef[[2]]) ) ^ 2
      }
      var_normal <- var_normal / sum(w_2)
      
      ll <- c(ll, loglik_MR(x, y, coef, knots, pi, d, p, var_normal, phi))
      
      ## update difference
      diff <- abs(diff(tail(ll, 2)))
    }
    
    return(list("coef" = coef, "degrees" = d, "knots" = p, "log_lik" = ll, "max_it_reached" = FALSE))
  }
  
  
  ## we only consider durations >= 19min as values before that does not have mixture
  ## auxiliary function 3: choose initial coefficients
  initial_coef <- function(knots, d, p) {
    init <- list("gamma" = c(), "normal" = c())
    mean_gamma <- df_gamma_normal_2$mean_gamma[19:max_dur]
    mean_normal <- df_gamma_normal_2$mean_normal[19:max_dur]
    sizes <- df_gamma_normal_2$size[19:max_dur]
    # duration
    dur <- 19:max_dur
    formula_1 <- "mean_gamma ~ dur"
    for (i in 2:d[1]) {
      formula_1 <- paste(formula_1, " + I(dur ^ (", toString(i), "))", sep="")
    }
    if (p[1] > 0 ) {
    for (i in 1:p[1]) {
        formula_1 <- paste(formula_1, 
                             " + I(ifelse(dur > knots[[1]][", toString(i), "], dur - knots[[1]][", toString(i), "], 0) ^ (d[1]))",
                             sep="")
      }
    }
    formula_1 <- formula(formula_1)
    
    lmod_1 <- lm(formula_1, weights = sizes)
    init[[1]] <- lmod_1$coefficients
    
    plot(dur, mean_gamma, xlim = c(0, 180), xaxp = c(0, 180, 6))
    lines(dur, lmod_1$fitted.values, col="red")
    
    formula_2 <- "mean_normal ~ dur"
    for (i in 2:d[2]) {
      formula_2 <- paste(formula_2, " + I(dur ^ (", toString(i), "))", sep="")
    }
    if (p[2] > 0) {
      for (i in 1:p[2]) {
        formula_2 <- paste(formula_2, 
                           " + I(ifelse(dur > knots[[2]][", toString(i), "], dur - knots[[2]][", toString(i), "], 0) ^ (d[2]))",
                           sep="")
      }
    }
    formula_2 <- formula(formula_2)
    lmod_2 <- lm(formula_2, weights = sizes)
    init[[2]] <- lmod_2$coefficients
    plot(dur, mean_normal, xlim = c(0, 180), xaxp = c(0, 180, 6))
    lines(dur, lmod_2$fitted.values, col="red")
    
    return(init)
  }
  
  
  ## test the algorithm on a small set
  test_proportion <- 0.005
  total_len <- length(df_one_side_trimmed$duration)
  upper_lim <- ceiling(total_len * test_proportion)
  test_x <- df_one_side_trimmed$duration[1:upper_lim]
  test_y <- df_one_side_trimmed$dist[1:upper_lim]
  test_df <- data.frame("x" = test_x, "y" = test_y)
  test_df <- test_df[(test_df$x) > 18 & (test_df$y > 0), ]
  
  ## trial 1 (knots decided from the graph)
  knots_1 <- list("gamma" = c(24, 38, 60), "normal" = c(38, 60))
  init_1 <- initial_coef(knots_1, d = c(3, 3), p = c(3, 2))
  test_result_1 <- MR_gamma_normal(test_df$x, test_df$y, init_1, knots_1, pi = 0.5, d = c(3, 3), p = c(3, 2), tol = 1e-5, max_it = 250)
  
  ## trial 2 (knots are exactly 30min apart)
  knots_2 <- list("gamma" = seq(from = 30, to = 150, by = 30), "normal" = seq(from = 30, to = 150, by = 30))
  init_2 <- initial_coef(knots_2, d = c(3, 3), p = c(5, 5))
  test_result_2 <- MR_gamma_normal(test_df$x, test_df$y, init_2, knots_2, pi = 0.5, d = c(3, 3), p = c(5, 5), tol = 1e-5, max_it = 250)
  
  ## trial 3 (no knot)
  knots_3 <- list("gamma" = c(), "normal" = c())
  init_3 <- initial_coef(knots_3, d = c(3, 3), p = c(0, 0))
  test_result_3 <- MR_gamma_normal(test_df$x, test_df$y, init_3, knots_3, pi = 0.5, d = c(3, 3), p = c(0, 0), tol = 1e-5, max_it = 250)
  
  ## all fails. Try the data of the second week
  ## for accuracy, we rather focus on the data between 19min and 60min
  week_2_start <- df$start_time[1] + 3600 * 24 * 7
  week_2_end <- week_2_start + 3600 * 24 * 7
  df_week2_MR <- df_one_side_trimmed[(df_one_side_trimmed$duration > 18) & (df_one_side_trimmed$duration <= 60) &
                                       (df_one_side_trimmed$dist > 0) & 
                                       (df_one_side_trimmed$start_time >= week_2_start) & (df_one_side_trimmed$start_time <= week_2_end), ]
  n <- length(df_week2_MR$dist)
  
  
  train_proportion_MR <- 0.7
  train_length_MR <- ceiling(n * train_proportion_MR)
  dur_train_MR <- df_week2_MR$duration[1:train_length_MR]
  dist_train_MR <- df_week2_MR$dist[1:train_length_MR]
  dur_test_MR <- df_week2_MR$duration[(train_length_MR+1):n]
  dist_test_MR <- df_week2_MR$dist[(train_length_MR+1):n]
  
  ## knots decided from the graph
  knots_1 <- list("gamma" = c(24, 38), "normal" = c(38))
  init_1 <- initial_coef(knots_1, d = c(3, 3), p = c(2, 1))
  MR_result_1 <- MR_gamma_normal(dur_train_MR, dist_train_MR, init_1, knots_1, pi = 0.5, d = c(3, 3), p = c(2, 1), tol = 1e-5, max_it = 250)
  
  ## knots are exactly 30min apart
  knots_2 <- list("gamma" = c(30), "normal" = c(30))
  init_2 <- initial_coef(knots_1, d = c(3, 3), p = c(1, 1))
  MR_result_2 <- MR_gamma_normal(dur_train_MR, dist_train_MR, init_2, knots_2, pi = 0.5, d = c(3, 3), p = c(1, 1), tol = 1e-5, max_it = 250)
  
  ## no knot
  MR_result_3 <- MR_gamma_normal(dur_train_MR, dist_train_MR, init_3, knots_3, pi = 0.5, d = c(3, 3), p = c(0, 0), tol = 1e-5, max_it = 250)
  
  ## try whole data set 
  df_whole_MR <- df_one_side_trimmed[(df_one_side_trimmed$duration > 18) & (df_one_side_trimmed$duration <= 60) &
                                       (df_one_side_trimmed$dist > 0), ]
  n <- length(df_whole_MR$dist)
  
  train_proportion_MR <- 0.7
  train_length_MR <- ceiling(n * train_proportion_MR)
  dur_train_MR <- df_whole_MR$duration[1:train_length_MR]
  dist_train_MR <- df_whole_MR$dist[1:train_length_MR]
  dur_test_MR <- df_whole_MR$duration[(train_length_MR+1):n]
  dist_test_MR <- df_whole_MR$dist[(train_length_MR+1):n]
  
  ## knots decided from the graph
  MR_whole_result_1 <- MR_gamma_normal(dur_train_MR, dist_train_MR, init_1, knots_1, pi = 0.5, d = c(3, 3), p = c(2, 1), tol = 1e-5, max_it = 250)
  ## knots are exactly 30min apart
  MR__whole_result_2 <- MR_gamma_normal(dur_train_MR, dist_train_MR, init_2, knots_2, pi = 0.5, d = c(3, 3), p = c(1, 1), tol = 1e-5, max_it = 250)
  ## no knot
  MR__whole_result_3 <- MR_gamma_normal(dur_train_MR, dist_train_MR, init_3, knots_3, pi = 0.5, d = c(3, 3), p = c(0, 0), tol = 1e-5, max_it = 250)
  
  ## all fails, we have to make pi depend on d(duration). 
}


###### EM, but parameters allowed to depend on d ###### 
{
  ## the main idea is optimisation on the log likelihood
  ## coef is a list of five rows: logit(pi), ln(mean), ln(var), ln(alpha), ln(beta)
  ## knots: same as above
  ## d has length 5
  ## p has length 5
  par_opt <- function(x, coef, knots, d, p) {
    if ((length(which(is.nan(x))) > 0) || (length(which(is.na(x))) > 0)) {
      stop("NaN or NA in x")
    }
    
    n <- length(x)
    par <- matrix(0, n, 5)
    for (j in 1:5) {
      for (i in 0:d[j]) {
        par[, j] <- par[, j] + coef[[j]][i+1] * (x^i)
      }
      if (p[j] > 0) {
        for (i in 1:p[j]) {
          par[, j] <- par[, j] + coef[[j]][d[j] + i + 1] * (ifelse(x > knots[[j]][i], x - knots[[j]][i], 0))^d[j]
        }
      }
    }
    
    for (j in 1:5) {
      for (i in 1:n) {
        if (is.nan(par[i, j]) || is.na(par[i, j])) {
          stop(paste("some parameter(s) is NaN or NA: ", toString(i), ", ", toString(j), sep=""))
        }
        if ((j == 5) & ((exp(par[i, j]) == 0) || (is.infinite(exp(par[i, j]))))) {
          print(paste("ln(beta): ", toString(par[i, j])))
          stop(paste("some beta(s) is 0 or Inf: ", toString(i), ", ", toString(j), sep=""))
        }
      }
    }
    
    
    return(par)
  }
  
  
  ## when we only need to find pi. 
  ## coef, knots are vectors
  ## d, p are scalars
  pi_opt <- function(x, coef, knots, d, p) {
    n <- length(x)
    par <- rep(0, times = n)
    for (i in 0:d) {
      par <- par + coef[i+1] * (x^i)
    }
    if (p > 0) {
      for (i in 1:p) {
        par <- par + coef[d + i + 1] * (ifelse(x > knots[i], x - knots[i], 0))^d
      }
    }
  
    return(inv.logit(par))
  }
  
  ## the main code
  ## again, init consists of the initial coefficients. 
  EM_opt <- function(x, y, init, knots, d, p, tol=1e-5, max_it = 250) {
    library("boot")
    
    ###### basic checks #######
    if (length(x) != length(y)) {
      stop("length of x, y should be the same")
    }
    
    for (j in 1:5) {
      if (length(knots[[j]]) != p[j]) {
        stop("please check that length of knots is p")
      }
      if (length(init[[j]]) != d[j] + p[j] + 1) {
        stop("please make sure length of init is d + p + 1")
      }
    }
    
    
    ######### basic set up ########
    coef <- init
    pi <- pi_opt(x, coef[[1]], knots[[1]], d[1], p[1])
    
    ## for converting coef between list and vector
    bound1 <- d[1] + p[1] + 1
    bound2 <- bound1 + d[2] + p[2] + 1
    bound3 <- bound2 + d[3] + p[3] + 1
    bound4 <- bound3 + d[4] + p[4] + 1
    bound5 <- bound4 + d[5] + p[5] + 1
    
    ###### inner auxiliary function: ###### 
    ## used for maximisation, and convergence criteria
    ## beta does NOT include coefficients for logit(pi)
    loglik_opt <- function(beta) {
      if (length(roll_x) != length(roll_y)) {
        stop("length of x, y in loglik_opt differs")
      }
      
      if (length(beta) != (bound5 - bound1)) {
        stop("length of input is wrong")
      }
      beta <- c(coef[[1]], beta)
      
      if ((length(which(is.nan(beta))) > 0) || (length(which(is.na(beta))) > 0)) {
        stop("NaN or NA in coefficients")
      }

      beta2 <- beta[(bound1+1):bound2]
      beta3 <- beta[(bound2+1):bound3]
      beta4 <- beta[(bound3+1):bound4]
      beta5 <- beta[(bound4+1):bound5]
      
      coef_rebuild <- list(coef[[1]], beta2, beta3, beta4, beta5)
      
      par <- par_opt(roll_x, coef_rebuild, knots, d, p)
      
      loglik <- sum( log( inv.logit(par[,1]) * dgamma(roll_y, shape = exp(par[, 4]), scale = exp(par[, 5])) + (1-inv.logit(par[,1])) * dnorm(roll_y, mean =exp(par[, 2]), sd = sqrt(exp(par[, 3]))) ))
      
      return(loglik)
    }
    
    
    #### basic set up ####
    n <- length(x)
    
    ## diff measures difference between log likelihoods of two iterations
    diff <- 1 
    
    ## r is number of iterations 
    r <- 0
    
    ## store log-likelihoods
    initial_coef <- c(coef[[2]], coef[[3]], coef[[4]], coef[[5]])
    ## roll_x, y are for loglik_opt calculations
    roll_x <- x
    roll_y <- y
    ll <- c(loglik_opt(initial_coef))
    
    #### iterations ####
    while (diff > tol) {
      if (r >= max_it) {
        print(paste("maximum iterations:", toString(max_it), "is reached", sep=" "))
        return(list("coef" = coef, "degrees" = d, "knots" = knots, "log_lik" = ll, "max_it_reached" = TRUE))
      }
      
      r = r + 1
      
      print(paste("iteration: ", toString(r), sep=""))
      
      #### E-step ####
      par <- par_opt(x, coef, knots, d, p)
      pi <- inv.logit(par[,1])
      # w_1 is the vector of probabilities of each data point y belonging to gamma component conditioned on x
      w_1 <- c()
      w_1 <- (pi * dgamma(y, shape=exp(par[, 4]), scale=exp(par[,5]))) / 
        (pi * dgamma(y, shape=exp(par[, 4]), scale=exp(par[,5])) + (1-pi) * dnorm(y, mean=exp(par[,2]), sd = sqrt(exp(par[,3]))))
      
      nan_index <- which(is.nan(w_1))
      if (length(nan_index) != 0) {
        print("------ parameters producing NaN w_1 ---------")
        print(exp(par[nan_index, 2:5]))
        print("------ corresponding pi producing NaN w_1 ---------")
        print(pi[nan_index])
        stop("NaNs produced in w_1")
      }
      
      if ((length(w_1[w_1 == 0]) == length(w_1)) || (length(w_1[w_1 == 1]) == length(w_1))) {
        stop("fails, one component eliminated")
      }
      
      ## prop_0 <- length(w_1[w_1 == 0]) /  length(w_1)
      ## prop_1 <- length(w_1[w_1 == 1]) / length(w_1)
      ## print(paste("proportion of responsibilities being 1:", toString(prop_0)))
      ## print(paste("proportion of responsibilities being 0:", toString(prop_1)))
      
      
      ## 0 or 1s in w_1 causes problem to logit, so replace them with small deviation
      w_1 <- replace(w_1, w_1 == 0, 1e-4)
      w_1 <- replace(w_1, w_1 == 1, 1 - 1e-4)
      
      # probabilities of belonging to normal component
      w_2 <- 1 - w_1
      
      #### M-step #### 
      # update coefficients for pi 
      logit_pi <- logit(w_1)
      formula <- paste("logit_pi", " ~ x", sep="")
      for (i in 2:d[1]) {
        formula <- paste(formula, " + I(x ^ (", toString(i), "))", sep="")
      }
      if (p[1] > 0 ) {
        for (i in 1:p[1]) {
          formula <- paste(formula, 
                           " + I(ifelse(x > knots[[1]][", toString(i), "], x - knots[[1]][", toString(i), "], 0) ^ (d[1]))",
                           sep="")
        }
      }
      formula <- formula(formula)
      coef[[1]] <- lm(formula)$coefficients
      
      if (length(which(is.nan(coef[[1]]))) > 0) {
        stop("NaN produced when pi is updated")
      }
      
      if (length(which(is.na(coef[[1]]))) > 0) {
        stop("NA produced when pi is updated")
      }
      
      pi <- pi_opt(x, coef[[1]], knots[[1]], d[1], p[1])
      extreme_pi <- length(pi[pi < 1e-3 || (1 - pi) < 1e-3]) / length(pi)
      semi_extreme_pi <- length(pi[pi < 1e-2 || (1 - pi) < 1e-2]) / length(pi)
      print(paste("proportion of extreme pi: ", toString(extreme_pi)))
      print(paste("proportion of relatively extreme pi: ", toString(semi_extreme_pi)))
      
      print("pi updated")
                
      ## update coefficients for other parameters
      initial_coef <- c(coef[[2]], coef[[3]], coef[[4]], coef[[5]])
      # perform stochastic optimisation, take a small sample
      data_df <- data.frame("x" = x, "y" = y)
      n_sample <- 1000
      sample_df <- data_df[sample(nrow(data_df), n_sample), ]
      roll_x <- sample_df$x
      roll_y <- sample_df$y
      opt_result <- optim(initial_coef, loglik_opt, method = "Nelder-Mead", control = list(fnscale = -1))
      print("optimisation finished")
      placeholder <- rep(0, times = bound1)
      opt_coef <- c(placeholder, opt_result$par)
      coef[[2]] <- opt_coef[(bound1+1):bound2]
      coef[[3]] <- opt_coef[(bound2+1):bound3]
      coef[[4]] <- opt_coef[(bound3+1):bound4]
      coef[[5]] <- opt_coef[(bound4+1):bound5]
      
      roll_x <- x
      roll_y <- y
      coef_vec <- c(coef[[2]], coef[[3]], coef[[4]], coef[[5]])
      new_ll <- loglik_opt(coef_vec)
      if (is.nan(new_ll)) {
        stop("new loglik is NaN")
      }
      ll <- c(ll, new_ll)
      
      ## update difference
      diff <- abs(diff(tail(ll, 2)))
    }
    
    return(list("coef" = coef, "degrees" = d, "knots" = knots, "log_lik" = ll, "max_it_reached" = FALSE))
  }
  
  
  
  ## initial value chosen to be the fit from individual models
  initial_coef_opt <- function(knots, d, p) {
    init <- vector(mode = "list", length = 5)
    pi <- df_gamma_normal_2$proportion[19:max_dur]
    logit_pi <- logit(pi)
    alpha <- (df_gamma_normal_2$mean_gamma[19:max_dur])^2 / df_gamma_normal_2$var_gamma[19:max_dur]
    log_alpha <- log(alpha)
    beta <- df_gamma_normal_2$var_gamma[19:max_dur] / df_gamma_normal_2$mean_gamma[19:max_dur]
    log_beta <- log(beta)
    mu <- df_gamma_normal_2$mean_normal[19:max_dur]
    log_mu <- log(mu)
    var <- df_gamma_normal_2$var_normal[19:max_dur]
    log_var <- log(var)
    sizes <- df_gamma_normal_2$size[19:max_dur]
    # duration
    dur <- 19:max_dur
    
    par_individual <- list(logit_pi, log_mu, log_var, log_alpha, log_beta)
    
    par_names <- c("logit_pi", "log_mu", "log_var", "log_alpha", "log_beta")
    
    for (j in 1:5) {
      formula <- paste(par_names[j], " ~ dur", sep="")
      for (i in 2:d[j]) {
        formula <- paste(formula, " + I(dur ^ (", toString(i), "))", sep="")
      }
      if (p[j] > 0 ) {
        for (i in 1:p[j]) {
          formula <- paste(formula, 
                             " + I(ifelse(dur > knots[[j]][", toString(i), "], dur - knots[[j]][", toString(i), "], 0) ^ (d[j]))",
                             sep="")
        }
      }
      formula <- formula(formula)
      
      lmod <- lm(formula, weights = sizes)
      init[[j]] <- lmod$coefficients
      
      plot(dur, par_individual[[j]], xlim = c(0, 180), xaxp = c(0, 180, 6), xlab="duration/min", ylab=par_names[j])
      lines(dur, lmod$fitted.values, col="red")
    }
    
    return(init)
  }
  
  
  #### test with week 2 data again (ONLY duration between 19-60 taken)
  train_proportion <- 0.7
  n <- length(df_week2_MR$duration)
  train_length_week_2 <- ceiling(n * train_proportion)
  dur_train_week2 <- df_week2_MR$duration[1:train_length_week_2]
  dist_train_week2 <- df_week2_MR$dist[1:train_length_week_2]
  dur_test_MR_week2 <- df_week2_MR$duration[(train_length_week_2+1):n]
  dist_test_MR_week2 <- df_week2_MR$dist[(train_length_week_2+1):n]
  
  ## model 1: no knot
  knots_1 <- list(c(), c(), c(), c(), c())
  init_1 <- initial_coef_opt(knots, d = c(3, 3, 3, 3, 3), p = c(0, 0, 0, 0, 0))
  EM_opt_result_1 <- EM_opt(dur_train_week2, dist_train_week2, init_1, knots_1, d = c(3, 3, 3, 3, 3), p = c(0, 0, 0, 0, 0), tol = 1e-3)
  # doesn't converge
  plot(1:251, EM_opt_result_1$log_lik, xlab="iteration", ylab="log likelihood", pch=20)
  
  ## model 2: equal distant knots
  knots_2 <- list(c(30), c(30), c(30), c(30), c(30))
  init_2 <- initial_coef_opt(knots_2, d = c(3, 3, 3, 3, 3), p = c(1, 1, 1, 1, 1))
  EM_opt_result_2 <- EM_opt(dur_train_week2, dist_train_week2, init_2, knots_2, d = c(3, 3, 3, 3, 3),  p = c(1, 1, 1, 1, 1), tol = 1e-3)
  
  ## model 3: knots from graphs
  knots_3 <- list(c(23, 38), c(37), c(37), c(24, 37), c(37))
  init_3 <- initial_coef_opt(knots_3, d = c(3, 3, 3, 3, 3), p = c(2, 1, 1, 2, 1))
  EM_opt_result_3 <- EM_opt(dur_train_week2, dist_train_week2, init_3, knots_3, d = c(3, 3, 3, 3, 3), p = c(2, 1, 1, 2, 1), tol = 1e-3)
  

}


#### optimisation with beta(scale of gamma component) fixed ###### 
{
  ## all of above fails, Inf mean, var, alpha, beta values occurs, so make them all fixed in this section 
  
  ## pi_opt is the same as before
  
  ## the main code
  ## d, p are scalars this time
  ## knots, init are a vectors. Init represents initial value of coefficients for pi
  ## again, init consists of the initial coefficients, but only for mixing proportion
  ## par is a matrix of the parameters taken from individual models
    # columns of par: mean, var, alpha, beta
    # par remains fixed
  EM_opt_2 <- function(x, y, init, par, knots, d, p, tol=1e-5, max_it = 250) {
    library("boot")
    
    ###### basic checks #######
    if (length(x) != length(y)) {
      stop("length of x, y should be the same")
    }
    
    if (length(x) != length(par[,1])) {
      stop("length of par, x should be the same")
    }
    
    if (length(knots) != p) {
        stop("please check that length of knots is p")
    }
    if (length(init) != d + p + 1) {
      stop("please make sure length of init is d + p + 1")
    }
    
    
    ######### basic set up ########
    coef <- init
    pi <- pi_opt(x, coef, knots, d, p)
    n <- length(x)
    diff <- 1 
    ## r is number of iterations 
    r <- 0
    
    ###### inner auxiliary function: ###### 
    ## used for maximisation, and convergence criteria
    ## beta does NOT include coefficients for logit(pi)
    loglik_opt <- function(pi) {
      gamma_comp_pdf <- dgamma(y, shape = par[,3], scale = par[,4])
      if (length(gamma_comp_pdf) != length(x)) {
        stop("length of gamma_comp_pdf wrong")
      }
      loglik <- sum( log( pi * gamma_comp_pdf + (1-pi) * dnorm(y, mean =par[, 1], sd = sqrt(par[, 2])) ) )
      
      return(loglik)
    }
    
    ll <- c(loglik_opt(pi))
    
    #### iterations ####
    while (diff > tol) {
      if (r >= max_it) {
        print(paste("maximum iterations:", toString(max_it), "is reached", sep=" "))
        return(list("coef" = coef, "degrees" = d, "knots" = knots, "log_lik" = ll, "max_it_reached" = TRUE))
      }
      
      r = r + 1
      
      #### E-step ####
      # w_1 is the vector of probabilities of each data point y belonging to gamma component conditioned on x
      w_1 <- c()
      w_1 <- (pi * dgamma(y, shape=par[,3], scale=par[,4])) / 
        (pi * dgamma(y, shape=par[,3], scale=par[,4]) + (1-pi) * dnorm(y, mean=par[,1], sd = sqrt(par[,2])))
      
      nan_index <- which(is.nan(w_1))
      if (length(nan_index) != 0) {
        stop("NaNs produced in w_1")
      }
      
      if ((length(w_1[w_1 == 0]) == length(w_1)) || (length(w_1[w_1 == 1]) == length(w_1))) {
        stop("fails, one component eliminated")
      }
      
      
      ## 0 or 1s in w_1 causes problem to logit, so replace them with small deviation
      w_1 <- replace(w_1, w_1 == 0, 1e-4)
      w_1 <- replace(w_1, w_1 == 1, 1 - 1e-4)
      
      # probabilities of belonging to normal component
      w_2 <- 1 - w_1
      
      #### M-step #### 
      # update coefficients for logit(pi)
      logit_pi <- logit(w_1)
      formula <- paste("logit_pi", " ~ x", sep="")
      for (i in 2:d) {
        formula <- paste(formula, " + I(x ^ (", toString(i), "))", sep="")
      }
      if (p > 0 ) {
        for (i in 1:p) {
          formula <- paste(formula, 
                           " + I(ifelse(x > knots[", toString(i), "], x - knots[", toString(i), "], 0) ^ d)",
                           sep="")
        }
      }
      formula <- formula(formula)
      coef <- lm(formula)$coefficients
      
      if (length(which(is.nan(coef))) > 0) {
        stop("NaN produced when pi is updated")
      }
      
      if (length(which(is.na(coef))) > 0) {
        stop("NA produced when pi is updated")
      }
      
      pi <- pi_opt(x, coef, knots, d, p)
      extreme_pi <- length(pi[pi < 1e-3 || (1 - pi) < 1e-3]) / length(pi)
      semi_extreme_pi <- length(pi[pi < 1e-2 || (1 - pi) < 1e-2]) / length(pi)

      
      new_ll <- loglik_opt(pi)
      if (is.nan(new_ll)) {
        stop("new loglik is NaN")
      }
      ll <- c(ll, new_ll)
      
      ## update difference
      diff <- abs(diff(tail(ll, 2)))
    }
    
    return(list("coef" = coef, "degrees" = d, "knots" = knots, "log_lik" = ll, "max_it_reached" = FALSE))
  }
  
  
  ## check if there is any NaN, Inf, NA in vec
  irregular_data_check <- function(vec) {
    if (length(which(is.infinite(vec))) > 0) {
      print("Waning: Inf in vec")
    } else if (length(which(is.na(vec))) > 0) {
      print("Waning: NA in vec")
    } else if (length(which(is.nan(vec))) > 0) {
      print("Waning: NaN in vec")
    } else {
      print("No abnormal data detected")
    }
  }
  
  # fixed means these parameters are fixed
  fixed_par_1 <- exp(par_opt(dur_train_week2, init_1, knots_1, d = c(3, 3, 3, 3, 3), p = c(0, 0, 0, 0, 0))[,2:5])
  init_pi_coef_1 <- init_1[[1]]
  EM_opt2_result_1 <- EM_opt_2(dur_train_week2, dist_train_week2, init_pi_coef_1, fixed_par_1, knots_1[[1]], d = 3, p = 0, tol = 1e-3)
  plot((1:251), EM_opt2_result_1$log_lik)
  
  fixed_par_2 <- exp(par_opt(dur_train_week2, init_2, knots_2, d = c(3, 3, 3, 3, 3), p = c(1, 1, 1, 1, 1))[,2:5])
  init_pi_coef_2 <- init_2[[1]]
  EM_opt2_result_2 <- EM_opt_2(dur_train_week2, dist_train_week2, init_pi_coef_2, fixed_par_2, knots_2[[1]], d = 3,  p = 1, tol = 1e-3)
  plot((1:251), EM_opt2_result_2$log_lik)
  
  fixed_par_3 <- exp(par_opt(dur_train_week2, init_3, knots_3, d = c(3, 3, 3, 3, 3), p = c(2, 1, 1, 2, 1))[,2:5])
  init_pi_coef_3 <- init_3[[1]]
  EM_opt2_result_3 <- EM_opt_2(dur_train_week2, dist_train_week2, init_pi_coef_3, fixed_par_3, knots_3[[1]], d = 3, p = 2, tol = 1e-3)
  plot((1:251), EM_opt2_result_3$log_lik)
  
  ## fails, doesn't converge. 
}

####### try coordinate ascent method instead of optim() ########
{
  # f is the k-dimensional function to optimise (but input to f must be a single vector with length k!)
  # tol is used as convergence criteria
  # interval is a list of length k, each entry is a vectors of length 2, specifying the interval to search for minima/maxima
  optimize_coord <- function(f, init, interval, tol = 1e-3, maximum = FALSE) {
    input_test <- try({
      f(init)
    }, silent = TRUE)
    if (class(input_test) == "try-error") {
      print(geterrmessage())
      stop("function cannot be evaluated at init")
    }
    
    if (length(f(init)) > 1) {
      stop("f must have 1-dimensional output")
    }
    
    ## number of coordinates
    k <- length(init)
    
    if (length(interval) != k) {
      stop("length of interval should be the same as input of f")
    }
    
    ## euclidean distance of points between two iterations 
    diff <- 1
    ## current position
    curr <- init
    ## number of iterations 
    r <- 0
    
    repeat {
      r <- r + 1
      prev <- curr
      
      for (i in 1:k) {
        fi <- function(xi) {
          if (length(xi) > 1) {
            stop("xi should be a scalar")
          }
          point_copy <- curr
          point_copy[i] <- xi
          return(f(point_copy))
        }
        
        opt_result <- optimize(fi, interval[[i]], maximum = maximum)
        opt_xi <- opt_result$minimum
        
        curr[i] <- opt_xi
      }
      
      diff <- sqrt(sum((curr - prev)^2))
      
      if (diff < tol) {
        break
      }
    }
    
    fval <- f(curr)
    
    return(list("coordinates" = curr, "iterations" = r, "fval" = fval, "diff" = diff))
  }
}



##### distance conditioned on starting station #####
{
  ## plot some examples
  ## this number includes those abandoned station ids
  n_stations <- tail(stations$Station.Id, 1)
  dist_by_start <-vector(mode = "list", length = n_stations)
  for (i in 1:n_stations) {
    dist_by_start[[i]] <- df_one_side_trimmed[df_one_side_trimmed$start_id == i,]$dist
  }
  
  ## exclude zero distance
  plot_dist_hist_by_start <- function(start_id, breaks = 100, xlim = c(0, 10)) {
    library('modeest')
    if (length(dist_by_start[[start_id]]) == 0) {
      stop("No data")
    }
    station_name <- stations$StationName[stations$Station.Id == start_id]
    
    dists <- dist_by_start[[start_id]]
    
    dists <- dists[dists != 0]
    
    ## rounding and count
    dists <- round(dists, digits = 3)
    dists_count <- as.data.frame(table(dists))
    dists_count <- dists_count[order(dists_count$Freq),]
    print(tail(dists_count, 30))
    
    mode_dist <- tail(dists_count$dists, 1)
    n_occur <- tail(dists_count$Freq, 1)
    print(paste("mode distance: ", toString(mode_dist), "km with occurrence ", toString(n_occur), sep=""))
    
    hist(dists, breaks = breaks, xlab = "distance/km", xlim = xlim, main = paste("distances of journeys starting at \n", station_name, sep=""))
    
  }
  
  plot_dist_hist_by_start(785, breaks=200)
  plot_dist_hist_by_start(191, breaks=200, xlim = c(0, 6))
  plot_dist_hist_by_start(111)
  plot_dist_hist_by_start(248)
  plot_dist_hist_by_start(789)
  plot_dist_hist_by_start(307)
  plot_dist_hist_by_start(303)
}

