#' Four Point SSID Filter
#'
#'identifies probable SS in a PV, for use in either process engineering or economic data
#' @param data Vector containing PV to be analyzed
#' @param tcrit_u The upper lim for the t-statistic, probable TS when returned t exceeds
#' @param tcrit_l The lower lim for the t-statistic, probable SS when returned t is below
#' @param n The window size for analysis
#' @param ewma The exponential weighted average filter value (typ. 0.1)
#' @param stp Step size (typ. 1)
#' @return Dataframe containing T-statistics and SS indicator value
#' @examples 
#' SSID_4ptFilter <- (data = ProcessData$Variable, tcrit_u = 3.2, tcrit_l = 1.0, n = 10, ewma = 0.07, stp = 1);
#' @export


SSID_4ptFilter <- function(data, tcrit_u, tcrit_l, n, ewma, stp) {
  
  # Title: 4 Points Filter Function
  # Auth: Max Horner
  # Description: Statistical filter for SS and TS indication of PV. 
  # Uses 4 corner approach where 4 points in a moving window are assessed
  # by taking the ratio of the difference between the min/max filtered values
  # and the filtered variance for the window. 4 points selected as a balance for
  # computing performance vs ability to identify periodicity in data. 
  # References: Dr. R Rhinehard 2019, OSU
  
  # Initialize required constants, vectors, and variables.
  progress(0/100, max.value = 100, progress.bar=TRUE)
  SSData <- setNames(data.frame(matrix(ncol = 11, nrow = length(data))), c("data","tstat", "ss","iread1","iread2","iread3","iread4","y1filt","y2filt","y3filt","y4filt"))
  
  if(stp < 1) {stp <- 1}
  y <- rep(0,n)                     # Init working window to length of window as determined by user. 
  yold <- 0                         # Init n-1 y variable as 0
  data_var<- 0                      # Init variance as 0
  cewma <- 1- ewma
  g <- ((5^0.5-1)/2)
  random <- 0
  # 
  iput <- 1                         # Write data for placing new data in array, read num for most recent
  iread4 <- n                       # Start num for oldest dataset
  iread3 <- as.integer((1-g)*n)     # Start num for third dataset
  iread2 <- as.integer((g*n))       # Start num for second data set
  
  y1filt <- 0                       # Initialize filtered variables to 0
  y2filt <- 0
  y3filt <- 0 
  y4filt <- 0
  
  SS <- 0.5                         # Initialize SSI as indeterminate
  noise_ampl <- 0                   # Add in later as needed, noise addition is currently unused
  
  for(i in 2:length(data)) {                                        # Loop through dataset, filter, data, find SS
    progress((i/length(data))*100,
             max.value = (length(data)/length(data))*100,
             progress.bar=TRUE)                                     # Progress bar so user doesn't rage quit because it takes ~10-15 minutes to process large datasets
    
    y[iput]       <- data[i] + noise_ampl*random*0                  # Pull PV from input dataset, noise amplitude currently ignored
    SSData$data[i] <- data[i]                                       # Add PV to output dataframe for debugging
    
    data_var <- (0.05) * 0.5 * (y[iput] - yold)^2 + 0.95*data_var   # Calculate filtered variance, filter value of 0.05 seleced based on advice from Russell Rhinehart 2019 (OK State) 
    yold          <- y[iput]                                        # Add current PV to yold for next variance calc
    
    y1filt <- ewma * y[iput] + cewma * y1filt                       # Calculate EWMA value at first index
    y2filt <- ewma * y[iread2] + cewma * y2filt                     # Calculate EWMA value at second index
    y3filt <- ewma * y[iread3] + cewma * y3filt                     # Calculate EWMA value at third index
    y4filt <- ewma * y[iread4] + cewma * y4filt                     # Calculate EWMA value at fourth index 
    
    SSData$y1filt[i] <- y1filt                                      # Add PV to output dataframe for debugging
    SSData$y2filt[i] <- y2filt                                      # Add PV to output dataframe for debugging
    SSData$y3filt[i] <- y3filt                                      # Add PV to output dataframe for debugging
    SSData$y4filt[i] <- y4filt                                      # Add PV to output dataframe for debugging
    
    SSData$iread1[i] <- iput                                        # Add PV to output dataframe for debugging
    SSData$iread2[i] <- iread2                                      # Add PV to output dataframe for debugging
    SSData$iread3[i] <- iread3                                      # Add PV to output dataframe for debugging
    SSData$iread4[i] <- iread4                                      # Add PV to output dataframe for debugging
    
    iput   <- iput   + 1                                            # Increment each corner index
    iread2 <- iread2 + 1
    iread3 <- iread3 + 1
    iread4 <- iread4 + 1
    
    if(iput > n) { iput <- 1 }                                      # When an index exceeds the window size, re-initialize the index
    if(iread2 > n) { iread2 <- 1 }
    if(iread3 > n) { iread3 <- 1 }
    if(iread4 > n) { iread4 <- 1 }
    
    maxfilt <- y1filt                                               # Initialize variables to hold min and max filtered values
    minfilt <- y1filt
    
    if (y2filt > maxfilt) { maxfilt = y2filt }                      # Identify minimum/maximum filtered values
    if (y2filt < minfilt) { minfilt = y2filt }
    if (y3filt > maxfilt) { maxfilt = y3filt }
    if (y3filt < minfilt) { minfilt = y3filt }
    if (y4filt > maxfilt) { maxfilt = y4filt }
    if (y4filt < minfilt) { minfilt = y4filt }
    
    if (data_var < 0.01) {data_var <- 0.01}                         # Set minimum variance of 0.01
    t <- (maxfilt - minfilt) / (data_var)^0.5                       # Calculate T-statistic
    
    if (t > 5) {t <- 5}                                             # If T-statistic exceeds 5, force to 5
    
    for (ii in i:(i+stp)) {                                         # Add T-statistic to output dataframe for later troubleshooting
      if(ii<=length(data)) {SSData$tstat[ii] <- t}
    }
    if (i >= n ) {                                                  # If the T-statistic is lower than the threshold value, probable SS
      if (t <= tcrit_l) { SS = 1 }
      if (t > tcrit_u)  { SS = 0 }                                  # If the T-statistic exceeds the upper threshold variable, probable TS
    }
    
    for (ii in i:(i+stp)) { 
      if(ii<=length(data)) {SSData$ss[ii] <- SS }                   # Add steady-state indicator to output dataframe
    }
  }
  
  return(SSData)                                                    # Return steady-state dataframe
}