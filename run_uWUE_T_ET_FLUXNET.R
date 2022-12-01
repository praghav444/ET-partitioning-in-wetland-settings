#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Water Flux Partitioning using underlying Water Use Efficiency (uWUE) Method
# Writen by: Pushpendra Raghav on July 27, 2021
# @ppushpendra@crimson.ua.edu
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# uWUE <- GPP*sqr(VPD)/ET      -------(1)
# T:ET <- uWUEa/uWUEp
# where uWUEa is Apparent uWUE [estimated as the linear regression slope from a 
# moving window spanning either one or eight days (depending on desired 
# smoothening or data availability) or directly from Eqn 1 when estimating at 
# 30 minutes temporal resolution]
# uWUEp is potential uWUE [calculated at annual or seasonal scale using 95th 
# percentile regression between GPP*sqrt(VPD) and ET]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Some important functions

# LE to ET
LE2ET <- function(LE, Ta){
  "
  Convert LE to ET
  
  Parameters
  ---------------
  LE: Latent Heat flux (W m-2)
  Ta: Air temperature (deg C)
  
  Returns
  --------------
  ET: Evapotranspiration (mm s-1)
  "
  lambda <- (2.501 - 0.00237*Ta)*1e6  # Latent heat of vaporization (J kg-1)
  ET <- LE/lambda
  return(ET)
}


calc_PET_PT <- function(Ta, Pa, Rn, G=NA, S=NA, alpha=NA){
  "
  Calculates potential evapotranspiration (PET) according to Priestley & Taylor 1972
  LE_pot = alpha * delta * (Rn - G)) / (delta + gamma)}
  
  Parameters
    ----------
    Ta : Air temperature (deg C)
    Pa : Atmospheric pressure (kPa)
    Rn : Net radiation (W m-2)
    G : Ground heat flux (W m-2); optional
    S : Sum of all storage fluxes (W m-2); optional
    alpha : Priestley-Taylor coefficient (default = 1.26)
    
  Returns
    -------
    PET : Potential evapotranspiration (kg m-2 s-1 or mm/s)
    LE_pot : Potential latent heat flux (W m-2)
    
  
  References
    ----------
    - Priestley, C.H.B., Taylor, R.J., 1972: On the assessment of surface heat flux
      and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.
  "
  G[is.na(G)] <- 0  # Set G to 0 if not provided
  S[is.na(S)] <- 0  # Set S to 0 if not provided
  lambda <- (2.501 - 0.00237*Ta)*1e6  # Latent heat of vaporization (J kg-1)
  Cp <- 1004.834   # specific heat of air for constant pressure (J K-1 kg-1)
  eps <- 0.622     # ratio of the molecular weight of water vapor to dry air (=Mw/Md)
  gamma  <- (Cp * Pa) / (eps * lambda)  # Psychrometric constant (kPa K-1)
  esat <- 0.6108*exp(17.27*Ta/(Ta+237.3))*1e3 # Saturation vapor pressure (Pa)
  delta <- 4098*esat/(237.3+Ta)^2*1e-3 # Slope of the saturation vapor pressure curve (kPa K-1)
  alpha[is.na(alpha)] <- 1.26
  LE_pot = (alpha * delta * (Rn - G - S)) / (delta + gamma)
  PET = LE2ET(LE_pot, Ta)
  return(data.frame(LE_pot, PET))
}

quantile_reg <- function(x, y , PolyDeg=1, tau=0.95, weights){
  "
  Quantile regression
    Fits a polynomial function (of degree PolyDeg) using quantile regression based on a percentile (tau).
    Based on script by Dr. Phillip M. Feldman, and based on method by Koenker, Roger, and
    Gilbert Bassett Jr. Regression Quantiles. Econometrica: Journal of the Econometric Society, 1978, 33-50.
  
  Parameters
    ----------
    x : independent variable
    y : dependent variable
    PolyDeg : Degree of polynomial function
    tau : Percentile for the data to fit to [0-1]
    weights : Vector to weight each point, must be same size as x
    
  Returns
    -------
    The resulting parameters in order of degree from high to low
  
  "
  model <- function(x, beta){
    "
    This example defines the model as a polynomial, where the coefficients of the
    polynomial are passed via `beta`.
    "
    library(signal)
    if(PolyDeg==0){
      return(x*beta)
    } else{
      return(polyval(beta, x))
    }
  }
  N_coefficients <- PolyDeg+1
  
  tilted_abs <- function(tau, x, weights){
    "
     The tilted absolute value function is used in quantile regression.
      INPUTS
       tau: This parameter is a probability, and thus takes values between 0 and 1.
       x: This parameter represents a value of the independent variable, and in
       general takes any real value (float).
    "
    return (weights * x * (tau - (x < 0)))
  }
  objective <- function(beta, tau, weights){
    "
    The objective function to be minimized is the sum of the tilted absolute
    values of the differences between the observations and the model.
    "
    return(sum(tilted_abs(tau, y - model(x, beta), weights)))
  }
  # Build weights if they don't exits:
  if(is.na(weights)){
    weights <- rep(1,length(x))
  }
  # Define starting point for optimization:
  beta_0 <- rep(0, N_coefficients)  
  if(N_coefficients >= 2){
    beta_0[1] <- 1.0
  }
  # `beta_hat[i]` will store the parameter estimates for the quantile
  # corresponding to `fractions[i]`:
  library(stats)
  beta_hat <- optim(beta_0, objective,weights=weights, tau=tau)
  beta_hat <- beta_hat$par
  return(beta_hat)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Data Needed: GPP, VPD, ET, and other supporting variables (if available)
#-------------------------------------------------------------------------------
library(openxlsx)
library(lubridate)
library(dplyr)
library(data.table)
library(quantreg)
library(bigleaf)
#-------------------------------------------------------------------------------
# Creating a `dataframe` having all the required variables
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#files = list.files('/Users/raghav/Library/CloudStorage/Box-Box/SW_FluxNet_NEW/Input_Data/', pattern = "*.csv")
#sites = substr(files,1,6)
#sites <- sites[1:length(sites)]
site <- "US-Sne"
for (i in 1:1){
  ds <- read.csv(paste0("C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/AmeriFlux_Wetland_Data/AMF_US-Sne_FLUXNET_FULLSET_2016-2020_3-5/AMF_US-Sne_FLUXNET_FULLSET_HH_2016-2020_3-5.csv"))
  ds[ds==-9999] <- NA
  ds$DateTime <- strptime(ds$TIMESTAMP_START, format = "%Y%m%d%H%M")
  # Remove duplicates
  ds <- ds[!duplicated(ds[c('DateTime')]),]
  ds <- ds[order(ds$DateTime),]
  ds <- ds[!(is.na(ds$DateTime)),]
  #---First ensure that DateTime is continuous with no gaps in between
  nStepsPerDay <- 24*60/as.numeric(difftime(ds$DateTime[2], ds$DateTime[1], units = 'mins'))
  if(nStepsPerDay == 48){
    time_start <- "00:00:00"
    time_end <- "23:30:00"
  } else {
    time_start <- "00:00:00"
    time_end <- "23:00:00"
  }
  
  date_start <- as.POSIXct(paste0(date(ds$DateTime[1]), time_start))
  date_end <- as.POSIXct(paste0(date(ds$DateTime[nrow(ds)]), time_end))
  seq_DateTime <- seq(from=date_start, by=24*60/nStepsPerDay*60, to=date_end)
  temp <- data.frame(DateTime = seq_DateTime)
  ds <- left_join(temp,ds,by='DateTime')
  if(nrow(ds)%%nStepsPerDay != 0){
    stop('---Killing the processing---Please check the length of the data frame (should be multiple of 24 or 48)----')
  } else{
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Add PET (using P-T 1972)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    ds$LE_pot <- calc_PET_PT(Ta = ds$TA_F_MDS, Pa = ds$PA, Rn = ds$NETRAD, G = ds$G_F_MDS, alpha=1.26)$LE_pot
    ds$PET <- calc_PET_PT(Ta = ds$TA_F_MDS, Pa = ds$PA, Rn = ds$NETRAD, G = ds$G_F_MDS, alpha=1.26)$PET
    
    #----Initial Checks----
    plot(ds$LE_pot, type='l',ylim=c(0,1000))
    lines(ds$LE_F_MDS, type='l', col='red')
    legend('top', legend=c("LE_pot", "LE"),
           col=c("black", "red"), lty=1:1, cex=1.2)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # ----Data screening and quality control----
    # --> (1) Use only half hourly data with high confidence, i.e., original or most reliable data acc. to quality flags (Quality mask)
    # --> (2) only daylight data (i.e., when SW_IN_POT > 0)  with positive Rn, GPP, ET, and VPD (Zero mask)
    # --> (3) Select data during the growing season: i.e., data for days when when average half-hourly 
    # GPP was at least 10% of the 95th percentile of all the half-hourly GPP for the site (SeasonMask)
    # --> (4) exclude data from rainy days and several days that followed rainy days
    # as follows: Two dry days following a rainy day will be excluded when P >= 2*PET, 
    # otherwise one dry day will be excluded when P >= PET. Only the rainy days will be
    # excluded when P < PET for the day. (RainMask)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # 1. ----Quality Mask----
    ds$QualityMask <- TRUE
    ds$NEE_VUT_REF_QC[is.na(ds$NEE_VUT_REF_QC)] <- 0
    ds$LE_F_MDS_QC[is.na(ds$LE_F_MDS_QC)] <- 0
    ds$VPD_F_MDS_QC[is.na(ds$VPD_F_MDS_QC)] <- 0
    ds$TA_F_MDS_QC[is.na(ds$TA_F_MDS_QC)] <- 0
    ds$QualityMask[ds$NEE_VUT_REF_QC >= 2 | ds$NEE_VUT_REF_QC < 0 | !is.finite(ds$NEE_VUT_REF_QC) | !is.finite(ds$GPP_NT_VUT_REF) |
                     ds$LE_F_MDS_QC >= 2 | ds$LE_F_MDS_QC < 0 | !is.finite(ds$LE_F_MDS_QC) | !is.finite(ds$LE_F_MDS) |
                     ds$VPD_F_MDS_QC >= 2 | ds$VPD_F_MDS_QC < 0 | !is.finite(ds$VPD_F_MDS_QC) | !is.finite(ds$VPD_F_MDS) |
                     ds$TA_F_MDS_QC >= 2 | ds$TA_F_MDS_QC < 0 | !is.finite(ds$TA_F_MDS_QC) | !is.finite(ds$TA_F_MDS) |
                     !is.finite(ds$NETRAD)] <- FALSE
    
    # 2. ----Zero Mask ----
    ds$ZeroMask <- TRUE
    ds$ZeroMask[ds$SW_IN_POT <= 0 | ds$LE_F_MDS < 0 | ds$NETRAD < 0 | ds$GPP_NT_VUT_REF < 0 | ds$VPD_F_MDS < 0] <- FALSE
    
    # 3. ----Season Mask----
    # ----RAGHAV----
    #temp <- data.frame(DateTime=ds$DateTime, Ta = ds$TA_F_MDS, SW_IN=ds$SW_IN_F_MDS, Prcp = ds$P)
    #temp <- temp %>% group_by(Month = month(DateTime)) %>%
    #  summarise_all(funs(mean(.,na.rm=T)))
    #temp$Prcp <- temp$Prcp*30*24*3600  # mm/month
    #----
    ds$SeasonMask <- TRUE
    ds <- ds %>% 
      group_by(Date = date(DateTime)) %>% 
      mutate(GPP_day = mean(GPP_NT_VUT_REF, na.rm=T))
    ds$SeasonMask[ds$GPP_day <= 0.10*quantile(ds$GPP_day,0.95, na.rm=T)] <- FALSE
    
    # 4. ----Precipitation Mask----
    ds$PET <- ds$PET*30*60
    ds$PET[is.na(ds$P)] <- NA
    ds$P[is.na(ds$PET)] <- NA
    ds$PrecipMask <- TRUE
    ds <- ds %>% 
      group_by(Date = date(DateTime)) %>% 
      mutate(P_daily = mean(P, na.rm=T)) %>%  # P: kg m2 s-1
      mutate(PET_daily = mean(PET, na.rm=T))  # PET: kg m2 s-1
    ds$PrecipMask[ds$P_daily > 0] <- FALSE
    ds$DateTime <- as.POSIXct(ds$DateTime)  # <------<----- make POSIXct DateTime
    ds <- ds %>% as.data.table() # to simplify operations on ds
    ds <- ds %>% 
      mutate(PrecipMask = case_when(Date %in% c(ds[P_daily > PET_daily, Date],ds[P_daily > PET_daily, Date]+1) ~ FALSE,
                                    TRUE ~PrecipMask ))
    ds <- ds %>% 
      mutate(PrecipMask = case_when(Date %in% c(ds[P_daily > 2*PET_daily, Date],ds[P_daily > 2*PET_daily, Date]+1, 
                                                ds[P_daily > 2*PET_daily, Date]+2) ~ FALSE, TRUE ~PrecipMask ))
    
    #-------------------------------------------------------------------------------
    ds$uWUEa_Mask   = ds$ZeroMask & ds$QualityMask
    ds$uWUEp_Mask   = ds$ZeroMask & ds$QualityMask & ds$PrecipMask & ds$SeasonMask
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Main
    #-------------------------------------------------------------------------------
    ds$GPP <- ds$GPP_NT_VUT_REF*12.001/1e6 * 24*3600 # GPP from umolC m-2 s-1 to gC m-2 day-1
    ds$GPP_mul_sqrt_VPD <- ds$GPP*sqrt(ds$VPD_F_MDS)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    ds$ET <- LE.to.ET(ds$LE_F_MDS,ds$TA_F_MDS)
    ds$ET <- ds$ET * 24*3600  # mm s-1 to mm day-1 or kgH2O m-2 day-1
    ds <- ds[!(is.na(DateTime)),]
    ds$Date <- date(ds$DateTime)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Visualization (Potential underlying water use efficiency (uWUEp))
    #-------------------------------------------------------------------------------
    # Using all the data at a site
    plot(ds$ET[ds$uWUEp_Mask==TRUE], ds$GPP_mul_sqrt_VPD[ds$uWUEp_Mask==TRUE], xlab = expression(ET ~ (kg ~ H[2]*O~m^-2 ~ day^-1)), ylab="", 
         cex.lab=1.5, cex.axis=1.2, xlim = c(0,30), ylim=c(0,250), pch=20)
    title(main=paste0(site), ylab=expression(GPP%.%VPD^0.5~(gC~hPa^0.5~m^-2~ day^-1)), line=2.2, cex.lab=1.5)
    # Fit 95th Quantile regression
    temp_ds <- data.frame(x=ds$ET[ds$uWUEp_Mask==TRUE], y=ds$GPP_mul_sqrt_VPD[ds$uWUEp_Mask==TRUE])
    temp_ds <- temp_ds[complete.cases(temp_ds), ]
    qt_fit <- quantile_reg(temp_ds$x, temp_ds$y,  PolyDeg=0, tau=0.95, weights=NA)
    
    uWUEp_all <- qt_fit[[1]]
    print(uWUEp_all)
    x <- seq(0,16, 0.01)
    y <- uWUEp_all*x #+ qt_fit[[1]][1]
    par(new=TRUE)
    lines(x,y,lty=2,col="red",lwd=3)
    text(10,230, bquote(~ uWUE[p] == .(round(uWUEp_all,2)) ~ (gC%.%hPa^0.5/kgH[2]~O)))
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # ----uWUEp: year-wise----
    years <- unique(year(ds$DateTime))
    uWUEps <- NULL
    for (yr in years) {
      temp_ds <- data.frame(x=ds$ET[year(ds$DateTime)==yr & ds$uWUEp_Mask==TRUE], 
                            y=ds$GPP_mul_sqrt_VPD[year(ds$DateTime)==yr & ds$uWUEp_Mask==TRUE])
      temp_ds <- temp_ds[complete.cases(temp_ds), ]
      qt_fit <- quantile_reg(temp_ds$x, temp_ds$y,  PolyDeg=0, tau=0.95, weights=NA)
      uWUEp <- qt_fit[[1]]
      temp <- data.frame(Year=yr, uWUEp=uWUEp)
      uWUEps <- rbind(uWUEps,temp)
    }
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # ----T:ET at half-hourly scale----
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    ds$uWUEp <- NA
    ds$Year <- year(ds$DateTime)
    for( yr in years){
      ds$uWUEp[ds$Year==yr] <- uWUEps$uWUEp[uWUEps$Year==yr]
    }
    #--- Apparent uWUE ----
    ds$uWUEa <- ds$GPP_mul_sqrt_VPD/ds$ET
    ds$T_ET <- ds$uWUEa/ds$uWUEp
    ds$T_ET[!is.finite(ds$T_ET)] <- NA
    ds$T_ET[ds$T_ET<0] <- NA
    ds$T_ET[ds$T_ET>1] <- NA
    plot(ds$VPD_F_MDS, ds$T_ET, ylim=c(0,1), pch=20)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Save output (half hourly)
    write.csv(ds, paste0('C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/Scripts/Model_uWUE/Output/',
                         site,'_uWUE_output_30min.csv'), row.names = FALSE)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # ----T:ET at daily scale----
    # uWUEa will be estimated as the linear regression slope from a moving window spanning one day
    # T:ET will be estimated only for days when there are at least 10 effective entries
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    temp <- ds
    temp <- temp %>% as.data.table() # to simplify operations on ds
    temp$Day_Mask <- is.finite(ds$ET) & is.finite(ds$GPP_mul_sqrt_VPD)
    temp <- temp %>% 
      group_by(Date) %>%
      mutate(ET_vals_count = case_when(Day_Mask==TRUE ~ sum(!(is.na(ET)) & Day_Mask==TRUE))) %>%  # Count number of effective entries of ET data during DayTime
      mutate(GPP_mul_sqrt_VPD_count = case_when(Day_Mask==TRUE ~ sum(!(is.na(GPP_mul_sqrt_VPD)) & Day_Mask==TRUE))) # Count number of effective entries of GPP*sqrt(VPD) data during DayTime
    
    temp$Day_Mask[(temp$ET_vals_count < 5 | temp$GPP_mul_sqrt_VPD_count < 5)] <- FALSE  # 5: Atleast 5 effective entries during a day
    
    temp <- temp[!(temp$Day_Mask==FALSE),]
    #-------------------------------------------------------------------------------
    # Linear Regression for uWUEa for each day
    temp$uWUEa <- NA
    for(i in 1:length(unique(temp$Date))){
      ET <- temp$ET[temp$Date==unique(temp$Date)[i]]
      GPP_mul_sqrt_VPD <- temp$GPP_mul_sqrt_VPD[temp$Date==unique(temp$Date)[i]]
      if(sum(!(is.na(ET)))>1 & sum(!(is.na(GPP_mul_sqrt_VPD)))>1){
        temp.ds <- data.frame(x=ET, y=GPP_mul_sqrt_VPD)
        temp.ds <- temp.ds[complete.cases(temp.ds),]
        linearMod <- lm(y ~ x, data=temp.ds)
        uWUEa <- linearMod["coefficients"][[1]][2]
        temp$uWUEa[temp$Date==unique(temp$Date)[i]] <- uWUEa
      }
    }
    temp$uWUEp <- NA
    temp$Year <- year(temp$DateTime)
    for( yr in years){
      temp$uWUEp[temp$Year==yr] <- uWUEps$uWUEp[uWUEps$Year==yr]
    }
    temp$T_ET <- temp$uWUEa/temp$uWUEp
    temp$T_ET[temp$T_ET<0 | temp$T_ET>1] <- NA
    temp <- temp %>% group_by(Date=date(DateTime)) %>%
      summarise_each(funs(mean(., na.rm=T)))
    write.csv(temp, paste0('C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/Scripts/Model_uWUE/Output/',
                           site,'_uWUE_output_daily.csv'), row.names = FALSE)
    plot(temp$Date,temp$T_ET, pch=20)
    plot(temp$VPD_F_MDS,temp$T_ET)
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # ----T:ET at daily scale using 8-day moving window----
    # uWUEa will be estimated as the linear regression slope from a moving window spanning 8 day
    # T:ET will be estimated only for days when there are at least 80 effective entries
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    temp <- ds
    temp$uWUEa <- NA
    for(i in 1:length(unique(temp$Date))){
      if(i <= 5){
        istart = 1
        iend = (i+4)*48
      } else if (i > length(unique(temp$Date))-4){
        istart = (i-5)*48
        iend = length(unique(temp$Date))*48
      } else{
        istart = (i-5)*48
        iend = (i+4)*48
      }
      ET <- temp$ET[istart:iend]
      GPP_mul_sqrt_VPD <- temp$GPP_mul_sqrt_VPD[istart:iend]
      if(sum(!(is.na(ET)))>80 & sum(!(is.na(GPP_mul_sqrt_VPD)))>80){ # At least 80 effective entries over 8 days
        temp.ds <- data.frame(x=ET, y=GPP_mul_sqrt_VPD)
        temp.ds <- temp.ds[complete.cases(temp.ds),]
        linearMod <- lm(y ~ x + 0, data=temp.ds)  # With zero intercept
        uWUEa <- linearMod["coefficients"][[1]]
        temp$uWUEa[temp$Date==unique(temp$Date)[i]] <- uWUEa
      }
    }
    temp$uWUEp <- NA
    temp$Year <- year(temp$DateTime)
    for( yr in years){
      temp$uWUEp[temp$Year==yr] <- uWUEps$uWUEp[uWUEps$Year==yr]
    }
    temp$T_ET <- temp$uWUEa/temp$uWUEp
    temp$T_ET[temp$T_ET<0 | temp$T_ET>1] <- NA
    temp <- temp %>% group_by(Date=date(DateTime)) %>%
      summarise_each(funs(mean(., na.rm=T)))
    write.csv(temp, paste0('C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/Scripts/Model_uWUE/Output/',
                           site,'_uWUE_output_daily_8day.csv'), row.names = FALSE)
    plot(temp$Date,temp$T_ET, type='l',pch=20)
    #-------------------------------------------------------------------------------
  }
}