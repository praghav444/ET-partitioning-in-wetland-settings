library(bigleaf)
library(dplyr)
library(lubridate)
site <- "US-Myb"
ds <- read.csv(paste0("C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/AmeriFlux_Wetland_Data/AMF_US-Myb_FLUXNET_FULLSET_2010-2021_3-5/AMF_US-Myb_FLUXNET_FULLSET_HH_2010-2021_3-5.csv"))
ds[ds==-9999] <- NA
ds$DateTime <- strptime(ds$TIMESTAMP_START, format = "%Y%m%d%H%M")
ds$ET <- LE.to.ET(ds$LE_F_MDS,ds$TA_F_MDS)
ds <- data.frame(DateTime=ds$DateTime,ET=ds$ET,GPP=ds$GPP_NT_VUT_USTAR50)
ds <- ds[complete.cases(ds),]
ds$DateTime <- strptime(ds$DateTime, format = "%Y-%m-%d %H:%M:%S")

#------Let's obtain site specific marginal water use efficiency; m; (ET = m*GPP + E')-----using mean monthly data of ET and GPP----
data_month <- ds %>% group_by(Year=year(DateTime),Month=month(DateTime)) %>%
  summarise_all(funs(mean(.,na.rm=T)))
data_month$ET <- data_month$ET*30*24*3600   # mm/month
data_month$GPP <- umolCO2.to.gC(data_month$GPP)*30    # gC m-2 month-1

lm_fit <- lm(ET~GPP, data=data_month)
m_site <- as.numeric(lm_fit$coefficients[2])

#---Now once we have site specific m; T can be calculated as: T = m*GPP----
ds$ET <- ds$ET*30*24*3600    # mm month-1
ds$GPP <- umolCO2.to.gC(ds$GPP)*30    # gC m-2 month-1
ds$Tr <-  m_site*ds$GPP
# Now to calculate monthly T_ET; calculate monthly mean T and monthly mean ET
ds <- ds %>% group_by(Year=year(DateTime),Month=month(DateTime)) %>%
  summarise_all(funs(mean(.,na.rm=T))) 
ds$T_ET <- ds$Tr/ds$ET
ds$T_ET[ds$T_ET<0 | ds$T_ET>1] <- NA
plot(ds$GPP,ds$ET)

#---save monthly T_ET data----
write.csv(ds,paste0('C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/Scripts/Model_Scott_Biederman/Output/',site,'_output_Scott.csv'), row.names = FALSE)
