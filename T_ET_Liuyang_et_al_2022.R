#----Evapotranspiration Partitioning Based on Leaf and Ecosystem Water Use Efficiency----
#@Ref: Yu, L., Zhou, S., Zhao, X., Gao, X., Jiang, K., Zhang, B., et al. (2022). 
# Evapotranspiration partitioning based on leaf and ecosystem water use efficiency. Water Resources Research, 58, e2021WR030629
# https://doi.org/10.1029/2021WR030629
site <- "US-Myb"
g1s <- c(7.5,8.5)
for(g1 in g1s){
  df <-  read.csv(paste0("C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/AmeriFlux_Wetland_Data/AMF_US-Myb_FLUXNET_FULLSET_2010-2021_3-5/AMF_US-Myb_FLUXNET_FULLSET_HH_2010-2021_3-5.csv"))
  df[df==-9999] <- NA
  ET <- df$LE_F_MDS/(2.501-0.002361*(df$TA_F_MDS))/1000000
  WUE_eco <- df$GPP_NT_VUT_USTAR50/ET/1000000*12
 # g1 <- 3.5    # <----- Change here (Critical parameter)
  WUE_leaf <- (df$CO2_F_MDS*df$PA)/(1.6*100*(df$VPD_F_MDS+g1*(df$VPD_F_MDS)^0.5*10^0.5))
  T_ET <- WUE_eco/WUE_leaf
  T_ET[T_ET < 0 | T_ET > 1] <- NA
  df$T_ET_Liuyang_2022 <- T_ET
  plot(df$T_ET_Liuyang_2022)
  
  df[is.na(df)] <- -9999
  write.csv(df,paste0('C://Users/ppushpendra/OneDrive - The University of Alabama/Manuscript_PhD/Next_Tentative_Project/Scripts/Model_Yu22/Output/',site,'_',g1,'_output_Yu22.csv'),row.names = FALSE)
}
