# Chapter 1 Article code

#-------LIBRARIES--------
library(dismo)
library(tidyverse)
library(mgcv)
library(mgcViz)
library(RColorBrewer)
library(sdm)
library(sf)
library(caret)
library(Rcpp)
library(performance)
library(pROC)
library(yardstick)
library(ROCR)
library(gamm4)
library(readxl)
library(scoring)
library(multcompView)
library(gridExtra)
library(grid)
library(ggpubr)
library(dismo)
library(rJava)
library(ggrepel)




#-----DATA------

##-----General albo data----
load("~/Desktop/R/SDM_LI_ALB_FV/Data/AlboTraining.RData")
alboraw <- albotraining
rm(albotraining)
alboraw$Occurrence <- alboraw$Occurrence %>% as.logical()

# collapse down land cover 
alboraw$Type_.25km2 <- alboraw$Type_.25km2 %>% 
  factor(levels = c("Developed - Open","Developed -Low","Developed - Medium",
                    "Developed - High","Cultivated","Forest","Herbaceous","Water",
                    "Wetlands","Shrub", "Barren"), 
         labels = c("Altered","Developed - Low","Developed - Med/High",
                    "Developed - Med/High","Altered","Natural","Natural","Natural",
                    "Natural", "Natural", "Natural"))

albo <- alboraw %>% dplyr::select(Occurrence,Type,Year,Month,
                                  Imp_.25km2,EVI_.25km2,DLST_1km2,NLST_1km2,
                                  Type_.25km2,Longitude,Latitude,Site, Mosquitoes, Township)
# Remove NAs
albo <- albo %>% drop_na(Type,Year,Month,Imp_.25km2,EVI_.25km2,DLST_1km2,NLST_1km2,Type_.25km2,Longitude,Latitude,Site)

# alboME (maxent) make (with type_.25km2 made into numbers for MaxEnt)
alboME <- albo
alboME <- alboME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                 "Natural"=1,
                                                 "Altered"=2,
                                                 "Developed - Low"=3,
                                                 "Developed - Med/High"=4))

##-----General prediction data----
# Load in Prediction dataframe 
load("~/Desktop/R/SDM_LI_ALB_FV/Output/LSTmosquitoPredict.RData")

# add in a type column for to control for TYPE of trap 
LST_Predict$Type <- 1
LST_Predict$Month <- LST_Predict$Month %>% as.numeric()
LST_Predict$Type <- LST_Predict$Type %>% as.factor()

# Remove NAs for DLST and NLST and Type_Site
LST_Predict <- LST_Predict %>% 
  drop_na(Type,Year,Month,Imp_.25km2,EVI_.25km2,DLST_1km2,NLST_1km2,Type_.25km2,coords.x1.x,coords.x2.x,Site)
LST_Predict$Type<-LST_Predict$Type%>%factor(levels=c("1","2","3","4"),labels=c("1","0","BG Sentinel","Mosquito Magnet"))

# I added in shrub and barren!
LST_Predict$Type_.25km2 <- LST_Predict$Type_.25km2 %>% 
  factor(levels = c("Developed - Open","Developed -Low","Developed - Medium",
                    "Developed - High","Cultivated","Forest","Herbaceous","Water",
                    "Wetlands","Shrub", "Barren"), 
         labels = c("Altered","Developed - Low","Developed - Med/High",
                    "Developed - Med/High","Altered","Natural","Natural","Natural",
                    "Natural", "Natural", "Natural"))

LST_Predict <- LST_Predict %>% dplyr::select( Type,Year,Month,Imp_.25km2,EVI_.25km2,DLST_1km2,NLST_1km2,Type_.25km2,coords.x1.x,coords.x2.x,Site)
# averaging lst_predict for each month per site per year
LST_Predict <- LST_Predict %>% 
  group_by(Type, Year, Month, Site, coords.x1.x, coords.x2.x) %>% 
  summarise(Imp_.25km2 = mean(Imp_.25km2), 
            EVI_.25km2 = mean(EVI_.25km2),
            DLST_1km2 = mean(DLST_1km2),
            NLST_1km2 = mean(NLST_1km2),
            Type_.25km2 = first(Type_.25km2)) %>% 
  ungroup()%>%
  data.frame()

LST_Predict <- LST_Predict %>% filter(Year == 2023)

# Make a copy of the data that includes lambda 
LST_Predict_Lambda <- LST_Predict


###-----w/o lambda data for models----
##SDMdata
# 1km level SDMdataset because all variables were aggregated at the same resolution (1km is finest LST res)
sdmAlbo <- sdmData(Occurrence ~ 
                     Type +
                     Year + 
                     Month +
                     Imp_.25km2 + 
                     EVI_.25km2 + 
                     DLST_1km2 +
                     NLST_1km2 +
                     Type_.25km2 +
                     Longitude + 
                     Latitude +
                     Site+
                     coords(Longitude+Latitude), 
                   train=data.frame(albo))

sdmAlboMaxEnt <- sdmData(Occurrence ~
                           Type +
                           Year +
                           Month +
                           Imp_.25km2 +
                           EVI_.25km2 +
                           DLST_1km2 +
                           NLST_1km2 +
                           Type_.25km2 +
                           coords(Longitude+Latitude),
                         train=data.frame(alboME))

###-----w/ lambda data for models----- 

#create functions for mosquito traits
#EFD (eggs laid per female per gonotrophic cycle (number/female))
EFD <- function(T){
  pmax((0.0488)*T*(T-8.02)*(35.65-T)^(1/2),0, na.rm=T)
}
# pEA (mosquito egg-to-adult survival probability)
PEA <- function(T){
  pmax((-0.00361)*(T-9.04)*(T-39.33),0,na.rm=T)
}
# MDR
MDR <- function(T){
  pmax((0.0000638)*T*(T-8.60)*(39.66-T)^(1/2),0,na.rm=T)
}
# mu = 1/ls mosquito adult lifespan (days) 
ls <- function(T){
  pmax((-1.43)*(T-13.41)*(T-31.51),0,na.rm=T)
}
# lambda 
Lambda<-function(MDR, EFD, PEA, ls){
  pmax((EFD*PEA*MDR)/(1/ls),0,na.rm=T)
}

####-----Albo----
# Add lambda to data
albomechanistic <- albo
albomechanistic <- albomechanistic %>% mutate(LST_Mean = ((DLST_1km2+NLST_1km2)/2))
albomechanistic <- albomechanistic %>% mutate(MDR = MDR(LST_Mean),
                                              EFD = EFD(LST_Mean),
                                              PEA = PEA(LST_Mean),
                                              ls = ls(LST_Mean))
albomechanistic <- albomechanistic %>% mutate(Lambda = Lambda(MDR = MDR,EFD = EFD, PEA = PEA, ls = ls))
albomechanistic[,"Lambda"] <- (albomechanistic$Lambda - min(albomechanistic$Lambda))/(max(albomechanistic$Lambda)-min(albomechanistic$Lambda))

# remove extra columns not used in prediction dataframe
albomechanistic <- albomechanistic %>%
  dplyr::select(-LST_Mean, -MDR, -EFD, -PEA, -ls, -Mosquitoes)

# alboME make (with type_.25km2 made into numbers)
alboME_lambda <- albomechanistic
alboME_lambda <- alboME_lambda %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                               "Natural"=1,
                                                               "Altered"=2,
                                                               "Developed - Low"=3,
                                                               "Developed - Med/High"=4))

##SDMdata
# 1km level SDMdataset because all variables were aggregated at the same resolution (1km is finest LST res)
sdmAlbo_Lambda <- sdmData(Occurrence ~ 
                            Type +
                            Year + 
                            Month +
                            Imp_.25km2 + 
                            EVI_.25km2 + 
                            DLST_1km2 +
                            NLST_1km2 +
                            Type_.25km2 +
                            Longitude + 
                            Latitude +
                            Site+
                            Lambda +
                            coords(Longitude+Latitude), 
                          train=data.frame(albomechanistic))

sdmAlboMaxEnt_Lambda <- sdmData(Occurrence ~
                                  Type +
                                  Year +
                                  Month +
                                  Imp_.25km2 +
                                  EVI_.25km2 +
                                  DLST_1km2 +
                                  NLST_1km2 +
                                  Type_.25km2 +
                                  Lambda +
                                  coords(Longitude+Latitude),
                                train=data.frame(alboME_lambda))

####-----Predict----
LST_Predict_Lambda <- LST_Predict_Lambda %>% mutate(LST_Mean = ((DLST_1km2+NLST_1km2)/2))

LST_Predict_Lambda <- LST_Predict_Lambda %>% mutate(MDR = MDR(LST_Mean),
                                                    EFD = EFD(LST_Mean),
                                                    PEA = PEA(LST_Mean),
                                                    ls = ls(LST_Mean))
LST_Predict_Lambda <- LST_Predict_Lambda %>% mutate(Lambda = Lambda(MDR = MDR,EFD = EFD, PEA = PEA, ls = ls))
LST_Predict_Lambda[,"Lambda"] <- (LST_Predict_Lambda$Lambda - min(LST_Predict_Lambda$Lambda))/(max(LST_Predict_Lambda$Lambda)-min(LST_Predict_Lambda$Lambda))

LST_Predict_Lambda <- LST_Predict_Lambda %>%
  dplyr::select(-LST_Mean, -MDR, -EFD, -PEA, -ls)


#-----MODELS w/o lambda -----

##-----------Regression----------
###-----------GAM----------
if (file.exists("./Data/SDM/AlboGAM.RData")) {
  # If the file exists, load it
  load("./Data/SDM/AlboGAM.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  AlboGAM <- gam(Occurrence ~
                   Type +
                   Year +
                   Month +
                   s(Imp_.25km2,k=5) + 
                   s(EVI_.25km2,k=5) +
                   s(DLST_1km2,k=5) +
                   s(NLST_1km2,k=5) +
                   Type_.25km2 +
                   s(Site,bs = "re")
                 ,
                 data = albo,
                 family = binomial)
  summary(AlboGAM)
  save(AlboGAM, file = "./Data/SDM/AlboGAM.RData")
  
  print("File created and saved successfully.")
}


##-----------Machine-learning----------
###-----------RandomForest----------
if (file.exists("./Data/SDM/AlboRF.RData")) {
  # If the file exists, load it
  load("./Data/SDM/AlboRF.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  AlboRF <- randomForest(Occurrence ~ 
                           Type +
                           Year + 
                           Month +
                           Imp_.25km2 + 
                           EVI_.25km2 + 
                           DLST_1km2 +
                           NLST_1km2 +
                           Type_.25km2, 
                         data = albo)
  summary(AlboRF)
  
  save(AlboRF, file = "./Data/SDM/AlboRF.RData")
  
  print("File created and saved successfully.")
}

###-----------BRT----------
if (file.exists("./Data/SDM/AlboBRT.RData")) {
  # If the file exists, load it
  load("./Data/SDM/AlboBRT.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  AlboBRT <- sdm(Occurrence ~ 
                   Type +
                   Year + 
                   Month +
                   Imp_.25km2 + 
                   EVI_.25km2 + 
                   DLST_1km2 +
                   NLST_1km2 +
                   Type_.25km2, 
                 data = sdmAlbo,
                 methods = c('brt'))
  summary(AlboBRT)
  
  save(AlboBRT, file = "./Data/SDM/AlboBRT.RData")
  
  print("File created and saved successfully.")
}

###-----------MaxEnt----------
if (file.exists("./Data/SDM/AlboMaxEnt.RData")) {
  # If the file exists, load it
  load("./Data/SDM/AlboMaxEnt.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  AlboMaxEnt <- sdm(Occurrence ~ 
                      Type +
                      Year + 
                      Month +
                      Imp_.25km2 + 
                      EVI_.25km2 + 
                      DLST_1km2 +
                      NLST_1km2+
                      Type_.25km2, 
                    data = sdmAlboMaxEnt,
                    methods = c('maxent'))
  summary(AlboMaxEnt)
  
  save(AlboMaxEnt, file = "./Data/SDM/AlboMaxEnt.RData")
  
  print("File created and saved successfully.")
}


#-------PREDICTIONS-------

##-----PREDICTIONS w/o lambda -----
# BRT
LST_Predict <- data.frame(LST_Predict,
                          SuitabilityBRT=predict(AlboBRT, newdata = LST_Predict, method = "brt", mean=T))
LST_Predict <- LST_Predict %>% rename(SuitabilityBRT=id_1.sp_1.m_brt)
# following function is for scaling suitability values as they are only meaningful in relationship to each oth and its good viewing to see them all on the same scale 
LST_Predict$SuitabilityBRT <- (LST_Predict$SuitabilityBRT - min(LST_Predict$SuitabilityBRT))/(max(LST_Predict$SuitabilityBRT)-min(LST_Predict$SuitabilityBRT))

# random forest (using random forest package)
SuitabilityRF <- predict(AlboRF, newdata = LST_Predict)
LST_Predict <- data.frame(LST_Predict,
                          SuitabilityRF=SuitabilityRF)
LST_Predict$SuitabilityRF <- (LST_Predict$SuitabilityRF - min(LST_Predict$SuitabilityRF))/(max(LST_Predict$SuitabilityRF)-min(LST_Predict$SuitabilityRF))

# GAM
SuitabilityGAM <- predict(AlboGAM, newdata = LST_Predict,exclude="s(Site)",type="response")
LST_Predict <- data.frame(LST_Predict,
                          SuitabilityGAM=SuitabilityGAM)
LST_Predict$SuitabilityGAM <- (LST_Predict$SuitabilityGAM - min(LST_Predict$SuitabilityGAM))/(max(LST_Predict$SuitabilityGAM)-min(LST_Predict$SuitabilityGAM))

# MAXENT (had to get rid of Type_.25km2)
LST_PredictMaxEnt <- LST_Predict %>%  dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
LST_PredictMaxEnt <- LST_PredictMaxEnt %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                       "Natural"=1,
                                                                       "Altered"=2,
                                                                       "Developed - Low"=3,
                                                                       "Developed - Med/High"=4))
SuitabilityMaxEnt <- predict(AlboMaxEnt,LST_PredictMaxEnt, methods='maxent')
LST_Predict <- data.frame(LST_Predict,
                          SuitabilityMaxEnt=SuitabilityMaxEnt)
LST_Predict <- LST_Predict %>% rename(SuitabilityMaxEnt=id_1.sp_1.m_maxent)
LST_Predict$SuitabilityMaxEnt <- (LST_Predict$SuitabilityMaxEnt - min(LST_Predict$SuitabilityMaxEnt))/(max(LST_Predict$SuitabilityMaxEnt)-min(LST_Predict$SuitabilityMaxEnt))

rm(SuitabilityGAM, SuitabilityMaxEnt, SuitabilityRF)


#-------CROSS-VALIDATION-------
# ---------Township Subsampling (10 townships, hold one out for testing)---------
if (file.exists("/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_township.RData")) {
  # If the file exists, load it
  load("/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_township.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  Sys.time()
  fit_models_township<-lapply(1:10,function(i){
    # split albo into training and testing data (for GAM and RF) - by site (not datapoint)
    township <- data.frame(unique(albo$Township))
    
    testing_dataset_albo <- albo %>% filter(Township%in%township$unique.albo.Township.[i])
    training_dataset_albo <- albo %>% filter(!Township%in%township$unique.albo.Township.[i])
    
    
    # make datasets for BRT based on data partitioned above (can use same testing data as RF and GAM)
    training_dataset_albo_BRT <- sdmData(Occurrence ~ 
                                           Type +
                                           Year + 
                                           Month +
                                           Imp_.25km2 + 
                                           EVI_.25km2 + 
                                           DLST_1km2 +
                                           NLST_1km2 +
                                           Type_.25km2 +
                                           Longitude + 
                                           Latitude +
                                           Site+
                                           coords(Longitude+Latitude), 
                                         train=data.frame(training_dataset_albo))
    
    # Make datasets for MaxEnt Based on Data partitioned above
    # alboME make (with type_.25km2 made into numbers)
    training_dataset_albo_ME <- training_dataset_albo
    training_dataset_albo_ME <- training_dataset_albo_ME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                                         "Natural"=1,
                                                                                         "Altered"=2,
                                                                                         "Developed - Low"=3,
                                                                                         "Developed - Med/High"=4))
    training_dataset_albo_ME_m <- sdmData(Occurrence ~
                                            Type +
                                            Year +
                                            Month +
                                            Imp_.25km2 +
                                            EVI_.25km2 +
                                            DLST_1km2 +
                                            NLST_1km2 +
                                            Type_.25km2 +
                                            coords(Longitude+Latitude),
                                          train=data.frame(training_dataset_albo_ME))
    
    testing_dataset_albo_ME <- testing_dataset_albo
    testing_dataset_albo_ME <- testing_dataset_albo_ME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                                       "Natural"=1,
                                                                                       "Altered"=2,
                                                                                       "Developed - Low"=3,
                                                                                       "Developed - Med/High"=4))
    testing_dataset_albo_ME <- testing_dataset_albo_ME %>% 
      dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
    
    training_dataset_albo_ME <- training_dataset_albo_ME %>%
      dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
    
    tryCatch({ ##This is done so that it doesn't stop everytime something throws an error
      print(paste0(i,": Fitting Models"))
      # GAM and Random Forest
      Albo_GAM_.25 <- gam(Occurrence ~
                            Type +
                            Year +
                            Month +
                            s(Imp_.25km2,k=5) + 
                            s(EVI_.25km2,k=5) +
                            s(DLST_1km2,k=5) +
                            s(NLST_1km2,k=5) +
                            Type_.25km2 +
                            s(Site,bs = "re"),
                          data = training_dataset_albo,
                          family = binomial)
      
      Albo_RF_.25 <- randomForest(Occurrence ~ 
                                    Type +
                                    Year + 
                                    Month +
                                    Imp_.25km2 + 
                                    EVI_.25km2 + 
                                    DLST_1km2 +
                                    NLST_1km2 +
                                    Type_.25km2, 
                                  data = training_dataset_albo)
      
      # Boosted Regression Tree (BRT)
      Albo_BRT_.25 <- sdm(Occurrence ~ 
                            Type +
                            Year + 
                            Month +
                            Imp_.25km2 + 
                            EVI_.25km2 + 
                            DLST_1km2 +
                            NLST_1km2 +
                            Type_.25km2, 
                          data = training_dataset_albo_BRT,
                          methods = c('brt'))
      
      # MaxEnt
      Albo_MaxEnt_.25 <- sdm(Occurrence ~ 
                               Type +
                               Year + 
                               Month +
                               Imp_.25km2 + 
                               EVI_.25km2 + 
                               DLST_1km2 +
                               NLST_1km2+
                               Type_.25km2, 
                             data = training_dataset_albo_ME_m,
                             methods = c('maxent'))
      
      #Predict outcome for the testing data
      print(paste0(i,": Predictions"))
      
      pred_Albo_GAM_.25_testing <- predict(Albo_GAM_.25, newdata = testing_dataset_albo,exclude="s(Site)",type="response")
      pred_Albo_GAM_.25_training <- predict(Albo_GAM_.25, type="response",keep_prediction_data=T)
      
      pred_Albo_RF_.25_testing <- predict(Albo_RF_.25, newdata = testing_dataset_albo)
      pred_Albo_RF_.25_training <- predict(Albo_RF_.25)
      
      pred_Albo_BRT_.25_testing <- predict(Albo_BRT_.25, newdata = testing_dataset_albo, method = "brt", mean=T)%>% data.frame() %>% rename(pred_Albo_BRT_.25_testing=id_1.sp_1.m_brt)
      pred_Albo_BRT_.25_training <- predict(Albo_BRT_.25, newdata = training_dataset_albo, method = "brt", mean=T)%>% data.frame() %>% rename(pred_Albo_BRT_.25_training=id_1.sp_1.m_brt)
      
      pred_Albo_MaxEnt_.25_testing <- predict(Albo_MaxEnt_.25, newdata=testing_dataset_albo_ME, method='maxent')%>% data.frame() %>% rename(pred_Albo_MaxEnt_.25_testing=id_1.sp_1.m_maxent)
      pred_Albo_MaxEnt_.25_training <- predict(Albo_MaxEnt_.25, newdata=training_dataset_albo_ME, method='maxent')%>% data.frame() %>% rename(pred_Albo_MaxEnt_.25_training=id_1.sp_1.m_maxent)
      
      #Create dataframes of predictions and test data and save into respective lists
      #Names columns in dataframes with run number
      pred_df_.25 <- data.frame(Site=testing_dataset_albo$Site, Index= 1:dim(pred_Albo_GAM_.25_testing)[1],
                                Type="OOS",
                                pred_Albo_GAM_.25=pred_Albo_GAM_.25_testing,
                                pred_Albo_RF_.25=pred_Albo_RF_.25_testing,
                                pred_Albo_BRT_.25=pred_Albo_BRT_.25_testing,
                                pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_testing,
                                Occurrence=testing_dataset_albo$Occurrence) %>% rename(pred_Albo_BRT_.25=pred_Albo_BRT_.25_testing, pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_testing)
      
      pred_df_.25_training <- data.frame(Site=Albo_GAM_.25$model$Site, Index= 1:dim(pred_Albo_GAM_.25_training)[1],
                                         Type="WS",
                                         pred_Albo_GAM_.25=pred_Albo_GAM_.25_training,
                                         pred_Albo_RF_.25=pred_Albo_RF_.25_training,
                                         pred_Albo_BRT_.25=pred_Albo_BRT_.25_training,
                                         pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_training,
                                         Occurrence=Albo_GAM_.25$model$Occurrence)%>% rename(pred_Albo_BRT_.25=pred_Albo_BRT_.25_training, pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_training)
      
      pred_df_.25 <- rbind(pred_df_.25,pred_df_.25_training)
      
      # Scale suitability
      pred_df_.25[,"pred_Albo_GAM_.25"] <- (pred_df_.25$pred_Albo_GAM_.25 - min(pred_df_.25$pred_Albo_GAM_.25))/(max(pred_df_.25$pred_Albo_GAM_.25)-min(pred_df_.25$pred_Albo_GAM_.25))
      pred_df_.25[,"pred_Albo_RF_.25"] <- (pred_df_.25$pred_Albo_RF_.25 - min(pred_df_.25$pred_Albo_RF_.25))/(max(pred_df_.25$pred_Albo_RF_.25)-min(pred_df_.25$pred_Albo_RF_.25))
      pred_df_.25[,"pred_Albo_BRT_.25"] <- (pred_df_.25$pred_Albo_BRT_.25 - min(pred_df_.25$pred_Albo_BRT_.25))/(max(pred_df_.25$pred_Albo_BRT_.25)-min(pred_df_.25$pred_Albo_BRT_.25))
      pred_df_.25[,"pred_Albo_MaxEnt_.25"] <- (pred_df_.25$pred_Albo_MaxEnt_.25 - min(pred_df_.25$pred_Albo_MaxEnt_.25))/(max(pred_df_.25$pred_Albo_MaxEnt_.25)-min(pred_df_.25$pred_Albo_MaxEnt_.25))
      
      
      pred_df_.25 <- pred_df_.25 %>% 
        mutate(Run = i)
      
      #Create prediction objects for performance
      temp <- pred_df_.25%>%
        filter(Type=="OOS")
      
      Albo_GAM_.25_pred_test<-prediction(temp$pred_Albo_GAM_.25,temp$Occurrence)
      Albo_RF_.25_pred_test <- prediction(temp$pred_Albo_RF_.25,temp$Occurrence)
      Albo_BRT_.25_pred_test <- prediction(temp$pred_Albo_BRT_.25,temp$Occurrence)
      Albo_MaxEnt_.25_pred_test <- prediction(temp$pred_Albo_MaxEnt_.25,temp$Occurrence)
      
      temp <- pred_df_.25%>%
        filter(Type=="WS")
      Albo_GAM_.25_pred_train<-prediction(temp$pred_Albo_GAM_.25,temp$Occurrence)
      Albo_RF_.25_pred_train <- prediction(temp$pred_Albo_RF_.25,temp$Occurrence)
      Albo_BRT_.25_pred_train <- prediction(temp$pred_Albo_BRT_.25,temp$Occurrence)
      Albo_MaxEnt_.25_pred_train <- prediction(temp$pred_Albo_MaxEnt_.25,temp$Occurrence)
      
      #Return a list of needed objects
      r_list_.25<-list(pred_df_.25,
                       Albo_GAM_.25_pred_test,
                       Albo_RF_.25_pred_test,
                       Albo_BRT_.25_pred_test,
                       Albo_MaxEnt_.25_pred_test,
                       Albo_GAM_.25_pred_train,
                       Albo_RF_.25_pred_train,
                       Albo_BRT_.25_pred_train,
                       Albo_MaxEnt_.25_pred_train
      )
      # names(r_list_.25)<-c("Pred_DF",
      #                      "GAM_Pred","RF_Pred","BRT_Pred","MaxEnt_Pred")
      names(r_list_.25)<-c(paste0("Pred_DF: ",i),
                           paste0("GAM_Pred_test: ",i),
                           paste0("RF_Pred_test: ",i),
                           paste0("BRT_Pred_test: ",i),
                           paste0("MaxEnt_Pred_test: ",i),
                           paste0("GAM_Pred_train: ",i),
                           paste0("RF_Pred_train: ",i),
                           paste0("BRT_Pred_train: ",i),
                           paste0("MaxEnt_Pred_train: ",i)
      )
      
      
      return(r_list_.25)
    },
    error=function(cond){ #This tells it what to do if there is an error
      return(NULL)        #This just says to return NULL and go to the next iteration
    })
  })%>%unlist(recursive = F)
  Sys.time()
  
  save(fit_models_township, file="/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_township.RData")
  
  
  print("File created and saved successfully.")
}


# ---------Random Subsampling (80/20%, train/test)---------

#Set number of model runs we want
runs=50
# Set seed for reproducibility of randomness  
set.seed(6)

if (file.exists("/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_random.RData")) {
  # If the file exists, load it
  load("/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_random.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  Sys.time()
  fit_models<-lapply(1:runs,function(i){
    # split albo into training and testing data (for GAM and RF) - by site (not datapoint)
    training_sites <- slice_sample(data.frame(unique(albo$Site)),prop=.8) 
    training_dataset_albo <- albo %>% filter(Site%in%training_sites$unique.albo.Site.)
    testing_dataset_albo <- albo %>% filter(!Site%in%training_sites$unique.albo.Site.)
    
    # make datasets for BRT based on data partitioned above (can use same testing data as RF and GAM)
    training_dataset_albo_BRT <- sdmData(Occurrence ~
                                           Type +
                                           Year +
                                           Month +
                                           Imp_.25km2 +
                                           EVI_.25km2 +
                                           DLST_1km2 +
                                           NLST_1km2 +
                                           Type_.25km2 +
                                           Longitude +
                                           Latitude +
                                           Site+
                                           coords(Longitude+Latitude),
                                         train=data.frame(training_dataset_albo))
    
    # Make datasets for MaxEnt Based on Data partitioned above
    # alboME make (with type_.25km2 made into numbers)
    training_dataset_albo_ME <- training_dataset_albo
    training_dataset_albo_ME <- training_dataset_albo_ME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                                         "Natural"=1,
                                                                                         "Altered"=2,
                                                                                         "Developed - Low"=3,
                                                                                         "Developed - Med/High"=4))
    training_dataset_albo_ME_m <- sdmData(Occurrence ~
                                            Type +
                                            Year +
                                            Month +
                                            Imp_.25km2 +
                                            EVI_.25km2 +
                                            DLST_1km2 +
                                            NLST_1km2 +
                                            Type_.25km2 +
                                            coords(Longitude+Latitude),
                                          train=data.frame(training_dataset_albo_ME))
    
    testing_dataset_albo_ME <- testing_dataset_albo
    testing_dataset_albo_ME <- testing_dataset_albo_ME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                                       "Natural"=1,
                                                                                       "Altered"=2,
                                                                                       "Developed - Low"=3,
                                                                                       "Developed - Med/High"=4))
    testing_dataset_albo_ME <- testing_dataset_albo_ME %>%
      dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
    
    training_dataset_albo_ME <- training_dataset_albo_ME %>%
      dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
    
    tryCatch({ ##This is done so that it doesn't stop everytime something throws an error
      print(paste0(i,": Fitting Models"))
      # GAM and Random Forest
      Albo_GAM_.25 <- gam(Occurrence ~
                            Type +
                            Year +
                            Month +
                            s(Imp_.25km2,k=5) + 
                            s(EVI_.25km2,k=5) +
                            s(DLST_1km2,k=5) +
                            s(NLST_1km2,k=5) +
                            Type_.25km2 +
                            s(Site,bs = "re"),
                          data = training_dataset_albo,
                          family = binomial)
      
      Albo_RF_.25 <- randomForest(Occurrence ~ 
                                    Type +
                                    Year + 
                                    Month +
                                    Imp_.25km2 + 
                                    EVI_.25km2 + 
                                    DLST_1km2 +
                                    NLST_1km2 +
                                    Type_.25km2, 
                                  data = training_dataset_albo)
      
      # Boosted Regression Tree (BRT)
      Albo_BRT_.25 <- sdm(Occurrence ~ 
                            Type +
                            Year + 
                            Month +
                            Imp_.25km2 + 
                            EVI_.25km2 + 
                            DLST_1km2 +
                            NLST_1km2 +
                            Type_.25km2, 
                          data = training_dataset_albo_BRT,
                          methods = c('brt'))
      
      # MaxEnt
      Albo_MaxEnt_.25 <- sdm(Occurrence ~ 
                               Type +
                               Year + 
                               Month +
                               Imp_.25km2 + 
                               EVI_.25km2 + 
                               DLST_1km2 +
                               NLST_1km2+
                               Type_.25km2, 
                             data = training_dataset_albo_ME_m,
                             methods = c('maxent'))
      
      
      #Predict outcome for the testing data
      print(paste0(i,": Predictions"))
      
      pred_Albo_GAM_.25_testing <- predict(Albo_GAM_.25, newdata = testing_dataset_albo,exclude="s(Site)",type="response")
      pred_Albo_GAM_.25_training <- predict(Albo_GAM_.25, type="response",keep_prediction_data=T)
      
      pred_Albo_RF_.25_testing <- predict(Albo_RF_.25, newdata = testing_dataset_albo)
      pred_Albo_RF_.25_training <- predict(Albo_RF_.25)
      
      pred_Albo_BRT_.25_testing <- predict(Albo_BRT_.25, newdata = testing_dataset_albo, method = "brt", mean=T)%>% data.frame() %>% rename(pred_Albo_BRT_.25_testing=id_1.sp_1.m_brt)
      pred_Albo_BRT_.25_training <- predict(Albo_BRT_.25, newdata = training_dataset_albo, method = "brt", mean=T)%>% data.frame() %>% rename(pred_Albo_BRT_.25_training=id_1.sp_1.m_brt)
      
      pred_Albo_MaxEnt_.25_testing <- predict(Albo_MaxEnt_.25, newdata=testing_dataset_albo_ME, method='maxent')%>% data.frame() %>% rename(pred_Albo_MaxEnt_.25_testing=id_1.sp_1.m_maxent)
      pred_Albo_MaxEnt_.25_training <- predict(Albo_MaxEnt_.25, newdata=training_dataset_albo_ME, method='maxent')%>% data.frame() %>% rename(pred_Albo_MaxEnt_.25_training=id_1.sp_1.m_maxent)
      
      #Create dataframes of predictions and test data and save into respective lists
      #Names columns in dataframes with run number
      pred_df_.25 <- data.frame(Site=testing_dataset_albo$Site, Index= 1:dim(pred_Albo_GAM_.25_testing)[1],
                                Type="OOS",
                                pred_Albo_GAM_.25=pred_Albo_GAM_.25_testing,
                                pred_Albo_RF_.25=pred_Albo_RF_.25_testing,
                                pred_Albo_BRT_.25=pred_Albo_BRT_.25_testing,
                                pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_testing,
                                Occurrence=testing_dataset_albo$Occurrence) %>% rename(pred_Albo_BRT_.25=pred_Albo_BRT_.25_testing, pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_testing)
      
      pred_df_.25_training <- data.frame(Site=Albo_GAM_.25$model$Site, Index= 1:dim(pred_Albo_GAM_.25_training)[1],
                                         Type="WS",
                                         pred_Albo_GAM_.25=pred_Albo_GAM_.25_training,
                                         pred_Albo_RF_.25=pred_Albo_RF_.25_training,
                                         pred_Albo_BRT_.25=pred_Albo_BRT_.25_training,
                                         pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_training,
                                         Occurrence=Albo_GAM_.25$model$Occurrence)%>% rename(pred_Albo_BRT_.25=pred_Albo_BRT_.25_training, pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_training)
      
      pred_df_.25 <- rbind(pred_df_.25,pred_df_.25_training)
      
      # Scale suitability
      pred_df_.25[,"pred_Albo_GAM_.25"] <- (pred_df_.25$pred_Albo_GAM_.25 - min(pred_df_.25$pred_Albo_GAM_.25))/(max(pred_df_.25$pred_Albo_GAM_.25)-min(pred_df_.25$pred_Albo_GAM_.25))
      pred_df_.25[,"pred_Albo_RF_.25"] <- (pred_df_.25$pred_Albo_RF_.25 - min(pred_df_.25$pred_Albo_RF_.25))/(max(pred_df_.25$pred_Albo_RF_.25)-min(pred_df_.25$pred_Albo_RF_.25))
      pred_df_.25[,"pred_Albo_BRT_.25"] <- (pred_df_.25$pred_Albo_BRT_.25 - min(pred_df_.25$pred_Albo_BRT_.25))/(max(pred_df_.25$pred_Albo_BRT_.25)-min(pred_df_.25$pred_Albo_BRT_.25))
      pred_df_.25[,"pred_Albo_MaxEnt_.25"] <- (pred_df_.25$pred_Albo_MaxEnt_.25 - min(pred_df_.25$pred_Albo_MaxEnt_.25))/(max(pred_df_.25$pred_Albo_MaxEnt_.25)-min(pred_df_.25$pred_Albo_MaxEnt_.25))
      
      
      pred_df_.25 <- pred_df_.25 %>% 
        mutate(Run = i)
      
      #Create prediction objects for performance
      temp <- pred_df_.25%>%
        filter(Type=="OOS")
      
      Albo_GAM_.25_pred_test<-prediction(temp$pred_Albo_GAM_.25,temp$Occurrence)
      Albo_RF_.25_pred_test <- prediction(temp$pred_Albo_RF_.25,temp$Occurrence)
      Albo_BRT_.25_pred_test <- prediction(temp$pred_Albo_BRT_.25,temp$Occurrence)
      Albo_MaxEnt_.25_pred_test <- prediction(temp$pred_Albo_MaxEnt_.25,temp$Occurrence)
      
      temp <- pred_df_.25%>%
        filter(Type=="WS")
      Albo_GAM_.25_pred_train<-prediction(temp$pred_Albo_GAM_.25,temp$Occurrence)
      Albo_RF_.25_pred_train <- prediction(temp$pred_Albo_RF_.25,temp$Occurrence)
      Albo_BRT_.25_pred_train <- prediction(temp$pred_Albo_BRT_.25,temp$Occurrence)
      Albo_MaxEnt_.25_pred_train <- prediction(temp$pred_Albo_MaxEnt_.25,temp$Occurrence)
      
      #Return a list of needed objects
      r_list_.25<-list(pred_df_.25,
                       Albo_GAM_.25_pred_test,
                       Albo_RF_.25_pred_test,
                       Albo_BRT_.25_pred_test,
                       Albo_MaxEnt_.25_pred_test,
                       Albo_GAM_.25_pred_train,
                       Albo_RF_.25_pred_train,
                       Albo_BRT_.25_pred_train,
                       Albo_MaxEnt_.25_pred_train
      )
      # names(r_list_.25)<-c("Pred_DF",
      #                      "GAM_Pred","RF_Pred","BRT_Pred","MaxEnt_Pred")
      names(r_list_.25)<-c(paste0("Pred_DF: ",i),
                           paste0("GAM_Pred_test: ",i),
                           paste0("RF_Pred_test: ",i),
                           paste0("BRT_Pred_test: ",i),
                           paste0("MaxEnt_Pred_test: ",i),
                           paste0("GAM_Pred_train: ",i),
                           paste0("RF_Pred_train: ",i),
                           paste0("BRT_Pred_train: ",i),
                           paste0("MaxEnt_Pred_train: ",i)
      )
      
      
      return(r_list_.25)
    },
    error=function(cond){ #This tells it what to do if there is an error
      return(NULL)        #This just says to return NULL and go to the next iteration
    })
  })%>%unlist(recursive = F)
  Sys.time()
  save(fit_models, file="/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_random.RData")
  
  print("File created and saved successfully.")
}
rm(runs)

# ---------AUC values---------

##---------AUC Township Sampling---------
AUCs_township<-lapply(seq(2,(90),by=9),function(i){
  # Test
  temp_GAM_test<-performance(fit_models_township[[i]],measure = "auc", x.measure = "cutoff")
  AUC_GAM_test<-temp_GAM_test@y.values[[1]]
  temp_RF_test<-performance(fit_models_township[[i+1]],measure = "auc", x.measure = "cutoff")
  AUC_RF_test<-temp_RF_test@y.values[[1]]
  temp_BRT_test<-performance(fit_models_township[[i+2]],measure = "auc", x.measure = "cutoff")
  AUC_BRT_test<-temp_BRT_test@y.values[[1]]
  temp_MaxEnt_test<-performance(fit_models_township[[i+3]],measure = "auc", x.measure = "cutoff")
  AUC_MaxEnt_test<-temp_MaxEnt_test@y.values[[1]]
  # train
  temp_GAM_train<-performance(fit_models_township[[i+4]],measure = "auc", x.measure = "cutoff")
  AUC_GAM_train<-temp_GAM_train@y.values[[1]]
  temp_RF_train<-performance(fit_models_township[[i+5]],measure = "auc", x.measure = "cutoff")
  AUC_RF_train<-temp_RF_train@y.values[[1]]
  temp_BRT_train<-performance(fit_models_township[[i+6]],measure = "auc", x.measure = "cutoff")
  AUC_BRT_train<-temp_BRT_train@y.values[[1]]
  temp_MaxEnt_train<-performance(fit_models_township[[i+7]],measure = "auc", x.measure = "cutoff")
  AUC_MaxEnt_train<-temp_MaxEnt_train@y.values[[1]]
  temp<-data.frame(GAM_test=AUC_GAM_test,
                   RF_test=AUC_RF_test,
                   BRT_test=AUC_BRT_test, 
                   MaxEnt_test=AUC_MaxEnt_test,
                   GAM_train=AUC_GAM_train,
                   RF_train=AUC_RF_train,
                   BRT_train=AUC_BRT_train, 
                   MaxEnt_train=AUC_MaxEnt_train)%>%
    mutate(Best=names(which.max(.[,1:4])))    ##This will determine which model has the best AUC for each iter
  return(temp)
})%>%
  bind_rows()
AUCs_township$Best<-AUCs_township$Best%>%factor()

##---------AUC Random Sampling---------
AUCs_random<-lapply(seq(2,(450),by=9),function(i){
  # Test
  temp_GAM_test<-performance(fit_models[[i]],measure = "auc", x.measure = "cutoff")
  AUC_GAM_test<-temp_GAM_test@y.values[[1]]
  temp_RF_test<-performance(fit_models[[i+1]],measure = "auc", x.measure = "cutoff")
  AUC_RF_test<-temp_RF_test@y.values[[1]]
  temp_BRT_test<-performance(fit_models[[i+2]],measure = "auc", x.measure = "cutoff")
  AUC_BRT_test<-temp_BRT_test@y.values[[1]]
  temp_MaxEnt_test<-performance(fit_models[[i+3]],measure = "auc", x.measure = "cutoff")
  AUC_MaxEnt_test<-temp_MaxEnt_test@y.values[[1]]
  # train
  temp_GAM_train<-performance(fit_models[[i+4]],measure = "auc", x.measure = "cutoff")
  AUC_GAM_train<-temp_GAM_train@y.values[[1]]
  temp_RF_train<-performance(fit_models[[i+5]],measure = "auc", x.measure = "cutoff")
  AUC_RF_train<-temp_RF_train@y.values[[1]]
  temp_BRT_train<-performance(fit_models[[i+6]],measure = "auc", x.measure = "cutoff")
  AUC_BRT_train<-temp_BRT_train@y.values[[1]]
  temp_MaxEnt_train<-performance(fit_models[[i+7]],measure = "auc", x.measure = "cutoff")
  AUC_MaxEnt_train<-temp_MaxEnt_train@y.values[[1]]
  temp<-data.frame(GAM_test=AUC_GAM_test,
                   RF_test=AUC_RF_test,
                   BRT_test=AUC_BRT_test, 
                   MaxEnt_test=AUC_MaxEnt_test,
                   GAM_train=AUC_GAM_train,
                   RF_train=AUC_RF_train,
                   BRT_train=AUC_BRT_train, 
                   MaxEnt_train=AUC_MaxEnt_train)%>%
    mutate(Best=names(which.max(.[,1:4])))    ##This will determine which model has the best AUC for each iter
  return(temp)
})%>%
  bind_rows()
AUCs_random$Best<-AUCs_random$Best%>%factor()

#--------ROC Datasets--------


##------ROC Township Sampling----------------
ROCs_township<-lapply(seq(2,(90),by=9),function(i){
  # test
  temp_GAM_test<-performance(fit_models_township[[i]],measure = "tpr", x.measure = "fpr")
  temp_RF_test<-performance(fit_models_township[[i+1]],measure = "tpr", x.measure = "fpr")
  temp_BRT_test<-performance(fit_models_township[[i+2]],measure = "tpr", x.measure = "fpr")
  temp_MaxEnt_test<-performance(fit_models_township[[i+3]],measure = "tpr", x.measure = "fpr")
  # train
  temp_GAM_train<-performance(fit_models_township[[i+4]],measure = "tpr", x.measure = "fpr")
  temp_RF_train<-performance(fit_models_township[[i+5]],measure = "tpr", x.measure = "fpr")
  temp_BRT_train<-performance(fit_models_township[[i+6]],measure = "tpr", x.measure = "fpr")
  temp_MaxEnt_train<-performance(fit_models_township[[i+7]],measure = "tpr", x.measure = "fpr")
  
  # test
  temp_GAM_df_test<-data.frame(Run=i/5,Model="GAM",FPR=temp_GAM_test@x.values[[1]],TPR=temp_GAM_test@y.values[[1]]) %>% mutate(Data = "Test")
  temp_RF_df_test<-data.frame(Run=i/5,Model="RF",FPR=temp_RF_test@x.values[[1]],TPR=temp_RF_test@y.values[[1]])%>% mutate(Data = "Test")
  temp_BRT_df_test<-data.frame(Run=i/5,Model="BRT",FPR=temp_BRT_test@x.values[[1]],TPR=temp_BRT_test@y.values[[1]])%>% mutate(Data = "Test")
  temp_MaxEnt_df_test<-data.frame(Run=i/5,Model="MaxEnt",FPR=temp_MaxEnt_test@x.values[[1]],TPR=temp_MaxEnt_test@y.values[[1]])%>% mutate(Data = "Test")
  # train
  temp_GAM_df_train<-data.frame(Run=i/5,Model="GAM",FPR=temp_GAM_train@x.values[[1]],TPR=temp_GAM_train@y.values[[1]])%>% mutate(Data = "Train")
  temp_RF_df_train<-data.frame(Run=i/5,Model="RF",FPR=temp_RF_train@x.values[[1]],TPR=temp_RF_train@y.values[[1]])%>% mutate(Data = "Train")
  temp_BRT_df_train<-data.frame(Run=i/5,Model="BRT",FPR=temp_BRT_train@x.values[[1]],TPR=temp_BRT_train@y.values[[1]])%>% mutate(Data = "Train")
  temp_MaxEnt_df_train<-data.frame(Run=i/5,Model="MaxEnt",FPR=temp_MaxEnt_train@x.values[[1]],TPR=temp_MaxEnt_train@y.values[[1]])%>% mutate(Data = "Train")
  
  
  return(bind_rows(temp_GAM_df_test,temp_RF_df_test,temp_BRT_df_test,temp_MaxEnt_df_test,
                   temp_GAM_df_train,temp_RF_df_train,temp_BRT_df_train,temp_MaxEnt_df_train))
})%>%
  bind_rows()

ROCs_township$Run<-factor(ROCs_township$Run)
ROCs_township$Model<-factor(ROCs_township$Model,
                            levels=c("GAM","RF","BRT", "MaxEnt"),
                            labels = c("GAM","RF","BRT", "MaxEnt"))
ROCs_township$Data<-factor(ROCs_township$Data)
ROCs_township<-ROCs_township%>%mutate(Buffer=".25km2")
ROCs_township<-ROCs_township%>%mutate(Lambda="No")


##------ROC random Sampling----------------
ROCs_random<-lapply(seq(2,(450),by=9),function(i){
  # test
  temp_GAM_test<-performance(fit_models[[i]],measure = "tpr", x.measure = "fpr")
  temp_RF_test<-performance(fit_models[[i+1]],measure = "tpr", x.measure = "fpr")
  temp_BRT_test<-performance(fit_models[[i+2]],measure = "tpr", x.measure = "fpr")
  temp_MaxEnt_test<-performance(fit_models[[i+3]],measure = "tpr", x.measure = "fpr")
  # train
  temp_GAM_train<-performance(fit_models[[i+4]],measure = "tpr", x.measure = "fpr")
  temp_RF_train<-performance(fit_models[[i+5]],measure = "tpr", x.measure = "fpr")
  temp_BRT_train<-performance(fit_models[[i+6]],measure = "tpr", x.measure = "fpr")
  temp_MaxEnt_train<-performance(fit_models[[i+7]],measure = "tpr", x.measure = "fpr")
  
  # test
  temp_GAM_df_test<-data.frame(Run=i/5,Model="GAM",FPR=temp_GAM_test@x.values[[1]],TPR=temp_GAM_test@y.values[[1]]) %>% mutate(Data = "Test")
  temp_RF_df_test<-data.frame(Run=i/5,Model="RF",FPR=temp_RF_test@x.values[[1]],TPR=temp_RF_test@y.values[[1]])%>% mutate(Data = "Test")
  temp_BRT_df_test<-data.frame(Run=i/5,Model="BRT",FPR=temp_BRT_test@x.values[[1]],TPR=temp_BRT_test@y.values[[1]])%>% mutate(Data = "Test")
  temp_MaxEnt_df_test<-data.frame(Run=i/5,Model="MaxEnt",FPR=temp_MaxEnt_test@x.values[[1]],TPR=temp_MaxEnt_test@y.values[[1]])%>% mutate(Data = "Test")
  # train
  temp_GAM_df_train<-data.frame(Run=i/5,Model="GAM",FPR=temp_GAM_train@x.values[[1]],TPR=temp_GAM_train@y.values[[1]])%>% mutate(Data = "Train")
  temp_RF_df_train<-data.frame(Run=i/5,Model="RF",FPR=temp_RF_train@x.values[[1]],TPR=temp_RF_train@y.values[[1]])%>% mutate(Data = "Train")
  temp_BRT_df_train<-data.frame(Run=i/5,Model="BRT",FPR=temp_BRT_train@x.values[[1]],TPR=temp_BRT_train@y.values[[1]])%>% mutate(Data = "Train")
  temp_MaxEnt_df_train<-data.frame(Run=i/5,Model="MaxEnt",FPR=temp_MaxEnt_train@x.values[[1]],TPR=temp_MaxEnt_train@y.values[[1]])%>% mutate(Data = "Train")
  
  
  return(bind_rows(temp_GAM_df_test,temp_RF_df_test,temp_BRT_df_test,temp_MaxEnt_df_test,
                   temp_GAM_df_train,temp_RF_df_train,temp_BRT_df_train,temp_MaxEnt_df_train))
})%>%
  bind_rows()

ROCs_random$Run<-factor(ROCs_random$Run)
ROCs_random$Model<-factor(ROCs_random$Model,
                          levels=c("GAM","RF","BRT", "MaxEnt"),
                          labels = c("GAM","RF","BRT", "MaxEnt"))
ROCs_random$Data<-factor(ROCs_random$Data)
ROCs_random<-ROCs_random%>%mutate(Buffer=".25km2")
ROCs_random<-ROCs_random%>%mutate(Lambda="No")



#------Mean AUC ----------------

##-----Mean AUC township sampling -----
meanAUC_township<-data.frame(Model=c("GAM", "RF", "BRT", "MaxEnt", "GAM", "RF", "BRT", "MaxEnt"),
                             Data = c("Test","Test","Test","Test","Train","Train","Train","Train"),
                             meanAUC=c(mean(AUCs_township$GAM_test,na.rm = T),
                                       mean(AUCs_township$RF_test,na.rm = T),
                                       mean(AUCs_township$BRT_test,na.rm = T),
                                       mean(AUCs_township$MaxEnt_test,na.rm=T),
                                       mean(AUCs_township$GAM_train,na.rm=T),
                                       mean(AUCs_township$RF_train,na.rm=T),
                                       mean(AUCs_township$BRT_train,na.rm=T),
                                       mean(AUCs_township$MaxEnt_train,na.rm=T))
)
meanAUC_township <- meanAUC_township %>% mutate(Buffer = ".25km2")
meanAUC_township <- meanAUC_township %>% mutate(Lambda = "No")
meanAUC_township$Model<-factor(meanAUC_township$Model,levels=c("GAM", "RF", "BRT", "MaxEnt"))



##-----Mean AUC random sampling -----
meanAUC_random<-data.frame(Model=c("GAM", "RF", "BRT", "MaxEnt", "GAM", "RF", "BRT", "MaxEnt"),
                           Data = c("Test","Test","Test","Test","Train","Train","Train","Train"),
                           meanAUC=c(mean(AUCs_random$GAM_test,na.rm = T),
                                     mean(AUCs_random$RF_test,na.rm = T),
                                     mean(AUCs_random$BRT_test,na.rm = T),
                                     mean(AUCs_random$MaxEnt_test,na.rm=T),
                                     mean(AUCs_random$GAM_train,na.rm=T),
                                     mean(AUCs_random$RF_train,na.rm=T),
                                     mean(AUCs_random$BRT_train,na.rm=T),
                                     mean(AUCs_random$MaxEnt_train,na.rm=T))
)
meanAUC_random <- meanAUC_random %>% mutate(Buffer = ".25km2")
meanAUC_random <- meanAUC_random %>% mutate(Lambda = "No")
meanAUC_random$Model<-factor(meanAUC_random$Model,levels=c("GAM", "RF", "BRT", "MaxEnt"))




#----------Mechanistic Model---------
## ---------Mechanistic random sampling---------
set.seed(6)
runs = 50

Sys.time()
fit_models_lambda<-lapply(1:runs,function(i){
  # split albo into training and testing data (for GAM and RF) - by site (not datapoint)
  training_sites <- slice_sample(data.frame(unique(albomechanistic$Site)),prop=.8) 
  training_dataset_albo <- albomechanistic %>% filter(Site%in%training_sites$unique.albomechanistic.Site.)
  testing_dataset_albo <- albomechanistic %>% filter(!Site%in%training_sites$unique.albomechanistic.Site.)
  
  tryCatch({ ##This is done so that it doesn't stop everytime something throws an error
    print(paste0(i,": Fitting Models"))
    
    Albo_Lambda <- glm(Occurrence ~
                         Lambda,
                       data = training_dataset_albo,
                       family = binomial)
    
    #Predict outcome for the testing data
    print(paste0(i,": Predictions"))
    
    pred_Albo_Lambda_testing <- predict(Albo_Lambda, newdata = testing_dataset_albo,type="response") %>% data.frame
    pred_Albo_Lambda_training <- predict(Albo_Lambda, type="response",keep_prediction_data=T)%>% data.frame
    
    #Create dataframes of predictions and test data and save into respective lists
    #Names columns in dataframes with run number
    pred_df <- data.frame(Site=testing_dataset_albo$Site, Index= 1:dim(pred_Albo_Lambda_testing)[1],
                          Type="OOS",
                          pred_Albo_Lambda=pred_Albo_Lambda_testing,
                          Occurrence=testing_dataset_albo$Occurrence) 
    
    pred_df_training <- data.frame(Site=training_dataset_albo$Site, Index= 1:dim(pred_Albo_Lambda_training)[1],
                                   Type="WS",
                                   pred_Albo_Lambda=pred_Albo_Lambda_training,
                                   Occurrence=training_dataset_albo$Occurrence)
    
    pred_df <- rbind(pred_df,pred_df_training) %>% rename(pred_Albo_Lambda =".")
    
    # Scale suitability
    pred_df[,"pred_Albo_Lambda"] <- (pred_df$pred_Albo_Lambda - min(pred_df$pred_Albo_Lambda))/(max(pred_df$pred_Albo_Lambda)-min(pred_df$pred_Albo_Lambda))
    
    pred_df <- pred_df %>% 
      mutate(Run = i)
    
    #Create prediction objects for performance
    temp <- pred_df%>%
      filter(Type=="OOS")
    Albo_Lambda_pred_test<-prediction(temp$pred_Albo_Lambda,temp$Occurrence)
    
    temp <- pred_df%>%
      filter(Type=="WS")
    Albo_Lambda_pred_train<-prediction(temp$pred_Albo_Lambda,temp$Occurrence)
    
    
    #Return a list of needed objects
    r_list<-list(pred_df,
                 Albo_Lambda_pred_test,
                 Albo_Lambda_pred_train)
    # names(r_list_1)<-c("Pred_DF",
    #                      "GAM_Pred","RF_Pred","BRT_Pred","MaxEnt_Pred")
    names(r_list)<-c(paste0("Pred_DF: ",i),
                     paste0("Lambda_Pred_test: ",i),
                     paste0("Lambda_Pred_train: ",i))
    
    
    return(r_list)
  },
  error=function(cond){ #This tells it what to do if there is an error
    return(NULL)        #This just says to return NULL and go to the next iteration
  })
})%>%unlist(recursive = F)
Sys.time()
save(fit_models_lambda, file="/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/SDM_lambda.RData")

##---------AUC Random Sampling---------
AUCs_Lambda<-lapply(seq(2,(150),by=3),function(i){
  # Test
  temp_Lambda_test<-performance(fit_models_lambda[[i]],measure = "auc", x.measure = "cutoff")
  AUC_Lambda_test<-temp_Lambda_test@y.values[[1]]
  temp_Lambda_train<-performance(fit_models_lambda[[i+1]],measure = "auc", x.measure = "cutoff")
  AUC_Lambda_train<-temp_Lambda_train@y.values[[1]]
  
  temp<-data.frame(Lambda_test=AUC_Lambda_test,
                   Lambda_train=AUC_Lambda_train)
  return(temp)
})%>%
  bind_rows()

##------ROC Random Sampling----------------
ROCs_Lambda<-lapply(seq(2,(150),by=3),function(i){
  # test
  temp_Lambda_test<-performance(fit_models_lambda[[i]],measure = "tpr", x.measure = "fpr")
  # train
  temp_Lambda_train<-performance(fit_models_lambda[[i+1]],measure = "tpr", x.measure = "fpr")
  
  # test
  temp_Lambda_df_test<-data.frame(Run=i/5,Model="Lambda",FPR=temp_Lambda_test@x.values[[1]],TPR=temp_Lambda_test@y.values[[1]]) %>% mutate(Data = "Test")
  # train
  temp_Lambda_df_train<-data.frame(Run=i/5,Model="Lambda",FPR=temp_Lambda_train@x.values[[1]],TPR=temp_Lambda_train@y.values[[1]])%>% mutate(Data = "Train")
  
  return(bind_rows(temp_Lambda_df_test,
                   temp_Lambda_df_train))
})%>%
  bind_rows()

ROCs_Lambda$Run<-factor(ROCs_Lambda$Run)
ROCs_Lambda$Data<-factor(ROCs_Lambda$Data)

meanAUC_lambda <- data.frame(Model = c("Lambda"),
                             Data = c("Test", "Train"),
                             meanAUC=c(mean(AUCs_Lambda$Lambda_test,na.rm = T),
                                       mean(AUCs_Lambda$Lambda_train,na.rm = T)))



#----------Suitability thresholds---------
##--------Random Sampling-------
thresholds<-lapply(seq(2,(450),by=9),function(i){
  # test
  temp_GAM_test<-performance(fit_models[[i]],measure = "tpr", x.measure = "fpr")
  temp_RF_test<-performance(fit_models[[i+1]],measure = "tpr", x.measure = "fpr")
  temp_BRT_test<-performance(fit_models[[i+2]],measure = "tpr", x.measure = "fpr")
  temp_MaxEnt_test<-performance(fit_models[[i+3]],measure = "tpr", x.measure = "fpr")
  
  # test
  temp_GAM_df_test<-data.frame(Run=i/5,Model="GAM",FPR=temp_GAM_test@x.values[[1]],TPR=temp_GAM_test@y.values[[1]], Threshold=temp_GAM_test@alpha.values[[1]]) %>% mutate(Data = "Test")
  temp_GAM_df_test <- temp_GAM_df_test %>% filter(TPR>=0.8)
  result_GAM <- temp_GAM_df_test %>%
    arrange(FPR) %>%
    slice(1) 
  temp_RF_df_test<-data.frame(Run=i/5,Model="RF",FPR=temp_RF_test@x.values[[1]],TPR=temp_RF_test@y.values[[1]],Threshold=temp_RF_test@alpha.values[[1]])%>% mutate(Data = "Test")
  temp_RF_df_test <- temp_RF_df_test %>% filter(TPR>=0.8)
  result_RF <- temp_RF_df_test %>%
    arrange(FPR) %>%
    slice(1) 
  temp_BRT_df_test<-data.frame(Run=i/5,Model="BRT",FPR=temp_BRT_test@x.values[[1]],TPR=temp_BRT_test@y.values[[1]],Threshold=temp_BRT_test@alpha.values[[1]])%>% mutate(Data = "Test")
  temp_BRT_df_test <- temp_BRT_df_test %>% filter(TPR>=0.8)
  result_BRT <- temp_BRT_df_test %>%
    arrange(FPR) %>%
    slice(1) 
  temp_MaxEnt_df_test<-data.frame(Run=i/5,Model="MaxEnt",FPR=temp_MaxEnt_test@x.values[[1]],TPR=temp_MaxEnt_test@y.values[[1]],Threshold=temp_MaxEnt_test@alpha.values[[1]])%>% mutate(Data = "Test")
  temp_MaxEnt_df_test <- temp_MaxEnt_df_test %>% filter(TPR>=0.8)
  result_MaxEnt <- temp_MaxEnt_df_test %>%
    arrange(FPR) %>%
    slice(1) 
  
  return(bind_rows(result_GAM,result_RF,result_BRT,result_MaxEnt))
})%>%
  bind_rows()

thresholds <- thresholds %>% 
  group_by(Model) %>% 
  summarise(meanThreshold = mean(Threshold)) %>% 
  mutate(Buffer = ".25km2") %>% 
  mutate(Lambda = "No")

#------Thresholds for P/A------
BRTThreshold <- thresholds[[1,2]]
GAMThreshold <- thresholds[[2,2]]
MaxEntThreshold <- thresholds[[3,2]]
RFThreshold <- thresholds[[4,2]]

#make new Y/N column based on threshold 
LST_Predict <- LST_Predict %>%
  mutate(OccurrenceGAM = (SuitabilityGAM >= GAMThreshold),
         OccurrenceBRT = (SuitabilityBRT >= BRTThreshold),
         OccurrenceRF = (SuitabilityRF >= RFThreshold),
         OccurrenceMaxEnt = (SuitabilityMaxEnt >= MaxEntThreshold))
rm(BRTThreshold, GAMThreshold, MaxEntThreshold, RFThreshold)


#-----------Get Disagreement----------
LST_Predict <- LST_Predict %>% 
  mutate(OccurrenceGAM = as.integer(OccurrenceGAM),
         OccurrenceBRT = as.integer(OccurrenceBRT),
         OccurrenceRF = as.integer(OccurrenceRF),
         OccurrenceMaxEnt = as.integer(OccurrenceMaxEnt))

LST_Predict <- LST_Predict %>% 
  mutate(Concur = (abs(OccurrenceBRT - OccurrenceGAM)+
                     abs(OccurrenceBRT - OccurrenceMaxEnt)+
                     abs(OccurrenceBRT - OccurrenceRF)+
                     abs(OccurrenceGAM - OccurrenceMaxEnt)+
                     abs(OccurrenceGAM - OccurrenceRF)+
                     abs(OccurrenceMaxEnt - OccurrenceRF)))

#-----------Select Field Sites from July-September----------
# Filter for july-Septemeber
new_JAS <- LST_Predict %>% 
  filter(Month %in% c(7, 8, 9)) 

# make a dataset with total disagreement values (higher value = more disagreement)
JAS_disagreement <- new_JAS %>% 
  group_by(Site,coords.x1.x, coords.x2.x) %>% 
  summarize(TotalConcur = sum(Concur)) %>% 
  ungroup()

#We take the top disagreement sites (where disagreement is greater than 9)
site_candidates_JAS <- JAS_disagreement %>%
  filter(TotalConcur >= 9)
site_candidates_JAS$Site <- as.factor(site_candidates_JAS$Site)

# Now, I will sort through this site list for 30 sites that have high 
# disagreement but also are accessible for sampling by using google maps to view terrain
# and also discussing with collaborators in Suffolk County

# The final Sites chosen:
sites_list_AVS <- read_excel("Data/sites_list_AVS.xlsx")

# -------LOG-SCORING-----
# load in field collection data
field_mosquitoes2024_AVS <- read_excel("Data/field_mosquitoes2024_AVS.xlsx")
View(field_mosquitoes2024_AVS)

scoring_df <- new_JAS %>% filter(Site %in% sites_list_AVS$Site)
scoring_df <- scoring_df %>% dplyr::select(Month, Site, SuitabilityBRT, SuitabilityRF, SuitabilityGAM, SuitabilityMaxEnt)

# remove sites where traps were impacted and sampling did not occur (3 in total)
fieldmosq <- field_mosquitoes2024_AVS %>% filter(Trap_Status == "Y")

# Mark sites as positive for female albopictus based on total caught over sampling period
fieldmosq <- fieldmosq %>% dplyr::select(Site, Bag, Month, ALB_F)
fieldmosq <- fieldmosq %>% 
  group_by(Site, Month) %>% 
  summarize(Total_F = sum(ALB_F))
fieldmosq <- fieldmosq %>% 
  mutate(Occurrence_ALBF = (Total_F >= 1))
fieldmosq$Occurrence_ALBF <- fieldmosq$Occurrence_ALBF %>% 
  as.numeric()

# Use dplyr to mutate and map month abbreviations to numbers to match scoring dataframe
month_map <- c("July" = 7, "Aug" = 8, "Sept" = 9)
fieldmosq <- fieldmosq %>%
  mutate(Month = month_map[Month])

# combine field data and model pred
field_pred <- left_join(scoring_df, fieldmosq, by = c("Month", "Site"))

# remove any rows with NAs in the occurrence column 
field_pred <- field_pred %>% drop_na(Occurrence_ALBF)

# Calculate log scores 
GAM_LS <- logscore(Occurrence_ALBF~SuitabilityGAM, data = field_pred, group = "Month")
field_pred$GAM_LS <- GAM_LS$rawscores

BRT_LS <- logscore(Occurrence_ALBF~SuitabilityBRT, data = field_pred, group = "Month")
field_pred$BRT_LS <- BRT_LS$rawscores

RF_LS <- logscore(Occurrence_ALBF~SuitabilityRF, data = field_pred, group = "Month")
field_pred$RF_LS <- RF_LS$rawscores

MaxEnt_LS <- logscore(Occurrence_ALBF~SuitabilityMaxEnt, data = field_pred, group = "Month")
field_pred$MaxEnt_LS <- MaxEnt_LS$rawscores

# Compare log scores
# make data long
field_pred_long <- field_pred %>%
  pivot_longer(
    cols = c(GAM_LS, BRT_LS,RF_LS,MaxEnt_LS),
    names_to = "Model",
    values_to = "LogScore"
  )
field_pred_long <- field_pred_long %>% mutate(Model = recode(Model,
                                                             "GAM_LS"="GAM",
                                                             "MaxEnt_LS"="MaxEnt",
                                                             "RF_LS" = "Random Forest",
                                                             "BRT_LS" = "BRT"))

field_pred_long$Month <- field_pred_long$Month %>% as.factor()

Model1 <- aov(LogScore~Model*Month, data = field_pred_long)
summary(Model1)

# Tukey because there is a significant effect of the interaction of model and month
Tukey1 <- TukeyHSD(Model1, "Model:Month")
Tukey1

# visualize this mess 
# ggplot(field_pred_long, aes(x = Month, y = LogScore, fill = Model)) +
#   geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
#   labs(title = "Log-Transformed Scores by Month and Model",
#        x = "Month",
#        y = "Mean Log Score") +
#   theme_minimal()

# Extract p-values and determine significance groups
tukey1_table <- as.data.frame(Tukey1$`Model:Month`)
tukey1_table$pair <- rownames(tukey1_table)
tukey1_table <- tukey1_table %>% 
  mutate(significant = ifelse(`p adj` < 0.05, "*", ""))

data_summary <- field_pred_long %>%
  group_by(Model, Month) %>%
  summarize(mean_value = mean(LogScore), .groups = 'drop') %>%
  mutate(interaction_term = paste(Model, Month, sep = ":"))

p_values <- tukey1_table$`p adj`
names(p_values) <- tukey1_table$pair
letters <- multcompLetters(p_values)$Letters

data_summary <- data_summary %>%
  mutate(groups = letters[as.character(interaction_term)])


LogA <- ggplot(data_summary, aes(x = Month, y = -mean_value, color = Model, group = Model)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = groups, 
                y = -mean_value + ifelse(Model == "BRT", 0.03,
                                         ifelse(Model == "GAM", 0.02,
                                                ifelse(Model == "MaxEnt", -0.018, -0.03)))),size = 5)+
  labs(x = "Month",
       y = "Model Performance (negative log score)") +
  scale_x_discrete(labels = c("7" = "July", "8" = "August","9" = "September")) +
  scale_color_viridis_d(option = "viridis", labels = c("BRT", "GAM", "MaxEnt", "RF")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))


# how did models mispredict? were they biased towards over or underpredicting?


over_under <- new_JAS %>% filter(Site %in% sites_list_AVS$Site)
over_under <- left_join(over_under, fieldmosq, by = c("Month", "Site"))

over_under <- over_under %>% 
  mutate(BiasGAM = (OccurrenceGAM - Occurrence_ALBF),
         BiasRF=(OccurrenceRF - Occurrence_ALBF),
         BiasMaxEnt=(OccurrenceMaxEnt - Occurrence_ALBF),
         BiasBRT=(OccurrenceBRT - Occurrence_ALBF))

over_under_long <- over_under %>% 
  pivot_longer(cols = c(BiasGAM, BiasRF,BiasMaxEnt,BiasBRT),
               names_to = "Model",
               values_to = "Bias")

over_under_long$Month <- as.factor(over_under_long$Month)
levels(over_under_long$Month) <- c("July", "August", "September")

over_under_long$Model <- as.factor(over_under_long$Model)
levels(over_under_long$Model) <- c("BRT", "GAM", "MaxEnt", "RF")


LogB <- ggplot(over_under_long, aes(x = Month, y = Bias, fill = Model)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9), color = "black") +
  labs(x = "Month",
       y = "Model Bias") +
  scale_fill_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))

ggarrange(LogA, LogB, ncol = 2, labels = c("a)", "b)"))




# map of over and underpredicting 
bias_landuse_map <- LST_Predict %>% select(Year, Month, Site, coords.x1.x, coords.x2.x, Imp_.25km2, EVI_.25km2, DLST_1km2, NLST_1km2, Type_.25km2)
bias_landuse_map <- bias_landuse_map %>% filter(Year == "2023", Month == 7) 
#^^^^^land use does not change through a season but other variables do! Be cautious with mapping this!!!

spatial_bias_landuse_map <- st_as_sf(bias_landuse_map, coords = c("coords.x1.x","coords.x2.x"),crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))


boundary_AEA<-st_read("./Data/Remote Sensed Data/NYS_Civil_Boundaries/Counties_Shoreline.shp",quiet=T)%>%
  filter(NAME=="Suffolk")%>%
  st_transform(crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))%>%
  dplyr::select(geometry)%>%
  st_union()%>%
  as_Spatial()
boundary_AEA <- st_as_sf(boundary_AEA)

over_under_long_biasmap <- over_under_long %>% 
  group_by(Site, coords.x1.x, coords.x2.x, Model) %>% 
  summarize(Bias = mean(Bias)) %>% 
  ungroup()

over_under_long_biasmap <- over_under_long_biasmap %>% 
  mutate(bias_word = case_when(
    Bias > 0 ~ "overpredicting",
    Bias < 0 ~ "underpredicting",
    Bias == 0 ~ "accurate",
    TRUE ~ NA_character_  # Handles any other unexpected cases
  ))

ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_bias_landuse_map, aes(color = Type_.25km2), size = 1.5) +
  geom_point(data = over_under_long_biasmap, aes(x = coords.x1.x, y = coords.x2.x, fill = bias_word), color = "black", shape = 21, size = 3)+
  facet_wrap(Model~.)

# -----INVASION------

# Invasion looping (continuously add years of data then train the models)
if (file.exists("/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/invasion_series.RData")) {
  # If the file exists, load it
  load("/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/invasion_series.RData")
  print("File loaded successfully.")
} else {
  # If the file does not exist, run code to create the file
  print("File does not exist. Running code to create file.")
  
  fit_models_invasion<-lapply(1:16,function(i){
    # filter the full albo training dataset by year, adding a year each loop
    training_dataset_albo <- albo %>% filter(Year <= 2007+i)
    testing_dataset_albo <- LST_Predict %>% filter(Month %in% c(7, 8, 9))
    testing_dataset_albo <- testing_dataset_albo %>% filter(Site %in% sites_list_AVS$Site)
    
    # make datasets for BRT based on data partitioned above (can use same testing data as RF and GAM)
    training_dataset_albo_BRT <- sdmData(Occurrence ~ 
                                           Type +
                                           Year + 
                                           Month +
                                           Imp_.25km2 + 
                                           EVI_.25km2 + 
                                           DLST_1km2 +
                                           NLST_1km2 +
                                           Type_.25km2 +
                                           Longitude + 
                                           Latitude +
                                           Site+
                                           coords(Longitude+Latitude), 
                                         train=data.frame(training_dataset_albo))
    
    # Make datasets for MaxEnt Based on Data partitioned above
    # alboME make (with type_.25km2 made into numbers)
    training_dataset_albo_ME <- training_dataset_albo
    training_dataset_albo_ME <- training_dataset_albo_ME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                                         "Natural"=1,
                                                                                         "Altered"=2,
                                                                                         "Developed - Low"=3,
                                                                                         "Developed - Med/High"=4))
    training_dataset_albo_ME_m <- sdmData(Occurrence ~
                                            Type +
                                            Year +
                                            Month +
                                            Imp_.25km2 +
                                            EVI_.25km2 +
                                            DLST_1km2 +
                                            NLST_1km2 +
                                            Type_.25km2 +
                                            coords(Longitude+Latitude),
                                          train=data.frame(training_dataset_albo_ME))
    
    training_dataset_albo_ME <- training_dataset_albo_ME %>%
      dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
    
    testing_dataset_albo_ME <- testing_dataset_albo
    testing_dataset_albo_ME <- testing_dataset_albo_ME %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                                       "Natural"=1,
                                                                                       "Altered"=2,
                                                                                       "Developed - Low"=3,
                                                                                       "Developed - Med/High"=4))
    testing_dataset_albo_ME <- testing_dataset_albo_ME %>% 
      dplyr::select(Type,Year, Month,Imp_.25km2, EVI_.25km2,DLST_1km2,NLST_1km2, Type_.25km2)
    
    
    tryCatch({ ##This is done so that it doesn't stop everytime something throws an error
      print(paste0(i,": Fitting Models"))
      # GAM and Random Forest
      Albo_GAM_.25 <- gam(Occurrence ~
                            Type +
                            Year +
                            Month +
                            s(Imp_.25km2,k=5) + 
                            s(EVI_.25km2,k=5) +
                            s(DLST_1km2,k=5) +
                            s(NLST_1km2,k=5) +
                            Type_.25km2 +
                            s(Site,bs = "re"),
                          data = training_dataset_albo,
                          family = binomial)
      
      Albo_RF_.25 <- randomForest(Occurrence ~ 
                                    Type +
                                    Year + 
                                    Month +
                                    Imp_.25km2 + 
                                    EVI_.25km2 + 
                                    DLST_1km2 +
                                    NLST_1km2 +
                                    Type_.25km2, 
                                  data = training_dataset_albo)
      
      # Boosted Regression Tree (BRT)
      Albo_BRT_.25 <- sdm(Occurrence ~ 
                            Type +
                            Year + 
                            Month +
                            Imp_.25km2 + 
                            EVI_.25km2 + 
                            DLST_1km2 +
                            NLST_1km2 +
                            Type_.25km2, 
                          data = training_dataset_albo_BRT,
                          methods = c('brt'))
      
      # MaxEnt
      Albo_MaxEnt_.25 <- sdm(Occurrence ~ 
                               Type +
                               Year + 
                               Month +
                               Imp_.25km2 + 
                               EVI_.25km2 + 
                               DLST_1km2 +
                               NLST_1km2+
                               Type_.25km2, 
                             data = training_dataset_albo_ME_m,
                             methods = c('maxent'))
      
      #Predict outcome for the testing data
      print(paste0(i,": Predictions"))
      
      pred_Albo_GAM_.25_testing <- predict(Albo_GAM_.25, newdata = testing_dataset_albo,exclude="s(Site)",type="response")
      # pred_Albo_GAM_.25_training <- predict(Albo_GAM_.25, type="response",keep_prediction_data=T)
      
      pred_Albo_RF_.25_testing <- predict(Albo_RF_.25, newdata = testing_dataset_albo)
      # pred_Albo_RF_.25_training <- predict(Albo_RF_.25)
      
      pred_Albo_BRT_.25_testing <- predict(Albo_BRT_.25, newdata = testing_dataset_albo, method = "brt", mean=T)%>% data.frame() %>% rename(pred_Albo_BRT_.25_testing=id_1.sp_1.m_brt)
      # pred_Albo_BRT_.25_training <- predict(Albo_BRT_.25, newdata = training_dataset_albo, method = "brt", mean=T)%>% data.frame() %>% rename(pred_Albo_BRT_.25_training=id_1.sp_1.m_brt)
      
      pred_Albo_MaxEnt_.25_testing <- predict(Albo_MaxEnt_.25, newdata=testing_dataset_albo_ME, method='maxent')%>% data.frame() %>% rename(pred_Albo_MaxEnt_.25_testing=id_1.sp_1.m_maxent)
      # pred_Albo_MaxEnt_.25_training <- predict(Albo_MaxEnt_.25, newdata=training_dataset_albo_ME, method='maxent')%>% data.frame() %>% rename(pred_Albo_MaxEnt_.25_training=id_1.sp_1.m_maxent)
      
      #Create dataframes of predictions and test data and save into respective lists
      #Names columns in dataframes with run number
      pred_df_.25 <- data.frame(Site=testing_dataset_albo$Site,
                                Type="OOS",
                                Month = testing_dataset_albo$Month,
                                SuitabilityGAM=pred_Albo_GAM_.25_testing,
                                SuitabilityRF=pred_Albo_RF_.25_testing,
                                SuitabilityBRT=pred_Albo_BRT_.25_testing,
                                SuitabilityMaxEnt=pred_Albo_MaxEnt_.25_testing) %>% 
        rename(SuitabilityBRT=pred_Albo_BRT_.25_testing, SuitabilityMaxEnt=pred_Albo_MaxEnt_.25_testing)
      
      # pred_df_.25_training <- data.frame(Site=Albo_GAM_.25$model$Site, Index= 1:dim(pred_Albo_GAM_.25_training)[1],
      #                                    Type="WS",
      #                                    pred_Albo_GAM_.25=pred_Albo_GAM_.25_training,
      #                                    pred_Albo_RF_.25=pred_Albo_RF_.25_training,
      #                                    pred_Albo_BRT_.25=pred_Albo_BRT_.25_training,
      #                                    pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_training,
      #                                    Occurrence=Albo_GAM_.25$model$Occurrence)%>% rename(pred_Albo_BRT_.25=pred_Albo_BRT_.25_training, pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_training)
      # 
      # pred_df_.25 <- rbind(pred_df_.25,pred_df_.25_training)
      
      
      # Scale suitability
      pred_df_.25[,"SuitabilityGAM"] <- (pred_df_.25$SuitabilityGAM - min(pred_df_.25$SuitabilityGAM))/(max(pred_df_.25$SuitabilityGAM)-min(pred_df_.25$SuitabilityGAM))
      pred_df_.25[,"SuitabilityRF"] <- (pred_df_.25$SuitabilityRF - min(pred_df_.25$SuitabilityRF))/(max(pred_df_.25$SuitabilityRF)-min(pred_df_.25$SuitabilityRF))
      pred_df_.25[,"SuitabilityBRT"] <- (pred_df_.25$SuitabilityBRT - min(pred_df_.25$SuitabilityBRT))/(max(pred_df_.25$SuitabilityBRT)-min(pred_df_.25$SuitabilityBRT))
      pred_df_.25[,"SuitabilityMaxEnt"] <- (pred_df_.25$SuitabilityMaxEnt - min(pred_df_.25$SuitabilityMaxEnt))/(max(pred_df_.25$SuitabilityMaxEnt)-min(pred_df_.25$SuitabilityMaxEnt))
      
      
      pred_df_.25 <- pred_df_.25 %>% 
        mutate(TrainingYrs = 2007+i)
      
      
      
      # log scoring 
      field_pred_temp <- left_join(pred_df_.25, fieldmosq, by = c("Month", "Site"))
      field_pred_temp <- field_pred_temp %>% drop_na(Occurrence_ALBF)
      GAM_LS_temp <- logscore(Occurrence_ALBF~SuitabilityGAM, data = field_pred_temp, group = "Month")
      field_pred_temp$GAM_LS <- GAM_LS_temp$rawscores
      
      BRT_LS_temp <- logscore(Occurrence_ALBF~SuitabilityBRT, data = field_pred_temp, group = "Month")
      field_pred_temp$BRT_LS <- BRT_LS_temp$rawscores
      
      RF_LS_temp <- logscore(Occurrence_ALBF~SuitabilityRF, data = field_pred_temp, group = "Month")
      field_pred_temp$RF_LS <- RF_LS_temp$rawscores
      
      MaxEnt_LS_temp <- logscore(Occurrence_ALBF~SuitabilityMaxEnt, data = field_pred_temp, group = "Month")
      field_pred_temp$MaxEnt_LS <- MaxEnt_LS_temp$rawscores
      
      
      
      
      # #Create prediction objects for performance
      # temp <- pred_df_.25%>%
      #   filter(Type=="OOS")
      # 
      # Albo_GAM_.25_pred_test<-prediction(temp$pred_Albo_GAM_.25,temp$Occurrence)
      # Albo_RF_.25_pred_test <- prediction(temp$pred_Albo_RF_.25,temp$Occurrence)
      # Albo_BRT_.25_pred_test <- prediction(temp$pred_Albo_BRT_.25,temp$Occurrence)
      # Albo_MaxEnt_.25_pred_test <- prediction(temp$pred_Albo_MaxEnt_.25,temp$Occurrence)
      # 
      # temp <- pred_df_.25%>%
      #   filter(Type=="WS")
      # Albo_GAM_.25_pred_train<-prediction(temp$pred_Albo_GAM_.25,temp$Occurrence)
      # Albo_RF_.25_pred_train <- prediction(temp$pred_Albo_RF_.25,temp$Occurrence)
      # Albo_BRT_.25_pred_train <- prediction(temp$pred_Albo_BRT_.25,temp$Occurrence)
      # Albo_MaxEnt_.25_pred_train <- prediction(temp$pred_Albo_MaxEnt_.25,temp$Occurrence)
      
      #Return a list of needed objects
      r_list_.25<-list(field_pred_temp)
      
      
      # names(r_list_.25)<-c("Pred_DF",
      #                      "GAM_Pred","RF_Pred","BRT_Pred","MaxEnt_Pred")
      names(r_list_.25)<-c(paste0("Pred_DF: ",2007+i))
      
      
      return(r_list_.25)
    },
    error=function(cond){ #This tells it what to do if there is an error
      return(NULL)        #This just says to return NULL and go to the next iteration
    })
  })%>%unlist(recursive = F)
  
  save(fit_models_invasion, file = "/Users/avs79/Desktop/R/SDM_LI_ALB_FV/Data/SDM/invasion_series.RData")
  
  print("File created and saved successfully.")
}




invasion_scores <- do.call(rbind, fit_models_invasion)

invasion_scores_long <- invasion_scores %>% pivot_longer(cols = c(GAM_LS, BRT_LS,RF_LS,MaxEnt_LS),
                                                         names_to = "Model",
                                                         values_to = "LogScore")

invasion_scores_long <- invasion_scores_long %>% mutate(Model = recode(Model,
                                                                       "GAM_LS"="GAM",
                                                                       "MaxEnt_LS"="MaxEnt",
                                                                       "RF_LS" = "Random Forest",
                                                                       "BRT_LS" = "BRT"))

invasion_scores_long$Month <- invasion_scores_long$Month %>% as.factor()

# truncade scores of inf to be 10, and scores of 0 to be 0.00001
invasion_scores_long$LogScore <- ifelse(invasion_scores_long$LogScore >= 10, 10, invasion_scores_long$LogScore)
invasion_scores_long$LogScore <- ifelse(invasion_scores_long$LogScore == 0, 0.00001, invasion_scores_long$LogScore)

# plot
invasion_summary <- invasion_scores_long %>%
  group_by(Model, TrainingYrs) %>%
  summarise(
    mean = mean(LogScore),
    se = sd(LogScore)/sqrt(n()),
    .groups = 'drop')


# plot
ggplot()+
  geom_errorbar(data = invasion_summary, aes(ymin = -(mean - se), ymax = -(mean + se), x = TrainingYrs, group = Model, color = Model), width = 0.2, stat = "identity")+
  geom_line(data = invasion_summary, aes(x = TrainingYrs, y = -mean, group = Model, color = Model))+
  theme_minimal()+
  scale_color_viridis_d(option = "viridis", labels = c("BRT","GAM", "MaxEnt", "RF"))+
  labs(x = "Training data", y = "Model Performance (negative log score)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))+
  scale_x_continuous(breaks = invasion_summary$TrainingYrs, labels = invasion_summary$TrainingYrs)
  


#---------Plots---------

##------field Sites-----
sites <- read_excel("Data/final_SITES.xlsx")
sites <- sites %>% mutate(legend = "Field Site")

boundary_AEA<-st_read("./Data/Remote Sensed Data/NYS_Civil_Boundaries/Counties_Shoreline.shp",quiet=T)%>%
  filter(NAME=="Suffolk")%>%
  st_transform(crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))%>%
  dplyr::select(geometry)%>%
  st_union()%>%
  as_Spatial()
boundary_AEA <- st_as_sf(boundary_AEA)

ggplot()+
  geom_sf(data=boundary_AEA,aes(geometry=geometry), fill="grey90",color="black")+
  geom_point(data=sites, 
             aes(x = LongNEW, y = LatNEW, color = legend), size = 1.5) + 
  labs(x = "Longitude", y = "Latitude", title = "2024 Field Sites")+
  theme_minimal() + 
  scale_color_manual(values = "black", name = NULL)+
  theme(legend.position = c(.9, .1))

##--------Visualizing invasion---------

first_ALB <- alboraw %>%
  filter(Mosquitoes > 0) %>% 
  group_by(Site) %>%                  # Group the data by site
  arrange(Site, Date) %>%             # Arrange by site and date
  slice_min(n = 1, order_by = Date)

first_ALB$Year <- as.factor(first_ALB$Year)

ggplot()+
  geom_sf(data=boundary_AEA,aes(geometry=geometry), fill="grey90",color="black")+
  geom_point(data=first_ALB, 
             aes(x = Longitude.y, y = Latitude.y, color = Year), size = 3) + 
  scale_color_viridis_d(option = "viridis")+
  labs(x = "Longitude", y = "Latitude")+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.7, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )


##-----maps for  suitability across models in prediction dataset----
# make data long 
LST_Predict_long <- LST_Predict %>%
  pivot_longer(
    cols = c(SuitabilityGAM, SuitabilityBRT,SuitabilityRF,SuitabilityMaxEnt),
    names_to = "suitability_type",
    values_to = "suitability_value") %>% ungroup()

LST_Predict_long$Month <- as.factor(LST_Predict_long$Month)
LST_Predict_long$suitability_type <- as.factor(LST_Predict_long$suitability_type)

levels(LST_Predict_long$Month) <- c("May", "June", "July", "August", "September", "October")
levels(LST_Predict_long$suitability_type) <- c("BRT", "GAM", "MaxEnt", "RF")

LST_Predict_long <- LST_Predict_long %>% 
  group_by(Site, Month, suitability_type, coords.x1.x, coords.x2.x) %>% 
  summarise(suitability = mean(suitability_value)) %>% ungroup()

spatial_LST_Predict_long <- st_as_sf(LST_Predict_long, coords = c("coords.x1.x","coords.x2.x"),crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# A
ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_long %>% filter(suitability_type == "GAM"), aes(color = suitability), size = .75) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Habitat Suitability",        # Title for the legend
                        breaks = c(0.1, .45, 0.9),
                        labels = c("Low", "Medium", "High")) +  # Suitability color scale
  theme_minimal() +
  facet_wrap(~Month)+
  labs(title = "Suitability Map A", x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right")
# B
ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_long %>% filter(suitability_type == "RF"), aes(color = suitability), size = .75) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Habitat Suitability",        # Title for the legend
                        breaks = c(0.1, .45, 0.9),
                        labels = c("Low", "Medium", "High")) +  # Suitability color scale
  theme_minimal() +
  facet_wrap(~Month)+
  labs(title = "Suitability Map B", x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right")
# C
ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_long %>% filter(suitability_type == "BRT"), aes(color = suitability), size = .75) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Habitat Suitability",        # Title for the legend
                        breaks = c(0.1, .45, 0.9),
                        labels = c("Low", "Medium", "High")) +  # Suitability color scale
  theme_minimal() +
  facet_wrap(~Month)+
  labs(title = "Suitability Map C", x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right")
# D
ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_long %>% filter(suitability_type == "MaxEnt"), aes(color = suitability), size = .75) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Habitat Suitability",        # Title for the legend
                        breaks = c(0.1, .45, 0.9),
                        labels = c("Low", "Medium", "High")) +  # Suitability color scale
  theme_minimal() +
  facet_wrap(~Month)+
  labs(title = "Suitability Map D", x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right")

# average across months 
ave_LST_Predict_long <- LST_Predict_long %>% 
  group_by(Site, suitability_type, coords.x1.x, coords.x2.x) %>% 
  summarise(suitability = mean(suitability)) %>% ungroup()

spatial_ave_LST_Predict_long <- st_as_sf(ave_LST_Predict_long, coords = c("coords.x1.x","coords.x2.x"),crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_ave_LST_Predict_long, aes(color = suitability), size = 1) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Habitat Suitability",        # Title for the legend
                        breaks = c(0.1, .45, 0.8),
                        labels = c("Low", "Medium", "High")) +  # Suitability color scale
  theme_minimal() +
  facet_wrap(~suitability_type)+
  labs(x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 15))+
  annotation_custom(
          grob = textGrob("Source: gis.ny.gov", x = 0.5, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)


# graphical abstract figure
ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_ave_LST_Predict_long, aes(color = suitability), size = 1.5) +
  scale_color_viridis_c(
    option = "plasma",
    direction = -1,
    name = "Habitat Suitability",
    breaks = c(0.1, 0.45, 0.8),
    labels = c("Low", "Medium", "High"),
    guide = guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5)
  ) +
  theme_minimal() +
  facet_wrap(~suitability_type, ncol = 1) +
  labs(title = "Suitability Map", x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


##-----predictor maps----
boundary_AEA<-st_read("./Data/Remote Sensed Data/NYS_Civil_Boundaries/Counties_Shoreline.shp",quiet=T)%>%
  filter(NAME=="Suffolk")%>%
  st_transform(crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))%>%
  dplyr::select(geometry)%>%
  st_union()%>%
  as_Spatial()
boundary_AEA <- st_as_sf(boundary_AEA)

LST_Predict_predictors <- LST_Predict %>% 
  group_by(Site, coords.x1.x, coords.x2.x) %>% 
  summarise(ave_imp = mean(Imp_.25km2),
            ave_DLST = mean(DLST_1km2),
            ave_NLST = mean(NLST_1km2),
            ave_EVI = (mean(EVI_.25km2)*0.0001),
            ave_LC = first(Type_.25km2)) %>% 
  ungroup()

spatial_LST_Predict_predictors <- st_as_sf(LST_Predict_predictors, coords = c("coords.x1.x","coords.x2.x"),crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# imp
IMP <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_imp), size = 1) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "% Impervious Surface") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
# DLST
DLST <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_DLST), size = 1) +
  scale_color_viridis_c(option = "inferno", direction = 1, name = "Day Surface Temperature") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
# NLST
NLST <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_NLST), size = 1) +
  scale_color_viridis_c(option = "inferno", direction = 1, name = "Night Surface Temperature") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.6, y = 0.1, hjust = 0, gp = gpar(fontsize = 8)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
# EVI
EVI <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_EVI), size = 1) +
  scale_color_viridis_c(option = "mako", direction = 1, name = "EVI                                       ") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
# Land cover
LC <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_LC), size = 1) +
  # scale_color_viridis(option = "plasma", name = "Landcover") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Land Cover") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

ggarrange(IMP, LC, EVI, DLST, NLST, ncol = 1, nrow = 5, labels = c("a)","b)","c)","d)","e)"))


##-----suitability estimates across months and models----
model2 <- aov(suitability~suitability_type*Month, data = LST_Predict_long)
summary(model2)

tukey2 <- TukeyHSD(model2, "suitability_type:Month")
tukey2

# Extract p-values and determine significance groups
tukey2_table <- as.data.frame(tukey2$`suitability_type:Month`)
tukey2_table$pair <- rownames(tukey2_table)
tukey2_table <- tukey2_table %>% 
  mutate(significant = ifelse(`p adj` < 0.05, "*", ""))

data2_summary <- LST_Predict_long %>%
  group_by(suitability_type, Month) %>%
  summarize(mean_value = mean(suitability), .groups = 'drop') %>%
  mutate(interaction_term = paste(suitability_type, Month, sep = ":"))

p_values2 <- tukey2_table$`p adj`
names(p_values2) <- tukey2_table$pair
letters2 <- multcompLetters(p_values2)$Letters

data2_summary <- data2_summary %>%
  mutate(groups = letters2[as.character(interaction_term)])

ggplot(data2_summary, aes(x = Month, y = mean_value, fill = suitability_type))+
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9), color = "black")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
    )+
  scale_fill_viridis_d(option = "viridis")+
  geom_text(aes(label = groups, group = suitability_type), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 5)+
  labs( x = "Month", y = "Average Habitat Suitability", fill = "Model") 


##---------AUC---------

# random sampling
ggplot(ROCs_random,aes(x=FPR, y=TPR, color=Run, linetype = Data))+
  geom_line(alpha=.5,linewidth=.25)+
  geom_smooth(color="black")+
  geom_abline(slope=1,intercept = 0,linetype="dotted")+
  geom_label(data=meanAUC_random %>% filter(Data == "Test"),aes(x=.7,y=0.05,label=paste("Mean Testing AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  geom_label(data=meanAUC_random %>% filter(Data == "Train"),aes(x=.7,y=0.15,label=paste("Mean Training AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  scale_color_discrete(guide=NULL)+
  facet_wrap(Model~.)+
  xlab("False Positivity Rate")+
  ylab("True Positivity Rate")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 15))

# township sampling
ggplot(ROCs_township,aes(x=FPR, y=TPR, color=Run, linetype = Data))+
  geom_line(alpha=.5,linewidth=.25)+
  geom_smooth(color="black")+
  geom_abline(slope=1,intercept = 0,linetype="dotted")+
  geom_label(data=meanAUC_township %>% filter(Data == "Test"),aes(x=.7,y=0.05,label=paste("Mean Testing AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  geom_label(data=meanAUC_township %>% filter(Data == "Train"),aes(x=.7,y=0.15,label=paste("Mean Training AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  scale_color_discrete(guide=NULL)+
  facet_wrap(Model~.)+
  xlab("False Positivity Rate")+
  ylab("True Positivity Rate")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 15))

# lambda with random sampling 
ggplot(ROCs_Lambda,aes(x=FPR, y=TPR, color=Run, linetype = Data))+
  geom_line(alpha=.5,linewidth=.25)+
  geom_smooth(color="black")+
  geom_abline(slope=1,intercept = 0,linetype="dotted")+
  geom_label(data=meanAUC_lambda %>% filter(Data == "Test"),aes(x=.7,y=0.05,label=paste("Mean Testing AUC = ",round(meanAUC,digits = 3))),
             color="black",size=6)+
  geom_label(data=meanAUC_lambda %>% filter(Data == "Train"),aes(x=.7,y=0.15,label=paste("Mean Training AUC = ",round(meanAUC,digits = 3))),
             color="black",size=6)+
  scale_color_discrete(guide=NULL)+
  xlab("False Positivity Rate")+
  ylab("True Positivity Rate")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 15))

##------maps of disagreement across months or totalled with one map with FIELD SITE LOCATIONS------
spatial_JAS_disagreement <- st_as_sf(JAS_disagreement, coords = c("coords.x1.x","coords.x2.x"),crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

spatial_new_JAS <- st_as_sf(new_JAS, coords = c("coords.x1.x","coords.x2.x"),crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_JAS_disagreement, aes(color = TotalConcur), size = 2.25) +
  scale_color_distiller(limits = c(-1,13), palette = "YlOrRd", direction=+1,
                       breaks = c(0, 6, 12),
                       labels = c("Low", "Medium", "High")) +
  labs(color = "Disagreement", fill = NULL, x = "Longitude", y = "Latitude") +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))+
  geom_point(data = sites, 
             aes(x = LongNEW, y = LatNEW, fill = legend), size = 2)+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.7, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

# what land uses did models disagree the most in across july-sept?

landuse_disagreement <- new_JAS %>% group_by(Site) %>% summarise(Disagreement = sum(Concur),
                                                                 Type_.25km2 = first(Type_.25km2))
landuse_disagreement <- landuse_disagreement %>% mutate(Type_.25km2 = recode(Type_.25km2,
                                                                     "Natural" = "Natural",
                                                                     "Altered" = "Altered",
                                                                     "Developed - Med/High" = "DevMed/High",
                                                                     "Developed - Low" = "DevLow"))


landuse_bias <- landuse_disagreement %>% group_by(Type_.25km2) %>% 
  summarise(Count = n(), 
            Percentage = (n()/922)*100)

model3 <- aov(Disagreement ~ Type_.25km2, data = landuse_disagreement)
summary(model3)
tukey3 <- TukeyHSD(model3)


# Extract p-values and determine significance groups
tukey3_table <- as.data.frame(tukey3$Type_.25km2)
tukey3_table$pair <- rownames(tukey3_table)
tukey3_table <- tukey3_table %>% 
  mutate(significant = ifelse(`p adj` < 0.05, "*", ""))

data3_summary <- landuse_disagreement %>%
  group_by(Type_.25km2) %>%
  summarize(mean_value = mean(Disagreement), .groups = 'drop')

p_values3 <- tukey3_table$`p adj`
names(p_values3) <- tukey3_table$pair
letters3 <- multcompLetters(p_values3)$Letters

data3_summary <- data3_summary %>%
  mutate(groups = letters3[as.character(Type_.25km2)])

ggplot(data3_summary, aes(x = Type_.25km2, y = mean_value))+
  geom_bar(stat = "summary", fun = "mean", fill = "#472F7D", color = "black")+
  theme_minimal()+
  geom_text(aes(label = groups), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 5)+
  labs(x = "Landuse", y = "Disagreement")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))+
  scale_x_discrete(labels = c("Altered", "Developed - Low", "Developed - Medium/High", "Natural"))







## ------varimp------
variables <- c("Imp_.25km2", "EVI_.25km2", "DLST_1km2", "NLST_1km2", "Month", "Year", "Type_.25km2")

strings <- c("Altered", "Natural", "Developed - Med/High", "Developed - Low")

# Use lapply to iterate over variable names, not indices
pdp <- lapply(variables, function(var_name){
  # Calculate feature grid for each variable based on its actual range in 'albo'
  feature_grid <- seq(min(albo[[var_name]], na.rm = TRUE), max(albo[[var_name]], na.rm = TRUE), length.out = 100)
  feature_grid_type <- sample(strings, size = 100, replace = TRUE)
  
  
  
  # prediction data for all models except maxent
  prediction_data <- data.frame(
    Imp_.25km2 = if(var_name != "Imp_.25km2") mean(albo$Imp_.25km2, na.rm = TRUE) else feature_grid,
    EVI_.25km2 = if(var_name != "EVI_.25km2") mean(albo$EVI_.25km2, na.rm = TRUE) else feature_grid,
    DLST_1km2 = if(var_name != "DLST_1km2") mean(albo$DLST_1km2, na.rm = TRUE) else feature_grid,
    NLST_1km2 = if(var_name != "NLST_1km2") mean(albo$NLST_1km2, na.rm = TRUE) else feature_grid,
    Type_.25km2  = if(var_name != "Type_.25km2") "Developed - Low" else feature_grid_type,
    Type = 1,
    Year = if(var_name != "Year") mean(albo$Year, na.rm = TRUE) else feature_grid,
    Month = if(var_name != "Month") mean(albo$Month, na.rm = TRUE) else feature_grid,
    Site = 1)
  prediction_data$Type <- as.factor(prediction_data$Type)
  prediction_data$Type <- factor(prediction_data$Type, levels = c(1, 0, "BG Sentinel", "Mosquito Magnet"), 
                                 labels = c("1","0", "BG Sentinel", "Mosquito Magnet"))
  prediction_data$Type_.25km2<- as.factor(prediction_data$Type_.25km2)
  prediction_data$Type_.25km2 <- factor(prediction_data$Type_.25km2, 
                                        levels = c("Altered","Developed - Low","Developed - Med/High","Natural"),
                                        labels = c("Altered","Developed - Low","Developed - Med/High","Natural"),
                                        ordered = TRUE)
  
  numbers <- c(2,1,4,3)
  
  feature_grid_type_ME <- numbers[match(feature_grid_type, strings)]
  
  # prediction data for maxent
  prediction_data_ME <- data.frame(
    Imp_.25km2 = if(var_name != "Imp_.25km2") mean(alboME$Imp_.25km2, na.rm = TRUE) else feature_grid,
    EVI_.25km2 = if(var_name != "EVI_.25km2") mean(alboME$EVI_.25km2, na.rm = TRUE) else feature_grid,
    DLST_1km2 = if(var_name != "DLST_1km2") mean(alboME$DLST_1km2, na.rm = TRUE) else feature_grid,
    NLST_1km2 = if(var_name != "NLST_1km2") mean(alboME$NLST_1km2, na.rm = TRUE) else feature_grid,
    Type_.25km2  = if(var_name != "Type_.25km2") 3 else feature_grid_type_ME,
    Type = 1,
    Year = if(var_name != "Year") mean(alboME$Year, na.rm = TRUE) else feature_grid,
    Month = if(var_name != "Month") mean(alboME$Month, na.rm = TRUE) else feature_grid)
  prediction_data_ME$Type <- as.factor(prediction_data_ME$Type)
  prediction_data_ME$Type <- factor(prediction_data_ME$Type, levels = c(1, 0, "BG Sentinel", "Mosquito Magnet"), 
                                    labels = c("1","0", "BG Sentinel", "Mosquito Magnet"))
  
  # Predict using the model
  predictions_GAM <- predict(AlboGAM, newdata = prediction_data,exclude="s(Site)")%>% data.frame()
  predictions_RF <- predict(AlboRF, newdata = prediction_data)%>% data.frame()
  predictions_BRT <- predict(AlboBRT, newdata = prediction_data, method = "brt", mean=T)%>% data.frame()
  predictions_MaxEnt <- predict(AlboMaxEnt, newdata=prediction_data_ME, method='maxent')%>% data.frame() 
  
  names(predictions_GAM)[names(predictions_GAM) == "."] <- "GAM"
  names(predictions_RF)[names(predictions_RF) == "."] <- "RF"
  names(predictions_BRT)[names(predictions_BRT) == "id_1.sp_1.m_brt"] <- "BRT"
  names(predictions_MaxEnt)[names(predictions_MaxEnt) == "id_1.sp_1.m_maxent"] <- "MaxEnt"
  
  dataframe <- cbind(feature_grid, feature_grid_type, predictions_GAM, predictions_RF, predictions_BRT, predictions_MaxEnt)
  names(dataframe)[names(dataframe) == "feature_grid"] <- var_name
  
  # # scaling 
  # dataframe[,"GAM"] <- (dataframe$GAM - min(dataframe$GAM))/(max(dataframe$GAM)-min(dataframe$GAM))
  # dataframe[,"RF"] <- (dataframe$RF - min(dataframe$RF))/(max(dataframe$RF)-min(dataframe$RF))
  # dataframe[,"BRT"] <- (dataframe$BRT - min(dataframe$BRT))/(max(dataframe$BRT)-min(dataframe$BRT))
  # dataframe[,"MaxEnt"] <- (dataframe$MaxEnt - min(dataframe$MaxEnt))/(max(dataframe$MaxEnt)-min(dataframe$MaxEnt))
  
  dataframe <- dataframe %>% pivot_longer(
    cols = c(GAM, BRT,RF,MaxEnt),
    names_to = "Model",
    values_to = "Suitability")
  
  
  # pred_df_.25 <- data.frame(Variable= feature_grid, 
  #                           Type= var_name,
  #                           GAM=predictions_GAM,
  #                           # MaxEnt=predictions_MaxEnt,
  #                           BRT=predictions_BRT,
  #                           MaxEnt=predictions_RF) 
  # # %>% rename(pred_Albo_BRT_.25=pred_Albo_BRT_.25_testing, pred_Albo_MaxEnt_.25=pred_Albo_MaxEnt_.25_testing)
  # 
  # 
  # rlist <- list(dataframe)
  # names(rlist)<-c(paste0(var_name))
  # 
  return(dataframe)
  
})


# Plot imp
pdp_Imp <- pdp[[1]]
pdp_EVI <- pdp[[2]]
pdp_DLST <- pdp[[3]]
pdp_NLST <- pdp[[4]]
pdp_Month <- pdp[[5]]
pdp_Year <- pdp[[6]]
pdp_Type <- pdp[[7]]


pdpIMP <- ggplot(pdp_Imp, aes(x = Imp_.25km2, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "Impervious Surface", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
pdpEVI <- ggplot(pdp_EVI, aes(x = EVI_.25km2*0.0001, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "EVI", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
pdpDLST <- ggplot(pdp_DLST, aes(x = DLST_1km2, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "Day Land Surface Temperature", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
pdpNLST <- ggplot(pdp_NLST, aes(x = NLST_1km2, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "Night Land Surface Temperature", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
pdpMonth <- ggplot(pdp_Month, aes(x = Month, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "Month", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
pdpYear <- ggplot(pdp_Year, aes(x = Year, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "Year", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
pdpType <- ggplot(pdp_Type, aes(x = feature_grid_type, y = Suitability, color = Model))+
  geom_point()+
  labs(x = "Land Cover Type", y = "Suitability")+
  scale_color_viridis_d(option = "viridis")+
  theme_minimal()+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 25))


ggarrange(pdpIMP, pdpEVI, pdpDLST, pdpNLST, pdpMonth, pdpYear, pdpType, ncol = 3, nrow =3, labels = c("a)","b)","c)","d)","e)","f)","g)"), common.legend = TRUE)


## --------site sampling------

# plot all unique sites colored by the frequncy of years used
sampling_ALB <- alboraw %>%
  group_by(Site, Township, Latitude.y, Longitude.y) %>%                  # Group the data by site
  summarize(n_years = n_distinct(Year), , .groups = "drop")

site_summary <- sampling_ALB %>%
  count(n_years, name = "n_sites") %>%
  mutate(percent = 100 * n_sites / sum(n_sites))

site_summary

# plotting
boundary_AEA<-st_read("./Data/Remote Sensed Data/NYS_Civil_Boundaries/Counties_Shoreline.shp",quiet=T)%>%
  filter(NAME=="Suffolk")%>%
  st_transform(crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))%>%
  dplyr::select(geometry)%>%
  st_union()%>%
  as_Spatial()
boundary_AEA <- st_as_sf(boundary_AEA)

ggplot()+
  geom_sf(data=boundary_AEA,aes(geometry=geometry), fill="grey90",color="black")+
  geom_point(data=sampling_ALB,
             aes(x = Longitude.y, y = Latitude.y, color = n_years), size = 1) +
  scale_color_distiller(limits = c(0,17), palette = "YlOrRd", direction=+1,
                        breaks = c(1, 4, 8, 12, 16),
                        labels = c("1", "4", "8", "12", "16")) +
  labs(x = "Longitude", y = "Latitude", color = "Sampling Years")+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.7, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
