
################################################################################
## Title: Community characteristics of COVID-19 vaccine hesitancy in England: ## 
##                  A nationwide cross-sectional study                        ##
##                                                                            ##
## Authors: Bucyibaruta, G. - Blangiardo, M. -  Konstantinoudis, G.           ##
##                                                                            ##
##                                                                            ##
##                                                                            ##
################################################################################

# Loading packages

library(readODS);library(sf);library(leaflet);library(RColorBrewer)
library(tidyverse)
library(xtable)
library(classInt)
library(viridis)
library(gridExtra)
library(Hmisc)
library(ggplot2)
library(INLA)
library(spdep)
library(ggplot2)
library(patchwork)
library(ghibli)
library(ggplotify)
library(ggtext)
#--------------------------------------------------------------

###--------- Get MSOA boundaries -----------------###
#------- Download the MSOA Shapefile (polygon layer)-------#

url.shp <- "https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/England_msoa_2011_sgen_clipped.zip"
download.file(url.shp, destfile = "temp.zip", method = "auto" , quiet=FALSE)
unzip("temp.zip")
# Delete the original zip file
unlink("temp.zip")

# Read vac_msoa
msoa_poly <- st_read("england_msoa_2011_sgen_clipped.shp")
# Fix projection
msoa_poly <- st_transform(msoa_poly, 4326)
sf::sf_use_s2(FALSE)


###-------------- Loading  the vaccine  data------------###

data_Dec <- read.csv("msoa_2022-01-02.csv",header = T)


names(data_Dec)[11] <- "PopVac"
names(data_Dec)[12] <- "Cum_firstDose"
names(data_Dec)[13] <- "First_dosePct"

##---Linking vaccine data with shapefile
vac_poly_Dec <- merge(msoa_poly, data_Dec[,c(7,11,12,13)],  by.x = "code", by.y = "areaCode", all.x = T)

data_Vaccine <- read.csv("msoa_2022-01-02.csv",header = T)


names(data_Vaccine)[6] <- "PopVac"
names(data_Vaccine)[7] <- "Cum_firstDose"
names(data_Vaccine)[8] <- "First_dosePct"

vac_poly_Dec <- merge(msoa_poly, data_Vaccine[,c(3,6,7,8)],  by.x = "code", by.y = "areaCode", all.x = T)

##################-- Loading Covariate in orginal scale--------#################

datav_og <- read.csv("DataVac_OG.csv", header = T)

datav_og <- datav_og[,-c(1,3,4,5,6,12,16,24,29,30)]

###-----Calculating min, quintiles and max 

round(quantile(datav_og$bme,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$pctYoungAge,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$PctOld,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$Pct_LAB,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$Pct_con,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$Pct_Leave,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$DeathRates,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$Asthma,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$blood_pre,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$Diabetes,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)
round(quantile(datav_og$Depression,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)



datav_og$rateCat <- factor(datav_og$rateCat)
datav_og$rateCat <- relevel(datav_og$rateCat,ref = "low")
datav_og$imd_decile_msoaLast <- factor(datav_og$imd_decile_msoaLast)
datav_og$B_RUC11 <- factor(datav_og$B_RUC11)
datav_og$B_RUC11 <- relevel(datav_og$B_RUC11,ref = "PU")
datav_og <- fastDummies::dummy_cols(datav_og,
                                    select_columns = "B_RUC11",
                                    remove_first_dummy = TRUE)
vac_poly_Dec_og <- merge(vac_poly_Dec, datav_og,  by.x = "code", by.y = "areaCode", all.x = T)


######-----Convert covariates into categories-------######

vac_poly_Dec_og$bmeG <- cut2(vac_poly_Dec_og$bme,g=5)
levels(vac_poly_Dec_og$bmeG) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$leave <- cut2(vac_poly_Dec_og$Pct_Leave,g=5)
levels(vac_poly_Dec_og$leave) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$asthma <- cut2(vac_poly_Dec_og$Asthma,g=5)
levels(vac_poly_Dec_og$asthma) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$diabetes <- cut2(vac_poly_Dec_og$Diabetes,g=5)
levels(vac_poly_Dec_og$diabetes) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$bp <- cut2(vac_poly_Dec_og$blood_pre,g=5)
levels(vac_poly_Dec_og$bp) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$depression <- cut2(vac_poly_Dec_og$Depression,g=5)
levels(vac_poly_Dec_og$depression) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$Cons <- cut2(vac_poly_Dec_og$Pct_con,g=5)
levels(vac_poly_Dec_og$Cons) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$LAB <- cut2(vac_poly_Dec_og$Pct_LAB,g=5)
levels(vac_poly_Dec_og$LAB) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$death <- cut2(vac_poly_Dec_og$rates,g=5)
levels(vac_poly_Dec_og$death) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$young <- cut2(vac_poly_Dec_og$pctYoungAge,g=5)
levels(vac_poly_Dec_og$young) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$old <- cut2(vac_poly_Dec_og$PctOld,g=5)
levels(vac_poly_Dec_og$old) <- c("1(low)", "2","3","4","5")

vac_poly_Dec_og$imd_rank_msoaLast <-round(rank(-vac_poly_Dec_og$imd_medianS_msoa))
vac_poly_Dec_og$imd_quintile <- cut2(vac_poly_Dec_og$imd_rank_msoaLast, g=5)
levels(vac_poly_Dec_og$imd_quintile) <- 
  c("low", "second","third","fourth","high")
vac_poly_Dec_og$imd_quintile <- relevel(vac_poly_Dec_og$imd_quintile,ref = "high")


vac_poly_Dec_og$urbanicity <- vac_poly_Dec_og$B_RUC11

levels(vac_poly_Category$urbanicity)<  c("Urban","Rural","Semi")
#vac_poly_Category <- vac_poly_Dec_og

##--------- Subseting  categorical data to use in model------------##

vac_poly_Category <-vac_poly_Dec_og[,c(1,2,3,4,5,46:60)]

##---Reverse the order of levels in IMD
vac_poly_Category$imd_quintileNew <- ifelse(vac_poly_Category$imd_quintile=="low", "5",
                                            ifelse(vac_poly_Category$imd_quintile=="second", "4",
                                                   ifelse(vac_poly_Category$imd_quintile=="third","3",
                                                          ifelse(vac_poly_Category$imd_quintile=="fourth","2",
                                                                 "1(least)"))))
vac_poly_Category$imd_quintileNew <- factor(vac_poly_Category$imd_quintileNew, levels = c("1(least)", "2","3","4","5"))


####------ Laoding vaccine access data-------------######

Vac_Acces <- read.csv("Weighted_VAc_Acc_Data.csv")

round(quantile(Vac_Acces$WeightVacAcc,probs = c(0,0.2, 0.4, 0.6,0.8,1)),2)

vac_poly_Category <- merge(vac_poly_Category, Vac_Acces[,c(2,7,8)],  by.x = "code", by.y = "code", all.x = T)

vac_poly_Category$WVA <- cut2(vac_poly_Category$WeightVacAcc, g=5)
levels(vac_poly_Category$WVA) <- c("1(low)", "2","3","4","5")



############################################################################################
###############################----Data exploration----------#############################
###########################################################################################

vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = bmeG), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("BME")->Bme
#--------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = imd_quintileNew), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title = element_text(face="bold"),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  labs(fill='IMD')+
  ggtitle("IMD")->imd
#-------------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = young), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Young")->Young

#---------------------

vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = old), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Old")->Old

#---------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = death), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Mortality")->Death

#------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = urbanicity), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Urbanicity")->Urban

#----------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = LAB), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Labour")->lab

#----------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = Cons), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Conservatives")->cons

#------------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = leave), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Brexit")->Leave

#----------------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = asthma), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Asthma")->Asthma

#------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = bp), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Blood Pressure")->BP

#--------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = diabetes), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Diabetes")->Diab

#----------------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = depression), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Depression")->Depres

#-----------
vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = WVA), color = NA, show.legend = TRUE) +
  theme_void() +
  #scale_fill_viridis_d(name = "", alpha = 0.8,begin = 0,
  #end = 1,
  #direction = 1,
  #option = "D",
  #aesthetics = "colour")  +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.5,end=1)+
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Access")->access

Socio=((Bme|Old|Young)/(imd|Urban))

x11()

Socio +plot_layout(guides = "collect")+
  plot_annotation(title = "Socio-Demographics and urbanicity",
                  theme = theme(plot.title = element_text(face = "bold",size = 14)))

ggsave("Socio-Demographics.pdf",width=10,height =8)
#-----------------------------

poli_acc=((lab|cons)/(Leave|access))

x11()
poli_acc +plot_layout(guides = "collect")+
  plot_annotation(title = "Politics view and Access",
                  theme = theme(plot.title = element_text(face = "bold",size = 14)))

ggsave("Politics and Access.pdf",width=10,height =8)

#-----------------------------

Health=((Death|Asthma|BP)/(Diab|Depres))

x11()
Health +plot_layout(guides = "collect")+
  plot_annotation(title = "COVID-19 awareness and targeting of high risk groups",
                  theme = theme(plot.title = element_text(face = "bold",size = 14)))

ggsave("health-conditions.pdf",width=10,height =8)


##########################################################################
###-----------Correlation between covariates---------------##############
#########################################################################


DD <- data.frame(vac_poly_Category$bmeG,
                 vac_poly_Category$imd_quintileNew,
                 vac_poly_Category$young,
                 vac_poly_Category$old,
                 vac_poly_Category$death,
                 vac_poly_Category$urbanicity,
                 vac_poly_Category$LAB,
                 vac_poly_Category$Cons,
                 vac_poly_Category$leave,
                 vac_poly_Category$asthma,
                 vac_poly_Category$bp,
                 vac_poly_Category$diabetes,
                 vac_poly_Category$depression,
                 vac_poly_Category$WVA)

colnames(DD) <- c("BME", "IMD", "Young",
                  "Old", "Mortality", "Urbanicity",
                  "Labour", "Conservatives", "Brexit",
                  "Asthma", "Bood pressure","Diabetes",
                  "Depression", "Access")

levels(DD$BME) <- c(1,2,3,4,5)
levels(DD$IMD) <- c(1,2,3,4,5)
levels(DD$Young) <- c(1,2,3,4,5)
levels(DD$ Old) <- c(1,2,3,4,5)
levels(DD$Mortality) <- c(1,2,3,4,5)
levels(DD$Urbanicity) <- c(1,2,3)
levels(DD$Labour) <- c(1,2,3,4,5)
levels(DD$Conservatives) <- c(1,2,3,4,5)
levels(DD$Brexit) <- c(1,2,3,4,5)
levels(DD$Asthma) <- c(1,2,3,4,5)
levels(DD$`Bood pressure`) <- c(1,2,3,4,5)
levels(DD$Diabetes) <- c(1,2,3,4,5)
levels(DD$Depression) <- c(1,2,3,4,5)
levels(DD$Access) <- c(1,2,3,4,5)

library("varhandle")
DD=unfactor(DD)

cov_cor=round(cor(DD, method="kendall", use="pairwise"),2)
cov_cor

#----Using ggcorrplot (plot to save)

library(ggcorrplot)

# plotting corr heatmap
x11()
ggcorrplot(cov_cor,lab = TRUE,colors = c("red", "white", "blue"))

ggsave("Heatmap_cor3.pdf",width = 10,height = 10)

#########################################################################
###-------------Fitting model in INLA--------------------------#########
########################################################################

vac_poly_Category$reid <- 1:nrow(vac_poly_Category)

nb <- poly2nb(vac_poly_Category)

nb2INLA("map_adj",nb)

# Graph
g <- inla.read.graph(filename="map_adj")

# With spatial (all categorical)
formulaNewFirstG <- Cum_firstDose~ 1+bmeG+imd_quintileNew+young+old+
  urbanicity +LAB+Cons+leave+death+
  asthma+bp+diabetes+depression+WVA+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstG = inla(formulaNewFirstG, family = "binomial", 
                          data = vac_poly_Category, Ntrials=PopVac,
                          control.predictor = list(compute=TRUE),
                          control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

summary(vac_ful3New2FirstG)


###--------Variability explained by covariance----------###
Data_Dummy <- fastDummies::dummy_cols(vac_poly_Category[,c(6:16,19,20,24)],
                                      select_columns = c("bmeG","leave","asthma","diabetes","bp","depression",
                                                         "Cons","LAB","death","young","old","imd_quintileNew","urbanicity","WVA"),
                                      remove_first_dummy = TRUE)
data_og <- Data_Dummy[,16:69]

set.seed(123)
samplesSpatial = inla.posterior.sample(10000,vac_ful3New2FirstG)



funtotalCovSpatial = function(...) {
  c(Intercept+bmeGsecond*data_og$bmeG_second+bmeGthird*data_og$bmeG_third+bmeGfourth*data_og$bmeG_fourth+
      imd_quintileNewsecond*data_og$imd_quintileNew_second+imd_quintileNewthird*data_og$imd_quintileNew_third+
      imd_quintileNewfourth*data_og$imd_quintileNew_fourth+imd_quintileNewhigh*data_og$imd_quintileNew_high+
      leavesecond*data_og$leave_second+leavethird*data_og$leave_third+leavefourth*data_og$leave_fourth+leavehigh*data_og$leave_high+
      asthmasecond*data_og$asthma_second+asthmathird*data_og$asthma_third+asthmafourth*data_og$asthma_fourth+asthmahigh*data_og$asthma_high+
      diabetessecond*data_og$diabetes_second+diabetesthird*data_og$diabetes_third+diabetesfourth*data_og$diabetes_fourth+diabeteshigh*data_og$diabetes_high+
      bpsecond*data_og$bp_second+bpthird*data_og$bp_third+bpfourth*data_og$bp_fourth+bphigh*data_og$bp_high+
      depressionsecond*data_og$depression_second+depressionthird*data_og$depression_third+depressionfourth*data_og$depression_fourth+depressionhigh*data_og$depression_high+
      Conssecond*data_og$Cons_second+Consthird*data_og$Cons_third+Consfourth*data_og$Cons_fourth+Conshigh*data_og$Cons_high+
      LABsecond*data_og$LAB_second+LABthird*data_og$LAB_third+LABfourth*data_og$LAB_fourth+LABhigh*data_og$LAB_high+
      deathsecond*data_og$death_second+deaththird*data_og$death_third+deathfourth*data_og$death_fourth+deathhigh*data_og$death_high+
      youngsecond*data_og$young_second+youngthird*data_og$young_third+youngfourth*data_og$young_fourth+younghigh*data_og$young_high+
      urbanicityPR*data_og$urbanicity_PR+urbanicityUR*data_og$urbanicity_UR+
      oldsecond*data_og$old_second+oldthird*data_og$old_third+oldfourth*data_og$old_fourth+oldhigh*data_og$old_high+
      WVAsecond*data_og$WVA_second+WVAthird*data_og$WVA_third+WVAfourth*data_og$WVA_fourth+
      WVAhigh*data_og$WVA_high+
      reid[1:6791])
}

funtotalCov = function(...) {
  c(Intercept+
      bmeGsecond*data_og$bmeG_second+bmeGthird*data_og$bmeG_third+bmeGfourth*data_og$bmeG_fourth+
      imd_quintileNewsecond*data_og$imd_quintileNew_second+imd_quintileNewthird*data_og$imd_quintileNew_third+
      imd_quintileNewfourth*data_og$imd_quintileNew_fourth+imd_quintileNewhigh*data_og$imd_quintileNew_high+
      leavesecond*data_og$leave_second+leavethird*data_og$leave_third+leavefourth*data_og$leave_fourth+leavehigh*data_og$leave_high+
      asthmasecond*data_og$asthma_second+asthmathird*data_og$asthma_third+asthmafourth*data_og$asthma_fourth+asthmahigh*data_og$asthma_high+
      diabetessecond*data_og$diabetes_second+diabetesthird*data_og$diabetes_third+diabetesfourth*data_og$diabetes_fourth+diabeteshigh*data_og$diabetes_high+
      bpsecond*data_og$bp_second+bpthird*data_og$bp_third+bpfourth*data_og$bp_fourth+bphigh*data_og$bp_high+
      depressionsecond*data_og$depression_second+depressionthird*data_og$depression_third+depressionfourth*data_og$depression_fourth+depressionhigh*data_og$depression_high+
      Conssecond*data_og$Cons_second+Consthird*data_og$Cons_third+Consfourth*data_og$Cons_fourth+Conshigh*data_og$Cons_high+
      LABsecond*data_og$LAB_second+LABthird*data_og$LAB_third+LABfourth*data_og$LAB_fourth+LABhigh*data_og$LAB_high+
      deathsecond*data_og$death_second+deaththird*data_og$death_third+deathfourth*data_og$death_fourth+deathhigh*data_og$death_high+
      youngsecond*data_og$young_second+youngthird*data_og$young_third+youngfourth*data_og$young_fourth+younghigh*data_og$young_high+
      urbanicityPR*data_og$urbanicity_PR+urbanicityUR*data_og$urbanicity_UR+
      oldsecond*data_og$old_second+oldthird*data_og$old_third+oldfourth*data_og$old_fourth+oldhigh*data_og$old_high
    +WVAsecond*data_og$WVA_second+WVAthird*data_og$WVA_third+WVAfourth*data_og$WVA_fourth+
      WVAhigh*data_og$WVA_high)
}

fcov = inla.posterior.sample.eval(funtotalCov,samplesSpatial)
fcovSpatial = inla.posterior.sample.eval(funtotalCovSpatial, samplesSpatial)


fftotalc <- as.matrix(fcov)
ff_covcs <- as.matrix(fcovSpatial)


# Variance due to covariates
New_sample=apply(fftotalc, 2, var)/apply(ff_covcs, 2, var)

M=mean(New_sample)
SD=sd(New_sample)
M-1.96*SD
M+1.96*SD

##-------------------------------------------------------------------###

##-----------Extract Odd Ratios ------------------------##
#names of levels of covariates were changed in numbers
tab_spatial <- exp(vac_ful3New2FirstG$summary.fixed[
  c("bmeGsecond","bmeGthird","bmeGfourth","bmeGhigh",
    "imd_quintileNewsecond","imd_quintileNewthird","imd_quintileNewfourth","imd_quintileNewhigh",
    "youngsecond","youngthird","youngfourth","younghigh",
    "oldsecond","oldthird","oldfourth","oldhigh",
    "deathsecond","deaththird","deathfourth","deathhigh",
    "urbanicityPR","urbanicityUR",
    "LABsecond","LABthird","LABfourth","LABhigh",
    "Conssecond","Consthird","Consfourth","Conshigh",
    "leavesecond","leavethird","leavefourth","leavehigh",
    "asthmasecond","asthmathird","asthmafourth","asthmahigh",
    "bpsecond","bpthird","bpfourth","bphigh",
    "diabetessecond","diabetesthird","diabetesfourth","diabeteshigh",
    "depressionsecond","depressionthird","depressionfourth","depressionhigh",
    "WVAsecond","WVAthird","WVAfourth", "WVAhigh"
  ),
  c("0.5quant", "0.025quant", "0.975quant")])
tab_spatial <- round(tab_spatial, digits = 3)

##----------------------------------------------------------------------

########################################################################
###-------Customizing Caterpillar plots for OR-----------------#########
########################################################################


New= c("1 (low)","2","3","4","5",
       "1 (least)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "Urban","Rural","Semi",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5"
)
A <- c(1,6,11,16,21,24,29,34,39,44,49,54,59,64)
B <- c(2:5,7:10,12:15,17:20,22:23,25:28,30:33,35:38,40:43,45:48,50:53,55:58,60:63,65:68)
N<- rep(1,14)
M= exp(vac_ful3New2FirstG$summary.fixed$`0.5quant`)[-1]
Lower = exp(vac_ful3New2FirstG$summary.fixed$`0.025quant`)[-1]
Upper= exp(vac_ful3New2FirstG$summary.fixed$`0.975quant`)[-1]
MN <- c()
MN[A]=N
MN[B]=M
LN <- c()
LN[A]=N
LN[B]=Lower
UN <- c()
UN[A]=N
UN[B]=Upper
df <- data.frame(index=c(1:23, 25:39, 41:65, 67:71),variable=New,Mean=MN,l=LN,u=UN)
length(c(1:28, 30:45, 47:67, 69:76))
df$class <- ifelse(df$index <24, "Socio_demographic and urbanicity",
                   ifelse(df$index >= 24 & df$index <= 40, "Political_belief",
                          ifelse(df$index > 40 & df$index <= 65,
                                 "COVID-19 awareness and targeting of high risk groups",
                                 "Vaccine_Accessibility")))


df$class <- factor(df$class,levels = c("Socio_demographic and urbanicity","Political_belief",
                                       "COVID-19 awareness and targeting of high risk groups",
                                       "Vaccine_Accessibility" ))

breaks = c(0.5,0.6,0.7,0.8,0.9,1,1.1, 1.2,1.3,1.4, 1.5)
labels = as.character(breaks)


pm <- ggplot(data=df)+
  geom_vline(xintercept=1, linetype="dashed") +
  geom_errorbar(aes(y=index, x=Mean, xmin=l, xmax=u))+
  #scale_color_manual("class",values=c( "#CC79A7","#009E73","#0072B2","#D55E00"))+
  #this adds the effect sizes to the plot
  geom_point(aes(x=Mean,y=index))+
  #adds the CIs
  scale_y_continuous(name = "", breaks=c(1:23, 25:39, 41:65, 67:71),
                     labels = df$variable, trans="reverse",
                     expand = c(0, 0)) +
  #adding a vertical line at the effect = 0 mark
  #thematic stuff
  ggtitle("")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.y = element_markdown(size = rel(1.2)),
    axis.title.y = element_text(size = rel(1.8) ),
    panel.spacing.y  = unit(2, "lines"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, 1, 0, 1), "cm")
  ) +
  xlab("Odds ratio")  +
  geom_hline(yintercept = c(5.5, 10.5, 15.5, 20.5, 24,
                            29.5, 34.5, 40, 45.5,
                            50.5, 55.5, 60.5, 66), col = "grey85")+
  geom_hline(yintercept = c(24, 40, 66), col = "black") +
  geom_vline(xintercept = breaks,
             col = "grey85",linetype="dashed") +
  scale_x_continuous(labels = labels,breaks = breaks)+
  annotate("text",
           x = rep(0.339, lenght.out = 14),angle = 0,
           #y = c(69, 64, 57, 54, 48, 44, 40, 35, 30, 24, 19, 14, 9, 3),
           y=c(3,8,13,18,22,27,32,37,43,48,53,58,63,69),
           label = c("BME", "IMD", "Young",
                     "Old", "Urbanicity",
                     "Labour", "Conservatives", "Brexit", "Mortality",
                     "Asthma", "Bood pressure","Diabetes",
                     "Depression", "Access")) +
  annotate("text",
           x = rep(1.60, lenght.out =3 ),angle = 3*90,
           y = c(10, 33, 53, 69),
           label = c("Socio-demographics\n and urbanicity",
                     "Political opinions",
                     "COVID-19 awareness and \n targeting of high risk groups",
                     "Access"), fontface = "bold", size =4
           #,colour=c( "#CC79A7","#009E73","#0072B2","#D55E00")
  ) +
  geom_vline(xintercept=1, linetype="dashed",col="red") +
  coord_cartesian(xlim = c(0.5, 1.52), clip = "off")
x11()
pm

ggsave("CovariatesEstimatesNew5.pdf", width = 10, height = 14)

#########################################################################
###----Mapping posterior median and spatial residuals (relative risk)--##
########################################################################

pred <-vac_ful3New2FirstG$summary.fitted.values$`0.5quant`
predrnk <- rank(pred)
low <- vac_ful3New2FirstG$summary.fitted.values$`0.025quant`
high <- vac_ful3New2FirstG$summary.fitted.values$`0.975quant`

vac_poly_Category$predictor<- pred*100

vac_poly_Category$pred_low<- low*100
vac_poly_Category$pred_upper<- high*100
vac_poly_Category$groupPred <- cut2(vac_poly_Category$predictor,g=5)

Q_prob <- quantile(vac_poly_Category$predictor,probs = c(0.2,0.4,0.6,0.8))

######-----vaccine Uptake and large cities--------------------##########
cities <- read_sf("Major_Towns_and_Cities_(December_2015)_Boundaries_V2.shp")
cities <- st_transform(cities, 4326)

London <- cities[cities$TCITY15NM %in% "London",]           
Birmingham <- cities[cities$TCITY15NM %in% "Birmingham",]   
Liverpool <- cities[cities$TCITY15NM %in% "Liverpool",]    
Bristol <- cities[cities$TCITY15NM %in% "Bristol",]   

################################################################
######------------Relative ratio---------------------##########
###############################################################

Bris <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Bristol",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = group), col = NA, show.legend = F) +
  geom_sf(data = Bristol, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  #scale_fill_viridis_d(option = "inferno",direction = -1) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end=1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("Bristol")

L <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "London",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = group), col = NA, show.legend = F) +
  geom_sf(data = London, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                            "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end=1) +
  #scale_fill_viridis_d(option = "inferno",direction = -1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("London")

Bir <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Birmingham",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = group), col = NA, show.legend = F) +
  geom_sf(data = Birmingham, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end=1) +
  #scale_fill_viridis_d(option = "inferno",direction = -1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("Birmingham")

Liv <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Liverpool",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = group), col = NA, show.legend = F) +
  geom_sf(data = Liverpool, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  scale_fill_viridis_d(option="viridis",direction = -1,begin = 0.2,end = 1) +
  #scale_fill_viridis_d(option = "inferno",direction = -1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("Liverpool")

vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = group), color = NA, show.legend = TRUE) +
  geom_sf(data = London, colour = "red", fill = NA, size = .1) + 
  geom_sf(data = Birmingham, colour = "red", fill = NA, size = .1) + 
  geom_sf(data = Liverpool, colour = "red", fill = NA, size = .1) + 
  geom_sf(data = Bristol, colour = "red", fill = NA, size = .1) + 
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                            "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end=1)+
  #scale_fill_viridis_d(option = "inferno",direction = -1) +
  theme_void() +
  # theme_bw() (to check!)
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        legend.text = element_text(size=14),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0) # Left margin
        ,plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("England") -> a

p_all=a|((L|Bir)/(Liv|Bris))
p_all +plot_layout(guides = "collect")+
  plot_annotation(title = "")

ggsave("VaccUptake_CitiesViridisN.pdf",width=10,height =8)

################################################################
######------------Relative ratio---------------------##########
###############################################################
ff = function(x){exp(x)}
vac_poly_Category$bym2RR <- round(ff(vac_ful3New2FirstG$summary.random$reid$'0.5quant'[1:6791]),2)
vac_poly_Category$bym2GR <-cut2(vac_poly_Category$bym2RR,g=6)


L_rr <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "London",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), col = NA, show.legend = F) +
  geom_sf(data = London, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" ))+
  #scale_fill_viridis_d(option ="G",begin=1, end=0) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  #scale_fill_viridis_d(option = "inferno",direction = -1,begin = 0.5,end = 1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("London")
L_rr

Bir_rr <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Birmingham",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), col = NA, show.legend = F) +
  geom_sf(data = Birmingham, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" ))+
  #scale_fill_viridis_d(option ="G",begin=1, end=0) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  #scale_fill_viridis_d(option = "inferno",direction = -1,begin = 0.5,end = 1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("Birmingham")
Bir_rr


Liv_rr <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Liverpool",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), col = NA, show.legend = F) +
  geom_sf(data = Liverpool, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" ))+
  #scale_fill_viridis_d(option ="G",begin=1, end=0) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  #scale_fill_viridis_d(option = "inferno",direction = -1,begin = 0.5,end=1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("Liverpool")
Liv_rr

Man_rr <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Manchester",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), col = NA, show.legend = F) +
  geom_sf(data = Manchester, fill = NA,col=NA) +
  scale_fill_manual(values = c("#FC4E07","#0072B2" , 
                               "#00AFBB","#F0E442", "#CC79A7","#999999" ))+
  #scale_fill_viridis_d(option ="G",begin=1, end=0) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Manchester")
Man_rr

Lds_rr <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Leeds",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), col = NA, show.legend = F) +
  geom_sf(data = Leeds, fill = NA,col=NA) +
  scale_fill_manual(values = c("#FC4E07","#0072B2" , 
                               "#00AFBB","#F0E442", "#CC79A7","#999999" ))+
  #scale_fill_viridis_d(option ="G",begin=1, end=0) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Leeds")

Bris_rr <- st_intersection(vac_poly_Category, cities[cities$TCITY15NM %in% "Bristol",]) %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), col = NA, show.legend = F) +
  geom_sf(data = Bristol, fill = NA,col=NA) +
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                             "#00AFBB","#F0E442", "#CC79A7","#999999" ))+
  #scale_fill_viridis_d(option ="G",begin=1, end=0) +
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  #scale_fill_viridis_d(option = "inferno",direction = -1,begin = 0.5,end = 1) +
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("Bristol")
Bris_rr

vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = bym2GR), color = NA, show.legend = TRUE) +
  geom_sf(data = London, colour = "red", fill = NA, size = .1) + 
  geom_sf(data = Birmingham, colour = "red", fill = NA, size = .1) + 
  geom_sf(data = Liverpool, colour = "red", fill = NA, size = .1) + 
  geom_sf(data = Bristol, colour = "red", fill = NA, size = .1) + 
  #scale_fill_manual(values = c("#FC4E07","#0072B2" , 
  #                            "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  #scale_fill_viridis_d(option ="G",begin=1, end=0)+
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  #scale_fill_viridis_d(option = "inferno",direction = -1,begin = 0.5,end = 1) +
  theme_void() +
  # theme_bw() (to check!)
  theme(legend.position = c(0.92, 0.85),
        legend.title =  element_blank(),
        legend.text = element_text(size=14),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0) # Left margin
        ,plot.title = element_text(hjust = 0.5,face = "bold"))+
  ggtitle("England") -> a_rr
a_rr


p_all_rr=a_rr|((L_rr|Bir_rr)/(Liv_rr|Bris_rr))

x11()

p_all_rr +plot_layout(guides = "collect")+
  plot_annotation(title = "")

ggsave("RelativeRatio_CitiesVirN.pdf",width=10,height = 8)

##########################################################################
###########----Heatmap of profiles and Poesterior probabilities-----#####
#########################################################################

data1<- vac_poly_Category
data1$pred <- pred
data1$predrnk <- predrnk
data1$low <- low


data1$bmeG <- as.numeric(data1$bmeG)
data1$imd_quintile <- as.numeric(data1$imd_quintile)
data1$deathNew <- as.numeric(data1$deathNew)
data1$young <- as.numeric(data1$young)
data1$old <- as.numeric(data1$old)
data1$depression <- as.numeric(data1$depression)
data1$urbanicity <- as.numeric(data1$urbanicity)
data1$Cons <- as.numeric(data1$Cons)
data1$LAB <- as.numeric(data1$LAB)
data1$leave <- as.numeric(data1$leave)
data1$asthma <- as.numeric(data1$asthma)
data1$diabetes <- as.numeric(data1$diabetes)
data1$bp <- as.numeric(data1$bp)
data1$WVacAcc <- as.numeric(data1$WVacAcc)

data1$groupPred <- cut2(data1$pred,g=5)


data1$geometry <- NULL
library(tibble)
library(dplyr)
library(ggplot2)

data1<-as_tibble(data1)
data1 <- data1[,-11]

data1$imd_quintileNew <- ifelse(data1$imd_quintile==1, 5,
                                ifelse(data1$imd_quintile==2, 4,
                                       ifelse(data1$imd_quintile==3,3,
                                              ifelse(data1$imd_quintile==2,4,
                                                     1))))

Data2<-data1 %>%
  group_by(bmeG,imd_quintileNew,deathNew,young,old,urbanicity,
           LAB,Cons,leave,asthma, bp,depression, diabetes,WVacAcc) %>%
  summarise(mean=mean(pred),sd=sd(pred),
            q1=mean(low,0.025),
            q3=mean(high,0.975),
            n=n())%>%
  arrange(mean)%>% 
  rename(bme=bmeG,IMD=imd_quintileNew, Mortality=deathNew,YoungAge=young
         ,OldAge=old,urbanicity=urbanicity,LAB=LAB,
         Cons=Cons,EURef=leave,asthma=asthma, bloodPressure=bp,depression=depression,
         diabetes=diabetes,vac_accessibility=WVacAcc
  ) %>%
  ungroup %>%   
  mutate(ID=row_number()) %>%
  tidyr::gather(Covariates,Value,1:13) %>%
  arrange(mean)



## Selecting MASOAs with lower and high probabilities

LH_msoaData <- data1[data1$groupPred == "[0.376,0.724)" | data1$groupPred == "[0.894,0.940]",]

LH_msoaDataLow <- data1[data1$pred <0.724,]

LH_msoaDatahigh <- data1[data1$pred >= 0.894,]

##-----All profiles
Data2<-LH_msoaData %>%
  group_by(bmeG,imd_quintileNew,young,old,deathNew,urbanicity,
           LAB,Cons,leave,asthma, bp, diabetes,depression,WVacAcc) %>%
  summarise(mean=mean(pred),sd=sd(pred),
            q1=mean(low,0.025),
            q3=mean(high,0.975),
            n=n())%>%
  arrange(mean)%>% 
  rename(bme=bmeG,IMD=imd_quintileNew, Mortality=deathNew,YoungAge=young
         ,OldAge=old,urbanicity=urbanicity,LAB=LAB,
         Cons=Cons,EURef=leave,asthma=asthma, bloodPressure=bp,depression=depression,
         diabetes=diabetes,vac_accessibility=WVacAcc
  ) %>%
  ungroup %>%   
  mutate(ID=row_number()) %>%
  tidyr::gather(Covariates,Value,1:14) %>%
  arrange(mean)

###------Low profiles
Data2low<-LH_msoaDataLow %>%
  group_by(bmeG,imd_quintileNew,young,old,deathNew,urbanicity,
           LAB,Cons,leave,asthma, bp, diabetes,depression,WVacAcc) %>%
  summarise(mean=mean(pred),sd=sd(pred),
            q1=mean(low,0.025),
            q3=mean(high,0.975),
            n=n())%>%
  arrange(mean)%>% 
  rename(bme=bmeG,IMD=imd_quintileNew, Mortality=deathNew,YoungAge=young
         ,OldAge=old,urbanicity=urbanicity,LAB=LAB,
         Cons=Cons,EURef=leave,asthma=asthma, bloodPressure=bp,depression=depression,
         diabetes=diabetes,vac_accessibility=WVacAcc
  ) %>%
  ungroup %>%   
  mutate(ID=row_number()) %>%
  tidyr::gather(Covariates,Value,1:14) %>%
  arrange(mean)

####-----High profiles
Data2high<-LH_msoaDatahigh %>%
  group_by(bmeG,imd_quintileNew,young,old,deathNew,urbanicity,
           LAB,Cons,leave,asthma, bp, diabetes,depression,WVacAcc) %>%
  summarise(mean=mean(pred),sd=sd(pred),
            q1=mean(low,0.025),
            q3=mean(high,0.975),
            n=n())%>%
  arrange(mean)%>% 
  rename(BME=bmeG,IMD=imd_quintileNew, Young= young, Mortality=deathNew
         ,OldAge=old,urbanicity=urbanicity,LAB=LAB,
         Cons=Cons,EURef=leave,asthma=asthma, bloodPressure=bp,depression=depression,
         diabetes=diabetes,vac_accessibility=WVacAcc
  ) %>%
  ungroup %>%   
  mutate(ID=row_number()) %>%
  tidyr::gather(Covariates,Value,1:14)%>% 
  arrange(mean)


library(patchwork)
library(forcats)
heatmap<- ggplot(Data2,aes(x=fct_inorder(Covariates),y=ID,fill=Value))+
  geom_tile()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank()) +
  ylab(label = "MSOA") +
  xlab(label = "Profile")+
  scale_fill_gradient(low="#0073C2FF", high="#CD534CFF")+
  scale_x_discrete(
    label = c("BME", "IMD", "Young",
              "Old", "Mortality", "Urbanicity",
              "Labour", "Conservatives", "Brexit",
              "Asthma", "Bood pressure","Diabetes",
              "Depression", "Access"))
x11()
heatmap
ggsave("heatmap.tiff",width = 12,height = 8)

Prob_CI<-ggplot(Data2,aes(x=mean,y=ID)) +
  geom_linerange(aes(xmin=q1  , xmax=q3)) +
  geom_point(aes(x=mean,y=ID)) +
  ylab(label = "") + 
  scale_fill_gradient(low="#56B1F7", high="#132B43")+
  xlab(label = "Vaccine uptake") + 
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background =  element_blank())
Prob_CI

x11()
heatmap|Prob_CI

ggsave("Profile_heatMap.pdf",width = 12,height = 10)

###----Separating Profiles with low probabilities and ones with high probability

##----Low probabilities
heatmaplow<- ggplot(Data2low,aes(x=fct_inorder(Covariates),y=ID,fill=factor(Value)))+
  geom_tile()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank()) +
  ylab(label = "Area") +
  xlab(label = "Profile")+
  #scale_fill_gradient(low="#868686FF" , high="#CD534CFF")+
  #scale_colour_gradientn(colours = terrain.colors(5))+
  #scale_fill_gradient(low = "white", high = "blue")+
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  scale_x_discrete(
    label = c("BME", "IMD", "Young",
              "Old", "Mortality", "Urbanicity",
              "Labour", "Conservatives", "Brexit",
              "Asthma", "Bood pressure","Diabetes",
              "Depression", "Access"))

x11()
heatmaplow

Prob_CIlow<-ggplot(Data2low,aes(x=mean,y=ID)) +
  geom_linerange(aes(xmin=q1  , xmax=q3)) +
  geom_point(aes(x=mean,y=ID)) +
  ylab(label = "") + 
  #scale_fill_gradient(low="#56B1F7", high="#132B43")+
  #scale_colour_gradientn(colours = terrain.colors(5))+
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  xlab(label = "Vaccine uptake") + 
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background =  element_blank())

x11()
Prob_CIlow

heatmaplow|Prob_CIlow
ggsave("heatMap_LowerProbN.jpeg",width = 10,height = 8)

#-----Higher probabilities

heatmaphigh<- ggplot(Data2high,aes(x=fct_inorder(Covariates),y=ID,fill=factor(Value)))+
  geom_tile()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background =  element_blank()) +
  ylab(label = "Area") +
  xlab(label = "Profile")+
  #scale_fill_gradient(low="#56B1F7", high="#132B43")
  #scale_fill_gradient(low = "white", high = "black")+
  scale_fill_viridis_d(option = "viridis",direction = -1,begin = 0.2,end = 1)+
  scale_x_discrete(
    label = c("BME", "IMD", "Young",
              "Old", "Mortality", "Urbanicity",
              "Labour", "Conservatives", "Brexit",
              "Asthma", "Bood pressure","Diabetes",
              "Depression", "Access"))

x11()
heatmaphigh

Prob_CIhigh<-ggplot(Data2high,aes(x=mean,y=ID)) +
  geom_linerange(aes(xmin=q1, xmax=q3)) +
  geom_point(aes(x=mean,y=ID)) +
  ylab(label = "") + 
  scale_fill_gradient(low="#56B1F7", high="#132B43")+
  xlab(label = "Vaccine uptake") + 
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background =  element_blank())

x11()
Prob_CIhigh

heatmaphigh|Prob_CIhigh

ggsave("heatmap_HigherProbN.jpeg",width = 10,height = 8)

########################################################################
###--------Vaccination centers and vaccine access-------------##########
#######################################################################

vw=vac_poly_Category %>% 
  ggplot() + 
  geom_sf(aes(fill = AccWeighted), color = NA, show.legend = TRUE)+
  scale_fill_manual(values = c("#FC4E07","#0072B2" , 
                               "#00AFBB","#F0E442", "#CC79A7","#999999" )) +
  theme_void() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size=9,face = "bold",hjust = 0.5))+
  ggtitle(" Vaccination Access ")

x11()
vw
ggsave("weighted_Vac_Acc.tiff",width = 14,height = 10)


va=ggplot() + 
  geom_sf(data = msoa_poly, fill = NA) +
  geom_point(data = vac_access,
             aes(x=longitude,y=latitude,colour= site),size=1)+
  
  scale_colour_manual(values = c("#FC4E07","#0072B2" , "#F0E442","#CC79A7" ))+
  theme_void() +
  theme(legend.title =  element_blank(),
        plot.title = element_text(size=9, face = "bold",hjust = 0.5))+
  ggtitle("Viccination sites")
va
ggsave("vac_AccessSites.pdf",width = 10,height = 12)

acces=va|vw

x11()

acces +plot_layout(guides = "collect")+
  plot_annotation(title = "Vaccine Accessibility",
                  theme = theme(plot.title = element_text(face = "bold",size = 14)))

ggsave("VacAccesNew.pdf",width=12,height =14)


###########################################################################
############----Fitting Univariate (unadjusted) model-------##############
##########################################################################

#-------------------Univariate model--------------------------#
# With spatial
##---mbe------
formulaNewFirstGbme <- Cum_firstDose~ 1+bmeG+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGbme = inla(formulaNewFirstGbme, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                             control.predictor = list(compute=TRUE),
                             control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_bme <-data.frame(
  M= exp(vac_ful3New2FirstGbme$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGbme$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGbme$summary.fixed$`0.975quant`)[-1]
)

#----------------------------
# imd
formulaNewFirstGimd <- Cum_firstDose~ 1+imd_quintileNew+
  f(reid,model="bym2",graph = g,scale.model =TRUE )

vac_ful3New2FirstGimd = inla(formulaNewFirstGimd, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                             control.predictor = list(compute=TRUE),
                             control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_imd <-data.frame(
  M= exp(vac_ful3New2FirstGimd$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGimd$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGimd$summary.fixed$`0.975quant`)[-1]
)
#------------------------
# leave
formulaNewFirstGleave <- Cum_firstDose~ 1+leave+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGleave = inla(formulaNewFirstGleave, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                               control.predictor = list(compute=TRUE),
                               control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_leave <-data.frame(
  M= exp(vac_ful3New2FirstGleave$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGleave$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGleave$summary.fixed$`0.975quant`)[-1]
)
#--------------------------------------------
# asthma
formulaNewFirstGasthma <- Cum_firstDose~ 1+asthma+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGasthma = inla(formulaNewFirstGasthma, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                                control.predictor = list(compute=TRUE),
                                control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_asthma <-data.frame(
  M= exp(vac_ful3New2FirstGasthma$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGasthma$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGasthma$summary.fixed$`0.975quant`)[-1]
)

#----------------------------------
# diabetes
formulaNewFirstGdiabetes <- Cum_firstDose~ 1+diabetes+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGdiabetes = inla(formulaNewFirstGdiabetes, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                                  control.predictor = list(compute=TRUE),
                                  control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_diabetes <-data.frame(
  M= exp(vac_ful3New2FirstGdiabetes$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGdiabetes$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGdiabetes$summary.fixed$`0.975quant`)[-1]
)

#--------------------------------------
# bp
formulaNewFirstGbp <- Cum_firstDose~ 1+bp+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGbp = inla(formulaNewFirstGbp, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                            control.predictor = list(compute=TRUE),
                            control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_bp <-data.frame(
  M= exp(vac_ful3New2FirstGbp$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGbp$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGbp$summary.fixed$`0.975quant`)[-1]
)

#-------------------------------------------
# depression
formulaNewFirstGdepression <- Cum_firstDose~ 1+depression+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGdepression = inla(formulaNewFirstGdepression, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                                    control.predictor = list(compute=TRUE),
                                    control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_depression <-data.frame(
  M= exp(vac_ful3New2FirstGdepression$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGdepression$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGdepression$summary.fixed$`0.975quant`)[-1]
)
#-----------------------------------------
# Cons
formulaNewFirstGcons <- Cum_firstDose~ 1+Cons+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGcons = inla(formulaNewFirstGcons, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                              control.predictor = list(compute=TRUE),
                              control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_cons <-data.frame(
  M= exp(vac_ful3New2FirstGcons$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGcons$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGcons$summary.fixed$`0.975quant`)[-1]
)

#--------------------------------------
# LAB
formulaNewFirstGlab <- Cum_firstDose~ 1+LAB+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGlab = inla(formulaNewFirstGlab, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                             control.predictor = list(compute=TRUE),
                             control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_lab <-data.frame(
  M= exp(vac_ful3New2FirstGlab$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGlab$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGlab$summary.fixed$`0.975quant`)[-1]
)
#---------------------------------------
# death
formulaNewFirstGdeath <- Cum_firstDose~ 1+death+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGdeath = inla(formulaNewFirstGdeath, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                               control.predictor = list(compute=TRUE),
                               control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_death <-data.frame(
  M= exp(vac_ful3New2FirstGdeath$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGdeath$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGdeath$summary.fixed$`0.975quant`)[-1]
)

#-------------------------------------------
# young
formulaNewFirstGyoung <- Cum_firstDose~ 1+young+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGyoung = inla(formulaNewFirstGyoung, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                               control.predictor = list(compute=TRUE),
                               control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_young <-data.frame(
  M= exp(vac_ful3New2FirstGyoung$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGyoung$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGyoung$summary.fixed$`0.975quant`)[-1]
)
#----------------------------------------
# urbanicity
formulaNewFirstGur <- Cum_firstDose~ 1+urbanicity+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGur = inla(formulaNewFirstGur, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                            control.predictor = list(compute=TRUE),
                            control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_ur <-data.frame(
  M= exp(vac_ful3New2FirstGur$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGur$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGur$summary.fixed$`0.975quant`)[-1]
)
#---------------------------
# old
formulaNewFirstGold <- Cum_firstDose~ 1+old+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGold = inla(formulaNewFirstGold, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                             control.predictor = list(compute=TRUE),
                             control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_old <-data.frame(
  M= exp(vac_ful3New2FirstGold$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGold$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGold$summary.fixed$`0.975quant`)[-1]
)
#-----------------------
## WVA

formulaNewFirstGWVA <- Cum_firstDose~ 1+WVA+
  f(reid,model="bym2",graph = g,scale.model =TRUE )


vac_ful3New2FirstGWVA = inla(formulaNewFirstGWVA, family = "binomial", data = vac_poly_Category, Ntrials=PopVac,
                             control.predictor = list(compute=TRUE),
                             control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE,config = TRUE))

df_WVA <-data.frame(
  M= exp(vac_ful3New2FirstGWVA$summary.fixed$`0.5quant`)[-1],
  Lower = exp(vac_ful3New2FirstGWVA$summary.fixed$`0.025quant`)[-1],
  Upper= exp(vac_ful3New2FirstGWVA$summary.fixed$`0.975quant`)[-1]
)


df_univ <- rbind(df_bme,df_imd,df_young,df_old,df_ur,df_lab,df_cons,
                 df_leave,df_death,df_asthma,df_bp,df_diabetes,
                 df_depression,df_WVA)

New= c("1 (low)","2","3","4","5",
       "1 (least)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "Urban","Rural","Semi",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5",
       "1 (low)","2","3","4","5"
)
A <- c(1,6,11,16,21,24,29,34,39,44,49,54,59,64)
B <- c(2:5,7:10,12:15,17:20,22:23,25:28,30:33,35:38,40:43,45:48,50:53,55:58,60:63,65:68)
N<- rep(1,14)

Mu <- c()
Mu[A]<- N
Mu[B] <- df_univ$M

Lu <- c()
Lu[A]<- N
Lu[B] <- df_univ$Lower

Uu <- c()
Uu[A]<- N
Uu[B] <- df_univ$Upper

df_univNew <- data.frame(index=c(1:23, 25:39, 41:65, 67:71),variable=New,M=Mu,Lower=Lu,Upper=Uu)



df_univNew$class <- ifelse(df_univNew$index <24, "Socio_demographic and urbanicity",
                           ifelse(df_univNew$index >= 24 & df_univNew$index <= 40, "Political_belief",
                                  ifelse(df_univNew$index > 40 & df_univNew$index <= 65,
                                         "COVID-19 awareness and targeting of high risk groups",
                                         "Vaccine_Accessibility")))


df_univNew$class <- factor(df_univNew$class,levels = c("Socio_demographic and urbanicity","Political_belief",
                                                       "COVID-19 awareness and targeting of high risk groups",
                                                       "Vaccine_Accessibility" ))

breaks= seq(0.3,2.8,by=0.1)

labels = as.character(breaks)

pmu <- ggplot(data=df_univNew)+
  geom_vline(xintercept=1, linetype="dashed") +
  geom_errorbar(aes(y=index, x=M, xmin=Lower, xmax=Upper))+
  #scale_color_manual("class",values=c( "#CC79A7","#009E73","#0072B2","#D55E00"))+
  #this adds the effect sizes to the plot
  geom_point(aes(x=M,y=index))+
  #adds the CIs
  scale_y_continuous(name = "", breaks=c(1:23, 25:39, 41:65, 67:71),
                     labels = df_univNew$variable, trans="reverse",
                     expand = c(0, 0)) +
  #adding a vertical line at the effect = 0 mark
  #thematic stuff
  ggtitle("")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.y = element_markdown(size = rel(1.2)),
    axis.title.y = element_text(size = rel(2) ),
    panel.spacing.y  = unit(2, "lines"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0,1, 0, 1), "cm")
  ) +
  xlab("Odds ratio")  +
  geom_hline(yintercept = c(5.5, 10.5, 15.5, 20.5, 24,
                            29.5, 34.5, 40, 45.5,
                            50.5, 55.5, 60.5, 66), col = "grey85")+
  geom_hline(yintercept = c(24, 40, 66), col = "black") +
  geom_vline(xintercept = breaks,
             col = "grey85",linetype="dashed") +
  scale_x_continuous(labels = labels,breaks = breaks)+
  annotate("text",
           x = rep(-0.09, lenght.out = 14),angle = 0,
           #y = c(69, 64, 57, 54, 48, 44, 40, 35, 30, 24, 19, 14, 9, 3),
           y=c(3,8,13,18,22,27,32,37,43,48,53,58,63,69),
           label = c("BME", "IMD", "Young",
                     "Old", "Urbanicity",
                     "Labour", "Conservatives", "Brexit", "Mortality",
                     "Asthma", "Bood pressure","Diabetes",
                     "Depression", "Access")) +
  annotate("text",
           x = rep(2.99, lenght.out =4 ),angle = 3*90,
           y = c(10, 33, 53, 69),
           label = c("Socio-demographics\n and urbanicity",
                     "Political opinions",
                     "COVID-19 awareness and \n targeting of high risk groups",
                     "Access"), fontface = "bold", size =3.5
           #,colour=c( "#CC79A7","#009E73","#0072B2","#D55E00")
  ) +
  geom_vline(xintercept=1, linetype="dashed",col="red") +
  coord_cartesian(xlim = c(0.28, 2.8), clip = "off")

x11()
pmu

ggsave("CovariatesEstimatesUni3.pdf", width = 11, height = 14)
