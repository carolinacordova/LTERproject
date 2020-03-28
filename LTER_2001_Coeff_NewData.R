
setwd("C:/Users/scm_c/Dropbox/GLR_MSU_postdoc/LTER KBS Project/Deep_soilcores/Calculations/Final_Calculations")
rm(list=ls())

library(tidyverse)
library(dplyr,warn.conflicts = FALSE)

packageVersion("dplyr")

df <-read.csv("LTERData_2001_ReRanSamples.csv", comment.char='#')
df$station <- as.factor(df$station)
df$replicate <- as.factor(df$replicate)
df$core <- as.factor(df$core)

## Non-linear Functions ####
## Exponential decay, 3 parameters function ####
ExpoD3p <- function (X, # explanatory variable (e.g. depth)
                     yb = -0.07, #variable (e.g., %C, or %N) at the base of the profile,         
                     A=2.8, # Difference between variable at the surface minus at the base of the profile a= Yo - Yb
                     k = 0.03 # det. the steepness of the curve
){
  yb+A*exp(-k*X)
}


##Exponential rise function ####
ExpoRise <- function (X, #explanatory variably (depth)
                      Bdmin = 1.6, # initial or maximum value of Bulk density,
                      b=0.02 # rate constant det. the steepness of the curve)
){
  Bdmin*(1-exp(-b*X))
}



##GRaphs for all treatments T1-T7
dff <- df[df$treatment %in% c("1","CF","DF","SF"),]
df_2 <- dff %>% 
  drop_na()


# Soil D.w.
df_2%>%
  ggplot(aes(depth,cum_soil_dw))+
  geom_point(data=dff, aes(shape= replicate, color=core),size=3)+
  geom_smooth(aes(), method="lm",
              formula = y ~ 0+x)+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),expand=c(0,0),
                     name="Depth (cm)")+
  scale_y_continuous(limits= c(0,2250),breaks = seq(0,2250,250),expand=c(0,0),
                     name= "(kg/m2)")+
  facet_grid(.~treatment)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 Soil d.w.")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_cum_soildw.png",width = 16,height = 5,units = "in")

# Bulk Density Gravel free graph
df_2%>%
  ggplot(aes(depth, gavel_free_Bd))+
  geom_point(data=dff, aes(shape=replicate, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoRise(x,Bdmin,b), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(Bdmin=1.34, b=0.17),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0.5,2),breaks = seq(0.5,2,0.25),expand=c(0,0),
                     name= "(g/cm3)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),expand=c(0,0),
                     name="Depth (cm)")+
  facet_grid(.~treatment)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 gravel-free Bulk density")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_bd.png",width = 16,height = 5,units = "in")

# Soil %OC graph

df_2%>%
  ggplot(aes(depth,c_prct))+
  geom_point(data=dff, aes(shape= replicate, color=core),size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb,A,k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb=-0.07,A=2.8, k=0.03),
                                 control=nls.control(maxiter = 1000,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),expand=c(0,0),
                     name="Depth (cm)")+
  #  scale_y_continuous(limits= c(0,2.5),breaks = seq(0,2.5,0.5),expand=c(0,0),
  #                    name= "(%)")+
  facet_grid(.~treatment)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 Soil Organic Carbon")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2013_SOCprct_v2.png",width = 16,height = 5,units = "in")

# Soil %N graph

df_2%>%
  ggplot(aes(depth,n_prct))+
  geom_point(data=dff, aes(shape= replicate, color=core),size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb,A,k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb=-0.07,A=2.8, k=0.03),
                                 control=nls.control(maxiter = 1000,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),expand=c(0,0),
                     name="Depth (cm)")+
  # scale_y_continuous(limits= c(0,0.25),breaks = seq(0,0.35,0.05),expand=c(0,0),
  #                   name= "(%)")+
  facet_grid(.~treatment)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 Soil Organic Nitrogen")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_SONprct_v2.png",width = 16,height = 5,units = "in")

##____________Treatment 1 _____________####
##1) Cumulative SOIL D.W. -->Fitting linear regression to Cumulative Soil dw, intercept=0 
df1 <- filter(df_2,treatment== '1')

linearMod <- lm(cum_soil_dw ~ 0+depth, data=df1)  # build linear regression model on full data
summary(linearMod)

df1%>%
  ggplot(aes(depth,cum_soil_dw))+
  geom_point(data=df1, aes(shape= replicate, color=core),size=3)+
  geom_smooth(aes(), method="lm", 
              formula = y ~ 0+x)+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25), expand=c(0,0))+
  scale_y_continuous(limits= c(0,2250),breaks = seq(0,2250,250), expand=c(0,0))+
  labs(x="Depth (cm)", y=expression(Cumulative~Soil~dw.~(kg~m^{-2})))+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 T1 Cum. Soil d.w.")
ggsave("2001_T1_cumsoildw.png",width = 7,height = 5,units = "in")


### 2) Obtaining initial coefficients for Exponential Decay fitted to SOC% 
reg.ExpoD1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                  data=df1,
                  start=c(yb= 0.0565, A= 1.25, k=0.032),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %C graph per replicate 

df1%>%
  ggplot(aes(depth, c_prct))+
  geom_point(data=df1, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb,A,k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb=0.0557, A= 1.25, k=0.032),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,2),breaks = seq(0,2,0.25),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 T1 Soil organic carbon")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_T1_SOCprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients T1 SOC 

df1 %>% split(.$replicate)%>% 
  map( ~nls(c_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb=0.0557, A= 1.25, k=0.032)))%>% 
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- df1 %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.028, A= 1.098, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- df1 %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 4
df_r4 <- df1 %>%
  dplyr::filter(replicate==("4"))

reg.ExpoD_r4 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r4,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r4)

# Rep 5
df_r5 <- df1 %>%
  dplyr::filter(replicate==("5"))

reg.ExpoD_r5 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r5,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r5)

# Rep 6
df_r6 <- df1 %>%
  dplyr::filter(replicate==("6"))

reg.ExpoD_r6 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r6,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r6)

##3) GRAVEL-FREE BULK DENSITY fitting Exponential maximum rise function per replicate

##Obtaining initial coefficients for Exponential maximum rise fitted to Bulk density per Treatment
reg.ExpoRise1 <- nls(gavel_free_Bd ~ ExpoRise(depth,Bdmin,b),
                     data=df1,
                     start=c(Bdmin=1.57, b=0.25),
                     control=nls.control(maxiter= 100, tol=1e-3))

summary(reg.ExpoRise1)

# Bulk Density Gravel free graph
df1%>%
  ggplot(aes(depth, gavel_free_Bd))+
  geom_point(data=df1, aes(shape=replicate, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoRise(x,Bdmin,b), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(Bdmin=1.34, b=0.17),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0.5,2),breaks = seq(0.5,2,0.25),
                     name= "(g/cm3)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 T1 gravel-free Bulk density")

ggsave("2001_T1bd.png",width = 5,height = 5,units = "in")


### 4) Obtaining initial coefficients for Exponential Decay fitted to Soil %N ####
reg.ExpoD1 <- nls(n_prct  ~ ExpoD3p(depth,yb, A, k),
                  data=df1,
                  start=c(yb= 0.006, A= 0.103, k=0.023),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %N graph per replicate 

df1%>%
  ggplot(aes(depth, n_prct))+
  geom_point(data=df1, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb, A, k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.006, A= 0.103, k=0.023),
                                 control=nls.control(maxiter = 1000,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,0.25),breaks = seq(0,0.25,0.05),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 T1 Soil Nitrogen")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_T1_SONprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients T1 SON 

df1 %>% split(.$replicate)%>% 
  map( ~nls(n_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.006, A= 0.103, k=0.023)))%>%  
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- df1 %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- df1 %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 4
df_r4 <- df1 %>%
  dplyr::filter(replicate==("4"))

reg.ExpoD_r4 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r4,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r4)

# Rep 5
df_r5 <- df1 %>%
  dplyr::filter(replicate==("5"))

reg.ExpoD_r5 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r5,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r5)

# Rep 6
df_r6 <- df1 %>%
  dplyr::filter(replicate==("6"))

reg.ExpoD_r6 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r6,
                    start=c(yb= 0.038, A= 1.090, k=0.030),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r6)



##____________Treatment CF _____________####
##1) Cumulative SOIL D.W. -->Fitting linear regression to Cumulative Soil dw, intercept=0 ####
dfcf <- filter(df_2,treatment== 'CF')

linearMod <- lm(cum_soil_dw ~ 0+depth, data=dfcf)  # build linear regression model on full data
summary(linearMod)

dfcf%>%
  ggplot(aes(depth,cum_soil_dw))+
  geom_point(data=dfcf, aes(shape= replicate, color=core),size=3)+
  geom_smooth(aes(), method="lm", 
              formula = y ~ 0+x)+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25), expand=c(0,0))+
  scale_y_continuous(limits= c(0,2250),breaks = seq(0,2250,250), expand=c(0,0))+
  labs(x="Depth (cm)", y=expression(Cumulative~Soil~dw.~(kg~m^{-2})))+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 CF Cum. Soil d.w.")
ggsave("2001_CF_cumsoildw.png",width = 7,height = 5,units = "in")


### 2) Obtaining initial coefficients for Exponential Decay fitted to SOC% ####
reg.ExpoD1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                  data=dfcf,
                  start=c(yb= 0.025, A= 3.33, k=0.053),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %C graph per replicate 

dfcf%>%
  ggplot(aes(depth, c_prct))+
  geom_point(data=dfcf, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb,A,k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.025, A= 3.33, k=0.053),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,2),breaks = seq(0,2,0.25),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 CF Soil organic carbon")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_CF_SOCprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients T1 SOC 

dfcf %>% split(.$replicate)%>% 
  map( ~nls(c_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.025, A= 3.33, k=0.053)))%>% 
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- dfcf %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- dfcf %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 3
df_r3 <- dfcf %>%
  dplyr::filter(replicate==("3"))

reg.ExpoD_r4 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r3,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r4)

##3) GRAVEL-FREE BULK DENSITY fitting Exponential maximum rise function per replicate ####

##Obtaining initial coefficients for Exponential maximum rise fitted to Bulk density per Treatment
reg.ExpoRise1 <- nls(gavel_free_Bd ~ ExpoRise(depth,Bdmin,b),
                     data=dfcf,
                     start=c(Bdmin=1.57, b=0.25),
                     control=nls.control(maxiter= 100, tol=1e-3))

summary(reg.ExpoRise1)

# Bulk Density Gravel free graph
dfcf%>%
  ggplot(aes(depth, gavel_free_Bd))+
  geom_point(data=dfcf, aes(shape=replicate, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoRise(x,Bdmin,b), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(Bdmin=1.34, b=0.17),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0.5,2),breaks = seq(0.5,2,0.25),
                     name= "(g/cm3)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 CF gravel-free Bulk density")

ggsave("2001_CFbd.png",width = 5,height = 5,units = "in")


### 4) Obtaining initial coefficients for Exponential Decay fitted to Soil %N ####
reg.ExpoD1 <- nls(n_prct  ~ ExpoD3p(depth,yb, A, k),
                  data=dfcf,
                  start=c(yb= 0.008, A= 0.189, k=0.047),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %N graph per replicate 

dfcf%>%
  ggplot(aes(depth, n_prct))+
  geom_point(data=dfcf, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb, A, k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.008, A= 0.189, k=0.047),
                                 control=nls.control(maxiter = 1000,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,0.25),breaks = seq(0,0.25,0.05),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 CF Soil Nitrogen")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_CF_SONprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients T1 SON 

dfcf %>% split(.$replicate)%>% 
  map( ~nls(n_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.008, A= 0.189, k=0.047)))%>%  
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- dfcf %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.008, A= 0.189, k=0.047),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- dfcf %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.008, A= 0.189, k=0.047),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 3
df_r3 <- dfcf %>%
  dplyr::filter(replicate==("3"))

reg.ExpoD_r3 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r3,
                    start=c(yb= 0.008, A= 0.189, k=0.047),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r3)




##____________Treatment DF _____________####
##1) Cumulative SOIL D.W. -->Fitting linear regression to Cumulative Soil dw, intercept=0 ####
dfdf <- filter(df_2,treatment== 'DF')

linearMod <- lm(cum_soil_dw ~ 0+depth, data=dfdf)  # build linear regression model on full data
summary(linearMod)

dfdf%>%
  ggplot(aes(depth,cum_soil_dw))+
  geom_point(data=dfdf, aes(shape= replicate, color=core),size=3)+
  geom_smooth(aes(), method="lm", 
              formula = y ~ 0+x)+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25), expand=c(0,0))+
  scale_y_continuous(limits= c(0,2250),breaks = seq(0,2250,250), expand=c(0,0))+
  labs(x="Depth (cm)", y=expression(Cumulative~Soil~dw.~(kg~m^{-2})))+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 DF Cum. Soil d.w.")
ggsave("2001_DF_cumsoildw.png",width = 7,height = 5,units = "in")


### 2) Obtaining initial coefficients for Exponential Decay fitted to SOC% ####
reg.ExpoD1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                  data=dfdf,
                  start=c(yb= 0.005, A= 3.19, k=0.039),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %C graph per replicate 

dfdf%>%
  ggplot(aes(depth, c_prct))+
  geom_point(data=dfdf, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb,A,k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.005, A= 3.19, k=0.039),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,2),breaks = seq(0,2,0.25),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 DF Soil organic carbon")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_DF_SOCprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients SOC 

dfdf %>% split(.$replicate)%>% 
  map( ~nls(c_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.025, A= 3.33, k=0.053)))%>% 
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- dfdf %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- dfdf %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 3
df_r3 <- dfdf %>%
  dplyr::filter(replicate==("3"))

reg.ExpoD_r4 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r3,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r4)

##3) GRAVEL-FREE BULK DENSITY fitting Exponential maximum rise function per replicate ####
#Obtaining initial coefficients for Exponential maximum rise fitted to Bulk density per Treatment
reg.ExpoRise1 <- nls(gavel_free_Bd ~ ExpoRise(depth,Bdmin,b),
                     data=dfdf,
                     start=c(Bdmin=1.57, b=0.25),
                     control=nls.control(maxiter= 100, tol=1e-3))

summary(reg.ExpoRise1)

# Bulk Density Gravel free graph
dfdf%>%
  ggplot(aes(depth, gavel_free_Bd))+
  geom_point(data=dfdf, aes(shape=replicate, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoRise(x,Bdmin,b), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(Bdmin=1.34, b=0.17),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0.5,2),breaks = seq(0.5,2,0.25),
                     name= "(g/cm3)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 DF gravel-free Bulk density")

ggsave("2001_DFbd.png",width = 5,height = 5,units = "in")


### 4) Obtaining initial coefficients for Exponential Decay fitted to Soil %N ####
reg.ExpoD1 <- nls(n_prct  ~ ExpoD3p(depth,yb, A, k),
                  data=dfdf,
                  start=c(yb= 0.002, A= 0.223, k=0.034),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %N graph per replicate 

dfdf%>%
  ggplot(aes(depth, n_prct))+
  geom_point(data=dfdf, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb, A, k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.002, A= 0.223, k=0.034),
                                 control=nls.control(maxiter = 1000,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,0.25),breaks = seq(0,0.25,0.05),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 DF Soil Nitrogen")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_DF_SONprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients T1 SON 

dfdf %>% split(.$replicate)%>% 
  map( ~nls(n_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.002, A= 0.223, k=0.034)))%>%  
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- dfdf %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.002, A= 0.223, k=0.034),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- dfdf %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.002, A= 0.223, k=0.034),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 3
df_r3 <- dfdf %>%
  dplyr::filter(replicate==("3"))

reg.ExpoD_r3 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r3,
                    start=c(yb= 0.002, A= 0.223, k=0.034),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r3)

##____________Treatment SF _____________####
##1) Cumulative SOIL D.W. -->Fitting linear regression to Cumulative Soil dw, intercept=0 ####
dfsf <- filter(df_2,treatment== 'SF')

linearMod <- lm(cum_soil_dw ~ 0+depth, data=dfsf)  # build linear regression model on full data
summary(linearMod)

dfsf%>%
  ggplot(aes(depth,cum_soil_dw))+
  geom_point(data=dfsf, aes(shape= replicate, color=core),size=3)+
  geom_smooth(aes(), method="lm", 
              formula = y ~ 0+x)+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25), expand=c(0,0))+
  scale_y_continuous(limits= c(0,2250),breaks = seq(0,2250,250), expand=c(0,0))+
  labs(x="Depth (cm)", y=expression(Cumulative~Soil~dw.~(kg~m^{-2})))+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 SF Cum. Soil d.w.")
ggsave("2001_SF_cumsoildw.png",width = 7,height = 5,units = "in")


### 2) Obtaining initial coefficients for Exponential Decay fitted to SOC% ####
reg.ExpoD1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                  data=dfsf,
                  start=c(yb= 0.078, A= 3.12, k=0.073),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %C graph per replicate 

dfsf%>%
  ggplot(aes(depth, c_prct))+
  geom_point(data=dfsf, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb,A,k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.078, A= 3.12, k=0.073),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,2),breaks = seq(0,2,0.25),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 SF Soil organic carbon")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_SF_SOCprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients SOC 

dfsf %>% split(.$replicate)%>% 
  map( ~nls(c_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.025, A= 3.33, k=0.053)))%>% 
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- dfsf %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- dfsf %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 3
df_r3 <- dfsf %>%
  dplyr::filter(replicate==("3"))

reg.ExpoD_r4 <- nls(c_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r3,
                    start=c(yb= 0.025, A= 3.33, k=0.053),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r4)

##3) GRAVEL-FREE BULK DENSITY fitting Exponential maximum rise function per replicate ####
#Obtaining initial coefficients for Exponential maximum rise fitted to Bulk density per Treatment
reg.ExpoRise1 <- nls(gavel_free_Bd ~ ExpoRise(depth,Bdmin,b),
                     data=dfsf,
                     start=c(Bdmin=1.47, b=0.08),
                     control=nls.control(maxiter= 100, tol=1e-3))

summary(reg.ExpoRise1)

# Bulk Density Gravel free graph
dfsf%>%
  ggplot(aes(depth, gavel_free_Bd))+
  geom_point(data=dfsf, aes(shape=replicate, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoRise(x,Bdmin,b), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(Bdmin=1.47, b=0.08),
                                 control=nls.control(maxiter = 100,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0.5,2),breaks = seq(0.5,2,0.25),
                     name= "(g/cm3)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 SF gravel-free Bulk density")

ggsave("2001_SFbd.png",width = 5,height = 5,units = "in")


### 4) Obtaining initial coefficients for Exponential Decay fitted to Soil %N ####
reg.ExpoD1 <- nls(n_prct  ~ ExpoD3p(depth,yb, A, k),
                  data=dfsf,
                  start=c(yb= 0.017, A= 0.2, k=0.056),
                  control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD1)

#Soil Organic %N graph per replicate 

dfsf%>%
  ggplot(aes(depth, n_prct))+
  geom_point(data=dfsf, aes(shape=station, color=core), size=3)+
  geom_smooth(method = "nls",
              formula=y~ExpoD3p(x,yb, A, k), # this is an nls argument, 
              #but stat_smooth passes the parameter along
              method.args = list(start=c(yb= 0.017, A= 0.2, k=0.056),
                                 control=nls.control(maxiter = 1000,
                                                     tol = 1e-3)), 
              se=FALSE) +
  scale_y_continuous(limits= c(0,0.25),breaks = seq(0,0.25,0.05),
                     name= "(%)")+
  scale_x_continuous(limits= c(0,125),breaks = seq(0,125,25),
                     name="Depth (cm)")+
  facet_grid(.~replicate)+
  ggthemes::theme_base() + theme(plot.background = element_blank(),
                                 legend.position = "right")+ 
  labs(title="2001 SF Soil Nitrogen")+
  theme(panel.spacing = unit(1.1, "lines"))

ggsave("2001_SF_SONprct.png",width = 16,height = 4,units = "in")

## Getting Coefficients SON 

dfsf %>% split(.$replicate)%>% 
  map( ~nls(n_prct ~ yb + A *(exp(-k*depth)), data = ., start = list(yb= 0.002, A= 0.223, k=0.034)))%>%  
  map(summary)%>% 
  map("coefficients") 

## Per Rep 
# Rep 1
df_r1 <- dfsf %>%
  dplyr::filter(replicate==("1"))

reg.ExpoD_r1 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r1,
                    start=c(yb= 0.002, A= 0.223, k=0.034),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r1)

# Rep 2
df_r2 <- dfsf %>%
  dplyr::filter(replicate==("2"))

reg.ExpoD_r2 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r2,
                    start=c(yb= 0.002, A= 0.223, k=0.034),
                    control=nls.control(maxiter= 1000, tol=1e-3))

summary(reg.ExpoD_r2)

# Rep 3
df_r3 <- dfsf %>%
  dplyr::filter(replicate==("3"))

reg.ExpoD_r3 <- nls(n_prct  ~ ExpoD3p(depth,yb,A,k),
                    data=df_r3,
                    start=c(yb= 0.002, A= 0.223, k=0.034),
                    control=nls.control(maxiter= 1000, tol=1e-3))
summary(reg.ExpoD_r3)