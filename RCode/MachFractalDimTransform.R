
# Title:       Transform the D(l) curve in M(D), for Beattie et al. 2019b
# Author:      James Beattie
# Created:     08/07/18

# Imports
########################################################################################################################

library(tidyverse)
library(gridExtra)
library(zoo)
library(DataCombine)
library(magrittr)
library(ggridges)
library(broom)
library(VGAM)
library(nlstools)

########################################################################################################################


# Global Functions
########################################################################################################################

options(warn=-1)

trunc <- function(x, ..., prec = 0) {
     base::trunc(x * 10^prec, ...) / 10^prec
}

ExponentLabel <- function(l) { 
     l <- paste('10^',l) 
     l <- gsub(" ","", l)
     l <- gsub("e","10^", l) 
     parse(text=l) 
}


FittingBuilder <- function(data,lower,upper){
     
     # Change the fitting variable "Fitting" along the cascade
     
     
     for(i in 1:nrow(data)){
          if(data$Length[i] > lower & data$Length[i] < upper){
               data$Fitting[i] <- 'Fitting'
          }
     }
     return(data)
}


nlsFunc2 <- function(x,parms){
     
     # Functional form from Beattie et al. 2019
     
     beta1     <- parms[1]
     beta0     <- parms[2]
     fmax      <- parms[3]
     fmin      <- parms[4]
     
     f <- (fmax-fmin)*(1-erf(beta1*x+beta0)) / 2 + fmin
     
     return(f)
}


dD_errorFunction <- function(l,dt,parms){
     # This funciton calculates the error of the fit, assuming it is independent and Normal
     # Input:
     # l : length scales that are not log transformed
     # dt: time fluctuations at each l
     # parms: model coefficients as a data.frame() object
     # Ouput:
     # error structure in the fractal dimension
     
     # Calculate the parameters
     beta0     <- parms$estimate[1]
     dbeta0    <- parms$std.error[1]
     
     beta1     <- parms$estimate[2]
     dbeta1     <- parms$std.error[2]
     
     Dmax      <- 2
     
     Dmin      <- parms$estimate[3]    
     dDmin     <- parms$std.error[3]
     
     
     #dfdl      <- - ( (Dmax - Dmin)*exp(-( beta1*log(l)/log(10)  + beta0)^2  )*beta1 ) / ( sqrt(pi)*l*log(10) )
     
     # Calculate the derivatives.
     dfdMin    <- - ( 1/2 + 1/2*erf( -beta1*log10(l) + beta0) )
     
     dfdbeta0  <- ( (Dmax - Dmin)*exp(-( beta1*log10(l)  + beta0)^2  )  ) / ( sqrt(pi) )
     
     dfdbeta1  <- - ( (Dmax - Dmin)*exp(-( beta1*log10(l)  + beta0)^2  )*log10(l) ) / ( sqrt(pi) )
     
     dD        <- sqrt( (dfdbeta1*dbeta1)^2 + (dfdbeta0*dbeta0)^2 + (dfdMin*dDmin)^2 )
     
     # The addition of the temporal fluctuations and the error from the model
     dsigma    <- dD + dt
     
     return(dsigma)
}

FittingBuilderSonic<-function(data,lower,upper){
     sonicscale_up       <- 10^(-1.87 + 0.11)*10048
     sonicscale_down     <- 10^(-1.87 - 0.11)*10048
     
     data$Fitting   <- NA
     data$Scale     <- NA
     
     for(i in 1:nrow(data)){
          if(data$Length[i] <= lower){
               data$Scale[i] <- 'Not fitting'
          }
          else if(data$Length[i] > lower & data$Length[i] <= sonicscale_down){
               data$Scale[i] <- 'Fitting'
          }
          else if(data$Length[i] > sonicscale_down & data$Length[i] < sonicscale_up){
               data$Scale[i] <- 'Not fitting'
          }
          else if(data$Length[i] >= sonicscale_up & data$Length[i] < upper){
               data$Scale[i] <- 'Fitting'
          }
          else if(data$Length[i] >= upper){
               data$Scale[i] <- 'Not fitting'
          }
          
          
          if(data$Length[i] > lower & data$Length[i] <= sonicscale_down){
               data$Fitting[i] <- 'Fit1'
          }
          else if(data$Length[i] >= sonicscale_up & data$Length[i] < upper){
               data$Fitting[i] <- 'Fit2'
          }
     }
     return(data)
}


nlsFunc3_mod <- function(x,parms){
     
     # Define the actual model. Have to hard code the parameters to use the inverse function. 
     
     beta0 <- -0.2331088
     beta1 <- 1.1382839
     fmax <- 2
     fmin <- 1.5522796
     
     f <- (fmax-fmin)*(1 - erf(beta1*x+beta0)) / 2 + fmin
     
     return(f)
}


nlsFunc3_upper <- function(x){
     
     # The upper 1sigma error
     
     beta0 <- -0.149476
     beta1 <- 1.275321
     fmax <- 2
     fmin <- 1.686239
     
     f <- (fmax-fmin)*(1 - erf(beta1*x+beta0)) / 2 + fmin
     
     return(f)
}


nlsFunc3_lower <- function(x){
     
     # The lower 1sigma error
     
     beta0 <- -0.2689987
     beta1 <- 1.0910387
     fmax <- 2
     fmin <- 1.4202617
     
     f <- (fmax-fmin)*(1 - erf(beta1*x+beta0)) / 2 + fmin
     
     return(f)
}

########################################################################################################################


# Directories
########################################################################################################################

setwd('/Volumes/JamesBe/PeriodicBoundaryDatasets/PeriodicSingleMaxPixel/')


# Datasets
########################################################################################################################

Mach1SMP_1024       <- read_csv('Mach1_1024_MassLength_Periodic_SingleMaxPixel_OddSample.csv')
Mach4SMP_2512       <- read_csv('Mach4_2512_MassLength_Periodic_SingleMaxPixel_OddSample.csv')
Mach10SMP_1024      <- read_csv('Mach10_1024_MassLength_Periodic_SingleMaxPixel_OddSample.csv')
Mach20SMP_1024      <- read_csv('Mach20_1024_MassLength_Periodic_SingleMaxPixel_OddSample.csv')
Mach40SMP_1024      <- read_csv('Mach40_1024_MassLength_Periodic_SingleMaxPixel_OddSample.csv')
Mach100SMP_1024     <- read_csv('Mach100_1024_MassLength_Periodic_SingleMaxPixel_OddSample.csv')


# Organsise projection data into xy, xz, zy projections

Mach4SMP_2512$Projection <- NA

for(i in 1:nrow(Mach4SMP_2512)){
     if(Mach4SMP_2512$Filename[i] %>% grepl("coldens_x",.,fixed=TRUE) == TRUE){
          Mach4SMP_2512$Projection[i] <- 'x' 
          next
     } else if(Mach4SMP_2512$Filename[i] %>% grepl("coldens_y",.,fixed=TRUE) == TRUE){
          Mach4SMP_2512$Projection[i] <- 'y'
          next
     } else if(Mach4SMP_2512$Filename[i] %>% grepl("coldens_z",.,fixed=TRUE) == TRUE){
          Mach4SMP_2512$Projection[i] <- 'z'
          next
     }
}

# Transform all of the datasets into the Mach 4 frame.
########################################################################################################################

################################
# Parameters and Length Scales 
################################

# Important parameters
MachHighRes    <- 10048
MachLowRes     <- 1024
SonicScale     <- 120.9765 

# RMS Mach numbers calculated from simulation
M1   <- 1.009288
M4   <- 4.127861
M10  <- 10.16454
M20  <- 20.07651
M40  <- 40.20628
M100 <- 100

# Length scales in the Mach 4 frame
l_1       <- (M1/M4)^2
l_10      <- (M10/M4)^2
l_20      <- (M20/M4)^2
l_40      <- (M40/M4)^2
l_100     <- (M100/M4)^2

###################################################################
# Subset all of the data so the turnover time matches Mach4_2512
###################################################################

# Truncate to two orders of precision
Mach1SMP_1024       %<>% mutate(TurnoverTime = trunc(TimeStep*1*2,prec=2))
Mach4SMP_2512       %<>% filter(Projection == 'x') %>% mutate(TurnoverTime = trunc(TimeStep*4*2,prec=2))
Mach10SMP_1024      %<>% mutate(TurnoverTime  = trunc(TimeStep*10*2,prec=2))
Mach20SMP_1024      %<>% mutate(TurnoverTime  = trunc(TimeStep*20*2,prec=2))
Mach40SMP_1024      %<>% mutate(TurnoverTime  = trunc(TimeStep*40*2,prec=2))
Mach100SMP_1024     %<>% mutate(TurnoverTime = trunc(TimeStep*100*2,prec=2))

# Extract unique turnover times
Mach4SMP_2512  %<>% filter(TurnoverTime != 2.0)
uniqueTT       <- Mach4SMP_2512$TurnoverTime %>% unique()

# Matching the turnover times:
Mach1SMP_1024       <- filter(Mach1SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach10SMP_1024      <- filter(Mach10SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach20SMP_1024      <- filter(Mach20SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach40SMP_1024      <- filter(Mach40SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach100SMP_1024     <- filter(Mach100SMP_1024, TurnoverTime %in%  uniqueTT) 

##############################################################################################
# Transform all data into the M4 frame and take averages of the FD across all turnover times
##############################################################################################

Mach4Rel_2512 <- Mach4SMP_2512 %>% 
     filter(Projection == 'x') %>%
     select(X1,BoxCount,Filename,FractalDim,Length,Mass,Std) %>%
     mutate(LengthRel = log10(Length/2512)) %>% 
     group_by(LengthRel) %>%
     summarise(`Fractal Dimension` = mean(FractalDim),
               `std` = sd(FractalDim),
               Length=mean(Length),
               Fitting = 'Not Fitting') %>%
     mutate(`Mach Number` = 'Mach 4')

translation <- log10(1024/(2*10048)) - log10(SonicScale/10048)

Mach1Rel_1024 <- Mach1SMP_1024 %>% 
     mutate(LengthRel = log10(Length/10048) - translation) %>%
     group_by(LengthRel) %>%
     summarise(`Fractal Dimension` = mean(FractalDim),
               `std` = sd(FractalDim),
               Length=mean(Length),
               Fitting = 'Not Fitting') %>%
     mutate(`Mach Number` = 'Mach 1')

Mach10Rel_1024 <- Mach10SMP_1024 %>% 
     mutate(LengthRel = log10((Length/1024)*l_10)) %>%
     group_by(LengthRel) %>%
     summarise(`Fractal Dimension` = mean(FractalDim),
               `std` = sd(FractalDim),
               Length=mean(Length),
               Fitting = 'Not Fitting') %>%
     mutate(`Mach Number` = 'Mach 10')

Mach20Rel_1024 <- Mach20SMP_1024 %>% 
     mutate(LengthRel = log10((Length/1024)*l_20)) %>%      
     group_by(LengthRel) %>%
     summarise(`Fractal Dimension` = mean(FractalDim),
               `std` = sd(FractalDim),
               Length=mean(Length),
               Fitting = 'Not Fitting') %>%
     mutate(`Mach Number` = 'Mach 20')

Mach40Rel_1024 <- Mach40SMP_1024 %>% 
     mutate(LengthRel = log10((Length/1024)*l_40)) %>%      
     group_by(LengthRel) %>%
     summarise(`Fractal Dimension` = mean(FractalDim),
               `std` = sd(FractalDim),
               Length=mean(Length),
               Fitting = 'Not Fitting') %>%
     mutate(`Mach Number` = 'Mach 40')

Mach100Rel_1024 <- Mach100SMP_1024 %>% 
     mutate(LengthRel = log10((Length/1024)*l_100)) %>%      
     group_by(LengthRel) %>%
     summarise(`Fractal Dimension` = mean(FractalDim),
               `std` = sd(FractalDim),
               Length=mean(Length),
               Fitting = 'Not Fitting') %>%
     mutate(`Mach Number` = 'Mach 100')

########################################################################################################################


# Construct Fractal Dimension fits
########################################################################################################################


################################
# Construction of the model
################################

erf  <- function(x) 2 * pnorm(x * sqrt(2)) - 1    # Error function
fmax <- 2                                         # assumed Dmax for the plane

# create forumula for the non-linear regression
formulaExp     <- as.formula(FD  ~ ( (fmax-fmin)*(1-erf(beta1*LR+beta0)) / 2 + fmin ) )


# Fit data only to the arguments in these functions.
Mach1_RelFit_1024        <- FittingBuilder(Mach1Rel_1024,30,max(Mach1Rel_1024$Length)/2)
Mach4_RelFit_2512        <- FittingBuilder(Mach4Rel_2512,7.5,max(Mach4Rel_2512$Length)/2)
Mach10_RelFit_1024       <- FittingBuilder(Mach10Rel_1024,30,max(Mach10Rel_1024$Length)/2)
Mach20_RelFit_1024       <- FittingBuilder(Mach20Rel_1024,30,max(Mach20Rel_1024$Length)/2)
Mach40_RelFit_1024       <- FittingBuilder(Mach40Rel_1024,30,max(Mach40Rel_1024$Length)/2)
Mach100_RelFit_1024      <- FittingBuilder(Mach100Rel_1024,30,max(Mach100Rel_1024$Length)/2)
ModelDat                 <- rbind(Mach1_RelFit_1024,
                                  Mach4_RelFit_2512,
                                  Mach10_RelFit_1024,
                                  Mach20_RelFit_1024,
                                  Mach40_RelFit_1024,
                                  Mach100_RelFit_1024)

# Subset for fitting
ModelDatFit    <- dplyr::filter(ModelDat,Fitting == 'Fitting')
ModelDatNotFit <- dplyr::filter(ModelDat,Fitting == 'Not Fitting')

# Change the names of length_relative and fractal dimension for ease
names(ModelDatFit)[c(1,2)] <- c('LR','FD')

nls1       <- nls(formulaExp, start = list(beta0 = 0,beta1 = 1/2,fmin=1.4),data = ModelDatFit,weights=std)

sum(residuals(nls1)^2)

tidynls <- nls1 %>% tidy()

parms <- c(tidynls$estimate[1],
           tidynls$estimate[2],
           tidynls$estimate[3])


ModelDatFit$`Mach Number`          <- with(ModelDatFit,reorder(`Mach Number`,LR,function(x)-min(x)))
ModelDatFit$FractalDimensionError  <- dD_errorFunction(10^ModelDatFit$LR,ModelDatFit$std,tidynls)

UpperError <- data.frame(LR = ModelDatFit$LR,FD = ModelDatFit$FractalDimensionError + ModelDatFit$FD)
LowerError <- data.frame(LR = ModelDatFit$LR,FD = ModelDatFit$FD - ModelDatFit$FractalDimensionError)


nls1_Upper       <- nls(formulaExp, start = list(beta0 = 0.5,beta1 = 1/2,fmin=1.7),data =UpperError)

tidynls_Upper <- nls1_Upper %>% tidy()
parms_Upper <- c(tidynls_Upper$estimate[1],
           tidynls_Upper$estimate[2],
           tidynls_Upper$estimate[3])
nls1_Lower       <- nls(formulaExp, start = list(beta0 = 0.5,beta1 = 1/2,fmin=1.2),data =LowerError)
tidynls_Lower <- nls1_Lower %>% tidy()
parms_Lower <- c(tidynls_Lower$estimate[1],
           tidynls_Lower$estimate[2],
           tidynls_Lower$estimate[3])


ModelDatFit %<>% mutate(`Mach Number` = factor(`Mach Number`, levels=c("Mach 1","Mach 4","Mach 10", "Mach 20","Mach 40","Mach 100")) )


png("PosterFig3.png", width = 2000, height = 1000, units = 'px', res = 200)
ggplot(aes(x=LR,y=FD),data=ModelDatFit) +
     geom_point(aes(col=`Mach Number`)) +
     geom_ribbon(aes(ymax=FD + std,ymin=FD -std,fill=`Mach Number`),alpha=0.1) +
     geom_point(aes(x=LengthRel,y=`Fractal Dimension`),data=ModelDatNotFit,col='grey')+
     theme_bw() +
     #geom_hline(yintercept = 2,linetype=3,col='black') +
     #geom_hline(yintercept = 1.552280,linetype=3,col='black') +
     scale_y_continuous(breaks = round(seq(1.4,2, by = 0.1),1)) + 
     scale_x_continuous(breaks = round(seq(min(Mach1Rel_1024$LengthRel), 
                                           max(Mach40Rel_2048$LengthRel)+5, by = 0.5),1)-0.1,
                        labels=ExponentLabel) + 
     theme(axis.title = element_text(size = 16), 
           axis.text = element_text(size = 8), 
           axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15)) + 
     labs(x = expression(paste(l/L)),y = 'Fractal Dimension') + 
     theme(legend.text = element_text(size = 10), 
           legend.title = element_text(size = 12,face = "bold")) + 
     theme(legend.background = element_rect(fill = "white"), legend.position = c(0.9, 0.80)) +
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     stat_function(fun=nlsFunc2,args=list(parms=parms),col='blue',linetype=3) +
     stat_function(fun=nlsFunc2,args=list(parms=parms),col='black',size=1.1) +
     stat_function(fun=nlsFunc2,args=list(parms=parms_Lower),col='black',size=0.5,linetype=2) +
     stat_function(fun=nlsFunc2,args=list(parms=parms_Upper),col='black',size=0.5,linetype=2) +
     geom_segment(y = parms[3],yend = parms[3],x=1.5,xend=3.5,linetype=3) +
     geom_segment(y = 2,yend = 2,x=-4.5,xend=-2.5,linetype=3) 

dev.off()


# Matching Structure Functions
############################################################################################################

head(ord2_struc_tot)
ord2_struc_tot <- dplyr::select(ord2_struc,Grid,`Total_sqrt(SF(order=02))`,`Total_SD(order=02)`) %>% 
     filter(Grid<0.5 & Grid > 0)
names(ord2_struc_tot) <- c('Length','Mach Number','Sigma')
ord2_struc_tot %<>% mutate(Type = '2nd Order Structure Function',Length = log10(Length))

Mach4SF<-dplyr::select(Mach4_RelFit_2512, Length = LengthRel,`Mach Number` = `Fractal Dimension`,Sigma = std) %>%
     mutate(Type = 'Fractal Dimension')


SFandFD <-rbind(Mach4SF,ord2_struc_tot)

SFandFD
ggplot(aes(x=Length,y=`Mach Number`),data=SFandFD) +
     geom_point() +
     facet_wrap(~Type,nrow=2,scale='free_y') +
     theme_bw() +
     scale_x_continuous(labels=ExponentLabel) +
     geom_vline(xintercept = -1.87) +
     geom_vline(xintercept = -1.87 + 0.11) + 
     geom_vline(xintercept = -1.87 - 0.11) +
     geom_vline(xintercept = log10(30/10048))


# Plot the Mach 4 data for paper
########################################################################################################################

Mach4_RelFit_10048 <- FittingBuilderSonic(Mach4Rel_10048,30,10048/10)
Mach4_RelFit_10048 %>% head()


Fit1<-Mach4_RelFit_10048 %>% 
     filter(Fitting == 'Fit1') %>% 
     lm(log10(`Fractal Dimension`) ~ LengthRel,.) %>% 
     tidy()

Fit1$estimate[2]

Fit2<-Mach4_RelFit_10048 %>% 
     filter(Fitting == 'Fit2') %>% 
     lm(log10(`Fractal Dimension`) ~ LengthRel,.) %>% 
     tidy()


sonicscale_up       <- 10^(-1.87 + 0.11)*10048
sonicscale_down     <- 10^(-1.87 - 0.11)*10048

Mach4_RelFit_10048 %<>% mutate(Subsonic =  `Fractal Dimension`/(((10^(LengthRel))*10048)^Fit1$estimate[2]),
                               Supersonic =  `Fractal Dimension`/(((10^(LengthRel))*10048)^Fit2$estimate[2]))

a <- ggplot(aes(x=LengthRel,y=log10(Supersonic)),data=Mach4_RelFit_10048) +
     geom_point(aes(col=Scale)) + 
     annotate('text',x=-1.87,y=0.34,label='paste(l[s])%+-%sigma',
              size=8,col='blue',angle=90,parse=TRUE) +
     geom_vline(xintercept = -1.87,linetype =2,col='blue')+
     geom_vline(xintercept = -1.87-0.11,linetype =2,col='red')+
     geom_vline(xintercept = -1.87+0.11,linetype =2,col='red') +
     theme_bw() +
     theme(legend.position = 'none') +
     labs(x = expression(paste('log'['10']*(l/L))),
          y = expression(paste('log'['10']('Fractal Dimension'/(l/L)^-0.04))),
          title= 'Supersonic Power Law') +
     #geom_abline(intercept = Fit1$estimate[1],slope = Fit1$estimate[2]) +
     geom_abline(intercept = 0.366,slope = 0) +
     scale_x_continuous(labels=ExponentLabel) +
     scale_y_continuous(labels=ExponentLabel) +
     theme(axis.text = element_text(size = 12), 
           axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15),
           strip.text.x = element_text(size = 12),
           legend.text = element_text(size = 10), 
           legend.title = element_text(size = 12,face = "bold"),
           axis.title = element_text(size = 16)) #+

b<- ggplot(aes(x=LengthRel,y=log10(Subsonic)),data=Mach4_RelFit_10048) +
     geom_point(aes(col=Scale)) + 
     annotate('text',x=-1.87,y=0.3,label='paste(l[s])%+-%sigma',
              size=8,col='blue',angle=90,parse=TRUE) +
     geom_vline(xintercept = -1.87,linetype =2,col='blue')+
     geom_vline(xintercept = -1.87-0.11,linetype =2,col='red')+
     geom_vline(xintercept = -1.87+0.11,linetype =2,col='red') +
     theme_bw() +
     theme(legend.position = 'none') +
     labs(x = expression(paste('log'['10']*(l/L))),
          y = expression(paste('log'['10']('Fractal Dimension'/(l/L)^-0.01))), 
          title= 'Subsonic Power Law') +
     geom_abline(intercept = 0.309,slope = 0) +
     #geom_abline(intercept = Fit2$estimate[1],slope = Fit2$estimate[2]) +
     scale_x_continuous(labels=ExponentLabel) +
     scale_y_continuous(labels=ExponentLabel) +
     theme(axis.text = element_text(size = 12), 
           axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15),
           strip.text.x = element_text(size = 12),
           legend.text = element_text(size = 10), 
           legend.title = element_text(size = 12,face = "bold"),
           axis.title = element_text(size = 16)) #+

grid.arrange(a,b)

# Create the Mach number curves
##############################################################################################################

rmsM4 = 4.1
ModelDatFit_Mach         <- ModelDatFit %>% mutate(Mach = sqrt(10^LR)*rmsM4)
ModelDatNotFit_Mach      <- ModelDatNotFit %>% mutate(Mach = sqrt(10^LengthRel)*rmsM4)
ModelDatFit_Mach %<>% mutate(`Mach Number` = factor(`Mach Number`, levels=c("Mach 1","Mach 4","Mach 10", "Mach 20","Mach 40","Mach 100")) )

formulaExp2    <- as.formula(FD  ~ ( (fmax-fmin)*(1-VGAM::erf(beta1*log10(Mach)+beta0)) / 2 + fmin ) )
nls2           <- nls(formulaExp2, start = list(beta0 = -1,beta1 = 1/2,fmin=1),data = ModelDatFit_Mach,weights = std)
nls2_tidy      <- tidy(nls2) 

MachFractalDimError <- dD_errorFunction(ModelDatFit_Mach$Mach,ModelDatFit_Mach$std,nls2_tidy)

UpperError <- data.frame(Mach = ModelDatFit_Mach$Mach,FD = MachFractalDimError + ModelDatFit_Mach$FD)
LowerError <- data.frame(Mach = ModelDatFit_Mach$Mach,FD = ModelDatFit_Mach$FD - MachFractalDimError)

nls2_Upper     <- nls(formulaExp2, start = list(beta0 = -1,beta1 = 1/2,fmin=1),data = UpperError)
Upper_tidy     <- tidy(nls2_Upper) 
parms_Upper_Mach <- c(Upper_tidy$estimate[1],
                 Upper_tidy$estimate[2],
                 Upper_tidy$estimate[3])


nls2_Lower     <- nls(formulaExp2, start = list(beta0 = -1,beta1 = 1/2,fmin=1),data = LowerError)
Lower_tidy     <- tidy(nls2_Lower) 
parms_Lower_Mach <- c(Lower_tidy$estimate[1],
                 Lower_tidy$estimate[2],
                 Lower_tidy$estimate[3])


summary(nls2)
sum(residuals(nls2)^2)

tidynls2 <- nls2 %>% tidy()
tidynls2
parms2 <- c(tidynls2$estimate[1],
           tidynls2$estimate[2],
           tidynls2$estimate[3])

parms2

ggplot(aes(x=log10(Mach),y=FD),data=ModelDatFit_Mach) +
     geom_point(aes(col=`Mach Number`)) +
     geom_ribbon(aes(ymax=FD + std,ymin=FD -std,fill=`Mach Number`),alpha=0.1) +
     geom_point(aes(x=log10(Mach),y=`Fractal Dimension`),data=ModelDatNotFit_Mach,col='grey')+
     theme_bw() +
     scale_x_continuous(breaks = seq(-1.5, 2, by = 0.25),labels=ExponentLabel) +
     theme(axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15),
           axis.title = element_text(size = 16)) + 
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     labs(x = 'Mach Number', 
          y = 'Fractal Dimension') +
     theme(legend.text = element_text(size = 10), 
           legend.title = element_text(size = 12,face = "bold")) + 
     theme(legend.background = element_rect(fill = "white"), legend.position = c(0.9, 0.80)) +
     stat_function(fun=nlsFunc2,args=list(parms=parms2),col='black',size=1.1) + 
     stat_function(fun=nlsFunc2,args=list(parms=parms_Lower_Mach),col='black',size=0.5,linetype=2) +
     stat_function(fun=nlsFunc2,args=list(parms=parms_Upper_Mach),col='black',size=0.5,linetype=2) 


# Find the inverse of the current model
##############################################################################################################

# Invert everything using the inverse function.
invfunc_mod    <- GoFKernel::inverse(nlsFunc3_mod)
invfunc_upper  <- GoFKernel::inverse(nlsFunc3_upper)
invfunc_lower  <- GoFKernel::inverse(nlsFunc3_lower)

# Find the limits of the model, and the 1 sigma fluctuations.
ModelLimit <- nlsFunc3_mod(max(ModelDatFit_Mach$Mach))
UpperLimit <- nlsFunc3_upper(max(ModelDatFit_Mach$Mach))
LowerLimit <- nlsFunc3_lower(max(ModelDatFit_Mach$Mach))

solvec <- rep(0,length(seq(ModelLimit,2,0.001)))
solvec_L <- rep(0,length(seq(LowerLimit,2,0.001)))
solvec_U <- rep(0,length(seq(UpperLimit,2,0.0001)))

counter = 1
for (i in seq(ModelLimit,2,0.001)){
     solvec[counter] <- invfunc_mod(i)
     counter = counter + 1
}

counter = 1
for (i in seq(UpperLimit,2,0.0001)){
     solvec_U[counter] <- invfunc_upper(i)
     counter = counter + 1
}

counter = 1
for (i in seq(LowerLimit,2,0.001)){
     solvec_L[counter] <- invfunc_lower(i)
     counter = counter + 1
}





InverseModel <- data.frame(cbind(FD = seq(ModelLimit,2,0.001),Mach = solvec)) 
InverseLow <- data.frame(cbind(FD = seq(LowerLimit,2,0.001),Mach = solvec_L))
InverseHigh <- data.frame(cbind(FD = seq(UpperLimit,2,0.0001),Mach = solvec_U))

head(InverseHigh)



ggplot(aes(x=log10(Mach),y=FD),data=ModelDatFit_Mach) +
     #geom_errorbarh(aes(xmax=FD + std,xmin=FD -std,col=`Mach Number`),alpha=0.05) +
     geom_ribbon(aes(ymax=FD + std,ymin=FD -std,fill=`Mach Number`),alpha=0.1) +
     geom_point(aes(col=`Mach Number`)) +
     coord_flip() +
     geom_point(aes(x=log10(Mach),y=`Fractal Dimension`),data=ModelDatNotFit_Mach,col='grey')+
     geom_line(aes(x=Mach,y=FD),data=InverseModel,size=1.1) +
     geom_line(aes(x=Mach,y=FD),data=InverseLow,size=0.5,linetype=2) +
     geom_line(aes(x=Mach,y=FD),data=InverseHigh,size=0.5,linetype=2) +
     theme_bw() +
     scale_x_continuous(breaks = seq(-1.5, 2, by = 0.5),labels=ExponentLabel,limit=c(-1.5,2)) +
     scale_y_continuous(breaks = round(seq(1.4,2, by = 0.1),1)) + 
     theme(axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15),
           axis.title = element_text(size = 16)) + 
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     labs(y = 'Fractal Dimension', x = 'Mach Number') +
     theme(legend.text = element_text(size = 10), 
           legend.title = element_text(size = 12,face = "bold"),
                                       legend.background = element_rect(fill = "white"), 
                                       legend.position = c(0.9, 0.80)) 

 head(ModelDatFit)


# Larson's Paper
#####################################################################################################################

setwd('/Volumes/JamesBe/PeriodicBoundaryDatasets/')

Larson <- read_csv('LarsonTableObjects.csv')

Larson %<>% mutate(`Mach Number` = `v (km/s)`/0.2,
                   `Fractal Dim` = nlsFunc3_mod(log10(`Mach Number`)),
                   `Error` = nlsFunc3_upper(log10(`Mach Number`)) - nlsFunc3_mod(log10(`Mach Number`)))

Larson





file<-Mach4SMP_2512 %>% filter(Projection=='x') 
file$Filename %>% unique()


# Paper Figures: M(D)
##########################################################################################################


DatTransform <- function(data,L,MachNumber){
     data %>% mutate(LR = log10(Length/L)) %>%
          group_by(LR) %>%
          summarise(FD = mean(FractalDim),
                    `std` = sd(FractalDim),
                    Length = mean(Length)) %>%
          mutate(`Mach Number` = MachNumber) %>%
          return()
}

Mach4SMP_2512 %<>% filter(TurnoverTime != 2.0)
uniqueTT <- Mach4SMP_2512$TurnoverTime %>% unique()

Mach1SMP_1024       <- filter(Mach1SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach10SMP_1024      <- filter(Mach10SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach20SMP_1024      <- filter(Mach20SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach40SMP_1024      <- filter(Mach40SMP_1024, TurnoverTime %in%  uniqueTT) 
Mach100SMP_1024     <- filter(Mach100SMP_1024, TurnoverTime %in%  uniqueTT) 

M1   <- DatTransform(Mach1SMP_1024,1024,"Mach 1")
M4   <- DatTransform(Mach4SMP_2512,2512,"Mach 4")
M10  <- DatTransform(Mach10SMP_1024,1024,"Mach 10")
M20  <- DatTransform(Mach20SMP_1024,1024,"Mach 20")
M40  <- DatTransform(Mach40SMP_1024,1024,"Mach 40")
M100 <- DatTransform(Mach100SMP_1024,1024,"Mach 100")

M    <- rbind(M1,M4,M10,M20,M40,M100)
M    %<>% mutate(`Mach Number` = factor(`Mach Number`, levels=c("Mach 1","Mach 4","Mach 10", "Mach 20","Mach 40","Mach 100")) )


# Figure 2
ggplot(aes(x=LR,y=FD),data=M) +
     geom_point(aes(col=`Mach Number`)) +
     geom_ribbon(aes(ymax = FD + std,ymin= FD - std,fill=`Mach Number`),alpha=0.2) +
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     facet_wrap(~`Mach Number`,nrow=1) +
     theme_bw() +
     scale_y_continuous(breaks = round(seq(1.4,2, by = 0.1),1)) + 
     scale_x_continuous(breaks = round(seq(-3,0+5, by = 1),1),
                        labels=ExponentLabel) +
     theme(axis.title = element_text(size = 16), 
           axis.text = element_text(size = 8), 
           axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15),
           strip.text.x = element_text(size = 15,face="bold")) + 
     theme(legend.position = 0) 


# Figure 3 part 1
ggplot(aes(x=LR,y=FD),data=ModelDatFit) +
     geom_point(aes(col=`Mach Number`)) +
     geom_ribbon(aes(ymax=FD + std,ymin=FD -std,fill=`Mach Number`),alpha=0.1) +
     geom_point(aes(x=LengthRel,y=`Fractal Dimension`),data=ModelDatNotFit,col='grey')+
     theme_bw() +
     scale_y_continuous(breaks = round(seq(1.4,2, by = 0.1),1)) + 
     scale_x_continuous(breaks = round(seq(min(Mach1Rel_1024$LengthRel), 
                                           max(Mach40Rel_2048$LengthRel)+5, by = 0.5),1)-0.1,
                        labels=ExponentLabel) + 
     theme(axis.title = element_text(size = 16), 
           axis.text = element_text(size = 8), 
           axis.text.x = element_text(size = 12), 
           axis.text.y = element_text(size = 12), 
           plot.title = element_text(size = 15)) + 
     labs(x = expression(paste(l/L)),y = 'Fractal Dimension') + 
     theme(legend.position = 0) +
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     stat_function(fun=nlsFunc2,args=list(parms=parms),col='blue',linetype=3) +
     stat_function(fun=nlsFunc2,args=list(parms=parms),col='black',size=1.1) +
     stat_function(fun=nlsFunc2,args=list(parms=parms_Lower),col='black',size=0.5,linetype=2) +
     stat_function(fun=nlsFunc2,args=list(parms=parms_Upper),col='black',size=0.5,linetype=2) +
     geom_segment(y = parms[3],yend = parms[3],x=1.5,xend=3.5,linetype=3) +
     geom_segment(y = 2,yend = 2,x=-4.5,xend=-2.5,linetype=3) 

# Figure 3 part 2
ggplot(aes(x=log10(Mach),y=FD),data=ModelDatFit_Mach) +
     geom_ribbon(aes(ymax=FD + std,ymin=FD -std,fill=`Mach Number`),alpha=0.1) +
     geom_point(aes(col=`Mach Number`)) +
     coord_flip() +
     geom_point(aes(x=log10(Mach),y=`Fractal Dimension`),data=ModelDatNotFit_Mach,col='grey')+
     geom_line(aes(x=Mach,y=FD),data=InverseModel,size=2) +
     geom_line(aes(x=Mach,y=FD),data=InverseLow,size=0.5,linetype=2) +
     geom_line(aes(x=Mach,y=FD),data=InverseHigh,size=0.5,linetype=2) +
     theme_bw() +
     scale_x_continuous(breaks = seq(-1.5, 2, by = 0.5),labels=ExponentLabel,limit=c(-1.5,2)) +
     scale_y_continuous(breaks = round(seq(1.4,2, by = 0.1),1)) + 
     theme(axis.text.x = element_text(size = 15), 
           axis.text.y = element_text(size = 15), 
           plot.title = element_text(size = 15),
           axis.title = element_text(size = 16)) + 
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     labs(y = 'Fractal Dimension', x = 'Mach Number') +
     theme(legend.text = element_text(size = 16), 
           legend.title = element_text(size = 16,face = "bold"),
           legend.background = element_rect(fill = "white"), 
           legend.position = c(0.9, 0.90)) 

ggplot(aes(x=log10(Mach),y=FD),data=ModelDatFit_Mach) +
     coord_flip() +
     geom_line(aes(x=Mach,y=FD),data=InverseModel,size=1.1) +
     geom_line(aes(x=Mach,y=FD),data=InverseLow,size=0.5,linetype=2) +
     geom_line(aes(x=Mach,y=FD),data=InverseHigh,size=0.5,linetype=2) +
     theme_bw() +
     scale_x_continuous(breaks = seq(-1.5, 2, by = 0.5),labels=ExponentLabel,limit=c(-1.5,2)) +
     scale_y_continuous(breaks = round(seq(1.4,2, by = 0.1),1)) + 
     theme(axis.text.x = element_text(size = 15), 
           axis.text.y = element_text(size = 15), 
           plot.title = element_text(size = 15),
           axis.title = element_text(size = 16)) + 
     scale_color_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     scale_fill_manual(values=c("#2121D9", "#9999FF", "#21D921", "#D92121", "#FF9326", "#CCCC00")) +
     labs(y = 'Fractal Dimension', x = 'Mach Number') +
     theme(legend.text = element_text(size = 16), 
           legend.title = element_text(size = 16,face = "bold"),
           legend.background = element_rect(fill = "white"), 
           legend.position = c(0.9, 0.90)) 
