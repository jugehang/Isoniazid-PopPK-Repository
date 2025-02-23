## Model Repository of Isoniazid
## Author info.
#### Author: Graham Ju
#### Address: Central South University, Xiangya Hospital, Hunan Province, China
#### Email: jugehang@163.com
#### Date: 2023/11/21

## Code
#**START**#------------Preparatory work--------------##----
#setwd and library R packages
rm(list=ls())
#set working directory to current folder
curr.dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)

# load R packages
library(tidyverse) # for data visualisation and manipulation
library(mrgsolve) # for simulation
library(rxode2) # for simulation for complex model
library(cowplot) # for combine figures
library(ggplot2) 
library(ggpubr)
library(dplyr)
library(PopED)
library(ggsci)
library(gridExtra)

#create fold for figure output
output_dir <- "INH0001"
if (!file.exists(output_dir)) {dir.create (output_dir)}
#**END**#------------Preparatory work--------------#----
##------------First Part: Adult Model-----------#------------------
#**START**#------Definition of Typical Patients------
#### Body Weight: 70kg
#### Fat-free Mass: 44kg
#### NAT2 phenotype: RA, IA, SA (IA/RA: FA)
#### Other Covariates effects set to 0
#### Administration Methods: oral 300mg qd
#**END**#------Definition of Typical Patients------
#**Model 1**#----Abdelgawad, N. 2022----#------
#**Model 1**# Model Characteristics----
## Model Structure: 2 CMT with FO elimination and transit CMT absorption
## Model Parameters
#### CL (L/h) = NAT2*(BW/66)^0.75; NAT2_SA = 9.76, NAT2_FA = 25.5
#### Vc (L)   = 59.0*(BW/66)
#### Q  (L/h) = 1.43*(BW/66)^0.75
#### Vp (L)   = 30.7*(BW/66)
#### Ka (h-1) = 2.43
#### MTT (h)  = 0.442
#### NN       = 8 FIX
#### BIO      = 1 FIX
#### BSV (CV%) : CL  = 25.3%
#### BSV (CV%) : MTT = 122%
#### BSV (CV%) : BIO = 34.9%
#### Prop.err = 13.9%
#### Add.err  = 0.021 FIX
#**Model 1**# MrgSolve----
set.seed(112101)

cod1 <- '
$PARAM CL=9.76, VC=59.0, Q=1.43, VP=30.7, KA=2.43, MTT=0.442, NN=8, BIO=1,
$CMT DEPOT CENT PERI
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[300];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC Q VP KA MTT NN BIO
'

mod1 <- mcode("transit", cod1, atol=1e-8, rtol=1e-8, maxsteps=50000)
#**Model 1**# PopED----
#**Model 1**# -- Model structure definition function----
ff1 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   Q  = p[["Q"]], 
                   VP = p[["VP"]], 
                   KA = p[["KA"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod1, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 1**# -- parameter definition function----
fg1 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/66)^(0.75) * exp(b[1]),
    VC = bpop[2] * (a[3]/66),
    Q  = bpop[3] * (a[3]/66)^(0.75),
    VP = bpop[4] * (a[3]/66),
    KA = bpop[5],
    BIO = bpop[6] * exp(b[2]),
    MTT = bpop[7] * exp(b[3]),
    NN = bpop[8],
    DOSE = a[1],
    TAU = a[2],
    BW  = a[3]
  )
  return(parameters)
}

#**Model 1**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 1**#**SA**# Create PopED databases----
poped_db_model1 <- create.poped.database(
  ff_fun =ff1,
  fg_fun =fg1,
  fError_fun=feps1,
  bpop=c(CL=9.76,
         VC=59.0,
         Q =1.43,
         VP=30.7,
         KA=2.43,
         BIO=1,
         MTT=0.442,
         NN=8), 
  notfixed_bpop = c(1,0,0,0,0,1,1,0),
  d=c(CL=0.253^2,BIO=0.349^2,MTT=1.22^2),
  sigma=c(prop=0.139, add=0.021),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl1 <- plot_model_prediction(poped_db_model1,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model1: Abdelgawad N.(2022) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.9), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl1

jpeg(filename = paste0(output_dir,"/1.Abdelgawad2022_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl1)
dev.off()
#**Model 1**#**FA**# Create PopED databases----
poped_db_model1.1 <- create.poped.database(
  ff_fun =ff1,
  fg_fun =fg1,
  fError_fun=feps1,
  bpop=c(CL=25.5,
         VC=59.0,
         Q =1.43,
         VP=30.7,
         KA=2.43,
         BIO=1,
         MTT=0.442,
         NN=8), 
  notfixed_bpop = c(1,0,0,0,0,1,1,0),
  d=c(CL=0.253^2,BIO=0.349^2,MTT=1.22^2),
  sigma=c(prop=0.139, add=0.021),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl101 <- plot_model_prediction(poped_db_model1.1,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model1: Abdelgawad N.(2022) #NAT2_FA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=2.8), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl101

jpeg(filename = paste0(output_dir,"/1.Abdelgawad2022_Adult_FA.jpg"), 
     width=1600, height=1600, res=300)
print(pl101)
dev.off()

#**Model 2**#----Abdelwahab, M. T. 2020 #------
#**Model 2**# Model Characteristics----
## Model Structure: 2 CMT with FO elimination and transit CMT absorption
## Model Parameters
#### CL (L/h) = NAT2*(BW/65.3)^0.75; SA = 29.0, IA = 75.7, RA = 97.1
#### Vc (L)   = 130*(BW/65.3)
#### Q  (L/h) = 12.4*(BW/65.3)^0.75
#### Vp (L)   = 28.5*(BW/65.3)
#### Ka (h-1) = KTR
#### MTT (h)  = 1.21
#### NN       = 8.01
#### BIO      = 1 FIX
#### BSV (CV%) : CL  = 12.7%
#### BOV (CV%) : MTT = 56.7%
#### BOV (CV%) : BIO = 36.7%
#### Prop.err = 22.2%
#### Add.err  = 0.045
#**Model 2**# MrgSolve----
set.seed(112101)

cod2 <- '
$PARAM CL=29, VC=130, Q=12.4, VP=28.5, MTT=1.21, NN=8.01, BIO=1,
$CMT DEPOT CENT PERI
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[300];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KTR*DEPOT;
dxdt_CENT = KTR*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC Q VP MTT NN BIO
'

mod2 <- mcode("transit", cod2, atol=1e-8, rtol=1e-8, maxsteps=50000)
#**Model 2**# PopED----
#**Model 2**# -- Model structure definition function----
ff2 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   Q  = p[["Q"]], 
                   VP = p[["VP"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod2, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 2**# -- parameter definition function----
fg2 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/66)^(0.75) * exp(b[1]),
    VC = bpop[2] * (a[3]/66),
    Q  = bpop[3] * (a[3]/66)^(0.75),
    VP = bpop[4] * (a[3]/66),
    BIO = bpop[5] * exp(b[2]),
    MTT = bpop[6] * exp(b[3]),
    NN = bpop[7],
    DOSE = a[1],
    TAU = a[2],
    BW  = a[3]
  )
  return(parameters)
}

#**Model 2**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 2**#**SA**# Create PopED databases----
poped_db_model2 <- create.poped.database(
  ff_fun =ff2,
  fg_fun =fg2,
  fError_fun=feps1,
  bpop=c(CL=29,
         VC=130.0,
         Q =12.4,
         VP=28.5,
         BIO=1,
         MTT=1.21,
         NN=8.01), 
  notfixed_bpop = c(1,0,0,0,1,1,0),
  d=c(CL=0.127^2,BIO=0,MTT=0),
  sigma=c(prop=0.222, add=0.045),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl2 <- plot_model_prediction(poped_db_model2,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model2: Abdelwahab, M. T. (2020) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=1.65), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl2

jpeg(filename = paste0(output_dir,"/2.Abdelwahab2020_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl2)
dev.off()
#**Model 2**#**IA**# Create PopED databases----
poped_db_model201 <- create.poped.database(
  ff_fun =ff2,
  fg_fun =fg2,
  fError_fun=feps1,
  bpop=c(CL=75.7,
         VC=130.0,
         Q =12.4,
         VP=28.5,
         BIO=1,
         MTT=1.21,
         NN=8.01), 
  notfixed_bpop = c(1,0,0,0,1,1,0),
  d=c(CL=0.127^2,BIO=0,MTT=0),
  sigma=c(prop=0.222, add=0.045),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl201 <- plot_model_prediction(poped_db_model201,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model2: Abdelwahab, M. T. (2020) #NAT2_IA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=1.27), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl201

jpeg(filename = paste0(output_dir,"/2.Abdelwahab2020_Adult_IA.jpg"), 
     width=1600, height=1600, res=300)
print(pl201)
dev.off()

#**Model 2**#**RA**# Create PopED databases----
poped_db_model202 <- create.poped.database(
  ff_fun =ff2,
  fg_fun =fg2,
  fError_fun=feps1,
  bpop=c(CL=97.1,
         VC=130.0,
         Q =12.4,
         VP=28.5,
         BIO=1,
         MTT=1.21,
         NN=8.01), 
  notfixed_bpop = c(1,0,0,0,1,1,0),
  d=c(CL=0.127^2,BIO=0,MTT=0),
  sigma=c(prop=0.222, add=0.045),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl202 <- plot_model_prediction(poped_db_model202,
                               PI=F,
                               sample.times = T,
                               PRED = T) +
  labs(title = "Model2: Abdelwahab, M. T. (2020) #NAT2_RA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=1.13), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl202

jpeg(filename = paste0(output_dir,"/2.Abdelwahab2020_Adult_RA.jpg"), 
     width=1600, height=1600, res=300)
print(pl202)
dev.off()

#**Model 3**#----Chen B. 2022----#------
#**Model 3**# Model Characteristics----
## Model Structure: 1 CMT with FO absorption and FO elimination
## Model Parameters
#### CL (L/h) = 28.7*e^(-0.55*NAT2) ; SA=2, IA=1, RA=0
#### Vc (L)   = 54.1
#### Ka (h-1) = 3.91
#### BSV (CV%) : CL  = 30.7%
#### BSV (CV%) : Vc  = 19.4%
#### BSV (CV%) : Ka  = 55.2%
#### Prop.err = 33.3%

#**Model 3**# MrgSolve----
set.seed(112101)
cod3 <- '
$PARAM CL=28.7, VC=54.1, KA=3.91
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod3<- mcode("optim", cod3, atol=1e-8, rtol=1e-8, maxsteps=50000)

#**Model 3**# PopED----
#**Model 3**# -- Model structure definition function----
ff3 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod3, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 3**# -- parameter definition function----
fg3 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(-0.55*a[3]) * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    KA = bpop[3] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    NAT = a[3]
  )
  return(parameters)
}

#**Model 3**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 3**#**SA**# Create PopED databases----
poped_db_model3 <- create.poped.database(
  ff_fun =ff3,
  fg_fun =fg3,
  fError_fun=feps1,
  bpop=c(CL=28.7,
         VC=54.1,
         KA=3.91), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.307^2,VC=0.194^2,KA=0.552^2),
  sigma=c(prop=0.333, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,NAT2=2))

pl3 <- plot_model_prediction(poped_db_model3,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model3: Chen (2022) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=4.12), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl3

jpeg(filename = paste0(output_dir,"/3.Chen2022_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl3)
dev.off()
#**Model 3**#**IA**# Create PopED databases----
poped_db_model301 <- create.poped.database(
  ff_fun =ff3,
  fg_fun =fg3,
  fError_fun=feps1,
  bpop=c(CL=28.7,
         VC=54.1,
         KA=3.91), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.307^2,VC=0.194^2,KA=0.552^2),
  sigma=c(prop=0.333, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,NAT2=1))

pl301 <- plot_model_prediction(poped_db_model301,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model3: Chen (2022) #NAT2_IA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.25), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl301

jpeg(filename = paste0(output_dir,"/3.Chen2022_Adult_IA.jpg"), 
     width=1600, height=1600, res=300)
print(pl301)
dev.off()
#**Model 3**#**RA**# Create PopED databases----
poped_db_model302 <- create.poped.database(
  ff_fun =ff3,
  fg_fun =fg3,
  fError_fun=feps1,
  bpop=c(CL=28.7,
         VC=54.1,
         KA=3.91), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.307^2,VC=0.194^2,KA=0.552^2),
  sigma=c(prop=0.333, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,NAT2=0))

pl302 <- plot_model_prediction(poped_db_model302,
                               PI=F,
                               sample.times = T,
                               PRED = T) +
  labs(title = "Model3: Chen (2022) #NAT2_RA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=2.2), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl302

jpeg(filename = paste0(output_dir,"/3.Chen2022_Adult_RA.jpg"), 
     width=1600, height=1600, res=300)
print(pl302)
dev.off()

#**Model 4**#----Chirehwa 2019----#------
#**Model 4**# Model Characteristics----
## Model Structure: 2CMT with FO absorption with lag time and FO elimination
## Model Parameters
#### CL (L/h) = NAT2*(FFM/43.3)^0.75*(1+0.541)^EFV ; SA=6.64, IA=11.4, RA=35.9, EFV=0
#### Vc (L)   = 47.7*(FFM/43.3)
#### Q  (L/h) = 4.8*(FFM/43.3)^0.75
#### Vp (L)   = 8.14*(FFM/43.3)
#### Ka (h-1) = 1.59
#### Tlag (h) = 0.289
#### BSV (CV%) : CL  = 35.6%
#### BSV (CV%) : Vc  = 34.6%
#### BSV (CV%) : Q   = 18.9%
#### BSV (CV%) : Ka  = 36.7%
#### Prop.err = 11.2%
#### Add.err  = 0.02 FIX

#**Model 4**# MrgSolve----
set.seed(112101)
cod4 <- '
$PARAM CL=6.64, VC=47.7, KA=1.59, Q=4.8, VP=8.14, TLAG=0.289
$CMT DEPOT CENT PERI
$MAIN
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA Q VP
'
mod4 <- mcode("optim", cod4, atol=1e-8, rtol=1e-8, maxsteps=50000)

#**Model 4**# PopED----
#**Model 4**# -- Model structure definition function----
ff4 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]],
                   Q  = p[["Q"]], 
                   VP = p[["VP"]], 
                   TLAG=p[["TLAG"]])
  out <- mrgsim_q(mod4, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 4**# -- parameter definition function----
fg4 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/43.3)^(0.75) * exp(b[1]),
    VC = bpop[2] * (a[3]/43.3) * exp(b[2]),
    Q  = bpop[3] * (a[3]/43.3)^(0.75) * exp(b[3]),
    VP = bpop[4] * (a[3]/43.3),
    KA = bpop[5] * exp(b[2]) * exp(b[4]),
    TLAG=bpop[6],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3]
  )
  return(parameters)
}

#**Model 4**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 4**#**SA**# Create PopED databases----
poped_db_model4 <- create.poped.database(
  ff_fun =ff4,
  fg_fun =fg4,
  fError_fun=feps1,
  bpop=c(CL=6.64,
         VC=47.7,
         Q=4.8,
         VP=8.14,
         KA=1.59,
         TLAG=0.289), 
  notfixed_bpop = c(1,1,1,1,0,0),
  d=c(CL=0.356^2,VC=0.346^2,Q=0.189^2,KA=0.367^2),
  sigma=c(prop=0.112, add=0.02),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44))

pl4 <- plot_model_prediction(poped_db_model4,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model4: Chirehwa (2019) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=4.83), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl4

jpeg(filename = paste0(output_dir,"/4.Chirehwa2019_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl4)
dev.off()
#**Model 4**#**IA**# Create PopED databases----
poped_db_model401 <- create.poped.database(
  ff_fun =ff4,
  fg_fun =fg4,
  fError_fun=feps1,
  bpop=c(CL=11.4,
         VC=47.7,
         Q=4.8,
         VP=8.14,
         KA=1.59,
         TLAG=0.289), 
  notfixed_bpop = c(1,1,1,1,0,0),
  d=c(CL=0.356^2,VC=0.346^2,Q=0.189^2,KA=0.367^2),
  sigma=c(prop=0.112, add=0.02),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44))

pl401 <- plot_model_prediction(poped_db_model401,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model4: Chirehwa (2019) #NAT2_IA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=4.1), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl401

jpeg(filename = paste0(output_dir,"/4.Chirehwa2019_Adult_IA.jpg"), 
     width=1600, height=1600, res=300)
print(pl401)
dev.off()
#**Model 4**#**RA**# Create PopED databases----
poped_db_model402 <- create.poped.database(
  ff_fun =ff4,
  fg_fun =fg4,
  fError_fun=feps1,
  bpop=c(CL=35.9,
         VC=47.7,
         Q=4.8,
         VP=8.14,
         KA=1.59,
         TLAG=0.289), 
  notfixed_bpop = c(1,1,1,1,0,0),
  d=c(CL=0.356^2,VC=0.346^2,Q=0.189^2,KA=0.367^2),
  sigma=c(prop=0.112, add=0.02),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44))

pl402 <- plot_model_prediction(poped_db_model402,
                               PI=F,
                               sample.times = T,
                               PRED = T) +
  labs(title = "Model4: Chirehwa (2019) #NAT2_RA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=2.35), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl402

jpeg(filename = paste0(output_dir,"/4.Chirehwa2019_Adult_RA.jpg"), 
     width=1600, height=1600, res=300)
print(pl402)
dev.off()

#**Model 5**#----Cho.Y.S 2021----#------
#**Model 5**# Model Characteristics----
## Model Structure: 2CMT with absorption lag time and sequential ZO (D0) and FO absorption with FO elimination
## Model Parameters
#### CL (L/h) = 22.2*(1-NAT2)*(FFM/50)^0.75 ; SA=0.646, IA=0.274, RA=0
#### Vc (L)   = 16.5*(FFM/50)
#### Q  (L/h) = 18.4
#### Vp (L)   = 36.4
#### Ka (h-1) = 1.21
#### Tlag (h) = 0.02 FIX
#### D0 (h)   = 0.47
#### BSV (CV%) : CL  = 14%
#### BSV (CV%) : Vc  = 3%  FIX
#### BSV (CV%) : VP  = 15% FIX
#### BSV (CV%) : Ka  = 60% FIX
#### BSV (CV%) : Tlag= 22% FIX
#### BSV (CV%) : D0  = 20% FIX
#### Prop.err = 29.2%
#### Add.err  = 0.134

#**Model 5**# MrgSolve----
set.seed(112101)
cod5 <- '
$PARAM CL=22.2, VC=16.5, KA=1.21, Q=18.4, VP=36.4, TLAG=0.289, DUR=0.47
$CMT DEPOT CENT PERI
$MAIN
D_CENT = DUR;
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA Q VP
'
mod5 <- mcode("optim", cod5, atol=1e-8, rtol=1e-8, maxsteps=50000)

#**Model 5**# PopED----
#**Model 5**# -- Model structure definition function----
ff5 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]],
                   Q  = p[["Q"]], 
                   VP = p[["VP"]], 
                   TLAG=p[["TLAG"]],
                   DUR= p[["DUR"]])
  out <- mrgsim_q(mod5, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 5**# -- parameter definition function----
fg5 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/50)^(0.75) * (1-a[4]) * exp(b[1]),
    VC = bpop[2] * (a[3]/50) * exp(b[2]),
    Q  = bpop[3],
    VP = bpop[4] * exp(b[3]),
    KA = bpop[5] * exp(b[2]) * exp(b[4]),
    TLAG=bpop[6] * exp(b[5]),
    DUR= bpop[7] * exp(b[6]),
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    NAT2= a[4]
  )
  return(parameters)
}

#**Model 5**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 5**#**SA**# Create PopED databases----
poped_db_model5 <- create.poped.database(
  ff_fun =ff5,
  fg_fun =fg5,
  fError_fun=feps1,
  bpop=c(CL=22.2,
         VC=16.5,
         Q=18.4,
         VP=36.4,
         KA=1.21,
         TLAG=0.02,
         DUR=0.47), 
  notfixed_bpop = c(1,1,1,1,1,1,1),
  d=c(CL=0.14^2,VC=0.03^2,VP=0.15^2,KA=0.6^2,TLAG=0.22^2,DUR=0.2^2),
  sigma=c(prop=0.112, add=0.02),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44,NAT2=0.646))

pl5 <- plot_model_prediction(poped_db_model5,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model5: Cho Y.S. (2021) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=5.0), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl5

jpeg(filename = paste0(output_dir,"/5.ChoYS2021_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl5)
dev.off()
#**Model 5**#**IA**# Create PopED databases----
poped_db_model501 <- create.poped.database(
  ff_fun =ff5,
  fg_fun =fg5,
  fError_fun=feps1,
  bpop=c(CL=22.2,
         VC=16.5,
         Q=18.4,
         VP=36.4,
         KA=1.21,
         TLAG=0.02,
         DUR=0.47), 
  notfixed_bpop = c(1,1,1,1,1,1,1),
  d=c(CL=0.14^2,VC=0.03^2,VP=0.15^2,KA=0.6^2,TLAG=0.22^2,DUR=0.2^2),
  sigma=c(prop=0.112, add=0.02),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44,NAT2=0.274))

pl501 <- plot_model_prediction(poped_db_model501,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model5: Cho Y.S. (2021) #NAT2_IA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.25), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl501

jpeg(filename = paste0(output_dir,"/5.ChoYS2021_Adult_IA.jpg"), 
     width=1600, height=1600, res=300)
print(pl501)
dev.off()
#**Model 5**#**RA**# Create PopED databases----
poped_db_model502 <- create.poped.database(
  ff_fun =ff5,
  fg_fun =fg5,
  fError_fun=feps1,
  bpop=c(CL=22.2,
         VC=16.5,
         Q=18.4,
         VP=36.4,
         KA=1.21,
         TLAG=0.02,
         DUR=0.47), 
  notfixed_bpop = c(1,1,1,1,1,1,1),
  d=c(CL=0.14^2,VC=0.03^2,VP=0.15^2,KA=0.6^2,TLAG=0.22^2,DUR=0.2^2),
  sigma=c(prop=0.112, add=0.02),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44,NAT2=0))

pl502 <- plot_model_prediction(poped_db_model502,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model5: Cho Y.S. (2021) #NAT2_RA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=2.5), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl502

jpeg(filename = paste0(output_dir,"/5.ChoYS2021_Adult_RA.jpg"), 
     width=1600, height=1600, res=300)
print(pl502)
dev.off()

#**Model 6**#----Denti, P. 2015----#------
#**Model 6**# Model Characteristics----
## Model Structure: 2 CMT with  transit comp absorption
## Model Parameters
#### CL (L/h) = NAT2*(FFM/43.3)^0.75; SA = 15.1, FA = 26.1
#### Vc (L)   = 48.2*(FFM/43)
#### Q  (L/h) = 16.1*(FFM/43.3)^0.75
#### Vp (L)   = 16.5*(FFM/50)
#### MTT (h)  = 0.924
#### NN       = 2.73
#### BIO      = 1 FIX
#### BSV (CV%) : CL  = 30.7%
#### BOV (CV%) : MTT = 37.4%
#### BOV (CV%) : BIO = 12.8%
#### Prop.err = 13.3%
#### Add.err  = 0.0224
#**Model 6**# MrgSolve----
set.seed(112101)

cod6 <- '
$PARAM CL=15.1, VC=48.2, Q=16.1, VP=16.5, MTT=0.924, NN=2.73, BIO=1,
$CMT DEPOT CENT PERI
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[300];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KTR*DEPOT;
dxdt_CENT = KTR*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC Q VP MTT NN BIO
'

mod6 <- mcode("transit", cod6, atol=1e-8, rtol=1e-8, maxsteps=50000)
#**Model 6**# PopED----
#**Model 6**# -- Model structure definition function----
ff6 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   Q  = p[["Q"]], 
                   VP = p[["VP"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod6, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 6**# -- parameter definition function----
fg6 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/43)^(0.75) * exp(b[1]),
    VC = bpop[2] * (a[3]/43),
    Q  = bpop[3] * (a[3]/43)^(0.75),
    VP = bpop[4] * (a[3]/43),
    BIO = bpop[5] * exp(b[2]),
    MTT = bpop[6] * exp(b[3]),
    NN = bpop[7],
    DOSE = a[1],
    TAU = a[2],
    FMM = a[3]
  )
  return(parameters)
}

#**Model 6**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 6**#**SA**# Create PopED databases----
poped_db_model6 <- create.poped.database(
  ff_fun =ff6,
  fg_fun =fg6,
  fError_fun=feps1,
  bpop=c(CL=15.1,
         VC=48.2,
         Q =16.1,
         VP=16.5,
         BIO=1,
         MTT=0.924,
         NN=2.73), 
  notfixed_bpop = c(1,1,1,1,0,0,0,0),
  d=c(CL=0.307^2,BIO=0,MTT=0),
  sigma=c(prop=0.133, add=0.0224),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44))

pl6 <- plot_model_prediction(poped_db_model6,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model6: Denti P. (2015) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.6), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl6

jpeg(filename = paste0(output_dir,"/6.DentiP2015_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl6)
dev.off()
#**Model 6**#**FA**# Create PopED databases----
poped_db_model601 <- create.poped.database(
  ff_fun =ff6,
  fg_fun =fg6,
  fError_fun=feps1,
  bpop=c(CL=26.1,
         VC=48.2,
         Q =16.1,
         VP=16.5,
         BIO=1,
         MTT=0.924,
         NN=2.73), 
  notfixed_bpop = c(1,1,1,1,0,0,0,0),
  d=c(CL=0.307^2,BIO=0,MTT=0),
  sigma=c(prop=0.133, add=0.0224),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,FFM=44))

pl601 <- plot_model_prediction(poped_db_model601,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model6: Denti P. (2015) #NAT2_FA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.0), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl601

jpeg(filename = paste0(output_dir,"/6.DentiP2015_Adult_FA.jpg"), 
     width=1600, height=1600, res=300)
print(pl601)
dev.off()

#**Model 7**#----Gao Y. 2021----#------
#**Model 7**# Model Characteristics----
## Model Structure: 1 CMT with FO absorption and FO elimination
## Model Parameters
#### CL (L/h) = NAT2*(BW/50)^0.55 ; SA=12.6, IA=16.0, RA=30.6
#### Vc (L)   = 21.2
#### Q (L/h)  = 8.7
#### Vp (L)   = 125.8
#### Ka (h-1) = 0.68
#### BSV (CV%) : CL  = 60.9%
#### BSV (CV%) : Vc  = 21.7%
#### BSV (CV%) : Ka  = 23.6%
#### Add.err  = 0.178

#**Model 7**# MrgSolve----
set.seed(112101)
cod7 <- '
$PARAM CL=12.6, VC=21.2, KA=0.68, Q=8.7, VP=125.8
$CMT DEPOT CENT PERI
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA Q VP
'
mod7 <- mcode("optim", cod7, atol=1e-8, rtol=1e-8, maxsteps=50000)

#**Model 7**# PopED----
#**Model 7**# -- Model structure definition function----
ff7 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   Q  = p[["Q"]], 
                   VP = p[["VP"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod7, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

#**Model 7**# -- parameter definition function----
fg7 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/50)^0.55 * exp(b[1]),
    VC = bpop[2] * exp(b[2]) * exp(b[2]),
    KA = bpop[3] * exp(b[3]) * exp(b[3]),
    Q  = bpop[4],
    VP = bpop[5],
    DOSE = a[1],
    TAU = a[2],
    BW  = a[3]
  )
  return(parameters)
}

#**Model 7**# -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps7 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

#**Model 7**#**SA**# Create PopED databases----
poped_db_model7 <- create.poped.database(
  ff_fun =ff7,
  fg_fun =fg7,
  fError_fun=feps1,
  bpop=c(CL=12.6,
         VC=21.2,
         KA=0.68,
         Q =8.7,
         VP=125.8), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.609^2,VC=0.217^2,KA=0.236^2),
  sigma=c(prop=0.333, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl7 <- plot_model_prediction(poped_db_model7,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model7: Gao Y. (2021) #NAT2_SA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.56), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl7

jpeg(filename = paste0(output_dir,"/7.GaoY2021_Adult_SA.jpg"), 
     width=1600, height=1600, res=300)
print(pl7)
dev.off()
#**Model 7**#**IA**# Create PopED databases----
poped_db_model701 <- create.poped.database(
  ff_fun =ff7,
  fg_fun =fg7,
  fError_fun=feps1,
  bpop=c(CL=16.0,
         VC=21.2,
         KA=0.68,
         Q =8.7,
         VP=125.8), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.609^2,VC=0.217^2,KA=0.236^2),
  sigma=c(prop=0.333, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl701 <- plot_model_prediction(poped_db_model701,
                             PI=F,
                             sample.times = T,
                             PRED = T) +
  labs(title = "Model7: Gao Y. (2021) #NAT2_IA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=3.0), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl701

jpeg(filename = paste0(output_dir,"/7.GaoY2021_Adult_IA.jpg"), 
     width=1600, height=1600, res=300)
print(pl701)
dev.off()
#**Model 7**#**RA**# Create PopED databases----
poped_db_model702 <- create.poped.database(
  ff_fun =ff7,
  fg_fun =fg7,
  fError_fun=feps1,
  bpop=c(CL=30.6,
         VC=21.2,
         KA=0.68,
         Q =8.7,
         VP=125.8), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.609^2,VC=0.217^2,KA=0.236^2),
  sigma=c(prop=0.333, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=300,TAU=24,BW=70))

pl702 <- plot_model_prediction(poped_db_model702,
                               PI=F,
                               sample.times = T,
                               PRED = T) +
  labs(title = "Model7: Gao Y. (2021) #NAT2_RA",
       x = " ",
       y = " ") +
  annotate("rect", xmin = 168, xmax = 192,  ymin = 3, ymax = 6, alpha = 0.2, fill = "gray")+
  geom_line(aes(color="Predicted"), size=1.5)+
  scale_color_manual(values = "#3C5488FF")+
  scale_y_continuous(limits = c(0,7))+
  geom_point(aes(x=170, y=1.7), size=5, colour="darkred")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.title = element_text(size = 10, face = "bold.italic"))
pl702

jpeg(filename = paste0(output_dir,"/7.GaoY2021_Adult_RA.jpg"), 
     width=1600, height=1600, res=300)
print(pl702)
dev.off()




