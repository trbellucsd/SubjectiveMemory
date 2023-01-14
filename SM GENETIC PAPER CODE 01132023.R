
#library(devtools)
library(umx)
library (OpenMx)
#require(OpenMx)   #Loads OpenMx
library (haven)
library (dplyr)
library(car)
library(psych) 
#require(psych)   #Loads Psych package
#source("GenEpiHelperFunctions.R")
#source("polychoricMeansMatrix3.R")
library(dplyr)
library(haven)
library(lme4)
library(lmerTest)
library(readxl)
library(psych)
source("C:/Users/tyler/Dropbox/Lab Drive/SM and OM Twin Analyses/Twin Analysis/Older Files/GenEpiHelperFunctions.R")
source("C:/Users/tyler/Dropbox/Lab Drive/SM and OM Twin Analyses/Twin Analysis/Older Files/polychoricMeansMatrix3.R")
umx::umx_set_optimizer(opt="SLSQP")
mxOption(NULL,"mvnRelEps",0.0055)
mxOption(NULL, 'Number of Threads', parallel::detectCores())

# PULL FROM MASTER DATA. 
master<-read_sav("L:/VETSA Data/CurrentData/VETSA3_master_with_admin_RSI V3 05072021.sav")
names(master)
names(master)<- toupper(names(master))
selVars0<-c("VETSAID", "CASE", "TWIN", "ZYG2019", "RMCI_CONS_V1", "RMCI_CONS_V2",  "RMCI_CONS_V3", "NAS201TRAN","WNW", "CHARLSON_V1", "CHARLSON_V2", "CHARLSON_V3",
            "RATE_V1", "CONCERN_V1","RATE_V2", "CONCERN_V2", "RATE_V3", "CONCERN_V3", "Q13ISOHR",
            "SM_38", "SM_AVG","SM_V1", "SM_V2","SM_V3", "MEMORY_V1", "MEMORY_V2", "MEMORY_V3", "DEP_V1",
            "DEP_V2", "DEP_V3","AGE_V1", "AGE_V2", "AGE_V3", "NMPQNEM", "NSREACT", "STAITANXTOT_V2", "STAITANXTOT_V3", "STAISANXTOT_V2", 
            "STAISANXTOT_V3", "APOE4VSNON", "ADPGS_05", "ADPGS_1", "ADPGS_2", "ADPGS_3", "ADPGS_4", "ADPGS_5", 
            "KUNKLES1", "KUNKLES1_RES", "ADPGS_NOAPOE_05", "ADPGS_NOAPOE_1", "ADPGS_NOAPOE_2", "ADPGS_NOAPOE_3", "ADPGS_NOAPOE_4", "ADPGS_NOAPOE_5", 
            "NPGRS_S1", "NPGRS_S2", "NPGRS_S3", "NPGRS_S4", "NPGRS_S5", "NPGRS_S6", 
            "NPGRS_S7", "NPGRS_S8", "NPGRS_S2PT5", "DEPPGRS_S1", "DEPPGRS_S2", "DEPPGRS_S2.5",
            "DEPPGRS_S3", "DEPPGRS_S4", "DEPPGRS_S5", "DEPPGRS_S6", "DEPPGRS_S7", "DEPPGRS_S8", "PC1",
            "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", 
            "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "NPGRS_S3_RES",
            "DEPPGRS_S3_RES") 
SMdata<-master[selVars0]

#load("C:/Users/tyler/Downloads/Tyler08092022.RData")
#SMdata<-twins

#write.csv(SMdata, "L:SMdata.csv", row.names = FALSE)
#save(SMdata, file = "L:SMdata.RData")
#write.csv(twins, "L:OrigSMdata.csv", row.names = FALSE)
#describe(SMdata)
twins<-SMdata
#twins<-subset(SMdata, SMdata$RMCI_CONS_V1 != 1 | SMdata$RMCI_CONS_V1 !=2 | SMdata$RMCI_CONS_V1  != 3 | SMdata$RMCI_CONS_V1  != 4)
#load("C:/Users/tyler/Downloads/Tyler08092022.RData")
#names(newtwins)
##Only include data points before MCI dx 
#twins$SM_V1 = ifelse(twins$RMCI_CONS_V1>0,NA,twins$SM_V1)
#twins$SM_V2 = ifelse(twins$RMCI_CONS_V2>0,NA,twins$SM_V2)
#twins$SM_V3 = ifelse(twins$RMCI_CONS_V3>0,NA,twins$SM_V3)
#Only include data points with no baseline MCI 
#twins$SM_V1 = ifelse(twins$RMCI_CONS_V1>0,NA,twins$SM_V1)
#twins$SM_V2 = ifelse(twins$RMCI_CONS_V1>0,NA,twins$SM_V2)
#twins$SM_V3 = ifelse(twins$RMCI_CONS_V1>0,NA,twins$SM_V3)
#twins$SM_38 = scale(twins$SM_38)
#twins$SM_V1 = scale(twins$SM_V1)
#twins$SM_V2 = scale(twins$SM_V2)
#twins$SM_V3 = scale(twins$SM_V3)

#Use CSV file for MPlus, convert NAs to -99, remove header row. 
#Residualize variables 
#Model<-lm(y~ CASE)
#X<-Model$resid 
#View(twins)
names(twins)
names(twins)<- toupper(names(twins))
View(twins)


## Creating MZ and DZ data sets ##
twinA <- twins[twins$TWIN=="1",]
twinB <- twins[twins$TWIN=="2",]

df2 <- merge(twinA, twinB, by=c("CASE","ZYG2019"),all.x=TRUE, all.y=TRUE,suffixes=c("_T1","_T2"))
names(df2)

newtwins <- merge(twinA, twinB, by=c("CASE","ZYG2019"),all.x=TRUE, all.y=TRUE,suffixes=c("_T1","_T2"))
names(newtwins)
#save.image(file = "SM GENETIC PAPER MAIN MODELS.RData")
#load("/Volumes/ngillespie/Documents/work/projects/2016/2016. 8. VETSA/analyses/1. Tyler/subjective memory/SM GENETIC PAPER MAIN MODELS.RData")

df2<-newtwins
table(df2$SM_V1_T1)
# Convert twin file to long format

# Correlations (phenotypic)
 # Method 1 
 # Method 2 FIML

# Univariates
# SM_38
# SM_V1
# SM_V2
# SM_V3

# MULTIVARIATE ANALYSES
# Correlated factors
# Common pathway
# Independent Pathways
# Hybrid correlated factors & autoregression
# Autoregression
# Compare multivariate models

# MULTIVARIATE ANALYSES - VETSA waves only
 # Correlated factors - VETSA waves only
 # Common pathway - VETSA data only
 # Independent Pathways - VETSA data only
 # Autoregression - VETSA only 
 # Compare multivariate models - VETSA data only 



# Correlations (phenotypic)
source("/Volumes/ngillespie/Documents/work/scripts/OpenMx/correlations/polychoricMatrix3b.R")
source("/Volumes/ngillespie/Documents/work/scripts/OpenMx/correlations/polychoricMeansMatrix3.R")


# Convert twin file to long format
 df2_long = umx_wide2long(df2)

 head(df2_long)
 vars = c("SM_38","SM_V1","SM_V2","SM_V3")
 describe(df2_long[vars])
 #       vars    n mean sd median trimmed  mad   min  max range skew kurtosis   se
 # SM_38    1 1555    0  1  -0.36   -0.12 0.85 -0.94 3.59  4.52 0.85     0.01 0.03
 # SM_V1    2  520    0  1  -0.05   -0.05 0.99 -1.59 3.02  4.61 0.59     0.68 0.04
 # SM_V2    3 1199    0  1   0.16   -0.09 0.96 -1.34 3.14  4.48 0.60     0.27 0.03
 # SM_V3    4 1192    0  1   0.06   -0.06 0.97 -1.44 3.07  4.51 0.42    -0.02 0.03






 # Method 1 

 vars = c("SM_38","SM_V1","SM_V2","SM_V3")
 cor(df2_long[vars],use = "complete.obs")
  #           SM_38     SM_V1     SM_V2     SM_V3
  # SM_38 1.0000000 0.2505290 0.1883033 0.2384978
  # SM_V1 0.2505290 1.0000000 0.5277821 0.4893683
  # SM_V2 0.1883033 0.5277821 1.0000000 0.6214725
  # SM_V3 0.2384978 0.4893683 0.6214725 1.0000000

 vars = c("SM_38","SM_V1","SM_V2","SM_V3")
 psych::describe(df2_long[vars])
 # Pairwise estimates
 corr = polypairwise(df2_long[vars],printFit=T)
 corr$R
 #           SM_38     SM_V1     SM_V2     SM_V3
 # SM_38 1.0000000 0.2584764 0.2146674 0.2409272
 # SM_V1 0.2584764 1.0000000 0.5818389 0.5048188
 # SM_V2 0.2146674 0.5818389 1.0000000 0.6086149
 # SM_V3 0.2409272 0.5048188 0.6086149 1.0000000

 # Complete observations
 corr = polychoricMatrix(df2_long[vars],useMeans=T)	# polychoricMeansMatrix3.R
 corr$polychorics
  #           SM_38     SM_V1     SM_V2     SM_V3
  # SM_38 1.0000000 0.2708021 0.2161412 0.5738691
  # SM_V1 0.2708021 1.0000000 0.2379429 0.5071070
  # SM_V2 0.2161412 0.2379429 1.0000000 0.6068636
  # SM_V3 0.5738691 0.5071070 0.6068636 1.0000000


 # Method 2 FIML

 vars 		= c("SM_38","SM_V1","SM_V2","SM_V3")
 nv 			= length(vars)
 labVect 	= function(lab,nv) { paste(lab,1:nv, sep="") }
 labSymm 	= function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") }
 mlabs = labVect("u",4)
 elabs = labSymm("r",4)
 evals = round(cov(df2_long[vars],use="complete.obs"),2)
 mvals = round(colMeans(df2_long[vars],na.rm = T),2) 
 
  corrs_FIML = 
  mxModel( "phen_corrs",
  mxMatrix( name="Means", type = "Full", nrow= 1, ncol=4, free=T, values=mvals, labels=mlabs),   
  mxMatrix( name="E", 		type = "Symm", nrow=nv, ncol=nv, free=T, values=evals, labels=elabs),     
  mxAlgebra(name="expCov",	expression = E),  
  mxAlgebra(name="corr",	expression = cov2cor(expCov)),
  mxData( 	df2_long[vars], 	type ="raw"),
  mxFitFunctionML(),
  mxCI(c("corr[2,1]",
  		   "corr[3,1]","corr[3,2]",
  		   "corr[4,1]","corr[4,2]","corr[4,3]")),
  mxExpectationNormal(covariance="expCov", means="Means", dimnames=vars)) 
  omxGetParameters(corrs_FIML)
  
 	corrs_FIML_fit = mxTryHard(corrs_FIML, intervals=T)
  summary(corrs_FIML_fit)
	print(summary(corrs_FIML_fit)$CI)
  #                        lbound  estimate    ubound 
  # phen_corrs.corr[2,1] 0.1965146 0.2708023 0.3412139     
  # phen_corrs.corr[3,1] 0.1622351 0.2161412 0.2686852     
  # phen_corrs.corr[3,2] 0.5063437 0.5738690 0.6324446     
  # phen_corrs.corr[4,1] 0.1838338 0.2379430 0.2902017     
  # phen_corrs.corr[4,2] 0.4332665 0.5071069 0.5724679     
  # phen_corrs.corr[4,3] 0.5661217 0.6068636 0.6441154 

	FIML_corrs = corrs_FIML_fit$algebras$corr$result
	rownames(FIML_corrs) = vars
	colnames(FIML_corrs) = vars 	 
  #           SM_38     SM_V1     SM_V2     SM_V3
  # SM_38 1.0000000 0.2708023 0.2161412 0.2379430
  # SM_V1 0.2708023 1.0000000 0.5738690 0.5071069
  # SM_V2 0.2161412 0.5738690 1.0000000 0.6068636
  # SM_V3 0.2379430 0.5071069 0.6068636 1.0000000



#-------------------------------------------------------------------------------------------------------------------#
#TWIN MODEL
#-------------------------------------------------------------------------------------------------------------------#


# ------------------------------------------------------------------------------------------------------------------ #
# SM_38
df2<-newtwins
hist(df2$SM_38_T1)
hist(df2$SM_38_T2)
table(df2$SM_38_T1,df2$SM_38_T2)

selVars <- c("SM_38_T1",
             "SM_38_T2")
#df2$SM_38_T1 = car::recode(df2$SM_38_T1, "4:5=3") 
#df2[,c("SM_38_T1","SM_38_T2")] = mxFactor(df2[,c("SM_38_T1","SM_38_T2")], levels= c(1,2))						 

mzdata	<- subset(df2, ZYG2019==1,selVars); # dim(mzdata) # 481  2
dzdata	<- subset(df2, ZYG2019==2,selVars); # dim(dzdata) # 342  2

table(mzdata$SM_38_T1)
table(mzdata$SM_38_T2)

nv     	= length(selVars)/2     
ntv    	= nv*2    
aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 

# Generate Descriptive Statistics
mzcov <- cov(mzdata,use="complete")
dzcov <- cov(dzdata,use="complete")


# Set Starting Values
avals   <- sqrt(2*(mzcov[1,2]-dzcov[1,2]))
cvals 	<- sqrt(2*dzcov[1,2]-mzcov[1,2])
evals 	<- sqrt(mzcov[2,2]-mzcov[1,2])

uni_SM_38 = mxModel("ACE",
 mxModel("top",
 mxMatrix(	name="Mean", 		type="Full", nrow=1, ncol=nv, free=T, labels="mean1",values = 0),
 mxAlgebra(	name="expMean", 	cbind(Mean,Mean)), 
 mxMatrix(	name="A",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=avals, labels=aLabs, lbound=-50, ubound=50),
 mxMatrix(	name="C",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, labels=cLabs, lbound=-50, ubound=50),
 mxMatrix(	name="E",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=evals, labels=eLabs, lbound=-50, ubound=50),
 mxAlgebra(	name="expCovMZ",	expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
 mxAlgebra(	name="expCovDZ",	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
 mxAlgebra(	name="VC",			expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))),
 mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML()), 
 mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
 mxCI(c("top.VC[1,4]","top.VC[1,5]","top.VC[1,6]")),
 mxFitFunctionMultigroup(c("MZ","DZ"))  )
 omxGetParameters(uni_SM_38)
 ##mxCheckIdentification(uni_SM_38)

# ACE 
SM_38_ACE_fit 	<- mxTryHard( uni_SM_38,extraTries=35, greenOK=FALSE,checkHess=FALSE,fit2beat=Inf,intervals=T)
summary( SM_38_ACE_fit )  
summary(omxRunCI(SM_38_ACE_fit, optimizer = "SLSQP"),verbose=T) 
print(summary(SM_38_ACE_fit)$CI)
 #                 lbound   estimate    ubound note
 # top.VC[1,4] -0.1900626 0.09146928 0.3782659     
 # top.VC[1,5] -0.1173160 0.12682184 0.3585739     
 # top.VC[1,6]  0.6945211 0.78170888 0.8741812     

# AE model
SM_38_AE 				<- mxModel( SM_38_ACE_fit,name="AE")
SM_38_AE 				<- omxSetParameters( SM_38_AE, cLabs, free=F, values=0)
SM_38_AE_fit 			<- mxTryHard( 		 SM_38_AE, extraTries=25, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_38_AE_fit, verbose=F )
SM_38_AE_fit$top$algebras$VC$result
print(summary(SM_38_AE_fit)$CI)
SM_38_AE_fit$top$expCovMZ
SM_38_AE_fit$top$expCovDZ

# CE model
SM_38_CE 					<- mxModel( SM_38_ACE_fit,name="CE")
SM_38_CE 					<- omxSetParameters( SM_38_CE, aLabs, free=F, values=0)
SM_38_CE 					<- omxSetParameters( SM_38_CE, cLabs, free=T, lbound=-20, ubound=20, values=10)
SM_38_CE_fit 				<- mxTryHard( 		 SM_38_CE, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=F)
summary( SM_38_CE_fit ) 

# E model
SM_38_E 					<- mxModel( SM_38_ACE_fit,name="E")
SM_38_E 					<- omxSetParameters( SM_38_E, c(aLabs,cLabs),free=F,values=0)
SM_38_E_fit 				<- mxTryHard( 		 SM_38_E, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_38_E_fit ) 

# Comparisons
subs 						<- c( SM_38_AE_fit, SM_38_CE_fit, SM_38_E_fit )
comps						<- mxCompare( SM_38_ACE_fit,subs );comps
summary(SM_38_ACE_fit) 
summary(SM_38_AE_fit)
summary(SM_38_CE_fit)
summary(SM_38_E_fit)


#   base comparison ep minus2LL   df      AIC     diffLL diffdf            p
# 1  ACE       <NA>  4 4329.744 1533 1263.744         NA     NA           NA
# 2  ACE         AE  3 4330.805 1534 1262.805  1.0608632      1 3.030184e-01
# 3  ACE         CE  3 4330.146 1534 1262.146  0.4021928      1 5.259590e-01
# 4  ACE          E  2 4358.988 1535 1288.988 29.2436879      2 4.464922e-07




# SM_V1

selVars <- c("SM_V1_T1",
             "SM_V1_T2")
						 
#hist(df2$SM_V1_T1)
#hist(df2$SM_V1_T2)

mzdata	<- subset(df2, ZYG2019==1,selVars); # dim(mzdata) # 481  2
dzdata	<- subset(df2, ZYG2019==2,selVars); # dim(dzdata) # 342  2


nv     	= length(selVars)/2     
ntv    	= nv*2    
aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 

# Generate Descriptive Statistics
mzcov <- cov(mzdata,use="complete")
dzcov <- cov(dzdata,use="complete")


# Set Starting Values
avals      <- sqrt(2*(mzcov[1,2]-dzcov[1,2]))
cvals <- sqrt(2*dzcov[1,2]-mzcov[1,2])
evals <- sqrt(mzcov[2,2]-mzcov[1,2])


uni_SM_V1 = mxModel("ACE",
 mxModel("top",
 mxMatrix(	name="Mean", 		type="Full", nrow=1, ncol=nv, free=T, labels="mean1",values = 0),
 mxAlgebra(	name="expMean", 	cbind(Mean,Mean)), 
 mxMatrix(	name="A",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=avals, labels=aLabs, lbound=-50, ubound=50),
 mxMatrix(	name="C",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, labels=cLabs, lbound=-50, ubound=50),
 mxMatrix(	name="E",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=evals, labels=eLabs, lbound=-50, ubound=50),
 mxAlgebra(	name="expCovMZ",	expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
 mxAlgebra(	name="expCovDZ",	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
 mxAlgebra(	name="VC",			expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))),
 mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML()), 
 mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
 mxCI(c("top.VC[1,4]","top.VC[1,5]","top.VC[1,6]")),
 mxFitFunctionMultigroup(c("MZ","DZ"))  )
 omxGetParameters(uni_SM_V1)
 #mxCheckIdentification(uni_SM_V1) 

# ACE 
SM_V1_ACE_fit 	<- mxTryHard( uni_SM_V1)
summary( SM_V1_ACE_fit )  
summary(omxRunCI(SM_V1_ACE_fit, optimizer = "SLSQP"),verbose=T) 
print(summary(SM_V1_ACE_fit)$CI)

# AE model
SM_V1_AE 				<- mxModel( SM_V1_ACE_fit,name="AE")
SM_V1_AE 				<- omxSetParameters( SM_V1_AE, cLabs, free=F, values=0)
SM_V1_AE_fit 			<- mxTryHard( 		 SM_V1_AE, extraTries=25, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_V1_AE_fit, verbose=T )
SM_V1_AE_fit$top$algebras$VC$result
print(summary(SM_V1_AE_fit)$CI)


# CE model
SM_V1_CE 					<- mxModel( SM_V1_ACE_fit,name="CE")
SM_V1_CE 					<- omxSetParameters( SM_V1_CE, aLabs, free=F, values=0)
SM_V1_CE 					<- omxSetParameters( SM_V1_CE, cLabs, free=T, lbound=-20, ubound=20, values=10)
SM_V1_CE_fit 				<- mxTryHard( 		 SM_V1_CE, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=F)
summary( SM_V1_CE_fit ) 

# E model
SM_V1_E 					<- mxModel( SM_V1_ACE_fit,name="E")
SM_V1_E 					<- omxSetParameters( SM_V1_E, c(aLabs,cLabs),free=F,values=0)
SM_V1_E_fit 				<- mxTryHard( 		 SM_V1_E, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_V1_E_fit ) 

# Comparisons
subs 						<- c( SM_V1_AE_fit, SM_V1_CE_fit, SM_V1_E_fit )
comps						<- mxCompare( SM_V1_ACE_fit,subs );comps
#   base comparison ep minus2LL  df      AIC      diffLL diffdf            p
# 1  ACE       <NA>  4 1455.840 516 423.8399          NA     NA           NA
# 2  ACE         AE  3 1455.867 517 421.8667  0.02685081      1 0.8698396599
# 3  ACE         CE  3 1457.999 517 423.9986  2.15873372      1 0.1417614721
# 4  ACE          E  2 1473.695 518 437.6947 17.85477988      2 0.0001327039
summary(SM_V1_ACE_fit) 
summary(SM_V1_AE_fit)
summary(SM_V1_CE_fit)
summary(SM_V1_E_fit)
# SM_V2

selVars <- c("SM_V2_T1",
             "SM_V2_T2")

mzdata	<- subset(df2, ZYG2019==1,selVars); # dim(mzdata) # 481  2
dzdata	<- subset(df2, ZYG2019==2,selVars); # dim(dzdata) # 342  2


nv     	= length(selVars)/2     
ntv    	= nv*2    
aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 

# Generate Descriptive Statistics
mzcov <- cov(mzdata,use="complete")
dzcov <- cov(dzdata,use="complete")


# Set Starting Values
avals <- sqrt(2*(mzcov[1,2]-dzcov[1,2]))
cvals <- sqrt(2*dzcov[1,2]-mzcov[1,2])
evals <- sqrt(mzcov[2,2]-mzcov[1,2])


uni_SM_V2 = mxModel("ACE",
 mxModel("top",
 mxMatrix(	name="Mean", 		type="Full", nrow=1, ncol=nv, free=T, labels="mean1",values = 0),
 mxAlgebra(	name="expMean", 	cbind(Mean,Mean)), 
 mxMatrix(	name="A",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=avals, labels=aLabs, lbound=-50, ubound=50),
 mxMatrix(	name="C",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, labels=cLabs, lbound=-50, ubound=50),
 mxMatrix(	name="E",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=evals, labels=eLabs, lbound=-50, ubound=50),
 mxAlgebra(	name="expCovMZ",	expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
 mxAlgebra(	name="expCovDZ",	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
 mxAlgebra(	name="VC",			expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))),
 mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML()), 
 mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
 mxCI(c("top.VC[1,4]","top.VC[1,5]","top.VC[1,6]")),
 mxFitFunctionMultigroup(c("MZ","DZ"))  )
 omxGetParameters(uni_SM_V2)
 #mxCheckIdentification(uni_SM_V2) 

# ACE 
SM_V2_ACE_fit 	<- mxTryHard( uni_SM_V2,extraTries=35, greenOK=FALSE,checkHess=FALSE,fit2beat=Inf,intervals=T)
summary( SM_V2_ACE_fit )  
summary(omxRunCI(SM_V2_ACE_fit, optimizer = "SLSQP"),verbose=T) 
print(summary(SM_V2_ACE_fit)$CI)

# AE model
SM_V2_AE 				<- mxModel( SM_V2_ACE_fit,name="AE")
SM_V2_AE 				<- omxSetParameters( SM_V2_AE, cLabs, free=F, values=0)
SM_V2_AE_fit 			<- mxTryHard( 		 SM_V2_AE, extraTries=25, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_V2_AE_fit, verbose=T )
SM_V2_AE_fit$top$algebras$VC$result
print(summary(SM_V2_AE_fit)$CI)


# CE model
SM_V2_CE 					<- mxModel( SM_V2_ACE_fit,name="CE")
SM_V2_CE 					<- omxSetParameters( SM_V2_CE, aLabs, free=F, values=0)
SM_V2_CE 					<- omxSetParameters( SM_V2_CE, cLabs, free=T, lbound=-20, ubound=20, values=10)
SM_V2_CE_fit 	<- mxTryHard( 		 SM_V2_CE, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=F)
summary( SM_V2_CE_fit ) 

# E model
SM_V2_E 					<- mxModel( SM_V2_ACE_fit,name="E")
SM_V2_E 					<- omxSetParameters( SM_V2_E, c(aLabs,cLabs),free=F,values=0)
SM_V2_E_fit 				<- mxTryHard( 		 SM_V2_E, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_V2_E_fit ) 

# Comparisons
subs 						<- c( SM_V2_AE_fit, SM_V2_CE_fit, SM_V2_E_fit )
comps						<- mxCompare( SM_V2_ACE_fit,subs );comps

 #   base comparison ep minus2LL   df       AIC     diffLL diffdf            p
 # 1  ACE       <NA>  4 3371.061 1195  981.0609         NA     NA           NA
 # 2  ACE         AE  3 3371.336 1196  979.3361  0.2752145      1 5.998551e-01
 # 3  ACE         CE  3 3376.238 1196  984.2381  5.1772204      1 2.288491e-02
 # 4  ACE          E  2 3401.614 1197 1007.6142 30.5533190      2 2.319696e-07
summary(SM_V2_ACE_fit) 
summary(SM_V2_AE_fit)
summary(SM_V2_CE_fit)
summary(SM_V2_E_fit)

# SM_V3

selVars <- c("SM_V3_T1",
             "SM_V3_T2")

mzdata	<- subset(df2, ZYG2019==1,selVars); # dim(mzdata) # 481  2
dzdata	<- subset(df2, ZYG2019==2,selVars); # dim(dzdata) # 342  2


nv     	= length(selVars)/2     
ntv    	= nv*2    
aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 

# Generate Descriptive Statistics
mzcov <- cov(mzdata,use="complete")
dzcov <- cov(dzdata,use="complete")


# Set Starting Values
avals      <- sqrt(2*(mzcov[1,2]-dzcov[1,2]))
cvals <- sqrt(2*dzcov[1,2]-mzcov[1,2])
evals <- sqrt(mzcov[2,2]-mzcov[1,2])


uni_SM_V3 = mxModel("ACE",
 mxModel("top",
 mxMatrix(	name="Mean", 		type="Full", nrow=1, ncol=nv, free=T, labels="mean1",values = 0),
 mxAlgebra(	name="expMean", 	cbind(Mean,Mean)), 
 mxMatrix(	name="A",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=avals, labels=aLabs, lbound=-50, ubound=50),
 mxMatrix(	name="C",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, labels=cLabs, lbound=-50, ubound=50),
 mxMatrix(	name="E",			type="Symm", nrow=nv, ncol=nv, free=TRUE, values=evals, labels=eLabs, lbound=-50, ubound=50),
 mxAlgebra(	name="expCovMZ",	expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
 mxAlgebra(	name="expCovDZ",	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
 mxAlgebra(	name="VC",			expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))),
  mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML()), 
  mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
  mxCI(c("top.VC[1,4]","top.VC[1,5]","top.VC[1,6]")),
  mxFitFunctionMultigroup(c("MZ","DZ"))  )
# omxGetParameters(uni_MDS_51)

# ACE 
SM_V3_ACE_fit 	<- mxTryHard( uni_SM_V3,extraTries=35, greenOK=FALSE,checkHess=FALSE,fit2beat=Inf,intervals=T)
summary( SM_V3_ACE_fit )  
summary(omxRunCI(SM_V3_ACE_fit, optimizer = "SLSQP"),verbose=T) 
print(summary(SM_V3_ACE_fit)$CI)

# AE model
SM_V3_AE 				<- mxModel( SM_V3_ACE_fit,name="AE")
SM_V3_AE 				<- omxSetParameters( SM_V3_AE, cLabs, free=F, values=0)
SM_V3_AE_fit 			<- mxTryHard( 		 SM_V3_AE, extraTries=25, greenOK=TRUE,checkHess=FALSE,intervals=T)
SM_V3_AE <- mxAutoStart(SM_V3_AE)
summary( SM_V3_AE_fit 	<- mxTryHard( 		 SM_V3_AE, extraTries=25, greenOK=TRUE,checkHess=FALSE,intervals=T))

summary( SM_V3_AE_fit, verbose=T )
SM_V3_AE_fit$top$algebras$VC$result
print(summary(SM_V3_AE_fit)$CI)


# CE model
SM_V3_CE 					<- mxModel( SM_V3_ACE_fit,name="CE")
SM_V3_CE 					<- omxSetParameters( SM_V3_CE, aLabs, free=F, values=0)
SM_V3_CE 					<- omxSetParameters( SM_V3_CE, cLabs, free=T, lbound=-20, ubound=20, values=10)
SM_V3_CE_fit 				<- mxTryHard( 		 SM_V3_CE, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=F)
summary( SM_V3_CE_fit ) 

# E model
SM_V3_E 					<- mxModel( SM_V3_ACE_fit,name="E")
SM_V3_E 					<- omxSetParameters( SM_V3_E, c(aLabs,cLabs),free=F,values=0)
SM_V3_E_fit 				<- mxTryHard( 		 SM_V3_E, extraTries=5, greenOK=TRUE,checkHess=FALSE,intervals=T)
summary( SM_V3_E_fit ) 

# Comparisons
subs 						<- c( SM_V3_AE_fit, SM_V3_CE_fit, SM_V3_E_fit )
comps						<- mxCompare( SM_V3_ACE_fit,subs );comps
# base comparison ep  minus2LL   df       AIC     diffLL diffdf             p
# 1  ACE       <NA>  4 3332.2058 1188 3340.2058         NA     NA            NA
# 2  ACE         AE  3 3336.5440 1189 3342.5440  4.3382185      1 3.7265896e-02
# 3  ACE         CE  3 3348.9492 1189 3354.9492 16.7434383      1 4.2789992e-05
# 4  ACE          E  2 3381.3643 1190 3385.3643 49.1585151      2 2.1152616e-11
summary(SM_V3_ACE_fit) 
summary(SM_V3_AE_fit)
summary(SM_V3_CE_fit)
summary(SM_V3_E_fit)


#---------------------------------------------------------------------------------------------------#
# MULTIVARIATE ANALYSES
#---------------------------------------------------------------------------------------------------#

# Correlated factors


selVars    = c("SM_38_T1","SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
               "SM_38_T2","SM_V1_T2", "SM_V2_T2", "SM_V3_T2") 
nv        	= 4                   		
ntv       	= nv*2
mzdata		= subset(df2, ZYG2019==1,selVars) 
dzdata		= subset(df2, ZYG2019==2,selVars) 

psych::describe(mzdata)
psych::describe(dzdata)

nv     	= length(selVars)/2     
ntv    	= nv*2    
thVals 	= 0.5
aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 


multi_ace = mxModel("Correlated factors",
 mxModel("top",
 mxMatrix( name="Means", 		type="Full", nrow=1, ncol=nv, free=c(T,T,T,T), labels=c("m1","m2", "m3","m4")), 
 mxAlgebra(name="expMean", 		cbind(Means,Means)), # Mean matrices for twin 1 and twin 2
 mxMatrix( name="inc",    		type="Lower", nrow=1,  ncol=1,  free=F, values=1),
 mxMatrix( name="A",			type="Symm",  nrow=nv, ncol=nv, free=T, values=0, 	labels=aLabs, lbound= -20, ubound= 20),
 mxMatrix( name="C",			type="Symm",  nrow=nv, ncol=nv, free=T, values = 0, labels=cLabs, lbound= -25, ubound=20),
 mxMatrix( name="E",			type="Symm",  nrow=nv, ncol=nv, free=T, 			labels=eLabs, lbound= -20, ubound=20),
 mxAlgebra(name="expCovMZ",		expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
 mxAlgebra(name="expCovDZ",		expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
 mxAlgebra(name="corrP",		expression=cov2cor(A+C+E)),
 mxAlgebra(name="corrA",		expression=cov2cor(A)),
 mxAlgebra(name="corrE",		expression=cov2cor(E)),
   mxCI(c("VC[1,13]","VC[2,14]","VC[3,15]","VC[3,16]",
          "VC[1,17]","VC[2,18]","VC[3,19]","VC[3,20]",
          "VC[1,21]","VC[2,22]","VC[3,23]","VC[3,24]")),
 mxAlgebra(name="VC", expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)) )),
 mxModel("MZ", mxData(mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML() ),
 mxModel("DZ", mxData(dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML() ),
 mxFitFunctionMultigroup(c("MZ","DZ")))
omxGetParameters(multi_ace)
#mxCheckIdentification(multi_ace)

# Run ULS for better start values
 uls 				<- mxAutoStart( multi_ace )
 uls_fit			<- mxTryHard( uls )
 summary( uls_fit )

 	# uls_fit$top$algebras$expCovDZ
 	# uls_fit$top$algebras$expCovMZ
 	 round(multi_ace_fit$top$algebras$VC$result,2)

# ACE 
 multi_ace_fit 				<- mxTryHard(uls_fit)
 summary(omxRunCI(multi_ace_fit, optimizer = "SLSQP"),verbose=F) 
 summary( multi_ace_fit ) 
 print(summary( multi_ace_fit )$CI)

# AE model
 multi_ae 					<- mxModel( multi_ace_fit,name="AE")
 multi_ae 					<- omxSetParameters( multi_ae, label=cLabs,free=F,values= 0)
 multi_ae_fit				<- mxTryHard( multi_ae )
 summary( multi_ae_fit ) 
 round(multi_ae_fit$top$algebras$VC$result,2)
# CE model
 multi_ce 					<- mxModel( multi_ace_fit,name="CE")
 multi_ce 					<- omxSetParameters( multi_ce, label=aLabs,free=F,values= 0)
 multi_ce_fit				<- mxTryHard( multi_ce )
 summary( multi_ce_fit ) 

# E model
 multi_e 	 				<- mxModel( multi_ace_fit,name="E")
 multi_e 	 				<- omxSetParameters( multi_e, label=aLabs,free=F,values= 0)
 multi_e 					<- omxSetParameters( multi_e, label=cLabs,free=F,values= 0)
 multi_e_fit				<- mxTryHard( multi_e )
 summary( multi_e_fit ) 

# Compare models
 mxCompare(multi_ace_fit,c(multi_ae_fit, multi_ce_fit, multi_e_fit)) 

 #   base comparison ep minus2LL   df      AIC     diffLL diffdf            p
 # 1  ACE       <NA> 34 11788.83 4414 11856.83         NA     NA           NA
 # 2  ACE         AE 24 11796.63 4424 11844.63   7.803988     10 6.479763e-01
 # 3  ACE         CE 24 11812.58 4424 11860.58  23.753844     10 8.281267e-03
 # 4  ACE          E 14 11904.71 4434 11932.71 115.887029     20 1.638649e-15





# Common pathway

selVars    = c("SM_38_T1","SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
               "SM_38_T2","SM_V1_T2", "SM_V2_T2", "SM_V3_T2")

mzdata		<- subset(df2, ZYG2019==1, selVars); # dim(mzdata) # 475  10
dzdata		<- subset(df2, ZYG2019==2, selVars); # dim(dzdata) # 335  10

nv     	= length(selVars)/2     
ntv    	= nv*2    


nVariables	= 4
nFactors	= 1
#thVals 	= 0.5
#nth 		= 1
#svLTh  	= 1.1    # start value for first threshold
#svITh  	= 1.0    # start value for increments
#svTh   	= matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv)     # start value for thresholds
#round(auto_ACE_fit$matrices$thresh$values,2)
#Th_val		= c(1.2,1.1,1.1)
#Th_lb   	= matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)     # lower bounds for thresholds
#Th_lab  	= c(paste("th",1:nth,"v1",sep=""),paste("th",1:nth,"v2",sep=""),paste("th",1:nth,"v3",sep=""))
psi_lab_a	= c("ia11")
psi_lab_c	= c("ic11")
psi_lab_e	= c("ie11")
#round(t(diag2vec(auto_ACE_fit$matrices$psi_c$values)),2)
psi_a_val	= 0.7
psi_c_val	= 0.2
psi_e_val	= 0.2
res_lab_a	= c("res_a1","res_a2","res_a3","res_a4")
res_lab_c	= c("res_c1","res_c2","res_c3","res_c4")
res_lab_e	= c("res_e1","res_e2","res_e3","res_e4")
loadS 		= 0.8
loadF 		= T
lamba_lab	= c("f11","f21","f31","f41")

cp_ACE = mxModel("CP",
 mxModel("top",
 mxMatrix(name="Mean", 		type = "Full", nrow=1, ncol=nv, free = T, labels = c("m1","m2","m3","m4"), lbound = -5, ubound = 5 ),
 mxAlgebra(name="expMean", 	expression = cbind(Mean, Mean)),  
 mxMatrix(name="psi_a", 	type = "Symm", nrow = nFactors, ncol = nFactors, free = T, labels = psi_lab_a, values = 0, lbound = -5, ubound = 5),
 mxMatrix(name="psi_c", 	type = "Symm", nrow = nFactors, ncol = nFactors, free = T, labels = psi_lab_c, values = 0, lbound = -5, ubound = 5),
 mxMatrix(name="psi_e", 	type = "Symm", nrow = nFactors, ncol = nFactors, free = T, labels = psi_lab_e, values = 1, lbound = -5, ubound = 5),
 mxMatrix(name="lamba",		type = "Full", nrow = nv, ncol = nFactors, free = T, labels= lamba_lab, values = 1, lbound = 0, ubound = 20),
 mxMatrix(name="epsilon_a", type = "Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_a, values =   0, lbound = -5, ubound = 5), 
 mxMatrix(name="epsilon_c", type = "Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_c, values =   0, lbound = -5, ubound = 5), 
 mxMatrix(name="epsilon_e", type = "Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_e, values = 0.1, lbound = -5, ubound = 5), 
 mxAlgebra(name="A", 		expression = lamba %&% psi_a + epsilon_a),
 mxAlgebra(name="C", 		expression = lamba %&% psi_c + epsilon_c),
 mxAlgebra(name="E", 		expression = lamba %&% psi_e + epsilon_e),  
 mxAlgebra(name="expCovMZ", expression = rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
 mxAlgebra(name="expCovDZ", expression = rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
 mxAlgebra(name="corrP",	expression = cov2cor(A+C+E)),
 mxAlgebra(name="corrA", 	expression = cov2cor(A)),
 mxAlgebra(name="corrE", 	expression = cov2cor(E)),
# Standardization - constrain variance of CP = 1
mxMatrix(name="unitM",type="Unit", nrow=nFactors, ncol=nFactors),
mxConstraint(name="ConVar",expression=diag2vec(psi_a+psi_c+psi_e)==unitM),  
mxCI(c("psi_a","psi_e", "ia11","ie11","f11","f21","f31","f41")), 
#"ia11","ie11", "res_a1","res_a2","res_a3","res_a4","res_e1","res_e2","res_e3","res_e4", "f11","f21","f31","f41",
# "VC[1,13]","VC[2,14]","VC[3,15]","VC[4,16]", "VC[1,21]","VC[2,22]","VC[3,23]","VC[4,24]","corrA[2,1]"
mxAlgebra(name="VC",expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))), 
# Multiple groups
mxModel("MZ", mxData(mzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),   
mxModel("DZ", mxData(dzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),   
mxFitFunctionMultigroup(c("MZ","DZ")))
omxGetParameters(cp_ACE)
#mxCheckIdentification(cp_ACE) # = GOOD & identified

# ACE 
 cp_ACE_fit 		<- mxTryHard(cp_ACE)
 summary( cp_ACE_fit )  
 cp_ACE_fit
 cp_ACE_fit$top$algebras$VC$result 
 #summary(omxRunCI(cp_ACE_fit,optimizer = "SLSQP"),verbose=T)$CI	
 

# AE model
cp_1F_AE 			<- mxModel( cp_ACE_fit,name="AE")
cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label=psi_lab_c,free=F,values= 0)
cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label=res_lab_c,free=F,values= 0) 
cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label=c("f11","f21","f31","f41"),	free=T,values=0.5, ubound=10)  
#cp_1F_AE	 		<- omxSetParameters( cp_1F_AE, label=c(res_lab_a,res_lab_e),		free=T,values=0.5, ubound=20)    
#cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label="ia11",free=T,values=20,ubound=40) 
cp_1F_AE_fit		<- mxTryHard( 		 cp_1F_AE )
summary( cp_1F_AE_fit ) 

# CE model
cp_1F_CE 			<- mxModel( cp_ACE_fit,name="CE")
cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label=psi_lab_a,free=F,values= 0)
cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label=res_lab_a,free=F,values= 0) 
cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label=c("f11","f21","f31","f41"),	free=T,values=0.5, ubound=10)  
#cp_1F_CE	 		<- omxSetParameters( cp_1F_CE, label=c(res_lab_c,res_lab_e),		free=T,values=0.5, ubound=20)    
#cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label="ia11",free=T,values=20,ubound=40) 
cp_1F_CE_fit		<- mxTryHard( 		 cp_1F_CE )
summary( cp_1F_CE_fit ) 

# E model
cp_1F_E 			<- mxModel( cp_ACE_fit,name="E")
cp_1F_E 			<- omxSetParameters( cp_1F_E, label=psi_lab_a,free=F,values= 0)
cp_1F_E 			<- omxSetParameters( cp_1F_E, label=res_lab_a,free=F,values= 0) 
cp_1F_E 			<- omxSetParameters( cp_1F_E, label=psi_lab_c,free=F,values= 0)
cp_1F_E 			<- omxSetParameters( cp_1F_E, label=res_lab_c,free=F,values= 0) 
cp_1F_E 			<- omxSetParameters( cp_1F_E, label=c("f11","f21","f31","f41"),	free=T,values=0.5, ubound=10)  
cp_1F_E_fit			<- mxTryHard( 		 cp_1F_E )
summary( cp_1F_E_fit ) 

# Compare models
 mxCompare(cp_ACE_fit,c(cp_1F_AE_fit, cp_1F_CE_fit, cp_1F_E_fit)) 
 #   base comparison ep minus2LL   df      AIC     diffLL diffdf            p
 # 1  ACE       <NA> 23 11801.71 4426 11847.71         NA     NA           NA
 # 2  ACE         AE 18 11807.61 4431 11843.61   5.903541      5 3.157183e-01
 # 3  ACE         CE 18 11823.00 4431 11859.00  21.285104      5 7.154911e-04
 # 4  ACE          E 13 11911.06 4436 11937.06 109.353565     10 7.203119e-19






# Independent Pathways

selVars    = c("SM_38_T1","SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
               "SM_38_T2","SM_V1_T2", "SM_V2_T2", "SM_V3_T2")

mzdata		<- subset(df2, ZYG2019==1, selVars); # dim(mzdata) # 475  10
dzdata		<- subset(df2, ZYG2019==2, selVars); # dim(dzdata) # 335  10

nv     	= length(selVars)/2     
ntv    	= nv*2    
nVariables	= 4
nFactors	= 1
#thVals 	= 0.5
#nth 		= 1
#svLTh  	= 1.1    # start value for first threshold
#svITh  	= 1.0    # start value for increments
#svTh   	= matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv)     # start value for thresholds
#round(auto_ACE_fit$matrices$thresh$values,2)
#Th_val		= c(1.2,1.1,1.1)
#Th_lb   	= matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)     # lower bounds for thresholds
#Th_lab  	= c(paste("th",1:nth,"v1",sep=""),paste("th",1:nth,"v2",sep=""),paste("th",1:nth,"v3",sep=""))
psi_lab_a	= c("ia11")
psi_lab_c	= c("ic11")
psi_lab_e	= c("ie11")
#round(t(diag2vec(auto_ACE_fit$matrices$psi_c$values)),2)
psi_a_val	= 0.7
psi_c_val	= 0.2
psi_e_val	= 0.2
res_lab_a	= c("res_a1","res_a2","res_a3","res_a4")
res_lab_c	= c("res_c1","res_c2","res_c3","res_c4")
res_lab_e	= c("res_e1","res_e2","res_e3","res_e4")
loadS 		= 0.8
loadF 		= T

#lamba_lab_a <- matrix(c("f11_a",
#                        "f21_a",
#                        "f31_a",
#                        "f41_a"),nv,byrow = TRUE)
#
lamba_lab_a = c("f11_a","f21_a","f31_a","f41_a")
lamba_lab_c = c("f11_c","f21_c","f31_c","f41_c")
lamba_lab_e = c("f11_e","f21_e","f31_e","f41_e")

ip_ACE = mxModel("IP",
 mxModel("top",
 mxMatrix(name="Mean", 		type="Full", nrow=1, ncol=nv, free = T, labels = c("m1","m2","m3","m4"), lbound = -15, ubound = 5 ),
 mxAlgebra(name="expMean", 	expression= cbind(Mean, Mean)), 
 mxMatrix(name="lamba_a",	type="Full", nrow = nv, ncol = 1,  free = T, labels = lamba_lab_a, 	values = 0.3, lbound = -10, ubound = 15), #   
 mxMatrix(name="lamba_c",	type="Full", nrow = nv, ncol = 1,  free = T, labels = lamba_lab_c, 	values = 0.3, lbound = -10, ubound = 15), #   
 mxMatrix(name="lamba_e",	type="Full", nrow = nv, ncol = 1,  free = T, labels = lamba_lab_e, 	values = 0.3, lbound = -10, ubound = 15), # 
 mxMatrix(name="epsilon_a", type="Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_a, 	values = 0.3, lbound = -10, ubound = 15), #  
 mxMatrix(name="epsilon_c", type="Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_c, 	values = 0.3, lbound = -10, ubound = 15), #  
 mxMatrix(name="epsilon_e", type="Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_e, 	values = 0.3, lbound = -10, ubound = 15), #  
 mxAlgebra(name="A", 		expression= lamba_a %*% t(lamba_a) + epsilon_a %*% t(epsilon_a)),
 mxAlgebra(name="C", 		expression= lamba_c %*% t(lamba_c) + epsilon_c %*% t(epsilon_c)),
 mxAlgebra(name="E", 		expression= lamba_e %*% t(lamba_e) + epsilon_e %*% t(epsilon_e)),  
 mxAlgebra(name="expCovMZ", expression= rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
 mxAlgebra(name="expCovDZ", expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
 mxAlgebra(name="VC",expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))),
 mxModel("MZ",mxData(mzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
 mxModel("DZ",mxData(dzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
 mxFitFunctionMultigroup(c("MZ","DZ")))
 omxGetParameters(ip_ACE)
 #mxCheckIdentification(ip_ACE)

# ACE 
ip_ACE_fit 				<- mxTryHard( ip_ACE )
summary( ip_ACE_fit )  




# Hybrid correlated factors & autoregression


selVars    = c("SM_38_T1","SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
               "SM_38_T2","SM_V1_T2", "SM_V2_T2", "SM_V3_T2")

mzdata			= subset(df2, ZYG2019==1, selVars); # dim(mzdata) # 475  10
dzdata			= subset(df2, ZYG2019==2, selVars); # dim(dzdata) # 335  10
nv        	= 4                   		
ntv       	= nv*2
# mxOption(NULL,"mvnRelEps",0.0045)
nVariables	= nv
nFactors		= nv
labSymm   <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") }
psi_lab_a		= labSymm("a",4)
psi_lab_c		= labSymm("c",4)
psi_lab_e		= labSymm("e",4)

# To identify the model, constratin E effects from SM38 to be equal across time
psi_lab_e=matrix(c("e11","e21","e31","e41",
									 "e21","e22","e32","e42",
									 "e31","e32","e33","e43",
									 "e41","e42","e43","e44"),4,byrow=T)


res_lab_e		= c("res_e1","res_e2","res_e2","res_e2")

psi_free		= matrix(c(T,T,T,T,
											 T,T,F,F,
											 T,F,T,F,
											 T,F,F,T),4,byrow=T)
psi_vals		= matrix(c(0.6,0.6,0.6,0.6,
											 0.6,0.3,0.0,0.0,
											 0.6,0.0,0.3,0.0,
											 0.6,0.0,0.0,0.3),4,byrow=T)


psi_vals_e	= matrix(c(0.3,0.3,0.3,0.3,
											 0.3,0.2,0.0,0.0,
											 0.3,0.0,0.1,0.0,
											 0.3,0.0,0.0,0.1),4,byrow=T)



betaF   		= matrix(c(F,F,F,F, 
          	           F,F,F,F,  
          	           F,T,F,F,
          	           F,F,T,F),nv,byrow = TRUE)
b_lab 			= matrix(c(NA, NA, NA, NA, 
                       NA, NA, NA, NA, 
                       NA,"b2",NA, NA,
                       NA, NA,"b2",NA),nv,byrow = TRUE)
loadS 			= diag(nFactors)
loadF 			= F

# round(t(diag2vec(auto_PBA_ACE_fit$top$matrices$psi_a$values)),1)
psi_a_val	= c( 0.5)
psi_c_val	= c( 0.1)
psi_e_val	= c( 0.3)

# round(t(diag2vec(auto_PBA_ACE_fit$top$matrices$epsilon_a$values)),1)
# epsi_a_val	=c(-0.7,-0.2,-0.2)
# epsi_c_val	=c(0.2,-0.1,0.1)
epsi_e_val	= 0.8

hybrid = mxModel("Hybrid",
mxModel("top",
# Means & definition variables
 mxMatrix( name="Mean", 		type="Full", nrow=1, ncol=nv, free = T, labels = c("m1","m2","m3","m4"), values=0.1, lbound = -5, ubound = 5 ),
 mxAlgebra(name="expMean", expression= cbind(Mean,Mean)), 
# Psi - innovations
 mxMatrix(name="psi_a", 		type="Symm", nrow = 4, ncol = 4, free = psi_free, labels = psi_lab_a, values = psi_vals, lbound = -5, ubound = 5 ), 
 mxMatrix(name="psi_c", 		type="Symm", nrow = 4, ncol = 4, free = psi_free, labels = psi_lab_c, lbound = -5, ubound = 5 ), 
 mxMatrix(name="psi_e", 		type="Symm", nrow = 4, ncol = 4, free = psi_free, labels = psi_lab_e, values = psi_vals_e, lbound = -10, ubound = 5 ), 
# Beta - causal pathways
 mxMatrix(name="beta", 		type="Full", nrow = 4, ncol = 4, free = betaF, labels = b_lab, lbound = -2.5, ubound = 2.5 ),
 mxMatrix(name="I",				type="Iden", nrow = 4, ncol = 4),
# Lamba                    	
 mxMatrix(name="lamba",			type="Full", nrow = nv, ncol = nv, free = F, values = diag(nv)), 
# Epsilon - errors                                                                                   
 mxMatrix(name="epsilon_e", 	type="Diag", nrow = 4, ncol = 4, free = c(F,T,T,T), labels = res_lab_e, values = c(0,0.6,0.6,0.6)),  
# ExpCov = 
 mxAlgebra(name="A", 			expression= lamba %&% (solve(I-beta) %&% psi_a) ), 				# + epsilon_a
 mxAlgebra(name="C", 			expression= lamba %&% (solve(I-beta) %&% psi_c) ), 				# + epsilon_c
 mxAlgebra(name="E", 			expression= lamba %&% (solve(I-beta) %&% psi_e) + epsilon_e), 	# 
 mxAlgebra(name="expCovMZ", 	expression= rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
 mxAlgebra(name="expCovDZ", 	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
 mxAlgebra(name="corrP",		expression= cov2cor(A+C+E)),
 mxAlgebra(name="corrA", 		expression= cov2cor(A)),
 mxAlgebra(name="corrE", 		expression= cov2cor(E)), 
 mxAlgebra(name="VC",			expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv))), 
 mxCI(c("ia11","ia22","ia33","ia44","ic11","ic22","ic33","ic44","ie11","ie22","ie33","ie44","res_e", "b") )), 
 mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML() ),  
 mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML() ),
 mxFitFunctionMultigroup(c("MZ","DZ"))
 )
omxGetParameters(hybrid)
#mxCheckIdentification( hybrid_fit )

 # options(max.print = 3000)
 # auto_uls 				<- mxAutoStart(hybrid)
 # summary(auto_fit 	<- mxTryHard(auto_uls))
 # auto_fit$top$algebras$VC$result
 # 
 # auto_fit$top$algebras$corrP$result
 # auto_fit$top$algebras$corrA$result
 # auto_fitci<-summary(omxRunCI(auto_fit, optimizer = "SLSQP"),verbose=T)

# ACE 
 hybrid_fit 				<- mxTryHard( auto_fit )
 summary( hybrid_fit )
 hybrid_fit$top$matrices$psi_a
 hybrid_fit$top$matrices$beta$values





 
 
 # Autoregression
 
 selVars    = c("SM_38_T1","SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
                "SM_38_T2","SM_V1_T2", "SM_V2_T2", "SM_V3_T2")
 
 mzdata			= subset(df2, ZYG2019==1, selVars); # dim(mzdata) # 475  10
 dzdata			= subset(df2, ZYG2019==2, selVars); # dim(dzdata) # 335  10
 nv        	= 4                   		
 ntv       	= nv*2
 # mxOption(NULL,"mvnRelEps",0.0045)
 nVariables	= nv
 nFactors		= nv
 
 psi_lab_a		= c("ia11","ia22","ia33","ia44")
 psi_lab_c		= c("ic11","ic22","ic33","ic44")
 psi_lab_e		= c("ie11","ie22","ie33","ie44")
 
 # Residual unexplained by model. Residual at age 38 is different from residuals at VETSA 1-3
 res_lab_e		= c("res_e1","res_e2","res_e2","res_e2")
 
 # Beta regression from age 38 to VETSA is different.
 betaF   	= matrix(c(  F,F,F,F, 
                       T,F,F,F,  
                       F,T,F,F,
                       F,F,T,F),nv,byrow = TRUE)
 b_lab 		= matrix(c(NA, NA, NA, NA, 
                    "b1",NA, NA, NA, 
                    NA,"b2", NA, NA,
                    NA, NA,"b2", NA),nv,byrow = TRUE)
 loadS 		= diag(nFactors)
 loadF 		= F
 
 
 auto_ACE = mxModel("Auto",
                    mxModel("top",
                            # Means & definition variables
                            mxMatrix( name="Mean", 		type="Full", nrow=1, ncol=nv, free = T, labels = c("m1","m2","m3","m4"), values=0.1, lbound = -5, ubound = 5 ),
                            mxAlgebra(name="expMean", expression= cbind(Mean,Mean)), 
                            # Psi - innovations
                            mxMatrix(name="psi_a", 		type="Diag", nrow = 4, ncol = 4, free = T, labels = psi_lab_a, values = 0.5, lbound = -5, ubound = 5 ), 
                            mxMatrix(name="psi_c", 		type="Diag", nrow = 4, ncol = 4, free = T, labels = psi_lab_c, values = 0.5, lbound = -5, ubound = 5 ), 
                            mxMatrix(name="psi_e", 		type="Diag", nrow = 4, ncol = 4, free = T, labels = psi_lab_e, values = 0.5, lbound = -5, ubound = 5 ), 
                            # Beta - causal pathways
                            mxMatrix(name="beta", 		type="Full", nrow = 4, ncol = 4, free = betaF, labels = b_lab, lbound = 0, ubound =5 ),
                            mxMatrix(name="I",				type="Iden", nrow = 4, ncol = 4),
                            # Lamba                    	
                            mxMatrix(name="lamba",			type="Full", nrow = nv, ncol = nv, free = F, values = diag(nv)), 
                            # Epsilon - errors                                                                                   
                            mxMatrix(name="epsilon_e", 	type="Diag", nrow = 4, ncol = 4, free = T, labels = res_lab_e, values = 4),  
                            mxAlgebra(name="A", 			expression= lamba %&% (solve(I-beta) %&% psi_a) ), 				# + epsilon_a
                            mxAlgebra(name="C", 			expression= lamba %&% (solve(I-beta) %&% psi_c) ), 				# + epsilon_c
                            mxAlgebra(name="E", 			expression= lamba %&% (solve(I-beta) %&% psi_e) + epsilon_e), 	# 
                            mxAlgebra(name="expCovMZ", 	expression= rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
                            mxAlgebra(name="expCovDZ", 	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
                            mxAlgebra(name="corrP",		expression= cov2cor(A+C+E)),
                            mxAlgebra(name="corrA", 		expression= cov2cor(A)),
                            mxAlgebra(name="corrE", 		expression= cov2cor(E)), 
                            mxAlgebra(name="VC",			expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv))), 
                            mxCI(c("ia11","ia22","ia33","ia44","ie11","ie22","ie33","ie44","res_e1","res_e2","b1", "b2") )), 
                    mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML() ),  
                    mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars ), mxFitFunctionML() ),
                    mxFitFunctionMultigroup(c("MZ","DZ"))
 )
 omxGetParameters(auto_ACE)
 #mxCheckIdentification( auto_ACE_fit )
 
 # ACE 
 auto_ACE_fit = mxTryHard( auto_ACE )
 summary( auto_ACE_fit )
 
 # AE model
 auto_AE 			<- mxModel( auto_ACE_fit,name="AE")
 auto_AE 			<- omxSetParameters( auto_AE, label=psi_lab_c,free=F,values= 0)
 auto_AE_fit = mxTryHard( auto_AE )
 summary( auto_AE_fit ) 
 
 tmp = summary(omxRunCI(auto_AE_fit, optimizer = "SLSQP")) 
 tmp
 uls<-mxAutoStart(auto_AE)
 auto_AE_fit = mxTryHard( uls )
 summary( auto_AE_fit )
 # CE model
 auto_CE 			<- mxModel( auto_ACE_fit,name="CE")
 auto_CE 			<- omxSetParameters( auto_CE, label=psi_lab_a,free=F,values= 0)
 auto_CE_fit		<- mxTryHard( 		 auto_CE )
 summary( auto_CE_fit ) 
 
 # E model
 auto_E 			<- mxModel( auto_ACE_fit,name="E")
 auto_E 			<- omxSetParameters( auto_E, label=psi_lab_a,free=F,values= 0)
 auto_E 			<- omxSetParameters( auto_E, label=psi_lab_c,free=F,values= 0) 
 auto_E_fit			<- mxTryHard( 		 auto_E )
 summary( auto_E_fit ) 
 
 # Compare models
 mxCompare(auto_ACE_fit,c(auto_AE_fit, auto_CE_fit, auto_E_fit)) 
 #  base comparison ep minus2LL   df      AIC     diffLL diffdf            p
 # 1 Auto       <NA> 20 11804.89 4428 11844.89         NA     NA           NA
 # 2 Auto         AE 16 11806.24 4432 11838.24   1.342687      4 8.540925e-01
 # 3 Auto        AE2 15 11807.06 4433 11837.06   2.164886      5 8.258916e-01
 # 4 Auto         CE 16 11822.89 4432 11854.89  17.992994      4 1.237995e-03
 # 5 Auto          E 12 11909.14 4436 11933.14 104.247305      8 5.769858e-19


summary(cp_ACE_fit)
cp_AE_fit$top$VC
summary(auto_ACE_fit)
summary(auto_AE_fit)
summary(auto_CE_fit)
summary(auto_E_fit)

summary(auto_AE_fit)
auto_AE_fit$top$VC

# Compare multivariate models
#Table S5
 mxCompare(multi_ace_fit,c(cp_ACE_fit, ip_ACE_fit, auto_ACE_fit) )
 #                 base comparison ep minus2LL   df      AIC    diffLL diffdf         p
 # 1 Correlated factors       <NA> 34 11788.83 4414 11856.83        NA     NA        NA
 # 2 Correlated factors         CP 23 11801.71 4426 11847.71 12.885216     12 0.3774387
 # 3 Correlated factors         IP 28 11798.98 4420 11854.98 10.153453      6 0.1183368
 # 4 Correlated factors     Hybrid 27 11794.99 4421 11848.99  6.159983      7 0.5211983
 # 5 Correlated factors       Auto 20 11804.89 4428 11844.89 16.068461     14 0.3092114










#---------------------------------------------------------------------------------------------------#
# MULTIVARIATE ANALYSES - VETSA waves only
#---------------------------------------------------------------------------------------------------#

 # Correlated factors - VETSA waves only


 selVars    = c("SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
                "SM_V1_T2", "SM_V2_T2", "SM_V3_T2") 
 nv        	= 3                   		
 ntv       	= nv*2
 mzdata		= subset(newtwins, ZYG2019==1,selVars) 
 dzdata		= subset(newtwins, ZYG2019==2,selVars) 

 psych::describe(mzdata)
 psych::describe(dzdata)

 nv     	= length(selVars)/2     
 ntv    	= nv*2    
 thVals 	= 0.5
 aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
 cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
 eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 


 multi_ace = mxModel("ACE",
  mxModel("top",
  mxMatrix( name="Means", 		type="Full", nrow=1, ncol=nv, free=c(T,T,T), labels=c("m1","m2", "m3")), 
  mxAlgebra(name="expMean", 		cbind(Means,Means)), # Mean matrices for twin 1 and twin 2
  mxMatrix( name="inc",    		type="Lower", nrow=1,  ncol=1,  free=F, values=1),
  mxMatrix( name="A",			type="Symm",  nrow=nv, ncol=nv, free=T, values=0, 	labels=aLabs, lbound= -20, ubound= 20),
  mxMatrix( name="C",			type="Symm",  nrow=nv, ncol=nv, free=T, values = 0, labels=cLabs, lbound= -25, ubound=20),
  mxMatrix( name="E",			type="Symm",  nrow=nv, ncol=nv, free=T, 			labels=eLabs, lbound= -20, ubound=20),
  mxAlgebra(name="expCovMZ",		expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
  mxAlgebra(name="expCovDZ",		expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
  mxAlgebra(name="corrP",		expression=cov2cor(A+C+E)),
  mxAlgebra(name="corrA",		expression=cov2cor(A)),
  mxAlgebra(name="corrE",		expression=cov2cor(E)),
    mxCI(c("VC[1,10]","VC[2,11]","VC[3,12]",
           "VC[1,13]","VC[2,14]","VC[3,15]",
           "VC[1,16]","VC[2,17]","VC[3,18]")),
  mxAlgebra(name="VC", expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)) )),
  mxModel("MZ", mxData(mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML() ),
  mxModel("DZ", mxData(dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML() ),
  mxFitFunctionMultigroup(c("MZ","DZ")))
 omxGetParameters(multi_ace)
 #mxCheckIdentification(multi_ace)

 # Run ULS for better start values
  uls 				<- mxAutoStart( multi_ace )
  uls_fit			<- mxTryHard( uls )
  summary( uls_fit )

  	# uls_fit$top$algebras$expCovDZ
  	# uls_fit$top$algebras$expCovMZ
  	# round(multi_ace_fit$top$algebras$VC$result,2)

 # ACE 
  multi_ace_fit 				<- mxTryHard(uls_fit)
  summary(omxRunCI(multi_ace_fit, optimizer = "SLSQP"),verbose=F) 
  summary( multi_ace_fit ) 
  print(summary( multi_ace_fit )$CI)

 # AE model
  multi_ae 					<- mxModel( multi_ace_fit,name="AE")
  multi_ae 					<- omxSetParameters( multi_ae, label=cLabs,free=F,values= 0)
  multi_ae_fit				<- mxTryHard( multi_ae )
  summary( multi_ae_fit ) 

 # CE model
  multi_ce 					<- mxModel( multi_ace_fit,name="CE")
  multi_ce 					<- omxSetParameters( multi_ce, label=aLabs,free=F,values= 0)
  multi_ce_fit				<- mxTryHard( multi_ce )
  summary( multi_ce_fit ) 

 # E model
  multi_e 	 				<- mxModel( multi_ace_fit,name="E")
  multi_e 	 				<- omxSetParameters( multi_e, label=aLabs,free=F,values= 0)
  multi_e 					<- omxSetParameters( multi_e, label=cLabs,free=F,values= 0)
  multi_e_fit				<- mxTryHard( multi_e )
  summary( multi_e_fit ) 

 # Compare models
  mxCompare(multi_ace_fit,c(multi_ae_fit, multi_ce_fit, multi_e_fit)) 

  #   base comparison ep minus2LL   df      AIC    diffLL diffdf            p
  # 1  ACE       <NA> 21 7565.943 2890 7607.943        NA     NA           NA
  # 2  ACE         AE 15 7569.875 2896 7599.875  3.931924      6 6.858886e-01
  # 3  ACE         CE 15 7584.553 2896 7614.553 18.609495      6 4.876575e-03
  # 4  ACE          E  9 7648.184 2902 7666.184 82.240515     12 1.540054e-12





 # Common pathway - VETSA data only

 selVars    = c("SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
                "SM_V1_T2", "SM_V2_T2", "SM_V3_T2")

 mzdata		<- subset(newtwins, ZYG2019==1, selVars); # dim(mzdata) # 475  10
 dzdata		<- subset(newtwins, ZYG2019==2, selVars); # dim(dzdata) # 335  10

 nv     	= length(selVars)/2     
 ntv    	= nv*2    


 nVariables	= 3
 nFactors	= 1
 #thVals 	= 0.5
 #nth 		= 1
 #svLTh  	= 1.1    # start value for first threshold
 #svITh  	= 1.0    # start value for increments
 #svTh   	= matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv)     # start value for thresholds
 #round(auto_ACE_fit$matrices$thresh$values,2)
 #Th_val		= c(1.2,1.1,1.1)
 #Th_lb   	= matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)     # lower bounds for thresholds
 #Th_lab  	= c(paste("th",1:nth,"v1",sep=""),paste("th",1:nth,"v2",sep=""),paste("th",1:nth,"v3",sep=""))
 psi_lab_a	= c("ia11")
 psi_lab_c	= c("ic11")
 psi_lab_e	= c("ie11")
 #round(t(diag2vec(auto_ACE_fit$matrices$psi_c$values)),2)
 psi_a_val	= 0.7
 psi_c_val	= 0.2
 psi_e_val	= 0.2
 res_lab_a	= c("res_a1","res_a2","res_a3")
 res_lab_c	= c("res_c1","res_c2","res_c3")
 res_lab_e	= c("res_e1","res_e2","res_e3")
 loadS 		= 0.8
 loadF 		= T
 lamba_lab	= c("f11","f21","f31")

 cp_ACE = mxModel("ACE",
  mxModel("top",
  mxMatrix(name="Mean", 		type = "Full", nrow=1, ncol=nv, free = T, labels = c("m1","m2","m3"), lbound = -5, ubound = 5 ),
  mxAlgebra(name="expMean", 	expression = cbind(Mean, Mean)),  
  mxMatrix(name="psi_a", 	type = "Symm", nrow = nFactors, ncol = nFactors, free = T, labels = psi_lab_a, values = 0, lbound = -5, ubound = 5),
  mxMatrix(name="psi_c", 	type = "Symm", nrow = nFactors, ncol = nFactors, free = T, labels = psi_lab_c, values = 0, lbound = -5, ubound = 5),
  mxMatrix(name="psi_e", 	type = "Symm", nrow = nFactors, ncol = nFactors, free = T, labels = psi_lab_e, values = 1, lbound = -5, ubound = 5),
  mxMatrix(name="lamba",		type = "Full", nrow = nv, ncol = nFactors, free = T, labels= lamba_lab, values = 1, lbound = 0, ubound = 20),
  mxMatrix(name="epsilon_a", type = "Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_a, values =   0, lbound = -5, ubound = 5), 
  mxMatrix(name="epsilon_c", type = "Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_c, values =   0, lbound = -5, ubound = 5), 
  mxMatrix(name="epsilon_e", type = "Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_e, values = 0.1, lbound = -5, ubound = 5), 
  mxAlgebra(name="A", 		expression = lamba %&% psi_a + epsilon_a),
  mxAlgebra(name="C", 		expression = lamba %&% psi_c + epsilon_c),
  mxAlgebra(name="E", 		expression = lamba %&% psi_e + epsilon_e),  
  mxAlgebra(name="expCovMZ", expression = rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
  mxAlgebra(name="expCovDZ", expression = rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
  mxAlgebra(name="corrP",	expression = cov2cor(A+C+E)),
  mxAlgebra(name="corrA", 	expression = cov2cor(A)),
  mxAlgebra(name="corrE", 	expression = cov2cor(E)),
 # Standardization - constrain variance of CP = 1
 mxMatrix(name="unitM",type="Unit", nrow=nFactors, ncol=nFactors),
 mxConstraint(name="ConVar",expression=diag2vec(psi_a+psi_c+psi_e)==unitM),  
 mxCI(c("psi_a","psi_e", "ia11","ie11","f11","f21","f31")), 
 #"ia11","ie11", "res_a1","res_a2","res_a3","res_a4","res_e1","res_e2","res_e3","res_e4", "f11","f21","f31","f41",
 # "VC[1,13]","VC[2,14]","VC[3,15]","VC[4,16]", "VC[1,21]","VC[2,22]","VC[3,23]","VC[4,24]","corrA[2,1]"
 mxAlgebra(name="VC",expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))), 
 # Multiple groups
 mxModel("MZ", mxData(mzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),   
 mxModel("DZ", mxData(dzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),   
 mxFitFunctionMultigroup(c("MZ","DZ")))
 omxGetParameters(cp_ACE)
 #mxCheckIdentification(cp_ACE) # = GOOD & identified

 # ACE 
  cp_ACE_fit 		<- mxTryHard(cp_ACE)
  summary( cp_ACE_fit )  
  cp_ACE_fit
  cp_ACE_fit$top$algebras$VC$result 
  cp_ACE_fit_CIs <- summary(omxRunCI(cp_1F_fit,optimizer = "SLSQP"),verbose=T)$CI	
  cp_ACE_fit_CIs

 # AE model
 cp_1F_AE 			<- mxModel( cp_ACE_fit,name="AE")
 cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label=psi_lab_c,free=F,values= 0)
 cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label=res_lab_c,free=F,values= 0) 
 cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label=c("f11","f21","f31"),	free=T,values=0.5, ubound=10)  
 #cp_1F_AE	 		<- omxSetParameters( cp_1F_AE, label=c(res_lab_a,res_lab_e),		free=T,values=0.5, ubound=20)    
 #cp_1F_AE 			<- omxSetParameters( cp_1F_AE, label="ia11",free=T,values=20,ubound=40) 
 cp_1F_AE_fit		<- mxTryHard( 		 cp_1F_AE )
 summary( cp_1F_AE_fit ) 

 # CE model
 cp_1F_CE 			<- mxModel( cp_ACE_fit,name="CE")
 cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label=psi_lab_a,free=F,values= 0)
 cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label=res_lab_a,free=F,values= 0) 
 cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label=c("f11","f21","f31"),	free=T,values=0.5, ubound=10)  
 #cp_1F_CE	 		<- omxSetParameters( cp_1F_CE, label=c(res_lab_c,res_lab_e),		free=T,values=0.5, ubound=20)    
 #cp_1F_CE 			<- omxSetParameters( cp_1F_CE, label="ia11",free=T,values=20,ubound=40) 
 cp_1F_CE_fit		<- mxTryHard( 		 cp_1F_CE )
 summary( cp_1F_CE_fit ) 

 # E model
 cp_1F_E 			<- mxModel( cp_ACE_fit,name="E")
 cp_1F_E 			<- omxSetParameters( cp_1F_E, label=psi_lab_a,free=F,values= 0)
 cp_1F_E 			<- omxSetParameters( cp_1F_E, label=res_lab_a,free=F,values= 0) 
 cp_1F_E 			<- omxSetParameters( cp_1F_E, label=psi_lab_c,free=F,values= 0)
 cp_1F_E 			<- omxSetParameters( cp_1F_E, label=res_lab_c,free=F,values= 0) 
 cp_1F_E 			<- omxSetParameters( cp_1F_E, label=c("f11","f21","f31"),	free=T,values=0.5, ubound=10)  
 cp_1F_E_fit			<- mxTryHard( 		 cp_1F_E )
 summary( cp_1F_E_fit ) 

 # Compare models
  mxCompare(cp_ACE_fit,c(cp_1F_AE_fit, cp_1F_CE_fit, cp_1F_E_fit)) 
  #   base comparison ep minus2LL   df      AIC   diffLL diffdf            p
  # 1  ACE       <NA> 18 7567.589 2894 7603.589       NA     NA           NA
  # 2  ACE         AE 14 7570.721 2898 7598.721  3.13159      4 5.360506e-01
  # 3  ACE         CE 14 7585.386 2898 7613.386 17.79756      4 1.351729e-03
  # 4  ACE          E 10 7648.184 2902 7668.184 80.59501      8 3.710330e-14






 # Independent Pathways - VETSA data only

 selVars    = c("SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
                "SM_V1_T2", "SM_V2_T2", "SM_V3_T2")

 mzdata		<- subset(df2, ZYG2019==1, selVars); # dim(mzdata) # 475  10
 dzdata		<- subset(df2, ZYG2019==2, selVars); # dim(dzdata) # 335  10

 nv     	= length(selVars)/2     
 ntv    	= nv*2    
 nVariables	= 4
 nFactors	= 1
 #thVals 	= 0.5
 #nth 		= 1
 #svLTh  	= 1.1    # start value for first threshold
 #svITh  	= 1.0    # start value for increments
 #svTh   	= matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv)     # start value for thresholds
 #round(auto_ACE_fit$matrices$thresh$values,2)
 #Th_val		= c(1.2,1.1,1.1)
 #Th_lb   	= matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv)     # lower bounds for thresholds
 #Th_lab  	= c(paste("th",1:nth,"v1",sep=""),paste("th",1:nth,"v2",sep=""),paste("th",1:nth,"v3",sep=""))
 #round(t(diag2vec(auto_ACE_fit$matrices$psi_c$values)),2)
 psi_a_val	= 0.7
 psi_c_val	= 0.2
 psi_e_val	= 0.2
 res_lab_a	= c("res_a1","res_a2","res_a3")
 res_lab_c	= c("res_c1","res_c2","res_c3")
 res_lab_e	= c("res_e1","res_e2","res_e3")
 loadS 		= 0.8
 loadF 		= T
 lamba_lab_a	= c("f11_a","f21_a","f31_a")
 lamba_lab_c	= c("f11_c","f21_c","f31_c") 
 lamba_lab_e	= c("f11_e","f21_e","f31_e") 

 ip_ACE = mxModel("ACE",
  mxModel("top",
  mxMatrix(name="Mean", 		type="Full", nrow=1, ncol=nv, free = T, labels = c("m1","m2","m3"), lbound = -15, ubound = 5 ),
  mxAlgebra(name="expMean", 	expression= cbind(Mean, Mean)), 
  mxMatrix(name="lamba_a",	type="Full", nrow = nv, ncol = 1,  free = T, labels = lamba_lab_a, 	values = 0.3, lbound = -10, ubound = 15), #   
  mxMatrix(name="lamba_c",	type="Full", nrow = nv, ncol = 1,  free = T, labels = lamba_lab_c, 	values = 0.3, lbound = -10, ubound = 15), #   
  mxMatrix(name="lamba_e",	type="Full", nrow = nv, ncol = 1,  free = T, labels = lamba_lab_e, 	values = 0.3, lbound = -10, ubound = 15), # 
  mxMatrix(name="epsilon_a", type="Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_a, 	values = 0.3, lbound = -10, ubound = 15), #  
  mxMatrix(name="epsilon_c", type="Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_c, 	values = 0.3, lbound = -10, ubound = 15), #  
  mxMatrix(name="epsilon_e", type="Diag", nrow = nv, ncol = nv, free = T, labels = res_lab_e, 	values = 0.3, lbound = -10, ubound = 15), #  
  mxAlgebra(name="A", 		expression= lamba_a %*% t(lamba_a) + epsilon_a %*% t(epsilon_a)),
  mxAlgebra(name="C", 		expression= lamba_c %*% t(lamba_c) + epsilon_c %*% t(epsilon_c)),
  mxAlgebra(name="E", 		expression= lamba_e %*% t(lamba_e) + epsilon_e %*% t(epsilon_e)),  
  mxAlgebra(name="expCovMZ", expression= rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
  mxAlgebra(name="expCovDZ", expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
  mxAlgebra(name="VC",expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)))),
  mxModel("MZ",mxData(mzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
  mxModel("DZ",mxData(dzdata, type ="raw"), mxExpectationNormal(	covariance="top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
  mxFitFunctionMultigroup(c("MZ","DZ")))
  omxGetParameters(ip_ACE)
  #mxCheckIdentification(ip_ACE)

 # ACE 
 ip_ACE_fit 				<- mxTryHard( ip_ACE )
 summary( ip_ACE_fit )  




 # Autoregression - VETSA only

 selVars    = c("SM_V1_T1", "SM_V2_T1", "SM_V3_T1",
                "SM_V1_T2", "SM_V2_T2", "SM_V3_T2")

 mzdata			= subset(newtwins, ZYG2019==1, selVars); # dim(mzdata) # 475  10
 dzdata			= subset(newtwins, ZYG2019==2, selVars); # dim(dzdata) # 335  10
 nv        	= 3                  		
 ntv       	= nv*2

 # mxOption(NULL,"mvnRelEps",0.0045)
 nv     		= length(selVars)/2     
 ntv    		= nv*2    
 nVariables	= nv
 nFactors		= nv
 psi_lab_a	= c("ia11","ia22","ia33")
 psi_lab_c	= c("ic11","ic22","ic33")
 psi_lab_e	= c("ie11","ie22","ie33")
 res_lab_e	= "res_e"
 betaF   	= matrix(c(F,F,F, 
                     T,F,F,  
                     F,T,F),nv,byrow = TRUE)
 b_lab 		= matrix(c(NA, NA, NA,
                    "b", NA, NA,
                     NA,"b", NA),nv,byrow = TRUE)
 loadS 		= diag(nFactors)
 loadF 		= F
 psi_a_val	= c( 3, 5, 3)
 psi_c_val	= c(-1,-1,-1)
 psi_e_val	= c( 4,-1,-1)
 epsi_e_val	= 0.8

 auto_ACE = mxModel("ACE",
 mxModel("top",
 # Means & definition variables
 mxMatrix( name="Mean", 		type="Full", nrow=1, ncol=nv, free = c(T,T,T), labels = c("m1","m2","m3"), values=c(0.1,0.1,0.1), lbound = -5, ubound = 5 ),
 mxAlgebra(name="expMean", expression= cbind(Mean,Mean)), 
 # Psi - innovations
 mxMatrix(name="psi_a", 		type="Diag", nrow = 3, ncol = 3, free = T, labels = psi_lab_a, lbound = -5, ubound = 5 ), 
 mxMatrix(name="psi_c", 		type="Diag", nrow = 3, ncol = 3, free = T, labels = psi_lab_c, lbound = -5, ubound = 5 ),  
 mxMatrix(name="psi_e", 		type="Diag", nrow = 3, ncol = 3, free = T, labels = psi_lab_e, lbound = -5, ubound = 5 ), 
 # Epsilon - errors                                                                                   
 mxMatrix(name="epsilon_e", type="Diag", nrow = nVariables, ncol = nVariables, free = T, labels = res_lab_e, values = 0.5, lbound = -5, ubound = 5 ), 
 # Beta - causal pathways
 mxMatrix(name="beta", 		type="Full", nrow = nFactors, ncol = nFactors, free = betaF, labels = b_lab, lbound = -5, ubound = 5 ), 
 mxMatrix(name="I",				type="Iden", nrow = nFactors, ncol = nFactors),
 # Lamba                    	
 mxMatrix(name="lamba",		type="Full", nrow = nv, ncol = nv, free = F, values = diag(nv)), 
 # Expected variance & covariance
 mxAlgebra(name="A",				expression = (lamba %&% solve(I-beta)) %&% psi_a ), 
 mxAlgebra(name="C",				expression = (lamba %&% solve(I-beta)) %&% psi_c ), 
 mxAlgebra(name="E",				expression = (lamba %&% solve(I-beta)) %&% psi_e + epsilon_e),
 mxAlgebra(name="expCovMZ",expression= rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E))),
 mxAlgebra(name="expCovDZ",expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E))),
 # Standardization
 mxAlgebra(name="corrP", 	expression= cov2cor(A+C+E)),
 mxAlgebra(name="corrA", 	expression= cov2cor(A)),
 mxAlgebra(name="corrE", 	expression= cov2cor(E)), 
 mxAlgebra(name="VC",     expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)) )),
 # Groups
 mxModel("MZ", mxData(			mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()), 
 mxModel("DZ", mxData(			dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML()),
 mxFitFunctionMultigroup(c("MZ","DZ")))
 omxGetParameters(auto_ACE)
 #mxCheckIdentification( auto_ACE_fit )

  # options(max.print = 3000)
  # auto_uls 				<- mxAutoStart(auto_PBA)
  # summary(auto_fit 	<- mxTryHard(auto_uls))
  # auto_fit$top$algebras$VC$result
  # 
  # auto_fit$top$algebras$corrP$result
  # auto_fit$top$algebras$corrA$result
  # auto_fitci<-summary(omxRunCI(auto_fit, optimizer = "SLSQP"),verbose=T)

 # ACE 
  auto_ACE_fit 				<- mxTryHard( auto_ACE )
  summary( auto_ACE_fit )  
  auto_ACE_fit$top$corrA
  
  auto_ACE_fit$top$VC
  
  
  ############################################################
  #ACE MODELS FOR PHENOTYPIC GENETIC AND ENVIRON CORRELATIONS 
  ############################################################

  selVars    = c("MEMORY_V2_T1", "SM_38_T1",
                 "MEMORY_V2_T2","SM_38_T2")
  nv        	= 2                   		
  ntv       	= nv*2
  mzdata		= subset(newtwins, ZYG2019==1,selVars) 
  dzdata		= subset(newtwins, ZYG2019==2,selVars) 
  
  psych::describe(mzdata)
  psych::describe(dzdata)
  
  nv     	= length(selVars)/2     
  ntv    	= nv*2    
  thVals 	= 0.5
  aLabs  	= paste("a",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
  cLabs  	= paste("c",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")
  eLabs  	= paste("e",rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") 
  
  
  lbTh   <- matrix(rep(c(-3,(rep(0.001,nth-1))),nvo),nrow=nth,ncol=nvo)   # lower bounds for thresholds
  ubTh   <- matrix(rep(c(5,(rep(3,nth-1))),nvo),nrow=nth,ncol=nvo)     	# lower bounds for thresholds
  labTh  <- c(paste("v2","t",1:nth,sep=""))
  
  umx::umx_set_optimizer(opt="SLSQP")
  mxOption(NULL,"mvnRelEps",0.0055)
  mxOption(NULL, 'Number of Threads', parallel::detectCores())
  
  bivariate = mxModel("ACE_E",
                      mxModel("top",
                              mxMatrix(	name="Means", 		type="Full", nrow=1, ncol=nv, free=c(T,T), labels=c("m1","m2")), # Mean matrix for one twin. The 2nd variable is ordinal, mean of zero.
                              mxAlgebra(name="expMean", 	cbind(Means,Means)), # Mean matrices for twin 1 and twin 2
                              mxMatrix( name="inc",    		type="Lower", nrow=1, ncol=1, free=F, values=1),
                              mxMatrix(	name="A",			type="Symm", nrow=nv, ncol=nv, free=T, values=0, labels=aLabs, lbound= -20, ubound= 20), 		# Constrains element [2,2] to 1
                              mxMatrix(	name="C",			type="Symm", nrow=nv, ncol=nv, free=F, values = 0, labels=cLabs, lbound= -20, ubound=20), 		# Constrains element [2,2] to 1
                              mxMatrix(	name="E",			type="Symm", nrow=nv, ncol=nv, free=T, labels=eLabs,  lbound= -20, ubound=20),
                              mxAlgebra(name="expCovMZ",	expression= rbind( cbind(A+C+E, A+C  ), cbind(  A+C, A+C+E))),
                              mxAlgebra(name="expCovDZ",	expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind( 0.5%x%A+C, A+C+E) )),
                              mxAlgebra(	name="corrP",	expression=cov2cor(A+C+E)),
                              mxAlgebra(	name="corrA",	expression=cov2cor(A)),
                              mxAlgebra(	name="corrE",	expression=cov2cor(E)),
                               mxCI(c("VC[1,7]","VC[2,12]", "corrP[1,2]", "corrA[1,2]", "corrE[1,2]")),
                              mxAlgebra(name="VC", expression=cbind(A,C,E,A/(A+C+E),C/(A+C+E),E/(A+C+E)), dimnames=list(rep('VC',nv),rep(c('A','C','E','SA','SC','SE'),each=nv)) ) ),
                      
                      mxModel("MZ", mxData(mzdata, type ="raw"), mxExpectationNormal("top.expCovMZ", means="top.expMean", dimnames=selVars), mxFitFunctionML() ),
                      mxModel("DZ", mxData(dzdata, type ="raw"), mxExpectationNormal("top.expCovDZ", means="top.expMean", dimnames=selVars), mxFitFunctionML() ),
                      mxFitFunctionMultigroup(c("MZ","DZ")))
  
  # mxCheckIdentification(multi_PBA)
  
  # Run unweighted least squares to find better start values
  uls <- mxAutoStart(bivariate)
  summary(uls_fit_emp 	<- mxTryHard(uls, intervals=T))
  uls_fit_emp$top$corrE
  uls_fit_emp$top$corrP
  uls_fit_emp$top$corrA
  uls_fit_emp$top$VC
  #uls_fit_emp_AE$top$algebras$expCovDZ
  #summary(uls_fit_emp_E)
  #uls_fit_emp_AE$top$algebras$expCovMZ
  #uls_fit_empci <-summary(omxRunCI(uls_fit_emp, optimizer = "SLSQP"),verbose=T); uls_fit_empci
  #uls_fit_empci$top
  uls_fitci<-summary(omxRunCI(uls_fit_emp, optimizer = "SLSQP"),verbose=T)
  uls_fitci
  # AE model
  uls_fit_emp_AE 					<- mxModel( uls_fit_emp,name="AE")
  uls_fit_emp_AE 					<- omxSetParameters( uls_fit_emp_AE, label=cLabs,free=F,values= 0)
  summary( uls_fit_emp_AE_fit		<- mxTryHard(uls_fit_emp_AE, extraTries=50, greenOK=TRUE,checkHess=FALSE,intervals=F))
  summary( uls_fit_emp_AE_fit ) 
  uls_fit_emp_AE_fit$top$algebras$VC$result
  uls_fitci<-summary(omxRunCI(uls_fit_emp_AE_fit, optimizer = "SLSQP"),verbose=T)
  uls_fitci
  
  ############################################################
  #POLYCHOR CORRELATIONS ACCOUNTING FOR TWINNESS AND ZYGOSITY
  ############################################################
  
  #POLYCHOR WITH CONTINUOUS VARIABLES
  
  selVarsMZ	= c("SM_V3_T1","KUNKLES1_RES_T1",
                "SM_V3_T2")
  selVarsDZ	= c("SM_V3_T1","KUNKLES1_RES_T1",
                "SM_V3_T2","KUNKLES1_RES_T2")
  
  df$NEUROPGRS_3_T1 = ifelse(df$ZYG2019==1 & is.na(df$NEUROPGRS_3_T1),df$NEUROPGRS_3_T2,df$NEUROPGRS_3_T1)
  
  mzdata 	= subset(df2, ZYG2019==1, selVarsMZ)
  dzdata 	= subset(df2, ZYG2019==2, selVarsDZ)
  
  describe(mzdata)
  describe(dzdata)
  
  nv     	= length(selVarsMZ)/2     
  ntv    	= length(selVarsMZ)
  
  
  labelsMZ =c(
    "a", 	
    "c", "b", 	
    "d", "c", "a")
  
  labelsDZ =c(
    "a", 	
    "c", "b", 	
    "d", "f", "a", 
    "f", "e", "c", "b")
  
  
  freeMZ =
    c(
      T,
      T,T,	
      T,T,T)
  
  freeDZ =
    c(
      T,
      T,T,	
      T,T,T,	
      T,T,T,T)
  
  valsMZ =
    c(
      5,
      0.5,5,	
      0.5,0.5,5)
  
  valsDZ =
    c(
      5,
      0.5,5,	
      0.5,0.5,5,	
      0.5,0.5,0.5,5)
  
  temp = mxModel("corr",
                 mxModel("top",
                         mxMatrix( name="MeansMZ", 	type="Full", nrow=1, ncol=3, free=T, labels=c("u1","u2","u1")),
                         mxMatrix( name="MeansDZ", 	type="Full", nrow=1, ncol=4, free=T, labels=c("u1","u2","u1","u2")),  
                         mxMatrix( name="CovMZ",		type="Symm", nrow=3, ncol=3, free=freeMZ, labels=labelsMZ, values=valsMZ,byrow=T, lbound=-20,ubound=50),
                         mxMatrix( name="CovDZ",		type="Symm", nrow=4, ncol=4, free=freeDZ, labels=labelsDZ, values=valsDZ,byrow=T, lbound=-20,ubound=50),
                         mxAlgebra(name="corrMZ",	expression=cov2cor(CovMZ)),
                         mxAlgebra(name="corrDZ",	expression=cov2cor(CovDZ)),
                         mxCI(c("corrDZ[4,3]")) ),    
                 mxModel("MZ", mxData(		mzdata, type ="raw"),mxExpectationNormal("top.CovMZ", means="top.MeansMZ", dimnames=selVarsMZ), mxFitFunctionML() ),
                 mxModel("DZ", mxData(		dzdata, type ="raw"),mxExpectationNormal("top.CovDZ", means="top.MeansDZ", dimnames=selVarsDZ), mxFitFunctionML() ),
                 mxFitFunctionMultigroup(c("MZ","DZ")) )
  omxGetParameters(temp)
  
  SM_38_fit = mxTryHard(temp)   
  SM_38_fit = mxRun(SM_38_fit,intervals=T)     
  summary( SM_38_fit)
  
  SM_38_fit$top$CovMZ
  SM_38_fit$top$CovDZ  
  SM_38_fit$top$corrMZ
  SM_38_fit$top$corrDZ
  
  SM_38 = paste0(
    format(round(print(summary(SM_38_fit)$CI)[1,2],2),nsmall=2)," (",
    format(round(print(summary(SM_38_fit)$CI)[1,1],2),nsmall=2),", ",
    format(round(print(summary(SM_38_fit)$CI)[1,3],2),nsmall=2),")")
  SM_38
  
  #POLYCHOR WITH ORDINAL VARIABLES
  
  
  
  #Ordinal Variable Input 
  MZvarsO		<- c("APOE4VSNON_T1")
  MZnth       	<- 1             				 
  MZnvo       	<- 1                       
  MZlabTh  <- c(paste("t",1:MZnth,sep=""))
  mxFactor(df2[MZvarsO], levels=c(0,1))
  
  DZvarsO		<- c("APOE4VSNON_T1","APOE4VSNON_T2")
  DZnth       	<- 1             				 
  DZnvo       	<- 2                       
  DZlabTh  <- c(paste("t",1:DZnth,sep=""))
  df2[DZvarsO]<- mxFactor(df2[DZvarsO], levels=c(0,1)) 
  
  # SM_V3 & NEUROPGRS
  
  selVarsMZ	= c("SM_V3_T1","APOE4VSNON_T1",
                "SM_V3_T2")
  selVarsDZ	= c("SM_V3_T1","APOE4VSNON_T1",
                "SM_V3_T2","APOE4VSNON_T2")
  
  #df$NEUROPGRS_3_T1 = ifelse(df$ZYG2019==1 & is.na(df$NEUROPGRS_3_T1),df$NEUROPGRS_3_T2,df$NEUROPGRS_3_T1)
  
  mzdata 	= subset(df2, ZYG2019==1, selVarsMZ)
  dzdata 	= subset(df2, ZYG2019==2, selVarsDZ)
  
  describe(mzdata)
  describe(dzdata)
  
  nv     	= length(selVarsMZ)/2     
  ntv    	= length(selVarsMZ)
  
  labelsMZ =c(
    "a", 	
    "c", "b", 	
    "d", "c", "a")
  
  labelsDZ =c(
    "a", 	
    "c", "b", 	
    "d", "f", "a", 
    "f", "e", "c", "b")
  
  freeMZ =
    c(
      T,
      T,T,	
      T,T,T)
  
  freeDZ =
    c(
      T,
      T,T,	
      T,T,T,	
      T,T,T,T)
  
  valsMZ =
    c(
      5,
      0.5,5,	
      0.5,0.5,5)
  
  valsDZ =
    c(
      5,
      0.5,5,	
      0.5,0.5,5,	
      0.5,0.5,0.5,5)
  
  temp = mxModel("corr",
                 mxModel("top",
                         mxMatrix( name="MeansMZ", 	type="Full", nrow=1, ncol=3, free=T, labels=c("u1","u2","u1")),
                         mxMatrix( name="MeansDZ", 	type="Full", nrow=1, ncol=4, free=T, labels=c("u1","u2","u1","u2")),  
                         mxMatrix( name="MZth",			type="Full",  nrow=1, ncol=1, free=T, labels=MZlabTh), 	# Threshold matrix for one twin
                         mxMatrix( name="DZth",			type="Full",  nrow=1, ncol=2, free=T, labels=DZlabTh), 	# Threshold matrix for one twin
                         mxMatrix( name="CovMZ",		type="Symm", nrow=3, ncol=3, free=freeMZ, labels=labelsMZ, values=valsMZ,byrow=T, lbound=-100,ubound=100),
                         mxMatrix( name="CovDZ",		type="Symm", nrow=4, ncol=4, free=freeDZ, labels=labelsDZ, values=valsDZ,byrow=T, lbound=-100,ubound=100),
                         mxAlgebra(name="corrMZ",	expression=cov2cor(CovMZ)),
                         mxAlgebra(name="corrDZ",	expression=cov2cor(CovDZ)),
                         mxCI(c("corrDZ[4,3]")) ),    
                 mxModel("MZ", mxData(		mzdata, type ="raw"),mxExpectationNormal("top.CovMZ", means="top.MeansMZ", dimnames=selVarsMZ, thresholds="top.MZth", threshnames=MZvarsO), mxFitFunctionML() ),
                 mxModel("DZ", mxData(		dzdata, type ="raw"),mxExpectationNormal("top.CovDZ", means="top.MeansDZ", dimnames=selVarsDZ, thresholds="top.DZth", threshnames=DZvarsO), mxFitFunctionML() ),
                 mxFitFunctionMultigroup(c("MZ","DZ")) )
  omxGetParameters(temp)
  uls <- mxAutoStart(temp)
  
  SM_V3_fit = mxTryHard(uls)   
  SM_V3_fit = mxRun(SM_V3_fit,intervals=T)     
  summary( SM_V3_fit)
  
  SM_V3_fit$top$CovMZ
  SM_V3_fit$top$CovDZ  
  SM_V3_fit$top$corrMZ
  SM_V3_fit$top$corrDZ
  
  SM_V3 = paste0(
    format(round(print(summary(SM_V3_fit)$CI)[1,2],2),nsmall=2)," (",
    format(round(print(summary(SM_V3_fit)$CI)[1,1],2),nsmall=2),", ",
    format(round(print(summary(SM_V3_fit)$CI)[1,3],2),nsmall=2),")")
  SM_V3 	# "0.12 (0.04, 0.20)"
  
  