# Dissertation

#EVERYTHING 

library(survival)
library(simsurv)
library(gsDesign)

#Simulate Survival Data

#Function that simulates one arm of exponential survival data: control, treatment and historical data can be generated through this 
#by just inputting different argument values. The function just simulates survival data and gives it a label to say which arm it is from.

SimOneArmSurvivalData <- function(n, lambda, label){
  
  #n: number of observations to be simulated 
  #lambda: scale parameter of the exponential distribution we assume the data follows
  #label: the value to assign the observations to indicate which arm the data belongs to. 
  
  x = data.frame(id = 1:n) 
  
  SurvivalData <- as.data.frame(simsurv(dist = "exponential", lambdas = lambda, x=x)) 
  
  Label <- as.factor(rep(label,n)) 
  
  FinalData <- cbind(Arm = Label, SurvivalData) 
  
  return(FinalData)
  
}

#When data is simulated through this function, it has no maximum upper limit, and therefore no censored observations, 
#cos we aren't taking into account the end of the trial yet and drop out until later on. This output represents how long this person actually lives for. 


#We have event times when assuming all participants were recruited at the start of the trial, this is not the case, so we need to shift the event times
#To reflect that they've started at different points. 

#With knowledge of the start and end point of the recruitment/accrual period, randomly assign each observation a time of recruitment, 
#assume uniform distribution.

#Once we have their recruitment time, we want to find the calendar time of when we stopped observing that person, which will be the recruitment
#time plus the event time
#So now we want to see if their calendar time is valid. If the calendar time is longer than the full length of the study, because of their given 
#recruitment time, we set their actual event time for the largest it would end up being, which is the maximum time, minus their recruitment time
#If they did not reach the end of the study, keep the event time the way it is. 

#The event time to be used in the model will not be eventtime but will be the newly defined ActualEventTime. 
#And then recruitment time can be used to make a time variable in the power functions.  

#This means though that the historical data will also have to have a recruitment time. 

TESTSimOneArmSurvivalData <- SimOneArmSurvivalData(10, 0.01, 1)
TESTDataToBeChanged <- as.data.frame(cbind(TESTSimOneArmSurvivalData, recruitmenttime = runif(nrow(TESTSimOneArmSurvivalData), 0, 20)))
TESTDataToBeChanged

TESTCalendarTime <- TESTDataToBeChanged$eventtime + TESTDataToBeChanged$recruitmenttime
TESTDataToBeChanged2 <- as.data.frame(cbind(TESTDataToBeChanged, calendarTime = TESTCalendarTime))
TESTDataToBeChanged2


RecruitmentAdjustment <- function(SimulatedData, StartofRecruitment, EndofRecruitment, mfut){
  #SimulatedData should be of the form Arm, eventtime, status in order to work. 
  
  #Assign each patient a random (uniform distribution) recruitment time. 
  DataToBeChanged <- as.data.frame(cbind(SimulatedData, recruitmenttime = runif(nrow(SimulatedData), StartofRecruitment, EndofRecruitment)))
  #Use this to find the time the observation was taken in calendar time. 
  calendarTime <- DataToBeChanged$eventtime + DataToBeChanged$recruitmenttime
  DataToBeChanged2 <- as.data.frame(cbind(DataToBeChanged, calendarTime = calendarTime))
  
  ActualEventTime<- rep(NA, nrow(DataToBeChanged2)) #Create a variable that is the actual event time we will use. 
  
  for(i in 1:nrow(DataToBeChanged2)){
    if (DataToBeChanged2[i, ]$calendarTime >= mfut){ #If calendar time is more than the full length of the trial, event time needs to be changed
      ActualEventTime[i] = mfut - DataToBeChanged2[i, ]$recruitmenttime #They reached end of trial so event time will be time between recruitment and end 
      DataToBeChanged2[i, ]$status = rep(0, 1) #say its been censored.
    }
    else{ 
      ActualEventTime[i] = DataToBeChanged2[i, ]$eventtime #If calendar time is not full length of trial, ActualEventTime = eventtime. 
    }
    
  }
  return(as.data.frame(cbind(DataToBeChanged2, ActualEventTime)))
}



#Drop Out Function

#To be applied to survival data sets that do not yet take into account drop out during the trial, and returns a data set 
#That has some observations amended in a way that reflects a certain expected drop out rate throughout the trial. 

#Says that whether or not someone has dropped out is binomially distributed with a probability p of dropping out: Bin(p). 
#Using this, randomly decide whether each observation is from someone who has dropped out or not.

#For the ones who have, give them a new event time: eventtime x p, where p is a randomly generated number between 0 and 1.
#Will also need to say that the observation is now censored, and the event time is the time of their dropping out
#This results in event times that are correlated to the time that they dropped out, which is reasonable because the event time 
#a person is censored at when they drop out, is certainly related to when they have dropped out. 
#For example, someone would not be censored at 5 months when they didn't drop out until 10 months in. 

#When someone has not dropped out, their observation remains the same. 

DropOut <- function(DataWithoutDropOut, probdropout){
  
  #DataWithoutDropOut: A survival data set that does not yet take into account drop out 
  #probdropout: the expected drop out rate of the trial: eg 0.05 for 5% drop out rate
  
  DropoutIndicator <- rbinom(n = nrow(DataWithoutDropOut), size = 1, prob = probdropout) #Create dropout indicator
  DataWithIndicator <- as.data.frame(cbind(DataWithoutDropOut, DropoutIndicator)) #Add the indicator variable to select observations that have dropped out
  
  for(i in 1:nrow(DataWithIndicator)){
    if (DataWithIndicator[i, ]$DropoutIndicator == 1) { 
      p <- runif(1, 0, 1)
      DataWithIndicator[i, ]$ActualEventTime <- DataWithIndicator[i, ]$ActualEventTime*p #If they've dropped out give random event time less than the simulated one
      DataWithIndicator[i, ]$status <- rep(0, 1) #And consider it censored
    } else{
      DataWithIndicator[i, ]$ActualEventTime <- DataWithIndicator[i, ]$ActualEventTime #If they haven't dropped out it stays the same.
      DataWithIndicator[i, ]$status <- DataWithIndicator[i, ]$status
    }
  }
  return(DataWithIndicator)
}

#Use these three functions to simulate a data set to be used directly in power analysis, that has been adjusted for drop out and to
#Reflect participants being recruited at different times across a period. 

CompleteTrialSimulation <- function(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout){
  
  
  ConcurrentControlData <- SimOneArmSurvivalData(n1, lambda1, label1) #Simulate concurrent control survival data
  
  TreatmentData <- SimOneArmSurvivalData(n2, lambda2, label2) #Simulate treatment group survival data
  
  FullTrialData <- rbind(ConcurrentControlData, TreatmentData) #Create full survival data set of survival times  
  FullTrialData <- FullTrialData[ , -2] #Remove duplicate id column
  
  FullTrialDataRecruitmentAdj <- RecruitmentAdjustment(FullTrialData, StartofRecruitment, EndofRecruitment, mfut)
  
  FinalTrialData <- DropOut(FullTrialDataRecruitmentAdj, probdropout)
  
  return(FinalTrialData)
  
  
}


#Power Functions

#Power when using Separate Analysis

PowerSeparateAnalysis <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment) { 
  
  #For the use of CompleteTrialSimulation: 
  #n1 : number of participants in control group
  #n2 : number of participants in treatment group
  #lambda1 : exponential scale parameter for survival distribution of control group
  #lambda2 : exponential scale parameter for survival distribution of treatment group
  #mfut: max follow up time
  #label1: identifer to be given to control group
  #label2: identifier to be given to treatment group
  
  #For actual power calculation
  #M: Number of trials to be simulated
  #alpha : two sided significance level
  
  N <- rep(NA, M) #Creating a vector in which we can store indicator variable of whether or not the null hypothesis was rejected for each 
  #simulated trial
  
  for (i in 1:M) {
    #Simulate current data set
    CurrentData <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    ModelSeparate <- coxph(Surv(ActualEventTime, status) ~ Arm, data = CurrentData) #Use survival object to fit a cox proportional hazards model. 
    
    if (summary(ModelSeparate)$coefficients[,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if wasn't
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #Number of trials that rejected the null/total trials = power.
}



#Note that the functions using historical trial data only work if the TREATMENT of historical is given a value of 1 - as this is what becomes control 
#in the current trial, and needs to have a value of 1 to be recognized as a control observation. 


#Power when using Pooled Analysis to borrow historical data

#Simulate one historical data set to use in all current data set, then use as an argument in the function.

PowerPooledAnalysis <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, HistoricalTrial, probdropout, StartofRecruitment, EndofRecruitment) { 
  
  #To be used in CompleteTrialSimulation function
  #n1 : number of participants in control group
  #n2 : number of participants in treatment group
  #lambda1 : exponential scale parameter for survival distribution of control group
  #lambda2 : exponential scale parameter for survival distribution of treatment group
  #M: Number of trials to be simulated
  #mfut: max follow up time
  #label1: identifier to be given to control group
  #label2: identifier to be given to treatment group
  #HistoricalControlData: the data set to be incorporated, only the observations to be used in analysis should be in this data set
  #probdropout : drop out rate for current trial
  
  #For power calculation
  #M: Number of trials to be simulated
  #alpha : two sided significance level
  
  
  N <- rep(NA, M) #Creating a vector in which we can store indicator variable of whether or not the null hypothesis was rejected for each 
  #simulated trial
  
  for (i in 1:M) {
    
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout) #Simulate current trial 
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract only Treatment data from historical trial. 
    #This is because treatment in historical trial = control in current trial
    
    FinalTrialDataPooled <- rbind(CurrentTrial, HistoricalControlData) #Add historical data into data set for analysis
    
    ModelPooled <- coxph(Surv(ActualEventTime, status) ~ Arm, data = FinalTrialDataPooled) #Use survival object to fit a cox proportional hazards model. 
    
    if (summary(ModelPooled)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if wasn't
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #Number of trials that rejected the null/total trials = power. 
}





#Power when modelling the time as a step change in time. 

PowerStepTimeTrend <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment, 
                               HistoricalTrial) { 
  
  N <- rep(NA, M) #Creating a vector in which we can store indicator variable of whether or not the null hypothesis was rejected for each 
  #simulated trial
  
  for (i in 1:M) {
    
    #Simulate current data set
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract the treatment data from historical trial. 
    #This is because treatment in historical trial = control in current trial
    
    Current <- rep(1 ,nrow(CurrentTrial)) 
    Historical <- rep(0, nrow(HistoricalControlData)) 
    
    CurrentTrial2 <- as.data.frame(cbind(CurrentTrial, HistOrCurrent = Current))  #Assign current data value of 1
    HistoricalControlData2 <-as.data.frame(cbind(HistoricalControlData, HistOrCurrent = Historical)) #Assign historical data value of 0
    
    CombinedFinalTrialData <- rbind(CurrentTrial2, HistoricalControlData2) #Add historical data into data set for analysis
    
    #Fit model, dependent on arm and whether data is historical or current, which is used a step change in time
    ModelStep <- coxph(Surv(ActualEventTime, status) ~ Arm + HistOrCurrent , data = CombinedFinalTrialData)  
    
    if (summary(ModelStep)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if wasn't
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #Number of trials that rejected the null/total trials = power. 
}





#Power when modelling the time as a linear change in time 

PowerLinearTimeTrend <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment, HistoricalTrial,
                                 mfutHist, TimeBetweenTrials) { 
  
  N <- rep(NA, M) #Creating a vector in which we can store indicator variable of whether or not the null hypothesis was rejected for each 
  #simulated trial
  
  for (i in 1:M) {
    
    #Simulate current data set 
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract the treatment data from historical data set
    #This is because treatment in historical trial = control in current trial
    
    Current <- rep(1,nrow(CurrentTrial))
    Historical <- rep(0, nrow(HistoricalControlData))
    
    CurrentTrial2 <- as.data.frame(cbind(CurrentTrial, HistOrCurrent = Current)) #Assign a value of 1 to current data
    HistoricalControlData2 <-as.data.frame(cbind(HistoricalControlData, HistOrCurrent = Historical)) #Assign a value of 0 to historical data
    
    CombinedFinalTrialData <- rbind(CurrentTrial2, HistoricalControlData2) #Add historical data into data set for analysis
    
    RecruitmentFromBaseline <- rep(NA, n1+n2+nrow(HistoricalControlData2)) #Create vector to store recruitmentfrombaseline/create the variable
    
    for(j in 1:nrow(CombinedFinalTrialData)){
      
      if(CombinedFinalTrialData[j, ]$HistOrCurrent == 0){ #If the observation is historical 
        
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime#Baseline recruitment is just the recruitment time. 
      }
      else{ #When the observation is current
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime + mfutHist + TimeBetweenTrials #Baseline recruitment is recruitment time 
        #of current trial, plus time it has been since the start of the historical trial. 
      }
    }
    
    CombinedFinalTrialData2 <- cbind(CombinedFinalTrialData, RecruitmentFromBaseline)  #Add baseline recruitment time as a variable in the data set
    
    #Fit model dependent on arm and baseline recruitment time, which represents linear time
    ModelLinear <- coxph(Surv(ActualEventTime, status) ~ Arm + RecruitmentFromBaseline, data = CombinedFinalTrialData2) 
    
    if (summary(ModelLinear)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if wasn't
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #Number of trials that rejected the null/total trials = power. 
}



#Power when modelling the time as both linear and a step change 



PowerStepLinearTimeTrend <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment, HistoricalTrial,
                                     mfutHist, TimeBetweenTrials) { 
  
  N <- rep(NA, M) #Creating a vector in which we can store indicator variable of whether or not the null hypothesis was rejected for each 
  #simulated trial
  
  for (i in 1:M) {
    
    #Simulate current data set 
    
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract the treatment data from historical trial 
    #This is because treatment in historical trial = control in current trial
    
    Current <- rep(1,nrow(CurrentTrial))
    Historical <- rep(0, nrow(HistoricalControlData))
    
    CurrentTrial2 <- as.data.frame(cbind(CurrentTrial, HistOrCurrent = Current)) #Assign a value of 1 to all current observations
    HistoricalControlData2 <-as.data.frame(cbind(HistoricalControlData, HistOrCurrent = Historical)) #Assign a value of 0 to historical observations
    
    CombinedFinalTrialData <- rbind(CurrentTrial2, HistoricalControlData2) #Add historical data into data set for analysis
    
    RecruitmentFromBaseline <- rep(NA, n1+n2+nrow(HistoricalControlData2))
    
    for(j in 1:nrow(CombinedFinalTrialData)){
      
      if(CombinedFinalTrialData[j, ]$HistOrCurrent == 0){ #When data is historical
        
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime #baseline recruitment is just recruitment time
      }
      else{ #When data is current
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime + mfutHist + TimeBetweenTrials #Baseline recruitment is recruitment time 
        #of current trial, plus time it has been since the start of the historical trial
      }
    }
    
    CombinedFinalTrialData2 <- cbind(CombinedFinalTrialData, RecruitmentFromBaseline)
    
    #Fit model dependent on arm, steo time and linear time
    ModelStepLinear <- coxph(Surv(ActualEventTime, status) ~ Arm + RecruitmentFromBaseline + HistOrCurrent , data = CombinedFinalTrialData2) 
    
    
    if (summary(ModelStepLinear)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if wasn't
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #Number of trials that rejected the null/total trials = power. 
}


#TESTING 


#Historical Trial 

#Assume there was a previous clinical trial where the experimental group is the same treatment as the current control
#Assume they were looking for a difference of hazard ratio of 2, with lambda = 0.03 and 0.02 under the control and treatment arm respectively
#Choosing this because Myeloma trial is expecting lambda1 approx 0.021.  Want the control parameter 
#very similar to historical experimental arm parameter (same treatment), but a small change in parameter takes large sample
#size to detect and would make it computationally intensive, so make control in historical 0.03.
#Assume that the trial was half as long as Myeloma XI. 
#Using the sample size function, it tells us to use 328 participants in each arm, for witnessing 67 events.
#alpha, beta and drop out is same, as usually chosen from convention
nSurvival(lambda1 = 0.03 , lambda2 = 0.02 , Ts = 48 , Tr = 24, ratio = 1, alpha = 0.05 , beta = 0.2 , sided =2, eta = 0.05) 

set.seed(1234) #To obtain the same historical data set 
HistoricalTrial <- CompleteTrialSimulation(628, 0.03, 48, 0, 628, 0.02, 1, 0, 24, 0.05)
HistoricalTrial


#Current Trial Power Analysis

#The Myeloma trial non-intensive pathway needed 545 events observed, needing 787 participants 

perarmsample = 787/2

lambda1Myeloma = log(2)/33
lambda1Myeloma
lambda2Myeloma = log(2)/42
lambda2Myeloma

#Makes 394 in each arm

#lambda1 = 0.02100446
#Lambda2 = 0.0165035
#alpha = 0.05
#96 month maximum follow up
#48 month recruitment
#5% drop out

 
#Seperate Analysis 
set.seed(1234)
PowerSeparateAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48) #0.8

#For the methods that borrow historical data, we will be assuming they borrow all data available 
#from the previous trial, on the control treatment

#Pooled Analysis
set.seed(1234)
PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial, 0.05, 0, 48) #0.798

#Modelling time trend as a step 
set.seed(1234)
PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial) #0.801

#Modelling time trend as linear
set.seed(1234)
PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial, 70, 60) #0.819

#Modelling time trend as linear and step
set.seed(1234)
PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial, 70, 60) #0.805




#Type 1 errors. 
set.seed(1234)
Type1ErrorSeparate <- PowerSeparateAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48)

set.seed(1234)
Type1ErrorPooled <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial, 0.05, 0, 48)

set.seed(1234)
Type1ErrorTimeTrend <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial)

set.seed(1234)
Type1ErrorLinearTimeTrend <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial, 70, 60)

set.seed(1234)
Type1ErrorStepLinearTimeTrend <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial, 70, 60)

Type1ErrorSeparate #0.057
Type1ErrorPooled #0.068
Type1ErrorTimeTrend #0.057
Type1ErrorLinearTimeTrend #0.058
Type1ErrorStepLinearTimeTrend #0.055

#Sensitivity Analysis - HASN'T BEEN RUN

#Sensitivity Analysis when using pooling to borrow historical data. 
#Create a matrix of combinations of lambda and sample size we want to use for historical data. 

lambdaHistSet <- seq(0.01, 0.02, 0.005)
nHistSet <- seq(50, 60, 10)
nHistLambdaHistCombo1 <- cbind(expand.grid(nHist = nHistSet, LambdaHist = lambdaHistSet))
nHistLambdaHistCombo1
#Now we have a combination of lambdas and sample sizes we could use for historical borrowing. We want to test each combination of them. 

#For each row of the matrix of combinations, use the lambda and the sample size in that row to simulate a historical data set
#Then use this historical data set to calculate power when pooling 
#Obtain a vector of powers, can find which combination of lambda and n gives us highest and lowest power.


#CompleteTrialSimulation <- function(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)


SensitivityAnalysisPooled <- function(nControlHist, lambda1Hist, mfutHist, nHistLambdaHistCombo, StartofRecruitmentHist,EndofRecruitmentHist, probdropoutHist)
{
  
  N <- rep(NA, nrow(nHistLambdaHistCombo))
  
  for(i in 1:nrow(nHistLambdaHistCombo)){ 
    
    HistoricalTrialSensAn <- CompleteTrialSimulation(nControlHist, lambda1Hist, mfutHist, 0, nHistLambdaHistCombo[i,1], nHistLambdaHistCombo[i,2], 1, StartofRecruitmentHist, EndofRecruitmentHist, probdropoutHist)
    Power <- PowerPooledAnalysis(592, 592, 0.01050223, 0.008251752, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrialSensAn, 0.05, 0, 48)
    
    N[i] <- Power
  }
  
  return(N)}

SensitivityAnalysisStep <- function(nControlHist, lambda1Hist, mfutHist, nHistLambdaHistCombo, StartofRecruitmentHist,EndofRecruitmentHist, probdropoutHist)
{
  
  N <- rep(NA, nrow(nHistLambdaHistCombo))
  
  for(i in 1:nrow(nHistLambdaHistCombo)){ 
    
    HistoricalTrialSensAn <- CompleteTrialSimulation(nControlHist, lambda1Hist, mfutHist, 0, nHistLambdaHistCombo[i,1], nHistLambdaHistCombo[i,2], 1, StartofRecruitmentHist, EndofRecruitmentHist, probdropoutHist)
    Power <- PowerStepTimeTrend(592, 592, 0.01050223, 0.008251752, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialSensAn)
    
    N[i] <- Power
  }
  
  return(N)}

SensitivityAnalysisLinear <- function(nControlHist, lambda1Hist, mfutHist, nHistLambdaHistCombo, StartofRecruitmentHist,EndofRecruitmentHist, probdropoutHist)
{
  
  N <- rep(NA, nrow(nHistLambdaHistCombo))
  
  for(i in 1:nrow(nHistLambdaHistCombo)){ 
    
    HistoricalTrialSensAn <- CompleteTrialSimulation(nControlHist, lambda1Hist, mfutHist, 0, nHistLambdaHistCombo[i,1], nHistLambdaHistCombo[i,2], 1, StartofRecruitmentHist, EndofRecruitmentHist, probdropoutHist)
    Power <- PowerLinearTimeTrend(592, 592, 0.01050223, 0.008251752, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialSensAn, 70, 60)
    
    N[i] <- Power
  }
  
  return(N)}

SensitivityAnalysisLinearStep <- function(nControlHist, lambda1Hist, mfutHist, nHistLambdaHistCombo, StartofRecruitmentHist,EndofRecruitmentHist, probdropoutHist)
{
  
  N <- rep(NA, nrow(nHistLambdaHistCombo))
  
  for(i in 1:nrow(nHistLambdaHistCombo)){ 
    
    HistoricalTrialSensAn <- CompleteTrialSimulation(nControlHist, lambda1Hist, mfutHist, 0, nHistLambdaHistCombo[i,1], nHistLambdaHistCombo[i,2], 1, StartofRecruitmentHist, EndofRecruitmentHist, probdropoutHist)
    PowerStepLinearTimeTrend(592, 592, 0.01050223, 0.008251752, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialSensAn, 70, 60)
    
    N[i] <- Power
  }
  
  return(N)}


