
#Packages Used:
library(survival)
library(simsurv)
library(gsDesign)




# - SIMULATING SURVIVAL DATA - To be considered alongside Section 6.3 - 


#Function that simulates one arm of exponential survival data. There is no upper limit for the survival times when using this function, 
#as this is used to generate true event times. No censoring is applied until later steps. 

SimOneArmSurvivalData <- function(n, lambda, label){
  
  #n: number of observations to be simulated 
  #lambda: scale parameter of the exponential distribution survival times are assumed to follow
  #label: the value to assign the observations to indicate which arm the data belongs to 
  
  x = data.frame(id = 1:n)  
  
  SurvivalData <- as.data.frame(simsurv(dist = "exponential", lambdas = lambda, x=x)) #use simsurv function to simulate exponential data
  
  Label <- as.factor(rep(label,n)) #Assign label to indicate treatment arm
  
  FinalData <- cbind(Arm = Label, SurvivalData) #Combine survival object with arm label
  
  return(FinalData)
  
}



#Construct a function that takes a data set of true survival times, assigns a random recruitment time based on uniform distribution and adjusts their event times accordingly. Observations are #censored if the recruitment time results in the true event occuring outside of the trial time constraints. 

RecruitmentAdjustment <- function(SimulatedData, StartofRecruitment, EndofRecruitment, mfut){
  
  
  #Assign each patient a random (uniform distribution) recruitment time. 
  DataToBeChanged <- as.data.frame(cbind(SimulatedData, recruitmenttime = runif(nrow(SimulatedData), StartofRecruitment, EndofRecruitment)))
  
  #Use this to find the time the observation was taken in calendar time. 
  calendarTime <- DataToBeChanged$eventtime + DataToBeChanged$recruitmenttime
  
  DataToBeChanged2 <- as.data.frame(cbind(DataToBeChanged, calendarTime = calendarTime))
  
  ActualEventTime<- rep(NA, nrow(DataToBeChanged2)) #Create a vector to store actual event times after any censoring has been applied.
  
  for(i in 1:nrow(DataToBeChanged2)){
    if (DataToBeChanged2[i, ]$calendarTime >= mfut){ #If calendar time is more than the full length of the trial, event time is censored
      ActualEventTime[i] = mfut - DataToBeChanged2[i, ]$recruitmenttime #New event time will be time between recruitment and end 
      DataToBeChanged2[i, ]$status = rep(0, 1) #indicate this event time is censored
    }
    else{ 
      ActualEventTime[i] = DataToBeChanged2[i, ]$eventtime #If calendar time is not full length of trial, ActualEventTime = eventtime. 
    }
    
  }
  return(as.data.frame(cbind(DataToBeChanged2, ActualEventTime)))
}


#Construct a function to apply to a full survival data set that is otherwise finalised but no patients have dropped out 
#When there is no drop out assumed, this function does not need to be used 
#When some drop out is expected, this function alters the survival data set to reflect a certain incidence of drop out in the trial 

#Whether or not someone has dropped out is binomially distributed with a probability p of dropping out: Bin(p). 
#Use this to randomly decide whether each observation is from a patient who has dropped out, and randomly assign new event time no larger than their 
#simulated event time for patients considered to have dropped out.

DropOut <- function(DataWithoutDropOut, probdropout){ 
  
  #DataWithoutDropOut: A survival data set that does not yet take into account drop out 
  #probdropout: the expected drop out rate of the trial expressed as a decimal, for example: 0.05 for 5% drop out rate
  
  DropoutIndicator <- rbinom(n = nrow(DataWithoutDropOut), size = 1, prob = probdropout) #Create dropout indicator, binomially distributed with drop out rate p
  DataWithIndicator <- as.data.frame(cbind(DataWithoutDropOut, DropoutIndicator)) #Add the indicator variable to data set to assign observations from drop out
  
  for(i in 1:nrow(DataWithIndicator)){
    if (DataWithIndicator[i, ]$DropoutIndicator == 1) { 
      p <- runif(1, 0, 1) #Generate a random number, uniformly distributed for equal chance of any number between 0 and 1
      DataWithIndicator[i, ]$ActualEventTime <- DataWithIndicator[i, ]$ActualEventTime*p #If patients has dropped out assign random event time less than simulated time
      DataWithIndicator[i, ]$status <- rep(0, 1) #Indicate that this observation has is censored
    } else{
      DataWithIndicator[i, ]$ActualEventTime <- DataWithIndicator[i, ]$ActualEventTime #If patient hasn't dropped out event time remains the same
      DataWithIndicator[i, ]$status <- DataWithIndicator[i, ]$status #censoring status remains the same. 
    }
  }
  return(DataWithIndicator) #Returns data set with drop out taken into account. 
}

#Combine functions SimOneArmSurvivalData, RecruitmentAdjustment and DropOut into one function that simulates a finalised survival (exponential) data set

CompleteTrialSimulation <- function(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout){
  
  
  ConcurrentControlData <- SimOneArmSurvivalData(n1, lambda1, label1) #Simulate true control survival data and give control label
  
  TreatmentData <- SimOneArmSurvivalData(n2, lambda2, label2) #Simulate true treatment group survival dara and give treatment label
  
  FullTrialData <- rbind(ConcurrentControlData, TreatmentData) #Combine true control and treatment data into one data set
  
  FullTrialData <- FullTrialData[ , -2] #Remove duplicate id column
  
  FullTrialDataRecruitmentAdj <- RecruitmentAdjustment(FullTrialData, StartofRecruitment, EndofRecruitment, mfut) #Assign recruitment time and adjust event time
                                                                                                                  #censoring when appropriate
  
  FinalTrialData <- DropOut(FullTrialDataRecruitmentAdj, probdropout) #Apply drop out function to take into account a certain rate of patient drop out
  
  return(FinalTrialData) #return final survival data set
  
}



# - POWER FUNCTIONS- To be considered alongside Sections 6.2, 6.5.1, 6.5.2, and 6.6 - 

#Some notes on the following functions:   

#THE CHOICE OF LAMBDA1 AND LAMBDA2 DETERMINES WHETHER THIS FUNCTION CALCULATES POWER OR TYPE 1 ERROR RATE.
#When lambda1 = lambda2, type 1 error rate is calculated, as no true difference in survival times between treatments exists
#When lambda1 is not equal to lambda2, this function calculates power, as there is a true difference in survival times between treatment arms.

#Any functions using historical trial data should have as its input a historical data set where the TREATMENT group observations 
#have the same label as the CONTROL group observations of the current data set, given to be 1 here.
#This is what is required for the model to recognise both historical and current data as being control observations. 


#Power when using Separate Analysis

PowerSeparateAnalysis <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment) { 
  
  #n1 : number of participants in control group
  #n2 : number of participants in treatment group
  #lambda1 : exponential scale parameter for survival times in the control group
  #lambda2 : exponential scale parameter for survival times in the treatment group
  
  #mfut: Maximum follow up time = accrual period + minimum follow up 
  #label1: identifer to be given to control group
  #label2: identifier to be given to treatment group
  #probdropout: drop out rate 
  #StartofRecruitment: Usually 0 as this would be the start of the trial
  #EndofRecruitment: The length of the accrual period
  
  #Recruitment and follow up times are given in months throughout this analysis, but can be expressed as any unit of time as long as the  
  #hazard rates lambda1 and lambda2 have also been calculated using this same measurement.
  
  #M: Number of trials to be simulated
  #alpha : two sided significance level
  
  
  N <- rep(NA, M) #Create a vector in which an indicator variable of whether or not the null hypothesis was rejected for each simulated trial will be stored
  
  for (i in 1:M) {
  
    #Simulate a full current data set, representing a current survival trial
    CurrentData <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    #Use survival object to fit a cox proportional hazards model, conditional on the arm of the trial to estimate treatment effect
    ModelSeparate <- coxph(Surv(ActualEventTime, status) ~ Arm, data = CurrentData)  
    
    if (summary(ModelSeparate)$coefficients[,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if it was not
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #return proportion of trials that rejected the null hypothesis, out of those simulated.
}


#Power when using Pooled Analysis to borrow historical data

#Simulate one historical data set to use in all current data set, then use as an argument in the function.

PowerPooledAnalysis <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, HistoricalTrial, probdropout, StartofRecruitment, EndofRecruitment) { 
  
  #n1 : number of participants in control group
  #n2 : number of participants in treatment group
  #lambda1 : exponential scale parameter for survival times in the control group
  #lambda2 : exponential scale parameter for survival times in the treatment group
  #mfut: maximum follow up time = accrual period + minimum follow up
  #label1: identifier to be given to control group
  #label2: identifier to be given to treatment group
  
  #probdropout: drop out rate 
  #StartofRecruitment: Usually 0 as this would be the start of the trial
  #EndofRecruitment: The length of the accrual period
  
  #For power calculation
  #M: Number of trials to be simulated
  #alpha : two sided significance level
  
  #HistoricalTrial: the data set to be incorporated into analysis, only the observations to be used in analysis should be in this data set
  
  
  N <- rep(NA, M) #Create a vector in which an indicator variable of whether or not the null hypothesis was rejected for each simulated trial is stored
  
  for (i in 1:M) {
    
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout) #Simulate current trial 
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract only Treatment data from historical trial as:
                                                                                #treatment in historical trial = control in current trial = historical control data
    
    FinalTrialDataPooled <- rbind(CurrentTrial, HistoricalControlData) #Add historical data into data set for analysis
    
    #Use survival object to fit a cox proportional hazards model using the arm of the trial as a covariate.
    
    ModelPooled <- coxph(Surv(ActualEventTime, status) ~ Arm, data = FinalTrialDataPooled) 
    
    if (summary(ModelPooled)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if it was not
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #return proportion of trials that rejected the null hypothesis, out of those simulated.
}




#Power when modelling the time trends using a step function.

PowerStepTimeTrend <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment, 
                               HistoricalTrial) { 
  
  N <- rep(NA, M) #Create a vector in which an indicator variable of whether or not the null hypothesis was rejected for each simulated trial is stored.
  
  for (i in 1:M) {
    
    #Simulate current data set
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract the treatment data from historical trial as:
                                                                                #treatment in historical trial = control in current trial = historical control data
    
    
    Current <- rep(1 ,nrow(CurrentTrial))  #Generate vector of 1s to act as a label for current data
    Historical <- rep(0, nrow(HistoricalControlData)) #Generate a vector of 0s to act as a label for historical data
    
    CurrentTrial2 <- as.data.frame(cbind(CurrentTrial, HistOrCurrent = Current))  #Assign current data value of 1
    HistoricalControlData2 <-as.data.frame(cbind(HistoricalControlData, HistOrCurrent = Historical)) #Assign historical data value of 0
    

    CombinedFinalTrialData <- rbind(CurrentTrial2, HistoricalControlData2) #Add historical data into data set for analysis
    
    #Use survival object to fit a cox proportional hazards model, using the covariates of: the arm of the trial and whether or not the data is historical 
    #or current (which is binary and therefore is a step function with a step between 0 (historical) and 1(current))
    
    ModelStep <- coxph(Surv(ActualEventTime, status) ~ Arm + HistOrCurrent , data = CombinedFinalTrialData)  
    
    if (summary(ModelStep)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if it was not
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M)  #return proportion of trials that rejected the null hypothesis, out of those simulated.
}




#Power when modelling time trends using a linear function

PowerLinearTimeTrend <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment, HistoricalTrial,
                                 mfutHist, TimeBetweenTrials) { 
  
  N <- rep(NA, M) #Create a vector in which we can store indicator variable of whether or not the null hypothesis was rejected for each simulated trial is stored
  
  for (i in 1:M) {
    
    #Simulate current data set 
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract the treatment data from historical data set as:
                                                                                #treatment in historical trial = control in current trial = historical control data
    
    
    Current <- rep(1,nrow(CurrentTrial)) #Generate vector of 1s to act as a label for current data
    Historical <- rep(0, nrow(HistoricalControlData)) #Generate vector of 0s to act as a label for historical data
    
    CurrentTrial2 <- as.data.frame(cbind(CurrentTrial, HistOrCurrent = Current)) #Assign a value of 1 to current data
    HistoricalControlData2 <-as.data.frame(cbind(HistoricalControlData, HistOrCurrent = Historical)) #Assign a value of 0 to historical data
    
    CombinedFinalTrialData <- rbind(CurrentTrial2, HistoricalControlData2) #Add historical data into data set for analysis
    
    #Make variable to used as a linear time trend
    
    RecruitmentFromBaseline <- rep(NA, n1+n2+nrow(HistoricalControlData2)) #Create vector to store recruitmentfrombaseline times
    
    for(j in 1:nrow(CombinedFinalTrialData)){
      
      if(CombinedFinalTrialData[j, ]$HistOrCurrent == 0){ #If the observation is historical 
        
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime  #Baseline recruitment is the recruitment time. 
      }
      else{ 
      #When the observation is current, baseline recruitment is recruitment time of current trial, in addition to the time since the start of the historical trial.
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime + mfutHist + TimeBetweenTrials  
      }
    }
    
    CombinedFinalTrialData2 <- cbind(CombinedFinalTrialData, RecruitmentFromBaseline)  #Add baseline recruitment time to the data set
    
    #Use survival object to fit a cox proportional hazards model, using the covariates of: the arm of the trial and the baseline recruitment time (which acts as 
    linear time between the start of the historical trial and the end of current trial recruitment) 
    
    ModelLinear <- coxph(Surv(ActualEventTime, status) ~ Arm + RecruitmentFromBaseline, data = CombinedFinalTrialData2) 
    
    if (summary(ModelLinear)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if it was not
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #return proportion of trials that rejected the null hypothesis, out of those simulated. 
}



#Power when modelling the time trend using both a linear and step function


PowerStepLinearTimeTrend <- function(n1, n2, lambda1, lambda2, M, mfut, alpha, label1, label2, probdropout, StartofRecruitment, EndofRecruitment, HistoricalTrial,
                                     mfutHist, TimeBetweenTrials) { 
  
  N <- rep(NA, M) #Create a vector in which an indicator variable of whether or not the null hypothesis was rejected for each simulated trial is stored.
  
  for (i in 1:M) {
    
    #Simulate current data set 
    
    CurrentTrial <- CompleteTrialSimulation(n1, lambda1, mfut, label1, n2, lambda2, label2, StartofRecruitment, EndofRecruitment, probdropout)
    
    HistoricalControlData <- HistoricalTrial[which(HistoricalTrial$Arm == 1), ] #Extract the treatment data from historical trial as:
                                                                                ##treatment in historical trial = control in current trial = historical control data
    
    Current <- rep(1,nrow(CurrentTrial)) #Generate vector of 1s to act as a label for current data
    Historical <- rep(0, nrow(HistoricalControlData)) #Generate vector of 0s to act as a label for historical data
    
    CurrentTrial2 <- as.data.frame(cbind(CurrentTrial, HistOrCurrent = Current)) #Assign a value of 1 to all current observations
    HistoricalControlData2 <-as.data.frame(cbind(HistoricalControlData, HistOrCurrent = Historical)) #Assign a value of 0 to historical observations
    
    CombinedFinalTrialData <- rbind(CurrentTrial2, HistoricalControlData2) #Add historical data into data set for analysis
    
    #Make variable to used as a linear time trend
    
    RecruitmentFromBaseline <- rep(NA, n1+n2+nrow(HistoricalControlData2)) #Create vector to store recruitmentfrombaseline times
    
    for(j in 1:nrow(CombinedFinalTrialData)){
      
      if(CombinedFinalTrialData[j, ]$HistOrCurrent == 0){ #If the observation is historical 
        
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime #Baseline recruitment is the recruitment time.
      }
      else{ 
      #When the observation is current, baseline recruitment is recruitment time of current trial, in addition to the time since the start of the historical trial.
        RecruitmentFromBaseline[j] = CombinedFinalTrialData[j, ]$recruitmenttime + mfutHist + TimeBetweenTrials 
      }
    }
    
    CombinedFinalTrialData2 <- cbind(CombinedFinalTrialData, RecruitmentFromBaseline) #Add baseline recruitment time to the data set
    
    #Use survival object to fit a cox proportional hazards model, using the covariates of: the arm of the trial, the baseline recruitment time (which acts as 
    linear time between the start of the historical trial and the end of current trial recruitment) and whether or not the data is historical or current (which 
    #is binary and therefore is a step function with a step between 0 (historical) and 1(current))
    
    ModelStepLinear <- coxph(Surv(ActualEventTime, status) ~ Arm + RecruitmentFromBaseline + HistOrCurrent , data = CombinedFinalTrialData2) 
    
    
    if (summary(ModelStepLinear)$coefficients[1,5] < alpha) { #Assign a value of 1 to the vector N if the null hypothesis was rejected, and zero if it was not
      N[i] <- 1
    } else {
      N[i] <- 0 
    }
  }
  return(sum(N)/M) #return proportion of trials that rejected the null hypothesis, out of those simulated. 
}








# - MAIN DATA ANALYSIS - To be considered alongside Section 6.6.1

#Simulate the Historical Trial  

set.seed(1234)
HistoricalTrialMyeloma <- CompleteTrialSimulation(100, 0.03500743, 72, 0, 100, 0.02100446, 1, 0, 36, 0.05)

#Check that the historical data trial is powered to at least 80%
set.seed(1234)
PowerHistTrialMyeloma <- PowerSeparateAnalysis(100, 100, 0.03500743, 0.02100446, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.873

#Calculate the Monte Carlo Error estimate of this power
MCVarianceSeparateHistMyeloma <-(PowerHistTrialMyeloma*(1 - PowerHistTrialMyeloma))/1000
MCSDSeparateHistMyeloma <- sqrt(MCVarianceSeparateHistMyeloma) #0.01052953


#Separate Analysis of current data set 
set.seed(1234)
SeparatePowerHistMyeloma <- PowerSeparateAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48) #0.8 

#Calculate the Monte Carlo Error estimate of this power
MCVarianceSeparateMyeloma <-(SeparatePowerHistMyeloma*(1 - SeparatePowerHistMyeloma))/1000
MCSDSeparateMyeloma <- sqrt(MCVarianceSeparateMyeloma) #0.01264911



#Use methods of historical borrowing incorporate HistoricalTrialMyeloma into the analysis of the current trial

#Pooled Analysis
set.seed(1234)
PooledPowerHistMyeloma <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrialMyeloma, 0.05, 0, 48) #0.84

#Calculate the estimated Monte Carlo Error of this power
MCVariancePooledMyeloma <-(PooledPowerHistMyeloma*(1 - PooledPowerHistMyeloma))/1000
MCSDPooledMyeloma <- sqrt(MCVariancePooledMyeloma) #0.0115931


#Modelling time trend using a step function
set.seed(1234)
PowerStepHistMyeloma <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma) #0.801

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepMyeloma <-(PowerStepHistMyeloma*(1 - PowerStepHistMyeloma))/1000
MCSDStepMyeloma <- sqrt(MCVarianceStepMyeloma) #0.01262533



#Modelling time trend using a linear function
set.seed(1234)
PowerLinearHistMyeloma <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 70, 60) #0.821

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinearMyeloma <-(PowerLinearHistMyeloma*(1 - PowerLinearHistMyeloma))/1000
MCSDLinearMyeloma <- sqrt(MCVarianceLinearMyeloma) #0.01212266



#Modelling time trend using a linear and step function
set.seed(1234)
PowerStepLinearHistMyeloma <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 70, 60) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinearMyeloma <-(PowerStepLinearHistMyeloma*(1 - PowerStepLinearHistMyeloma))/1000
MCSDStepLinearMyeloma <- sqrt(MCVarianceStepLinearMyeloma) #0.01255325


#Type 1 Errors

#Seperate Analysis of current trial
set.seed(1234)
Type1ErrorSeparateMyeloma <- PowerSeparateAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1SeparateMyeloma <-(Type1ErrorSeparateMyeloma*(1 - Type1ErrorSeparateMyeloma))/1000
MCSDType1SeparateMyeloma <- sqrt(MCVarianceType1SeparateMyeloma) #0.007331507

#Pooled Analysis 
set.seed(1234)
Type1ErrorPooledMyeloma <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrialMyeloma, 0.05, 0, 48) #0.05

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1PooledMyeloma <-(Type1ErrorPooledMyeloma*(1 - Type1ErrorPooledMyeloma))/1000
MCSDType1PooledMyeloma <- sqrt(MCVarianceType1PooledMyeloma) #0.006892024

#Modelling time trends using a step function
set.seed(1234)
Type1ErrorTimeTrendMyeloma <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepMyeloma <-(Type1ErrorTimeTrendMyeloma*(1 - Type1ErrorTimeTrendMyeloma))/1000
MCSDType1StepMyeloma <- sqrt(MCVarianceType1StepMyeloma) #0.007331507

#Modelling time trends using a linear function
set.seed(1234)
Type1ErrorLinearTrendMyeloma <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 60) #0.056

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1LinearMyeloma <-(Type1ErrorLinearTrendMyeloma*(1 - Type1ErrorLinearTrendMyeloma))/1000
MCSDType1LinearMyeloma <- sqrt(MCVarianceType1LinearMyeloma) #0.007270763

#Modelling time trends using a step and linear function
set.seed(1234)
Type1ErrorStepLinearTrendMyeloma <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 60) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinearMyeloma <-(Type1ErrorStepLinearTrendMyeloma*(1 - Type1ErrorStepLinearTrendMyeloma))/1000
MCSDType1StepLinearMyeloma <- sqrt(MCVarianceType1StepLinearMyeloma) #0.007147307



#Summary

Methods <- c("Separate", "Pooling", "Step Time", "Linear Time", "Step Linear Time")

MainAnalysisPowers <- c(SeparatePowerHistMyeloma, PooledPowerHistMyeloma, PowerStepHistMyeloma, PowerLinearHistMyeloma, PowerStepLinearHistMyeloma)

MainAnalysisType1Errors <- c(Type1ErrorSeparateMyeloma, Type1ErrorPooledMyeloma,Type1ErrorTimeTrendMyeloma, Type1ErrorLinearTrendMyeloma, Type1ErrorStepLinearTrendMyeloma)

MainAnalysisPowersMCE <- c(MCSDSeparateMyeloma, MCSDPooledMyeloma, MCSDStepMyeloma,MCSDLinearMyeloma, MCSDStepLinearMyeloma)

MainAnalysisType1ErrorsMCE <- c(MCSDType1SeparateMyeloma,MCSDType1PooledMyeloma, MCSDType1StepMyeloma, MCSDType1LinearMyeloma, MCSDType1StepLinearMyeloma)

MainAnalysisData <- as.data.frame(cbind(Methods, MainAnalysisPowers, MainAnalysisType1Errors, MainAnalysisPowersMCE, MainAnalysisType1ErrorsMCE))

#Plotting the power and type 1 Error

#Powers
par(mfrow = c(1,2))

plot(MainAnalysisPowers, data = MainAnalysisData, xlab = "Method of Historical Borrowing", ylab = "Power", Main = "Graph Showing Change in Power for Methods of Borrowing the HistoricalTrialMyeloma Data Set", 
     type = "b", xaxt='n', main = "Change in Power when Borrowing A Historical 
     Data Set with The Same Underlying Parameter 
     Using Different Methods")
axis(1, at = 1:nrow(MainAnalysisData), labels = Methods)


#Type 1 Errors
plot(MainAnalysisType1Errors, data = MainAnalysisData, xlab = "Method of Historical Borrowing", ylab = "Type 1 Error", 
     type = "b", xaxt='n', main = "Change in Type 1 Error when Borrowing 
     A Historical Data Set with The Same 
     Underlying Parameter Using Different Methods")
axis(1, at = 1:nrow(MainAnalysisData), labels = Methods)

















# - MAIN DATA ANALYSIS - To be considered alongside Section 6.6.2



#Now, since all of those powers calculated are conditional on the historical data set, we should make some changes to the historical data set 
#And see how the results are affected

#First, lets make changes to the parameters, more specifically, how close lambda2 in the historical data set is to lambda1 in the current data set (For same treatment)
#The lambda in the current control is 0.02100446, so lets call this 0.021 correct to 3 dp. 
#We can test out what happens to the power of the method if we make the historical lambda 0.18, 0.19, 0.2, 0.22, These values have been chosen to see what a change in 0.001 and 0.002 does
#and what it does when its higher or lower
#More values can be tested but in this case it was too computationally intensive to do more than this. 


#So when the parameter is 0.19 or 0.2, the risk of treatment A was less than it is currently assumed. 
#We expect this to skew the parameter estimation of the control group to one closer to the treatment, and we may be less likely to conclude change in treatment. 
#This may result in a higher false positive rate, or alternatively a lower power. 

#When the parameter is 0.22 or 0.23, the risk of treatment A was more historically than is currently assumed. 
#We expect this to skew the parameter estimation of the control to one further away than the treatment, and may be more likely to conclude change in treatment 
#This should result in a higher power, because more trials are rejecting the null and we know that there is indeed a treatment difference 
#But, the rate of false positives, the type I error rate, may increase. Alpha should be bigger. 



#Make 4 more data sets, each with a different value for lambda2, and lambda1 changes as well in order to keep the 0.6 hazard ratio 

set.seed(1234)
HistoricalTrial0.018 <- CompleteTrialSimulation(100, 0.03, 72, 0, 100, 0.018, 1, 0, 36, 0.05)

set.seed(1234)
HistoricalTrial0.019 <- CompleteTrialSimulation(100, 0.03166667, 72, 0, 100, 0.019, 1, 0, 36, 0.05)

set.seed(1234)
HistoricalTrial0.02 <- CompleteTrialSimulation(100, 0.03333333, 72, 0, 100, 0.02, 1, 0, 36, 0.05)

set.seed(1234)
HistoricalTrial0.022 <- CompleteTrialSimulation(100, 0.03666667, 72, 0, 100, 0.022, 1, 0, 36, 0.05)

set.seed(1234)
HistoricalTrial0.023 <- CompleteTrialSimulation(100, 0.03833333, 72, 0, 100, 0.023, 1, 0, 36, 0.05)

#We should check that they appear to have enough power in them as individual trials:

#Separate Analysis

set.seed(1234)
SeparatePowerHist0.018 <- PowerSeparateAnalysis(100, 100, 0.03, 0.018, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.851

set.seed(1234)
SeparatePowerHist0.019 <- PowerSeparateAnalysis(100, 100, 0.03166667, 0.019, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.856

set.seed(1234)
SeparatePowerHist0.02 <- PowerSeparateAnalysis(100, 100, 0.03333333, 0.02, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.858

set.seed(1234)
PowerHistTrialMyeloma <- PowerSeparateAnalysis(100, 100, 0.03500743, 0.02100446, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.873

set.seed(1234)
SeparatePowerHist0.022 <- PowerSeparateAnalysis(100, 100, 0.03666667, 0.022, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.877

set.seed(1234)
SeparatePowerHist0.023 <- PowerSeparateAnalysis(100, 100, 0.03833333, 0.023, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.888


#So when we incorporate each of the historical data sets, we know we are incorporating data sets that had adequate power themselves. 

#Summary 

Histlambda0.018 <- 0.018
Histlambda0.019 <- 0.019
Histlambda0.02 <- 0.02
HistlambdaMyeloma <- 0.02100446
Histlambda0.022 <- 0.022
Histlambda0.023 <- 0.023

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)

SeparatePowersHist <- c(SeparatePowerHist0.018,SeparatePowerHist0.019, SeparatePowerHist0.02, PowerHistTrialMyeloma, SeparatePowerHist0.022, SeparatePowerHist0.023)

PowersHistorialData <- as.data.frame(cbind(Lambda2 = lambdas, "Power of Historical Trial" = SeparatePowersHist))

PowersHistorialData


#Applying the methods to the 4 new data sets with a different value for lambda

#POWER

#Let us start with the POOLED ANALYSIS. How does the power change when we include data from a historical data set that is not the same, but we treat it as the same. 

#Pooled Analysis when the data set is exactly the same gives us a power of 0.84. 

PooledPowerHistMyeoma #0.84

#We would expect to get lower powers when the underlying parameters are not the same. 

#Pooled Analysis
set.seed(1234)
PooledPowerHist0.018 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.018, 0.05, 0, 48) #0.763
set.seed(1234)
PooledPowerHist0.019 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.019, 0.05, 0, 48) #0.794
set.seed(1234)
PooledPowerHist0.02 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.02, 0.05, 0, 48) #0.827
set.seed(1234)
PooledPowerHist0.022 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.022, 0.05, 0, 48) #0.848
set.seed(1234)
PooledPowerHist0.023 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.023, 0.05, 0, 48) #0.881


par(mfrow = c(2,2))

#GRAPH OF THIS

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
PooledPowers <- c(PooledPowerHist0.018, PooledPowerHist0.019, PooledPowerHist0.02, PooledPowerHistMyeloma, PooledPowerHist0.022, PooledPowerHist0.023)

plot(lambdas, PooledPowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of the Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Using Pooled Analysis")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PooledPowerHistMyeloma, lty = "dotted", col = "purple")



#Now what happens when we model the change in time as a STEP FUNCTION. When the parameter was exactly the same, using this method kept the power roughly the same: 0.801 

PowerStepHistMyeloma  #0.801

#So what about when it is not exactly the same. 

set.seed(1234)
PowerStepHist0.018 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018) #0.803

set.seed(1234)
PowerStepHist0.019 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019) #0.802

set.seed(1234)
PowerStepHist0.02 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02) #0.802

#For these three we actually see an increase in power. This could be because the model is detecting a change where there actually is one and is modelling it

set.seed(1234)
PowerStepHist0.022 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022) #0.798

set.seed(1234)
PowerStepHist0.023 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023) #0.801

#This last one gets a slightly lower power. Not quite sure why. 

#Graph of this 

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
StepTimePowers <- c(PowerStepHist0.018, PowerStepHist0.019, PowerStepHist0.02, PowerStepHistMyeloma, PowerStepHist0.022, PowerStepHist0.023)

plot(lambdas, StepTimePowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of the Current Trial when Borrowing
     Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PowerStepHistMyeloma, lty = "dotted", col = "purple")


#Now lets try the LINEAR function: 

#When lambda2 is exactly the same, this gives a power of 0.821, 
PowerLinearHistMyeloma #0.821

#Let us see what happens when we change the parameters. 

set.seed(1234)
PowerLinearHist0.018 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 70, 60) #0.815

set.seed(1234)
PowerLinearHist0.019 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 70, 60) #0.816

set.seed(1234)
PowerLinearHist0.02 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 70, 60) #0.821

set.seed(1234)
PowerLinearHist0.022 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 70, 60) #0.821

set.seed(123)
PowerLinearHist0.023 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 70, 60) #0.83


#Do we get an increase or decrease from 0.821. Since we're incorporating data that has a step change and we're modelling it only in a linear way, might be a misspecified model. 
#And we might expect to get a smaller power. 


#Graph of this 

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
LinearTimePowers <- c(PowerLinearHist0.018, PowerLinearHist0.019, PowerLinearHist0.02, PowerLinearHistMyeloma, PowerLinearHist0.022, PowerLinearHist0.023)

plot(lambdas, LinearTimePowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of the Current Trial when Borrowing
     Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PowerLinearHistMyeloma, lty = "dotted", col = "purple")


#And finally, we can see what happens when we use both the step and the linear model. 

#When using the main trial we got a power of 

PowerStepLinearHistMyeloma #0.804

#What happens when we change lambda2

set.seed(1234)
PowerStepLinearHist0.018 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 70, 60) #0.804

set.seed(1234)
PowerStepLinearHist0.019 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 70, 60) #0.803

set.seed(1234)
PowerStepLinearHist0.02 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 70, 60) #0.805

set.seed(1234)
PowerStepLinearHist0.022 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 70, 60) #0.802

set.seed(1234)
PowerStepLinearHist0.023 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 70, 60) #0.802

#What happens to the power? 


#Graph of This
lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
StepLinearTimePowers <- c(PowerStepLinearHist0.018, PowerStepLinearHist0.019, PowerStepLinearHist0.02, PowerStepLinearHistMyeloma, PowerStepLinearHist0.022, PowerStepLinearHist0.023)

plot(lambdas, StepLinearTimePowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of Current Trial when Borrowing Historical 
Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step and Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PowerStepLinearHistMyeloma, lty = "dotted", col = "purple")

#Table of powers when using methods for varying lambdas 


lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)

PooledPowers <- c(PooledPowerHist0.018, PooledPowerHist0.019, PooledPowerHist0.02, PooledPowerHistMyeloma, PooledPowerHist0.022, PooledPowerHist0.023)
StepTimePowers <- c(PowerStepHist0.018, PowerStepHist0.019, PowerStepHist0.02, PowerStepHistMyeloma, PowerStepHist0.022, PowerStepHist0.023)
LinearTimePowers <- c(PowerLinearHist0.018, PowerLinearHist0.019, PowerLinearHist0.02, PowerLinearHistMyeloma, PowerLinearHist0.022, PowerLinearHist0.023)
StepLinearTimePowers <- c(PowerStepLinearHist0.018, PowerStepLinearHist0.019, PowerStepLinearHist0.02, PowerStepLinearHistMyeloma, PowerStepLinearHist0.022, PowerStepLinearHist0.023)

AllPowersLambda <- as.data.frame(cbind(lambdas, PooledPowers, StepTimePowers, LinearTimePowers, StepLinearTimePowers)) 
AllPowersLambda


#TYPE 1 ERROR RATE



#We also need to take into account the Type I error rates 

#Pooled 

#The usual error is 
Type1ErrorPooledMyeloma #0.05

#Now we find it for the other values of lambda
set.seed(1234)
Type1ErrorPooled0.018 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.018, 0.05, 0, 48) #0.057

set.seed(1234)
Type1ErrorPooled0.019 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.019, 0.05, 0, 48) #0.053

set.seed(1234)
Type1ErrorPooled0.02 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.02, 0.05, 0, 48) #0.049

set.seed(1234)
Type1ErrorPooled0.022 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.022, 0.05, 0, 48) #0.048

set.seed(1234)
Type1ErrorPooled0.023 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.023, 0.05, 0, 48) #0.057

par(mfrow = c(2,2))

#Graph of this

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
PooledType1Errors <- c(Type1ErrorPooled0.018, Type1ErrorPooled0.019, Type1ErrorPooled0.02, Type1ErrorPooledMyeloma,Type1ErrorPooled0.022,Type1ErrorPooled0.023)

plot(lambdas, PooledType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Using Pooled Analysis")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorPooledMyeloma, lty = "dotted", col = "purple")


#Step

#The usual one is 
Type1ErrorTimeTrendMyeloma #0.062

#Now we find it for the other values of lambda
set.seed(1234)
Type1ErrorTimeTrend0.018 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018) #0.057

set.seed(1234)
Type1ErrorTimeTrend0.019 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019) #0.057

set.seed(1234)
Type1ErrorTimeTrend0.02 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02) #0.057

set.seed(1234)
Type1ErrorTimeTrend0.022 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022) #0.057

set.seed(1234)
Type1ErrorTimeTrend0.023 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023) #0.057

#Graph of this:

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
StepTimeType1Errors <- c(Type1ErrorTimeTrend0.018, Type1ErrorTimeTrend0.019 , Type1ErrorTimeTrend0.02, Type1ErrorTimeTrendMyeloma, Type1ErrorTimeTrend0.022, Type1ErrorTimeTrend0.023)

plot(lambdas, StepTimeType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorTimeTrendMyeloma, lty = "dotted", col = "purple")


#Linear 

#The usual one is 
Type1ErrorLinearTrendMyeloma #0.056

#Now we find it for the other values of lambda
set.seed(1234)
Type1ErrorLinearTrend0.018 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 72, 60) #0.056

set.seed(1234)
Type1ErrorLinearTrend0.019 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 72, 60) #0.056

set.seed(1234)
Type1ErrorLinearTrend0.02 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 72, 60) #0.057

set.seed(1234)
Type1ErrorLinearTrend0.022 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 72, 60) #0.057

set.seed(1234)
Type1ErrorLinearTrend0.023 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 72, 60) #0.059

#Graph of this:

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
LinearTimeType1Errors <- c(Type1ErrorLinearTrend0.018, Type1ErrorLinearTrend0.019, Type1ErrorLinearTrend0.02, Type1ErrorLinearTrendMyeloma, Type1ErrorLinearTrend0.022, Type1ErrorLinearTrend0.023)

plot(lambdas, LinearTimeType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorLinearTrendMyeloma, lty = "dotted", col = "purple")


#Step and Linear

#The usual one is 
Type1ErrorStepLinearTrendMyeloma  #0.054

#Now we find it for the other values of lambda
set.seed(1234)
Type1ErrorStepLinearTrend0.018 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 72, 60) #0.055

set.seed(1234)
Type1ErrorStepLinearTrend0.019 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 72, 60) #0.056

set.seed(1234)
Type1ErrorStepLinearTrend0.02 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 72, 60) #0.054

set.seed(1234)
Type1ErrorStepLinearTrend0.022 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 72, 60) #0.054

set.seed(1234)
Type1ErrorStepLinearTrend0.023 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 72, 60) #0.54

#Graph of This

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
StepLinearTimeType1Errors <- c(Type1ErrorStepLinearTrend0.018, Type1ErrorStepLinearTrend0.019, Type1ErrorStepLinearTrend0.02, Type1ErrorStepLinearTrendMyeloma,Type1ErrorStepLinearTrend0.022, Type1ErrorStepLinearTrend0.023)

plot(lambdas, StepLinearTimeType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step and Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorStepLinearTrendMyeloma, lty = "dotted", col = "purple")






#GRAPHS FOR CHANGE IN HISTORICAL LAMBDA

#Summary

#Historical Lambdas
Histlambda0.018 <- 0.018
Histlambda0.019 <- 0.019
Histlambda0.02 <- 0.02
HistlambdaMyeloma <- 0.02100446
Histlambda0.022 <- 0.022
Histlambda0.023 <- 0.023

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)

#Power
PooledPowers <- c(PooledPowerHist0.018, PooledPowerHist0.019, PooledPowerHist0.02, PooledPowerHistMyeloma, PooledPowerHist0.022, PooledPowerHist0.023)

StepTimePowers <- c(PowerStepHist0.018, PowerStepHist0.019, PowerStepHist0.02, PowerStepHistMyeloma, PowerStepHist0.022, PowerStepHist0.023)

LinearTimePowers <- c(PowerLinearHist0.018, PowerLinearHist0.019, PowerLinearHist0.02, PowerLinearHistMyeloma, PowerLinearHist0.022, PowerLinearHist0.023)

StepLinearTimePowers <- c(PowerStepLinearHist0.018, PowerStepLinearHist0.019, PowerStepLinearHist0.02, PowerStepLinearHistMyeloma, PowerStepLinearHist0.022, PowerStepLinearHist0.023)

#Type 1 Error
PooledType1Errors <- c(Type1ErrorPooled0.018, Type1ErrorPooled0.019, Type1ErrorPooled0.02, Type1ErrorPooledMyeloma,Type1ErrorPooled0.022, Type1ErrorPooled0.023)

StepTimeType1Errors <- c(Type1ErrorTimeTrend0.018, Type1ErrorTimeTrend0.019 , Type1ErrorTimeTrend0.02, Type1ErrorTimeTrendMyeloma, Type1ErrorTimeTrend0.022, Type1ErrorTimeTrend0.023)

LinearTimeType1Errors <- c(Type1ErrorLinearTrend0.018, Type1ErrorLinearTrend0.019, Type1ErrorLinearTrend0.02, Type1ErrorLinearTrendMyeloma, Type1ErrorLinearTrend0.022, Type1ErrorLinearTrend0.023)

StepLinearTimeType1Errors <- c(Type1ErrorStepLinearTrend0.018, Type1ErrorStepLinearTrend0.019, Type1ErrorStepLinearTrend0.02, Type1ErrorStepLinearTrendMyeloma,Type1ErrorStepLinearTrend0.022, Type1ErrorStepLinearTrend0.023)


DataHistLambda <- as.data.frame(cbind(Lambda = lambdas, PoolPower = PooledPowers, PoolError =PooledType1Errors, StepTimePower = StepTimePowers, 
                        StepTimeError = StepTimeType1Errors, LinearTimePower = LinearTimePowers, LinearTimeError = LinearTimeType1Errors, 
                        StepLinearTimePower = StepLinearTimePowers, StepLinearTimeError = StepLinearTimeType1Errors))

DataHistLambda

par(mfrow = c(1,2))

#Plotting Power
Historicallambdas <- c(rep(lambdas, 4)) #x values of the points to plot

PowersLambdas <- c(PooledPowers, StepTimePowers, LinearTimePowers, StepLinearTimePowers) #y values of the points to plot

MethodsLambdas <- c(rep("Pooling", 5), rep("Step", 5),rep("Linear", 5),rep("StepLinear", 5)) #Methods to colour points by

PowerDataLambdas <- as.data.frame(cbind(HistoricalLambda = Historicallambdas, Power = PowersLambdas, Method = as.factor(MethodsLambdas)))


plot(DataLambdas$HistoricalLambda, DataLambdas$Power, col = DataLambdas$Method, pch = 4, xlab = "Historical Lambda2", ylab = "Power", 
     main = "Plot of Power of the Current Trial when Borrowing Historical 
     Data using Various Methods, from Historical 
     Data with Varying Values of Lambda2") #Plot all powers of all methods

#Put lines through points to show trend for each method
lines(lambdas, PooledPowers, col = "red")
lines(lambdas,StepTimePowers, col = "green")
lines(lambdas,LinearTimePowers, col = "black")
lines(lambdas,StepLinearTimePowers, col = "blue")

#Add lines indicating lambda used in main analysis
abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")

#Add lines indicating power when using lambda in main analysis
abline(h = PooledPowerHistMyeloma, lty = "dotted", col = "red")
abline(h = PowerStepHistMyeloma, lty = "dotted", col = "green")
abline(h = PowerLinearHistMyeloma, lty = "dotted", col = "black")
abline(h = PowerStepLinearHistMyeloma, lty = "dotted", col = "blue")

#Add legend

legend("bottomright", legend=c("Pooled", "Step Time", "Linear Time", "Step and Linear Time"),col= c("red", "green", "black", "blue"), lty=1, cex=0.8)





#Plotting Type 1 Errors

Historicallambdas <- c(rep(lambdas, 4)) #x values of points to plot

ErrorsLambdas <- c(PooledType1Errors, StepTimeType1Errors, LinearTimeType1Errors, StepLinearTimeType1Errors) 

MethodsLambdas <- c(rep("Pooling", 5), rep("Step", 5),rep("Linear", 5),rep("StepLinear", 5)) #Methods to colour points by

Data2Lambdas <- as.data.frame(cbind(HistoricalLambda = Historicallambdas, Type1Error = ErrorsLambdas, Method = as.factor(MethodsLambdas)))

plot(Data2Lambdas$HistoricalLambda, Data2Lambdas$Type1Error, col = Data2Lambdas$Method, pch = 4, 
     xlab = "Historical Lambda2", ylab = "Type 1 Error", main = "Plot of the Type 1 Error of the Current Trial when Borrowing 
     Historical Data UsingVarious Methods, from Historical 
     Data with Varying Values of Lambda2") #Plot all errors from all methods

#Put lines through points to show trend for each method
lines(lambdas, PooledType1Errors, col = "red")
lines(lambdas, StepTimeType1Errors, col = "green")
lines(lambdas, LinearTimeType1Errors, col = "black")
lines(lambdas, StepLinearTimeType1Errors, col = "blue")

#Add lines indicating lambda used in main analysis
abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")

#Add lines indicating power when using lambda in main analysis
abline(h = Type1ErrorPooledMyeloma, lty = "dotted", col = "red")
abline(h = Type1ErrorTimeTrendMyeloma, lty = "dotted", col = "green")
abline(h = TypeIErrorLinearTrendMyeloma, lty = "dotted", col = "black")
abline(h = TypeIErrorStepLinearTrendMyeloma, lty = "dotted", col = "blue")

#Add legend

legend("topleft", legend=c("Pooled", "Step Time", "Linear Time", "Step and Linear Time"),col= c("red", "green", "black", "blue"), lty=1, cex=0.8)


#Changing the time between historical and current trials 

#Now that we have seen how all of the methods work when there is a change in underlying parameter, let us look at how the power changes between more recent and older trials. 

#We have been assuming that the trial was done 5 years ago, and this gives us a power of 0.821

#Let us see what happens if it was done 1-4 years ago. #This time we don't need to change the data set
#We keep the data set the same, but just change the time between trials in both of the functions that allow us to vary this: 

#First, if we model the time trend as purely linear: 

set.seed(1234)
PowerLinear1Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.824

set.seed(1234)
PowerLinear2Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.821

set.seed(1234)
PowerLinear3Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.82

set.seed(1234)
PowerLinear4Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.82

set.seed(1234)
PowerLinear6Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.82

#We don't really see much change in the power, but we are seeing a very slight increase in power when they're more recent. 

par(mfrow = c(1,2))

#Put this on a Graph itself 

Years <- c(1, 2, 3, 4, 5, 6)
LinearPowersYears <- c(PowerLinear1Year, PowerLinear2Year, PowerLinear3Year, PowerLinear4Year, PowerLinearHistMyeloma, PowerLinear6Year)

plot(Years, LinearPowersYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Power", main = "Plot of Power of the Current Trial When Borrowing
   Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = PowerLinearHistMyeloma, lty = "dotted", col = "purple")

#Now let us see what changes when we do both a step and linear trend. Assuming that the trial was done 5 years ago, this gives us a power of 0.804. Is there any increase or decrease in that?

set.seed(1234)
PowerStepLinear1Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.804

set.seed(1234)
PowerStepLinear2Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.804

set.seed(1234)
PowerStepLinear3Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.804

set.seed(1234)
PowerStepLinear4Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.804

set.seed(1234)
PowerStepLinear6Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.804

#Do a graph of this

Years <- c(1, 2, 3, 4, 5, 6)
StepLinearPowersYears <- c(PowerStepLinear1Year, PowerStepLinear2Year, PowerStepLinear3Year, PowerStepLinear4Year, PowerStepLinearHistMyeloma, PowerStepLinear6Year)

plot(Years, StepLinearPowersYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Power", main = "Plot of Power of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear and Step Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = PowerStepLinearHistMyeloma, lty = "dotted", col = "purple")



#And we also need to take into account the Type 1 Error

#Linear 

set.seed(1234)
Type1ErrorLinearTimeTrend1Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.059

set.seed(1234)
Type1ErrorLinearTimeTrend2Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.058

set.seed(1234)
Type1ErrorLinearTimeTrend3Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.059

set.seed(1234)
Type1ErrorLinearTimeTrend4Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.059

set.seed(1234)
Type1ErrorLinearTimeTrend6Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.057

#Graph the Type 1 Errors

par(mfrow = c(1,2))

Years <- c(1, 2, 3, 4, 5, 6)
LinearType1ErrorsYears <- c(Type1ErrorLinearTimeTrend1Year, Type1ErrorLinearTimeTrend2Year, Type1ErrorLinearTimeTrend3Year, Type1ErrorLinearTimeTrend4Year, TypeIErrorLinearTrendMyeloma, Type1ErrorLinearTimeTrend6Year)

plot(Years, LinearType1ErrorsYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Type 1 Error", main = "Plot of Type 1 Error of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = TypeIErrorLinearTrendMyeloma, lty = "dotted", col = "purple")

#Step and Linear 


set.seed(1234)
Type1ErrorStepLinearTimeTrend1Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.054

set.seed(1234)
Type1ErrorStepLinearTimeTrend2Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.054

set.seed(1234)
Type1ErrorStepLinearTimeTrend3Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.054

set.seed(1234)
Type1ErrorStepLinearTimeTrend4Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.054

set.seed(1234)
Type1ErrorStepLinearTimeTrend6Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.054

#Graph the Type 1 Errors

Years <- c(1, 2, 3, 4, 5, 6)
StepLinearType1ErrorsYears <- c(Type1ErrorStepLinearTimeTrend1Year, Type1ErrorStepLinearTimeTrend2Year, Type1ErrorStepLinearTimeTrend3Year, Type1ErrorStepLinearTimeTrend4Year, TypeIErrorStepLinearTrendMyeloma,
                                Type1ErrorStepLinearTimeTrend6Year)

plot(Years, StepLinearType1ErrorsYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Type 1 Error", main = "Plot of Type 1 Error of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear and Step Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = TypeIErrorStepLinearTrendMyeloma, lty = "dotted", col = "purple")



#GRAPH FOR CHANGING TIME BETWEEN TRIALS 

#Summary

#Years
Years <- c(1, 2, 3, 4, 5, 6)

#Powers
LinearPowersYears <- c(PowerLinear1Year, PowerLinear2Year, PowerLinear3Year, PowerLinear4Year, PowerLinearHistMyeloma, PowerLinear6Year)
StepLinearPowersYears <- c(PowerStepLinear1Year, PowerStepLinear2Year, PowerStepLinear3Year, PowerStepLinear4Year, PowerStepLinearHistMyeloma, PowerStepLinear6Year)

#Type 1 Errors
LinearType1ErrorsYears <- c(Type1ErrorLinearTimeTrend1Year, Type1ErrorLinearTimeTrend2Year, Type1ErrorLinearTimeTrend3Year, Type1ErrorLinearTimeTrend4Year, TypeIErrorLinearTrendMyeloma, Type1ErrorLinearTimeTrend6Year)
StepLinearType1ErrorsYears <- c(Type1ErrorStepLinearTimeTrend1Year, Type1ErrorStepLinearTimeTrend2Year, Type1ErrorStepLinearTimeTrend3Year, Type1ErrorStepLinearTimeTrend4Year, TypeIErrorStepLinearTrendMyeloma, Type1ErrorStepLinearTimeTrend6Year)

DataHistTime <- as.data.frame(cbind(YearsBetweenTrials = Years, LinearPower = LinearPowersYears, LinearErrors = LinearType1ErrorsYears, 
                                    StepLinearPower = StepLinearPowersYears , StepLinearError = StepLinearType1ErrorsYears))

DataHistTime

par(mfrow = c(1,2))

#Graph for Powers

YearsxAxis <- c(rep(Years, 2))
PowersYears <- c(LinearPowersYears, StepLinearPowersYears)
MethodsYears <- c(rep("Linear", 5),rep("StepLinear", 5))
DataYears <- as.data.frame(cbind(Years = YearsxAxis, Power = PowersYears, Method = as.factor(MethodsYears)))



#Plot the Powers
plot(DataYears$Years, DataYears$Power, col = DataYears$Method, pch = 4, xlab = "Years Between Trials", ylab = "Power", 
     main = "Plot of Power of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends in Two Different Ways")

#Put line through points
lines(Years, LinearPowersYears, col = "red")
lines(Years, StepLinearPowersYears, col = "green")

#Add lines indicating Year used in main analysis
abline(v = 5, lty = "dotted", col = "purple")

#Add lines indicating power when using Year in main analysis
abline(h = PowerLinearHistMyeloma, lty = "dotted", col = "red")
abline(h = PowerStepLinearHistMyeloma, lty = "dotted", col = "green")

#Add legend

legend(1, 0.808, legend=c("Linear Time", "Step and Linear Time"), col= c("red", "green"), lty=1, cex=0.8)


#Graph for Errors

#Then plot the type 1 Errors

YearsxAxis <- c(rep(Years, 2))
ErrorsYears <- c(LinearType1ErrorsYears, StepLinearType1ErrorsYears)
MethodsYears <- c(rep("Linear", 5),rep("StepLinear", 5))
Data2Years <- as.data.frame(cbind(Years = YearsxAxis, Type1Error = ErrorsYears, Method = as.factor(MethodsYears)))
Data2Years

plot(Data2Years$Years, Data2Years$Type1Error, col = Data2Years$Method, pch = 4, xlab = "Years Between Trials", ylab = "Type 1 Error", 
     main = "Plot of Type 1 Error of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends in Two Different Ways")

lines(Years, LinearType1ErrorsYears, col = "red")
lines(Years, StepLinearType1ErrorsYears, col = "green")


#Add lines indicating year used in main analysis
abline(v = 5, lty = "dotted", col = "purple")

#Add lines indicating power when using year in main analysis
abline(h = TypeIErrorLinearTrendMyeloma, lty = "dotted", col = "red")
abline(h = TypeIErrorStepLinearTrendMyeloma, lty = "dotted", col = "green")

#Add legend

legend(1, 0.055, legend=c("Linear Time", "Step and Linear Time"),col= c("red", "green"), lty=1, cex=0.8)




