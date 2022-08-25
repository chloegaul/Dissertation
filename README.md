
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

















# - SENSITIVITY ANALYSIS - CHANGE IN HISTORICAL LAMBDA2 - To be considered alongside Section 6.6.2


#SIMULATING NEW HISTORICAL DATA SETS
#Simulate 5 more historical data sets, each with a different value for lambda2, changing lambda1 as needed to maintain a 0.6 hazard ratio 

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

#Check that each trial individually has a power of at least 80%, and estimate monte carlo error. 

#Separate Analysis

#0.018
set.seed(1234)
SeparatePowerHist0.018 <- PowerSeparateAnalysis(100, 100, 0.03, 0.018, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.851

MCVariancePowerSeparate0.018 <-(SeparatePowerHist0.018*(1 - SeparatePowerHist0.018))/1000
MCSDPowerSeparate0.018 <- sqrt(MCVariancePowerSeparate0.018) #0.01126051

#0.019
set.seed(1234)
SeparatePowerHist0.019 <- PowerSeparateAnalysis(100, 100, 0.03166667, 0.019, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.856

#Calculate the estimated Monte Carlo Error of this power
MCVariancePowerSeparate0.019 <-(SeparatePowerHist0.019*(1 - SeparatePowerHist0.019))/1000
MCSDPowerSeparate0.019 <- sqrt(MCVariancePowerSeparate0.019) #0.01110243

#0.02
set.seed(1234)
SeparatePowerHist0.02 <- PowerSeparateAnalysis(100, 100, 0.03333333, 0.02, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.858

#Calculate the estimated Monte Carlo Error of this power
MCVariancePowerSeparate0.02 <-(SeparatePowerHist0.02*(1 - SeparatePowerHist0.02))/1000
MCSDPowerSeparate0.02 <- sqrt(MCVariancePowerSeparate0.02) #0.01103793

#0.022
set.seed(1234)
SeparatePowerHist0.022 <- PowerSeparateAnalysis(100, 100, 0.03666667, 0.022, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.877

#Calculate the estimated Monte Carlo Error of this power
MCVariancePowerSeparate0.022 <-(SeparatePowerHist0.022*(1 - SeparatePowerHist0.022))/1000
MCSDPowerSeparate0.022 <- sqrt(MCVariancePowerSeparate0.022) #0.0103861

#0.023
set.seed(1234)
SeparatePowerHist0.023 <- PowerSeparateAnalysis(100, 100, 0.03833333, 0.023, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.888

#Calculate the estimated Monte Carlo Error of this power
MCVariancePowerSeparate0.023 <-(SeparatePowerHist0.023*(1 - SeparatePowerHist0.023))/1000
MCSDPowerSeparate0.023 <- sqrt(MCVariancePowerSeparate0.023) #0.009972763

#Summary 

Histlambda0.018 <- 0.018
Histlambda0.019 <- 0.019
Histlambda0.02 <- 0.02
HistlambdaMyeloma <- 0.02100446
Histlambda0.022 <- 0.022
Histlambda0.023 <- 0.023

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)

SeparatePowersHist <- c(SeparatePowerHist0.018,SeparatePowerHist0.019, SeparatePowerHist0.02, PowerHistTrialMyeloma, SeparatePowerHist0.022, SeparatePowerHist0.023)
SeparatePowersHistMCE <- c(MCSDPowerSeparate0.018, MCSDPowerSeparate0.019, MCSDPowerSeparate0.02,MCSDSeparateHistMyeloma, MCSDPowerSeparate0.022, MCSDPowerSeparate0.023)

PowersHistorialData <- as.data.frame(cbind(lambdas, SeparatePowersHist,SeparatePowersHistMCE))






#POWER - Sensitivity Analysis - Change in lambda

#Pooled Analysis

#0.018
set.seed(1234)
PooledPowerHist0.018 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.018, 0.05, 0, 48) #0.763

#Calculate the estimated Monte Carlo Error of this power
MCVariancePooled0.018 <-(PooledPowerHist0.018*(1 - PooledPowerHist0.018))/1000
MCSDPooled0.018 <- sqrt(MCVariancePooled0.018) #0.01344734

#0.019
set.seed(1234)
PooledPowerHist0.019 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.019, 0.05, 0, 48) #0.794

#Calculate the estimated Monte Carlo Error of this power
MCVariancePooled0.019 <-(PooledPowerHist0.019*(1 - PooledPowerHist0.019))/1000
MCSDPooled0.019 <- sqrt(MCVariancePooled0.019) #0.01278921

#0.02
set.seed(1234)
PooledPowerHist0.02 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.02, 0.05, 0, 48) #0.827

#Calculate the estimated Monte Carlo Error of this power
MCVariancePooled0.02 <-(PooledPowerHist0.02*(1 - PooledPowerHist0.02))/1000
MCSDPooled0.02 <- sqrt(MCVariancePooled0.02) #0.01196123

#0.022
set.seed(1234)
PooledPowerHist0.022 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.022, 0.05, 0, 48) #0.848

#Calculate the estimated Monte Carlo Error of this power
MCVariancePooled0.022 <-(PooledPowerHist0.022*(1 - PooledPowerHist0.022))/1000
MCSDPooled0.022 <- sqrt(MCVariancePooled0.022) #0.01135324

#0.023
set.seed(1234)
PooledPowerHist0.023 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.023, 0.05, 0, 48) #0.881

#Calculate the estimated Monte Carlo Error of this power
MCVariancePooled0.023 <-(PooledPowerHist0.023*(1 - PooledPowerHist0.023))/1000
MCSDPooled0.023 <- sqrt(MCVariancePooled0.023) #0.01023909




#Step Function

#0.018
set.seed(1234)
PowerStepHist0.018 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018) #0.803

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStep0.018 <-(PowerStepHist0.018*(1 - PowerStepHist0.018))/1000
MCSDStep0.018 <- sqrt(MCVarianceStep0.018) #0.0125774

#0.019
set.seed(1234)
PowerStepHist0.019 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019) #0.802

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStep0.019 <-(PowerStepHist0.019*(1 - PowerStepHist0.019))/1000
MCSDStep0.019 <- sqrt(MCVarianceStep0.019) #0.01260143

#0.02
set.seed(1234)
PowerStepHist0.02 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02) #0.802

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStep0.02 <-(PowerStepHist0.02*(1 - PowerStepHist0.02))/1000
MCSDStep0.02 <- sqrt(MCVarianceStep0.02) #0.01260143

#0.022
set.seed(1234)
PowerStepHist0.022 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022) #0.798

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStep0.022 <-(PowerStepHist0.022*(1 - PowerStepHist0.022))/1000
MCSDStep0.022 <- sqrt(MCVarianceStep0.022) #0.0126963

#0.023
set.seed(1234)
PowerStepHist0.023 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023) #0.801

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStep0.023 <-(PowerStepHist0.023*(1 - PowerStepHist0.023))/1000
MCSDStep0.023 <- sqrt(MCVarianceStep0.023) #0.01262533




#Linear Function: 

#0.018
set.seed(1234)
PowerLinearHist0.018 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 70, 60) #0.815

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear0.018 <-(PowerLinearHist0.018*(1 - PowerLinearHist0.018))/1000
MCSDLinear0.018 <- sqrt(MCVarianceLinear0.018) #0.01227905

#0.019
set.seed(1234)
PowerLinearHist0.019 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 70, 60) #0.816

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear0.019 <-(PowerLinearHist0.019*(1 - PowerLinearHist0.019))/1000
MCSDLinear0.019 <- sqrt(MCVarianceLinear0.019) #0.01225333

#0.02
set.seed(1234)
PowerLinearHist0.02 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 70, 60) #0.821

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear0.02 <-(PowerLinearHist0.02*(1 - PowerLinearHist0.02))/1000
MCSDLinear0.02 <- sqrt(MCVarianceLinear0.02) #0.01212266

#0.022
set.seed(1234)
PowerLinearHist0.022 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 70, 60) #0.821

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear0.022 <-(PowerLinearHist0.022*(1 - PowerLinearHist0.022))/1000
MCSDLinear0.022 <- sqrt(MCVarianceLinear0.022) #0.01212266

#0.023
set.seed(123)
PowerLinearHist0.023 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 70, 60) #0.83

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear0.023 <-(PowerLinearHist0.023*(1 - PowerLinearHist0.023))/1000
MCSDLinear0.023 <- sqrt(MCVarianceLinear0.023) #0.01187855




#Step and linear function 

#0.018
set.seed(1234)
PowerStepLinearHist0.018 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 70, 60) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear0.018 <-(PowerStepLinearHist0.018*(1 - PowerStepLinearHist0.018))/1000
MCSDStepLinear0.018 <- sqrt(MCVarianceStepLinear0.018) #0.01255325

#0.019
set.seed(1234)
PowerStepLinearHist0.019 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 70, 60) #0.803

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear0.019 <-(PowerStepLinearHist0.019*(1 - PowerStepLinearHist0.019))/1000
MCSDStepLinear0.019 <- sqrt(MCVarianceStepLinear0.019) #0.0125774

#0.02
set.seed(1234)
PowerStepLinearHist0.02 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 70, 60) #0.805

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear0.02 <-(PowerStepLinearHist0.02*(1 - PowerStepLinearHist0.02))/1000
MCSDStepLinear0.02 <- sqrt(MCVarianceStepLinear0.02) #0.01252897

#0.022
set.seed(1234)
PowerStepLinearHist0.022 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 70, 60) #0.802

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear0.022 <-(PowerStepLinearHist0.022*(1 - PowerStepLinearHist0.022))/1000
MCSDStepLinear0.022 <- sqrt(MCVarianceStepLinear0.022) #0.01260143

#0.023
set.seed(1234)
PowerStepLinearHist0.023 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 70, 60) #0.802

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear0.023 <-(PowerStepLinearHist0.023*(1 - PowerStepLinearHist0.023))/1000
MCSDStepLinear0.023 <- sqrt(MCVarianceStepLinear0.023) #0.01260143



#Summary of Power - Sensitivity Analysis

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)
PooledPowers <- c(PooledPowerHist0.018, PooledPowerHist0.019, PooledPowerHist0.02, PooledPowerHistMyeloma, PooledPowerHist0.022, PooledPowerHist0.023)
StepTimePowers <- c(PowerStepHist0.018, PowerStepHist0.019, PowerStepHist0.02, PowerStepHistMyeloma, PowerStepHist0.022, PowerStepHist0.023)
LinearTimePowers <- c(PowerLinearHist0.018, PowerLinearHist0.019, PowerLinearHist0.02, PowerLinearHistMyeloma, PowerLinearHist0.022, PowerLinearHist0.023)
StepLinearTimePowers <- c(PowerStepLinearHist0.018, PowerStepLinearHist0.019, PowerStepLinearHist0.02, PowerStepLinearHistMyeloma, PowerStepLinearHist0.022, PowerStepLinearHist0.023)
AllPowersLambda <- as.data.frame(cbind(lambdas, PooledPowers, StepTimePowers, LinearTimePowers, StepLinearTimePowers)) 



#Graph of Power - Sensitivity Analysis -  change in historical lambda2 

#Pooled

par(mfrow = c(2,2))


plot(lambdas, PooledPowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of the Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Using Pooled Analysis")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PooledPowerHistMyeloma, lty = "dotted", col = "purple")

#Step Function 

plot(lambdas, StepTimePowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of the Current Trial when Borrowing
     Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PowerStepHistMyeloma, lty = "dotted", col = "purple")


#Linear Function 

plot(lambdas, LinearTimePowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of the Current Trial when Borrowing
     Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PowerLinearHistMyeloma, lty = "dotted", col = "purple")

#Step and Linear function

plot(lambdas, StepLinearTimePowers, type = "b", xlab = "Value of Historical Lambda2", ylab = "Power", main = "Plot of Power of Current Trial when Borrowing Historical 
Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step and Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = PowerStepLinearHistMyeloma, lty = "dotted", col = "purple")



#TYPE 1 ERROR RATE - Sensitivity Analysis - Change in Historical lambda2


#We also need to take into account the Type I error rates 

#Pooled Analysis

#0.018
set.seed(1234)
Type1ErrorPooled0.018 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.018, 0.05, 0, 48) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Pooled0.018 <-(Type1ErrorPooled0.018*(1 - Type1ErrorPooled0.018))/1000
MCSDType1Pooled0.018 <- sqrt(MCVarianceType1Pooled0.018) #0.007331507

#0.019
set.seed(1234)
Type1ErrorPooled0.019 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.019, 0.05, 0, 48) #0.053

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Pooled0.019 <-(Type1ErrorPooled0.019*(1 - Type1ErrorPooled0.019))/1000
MCSDType1Pooled0.019 <- sqrt(MCVarianceType1Pooled0.019) #0.007084561

#0.02
set.seed(1234)
Type1ErrorPooled0.02 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.02, 0.05, 0, 48) #0.049

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Pooled0.02 <-(Type1ErrorPooled0.02*(1 - Type1ErrorPooled0.02))/1000
MCSDType1Pooled0.02 <- sqrt(MCVarianceType1Pooled0.02) #0.006826346

#0.022
set.seed(1234)
Type1ErrorPooled0.022 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.022, 0.05, 0, 48) #0.048

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Pooled0.022 <-(Type1ErrorPooled0.022*(1 - Type1ErrorPooled0.022))/1000
MCSDType1Pooled0.022 <- sqrt(MCVarianceType1Pooled0.022) #0.006759882

#0.023
set.seed(1234)
Type1ErrorPooled0.023 <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrial0.023, 0.05, 0, 48) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Pooled0.023 <-(Type1ErrorPooled0.023*(1 - Type1ErrorPooled0.023))/1000
MCSDType1Pooled0.023 <- sqrt(MCVarianceType1Pooled0.023) #0.007331507
 




#Modelling time trends as a step function 

#0.018
set.seed(1234)
Type1ErrorTimeTrend0.018 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Step0.018 <-(Type1ErrorTimeTrend0.018*(1 - Type1ErrorTimeTrend0.018))/1000
MCSDType1Step0.018 <- sqrt(MCVarianceType1Step0.018) #0.007331507

#0.019
set.seed(1234)
Type1ErrorTimeTrend0.019 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Step0.019 <-(Type1ErrorTimeTrend0.019*(1 - Type1ErrorTimeTrend0.019))/1000
MCSDType1Step0.019 <- sqrt(MCVarianceType1Step0.019) #0.007331507

#0.02
set.seed(1234)
Type1ErrorTimeTrend0.02 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Step0.02 <-(Type1ErrorTimeTrend0.02*(1 - Type1ErrorTimeTrend0.02))/1000
MCSDType1Step0.02 <- sqrt(MCVarianceType1Step0.02) #0.007331507

#0.022
set.seed(1234)
Type1ErrorTimeTrend0.022 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Step0.022 <-(Type1ErrorTimeTrend0.022*(1 - Type1ErrorTimeTrend0.022))/1000
MCSDType1Step0.022 <- sqrt(MCVarianceType1Step0.022) #0.007331507

#0.023
set.seed(1234)
Type1ErrorTimeTrend0.023 <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Step0.023 <-(Type1ErrorTimeTrend0.023*(1 - Type1ErrorTimeTrend0.023))/1000
MCSDType1Step0.023 <- sqrt(MCVarianceType1Step0.023) #0.007331507





#Modelling time trends as a linear function 

#0.018
set.seed(1234)
Type1ErrorLinearTrend0.018 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 72, 60) #0.056

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear0.018 <-(Type1ErrorLinearTrend0.018*(1 - Type1ErrorLinearTrend0.018))/1000
MCSDType1Linear0.018 <- sqrt(MCVarianceType1Linear0.018) #0.007270763

#0.019
set.seed(1234)
Type1ErrorLinearTrend0.019 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 72, 60) #0.056

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear0.019 <-(Type1ErrorLinearTrend0.019*(1 - Type1ErrorLinearTrend0.019))/1000
MCSDType1Linear0.019 <- sqrt(MCVarianceType1Linear0.019) #0.007270763

#0.02
set.seed(1234)
Type1ErrorLinearTrend0.02 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 72, 60) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear0.02 <-(Type1ErrorLinearTrend0.02*(1 - Type1ErrorLinearTrend0.02))/1000
MCSDType1Linear0.02 <- sqrt(MCVarianceType1Linear0.02) #0.007331507

#0.022
set.seed(1234)
Type1ErrorLinearTrend0.022 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 72, 60) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear0.022 <-(Type1ErrorLinearTrend0.022*(1 - Type1ErrorLinearTrend0.022))/1000
MCSDType1Linear0.022 <- sqrt(MCVarianceType1Linear0.022) #0.007331507

#0.023
set.seed(1234)
Type1ErrorLinearTrend0.023 <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 72, 60) #0.059

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear0.023 <-(Type1ErrorLinearTrend0.023*(1 - Type1ErrorLinearTrend0.023))/1000
MCSDType1Linear0.023 <- sqrt(MCVarianceType1Linear0.023) #0.007451107





#Modelling time trends as a step and linear function

#0.018
set.seed(1234)
Type1ErrorStepLinearTrend0.018 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.018, 72, 60) #0.055

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear0.018 <-(Type1ErrorStepLinearTrend0.018*(1 - Type1ErrorStepLinearTrend0.018))/1000
MCSDType1StepLinear0.018 <- sqrt(MCVarianceType1StepLinear0.018) #0.007209369

#0.019
set.seed(1234)
Type1ErrorStepLinearTrend0.019 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.019, 72, 60) #0.056

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear0.019 <-(Type1ErrorStepLinearTrend0.019*(1 - Type1ErrorStepLinearTrend0.019))/1000
MCSDType1StepLinear0.019 <- sqrt(MCVarianceType1StepLinear0.019) #0.007270763

#0.02
set.seed(1234)
Type1ErrorStepLinearTrend0.02 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.02, 72, 60) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear0.02 <-(Type1ErrorStepLinearTrend0.02*(1 - Type1ErrorStepLinearTrend0.02))/1000
MCSDType1StepLinear0.02 <- sqrt(MCVarianceType1StepLinear0.02) #0.007147307

#0.022
set.seed(1234)
Type1ErrorStepLinearTrend0.022 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.022, 72, 60) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear0.022 <-(Type1ErrorStepLinearTrend0.022*(1 - Type1ErrorStepLinearTrend0.022))/1000
MCSDType1StepLinear0.022 <- sqrt(MCVarianceType1StepLinear0.022) #0.007147307

#0.023
set.seed(1234)
Type1ErrorStepLinearTrend0.023 <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrial0.023, 72, 60) #0.54

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear0.023 <-(Type1ErrorStepLinearTrend0.023*(1 - Type1ErrorStepLinearTrend0.023))/1000
MCSDType1StepLinear0.023 <- sqrt(MCVarianceType1StepLinear0.023) #0.007147307



#Summary 

Histlambda0.018 <- 0.018
Histlambda0.019 <- 0.019
Histlambda0.02 <- 0.02
HistlambdaMyeloma <- 0.02100446
Histlambda0.022 <- 0.022
Histlambda0.023 <- 0.023

lambdas <- c(Histlambda0.018,Histlambda0.019, Histlambda0.02, HistlambdaMyeloma,Histlambda0.022, Histlambda0.023)

PooledType1Errors <- c(Type1ErrorPooled0.018, Type1ErrorPooled0.019, Type1ErrorPooled0.02, Type1ErrorPooledMyeloma,Type1ErrorPooled0.022,Type1ErrorPooled0.023)

StepTimeType1Errors <- c(Type1ErrorTimeTrend0.018, Type1ErrorTimeTrend0.019 , Type1ErrorTimeTrend0.02, Type1ErrorTimeTrendMyeloma, Type1ErrorTimeTrend0.022, Type1ErrorTimeTrend0.023)

LinearTimeType1Errors <- c(Type1ErrorLinearTrend0.018, Type1ErrorLinearTrend0.019, Type1ErrorLinearTrend0.02, Type1ErrorLinearTrendMyeloma, Type1ErrorLinearTrend0.022, Type1ErrorLinearTrend0.023)

StepLinearTimeType1Errors <- c(Type1ErrorStepLinearTrend0.018, Type1ErrorStepLinearTrend0.019, Type1ErrorStepLinearTrend0.02, Type1ErrorStepLinearTrendMyeloma,Type1ErrorStepLinearTrend0.022, Type1ErrorStepLinearTrend0.023)

AllErrorsLambda <- as.data.frame(cbind(lambdas, PooledType1Errors, StepTimeType1Errors, LinearTimeType1Errors, StepLinearTimeType1Errors)) 


#Graphs of Type 1 Errors - Sensitivity Analysis - Change in Historical Lambda2

par(mfrow = c(2,2))

#Pooled

plot(lambdas, PooledType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Using Pooled Analysis")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorPooledMyeloma, lty = "dotted", col = "purple")

#Step function:

plot(lambdas, StepTimeType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorTimeTrendMyeloma, lty = "dotted", col = "purple")

#Linear function:

plot(lambdas, LinearTimeType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorLinearTrendMyeloma, lty = "dotted", col = "purple")

#Step and Linear function

plot(lambdas, StepLinearTimeType1Errors, type = "b", xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of Type 1 Error of Current Trial when Borrowing 
Historical Data Sets of Varying Lambda2, 
     Modelling Time Trends with a Step and Linear Function")

abline(v = HistlambdaMyeloma, lty = "dotted", col = "purple")
abline(h = Type1ErrorStepLinearTrendMyeloma, lty = "dotted", col = "purple")



#Graphs of power and type 1 error for all methods across differing values of historical lambda2

AllDataHistLambda <- as.data.frame(cbind(AllPowersLambda, AllErrorsLambda))

par(mfrow = c(1,2))


#Plotting Power
Historicallambdas <- c(rep(lambdas, 4)) #x values of the points to plot

PowersLambdas <- c(PooledPowers, StepTimePowers, LinearTimePowers, StepLinearTimePowers) #y values of the points to plot

MethodsLambdas <- c(rep("Pooling", 5), rep("Step", 5),rep("Linear", 5),rep("StepLinear", 5)) #Methods to colour points by

PowerDataLambdas <- as.data.frame(cbind(HistoricalLambda = Historicallambdas, Power = PowersLambdas, Method = as.factor(MethodsLambdas)))

plot(DataLambdas$HistoricalLambda, DataLambdas$Power, col = DataLambdas$Method, pch = 4, xlab = "Value of Historical Lambda2", ylab = "Power", 
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

legend(0.019, 0.78, legend=c("Pooled", "Step Time", "Linear Time", "Step and Linear Time"),col= c("red", "green", "black", "blue"), lty=1, cex=0.8)





#Plotting Type 1 Errors

Historicallambdas <- c(rep(lambdas, 4)) #x values of points to plot

ErrorsLambdas <- c(PooledType1Errors, StepTimeType1Errors, LinearTimeType1Errors, StepLinearTimeType1Errors) 

MethodsLambdas <- c(rep("Pooling", 5), rep("Step", 5),rep("Linear", 5),rep("StepLinear", 5)) #Methods to colour points by

Data2Lambdas <- as.data.frame(cbind(HistoricalLambda = Historicallambdas, Type1Error = ErrorsLambdas, Method = as.factor(MethodsLambdas)))

plot(Data2Lambdas$HistoricalLambda, Data2Lambdas$Type1Error, col = Data2Lambdas$Method, pch = 4, 
     xlab = "Value of Historical Lambda2", ylab = "Type 1 Error", main = "Plot of the Type 1 Error of the Current Trial when Borrowing 
     Historical Data Using Various Methods, from Historical 
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

legend(0.018, 0.0593, legend=c("Pooled", "Step Time", "Linear Time", "Step and Linear Time"),col= c("red", "green", "black", "blue"), lty=1, cex=0.8)







##POWER - Sensitivity Analysis - Change in time between trials 


#Model time trends using a linear function

#1 Year
set.seed(1234)
PowerLinear1Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.824

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear1year <-(PowerLinear1Year*(1 - PowerLinear1Year))/1000
MCSDLinear1year <- sqrt(MCVarianceLinear1year) #0.01204259

#2 Years
set.seed(1234)
PowerLinear2Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.821

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear2year <-(PowerLinear2Year*(1 - PowerLinear2Year))/1000
MCSDLinear2year <- sqrt(MCVarianceLinear2year) #0.01212266

#3 Years
set.seed(1234)
PowerLinear3Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.82

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear3year <-(PowerLinear3Year*(1 - PowerLinear3Year))/1000
MCSDLinear3year <- sqrt(MCVarianceLinear3year) #0.01214907

#4 Years
set.seed(1234)
PowerLinear4Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.82

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear4year <-(PowerLinear4Year*(1 - PowerLinear4Year))/1000
MCSDLinear4year <- sqrt(MCVarianceLinear4year) #0.01214907

#6 Years
set.seed(1234)
PowerLinear6Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.82

#Calculate the estimated Monte Carlo Error of this power
MCVarianceLinear6year <-(PowerLinear6Year*(1 - PowerLinear6Year))/1000
MCSDLinear6year <- sqrt(MCVarianceLinear6year) #0.01214907



#Model time trends using a step and linear function

#1 Year
set.seed(1234)
PowerStepLinear1Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear1year <-(PowerStepLinear1Year*(1 - PowerStepLinear1Year))/1000
MCSDStepLinear1year <- sqrt(MCVarianceStepLinear1year) #0.01255325

#2 Years
set.seed(1234)
PowerStepLinear2Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear2year <-(PowerStepLinear2Year*(1 - PowerStepLinear2Year))/1000
MCSDStepLinear2year <- sqrt(MCVarianceStepLinear2year) #0.01255325

#3 Years
set.seed(1234)
PowerStepLinear3Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear3year <-(PowerStepLinear3Year*(1 - PowerStepLinear3Year))/1000
MCSDStepLinear3year <- sqrt(MCVarianceStepLinear3year) #0.01255325

#4 Years
set.seed(1234)
PowerStepLinear4Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear4year <-(PowerStepLinear4Year*(1 - PowerStepLinear4Year))/1000
MCSDStepLinear4year <- sqrt(MCVarianceStepLinear4year) #0.01255325

#6 Years
set.seed(1234)
PowerStepLinear6Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.804

#Calculate the estimated Monte Carlo Error of this power
MCVarianceStepLinear6year <-(PowerStepLinear6Year*(1 - PowerStepLinear6Year))/1000
MCSDStepLinear6year <- sqrt(MCVarianceStepLinear6year) #0.01255325


#Summary

Years <- c(1, 2, 3, 4, 5, 6)
LinearPowersYears <- c(PowerLinear1Year, PowerLinear2Year, PowerLinear3Year, PowerLinear4Year, PowerLinearHistMyeloma, PowerLinear6Year)
StepLinearPowersYears <- c(PowerStepLinear1Year, PowerStepLinear2Year, PowerStepLinear3Year, PowerStepLinear4Year, PowerStepLinearHistMyeloma, PowerStepLinear6Year)

#Graphs

par(mfrow = c(1,2))

#Linear

plot(Years, LinearPowersYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Power", main = "Plot of Power of the Current Trial When Borrowing
   Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = PowerLinearHistMyeloma, lty = "dotted", col = "purple")

#Step Linear

plot(Years, StepLinearPowersYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Power", main = "Plot of Power of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear and Step Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = PowerStepLinearHistMyeloma, lty = "dotted", col = "purple")



##TYPE 1 ERROR - Sensitivity Analysis - Change in time between trials 

#Modeling time trends using a Linear 

#1 Year
set.seed(1234)
Type1ErrorLinearTimeTrend1Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.059

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear1year <-(Type1ErrorLinearTimeTrend1Year*(1 - Type1ErrorLinearTimeTrend1Year))/1000
MCSDType1Linear1year <- sqrt(MCVarianceType1Linear1year) #0.007451107

#2 Years
set.seed(1234)
Type1ErrorLinearTimeTrend2Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.058

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear2year <-(Type1ErrorLinearTimeTrend2Year*(1 - Type1ErrorLinearTimeTrend2Year))/1000
MCSDType1Linear2year <- sqrt(MCVarianceType1Linear2year) #0.007391617

#3 Years
set.seed(1234)
Type1ErrorLinearTimeTrend3Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.059

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear3year <-(Type1ErrorLinearTimeTrend3Year*(1 - Type1ErrorLinearTimeTrend3Year))/1000
MCSDType1Linear3year <- sqrt(MCVarianceType1Linear3year) #0.007451107

#4 Years
set.seed(1234)
Type1ErrorLinearTimeTrend4Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.059

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear4year <-(Type1ErrorLinearTimeTrend4Year*(1 - Type1ErrorLinearTimeTrend4Year))/1000
MCSDType1Linear4year <- sqrt(MCVarianceType1Linear4year) #0.007451107

#6 Years
set.seed(1234)
Type1ErrorLinearTimeTrend6Year <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.057

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1Linear6year <-(Type1ErrorLinearTimeTrend6Year*(1 - Type1ErrorLinearTimeTrend6Year))/1000
MCSDType1Linear6year <- sqrt(MCVarianceType1Linear6year) #0.007331507



#Modelling time trends using a step and linear function

#1 Year
set.seed(1234)
Type1ErrorStepLinearTimeTrend1Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 12) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear1year <-(Type1ErrorStepLinearTimeTrend1Year*(1 - Type1ErrorStepLinearTimeTrend1Year))/1000
MCSDType1StepLinear1year <- sqrt(MCVarianceType1StepLinear1year) #0.007147307

#2 Years
set.seed(1234)
Type1ErrorStepLinearTimeTrend2Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 24) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear2year <-(Type1ErrorStepLinearTimeTrend2Year*(1 - Type1ErrorStepLinearTimeTrend2Year))/1000
MCSDType1StepLinear2year <- sqrt(MCVarianceType1StepLinear2year) #0.007147307

#3 Years
set.seed(1234)
Type1ErrorStepLinearTimeTrend3Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 36) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear3year <-(Type1ErrorStepLinearTimeTrend3Year*(1 - Type1ErrorStepLinearTimeTrend3Year))/1000
MCSDType1StepLinear3year <- sqrt(MCVarianceType1StepLinear3year) #0.007147307

#4 Years
set.seed(1234)
Type1ErrorStepLinearTimeTrend4Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 48) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear4year <-(Type1ErrorStepLinearTimeTrend4Year*(1 - Type1ErrorStepLinearTimeTrend4Year))/1000
MCSDType1StepLinear4year <- sqrt(MCVarianceType1StepLinear4year) #0.007147307

#6 Years
set.seed(1234)
Type1ErrorStepLinearTimeTrend6Year <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 72) #0.054

#Calculate the estimated Monte Carlo Error of this power
MCVarianceType1StepLinear6year <-(Type1ErrorStepLinearTimeTrend6Year*(1 - Type1ErrorStepLinearTimeTrend6Year))/1000
MCSDType1StepLinear6year <- sqrt(MCVarianceType1StepLinear6year) # 0.007147307



#Summary

Years <- c(1, 2, 3, 4, 5, 6)
LinearType1ErrorsYears <- c(Type1ErrorLinearTimeTrend1Year, Type1ErrorLinearTimeTrend2Year, Type1ErrorLinearTimeTrend3Year, Type1ErrorLinearTimeTrend4Year, TypeIErrorLinearTrendMyeloma, Type1ErrorLinearTimeTrend6Year)
StepLinearType1ErrorsYears <- c(Type1ErrorStepLinearTimeTrend1Year, Type1ErrorStepLinearTimeTrend2Year, Type1ErrorStepLinearTimeTrend3Year, Type1ErrorStepLinearTimeTrend4Year, TypeIErrorStepLinearTrendMyeloma,
                                Type1ErrorStepLinearTimeTrend6Year)

#Graph the Type 1 Errors

par(mfrow = c(1,2))

#Linear 

plot(Years, LinearType1ErrorsYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Type 1 Error", main = "Plot of Type 1 Error of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = TypeIErrorLinearTrendMyeloma, lty = "dotted", col = "purple")

#Step Linear

plot(Years, StepLinearType1ErrorsYears, type = "b", xlab = "Years Between Historical and Current Trial", ylab = "Type 1 Error", main = "Plot of Type 1 Error of the Current Trial When Borrowing
     Historical Data Occuring Various Years Ago, 
     Modelling Time Trends with a Linear and Step Function")

abline(v = 5, lty = "dotted", col = "purple")
abline(h = TypeIErrorStepLinearTrendMyeloma, lty = "dotted", col = "purple")





#Power and type 1 error for all times between trials using all methods

#Summary

Years <- c(1, 2, 3, 4, 5, 6)

LinearPowersYears <- c(PowerLinear1Year, PowerLinear2Year, PowerLinear3Year, PowerLinear4Year, PowerLinearHistMyeloma, PowerLinear6Year)

StepLinearPowersYears <- c(PowerStepLinear1Year, PowerStepLinear2Year, PowerStepLinear3Year, PowerStepLinear4Year, PowerStepLinearHistMyeloma, PowerStepLinear6Year)

LinearType1ErrorsYears <- c(Type1ErrorLinearTimeTrend1Year, Type1ErrorLinearTimeTrend2Year, Type1ErrorLinearTimeTrend3Year, Type1ErrorLinearTimeTrend4Year, TypeIErrorLinearTrendMyeloma, Type1ErrorLinearTimeTrend6Year)

StepLinearType1ErrorsYears <- c(Type1ErrorStepLinearTimeTrend1Year, Type1ErrorStepLinearTimeTrend2Year, Type1ErrorStepLinearTimeTrend3Year, Type1ErrorStepLinearTimeTrend4Year, TypeIErrorStepLinearTrendMyeloma, Type1ErrorStepLinearTimeTrend6Year)

DataHistTime <- as.data.frame(cbind(LinearPowersYears,StepLinearPowersYears, LinearType1ErrorsYears, StepLinearType1ErrorsYears))

#Graphs

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




