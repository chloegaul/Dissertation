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


#Main Data analysis without variations:

#Generate a Historical Trial 

#Since the current trial has a control lambda of 0.02100446, we want to test what happens to 
#The power when historical lambda for the same treatment is the same, lower and higher. 

#So the normal historical data set we have generated as:
set.seed(1234)
HistoricalTrialMyeloma <- CompleteTrialSimulation(100, 0.03500743, 72, 0, 100, 0.02100446, 1, 0, 36, 0.05)

#100 in each group because this ensures at least 80% power for all variations of lambda2 that we will try.
#It also just seems like that of a reasonably sized trial. 
#lambda2 = 0.02100446, same as Myeloma Trial
#So this represents what will happen under ideal circumstances. 
#Look for hazard ratio of 0.6, so this makes lambda1 = 0.03500743
#Make the historical trial about 35% of the length of Myeloma, which makes 72 months total, 36 of which are for accrual period. 
#Same drop out rate. 
#We're assuming that the historical trial took place 5 years before, which is important for the last 2 functions. 

#Let us check that the historical data set is powered in its own right
set.seed(1234)
PowerHistTrialMyeloma <- PowerSeparateAnalysis(100, 100, 0.03500743, 0.02100446, 1000, 72, 0.05, 1, 2, 0.05, 0, 36) #0.873

#So yes the trial we use is powered on its own. 


#Separate Analysis of current data set 

set.seed(1234)
SeparatePowerHistMyeloma <- PowerSeparateAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48) #0.8 
#So without incorpoating any historical data, using a cox proportional hazards model to test for this hazard ratio gives us a power of 0.8


#Now we see what it does when we incorporate some of the historical data in various ways

#Pooled Analysis
set.seed(1234)
PooledPowerHistMyeloma <- PowerPooledAnalysis(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrialMyeloma, 0.05, 0, 48) #0.84

#This increases the power, which is to be expected since we're essentially just using more data from a data set with the same underlying parameter 
#And it is therefore just like we have a higher sample size, which would mean a higher power. 


#Modelling time trend as a step 
set.seed(1234)
PowerStepHistMyeloma <- PowerStepTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma) #0.801
#Doesn't seem to have made much change. Ever so slight increase in power. 


#Modelling time trend as linear
set.seed(1234)
PowerLinearHistMyeloma <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 70, 60) #0.821
#We see an increase in power


#Modelling time trend as linear and step
set.seed(1234)
PowerStepLinearHistMyeloma <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.0165035, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 70, 60) #0.804
#We see no change in power. 


#Type 1 Errors

set.seed(1234)
Type1ErrorSeparateMyeloma <- PowerSeparateAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48) #0.057

set.seed(1234)
Type1ErrorPooledMyeloma <- PowerPooledAnalysis(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, HistoricalTrial = HistoricalTrialMyeloma, 0.05, 0, 48) #0.05

set.seed(1234)
Type1ErrorTimeTrendMyeloma <- PowerStepTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma) #0.057

set.seed(1234)
Type1ErrorLinearTrendMyeloma <- PowerLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 60) #0.056

set.seed(1234)
Type1ErrorStepLinearTrendMyeloma <- PowerStepLinearTimeTrend(394, 394, 0.02100446, 0.02100446, 1000, 96, 0.05, 1, 2, 0.05, 0, 48, HistoricalTrial = HistoricalTrialMyeloma, 72, 60) #0.054





#Summary

Methods <- c("Separate", "Pooling", "Step Time Trend", "Linear Time Trend", "Step Linear Time Trend")
MainAnalysisPowers <- c(SeparatePowerHistMyeloma, PooledPowerHistMyeloma, PowerStepHistMyeloma, PowerLinearHistMyeloma, PowerStepLinearHistMyeloma)
MainAnalysisType1Errors <- c(Type1ErrorSeparateMyeloma, Type1ErrorPooledMyeloma,Type1ErrorTimeTrendMyeloma, Type1ErrorLinearTrendMyeloma, Type1ErrorStepLinearTrendMyeloma)

MainAnalysisData <- as.data.frame(cbind(Methods, MainAnalysisPowers, MainAnalysisType1Errors))
MainAnalysisData

#Graph

par(mfrow = c(1,1))

plot(MainAnalysisPowers, data = MainAnalysisData, xlab = "Method of Historical Borrowing", ylab = "Power", Main = "Graph Showing Change in Power for Methods of Borrowing the HistoricalTrialMyeloma Data Set", 
     type = "b", xaxt='n', main = "Change in Power when Borrowing A Historical 
     Data Set with The Same Underlying Parameter 
     Using Different Methods")
axis(1, at = 1:nrow(MainAnalysisData), labels = Methods)

abline(h = SeparatePowerHistMyeloma, col = "blue")

legend(x= "topright", legend= "SeparateAnalysis", lty = 1, col="blue")




#Plotting the Type 1 Errors

plot(MainAnalysisType1Errors, data = MainAnalysisData, xlab = "Method of Historical Borrowing", ylab = "Power", 
     type = "b", xaxt='n', main = "Change in Type 1 Error when Borrowing A Historical Data Set with The Same Underlying Parameter Using Different Methods")
axis(1, at = 1:nrow(MainAnalysisData), labels = Methods)

abline(h = Type1ErrorSeparateMyeloma, col = "blue")

legend(x= "bottomright", legend= "SeparateAnalysis", lty = 1, col="blue")


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




