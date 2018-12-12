## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----libraries, message=FALSE, warning=FALSE, results=FALSE--------------
#Find out what your working directory currently is
getwd()

#if necessary, set your working directory to whatever folder contains the data and scripts using the 'setwd' function
#notice that "/" is ued to separate the directories NOT "\" so if you copy and paste the address make sure to change those!

#setwd("C:/desktop/whatever/whatever")


#install the required package. 


#install.packages("bcpa")

#load the package
require(bcpa)



## ----getdata, message=FALSE, warning=FALSE-------------------------------


#I normally start by reading in the complete data set and then selecting one animal
#Read in "All Animals" data 

Data<-read.csv("Cow  BCPA.csv", na.strings="NA")

#select individual animals by ID
#In this case the data set has only one animal(collar number 129) but it is good to see how it works
mydata<-Data[Data$ID==129,]

#View the data to see if it looks right
head(mydata)

#Slect the columns that are of use for this analysis
mydata<-mydata[,c("timestamp","Lat","Lon","ID")]


#Rename columns in dataframe, BCPA requires the column names below
names(mydata)<-c("Time", "Y", "X","ID" )

#Format "Time" as a POSIXct object to make it easier to deal with as a time series 
#the format given here must match exactly the format in the original data

Da<-strptime(mydata$Time,format='%d/%m/%Y %H:%M', tz="GMT")
Da<-as.POSIXct(Da)
mydata$Time<-Da

#Check that the "Time" column looks okay and the formatting and time zone are correct
head(mydata$Time)

#Read in some plotting functions that I have tweaked for our purposes
source("BCPA_Functions.R")

## ----timeseries----------------------------------------------------------
AnimalVT<-GetVT(mydata, units="hour")

head(AnimalVT)

## ----WindowSweep---------------------------------------------------------


system.time(Animal.ws<-WindowSweep(AnimalVT,"V*cos(Theta)",windowsize=12,  progress=FALSE,  plotme=FALSE, K=2, tau=FALSE))




## ----plotBCPA------------------------------------------------------------
#plot outputs as "flat", ie.: without any smoothing between each point. Time needs to be a POSIXct object
#clusterwidth tells the function to group change points that were detected within a certain distance from each other. default =1
#"rho.where" and "mu.where" are for placing the legend if you want one.
bcpa.plot(Animal.ws, type="flat",ylab="Persistence Velocity",xlab="Time", clusterwidth=1, rho.where='topright', mu.where="topright")


## ----PathPlotPersist-----------------------------------------------------
PathPlot.Persistence(mydata, Animal.ws, plotlegend=TRUE)

## ----changepoints--------------------------------------------------------
changes<-ChangePointSummary(Animal.ws, tau=FALSE)
#how many change points are there?
length(changes$breaks$size)
changes


## ----phases--------------------------------------------------------------

phases<-changes$phases

# the first phase has nothing to compare to
Start<-9999

#create a column showing how 'mu.hat' has changed between each phase
mu.diff<-  diff(changes$phases$mu.hat, lag=1)
                                     
mu.diff<-c(Start, mu.diff)

phases<- cbind(phases, mu.diff)

#do the same for 's.hat'
s.diff<- diff(phases$s.hat, lag=1)
                                   
s.diff<-c(Start, s.diff)

phases<- cbind(phases, s.diff)

#and finally for 'rho'

rho.diff<- diff(phases$s.hat, lag=1)
                                   
rho.diff<- c(Start, rho.diff)

phases<- cbind(phases, rho.diff)

#check the data frame
head(phases)


## ----changes-------------------------------------------------------------
breaks<-changes$breaks
breaks<-breaks[,2:5]
Start<-data.frame(middle=9999, size=9999,modelmode=9999,middle.POSIX=AnimalVT$T.POSIX[1])


breaks<-rbind(Start,breaks)


Breaks_and_Phases<-cbind(phases, breaks)
Breaks_and_Phases

## ----classify------------------------------------------------------------
phases$Behaviour<-0                        
phases$Behaviour<-ifelse(phases$s.diff==9999 & phases$mu.diff==9999 & phases$rho.diff==9999,"Start of Day",
                            ifelse(phases$mu.diff> (2*mean(phases$mu.hat)),"Direct Movement" ,     
                                   ifelse(phases$s.diff>0.002  & phases$mu.diff>0.0001, "Foraging and Moving", 
                                          ifelse(phases$s.hat>0.001  & phases$mu.hat>0 , "Slow Foraging",
                                                 
                                                 ifelse(phases$s.hat>0.001  & phases$mu.diff<0 , "Slow Foraging",    
                                                        
                                                        ifelse(phases$mu.hat<0.001 & phases$s.hat<0.001, "Sedentary","Unclassified"))))))


head(phases$Behaviour,20)   
phases$Behaviour<-as.character(phases$Behaviour)
phases$Behaviour[is.na(phases$Behaviour)]<-"NA"

#a further series of classifiers for those phases at the start of each day
#detecting changes is not possible because we don't know what came before
#so we classify these based on absolute values of the three parameters

phases$Behaviour[phases$Behaviour=="Start of Day"]<-ifelse(phases$mu.hat[phases$Behaviour=="Start of Day"] >(2*mean(phases$mu.hat)) ,"Direct Movement",
                                                                 ifelse(phases$s.hat[phases$Behaviour=="Start of Day"]>0.001  & phases$mu.hat[phases$Behaviour=="Start of Day"]>0 , "Slow Foraging",
                                                                        ifelse(phases$mu.hat[phases$Behaviour=="Start of Day"] >(1*mean(phases$mu.hat)) &  phases$s.hat[phases$Behaviour=="Start of Day"]<mean(phases$s.hat), "Foraging and Moving" ,  
                                                                               
                                                                               ifelse(phases$mu.hat[phases$Behaviour=="Start of Day"]<0.001 & phases$s.hat[phases$Behaviour=="Start of Day"]<0.001, "Sedentary", "Unclassified"))))    

#add the timestamps for each phase form the 'breaks and phases' object above
phases$Time<-Breaks_and_Phases$middle.POSIX

head(phases,20)


## ----phase_lookup--------------------------------------------------------
#A series of "if else" statements to select the correct behavioural phase and add it to the main data set                      
mydata$Behaviour<-0                        
mydata$Behaviour<-ifelse(mydata$Time<=phases$Time[1],phases$Behaviour[1],
                            ifelse(mydata$Time<=phases$Time[3]& mydata$Time >phases$Time[1],phases$Behaviour[2] ,     
                                   ifelse(mydata$Time<=phases$Time[4]& mydata$Time > phases$Time[3],phases$Behaviour[3], 
                                          ifelse(mydata$Time<=phases$Time[5]& mydata$Time >phases$Time[4],phases$Behaviour[4],
                                                 ifelse(mydata$Time>=phases$Time[5],phases$Behaviour[5], "Unclassified")))))

head(mydata)

#take the "Behaviour" column we just created and append it to the raw data from the begining

Data$Behaviour<-mydata$Behaviour

#view the results
head(Data)

#write out the results as a CSV file
write.csv(Data, "A Day In The Life Of A Cow.csv", row.names=FALSE)

