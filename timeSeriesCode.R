install.packages("forecast")
install.packages("zoo")
install.packages("tidyverse")
install.packages("xts")
install.packages("ggplot2")
install.packages("magrittr") 
library(magrittr)
library(tidyverse)
library(forecast)
library(zoo)
library(xts)
library(ggplot2)
#Naive Model
df <- read.table("C:/Users/c1826277/OneDrive - Cardiff University/Desktop/Time series and forecasting/TimeSeriesCourseworkData21_22.csv", sep = ",", header = T)
df <- na.omit(df)
df$Date <- as.Date(df$Date, format = "%d/%m/%Y")
firstYear <- df[1:275,]
secondYear <- df[276:641,]
thirdYear <- df[642:1006,]
fourthYear <- df[1007:1371,]
fifthYear <- df[1372:1461,]
allYear <- df[1:1461,]

scatter.smooth(x=df$Date, y=df$Number.of.patients, main="Date ~ Number of Patients")

#Calculating MAPE
forecastFirstYear <- c(NA, firstYear[,2][-length(firstYear[,2])])
mean(abs((firstYear[,2]-forecastFirstYear)/firstYear[,2]),na.rm=T)*100

forecastSecondYear <- c(NA, secondYear[,2][-length(secondYear[,2])])
mean(abs((secondYear[,2]-forecastSecondYear)/secondYear[,2]),na.rm=T)*100

forecastThirdYear <- c(NA, thirdYear[,2][-length(thirdYear[,2])])
mean(abs((thirdYear[,2]-forecastThirdYear)/thirdYear[,2]),na.rm=T)*100

forecastFourthtYear <- c(NA, fourthYear[,2][-length(fourthYear[,2])])
mean(abs((fourthYear[,2]-forecastFourthtYear)/fourthYear[,2]),na.rm=T)*100

forecastFifthYear <- c(NA, fifthYear[,2][-length(fifthYear[,2])])
mean(abs((fifthYear[,2]-forecastFifthYear)/fifthYear[,2]),na.rm=T)*100

forecastAllYears <- c(NA, df[,2][-length(df[,2])])
mean(abs((df[,2]-forecastAllYears)/df[,2]),na.rm=T)*100

#Naive Graphs
#Split up into separate years to make it easier to see
plot(firstYear[,1],firstYear[,2], type="l", col="green",main="Actual vs naive Forecasted People First Year",xlab="Period",ylab="population")
lines(firstYear[,1],forecastFirstYear, type="l",col="blue")
legend("topleft",legend=c("Actual","Forecasted"),col=c("green","blue"),lty=1)

plot(secondYear[,1],secondYear[,2], type="l", col="green",main="Actual vs naive Forecasted People Second Year",xlab="Period",ylab="population")
lines(secondYear[,1],forecastSecondYear, type="l",col="blue")
legend("topleft",legend=c("Actual","Forecasted"),col=c("green","blue"),lty=1)

plot(thirdYear[,1],thirdYear[,2], type="l", col="green",main="Actual vs naive Forecasted People third Year",xlab="Period",ylab="population")
lines(thirdYear[,1],forecastThirdYear, type="l",col="blue")
legend("topleft",legend=c("Actual","Forecasted"),col=c("green","blue"),lty=1)

plot(fourthYear[,1],fourthYear[,2], type="l", col="green",main="Actual vs naive Forecasted People fourth Year",xlab="Period",ylab="population")
lines(fourthYear[,1],forecastFourthtYear, type="l",col="blue")
legend("topleft",legend=c("Actual","Forecasted"),col=c("green","blue"),lty=1)

plot(fifthYear[,1],fifthYear[,2], type="l", col="green",main="Actual vs naive Forecasted People fifth Year",xlab="Period",ylab="population")
lines(fifthYear[,1],forecastFifthYear, type="l",col="blue")
legend("topleft",legend=c("Actual","Forecasted"),col=c("green","blue"),lty=1)

plot(allYear[,1],allYear[,2], type="l", col="green",main="Actual vs naive Forecasted People",xlab="Period",ylab="population")
lines(allYear[,1],forecastAllYears, type="l",col="blue")
legend("topleft",legend=c("Actual","Forecasted"),col=c("green","blue"),lty=1)

#ARIMA
train_sample_size <- floor(0.75*nrow(df))

set.seed(123)
train_ind = sample(seq_len(nrow(df)),size = train_sample_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]
mape <- function(actual,pred){
  mape <- mean(abs((actual - pred)/actual))*100
  return (mape)
}
train.ts <- xts(train[,-1],order.by = train[,1])
naive_mod <- naive(train.ts, h = 7)
summary(naive_mod)
test$naive = 18
mape(test$Number.of.patients,test$naive)

se_model <- ses(train.ts, h = 366)
summary(se_model)

df_fc = as.data.frame(se_model)
test$simplexp = df_fc$`Point Forecast`
mape(test$Number.of.patients, test$simplexp) 

holt_model <- holt(train.ts, h = 366)
summary(holt_model)
df_holt = as.data.frame(holt_model)
test$holt = df_holt$`Point Forecast`
mape(test$Number.of.patients, test$holt)

arima_model <- auto.arima(train.ts)
summary(arima_model)
fore_arima = forecast::forecast(arima_model, h=366)
df_arima = as.data.frame(fore_arima)
test$arima = df_arima$`Point Forecast`
mape(test$Number.of.patients,test$arima)  

#Moving Average
#ARIMA
install.packages("urca")
install.packages("aTSA")
library(urca)
library(aTSA)
train.ts %>% ur.kpss() %>% summary()

train$Date <- as.numeric(train$Date)

train$Number.of.patients = tsclean(train.ts)

cbind("Original" = train.ts,
      "Logs" = log(train.ts),
      "Seasonally\n differenced logs" =
        diff(log(train.ts),12),
      "Doubly\n differenced logs" =
        diff(diff(log(train.ts),12),1)) %>%
  autoplot(facets=TRUE) +
  xlab("Year") + ylab("") +
  ggtitle("Number of patients over Years")
stationary <- diff(log(train.ts),12)
stationary %>% ur.kpss() %>% summary()

Acf(stationary,main="")
Pacf(stationary,main="")

count_d1 = diff(stationary, difference =1)
plot(count_d1)

Acf(count_d1,main="")
Pacf(count_d1,main="")

auto.arima(stationary, seasonal = FALSE)
fit<-auto.arima(stationary, seasonal = FALSE)
tsdisplay(residuals(fit),lag.max=45,main = "(1,1,1) model residuals")

fit2=arima(stationary, order=c(1,1,7)) 
tsdisplay(residuals(fit),lag.max=20,main = "seasonal model residuals")

fcast <- forecast:: forecast(fit2,h=7)
plot(fcast)
predictions = list(preds = -0.1930222,  -0.3669292 , -0.4253246,  -0.5727103,   -0.2758548,  -0.5658944 ,  -0.1077718)



fit_with_seasonalty = auto.arima(stationary, seasonal=TRUE)
seas_fcast <- forecast:: forecast(fit_with_seasonalty, h=7)
plot(seas_fcast)

