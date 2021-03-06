######################################################################
######################################################################
##               S�ries Temporais - Trabalho final                  ##
##                                                                  ##
## Prof. Maria Eduarda da Rocha                                     ##
## Aluno: Rodrigo Salles                                            ##
######################################################################
######################################################################

### Bibliotecas utilizadas
require(fpp2)
require(xts)
library(ggplot2)
library(tseries)
library(forecast)
library(astsa)

### Carregar os dados ###
ipeadata <- read.csv("C:/Users/RodrigoSalles/Desktop/ICMS_TS/data_3")
View(ipeadata)

### Eliminar primeira coluna ###
ipeadata$X <- NULL
View(ipeadata) # Observar os dados prontos

### Grafico inicial ###
ipeadata.ts = ts(ipeadata$valor, frequency = 12, start = c(1993,1))
plot(ipeadata.ts, col=1, lwd=2, main = 'ICMS - Minas Gerais(01/1993 - 12/2015 ', xlab = 'Anos', ylab = 'Valores')
grid (NULL,NULL, lty = 6,  col = "cornsilk2") 

### An�lise explorat�ria dos dados ###
summary(ipeadata)
sd(ipeadata$valor)


### Estacionaridade ###
# ACF
acf(ipeadata$valor, main = 'Fun��o de Autocorrela��o')

### Decomposi��o da s�rie ###
ipeadata_dec =decompose(ipeadata.ts)
plot(ipeadata_dec)

### Augmented Dickey Fuller Test : raiz unit�ria ###
adf.test(ipeadata$valor) #p-value = 0.5561 - S�rie n�o estacion�ria

### Teste KPSS ###
kpss.test(ipeadata$valor) # p-value = 0.01 - S�rie n�o estacion�ria

### N� de diferencia��es ###
ndiffs(ipeadata$valor) # Uma diferencia��o necess�ria

### N� de diferencia��es sazonais ###
nsdiffs(ipeadata.ts) # N�o precisa de diferencia��o sazonal



### Tratamento para alcan�ar estacionaridade ###
# Aplicando uma diferencia��o
ipea_dif_1 = diff(ipeadata$valor, differences = 1)

# Verificar ACF
acf(ipea_dif_1, main = 'Fun��o de Autocorrela��o')
pacf(ipea_dif_1, main = 'Fun��o de Autocorrela��o Parcial')

ndiffs(ipea_dif_1) # N�o precisa de mais diferencia��es

# Augmented Dickey Fuller Test : raiz unit�ria , 1 dif
adf.test(ipea_dif_1) # p-value = 0.01, com uma diferencia��o a s�rie � estacion�ria.

# Teste KPSS
kpss.test(ipea_dif_1) # p-value = 0.1, com uma diferencia��o a s�rie � estacion�ria.


#### Modelos ####
# Modelo inicial
mod = auto.arima(ipeadata$valor)
mod # ARIMA(1,1,1) with drift

# Modelo Arima(1,1,1) com drift
mod = Arima(ipeadata.ts, order = c(1,1,1), include.drift = T)

# Verificar par�metros e res�duos do modelo arima(1,1,1)
sarima(ipeadata$valor, 1,1,1)

### Testar modelos ###
sarima(ipeadata$valor, 2,1,1 )
sarima(ipeadata$valor, 1,1,2 )
sarima(ipeadata$valor, 3,1,1 )
sarima(ipeadata$valor, 3,1,2)
sarima(ipeadata$valor, 3,1,3 )
sarima(ipeadata$valor, 1,1,2 )
sarima(ipeadata$valor, 1,1,3 )


# Modelo com melhor resultado
mod1 <- Arima(ipeadata.ts,order = c(3,1,3),include.drift=TRUE)

# Verificar par�metros e res�duos do modelo arima(3,1,3)
sarima(ipeadata$valor, 3,1,3 )


# Testar modelo com sazonaridade
mod_saz = sarima(ipeadata$valor, 1,1,1, 2,0,1,12)
mod_saz = Arima(ipeadata$valor,order=c(1,1,1), seasonal=c(2,0,1), include.drift = T)
mod_saz$ttable

# Testar modelo com transforma��o Box Cox
# Valor ideal de lambda
BoxCox.lambda(ipeadata.ts) # 0.5993851
# Verificar efeito da transforma��o
plot.ts(BoxCox(ipeadata.ts, lambda = 0.5993851))
# An�lise de res�duos
acf(box_trans$residuals)



#### Previs�es ####
# Previs�o dos pr�ximos 12 meses: mod
forecast(mod, h = 12, level = c(80,95))
mod_prev = forecast(mod, h = 12, level = c(80,95))
# Plotar gr�fico de previs�o mod
plot(forecast(mod, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')

# Previs�o dos pr�ximos 12 meses: mod1
forecast(mod1, h = 12, level = c(80,95))
# Plotar gr�fico de previs�o mod1
plot(forecast(mod1, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')

# Previs�o dos pr�ximos 12 meses: mod_saz
forecast(mod_saz, h = 12, level = c(80,95))
# Plotar gr�fico de previs�o mod1
plot(forecast(mod_saz, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')

# Previs�o dos pr�ximos 12 meses: box_trans
box_trans = forecast(auto.arima(ipeadata.ts, lambda=0.5993851, biasadj=TRUE),h=12)
forecast(box_trans, h = 12, level = c(80,95))
# Plotar gr�fico de previs�o mod1
plot(forecast(box_trans, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')



###### Redes Neurais Artificiais ######
fit <- nnetar(ipeadata.ts, lambda=0)
autoplot(forecast(fit,h=12), xlab = 'Anos', ylab = 'Valores', main = 'Previs�o - RNA - 2016')
fit

fcast <- forecast(fit, PI=TRUE, h=12)
autoplot(fcast, xlab = 'Anos', ylab = 'Valores', main = 'Previs�o - RNA - 2016')

# Verificar res�duos das RNAs
acf(fcast$residuals[which(!is.na(fcast$residuals))], main = 'Fun��o de Autocorrela��o RNA')



#### Combina��o de previsores ####

train <- ipeadata.ts

ETS <- forecast(ets(train), h=12)
ARIMA <- forecast(auto.arima(train, lambda=0, biasadj=TRUE),h=12)
STL <- stlf(train, lambda=0, h=12, biasadj=TRUE)
NNAR <- forecast(nnetar(train), h=12)
TBATS <- forecast(tbats(train, biasadj=TRUE), h=12)
Combination <- (ETS[["mean"]] + ARIMA[["mean"]] + NNAR[["mean"]] + TBATS[["mean"]])/4
Combination


# Verificar residuos dos previsores
par(mfrow=c(2,2))

acf(ETS$residuals)

acf(ARIMA$residuals)

acf(STL$residuals)

acf(TBATS$residuals)

acf(na.remove(fit$residuals))


