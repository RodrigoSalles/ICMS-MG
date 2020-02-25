######################################################################
######################################################################
##               Séries Temporais - Trabalho final                  ##
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

### Análise exploratória dos dados ###
summary(ipeadata)
sd(ipeadata$valor)


### Estacionaridade ###
# ACF
acf(ipeadata$valor, main = 'Função de Autocorrelação')

### Decomposição da série ###
ipeadata_dec =decompose(ipeadata.ts)
plot(ipeadata_dec)

### Augmented Dickey Fuller Test : raiz unitária ###
adf.test(ipeadata$valor) #p-value = 0.5561 - Série não estacionária

### Teste KPSS ###
kpss.test(ipeadata$valor) # p-value = 0.01 - Série não estacionária

### Nº de diferenciações ###
ndiffs(ipeadata$valor) # Uma diferenciação necessária

### Nº de diferenciações sazonais ###
nsdiffs(ipeadata.ts) # Não precisa de diferenciação sazonal



### Tratamento para alcançar estacionaridade ###
# Aplicando uma diferenciação
ipea_dif_1 = diff(ipeadata$valor, differences = 1)

# Verificar ACF
acf(ipea_dif_1, main = 'Função de Autocorrelação')
pacf(ipea_dif_1, main = 'Função de Autocorrelação Parcial')

ndiffs(ipea_dif_1) # Não precisa de mais diferenciações

# Augmented Dickey Fuller Test : raiz unitária , 1 dif
adf.test(ipea_dif_1) # p-value = 0.01, com uma diferenciação a série é estacionária.

# Teste KPSS
kpss.test(ipea_dif_1) # p-value = 0.1, com uma diferenciação a série é estacionária.


#### Modelos ####
# Modelo inicial
mod = auto.arima(ipeadata$valor)
mod # ARIMA(1,1,1) with drift

# Modelo Arima(1,1,1) com drift
mod = Arima(ipeadata.ts, order = c(1,1,1), include.drift = T)

# Verificar parâmetros e resíduos do modelo arima(1,1,1)
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

# Verificar parâmetros e resíduos do modelo arima(3,1,3)
sarima(ipeadata$valor, 3,1,3 )


# Testar modelo com sazonaridade
mod_saz = sarima(ipeadata$valor, 1,1,1, 2,0,1,12)
mod_saz = Arima(ipeadata$valor,order=c(1,1,1), seasonal=c(2,0,1), include.drift = T)
mod_saz$ttable

# Testar modelo com transformação Box Cox
# Valor ideal de lambda
BoxCox.lambda(ipeadata.ts) # 0.5993851
# Verificar efeito da transformação
plot.ts(BoxCox(ipeadata.ts, lambda = 0.5993851))
# Análise de resíduos
acf(box_trans$residuals)



#### Previsões ####
# Previsão dos próximos 12 meses: mod
forecast(mod, h = 12, level = c(80,95))
mod_prev = forecast(mod, h = 12, level = c(80,95))
# Plotar gráfico de previsão mod
plot(forecast(mod, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')

# Previsão dos próximos 12 meses: mod1
forecast(mod1, h = 12, level = c(80,95))
# Plotar gráfico de previsão mod1
plot(forecast(mod1, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')

# Previsão dos próximos 12 meses: mod_saz
forecast(mod_saz, h = 12, level = c(80,95))
# Plotar gráfico de previsão mod1
plot(forecast(mod_saz, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')

# Previsão dos próximos 12 meses: box_trans
box_trans = forecast(auto.arima(ipeadata.ts, lambda=0.5993851, biasadj=TRUE),h=12)
forecast(box_trans, h = 12, level = c(80,95))
# Plotar gráfico de previsão mod1
plot(forecast(box_trans, h = 12, level = c(80,95)),xlab = 'Anos', ylab = 'Valores')



###### Redes Neurais Artificiais ######
fit <- nnetar(ipeadata.ts, lambda=0)
autoplot(forecast(fit,h=12), xlab = 'Anos', ylab = 'Valores', main = 'Previsão - RNA - 2016')
fit

fcast <- forecast(fit, PI=TRUE, h=12)
autoplot(fcast, xlab = 'Anos', ylab = 'Valores', main = 'Previsão - RNA - 2016')

# Verificar resíduos das RNAs
acf(fcast$residuals[which(!is.na(fcast$residuals))], main = 'Função de Autocorrelação RNA')



#### Combinação de previsores ####

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


