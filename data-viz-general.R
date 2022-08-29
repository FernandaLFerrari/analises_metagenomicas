#!/bin/bash/R
#Desenvolvido por: Fernanda Luiza Ferrari
#Com base em: https://r-graph-gallery.com/


## Data Viz

## Bibliotecas

library(jcolors)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)


myPalette <- brewer.pal(10, "Pastel2") 

setwd("/home/fernanda/Desktop/")
getwd()

## Importando os dados
dados_metadados = read.csv("metadados_teste.csv", sep ="\t", row.names = 1)

## Definindo padrões de exibição

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(dados_metadados$IDADE, horizontal=TRUE, ylim=c(18,90),xaxt="n", 
        col=rgb(0.1,0.1,0.7,0.5) , frame=F)

par(mar=c(4, 3.1, 1.1, 2.1))
hist(dados_metadados$IDADE , breaks=30 , col=rgb(0.1,0.6,0.6,0.9) , border=F , 
     main="" , xlab="Idade Geral", xlim=c(18,90))  



## Filtrando apenas os que possuem IMC

dados_metadados$imc_verifica = ifelse(dados_metadados$imc>0,dados_metadados$imc,"n")
dados_metadados_fill=filter(dados_metadados,dados_metadados$imc_verifica != "n" )

# Plotando gráfico de Idade x Sexo

dados_metadados %>%
  ggplot( aes(x=sexo, y=IDADE, fill=sexo)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_ipsum_rc(grid="Y") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Idade por sexo") +
  xlab("") + 
  scale_color_ipsum() +
  scale_fill_ipsum() 


# Plotando gráfico de IMC x Sexo

dados_metadados_fill %>%
  ggplot( aes(x=sexo, y=imc, fill=sexo)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_ipsum_rc(grid="X")  +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("IMC por sexo") +
  xlab("") + 
  scale_color_ipsum() +
  scale_fill_ipsum() 


## Plotando informações sobre UF's


layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

df <- data.frame(
  name=c("RJ","SC","SP","NULL","PE", "PR", "RS", "SE") ,  
  value=c(8,9,7,2,1,2,3,1))

par(mar=c(0, 3.1, 1.1, 2.1))


ggplot(df, aes(x=factor(name, level=c("SC", "RJ", "SP", "RS", "PR", "PE", "SE", "NULL")), 
                        y= value, fill=name )) + 
  geom_col() +
  theme_ipsum_rc(grid="X")  +
  scale_color_ipsum() +
  scale_fill_ipsum()  +
  theme(legend.position="none", plot.title = element_text(size=10))+
  ggtitle("Amostras por UF") +
  xlab("UF") +
  ylab("Qtde. Amostras") 


par(mar=c(4, 3.1, 1.1, 2.1))
pie(c(8,9,7,2,1,2,3,1) , labels = c("RJ","SC","SP","NULL","PE", "PR", "RS", "SE"), 
    border="white", col=myPalette )


