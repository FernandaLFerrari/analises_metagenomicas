#!/bin/bash/R
#Desenvolvido por: Fernanda Luiza Ferrari
#Com base em: https://r-graph-gallery.com/

library(tidyverse)
library(splitstackshape)

## Importando dados e duplicando coluna

dados_abundancia_metaphlan = read.csv("teste-ok.csv", sep=",", strip.white = T, 
                                      stringsAsFactors = F)

dados_abundancia_metaphlan$duplicate_tax = dados_abundancia_metaphlan$categoria


## Separando coluna de taxonomia

dados_abundancia_metaphlan = cSplit(dados_abundancia_metaphlan, "duplicate_tax", 
                                    sep="__")

## Somando abundância

dados_abundancia_metaphlan$soma = rowSums(dados_abundancia_metaphlan[,3:32])

## Ordenando e pegando os 10 de cada taxon

## Phylum

df_phylum = dados_abundancia_metaphlan %>%
                        filter(duplicate_tax_1=="p")  


df_phylum = df_phylum[order(df_phylum$soma, decreasing = TRUE),]
df_phylum$abd_relativa = df_phylum$soma/sum(df_phylum$soma)*100
df_phylum= df_phylum[1:5,]


ggplot(df_phylum, aes(x=reorder(duplicate_tax_2, abd_relativa), y=abd_relativa,
                      fill= duplicate_tax_2)) + 
  geom_col() +  
  scale_fill_ipsum() +
  coord_flip() + 
  theme(plot.title = element_text(size=10)) +
  labs(x="Phylum", y="Abundância Relativa",
       title="Abundância de principais filos") + 
  theme_ipsum_rc(grid="X") +
  guides(color=FALSE, fill=FALSE ) 



## Genus

df_genus = dados_abundancia_metaphlan %>%
  filter(duplicate_tax_1=="g")  


df_genus = df_genus[order(df_genus$soma, decreasing = TRUE),]
df_genus$abd_relativa = df_genus$soma/sum(df_genus$soma)*100
df_genus= df_genus[1:9,]


ggplot(df_genus, aes(x=reorder(duplicate_tax_2, abd_relativa), y=abd_relativa,
                      fill= duplicate_tax_2)) + 
  geom_col() +  
  scale_fill_ipsum() +
  coord_flip() + 
  theme(plot.title = element_text(size=10)) +
  labs(x="Genus", y="Abundância Relativa",
       title="Abundância de principais gêneros") + 
  theme_ipsum_rc(grid="X") +
  guides(color=FALSE, fill=FALSE ) 



## Espécie

df_species = dados_abundancia_metaphlan %>%
  filter(duplicate_tax_1=="s")  


df_species = df_species[order(df_species$soma, decreasing = TRUE),]
df_species$abd_relativa = df_species$soma/sum(df_species$soma)*100
df_species= df_species[1:9,]


ggplot(df_species, aes(x=reorder(duplicate_tax_2, abd_relativa), y=abd_relativa,
                     fill= duplicate_tax_2)) + 
  geom_col() +  
  scale_fill_ipsum() +
  coord_flip() + 
  theme(plot.title = element_text(size=10)) +
  labs(x="Species", y="Abundância Relativa",
       title="Abundância de principais espécies") + 
  theme_ipsum_rc(grid="X") +
  guides(color=FALSE, fill=FALSE ) 



