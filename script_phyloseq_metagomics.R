#!/bin/bash/R
#Desenvolvido por: Fernanda Luiza Ferrari
#Com base em: 

## Análise metagênomica de dados de sequeciamento shotgun mapeados e 
## classificados com metaphlan e humann

## Bibliotecas necessárias
## tydeverse - pacote para manipoulação de dados, plotagem de dados, etc.
## jcolors - pacote para customização de cores nos gráficos
## plyloseq - pacote do Bioconductor para importar, armazenar, analisar e exibir
## dados complexos de sequenciamento filogenético que já foram agrupados em OTUs
## ape - biblioteca para Análises de Filogenética e Evolução
## dendextend - biblioteca para criação de dendogramas
## curatedMetagenomicData - pacote do bioconductor para lidar com dados 
## metagenomicos oriundos do metaphlan e do humann 
## RColorBrewer - pacote de palettas de cores



library(tidyverse)
library(jcolors)
library(phyloseq)
library(dendextend)
library(ape)
library(curatedMetagenomicData)
library(RColorBrewer)
library(microbiome)

## Setando e verificando o diretório

setwd("/home/fernanda/Documents/")
getwd()

## Função para transformar dados do metaphlan em dados para phyloseq
## Retirado de: https://gist.github.com/lwaldron/512d1925a8102e921f05c5b25de7ec94

metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

################################################################################

## Importando os dados
## dados_abudancia_metaphlan - arquivo csv com amostras em colunas e táxons em
## linhas.
## dados_metadados - arquivos com amostras em linhas (da mesma forma que estão 
## identificados no dados_abundancia_metaplhan) e informações relevantes em 
## colunas.

dados_abundancia_metaphlan = read.csv("teste.csv", sep=",", strip.white = T, 
                                      stringsAsFactors = F, row.names = 1)

dados_metadados = read.csv("metadados_teste.csv", sep ="\t", row.names = 1)

head(dados_abundancia_metaphlan)
head(dados_metadados)

## Somando colunas para verificar a abundancia, todas devem ter soma de 700 ou 
## próximo disso, pois trata-se de abundância relativa cada taxon (p, g, o, etc)
## deve somar 100, como trata-se de 7 ranks o total deve ficar próximo desse 
## valor.

colSums(dados_abundancia_metaphlan[,])

## Covertando tabela para um objeto phyloseq

abundancia_metaphlan_phyloseq = metaphlanToPhyloseq(dados_abundancia_metaphlan)
                                                    
## Criando componentes e objeto phyloseq
## Unidade Taxanômica Operacional (OTU)
## Dados de taxonomia (TAX)
## Metadso e informações sobre as amostras (SAM)

OTU = otu_table(abundancia_metaphlan_phyloseq)
TAX = tax_table(abundancia_metaphlan_phyloseq)
SAM = sample_data(dados_metadados, errorIfNULL = TRUE) ###estrutura taxonômica

head(OTU)
head(TAX)
head(SAM)

## Adicionando dados de amostras ao objeto phyloseq

abundancia_phyloseq = metaphlanToPhyloseq(abundancia_metaphlan_phyloseq, metadat = SAM)
abundancia_phyloseq

## Criando a árvore taxonomica

random_tree = rtree(ntaxa(abundancia_phyloseq), rooted=TRUE, 
                    tip.label=taxa_names(abundancia_phyloseq))

abundancia_phyloseq = merge_phyloseq(abundancia_phyloseq, random_tree)

################################################################################

## Visualização gráfica de abundância de filos e gêneros
## Transformando dados
## Para gênero trocar o fill = "Phylum" para fill = "Genus"

ps_rel_abund = phyloseq::transform_sample_counts(abundancia_phyloseq, function(x){x / sum(x)})

## Mantendo apenas os mais abundantes em todas as amostras
top_taxon = genefilter_sample(ps_rel_abund, filterfun_sample(topp(0.30)), A=5)
summary(keepotu)
ps_rel_abund_top = subset_taxa(ps_rel_abund, top_taxon)

phyloseq::plot_bar(ps_rel_abund_top, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", 
           position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  # facet_wrap(~ sexo, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(legend.text = element_text(size = 10), legend.title=element_text(size=12), 
        plot.title = element_text(size = 12), plot.subtitle = 
          element_text(size = 10,  color = "darkslategrey"), 
        axis.text.x = element_text(angle=90,
                                   hjust=1,  
                                   size=10), axis.text.y =
          element_text(angle=0, hjust=1, size=10))+
  scale_fill_hue(c=45, l=50)

## Salvando gráfico

ggsave("phylum-abundance.svg", width = 10, height = 10, dpi = 300)


## Visualização gráfica de prevalência de filos e gêneros
## Transformando dados
## Para gênero trocar o fill = "Phylum" para fill = "Genus"

ps_taxon <- phyloseq::tax_glom(abundancia_phyloseq, "Phylum")
phyloseq::taxa_names(ps_taxon) <- phyloseq::tax_table(ps_taxon)[, "Phylum"]


phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = sample, y = Abundance, fill = OTU)) +
  geom_boxplot(outlier.shape  = TRUE) +
  geom_col() +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  theme(legend.text = element_text(size = 10), 
        legend.title=element_text(size=12), 
        plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 10, 
                                     color = "darkslategrey"), 
        axis.text.x = element_text(angle=90,
                                   hjust=1, 
                                   size=8),         
        axis.text.y = element_text(angle=0, 
                                   hjust=1, 
                                   size=8)) +
  guides(color=FALSE, fill=FALSE )


## Salvando gráfico

ggsave("phylum-prevalence.svg", width = 10, height = 10, dpi = 300)


################################################################################

## Criando dendograma de amostras

ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)

## Calculando a distância

bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")

## Dendograma
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))

## Anotando o dendograma com base no metadados

meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c("F" = "red", "M" = "blue")
labels_colors(ward) <- colorCode[meta$sexo][order.dendrogram(ward)]


plot(ward)

################################################################################

## PCoA - análise de coordenadas principais baseada em distâncias
## Bray é uma forma de calcular as distâncias

ord_bray <- ordinate(abundancia_phyloseq, method="PCoA", distance=dist_bray)
samples_infos <- data.frame(sample_data(abundancia_phyloseq))

samples_infos$bray_cluster_2 <- factor(pam(dist_bray, k=2, cluster.only = T))
sample_data(abundancia_phyloseq) <- samples_infos


## Plotando visualização com dois clusters

pcoa_2_clusters = plot_ordination(abundancia_phyloseq, ord_bray, 
                                  color="bray_cluster_2") + 
  ggplot2::ggtitle("Bray-Curtis PCoA with 2 Clusters") 

pcoa_2_clusters

## Plotando visualização com base na abundância de Bacteroides e Prevotella
## Extraindo os vetores

pc1 <- ord_bray$vectors[,1]
pc2 <- ord_bray$vectors[,2]

## Calculando abundancia de Bacteroides e Prevotela

otu_tax <- attr(otu_table(abundancia_phyloseq), "dimnames")[[1]]
abundancia_phyloseq_bact <- otu_table(abundancia_phyloseq)[grep("s__Bacteroides", 
                                                       otu_tax),]
sum_bacteroides <- apply(abundancia_phyloseq_bact, 2, sum)


abundancia_phyloseq_prev = as.numeric(otu_table(abundancia_phyloseq))[grep("s__Prevotella", 
                                                                           otu_tax),]
sum_prevotella <- apply(abundancia_phyloseq_prev, 2, sum)


sample_data(df2) <- samp
plot_ordination(abundancia_phyloseq, ord_bray, color="Prevotella")


## Ordenando dataframe

abundancia_phyloseq_ordinate <- data.frame(pc1, pc2, bact = sum_bacteroides, 
                                           prev = sum_prevotella, 
                                           bray2 = as.numeric(samples_infos$bray_cluster_2)+20)

df_bact <- abundancia_phyloseq_ordinate[abundancia_phyloseq_ordinate$bray2==21,]
df_prev <- abundancia_phyloseq_ordinate[abundancia_phyloseq_ordinate$bray2==22,]


## Criando os plot separados e plotando junto

p <- ggplot() 

## Prevotella
p <- p +  geom_point(data= df_prev, aes(x=pc1, y=pc2, shape=factor(bray2), 
                                        fill=sum_prevotella), shape=21, size=4) +
  scale_fill_gradient(low="white", high="red",
                                       guide = guide_colorbar(title="Prevotella \n(cluster 2)"))


## Bacteroides
p <- p + geom_point(data=df_bact, aes(x=pc1, y=pc2, shape=factor(bray2), color=bact), 
                    shape=22, size=4) +
  scale_color_gradient(low="gray", high="green", 
                              guide=guide_colorbar(title="Bacteroides \n(cluster 1)"))

## Labels
p <- p + labs(x="Axis 1",y="Axis 2", title="PCoA by selected abundances")

p


## PCoA com base em grupos
## Para mudar o grupo alterar color = "sexo"

p <- ggplot(data=sample, aes(x=pc1, y=pc2, shape=sexo,
                             color=sexo)) 
p <- p + geom_point()
p <- p + labs(x="Axis 1", y="Axis 2", title="PCoA by group")

p

################################################################################

## Visualizações com heatmaps

## Mantém apenas taxa que estão entre os 5% mais abundantes em pelo menos 
## cinco amostras: (A=5)

keepotu = genefilter_sample(abundancia_phyloseq, filterfun_sample(topp(0.05)), A=5)
subset_taxa(abundancia_phyloseq, keepotu)
loman.filt = subset_taxa(abundancia_phyloseq, keepotu)



## Esta função suporta um grande número de métodos de distância e ordenação

heatmap = plot_heatmap(loman.filt, method="PCoA", distance="bray", sample.label = NULL, 
                       taxa.label = NULL, low = "#000033", high = "#66CCFF", 
                       na.value = "black", max.label = 250,  
                       title = "Distance Bray-Curtis") +
  theme(legend.text = element_text(size = 12), 
        legend.title=element_text(size=15), 
        plot.title = element_text(size = 15), 
        plot.subtitle = element_text(size = 10, 
                                     color = "darkslategrey"), 
        axis.text.x = element_text(angle=90, 
                                   hjust=1, 
                                   size=10),
        axis.text.y = element_text(size=11))



heatmap



################################################################################

## Histograma de espécies e gêneros

loman.sp = subset_taxa(df2, !is.na(Species))


par(mar = c(20, 4, 10, 0) + 0.25) 

histo_species= barplot(sort(taxa_sums(loman.sp), TRUE)[1:20] / nsamples(loman.sp),
            ylab = "Abundance", 
            las = 2,
            col=rgb(0.1,0.1,0.7,0.5) ) +
  title(main = "20 most abundant species") +
  theme(plot.title = element_text(size =8))

histo_species


loman.gn = subset_taxa(df2, !is.na(Genus))


par(mar = c(20, 4, 10, 0) + 0.25) 

histo_genus= barplot(sort(taxa_sums(loman.gn), TRUE)[1:20] / nsamples(loman.gn),
                    ylab = "Abundance", 
                    las = 2,
                    col=rgb(0.1,0.1,0.7,0.5) ) +
  title(main = "20 most abundant genus") +
  theme(plot.title = element_text(size =8))

histo_genus



################################################################################

## Diversidade Alpha e Beta

## A diversidade alfa (α-diversidade) é definida como a diversidade média de 
## espécies em diferentes locais ou habitats dentro de uma escala local.

## Shannon - quantifica a incerteza na previsão da identidade da espécie de 
## um indivíduo que é retirado aleatoriamente do conjunto de dados.
## Simpson -  mede a probabilidade de que dois indivíduos selecionados
## aleatoriamente de uma amostra pertençam à mesma espécie
## InvSimpson - o índice representa a probabilidade de que dois indivíduos 
## selecionados aleatoriamente de uma amostra pertençam a espécies diferentes.

## https://www.cd-genomics.com/microbioseq/the-use-and-types-of-alpha-diversity-metrics-in-microbial-ngs.html
## http://www.countrysideinfo.co.uk/simpsons.htm


alphas = c("Shannon", "Simpson", "InvSimpson")

## Com base em espécies
riqueza_alpha_loman.sp = estimate_richness(loman.sp, measures = alphas)

## Com base em gêneros

riqueza_alpha_loman.gn = estimate_richness(loman.gn, measures = alphas)


## Plotando


riqueza_alpha_loman.sp$samples = row.names(riqueza_alpha_loman.sp)
riqueza_alpha_loman.gn$samples = row.names(riqueza_alpha_loman.gn)


riqueza_alpha_loman.sp  %>% gather(key = metric, value = value, 
                                   c("Shannon", "Simpson", "InvSimpson")) %>%
  mutate(metric = factor(metric, 
                         levels = c("Shannon", "Simpson", "InvSimpson"))) %>%
  ggplot(aes(x = samples, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = samples), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") + 
  theme(legend.text = element_text(size = 7), legend.title=element_text(size=6), 
        plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 10, color = "darkslategrey"), 
        axis.text.x = element_text(angle=90, hjust=1, size=8),
        axis.text.y = element_text(angle=0, hjust=1,  size=8))




riqueza_alpha_loman.gn  %>% gather(key = metric, value = value, 
                                   c("Shannon", "Simpson", "InvSimpson")) %>%
  mutate(metric = factor(metric, 
                         levels = c("Shannon", "Simpson", "InvSimpson"))) %>%
  ggplot(aes(x = samples, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = samples), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") + 
  theme(legend.text = element_text(size = 7), legend.title=element_text(size=6), 
        plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 10, color = "darkslategrey"), 
        axis.text.x = element_text(angle=90, hjust=1, size=8),
        axis.text.y = element_text(angle=0, hjust=1,  size=8))


sample_data(abundancia_phyloseq)$shannon_sp = riqueza_alpha_loman.sp$Shannon
sample_data(abundancia_phyloseq)$shannon_gn = riqueza_alpha_loman.gn$Shannon

################################################################################

## Outros índices

## Riqueza (número de espécies)
richness <- richness(loman.sp)

## Domínio - refere-se à abundância das espécies mais abundantes
dominance <- dominance(loman.sp, index = "all")


## Raridade - quantifica a concentração de taxa raros ou de baixa abundância
rarity <- rarity(loman.sp, index = "all")


## Cobertura
coverage <- coverage(abundancia_phyloseq, threshold = 0.5)


## Core  

core_abundance <- core_abundance(loman.sp, detection = .1/100, prevalence = 50/100)


## Medida de diversidade da comunidade.

inequality <- inequality(abundancia_phyloseq)

## Uniformidade

evenness <- evenness(abundancia_phyloseq, "all")


################################################################################

## Testando diferenças na diversidade alfa com base em grupos
## É recomendado o teste não paramétrico de Kolmogorov-Smirnov para comparações
## de dois grupos quando não há covariáveis relevantes.


d <- meta(abundancia_phyloseq)
d$diversity <- sample_data(abundancia_phyloseq)$diversity
spl <- split(d$diversity, d$Grupos)

## Kolmogorov-Smironv test
pv <- ks.test(spl$F, spl$M,)$p.value
pv

# Adjust the p-value
padj <- p.adjust(pv)

padj

################################################################################

## Arvore taxonomica

p = plot_tree(abundancia_phyloseq, color="Samples", shape="seo", 
              label.tips="Genus", size="Abundance")

p




