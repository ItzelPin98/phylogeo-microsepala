## En este script vamos a correr un Analisis de Componentes Princiaples (PCA)
## y un Analisis de Variables Canonicas (CVA)


##cargar librerias
library(geomorph)
library(Morpho)
library(ggplot2)
library(dplyr)


# Paso 1: Cargar datos alineados y aplicarles la funcion gm.prcomp para reducir 
# la dimensionalidad con un PCA

ProcMicro 
colnames(ProcMicro)
pcaMicro <- gm.prcomp(ProcMicro$coords)

## veremos le resumen para determinar cuantos componentes vamos a retener
summary(pcaMicro)

#Importance of Components:
                           #Comp1       Comp2        Comp3        Comp4        Comp5
#Eigenvalues            0.002530256 0.001074614 0.0004041302 0.0003305103 0.0002828739
#Proportion of Variance 0.458906983 0.194900440 0.0732961931 0.0599439291 0.0513042097
#Cumulative Proportion  0.458906983 0.653807423 0.7271036162 0.7870475453 0.8383517549

names(pcaMicro)
names(ProcMicro)
rownames(pcaMicro$x)

## vemos los scores

pc_scores <- as.data.frame(pcaMicro$x[, 1:3])
pc_scores$extractn <- rownames(pc_scores)  # añadir columna de ID
colnames(pc_scores)[1:3] <- c("PC1", "PC2", "PC3")
rownames(pc_scores) <- NULL
pc_scores$n_ord <- rownames(pc_scores)


capture.output(summary(pcaMicro), file="./02_out/PCAmorfo_RADcolls_output_20250809.txt")
write.csv(pc_scores, "./02_out/PCAmorfo_RADcolls_loadings_20250809.csv")



##ordenamos las columnas
pc_scores <- pc_scores[, c("n_ord" , "extractn" , "PC1", "PC2", "PC3")]
colnames(pc_scores)

##Paso 2: Importar el archivo de metadata para remplazar el coll_name por los nombres de la tip


list.files("01_data")
grupo <- read.csv("./01_data/metadata_RADseq_micro_sin_outgroups_20250809.csv")
dim(grupo) #178 22
colnames(grupo)

##juntamos los meta

metadata <- merge(grupo, pc_scores, by = "extractn")
dim(metadata) #178 26
colnames(pca_scores_tips)
pca_scores_tips <- metadata[c("n_ord.x", "extractn", "tip_iapr_1",
                                     "cluster_adm","color_admix", "PC1",
                                     "PC2", "PC3")]

#guardamos los pca scores para hacer las distancias

write.csv(pca_scores_tips, "./02_out/DAT_morfoPCA_RADcolls_toDistMor_20250809.csv", row.names=F)

# Graficar

pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 1] <- "grupo_01"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 2] <- "grupo_02"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 3] <- "grupo_03"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 4] <- "grupo_04"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 5] <- "grupo_05"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 6] <- "grupo_06"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 7] <- "grupo_07"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 8] <- "grupo_08"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 9] <- "grupo_09"
pca_scores_tips$cluster_adm[pca_scores_tips$cluster_adm == 10] <- "grupo_10"

summary(pcaMicro)

#PC1 46%
#PC2 19%
#PC3 7%

pdf("./02_out/PLOT_PCAmorfo_20250809.pdf", width=14, height = 10)
ggplot(pca_scores_tips, aes(x = PC1, y = PC2, fill = cluster_adm)) + 
  scale_fill_manual(values = c("grupo_01" = "#9e0142", 
                               "grupo_02" = "#d53e4f", 
                               "grupo_03" = "#f46d43",
                               "grupo_04" = "#fdae61", 
                               "grupo_05" = "#fee088", 
                               "grupo_06" = "#e6f598",
                               "grupo_07" = "#abdda4", 
                               "grupo_08" = "#66c2a5", 
                               "grupo_09" = "#3288bd",
                               "grupo_10" = "#5E4FA2")) +
  xlab("PC 1: 46%") +
  ylab("PC 2: 19%") +
  geom_hline(yintercept = 0, lty = 2) +  
  geom_vline(xintercept = 0, lty = 2) +  
  geom_point(shape = 21,alpha = 3, size = 4, ) 


#dev.off()

pdf("./02_out/PLOT_PCAmorfoPC2vsPC3_20250809.pdf", width=14, height = 10)
ggplot(pca_scores_tips, aes(x = PC2, y = PC3, fill = cluster_adm)) + 
  scale_fill_manual(values = c("grupo_01" = "#9e0142", 
                               "grupo_02" = "#d53e4f", 
                               "grupo_03" = "#f46d43",
                               "grupo_04" = "#fdae61", 
                               "grupo_05" = "#fee088", 
                               "grupo_06" = "#e6f598",
                               "grupo_07" = "#abdda4", 
                               "grupo_08" = "#66c2a5", 
                               "grupo_09" = "#3288bd",
                               "grupo_10" = "#5E4FA2")) +
  xlab("PC 2: 19%") +
  ylab("PC 3: 7%") +
  geom_hline(yintercept = 0, lty = 2) +  
  geom_vline(xintercept = 0, lty = 2) +  
  geom_point(shape = 21,alpha = 3, size = 4, ) 


#dev.off()

pdf("./02_out/PLOT_PCAmorfoPC1vsPC3_20250809.pdf", width=14, height = 10)
ggplot(pca_scores_tips, aes(x = PC1, y = PC3, fill = cluster_adm)) + 
  scale_fill_manual(values = c("grupo_01" = "#9e0142", 
                               "grupo_02" = "#d53e4f", 
                               "grupo_03" = "#f46d43",
                               "grupo_04" = "#fdae61", 
                               "grupo_05" = "#fee088", 
                               "grupo_06" = "#e6f598",
                               "grupo_07" = "#abdda4", 
                               "grupo_08" = "#66c2a5", 
                               "grupo_09" = "#3288bd",
                               "grupo_10" = "#5E4FA2")) +
  xlab("PC 1: 46%") +
  ylab("PC 3: 7%") +
  geom_hline(yintercept = 0, lty = 2) +  
  geom_vline(xintercept = 0, lty = 2) +  
  geom_point(shape = 21,alpha = 3, size = 4, ) 



#dev.off()
