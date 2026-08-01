# cargar la librería que tiene las funciones de estadística genética
library(hierfstat)

# Mostrar en qué carpeta estás trabajando (esto te ayuda a ubicar rutas relativas)
getwd()


# 1. cargar matriz genotípica en formato 012

# Aquí cargas tu matriz 012 que VCFtools generó
# Ojo: el formato 012 tiene:
# - La primera columna: el ID numérico de SNPs
# - Las filas: los individuos
# - Los valores: 0 = homocigoto referencia, 1 = heterocigoto, 2 = homocigoto alternativo
# - -1 significa missing data

list.files("01_datos")

micro <- read.delim("./01_datos/microv3_DP20_MAF002_miss30_012.012", 
  header = FALSE,        # No tiene cabeceras
  row.names = 1,         # La primera columna son IDs de SNPs
  na.strings = -1        # Los -1 se interpretan como NA
)
dim(micro) #187 29669

# Cargar metadata de individuos
list.files("02_out/01_outs_251019")
meta <- read.table("./02_out/01_outs_251019/microPobs_251019.txt",header = TRUE, sep = "\t")
dim(meta) #178 6
table(meta$cluster_fst_NEW)
#Grupo_A Grupo_B Grupo_C Grupo_D Grupo_E Grupo_F Grupo_G Grupo_H Grupo_I Grupo_J 
#23      32       8       7      10       9       9      34       8      23 
#Grupo_K 
#15 

# Cargar nombres de individuos
# Este archivo .indv lo genera VCFtools
ind <- read.delim( "./01_datos/microv3_DP20_MAF002_miss30_012.012.indv", header = FALSE)
dim(ind) #187 1

# 2. filtramos los datos geneticos con base en los datos de meta

muestras_micro <- ind$V1           # Nombres de las columnas de micro
muestras_meta <- meta$extrac_ID     # Ajusta si la columna tiene otro nombre

length(muestras_micro)  # 187
length(muestras_meta)   # 178

muestras_comunes <- intersect(muestras_meta, muestras_micro)
length(muestras_comunes)  # Debería ser 178 si todo calza perfecto

rownames(micro) <- muestras_micro
micro_filtrado <- micro[muestras_comunes,]
dim(micro_filtrado)  # 178 29669

#3.unimos meta con la matriz genética

# añadimos la columna de población al inicio de la matriz genética
colnames(meta)
micro_filtrado <- cbind(Pop = meta$cluster_fst_NEW, micro_filtrado)

#Imputar los NAs antes del análisis
#En lugar de borrar individuos o SNPs, podemos rellenar los NAs con una estrategia razonable. Algunas opciones:
#Imputación por moda (valor más común por SNP):

micro_imputado <- micro_filtrado  # Asegúrate de inicializar el objeto

for (i in 1:nrow(micro_filtrado)) {
  snp <- micro_filtrado[i, ]
  
  # Excluir NA antes de calcular la moda
  snp_no_na <- snp[!is.na(snp)]
  
  if (length(snp_no_na) == 0) {
    # Si todos son NA, puedes imputar con 0 o NA (depende de lo que prefieras)
    snp[is.na(snp)] <- 0
  } else {
    # Calcular la moda (valor más frecuente)
    moda <- as.numeric(names(sort(table(snp_no_na), decreasing = TRUE)[1]))
    snp[is.na(snp)] <- moda
  }
  
  micro_imputado[i, ] <- snp
}

# 4. codificar alelos en formato hierfstat

# hierfstat espera dos dígitos por alelo:
# Ejemplo:
# - 44 = homocigoto alelo 4 (ref)
# - 45 = heterocigoto (alelos 4 y 5)
# - 55 = homocigoto alelo 5 (alt)
# Por eso hacemos esta recodificación:

micro_imputado[micro_imputado == 0] <- 44
micro_imputado[micro_imputado == 1] <- 45
micro_imputado[micro_imputado == 2] <- 55


# 5. calcular estadísticas básicas de diversidad genética

# Calcula:
# - Heterocigosis observada (Ho)
# - Heterocigosis esperada (Hs)
# - Fis (endogamia)

div.micro <- basic.stats(micro_imputado)

# extraemos la heterocigosis observada por SNP y población

Ho <- div.micro$Ho

# promedios de Ho por población)

Ho.means <- as.data.frame(colMeans(Ho))

# extraemos la heterocigosis esperada

Hs <- div.micro$Hs

# promedios de Hs por población

Hs.means <- as.data.frame(colMeans(Hs, na.rm = TRUE))
# na.rm = TRUE elimina NAs al calcular la media

# 6. calculamos Fst entre poblaciones

fst.micro <- pairwise.neifst(micro_imputado) # Esto devuelve una matriz de distancias Fst por pares de poblaciones

# 7. bootstrapping de Fis: Estima intervalos de confianza de Fis por población
fis <- boot.ppfis(micro_imputado, nboot = 100)
FIS <- fis$fis.ci
# FIS contiene el límite inferior y superior del intervalo de confianza de Fis


# 8. combinamos resultados en un único data frame

info.pobs <- as.data.frame(
  cbind(
    Ho.means,
    Hs.means,
    FIS
  )
)

colnames(info.pobs) <- c("Ho", "He", "FisLow", "FisHigh")


# 8. guuardamos los  resultados

write.table(
  info.pobs,
  "./02_out/01_outs_251019/DiversidadMicro178k10_pops_251019.txt",
  sep = "\t",
  quote = FALSE
)

write.csv(info.pobs, "./02_out/01_outs_251019/DiversidadMicro_k10_251019.csv")
write.csv(fst.micro, "./02_out/01_outs_251019/TAB_fstMicro_allSamples_251019.csv")
###graficas 

library(reshape2)
library(ggplot2)



fst.df <- melt(fst.micro)

# Graficar


pdf("./02_out/01_outs_251019/PLOT_heatmapFst_allSamplesk10_251020.pdf", width = 14, height = 10)
ggplot(fst.df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
 # scale_fill_gradient2(low = "white", high = "hotpink", mid = "blue",
  #                     midpoint = median(fst.df$value, na.rm = TRUE),
   #                    name = "Fst") +
scale_fill_viridis_c(option = "turbo", name = "Fst", na.value = "grey90") +
  geom_text(aes(label = round(value, 5)), size = 5, color = "white") +
  theme_minimal() +
  labs(title = "Heatmap de Fst entre poblaciones",
       x = "Población",
       y = "Población") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

########## otras formas de graficar con diferentes colores

library(ggplot2)
library(viridis)



ggplot(fst.df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
 scale_fill_viridis_c(option = "viridis", name = "Fst", na.value = "grey90") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Heatmap de Fst entre poblaciones",
       x = "Población",
       y = "Población") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 14, face = "bold")
  )
dev.off()


ggplot(fst.df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = heat.colors(11), name = "Fst", na.value = "grey90") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Heatmap de Fst entre poblaciones",
       x = "Población",
       y = "Población") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 14, face = "bold")
  )


##version simposio de estudiantes

# Cargar librerías necesarias
library(ggplot2)
library(reshape2)  # para usar melt()

# Leer el archivo CSV (ajusta el separador si es necesario)
fst <- read.csv("./02_out/01_outs_251019/TAB_fstMicro_allSamples_251019.csv", header = TRUE, row.names = 1, sep = ",")

# Convertir la matriz a formato largo
fst.df <- melt(as.matrix(fst), varnames = c("Var1", "Var2"), value.name = "value")

# Graficar el heatmap
pdf("./02_out/01_outs_251019/PLOT_heatmapFst_allSamples_251120.pdf", width = 14, height = 10)

ggplot(fst.df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  #scale_fill_viridis_c(option = "turbo", name = "Fst", na.value = "grey90") +
  geom_text(data = subset(fst.df, !is.na(value)),
            aes(label = round(value, 3)),
            size = 8, color = "black") +
  scale_fill_gradient(low = "white", high = "#FF6EB4", name = "Fst", na.value = "grey90") +
  
  theme_minimal() +
  labs(title = "Heatmap de Fst entre poblaciones",
       x = "Población",
       y = "Población") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))
dev.off()
