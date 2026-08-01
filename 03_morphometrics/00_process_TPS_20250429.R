### En este script vamos a hacer un cva y morfoespacio para las landmarks de 
### las hojas de Ipomoea microsepala, en total deberiamos tener 190 especimenes (incluyendo los outgroups),
### pero aun no los vamos a incluir, entonces son 159
## Primero vamos a cargar los datos, transformarlos de 2D a 3D para hacer el CVA

## definimos el escritorio
setwd("~/Desktop/Proyecto_microsepala/MORFO/01_landmarks")
list.files()

# Cargamos librerias necesarias
#install.packages(c("geomorph", "Morpho", "ggplot2", "dplyr"))
library(geomorph)
library(Morpho)
library(ggplot2)
library(dplyr)

# Paso 1: Leer archivo TPS con coordenadas
list.files("./01_data")
landmarks <- readland.tps("./01_data/MicroProjectSubsetRADFinal_20250502ord.TPS",specID = c("ID"), 
                           negNA = TRUE, readcurves = TRUE, warnmsg = FALSE)


dim(landmarks) #72 landmarks 2 dimensiones 178 especimenes
head(landmarks)

## Paso 2: Estimamos las landmarks pérdidas

landmarks <- estimate.missing(landmarks, method = c("TPS"))
head(landmarks)

## Usaremos define.sliders() para definir las curvas de sliding semilandmarks a cada lado de la hoja,

sliders <- define.sliders(3:72,nsliders = 181) #definimos que landmarks son de las curvas

# Paso 3: Superposición generalizada de Procrustes (GPA) para ver las formas superimpuestas
# Alinear landmarks y semilandmarks (ProcD = TRUE alinea semilandmarks)
ProcMicro <- gpagen(landmarks, curves = sliders, ProcD = TRUE)

## guardamos procrustes para TNT

writeland.tps(ProcMicro$coords, "./02_out/ProcMicroToTNT_20250809.TPS")

# Paso 4: Visualización de variaciones de forma

# Comparar la forma consenso con una especie/forma específica (aquí la 1)

mean_shape <- mshape(ProcMicro$coords)
plotRefToTarget(mean_shape, ProcMicro$coords[,,1], method = "vector")
plot(mean_shape)

## o podemos plotear todas las formas hiperimpuestas (guardamos)

pdf("./02_out/formas_hiperimpuestas_micro_20250809.pdf")

plotAllSpecimens(ProcMicro$coords, mean = TRUE) 

dev.off()

# Evaluar si hay relacion significativa entre la forma (landmarks)
# y el tamaño del centroide

# Paso 5: Cálculo del tamaño del centroide

centroid_size <- ProcMicro$Csize

# Paso 6: Regresión de la forma sobre el tamaño (efecto del tamaño en la forma)
# Esta regresión evalúa la alometría (cambio de forma con el tamaño)

regresion <- procD.lm(ProcMicro$coords ~ centroid_size)

# Paso 6: Resumen de resultados de la regresión
summary(regresion)

               # Df      SS        MS     Rsq      F        Z Pr(>F)
#centroid_size   1 0.00269 0.0026913 0.00276 0.4867 -0.53823  0.709
#Residuals     176 0.97323 0.0055297 0.99724                       
#Total         177 0.97592 

## estadisticamente

### el modelo no cuenta con evidencia estadistica significativa de que el tamaño del centroide,
### influye en la forma (p=0.70) valor Pr(>F)

### la explicación del tamaño sobre la forma es muy baja (R2 = 0.27%) valor Rsq,
### asi que no hay alometria clara
##El valor-p (0.70) indica que la relación podría haberse observado por azar
### el valor de F tampoco es elevado, lo que refuerza que la variacion explicada es baja

### biologicamente

### la forma no varia sistematicamente con el tamaño en la muestra
### se puede considerar que la forma es independiente del tamaño, para estos datos

#El tamaño del centroide no tiene un efecto significativo sobre la forma en este análisis.
#La forma parece no depender del tamaño en este conjunto de datos, 
#al menos no de manera lineal detectable con este modelo.

