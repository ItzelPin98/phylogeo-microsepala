
##Parte 1: graficar los cv errors para escoger el valor optimo de K

#########  Plot  CV errors admixture  ########
library(ggplot2)

###Load error data
k.error<- read.delim("./01_datos/CV_errors_131samples.txt", header = F, sep = ":")

# Extraer el número K usando expresiones regulares
k.error$K <- as.numeric(gsub(".*\\(K=(\\d+)\\).*", "\\1", k.error$V1))

# Ordenar por K
k.error <- k.error[order(k.error$K), ]

rownames(k.error)<- c("k=2", "k=3", "k=4", "k=5", "k=6", "k=7", "k=8", "k=9", "k=10", "k=11", "k=12", "k=13", "k=14", "k=15", "k=16", "k=17", "k=18","k=19", "k=20")

#plot K error
#pdf("./02_out/admixture_out_filtered_20250717/PLOT_CVadmix_250818.pdf")
e.plot<- ggplot(data=k.error, aes(x=2:20, y=V2)) + geom_point(color = "magenta") + geom_line(color = "magenta")
e.plot + xlab("k") + ylab("Error") + theme(axis.text=element_text(size=21), axis.title=element_text(size=21)) + theme_bw()
#dev.off()



# Leer nombres de muestra desde el archivo .fam
list.files()
fam <- read.table("micro_131samples.fam", stringsAsFactors = FALSE)
sample_names <- fam$V2  # La segunda columna son los IDs de muestra

# Leer los archivos .Q
files <- list.files(pattern = "*.Q")

# Extraer el número K desde el nombre de archivo, asumiendo que tienen el formato "K#.Q"
k_vals <- as.numeric(gsub("[^0-9]", "", basename(files)))

# Ordenar archivos y qlist por K
ordered_indices <- order(k_vals)
files <- files[ordered_indices]
qlist <- readQ(files)

# Agregar nombres de muestra
for (i in seq_along(qlist)) {
  rownames(qlist[[i]]) <- sample_names
}

# ==============================
# Ver a qué grupo pertenece cada muestra (por asignación mayoritaria)
# ==============================
# Elige un K (por ejemplo, el tercero = K=3 si files están en orden)
k_target <- 7
k_index <- which(sapply(qlist, ncol) == k_target)

# Extrae la matriz
qmatrix <- qlist[[k_index]]

# Asegúrate de que las filas tengan nombres de muestra
rownames(qmatrix) <- sample_names

max_assign <- apply(qmatrix, 1, which.max)
assignments <- data.frame(Sample = rownames(qmatrix),
                          AssignedCluster = max_assign)

### guardamos el archivo csv para ver a que cluster admixture pertenece cada muestra

write.csv(assignments,"asignacionClusters_K7admixture_250924.csv")
