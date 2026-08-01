# Leer el archivo línea por línea
list.files()
lines <- readLines("CV_errors.txt")

# Extraer K y error usando expresiones regulares
K_values <- as.numeric(sub(".*\\(K=([0-9]+)\\).*", "\\1", lines))
errors <- as.numeric(sub(".*: ", "", lines))

# Combinar en un data.frame
cv <- data.frame(K = K_values, error = errors)

# Ordenar por K
cv <- cv[order(cv$K), ]

# Mejor valor de K
best_k <- cv$K[which.min(cv$error)]
best_error <- min(cv$error)

# Mostrar
cat("El mejor valor de K es:", best_k, "con un error de:", best_error, "\n")

# Graficar
library(ggplot2)
ggplot(cv, aes(x = K, y = error)) +
  geom_line(color = "gray40") +
  geom_point(size = 3, color = "black") +
  geom_point(data = subset(cv, K == best_k), aes(x = K, y = error), 
             color = "red", size = 4) +
  geom_text(data = subset(cv, K == best_k),
            aes(label = paste0("K=", best_k, "\nError=", round(best_error, 5))),
            vjust = -1.5, color = "red") +
  labs(title = "Error de validación cruzada (ADMIXTURE)",
       x = "Número de grupos K",
       y = "CV error") +
  theme_minimal()

## guardar



# Cargar librerías
library(ggplot2)
library(pophelper)


# Leer nombres de muestra desde el archivo .fam
list.files()
fam <- read.table("microv3_allSamples_filtered.fam", stringsAsFactors = FALSE)
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
k_target <- 2
k_index <- which(sapply(qlist, ncol) == k_target)
print(paste("Índice para K=10:", k_index))


max_assign <- apply(qmatrix, 1, which.max)
assignments <- data.frame(Sample = rownames(qmatrix),
                          AssignedCluster = max_assign)

### guardamos el archivo csv para ver a que cluster admixture pertenece cada muestra

write.csv(assignments,"asignacionClusters_K10admixture_20250723.csv")
# ==============================
# Graficar Q plots automáticamente
# ==============================
plotQ(qlist,
      imgoutput = "join",
      exportpath = getwd(),
      sortind = "all",
      showindlab = TRUE,
      indlabsize = 1,
      sharedindlab = FALSE,  # <<<< necesario para "join" + "sortind = all"
      showlegend = TRUE,
      legendkeysize = 4,
      legendtextsize = 4,
      panelspacer = 0.2,
      height = 3,
      width = 8)

# Ver número de columnas por archivo .Q
sapply(qlist, ncol)

