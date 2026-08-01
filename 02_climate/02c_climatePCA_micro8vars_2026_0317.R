###  Do climate PCA


###############
#############  cargamos los datos de clima para I. microsepala 
###############


	list.files("01_data")
	
	data <- read.csv("../CLIMA_V02/02_out/DAT_clima_RADcolls_metapobs_2026_0317.csv")
	str(data)
	head(data)
	colnames(data)
	dim(data)	 # 178 27

  # le cambiamos el nombre a las variables para que sea mas facil la visualizacion, 
	# las primeras 19 son las de BIOCLIM y la 20 es la elevación
	colnames(data)
	colnames(data)[6] <- "BIO_01"
	colnames(data)[7] <- "BIO_02"
	colnames(data)[8] <- "BIO_03"
	colnames(data)[9] <- "BIO_04"
	colnames(data)[10] <- "BIO_05"
	colnames(data)[11] <- "BIO_06"
	colnames(data)[12] <- "BIO_07"
	colnames(data)[13] <- "BIO_08"
	colnames(data)[14] <- "BIO_09"
	colnames(data)[15] <- "BIO_10"
	colnames(data)[16] <- "BIO_11"
	colnames(data)[17] <- "BIO_12"
	colnames(data)[18] <- "BIO_13"
	colnames(data)[19] <- "BIO_14"
	colnames(data)[20] <- "BIO_15"
	colnames(data)[21] <- "BIO_16"
	colnames(data)[22] <- "BIO_17"
	colnames(data)[23] <- "BIO_18"
	colnames(data)[24] <- "BIO_19"
	colnames(data)[25] <- "BIO_20"
	


###############
#############  Calculamos el PCA para clima
###############

#####  ## seleccionamos solo las variables que obtuvimos del script de correlacion; es decir, las que no estan
	## correlacionadas entre si y las que mas contribuyen al PC1 de todas las variables
###para este caso son BIO_02, BIO_04, BIO_05, BIO_13, BIO_15, BIO_18, BIO_19 y BIO_20
  dim(data) #178 27
  head(data)
 	colnames(data)
 
 	dat4pca <-   data[, c(7,9,10,18,20,23,24,25)] ## solo las de interes que contribuyen más al pc1
	head(dat4pca)

 	
	table(is.na(dat4pca))  # all FALSE 1424
 	
###    PCA clima

	clim.pca <- princomp(dat4pca, cor=TRUE)
	summary(clim.pca)
#Importance of components:
#	                          Comp.1    Comp.2    Comp.3    Comp.4
#	Standard deviation     1.6335853 1.3882486 1.1982588 0.9271728
#	Proportion of Variance 0.3335751 0.2409043 0.1794780 0.1074562
#	Cumulative Proportion  0.3335751 0.5744794 0.7539574 0.8614136


	#pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/plot_PCAclimatePCs_micro8vars_20250817.pdf")
				plot(clim.pca, main="Componentes principales del clima")
			# dev.off()

			 
			
					#pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/plot_PCAclimateBiplot_micro8vars_20250817.pdf")
				biplot(clim.pca, cex=0.5)
			# dev.off()

	# str(clim.pca)
				
## vemos los scores
	
	pca.scores <- as.data.frame(clim.pca$scores)
	# str(pca.scores)
	# rownames(pca.scores)
	rownames(pca.scores) <- rownames(dat4pca)
	head(pca.scores)
	

		colnames(dat4pca)
		rownames(dat4pca)
		head(dat4pca)
		dim(dat4pca)  # 178 8
		
	
##############
################## unir datos con componentes principales
##############

	data$n_ord <- rownames(data)
	
		pca.scores$n_ord <- rownames(pca.scores)
	
	colnames(data)
	head(pca.scores)
	dat2 <- merge(data, pca.scores, by= "n_ord")
	head(dat2)
	dim(dat2)  # 178 35
	
	names(dat2)
 	#names(dat2)[28:35] <- gsub("Comp.", "clim.PC_0", names(dat2)[28:35])

	### reordenamos las columnas
	
	dat3 <- dat2[c(1,2,3,4,5,26,27,7,9,10,18,20,23,24,25,28,29,30)]	
	
	names(dat3)
	head(dat3)
	
	
##############
################## guardamoss
##############

	capture.output(summary(clim.pca), file="../CLIMA_V02/02_out/03_data_8_bioclim_250817/PCAclimate_micro8vars_RADcolls_20250817_output.txt")
	write.csv(clim.pca$sdev, "../CLIMA_V02/02_out/03_data_8_bioclim_250817/PCAclimate_micro8vars_RADcolls_20250817_eigenvals.csv")
	write.csv(loadings(clim.pca), "../CLIMA_V02/02_out/03_data_8_bioclim_250817/PCAclimate_micro8vars_RADcolls_20250817_loadings.csv")

	write.csv(dat3, "../CLIMA_V02/02_out/03_data_8_bioclim_250817/DAT_climPCA_micro8vars_RADcolls_20250817.csv", row.names=F)
	summary(clim.pca)


##############
################## PLOTS
##############
	  library(ggplot2)
	  library(gridExtra)
	  library(patchwork)
    library(cowplot)
  	library(ggrepel)

###########
############ PC 2 vs PC 1
###########
#PC1 = 33%
#PC2 = 24%
#PC3 = 18%

	unique(dat3$cladoNew_260317)	
	colnames(dat3)
	setdiff(unique(dat3$color_clade_260317),
	        names(tu_vector_de_colores))
	
	#creamos el factor
	grupo_factor <- as.factor(dat3$cluster_fst_NEW)
	
	table(grupo_factor)
	
	#ACE ACO ATN ISN IST PCE PCH PNO POO POR PSU 
	#9   9  34   23  8  32  15  23  10   7   8 
	
	grupo_factor <- factor(dat3$cluster_fst_NEW,
	                       levels = c("PNO", "PCE", "PSU", "POR", "POO",
	                                  "ACO", "ACE", "ATN", "ISN", "IST", "PCH"))
	
	pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/PLOT_climPCA_PC1vsPC2meta_2026_0318v3.pdf", width=14, height = 10)
	
	pc12 <- ggplot(data = dat3, aes(x = Comp.1, y = Comp.2, fill = grupo_factor)) +	
	  scale_fill_manual(values = c("PNO" = "#D53E4F",
	                               "PCE" = "#FDAE61",
	                               "PSU" = "#FFEC8B",
	                               "POR" = "#FDD835",
	                               "POO" = "#CDE85F",
	                               "ACO" = "#ABDDA4",
	                               "ACE" = "#8CCF88",
	                               "ATN" = "#247067",
	                               "ISN" = "#3399FF",
	                               "IST" = "#99CCFF",
	                               "PCH" = "#5E4FA2"))  +
	  xlab("Climate PC 1 (33%)") + 
	  ylab("Climate PC 2 (24%)") +
	  geom_hline(yintercept = 0, lty = 2) +  
	  geom_vline(xintercept = 0, lty = 2) +  
	  geom_jitter(aes(fill = grupo_factor), 
	              shape = 21, alpha = 0.8, size = 6, 
	              width = 0.5, height = 0.5) +
	  theme_minimal()
	
     # pc12
 	  dev.off()
	
###########
############ PC 3 vs PC 1
###########
 	  pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/PLOT_climPCA_PC1vsPC3_20250817.pdf", width=14, height = 10)
 	  
	pc31 <- ggplot(data = dat3, aes(x = Comp.1, y = Comp.3, fill = cluster_adm)) +	
	  scale_fill_manual(values = c("Cluster_01" = "#FDD835", 
	                               "Cluster_02" = "#FFA94D", 
	                               "Cluster_03" = "#3399FF",
	                               "Cluster_04" = "#247067", 
	                               "Cluster_05" = "#8B4513", 
	                               "Cluster_06" = "#CDE85F",
	                               "Cluster_07" = "#FF9BB2", 
	                               "Cluster_08" = "#5E4FA2", 
	                               "Cluster_09" = "#EE0000",
	                               "Cluster_10" = "#ABDDA4")) +
	  xlab("Climate PC 1 (33%)") + 
	  ylab("Climate PC 3 (18%)") +
	  geom_hline(yintercept = 0, lty = 2) +  
	  geom_vline(xintercept = 0, lty = 2) +  
	  geom_jitter(aes(fill = cluster_adm), 
	              shape = 21, alpha = 0.8, size = 4, 
	              width = 0.5, height = 0.5)  
	dev.off()
	
###########
############ PC 3 vs PC 2
###########
	pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/PLOT_climPCA_PC2vsPC3_20250817.pdf", width=14, height = 10)
	
	pc32 <- ggplot(data = dat3, aes(x = Comp.2, y = Comp.3, fill = cluster_adm)) +	
	  scale_fill_manual(values = c("Cluster_01" = "#FDD835", 
	                               "Cluster_02" = "#FFA94D", 
	                               "Cluster_03" = "#3399FF",
	                               "Cluster_04" = "#247067", 
	                               "Cluster_05" = "#8B4513", 
	                               "Cluster_06" = "#CDE85F",
	                               "Cluster_07" = "#FF9BB2", 
	                               "Cluster_08" = "#5E4FA2", 
	                               "Cluster_09" = "#EE0000",
	                               "Cluster_10" = "#ABDDA4")) +
	  xlab("Climate PC 2 (24%)") + 
	  ylab("Climate PC 3 (18%)") +
	  geom_hline(yintercept = 0, lty = 2) +  
	  geom_vline(xintercept = 0, lty = 2) +  
	  geom_jitter(aes(fill = cluster_adm), 
	              shape = 21, alpha = 0.8, size = 4, 
	              width = 0.5, height = 0.5)  
	
	

#######
##### save plots 
######


## PDF con legenda
	pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/PLOT_climPCA_micro178_8vars_251003.pdf", width=14, height = 10)
	plot_grid(pc12, pc31, pc32,  labels = c('A', 'B', 'C'), ncol = 2)
	# dev.off()
	
##PDFs individuales

	rm(list = ls())
###############
#############  #############    END
###############
	
	##this for poster
	pdf("../CLIMA_V02/02_out/03_data_8_bioclim_250817/PLOT_climPCA_PC1vsPC3poster_20250825.pdf", width=14, height = 10)
	
	p13 <- ggplot(data = dat2, aes(x = clim.PC_01, y = clim.PC_03, fill = cluster_adm)) +	
	  scale_fill_manual(values = c("Cluster_01" = "#FDD835", 
	                               "Cluster_02" = "#FFA94D", 
	                               "Cluster_03" = "#3399FF",
	                               "Cluster_04" = "#247067", 
	                               "Cluster_05" = "#8B4513", 
	                               "Cluster_06" = "#CDE85F",
	                               "Cluster_07" = "#FF9BB2", 
	                               "Cluster_08" = "#5E4FA2", 
	                               "Cluster_09" = "#EE0000",
	                               "Cluster_10" = "#ABDDA4")) +
	  xlab("Clima PC 1 (33%)") + 
	  ylab("Clima PC 3 (18%)") +
	  geom_hline(yintercept = 0, lty = 2) +  
	  geom_vline(xintercept = 0, lty = 2) +  
	  geom_point(alpha = 1, size = 10, shape = 21, color = "black") +
	  theme_minimal() +
	  theme(legend.position="none") +
	  theme(panel.grid = element_blank(), 
	        panel.border = element_rect(fill= "transparent"),
	        axis.title.x = element_text(size = 28),
	        axis.title.y = element_text(size = 28),
	        axis.text.x  = element_text(size = 26),
	        axis.text.y  = element_text(size = 26))
	
	# p13g
	dev.off()
	
