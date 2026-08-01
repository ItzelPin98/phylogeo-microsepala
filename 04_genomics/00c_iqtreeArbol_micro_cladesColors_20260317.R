### Script para hacer PGLS en atributos de rango y amp climatica
#1. hacer el matching de las especies con el arbol
#2. pruning del arbol
#3. visualizar traits en arbol
#4. PGLS


#cargar librerias
   
    library(ape)
    library (phytools)
	# require(plyr)
	library (caper)	    
    # library(phangorn)
    # library(rBt)
    
list.files("01_data")

## cargar árbol
		arbol <- read.tree("./01_data/iqtree_trees/modelGTRFIG4_micro.contree")
		arbol  # 187 tips and 185 internal nodes.
		str(arbol)
		plot(arbol)
# Ver etiquetas de nodos internos (bootstrap)
		arbol$node.label
		length(arbol$node.label)  # ¿Coincide con el número de nodos internos? si 185

## cargar datos con colores asignados
		list.files("01_data")
		dat <- read.csv("./01_data/metadata_RADseq_microsepalaNewColorsTree_2026_0317.csv")
		dat1 <- read.csv("./01_data/metadata_RADseq_microNewColorsTreeALLouts_2026_0318.csv")
		head(dat1)
		dim(dat1)  # 191 29
    colnames(dat1)

##eliminamos las que no tenemos en el arbol "eIm_021"   "eIm_083"   "eIm_142"  "eIm_144"  "reference"

# creamos un vector de IDs a eliminar
		muestras_a_eliminar <- c("eIm_021", "eIm_083", "eIm_142", "reference", "e144_gro_Zihuatanejo")
		
# eliminamos las filas filas cuyo extrac_ID esté en el vector
		dat <- dat[!(dat$extrac_ID %in% muestras_a_eliminar), ]
		dat1 <- dat1[!(dat1$extrac_ID %in% muestras_a_eliminar), ]
		
## comparar conjuntos de puntas de arbol y extraccion
		setdiff(arbol$tip.label, dat1$extrac_ID) # character(0)
		setdiff(dat$extrac_ID, arbol$tip.label) # character(0)
		setdiff(dat1$extrac_ID, arbol$tip.label) # character(0)
		
## ¿están en el mismo orden las puntas del arbol y los metadatos?
		identical(arbol$tip.label, dat1$extrac_ID)  # FALSE
	
## organizar datos en función de puntas del arbol

		dat2 <- dat1[match(arbol$tip.label, dat1$extrac_ID),] 
		dim(dat2) #187 29
		head(dat2)
		
		identical(arbol$tip.label, dat2$extrac_ID)  # TRUE


## ahora sí podemos sustituir las puntas del arbol

### puntas arbol
		arbol$tip.label <- dat2$tip_iapr_1

### enraizar con minutiflora
		arbol <- root (arbol, "bract_e002_mor")
plot(arbol)

#Podar tips indeseados
		puntas_a_eliminar <- c(
		  "bract_e002_mor", "bract_e002_mor", "bract_e002_mor",
		  "e144_gro_Zihuatanejo", "purga_e154_ver", "purga_e155_ver",
		  "dumosa_e156_ver", "dumosa_e157_ver")
		puntas_a_eliminar <- c("e144_gro_Zihuatanejo")
		
		# Eliminar duplicados si hay
		puntas_a_eliminar <- unique(puntas_a_eliminar)
		
		# Podar el árbol
		arbol2 <- drop.tip(arbol, puntas_a_eliminar)
		plot(arbol2)

		#Agregar coloración de tips
		colores<- c(
		  "PNO" = "#D53E4F",
		  "PCE" = "#FDAE61",
		  "PSU" = "#FFEC8B",
		  "POR" = "#FDD835",
		  "POO" = "#CDE85F",
		  "ACO" = "#ABDDA4",
		  "ACE" = "#7FFF00",
		  "ATN" = "#247067",
		  "ISN" = "#99CCFF",
		  "IST" = "#3399FF",
		  "PCH" = "#5E4FA2"
		)
#dat2$color_clade_260317

#		tip_colors <- setNames(dat2$color_clade_260317, dat2$tip_iapr_1)
#		colors_in_order <- tip_colors[arbol2$tip.label]
		#	which(is.na(colors_in_order) | colors_in_order == "")
#		colnames(dat2)
		# 1. Creamos un vector nombrado para asociar cada 'tip' con su 'cluster_fst'
		tip_clusters <- setNames(dat2$cluster_fst, dat2$tip_iapr_1)
		
		# 2. Ordenamos las categorías siguiendo estrictamente el orden de las puntas del árbol
		clusters_in_tree_order <- tip_clusters[arbol2$tip.label]
		
		# 3. Traducimos las categorías a sus colores correspondientes usando tu vector 'colores'
		colors_in_order <- colores[clusters_in_tree_order]
		
		# Opcional: Verificar si alguna punta se quedó sin color (NA)
		# which(is.na(colors_in_order))		
	#	 arbol2 <-ladderize(arbol2, r=F)
		
##		ver si tenemos info de longitudes de rama

		arbol$edge.length #si tenemos
		
		ultra_tree <- chronos(arbol2) 

		is.ultrametric(ultra_tree) #true

		which(colors_in_order == "")
### dibujar arbol
		par(mfrow=c(1,1), mar=c(1.5,1,0.5,0.5), mgp=c(0.5,0.5,0.5))
		#plot(arbol2)

## con puntas organizadas hacia la izq (laderizado; puntas en cursivas),	y los tips de colores			

		#pdf("./02_out/FIG_arbol_iqtreeMicro_colorsCladeNEW_2026_0317.pdf", height=15, width=8.5)		
		#pdf("./02_out/FIG_arbol_iqtreeMicro_CCallOuts_2026_0702.pdf", height=15, width=8.5)	
		pdf("./02_out/FIG_arbol_MicroALLouts2026_0702.pdf", height=15, width=8.5)
		
		plot(ladderize(ultra_tree,r=F), lwd=1.5, cex=0.4, font = 4,
		     label.offset=0.0009, tip.color = colors_in_order)
		add.scale.bar(-.002,-4, cex = 0.4, font = 2, lcol="gray30", lwd=1, col = "gray30")
		# Mostrar valores de bootstrap (opcional: filtrar por ≥70)
		bootstrap_vals <- as.numeric(arbol2$node.label)
		bootstrap_filtrados <- ifelse(bootstrap_vals >= 50, bootstrap_vals, "")
		nodelabels(bootstrap_filtrados, frame="none", cex=0.4, adj=c(1, -0.3))
		#nodelabels(arbol2$node.label, frame = "none", cex=0.4, adj=c(1, -0.3))
		
		dev.off() 

		colors = ["#D53E4F",
		          "#FDAE61",
		          "#FFEC8B",
		          "#FDD835",
		          "#CDE85F",
		          "#ABDDA4",
		          "#7FFF00",
		          "#247067",
		          "#3399FF",
		          "#99CCFF",
		          "#5E4FA2"]