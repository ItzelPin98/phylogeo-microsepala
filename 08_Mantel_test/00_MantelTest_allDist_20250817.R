######################### MANTEL TEST #######################################################


## Load library
library(vegan)
library(dplyr)

##############################################################################################
########################### dist object for genetic data #####################################
##############################################################################################
## Load distance matrix genetic (parasitic distance)
list.files("01_data")
genetic <- read.csv ("./01_data/micro_PatristicDistMatrix_20250613.csv", row.names = 1)
head(genetic)
colnames(genetic)
dim(genetic)  #178 178                  
all(rownames(genetic) == colnames(genetic))  # Debe ser TRUE
genetic_dist <- as.dist(as.matrix(genetic))

## Estimate distance matrix using euclidean method
#genetic.dist <- vegdist(genetic, method = "euclidean")


##############################################################################################
########################### dist object for climatic data ####################################
##############################################################################################


## Load climatic distance matrix
list.files("01_data")
y <- read.csv("./01_data/DAT_climPCA_micro8vars_RADcolls_20250817.csv")
colnames(y)
dim(y) #178 18

## Subset climatic data to keep PCAs
y1 <- subset(y, select=c(16:18))

rownames(y1) <- y$tip_iapr_1

## Estimate distance matrix using euclidean method
y2 <- vegdist(y1, method = "euclidean")

##############################################################################################
########################### dist object for geographic data ##################################
##############################################################################################

############## COORDINATES

list.files("01_data")
values <- read.csv("./01_data/metadata_RADseq_micro_sinOut_250813.csv")
head(values)
colnames(values) 
dim(values) #178 23

## Subset data to keep coordinates: Latitude Longitude
coords <- subset(values, select=c("latitud", "longitud"))
rownames(coords) <- values$tip_iapr_1

# Convertir columnas latitud y longitud a numéricas
coords$latitud <- as.numeric(coords$latitud)
coords$longitud <- as.numeric(coords$longitud)
str(coords)

## Estimate distance matrix using euclidean method
coords1 <- vegdist(coords, method = "euclidean")

##############################################################################################
########################### dist object for morphology data ##################################
##############################################################################################


list.files("01_data")
morfo <- read.csv("./01_data/DAT_morfoPCA_RADcolls_toDistMor_20250616.csv")
head(morfo)
colnames(morfo) 
dim(morfo) #178 10

## Subset morpho data to keep PCAs
morfo2 <- subset(morfo, select=c(4:6))

rownames(morfo2) <- morfo$tip_iapr_1

## Estimate distance matrix using euclidean method
morfo3 <- vegdist(morfo2, method = "euclidean")



##############################################################################################
#################################### MANTEL TEST #############################################
##############################################################################################

## Genetic vs climatic
mantel(genetic_dist, y2, method="pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.01471 
#Significance: 0.318

#Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0440 0.0611 0.0777 0.0892  

## Genetic vs coords1
mantel(genetic_dist, coords1, method="pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.0989 
#Significance: 0.008  

#Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0460 0.0575 0.0701 0.0962
#Permutation: free
#Number of permutations: 999

## Genetic vs morpho
mantel(genetic_dist, morfo3, method="pearson", permutations = 999, na.rm = TRUE )

#Mantel statistic r: 0.00379 
#Significance: 0.405 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0576 0.0811 0.1281 0.2006 



## Clima vs coords
mantel(y2, coords1, method="pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.01817 
#Significance: 0.209 

#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
# 0.0284 0.0383 0.0463 0.0612 
#Permutation: free
#Number of permutations: 999


## morpho vs clima

mantel(morfo3, y2, method="pearson", permutations = 999, na.rm = TRUE)

# Mantel statistic r: 0.003808 
# Significance: 0.394 

# Upper quantiles of permutations (null model):
# 90%    95%  97.5%    99% 
# 0.0410 0.0570 0.0802 0.0950 


## morpho vs geo
mantel(morfo3, coords1, method="pearson", permutations = 999, na.rm = TRUE)

# Mantel statistic r: 0.004809 
# Significance: 0.413 

#Upper quantiles of permutations (null model):
#90%    95%  97.5%    99%
#0.0443 0.0598 0.0681 0.0800 9


##############################################################################################
################################ PARTIAL MANTEL TEST #########################################
##############################################################################################


#Genetic vs. climatic controlled by coords1
mantel.partial(genetic_dist, y2, coords1, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.01298 
#Significance: 0.31 

#Upper quantiles of permutations (null model):
 # 90%    95%  97.5%    99% 
  #0.0433 0.0597 0.0730 0.0881 
#Permutation: free
#Number of permutations: 999


#Genetic vs. coords1 controlled by climatic
mantel.partial(genetic_dist, coords1, y2, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.09866 
#Significance: 0.005 

#Upper quantiles of permutations (null model):
#     90%    95%  97.5%    99% 
# 0.0426 0.0568 0.0679 0.0831
#Permutation: free
#Number of permutations: 999


#Genetic vs. morfo controlled by geo
mantel.partial(genetic_dist, morfo3,coords1, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.003331 
#Significance: 0.369 

#Genetic vs morfo controlled by clima
mantel.partial(genetic_dist, morfo3,y2, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.003734 
#Significance: 0.39 

#Genetic vs clima controlled by morfo
mantel.partial(genetic_dist, y2,morfo3, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.0147 
#Significance: 0.319 

#Clima vs morfo controlled by geo
mantel.partial(y2,morfo3, coords1, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.004741  
#Significance: 0.421 

#Morfo vs geo controlled by clima
mantel.partial(morfo3, coords1, y2, method = "pearson", permutations = 999, na.rm = TRUE)

#Mantel statistic r: 0.003721 
#Significance: 0.422 
