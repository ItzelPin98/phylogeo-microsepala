# Hacer mapa de datos


## 0. cargar librerías
## 1. cargar datos
## 2. procesameinto: revisar que no haya NA
## 3. hacer mapa
## 4. guardar mapa

#####################
## 0. cargar librerías
#####################
### instalar librerías
install.packages("sf")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("ggplot2")

### cargar librerías
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggplot2")


#####################
## 1. cargar datos
#####################

list.files()

dat <- read.csv("DAT_georef_all_curated_IAPR_20240518.csv")
dim(dat) # 435 13
colnames(dat)

unique(dat$stateProvince) #para saber que estados tenemos y determinar los colores (posteriorimente se cambiaran por sitio)
# Veracruz, Oaxaca, Guerrero, Michoacan, Colima, Jalisco, Nayarit, Sinaloa, Chiapas, Puebla

## seleccionamos solo los registros que tienen SI en la columna de retener

dat1 <- dat[dat$retener == "si",]
dim(dat1) # 206 13
head(dat1)

## guardamos estos datos curados

write.csv(dat1, "DAT_georef_all_20240518.csv")

####### mapa, sin usar sf,
world <- ne_countries(scale="medium", returnclass="sf")
class(world)

ggplot(data = world) +
  geom_sf(fill = "white", color = "gray") +
  geom_point(data = dat1, aes(x = decimalLongitude, y = decimalLatitude), shape = 21, size = 4, color = "black", fill = "purple") +
  labs(x = "", y = "") +
  #scale_color_manual(values=c("Veracruz" = "#90CFFF", "Oaxaca" = "#FF6A8B", "Guerrero" = "#C774E7", "Michoacan" = "#FF6AD5", "Colima" = "#20DE8B","Jalisco"= "#CCDE8B", "Nayarit" = "#C874AA", "Sinaloa" = "#FFDE8B", "Chiapas" = "#966BFF", "Puebla"= "#FFA88D")) + 
  #scale_fill_manual(values=c("Veracruz" = "#90CFFF", "Oaxaca" = "#FF6A8B", "Guerrero" = "#C774E7", "Michoacan" = "#FF6AD5", "Colima" = "#20DE8B","Jalisco"= "#CCDE8B", "Nayarit" = "#C874AA", "Sinaloa" = "#FFDE8B", "Chiapas" = "#966BFF", "Puebla"= "#FFA88D")) +
  coord_sf(xlim = c(-117.06203, -86.71082), ylim = c(14.61324, 32.61718), expand = FALSE)

#dev.off()

### guardar el mapa

pdf("MAP_microsepala_20240518.pdf", width=14, height = 10)
