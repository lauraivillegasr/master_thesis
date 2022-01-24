# World map is available in the maps package
library(maps)

# No margin
par(mar=c(0,0,0,0))

# World map
#map('world',
    #col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
    #mar=rep(0,4),border=0, ylim=c(-90,90) 
#)
locations_pi$θ_π
my_nem_map <- map_data("world")
library(dplyr)
locations_pi <- read.csv2("~/Documents/Master_Thesis/locations_pi3.csv")
library(mapproj)
library(viridis)
##version without circle legend, with strain names and fixed horizontal lines ###
##real p.davidi and JB051 isolation location provided###
p<-  ggplot() +
  coord_quickmap() +
  geom_polygon(data = my_nem_map, aes(x=long, y = lat, group = group), fill="grey43", alpha=0.3) +
  geom_point(data=locations_pi, aes(x=long, y=lat, size=θ_π, color=θ_π)) +
  geom_text_repel( data=locations_pi %>% arrange(θ_π) %>% tail(10), aes(x=long, y=lat, label=strain), size=2, fontface = "bold") +
  scale_size_continuous(range=c(1,10)) +
  scale_color_viridis(trans="log") +
  theme_void() +
  guides(size = FALSE)+
  ggtitle("θ_π for asexual and sexual nematode populations")+
  theme(plot.title = element_text(hjust = 0.5))
p <- p+ geom_text_repel (data=locations_pi %>% arrange(θ_π) %>% tail(10), aes(x=long, y=lat, label=strain), size=2, fontface = "bold")
p


