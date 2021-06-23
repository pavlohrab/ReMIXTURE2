#Geospatial plots

#Table of lat("y")/long("x") positions for each region, and their chosen colour for the plot (recommend hex format [https://www.google.com/search?q=color+picker]). Format (example):

centres <- data.table::fread("C://Users/Tim/Dropbox/remixture_data/geo_centres.csv",col.names=c("region","x","y","col"))
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
circle <- function(x=0,y=0,rad=1,n_pts=200){
  theta <- seq(from=0,to=((2*pi)-(2*pi)/n_pts),length.out=n_pts)
  data.table(
    x = y+sin(theta)*rad,
    y = x+cos(theta)*rad
  )
}

#Plot the view of the globe. Use me to optimise the angle your globe will take (see geom_sf() help for details)
p <- ggplot2::ggplot(data = world) +
  ggplot2::geom_sf(lwd=0.05) +
  ggplot2::coord_sf(crs = "+proj=laea +lat_0=30 +lon_0=0 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")
p

#create circle data
centres[,size:=counts[p1==region & p2==region]$count,by=region]
centres[,size:=size %>% scale_between(2,7)]
centres <- centres[!is.na(size)]

############################################################################################
## JUST FOR THIS BIT WE NEED TO EXPLICITLY EXPORT DATA.TABLE FROM DATA.TABLE::DATA.TABLE
cdat <- plyr::ldply(1:nrow(centres),function(i){
  c <- circle(centres$x[i],centres$y[i],centres$size[i])
  c[,region:=centres$region[i]]
}) %>% data.table::setDT()
############################################################################################

cdat[region=="Africa"]

#Check they plot ok. Use to adjust colours and sizes
#p + geom_spatial_polygon(data=cdat[c((.N-199):.N)],aes(x=x,y=y,fill=region,group=region),alpha=0.8) #May generate a harmless warning for leaving a stat_spatial_identity() argument default

#Normalise the nearest-neighbour counts (not strictly necessary but, for sanity purposes)
cnormed[,id:=1:.N]


cnormed <- centres[,.(p1=region,x1=x,y1=y)][cnormed,on="p1"]
cnormed <- centres[,.(p2=region,x2=x,y2=y)][cnormed,on="p2"]

#compile data for lines joining regions
ldat <- cnormed[p1 != p2][, data.table(
  region = rep(p1,2),
  x = c(x1,x2),
  y = c(y1,y2),
  count = as.numeric(rep(prop,2))
) , by="id" ]

#Do the plots
data.table::setkey(ldat,region)
for(i in unique(ldat$region)){
  P <- p +
    ggspatial::geom_spatial_path(data=ldat[region==i],ggspatial::aes(x=y,y=x,alpha=count,size=count),colour=centres[region==i]$col,lineend="round") +
    ggspatial::geom_spatial_polygon(data=cdat[region==i],ggspatial::aes(x=x,y=y,group=region),fill=centres[region==i]$col) +
    ggplot2::theme(legend.position = "none")
  print(P)
}

i <- "Africa"

#... or in case you'd like to save them
# dev.off()
# pdf("ReMIXTURE_geospatial.pdf",height=5,width=5,onefile = TRUE)
# for(i in unique(ldat$region)){
#   P <- p +
#     geom_spatial_path(data=ldat[region==i],aes(x=y,y=x,alpha=count,size=count),colour=centres[region==i]$col,lineend="round") +
#     geom_spatial_polygon(data=cdat[region==i],aes(x=x,y=y,group=region),fill=centres[region==i]$col) +
#     theme(legend.position = "none")
#   print(P)
# }
# dev.off()


#Plot superimposition of top two (or more) RGOs for each region
topn <- 2 #SET: How many top RGOs?
ldat_top <- ldat[,.SD[order(-count)][1:(2*topn)],by="region"]

P <- p
for(i in unique(ldat_top$region)){
  P <- P +
    ggspatial::geom_spatial_path(data=ldat_top[region==i],ggspatial::aes(x=y,y=x,alpha=count,size=count),lineend="round") +
    #geom_spatial_polygon(data=cdat[region==i],aes(x=x,y=y,group=region),fill="") +
    ggplot2::theme(legend.position = "none")
}
print(P)
#To plot :
# dev.off()
# pdf("ReMIXTURE_geospatial_jux.pdf",height=5,width=5,onefile = TRUE)
# print(P)
# dev.off()