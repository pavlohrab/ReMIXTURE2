# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' sayname
#'
#' Says the name of the intended future of this package
#'
#' @section Warning:
#' Does nothing interesting.
#'
#' @param x Anything printable.
#' @return Nothing. Just prints nothing interesting.
#' @examples
#' sayname("Humboldt")
#' @export
sayname <- function(x="nothing") {
  print(paste0("This is ReMIXTURE: ",x))
  print(head(melt(iris)))
}



























###################################################################################################################################


#require(R6)
#require(data.table)
#require(magrittr)
#require(plyr)
#require(ggplot2)
#require(data.table)
#require(sp)
#require(rgdal)
#require(ggspatial)
#require(rnaturalearth)
#require(pheatmap)
#require(rnaturalearthdata)
#require(colorspace)

`%>%` <- magrittr::`%>%`
`data.table` <- data.table::`data.table`



### IMPORTED FUNCTIONS #####################################################
ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }
nu <-function(x){
  unique(x) %>% length
}
scale_between <- function(x,lower,upper){
  if(all(x==mean(x,na.rm=T))) return(rep(mean(c(lower,upper),na.rm=T),length(x)))
  ( x - min(x,na.rm=T) ) / (max(x,na.rm=T)-min(x,na.rm=T)) * (upper-lower) + lower
}

replace_levels_with_colours <- function(x,palette="Berlin",alpha=1,fun="diverge_hcl",plot=FALSE,newplot=TRUE){
  #require(colorspace)
  n <- nu(x[!is.na(x)])
  cols <- match.fun(fun)(n,palette = palette,alpha = alpha)
  colvec <- swap( x , unique(x[!is.na(x)]) , cols , na.replacement = NA )
  if(plot==FALSE) {
    return(colvec)
  } else {
    # null_plot(y=1:length(cols),x=rep(1,length(cols)),xaxt="n",yaxt="n")
    # text(y=1:length(cols),x=rep(1,length(cols)),labels=unique(x),col=cols)
    if(newplot) {null_plot(x=0,y=0,xaxt="n",yaxt="n",bty="n")}
    legend(x="topleft",legend=unique(x[!is.na(x)]),fill=cols,text.col=cols)
  }
}
swap <- function(vec,matches,names,na.replacement=NA){
  orig_vec <- vec
  #if(sum(! matches %in% names ) > 0 ) { stop("Couldn't find all matches in names") }
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you old bison!") }
  if(is.factor(vec)) { levels(vec) <- c(levels(vec),names,na.replacement) }
  vec[is.na(orig_vec)] <- na.replacement
  plyr::l_ply( 1:length(matches) , function(n){
    vec[orig_vec==matches[n]] <<- names[n]
  })
  vec
}
null_plot <- function(x,y,xlab=NA,ylab=NA,...){
  plot(NULL,xlim=range(x,na.rm=T),ylim=range(y,na.rm=T),xlab=xlab,ylab=ylab,...)
}
############################################################################

###################################################################################################################################
##################################################### Main Class ##################################################################
###################################################################################################################################

ReMIXTURE <- R6::R6Class(

  ################# Public ################
  public = list(
    initialize = function(distance_matrix,info_table=NULL){ #constructor, overrides self$new
      browser()

      if( #lower triangular dm --- fill
        all(distance_matrix[lower.tri(distance_matrix,diag=F)]==0) & !all(distance_matrix[upper.tri(distance_matrix,diag=F)]==0)
      ){
        warning("Detected a probable triangular distance matrix as input. Zero entries in lower triangle will be filled based on the upper triangle")
        dm[lower.tri(dm)] <- dm[upper.tri(dm)]
      }

      if( #upper triangular dm --- fill
        !all(distance_matrix[lower.tri(distance_matrix,diag=F)]==0) & all(distance_matrix[upper.tri(distance_matrix,diag=F)]==0)
      ){
        warning("Detected a probable triangular distance matrix as input. Zero entries in upper triangle will be filled based on the lower triangle")
        dm[upper.tri(dm)] <- dm[lower.tri(dm)]
      }




      #call validators for dm and it if they exist
      validate_dm(distance_matrix)
      if( !is.null(info_table) ){
        validate_it(info_table)
      } else {
        warning("No info table provided. Must be inputted manually with $info_table() before $run() can be called.")
      }


    },



    run = function(iterations=1000,resample=F){
      #run the method to fill private$counts (define this somewhere else for clarity and call it here)
      # if resample==T, then run the resampling stuff too
    },




    plot_heatmap = function(){
      #produce plots
    },




    plot_maps = function(){
      #check
      #produce plots
    }
  ),




  ################# Private ################
  private = list(
    dm = matrix(), # a distance matrix with rownames and colnames giving regions
    it = data.table::data.table(), # an info table with columns "region", "lat" , "long" , and optionally "colour"
    iterations = NA_integer_, # a record of the number of iterations used for the
    validate_dm = function(in_dm){
      #check matrix is a legit distance matrix
      if( !is.matrix(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a matrix (class(in_dm)==",class((in_dm)),")") )
      }
      if( ncol(in_dm) != nrow(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a square matrix") )
      }

      #check zeroes on diagonal
      if( !all(in_dm[diag(in_dm)]==0) ){
        stop("Self-distance (i.e. distance matrix diagonals) should always be zero")
      }

      #check rows and columns are the same
      sapply(1:nrow(in_dm),function(r) { all(in_dm[r,]==in_dm[,r]) })

      #if matrix is triangular, fix it
      #check groups have decent numbers
      #check rowsnames/colnames exist and rownames==colnames
      if (is.null(colnames(in_dm) | is.null(rownames(in_dm)))){
        stop( "Column and row names of input matrix must provide region information" )
      }
      if( !all(colnames(in_dm) == colnames(in_dm)) ) {
        stop( "Column and row names of input matrix must be the same" )
      }

      return(TRUE)
    },
    validate_it = function(in_it){
      #check all columns "region", "lat" , "long" present and character/numeric/numeric
      #if colour not present, auto-fill and
      if( is.null(info_table$col) ){ # No colours provided --- assign!
        warning("No colour column in info_table provided. Colour will be manually added.")
        info_table[ , col := replace_levels_with_colours(region) ]
      }
    },
    raw_out = data.table::data.table(), #raw output from sampling
    counts = data.table::data.table() #(normalised, prefereably) count data from sampling
  ),





  ################# Active ################
  active = list( #functions that look like vars. mostly for getters and setters of privates, since they can perform checks
    distance_matrix = function(in_dm){
      if(missing(in_dm)){
        dm
      } else {
        warning("Distance matrix cannot be set after initialisation")
      }
    },
    info_table = function(in_it){
      if(missing(in_it)){
        return(it)
      } else { #validate and replace private$ct
        if( validate_it(in_it) ) {
          private$it <- in_it
        }
      }
    }
  )




)



#################### WORKSPACE ########################


dm <- readRDS("../../dmfiltered.Rds") #a matrix of similarities (e.g., IBS scores, where higher=more similar) between individuals, where the row and column names give the regions from which the sample comes. Naturally, the individuals should be in the same order along rows and columns. Also, importantly, the regions MUST be grouped together, e.g. Regions(rows) = A,A,A,C,C,D,D,D,D,D,B,B and NOT A,B,D,A,C,C,A, ...
dm[lower.tri(dm)] <- dm[upper.tri(dm)] #assure it's symmetrical


all(colnames(dm)==rownames(dm)) #should be true
gpcol <- colnames(dm)
gplist <- data.table::data.table(region=colnames(dm))[,.N,by=.(region)]

#index the positions of each region group
gplist$offset <- c(0,rle(gpcol)$lengths) %>% `[`(-length(.)) %>% cumsum %>% `+`(1)
gplist[,idx:=1:.N]

nits <- 10 #SET: how many iterations?
sampsize <- (min(table(gpcol)) * (2/3)) %>% round #SET: how many samples per iteration (from each region)

#set up some vectors to store info later
outsize <- nits * sampsize * nrow(gplist)
select <- vector(mode="integer",length=sampsize*nrow(gplist)) #to store a list of the randomly selected samples each iteration
rawoutput <- data.table::data.table( #to store raw output each iteration
  p1 = character(length=outsize),
  p2 = character(length=outsize),
  dist = numeric(length=outsize),
  iteration = integer(length=outsize)
)
insert <- 1 #a flag

#run the iterations
for(iteration in 1:nits){
  #fill the `select` vector
  #dev iteration = 1
  gplist[,{
    select[(sampsize*(idx-1)+1):((sampsize*(idx-1))+sampsize)] <<- sample(N,sampsize)-1+offset
  },by="idx"] %>% invisible

  #Find closest neighbours for the selected sample, store results in output table
  rnum <- 1
  #r = dm[select,select][1,]
  apply(dm[select,select],1,function(r){
    rawoutput$p1[insert] <<- colnames(dm)[select][rnum]
    rawoutput$p2[insert] <<- colnames(dm)[select][which(r==min(r))[1]]
    rawoutput$dist[insert] <<- min(r)[1]
    rawoutput$iteration[insert] <<- iteration
    rnum <<- rnum+1
    insert <<- insert+1
  }) %>% invisible

  ce("% complete: ",(insert/outsize)*100)
}

#summarise the output
counts <- rawoutput[ , .(count=.N) , by=.(p1,p2) ][ is.na(count) , count:=0 ]

#look at the raw counts ("summary") matrix to check it all seems to have gone ok
data.table::dcast(counts,formula=p1~p2,value.var="count")

data.table::setorder(rawoutput,p1,p2,-dist)
rawoutput[,idx:=1:.N,by=.(p1)]

#saveRDS(counts,"COUNTS.Rds") #it's a good idea to save the output, these runs take a while
#saveRDS(rawoutput,"RAW_OUTPUT.Rds") #it's a good idea to save the output, these runs take a while


#heatmap
cnormed <- data.table::copy(counts)[,prop:=count/sum(count),by=.(p1)]
cnormed[p1!=p2][order(prop)]
cm <- as.matrix(data.table::dcast(cnormed,formula=p1~p2,value.var="prop")[,-"p1"])
dim(cm)
rownames(cm) <- colnames(cm)
hmplot <- pheatmap::pheatmap(cm,cluster_rows = F,cluster_cols = F)
hmplot

#Run me to save it
# dev.off()
# pdf("ReMIXTURE_heatmap.pdf",height=7,width=7,onefile = TRUE)
# print(hmplot)
# dev.off()

#significance testing by resampling from the raw output
nits <- 50 #SET: How many samples
samplesize <- nu(rawoutput$iteration)*0.1 #SET: How many items to sample each time
nrowsit <- (nu(rawoutput$p1)**2)
nrowsout <- nrowsit*nits
#to store output
itcount <- data.table::data.table(
  p1=character(length=nrowsout),
  p2=character(length=nrowsout),
  count=numeric(length=nrowsout),
  resamp=numeric(length=nrowsout)
)

#perform resampling
for(it in 1:nits){
  #it <- 1
  ce("It: ",it)
  selectit <- sample(unique(rawoutput$iteration),samplesize)

  fill <- data.table::setDT(expand.grid(p1=unique(rawoutput$p1),p2=unique(rawoutput$p2)))
  insert <- rawoutput[ iteration %in% selectit , .(count=.N,resamp=it) , by=.(p1,p2) ]
  insert <- insert[fill,on=.(p1,p2)]
  insert[is.na(count),count:=0]
  insert[is.na(resamp),resamp:=it]

  itcount[(nrowsit*(it-1)+1):((nrowsit*(it-1))+nrow(insert))] <- insert
}

#summarise output
itcount[, pct:=(count/sum(count))*100 , by=.(resamp,p1) ]
itcount <- itcount[, .(sd_pct=sd(pct),mean_pct=mean(pct)) , by=.(p1,p2) ]
itcount[, description:=paste0( round(mean_pct-(2*sd_pct),digits=2)," (",round(mean_pct,digits=2),") ",round(mean_pct+(2*sd_pct),digits=2)  )  ]
itcount
#To save the result ...
#write.csv(itcount,"MIX_data_stdev.csv")




#Geospatial plots

#Table of lat("y")/long("x") positions for each region, and their chosen colour for the plot (recommend hex format [https://www.google.com/search?q=color+picker]). Format (example):

centres <- data.table::fread("../../geo_centres.csv",col.names=c("region","x","y","col"))
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

