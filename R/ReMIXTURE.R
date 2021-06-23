
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
      private$validate_dm(distance_matrix)
      if( !is.null(info_table) ){
        private$validate_it(info_table)
      } else {
        warning("No info table provided. Must be inputted manually with $info_table() before $run() can be called.")
      }

      private$dm <- distance_matrix


    },



    run = function(iterations=1000, resample=F){
      #run the method to fill private$counts (define this somewhere else for clarity and call it here)
      # if resample==T, then run the resampling stuff too
      gpcol <- colnames(private$dm)
      gplist <- data.table::data.table(region=colnames(private$dm))[,.N,by=.(region)]

      #index the positions of each region group
      gplist$offset <- c(0,rle(gpcol)$lengths) %>% `[`(-length(.)) %>% cumsum %>% `+`(1)
      gplist[,idx:=1:.N]

      sampsize <- (min(table(gpcol)) * (2/3)) %>% round #SET: how many samples per iteration (from each region)

      #set up some vectors to store info later
      outsize <- iterations * sampsize * nrow(gplist)
      select <- vector(mode="integer",length=sampsize*nrow(gplist)) #to store a list of the randomly selected samples each iteration
      private$raw_out <- data.table::data.table( #to store raw output each iteration
        p1 = character(length=outsize),
        p2 = character(length=outsize),
        dist = numeric(length=outsize),
        iteration = integer(length=outsize)
      )
      insert <- 1 #a flag

      #run the iterations
      for(iteration in 1:iterations){
        #fill the `select` vector
        #dev iteration = 1
        gplist[,{
          select[(sampsize*(idx-1)+1):((sampsize*(idx-1))+sampsize)] <<- sample(N,sampsize)-1+offset
        },by="idx"] %>% invisible

        #Find closest neighbours for the selected sample, store results in output table
        rnum <- 1
        #r = dm[select,select][1,]
        apply(dm[select,select],1,function(r){
          private$raw_out$p1[insert] <<- colnames(private$dm)[select][rnum]
          private$raw_out$p2[insert] <<- colnames(private$dm)[select][which(r==min(r))[1]]
          private$raw_out$dist[insert] <<- min(r)[1]
          private$raw_out$iteration[insert] <<- iteration
          rnum <<- rnum+1
          insert <<- insert+1
        }) %>% invisible

        ce("% complete: ",(insert/outsize)*100)
      }

      #summarise the output
      private$counts <- private$raw_out[ , .(count=.N) , by=.(p1,p2) ][ is.na(count) , count:=0 ]
      data.table::setorder(private$raw_out,p1,p2,-dist)
      private$raw_out[,idx:=1:.N,by=.(p1)]

      if (resample){
        samplesize <- nu(private$raw_out$iteration)*0.1 #SET: How many items to sample each time
        nrowsit <- (nu(private$raw_out$p1)**2)
        nrowsout <- nrowsit*iterations
        #to store output
        itcount <- data.table::data.table(
          p1=character(length=nrowsout),
          p2=character(length=nrowsout),
          count=numeric(length=nrowsout),
          resamp=numeric(length=nrowsout)
        )

        #perform resampling
        for(it in 1:iterations){
          #it <- 1
          ce("It: ",it)
          selectit <- sample(unique(private$raw_out$iteration),samplesize)

          fill <- data.table::setDT(expand.grid(p1=unique(private$raw_out$p1),p2=unique(private$raw_out$p2)))
          insert <- private$raw_out[ iteration %in% selectit , .(count=.N,resamp=it) , by=.(p1,p2) ]
          insert <- insert[fill,on=.(p1,p2)]
          insert[is.na(count),count:=0]
          insert[is.na(resamp),resamp:=it]

          itcount[(nrowsit*(it-1)+1):((nrowsit*(it-1))+nrow(insert))] <- insert
        }

        #summarise output
        itcount[, pct:=(count/sum(count))*100 , by=.(resamp,p1) ]
        itcount <- itcount[, .(sd_pct=sd(pct),mean_pct=mean(pct)) , by=.(p1,p2) ]
        itcount[, description:=paste0( round(mean_pct-(2*sd_pct),digits=2)," (",round(mean_pct,digits=2),") ",round(mean_pct+(2*sd_pct),digits=2)  )  ]
        private$resample <- itcount
      }
    },




    plot_heatmap = function(){
      #produce plots
      cnormed <- data.table::copy(private$counts)[,prop:=count/sum(count),by=.(p1)]
      cnormed[p1!=p2][order(prop)]
      cm <- as.matrix(data.table::dcast(cnormed,formula=p1~p2,value.var="prop")[,-"p1"])
      rownames(cm) <- colnames(cm)
      hmplot <- pheatmap::pheatmap(cm,cluster_rows = F,cluster_cols = F)
      return(hmplot)
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

      #check if there is NAs or Inf on a diag. Convert to zeroes
      if (all(is.na(diag(in_dm)))){
        diag(in_dm) <- 0
      } else if (all(is.infinite(diag(in_dm)))){
        diag(in_dm) <- 0
      }
      #check zeroes on diagonal
      if( !all(in_dm[diag(in_dm)]==0) ){
        stop("Self-distance (i.e. distance matrix diagonals) should always be zero")
      }

      #check rows and columns are the same
      sapply(1:nrow(in_dm),function(r) { all(in_dm[r,]==in_dm[,r]) })

      #check groups have decent numbers
      #check rowsnames/colnames exist and rownames==colnames
      if (is.null(colnames(in_dm)) | is.null(rownames(in_dm))){
        stop( "Column and row names of input matrix must provide region information" )
      }
      if( !all(colnames(in_dm) == colnames(in_dm)) ) {
        stop( "Column and row names of input matrix must be the same" )
      }

      return(TRUE)
    },
    validate_it = function(in_it){
      #check all columns "region", "x"(longitude) , "y"(latitude) present and character/numeric/numeric
      #if colour not present, auto-fill and
      if( is.null(in_it$col) ){ # No colours provided --- assign!
        warning("No colour column in info_table provided. Colour will be manually added.")
        in_it[ , col := replace_levels_with_colours(region) ]
      }
    },
    raw_out = data.table(), #raw output from sampling
    counts = data.table(), #(normalised, prefereably) count data from sampling
    resample = data.table() #(normalised, prefereably) count data from sampling
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



