#' ReMIXTURE
#'
#' Regionwise similarity analysis using a resampled nearest-neighbour method.
#'
#' @section Warning:
#' Under development.
#'
#' @return A ReMIXTURE class object.
#' @export
ReMIXTURE <- R6::R6Class(
  ################# Public ################
  public = list(
    #' @description
    #' Create a new ReMIXTURE object.
    #' @param distance_matrix An all-vs-all, full numeric distance matrix, with rownames and
    #'      colnames giving the region of origin of the corresponding individual.
    #' @param info_table A data.table rescribing the lat(y)/long(x)s of each region, with columns named "region", "x", "y", and optionally "col", to give a HEX colour to each region.
    #' @return a new `ReMIXTURE`` object.
    initialize = function(distance_matrix,info_table=NULL){ #constructor, overrides self$new
      #browser()

      #call validators for dm and it if they exist
      #browser()

      if( #lower triangular dm --- fill
        all(distance_matrix[lower.tri(distance_matrix,diag=F)]==0) & !all(distance_matrix[upper.tri(distance_matrix,diag=F)]==0)
      ){
        warning("Detected a probable triangular distance matrix as input. Zero entries in lower triangle will be filled based on the upper triangle")
        dm <- fill_lower_from_upper(distance_matrix)
      }

      if( #upper triangular dm --- fill
        !all(distance_matrix[lower.tri(distance_matrix,diag=F)]==0) & all(distance_matrix[upper.tri(distance_matrix,diag=F)]==0)
      ){
        warning("Detected a probable triangular distance matrix as input. Zero entries in upper triangle will be filled based on the lower triangle")
        distance_matrix <- fill_upper_from_lower(distance_matrix)
      }




      #call validators for dm and it if they exist
      private$validate_dm(distance_matrix)
      if( !is.null(info_table) ){
        private$validate_it(info_table)
        private$it <- info_table
        #if colour not present, auto-fill
        if( is.null(info_table$col) ){ # No colours provided --- assign!
          warning("No colour column in info_table provided. Colour will be manually added.")
          info_table[ , col := replace_levels_with_colours(region) ]
        }
      } else {
        warning("No info table provided. Must be inputted manually with $info_table() before $run() can be called.")
      }


      private$dm <- distance_matrix

    },


    #' @description
    #' Run the ReMIXTURE analysis. Requires the information table to have been provided upon initialisation or later with $info_table().
    #' @param iterations The number of samplings requested.
    #' @param resample If TRUE, will resample the iterations to establish variance in the results.
    #' @return A sense of profound satisfaction.
    run = function(iterations=1000, resample=F){
      #run the method to fill private$counts (define this somewhere else for clarity and call it here)
      # if resample==T, then run the resampling stuff too
      dm <- private$dm
      gpcol <- colnames(dm)
      gplist <- data.table::data.table(region=colnames(dm))[,.N,by=.(region)]

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
          private$raw_out$p1[insert] <<- colnames(dm)[select][rnum]
          private$raw_out$p2[insert] <<- colnames(dm)[select][which(r==min(r))[1]]
          private$raw_out$dist[insert] <<- min(r)[1]
          private$raw_out$iteration[insert] <<- iteration
          rnum <<- rnum+1
          insert <<- insert+1
        }) %>% invisible

        ce("% complete: ",round((insert/outsize)*100, 4))
      }

      #summarise the output
      private$counts <- private$raw_out[ , .(count=.N) , by=.(p1,p2) ][ is.na(count) , count:=0 ]

      data.table::setorder(private$raw_out,p1,p2,-dist)
      private$raw_out[,idx:=1:.N,by=.(p1)]


      if (resample){
        samplesize <- nu(private$raw_out$iteration)*0.9 #SET: How many items to sample each time
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


    #' @description
    #' Plot the heatmap of RGOs. Requires that $run() has been called.
    #' @return A ggplot2 plot object.
    plot_heatmap = function(){
      #produce plots
      cnormed <- data.table::copy(private$counts)[,prop:=count/sum(count),by=.(p1)]
      cnormed[p1!=p2][order(prop)]
      cm <- as.matrix(data.table::dcast(cnormed,formula=p1~p2,value.var="prop")[,-"p1"])
      rownames(cm) <- colnames(cm)
      hmplot <- pheatmap::pheatmap(cm,cluster_rows = F,cluster_cols = F)
      return(hmplot)
    },



    #' @description
    #' Plot all inter-region RGOs on a map. Requires that $run() has been called.
    #' @return A pheatmap::heatmap object
    plot_maps = function(){
      #check
      #produce plots
      cnormed <- data.table::copy(private$counts)[,prop:=count/sum(count),by=.(p1)]
      cnormed<-cnormed[p1!=p2][order(prop)]

      coords<-private$it
      coords[,size:=private$counts[p1==region & p2==region]$count,by=region]
      coords[,size:=size %>% scale_between(2,7)]
      coords <- coords[!is.na(size)]

      cnormed[,id:=1:.N]
      cnormed <- coords[,.(p1=region,x1=x,y1=y)][cnormed,on="p1"]
      cnormed <- coords[,.(p2=region,x2=x,y2=y)][cnormed,on="p2"]
      # Plot the view of the globe.
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      p <- ggplot2::ggplot(data = world) +
        ggplot2::geom_sf(lwd=0.05)

      map<-  p+
        ggplot2::geom_curve(data = cnormed, ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2, size = prop),
                   curvature = 0.5,
                   alpha = 0.5,
                   lineend = "round")+
        ggplot2::geom_point(ggplot2::aes(x=x,y=y, size = size, colour = col),data= coords)+
        ggplot2::scale_colour_manual(values = coords$col)+
        ggplot2::theme(legend.position = "none")
      return(map)
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

      #check rows and columns are the same in_dm <- N
      if ( !sapply(1:nrow(in_dm),function(r) { all(in_dm[r,]==in_dm[,r]) }) %>% all ){
        stop("Distance matrix is not diagonal")
      }

      #check groups have decent numbers
      #check rowsnames/colnames exist and rownames==colnames
      if (is.null(colnames(in_dm)) | is.null(rownames(in_dm))){
        stop( "Column and row names of input matrix must provide region information" )
      }
      if( !all(colnames(in_dm) == colnames(in_dm)) ) {
        stop( "Column and row names of input matrix must be the same" )
      }

    },
    validate_it = function(in_it){
      #check all columns "region", "x"(longitude) , "y"(latitude) present and character/numeric/numeric
      if( !data.table::is.data.table(in_it) ){
        stop("Info table must be a data.table")
      }
      if( any(!c("region","x","y") %in% colnames(in_it) ) ){
        stop("Info table must have all( c(\"regions\",\"x\",\"y\") %in% colnames(.) )")
      }
      if( !all(unique(colnames(private$dm)) %in% in_it$region) ){
        stop("All regions present in distance matrix must have entries in the info table.")
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
      } else { #validate and replace private$it
        private$validate_it(in_it)
        #if colour not present, auto-fill
        if( is.null(in_it$col) ){ # No colours provided --- assign!
          warning("No colour column in info_table provided. Colour will be manually added.")
          in_it[ , col := replace_levels_with_colours(region) ]
        }
        private$it <- in_it
      }
    }
  )




)

