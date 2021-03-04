#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(mclust))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))


#' Discretize function
#'
#' @description
#' 'discretize' takes an .txt input file containing a matrix of fpkm values and returns an output file with
#' the discretized values (either classes 0,1 or -1,0,1)
#'
#' @param file (Path) input .txt file containing fpkm values (Columns=Samples, Rows=Genes)
#' @param figflag Binary parameter for plotting (TRUE = plot output)
#' @param binaryflag Binary parameter for binary discretization (TRUE = 0,1 classes)
#' @param output (Path) name output file
#' @return output file containing a discrete value for every fpkm value of the input file (0,1 values for binaryflag = TRUE, else -1,0,1)
#'
#' @export
#'
discretize <- function(file, figflag, binaryflag, output){

 # read input file
 data <- read.table(file, header = TRUE)

 # path settings for output files

   # get path of input file
   target.directory <- dirname(dirname(normalizePath(file)))

   # set path for plot output files
   if (figflag){

    # set path
    plot.folder <- file.path(target.directory, "Plots") # folder name for plot outputs

    # if folder "Sample Plots" does not exist, create folder
    if (!file.exists(plot.folder)){
      dir.create(plot.folder, showWarnings = TRUE)
     }
   }

   # set path for output .txt file
   if (is.na(output)){

    output.folder <- file.path(target.directory, "DGE") # folder name for plot outputs

    if (!file.exists(output.folder)){
      dir.create(output.folder, showWarnings = TRUE)
     }
   }

   # fpkm values without headers
   fpkm <- as.matrix(data)

   signal <- log2(fpkm)  # log-transform the fpkm values
   signal[is.infinite(signal)] = -10000  # transform all -inf to -10000

   discretized.keep <- matrix(nrow = nrow(fpkm), ncol = ncol(fpkm))
   rownames(discretized.keep) <- rownames(data)
   colnames(discretized.keep) <- colnames(data)

 # discretize each sample
 for (j in 1 : dim(signal)[2]){  # for each sample

  # console output: status
  write(paste("Calculation of Sample:", colnames(signal)[j]), "")

  signal.sample <- signal[,j]  # save current sample for density estimation
  data.keep <- signal[,j]  # all values of current sample for discretization
  signal.sample <- signal.sample[signal.sample > -10000]  # remove sample values with low expression for plotting

  # signal density estimation
  est <- density(signal.sample, bw = "nrd0", adjust = 1, kernel="gaussian", n = 100)


  if (figflag){

   # plot signal density
   est.df <- data.frame(X = est$x, Y = est$y)

   p <- ggplot2::ggplot(est.df, mapping = aes(x=est$x, y=est$y))
   p <- p + geom_line(aes(y=est$y, color = "signal", linetype = "signal"), size = 1.2)
   p <- p + labs(title=paste("Density Plot for Sample:", colnames(signal)[j]),
                x ="log2(FPKM)", y = "density")

  }


   # gaussian mixture model with 2 mixture components

    signalfit <- mclust::densityMclust(signal.sample, G=2)

    # first gaussian component of mixture model
    norm1 <- signalfit$parameters$pro[1]*dnorm(est$x, signalfit$parameters$mean[1], sqrt(signalfit$parameters$variance$sigmasq[1]))

    # Second gaussian component of mixture model
    if (length(signalfit$parameters$variance$sigmasq) == 1){ # both components have same variance -> only one entry
      norm2 <- signalfit$parameters$pro[2]*dnorm(est$x, signalfit$parameters$mean[2], sqrt(signalfit$parameters$variance$sigmasq[1]))
    }else{
      norm2 <- signalfit$parameters$pro[2]*dnorm(est$x, signalfit$parameters$mean[2], sqrt(signalfit$parameters$variance$sigmasq[2]))
    }


   # plot gaussian mixture components
   if(figflag){

     p <- p + geom_line(aes(y=norm1, color = "gauss1", linetype = "gauss1"), size = 1.2)
     p <- p + geom_line(aes(y=norm2,  color = "gauss2", linetype = "gauss2"), size = 1.2)

   }

  # treshold calculation

   # starting point for index search: mu1 - 2*sigma1
    low.idx.gauss1 <- signalfit$parameters$mean[1]-2*sqrt(signalfit$parameters$variance$sigmasq[1])

   # binary index: intersection of two gaussian components
    intersection.idx.bin <- (norm1 < norm2) & (est$x >= low.idx.gauss1)  # y-value first gaussian is lower than second gaussian
   # intersection.idx2 <- (norm1 < est$y) & (est$x >= min(est$x[est$y > 0.001])) #(est$x >= low.idx.gauss1)  # y-value first gaussian is lower than signal-curve

   # maximum signal density estimation curve
    peak.max.y <- which.max(est$y)  # y-value
    peak.max.x <- est$x[peak.max.y]  # x-value
    # peak.max.idx <- c(peak.max.x,peak.max.y)

    if(binaryflag){ # calculation expression treshold

      # intersection two gaussian curves
      if (any(intersection.idx.bin)){
          intersection.point.bin <- min(est$x[intersection.idx.bin]) #leftmost intersection

         # lower bound for intersection point
         if(intersection.point.bin < min(est$x[est$y > 0.001])){
            intersection.point.bin <- peak.max.x
         }
      } else {
         intersection.point.bin <- peak.max.x
      }
    } else { # no binary flag: calculation inexpression and expression treshold

       # means of gaussian mixture components as tresholds
       intersection.point.inexp <- signalfit$parameters$mean[1]
       intersection.point.exp <- signalfit$parameters$mean[2]

       # lower bound inexpression treshold
       if(intersection.point.inexp < min(est$x[est$y > 0.001])){
          intersection.point.inexp <- min(est$x[est$y > 0.001])
       }
    }


    # plot tresholds
    if (figflag){
      if (binaryflag){
        p <- p + geom_line(aes(x =intersection.point.bin, linetype = "treshold", color = "treshold"), size = 1)
      }else{
        p <- p + geom_line(aes(x=intersection.point.inexp, color = "inexpression treshold", linetype = "inexpression treshold"), size = 1)
        p <- p + geom_line(aes(x=intersection.point.exp, color = "expression treshold",  linetype = "expression treshold"), size = 1)
      }
    }


   # set expression tresholds
    if(binaryflag){
      expression.treshold <- intersection.point.bin
    }else{
      inexpression.treshold <- intersection.point.inexp
      expression.treshold <- intersection.point.exp
    }

    # indices expression classes according to tresholds
    if (binaryflag){
      activated.exp <- which(data.keep >= expression.treshold)
      not.exp <- which(data.keep < expression.treshold)
    }else{
     activated.exp <- which(data.keep >= expression.treshold)
     not.exp <- which(data.keep <= inexpression.treshold)
     normal.exp <- which(data.keep > inexpression.treshold & data.keep < expression.treshold)
    }

    # log2(fpkm) values to discrete values
    if(binaryflag){
      data.keep[activated.exp] <- 1
      data.keep[not.exp] <- 0
    }else{
      data.keep[activated.exp] <- 1
      data.keep[not.exp] <- -1
      data.keep[normal.exp] <- 0
    }

    # save discretized gene values in output matrix
    discretized.keep[,j] <- data.keep

    # save plot
    if (figflag){

      # add legend
      if (binaryflag){
        legend.title <- "Legend"
        p <- p +
          scale_color_manual(name = legend.title,
                             values = c("treshold"="black", "gauss1"="blue", "gauss2"="red", "signal"="black")) +
          # labels = c("treshold", "gauss1", "gauss2", "signal")) +
          scale_linetype_manual(name = legend.title, values=c("treshold"="dotted",  "gauss1"="solid", "gauss2"="solid", "signal"="solid")) +
          theme(legend.direction = "vertical")
        #guides(linetype=guide_legend(keywidth = 3, keyheight = 1), colour=guide_legend(keywidth = 3, keyheight = 1))

      }else{
        legend.title <- "Legend"
        p <- p +
          scale_color_manual(name = legend.title,
                             # breaks = c("inexpression treshold", "expression treshold", "gauss1", "gauss2", "signal"),
                             values = c( "inexpression treshold" = "blue", "expression treshold" = "red", "gauss1" = "blue", "gauss2" = "red", "signal" = "black")) +
          scale_linetype_manual(name = legend.title, values = c("inexpression treshold" = "dotted", "expression treshold" = "dotted", "gauss1" = "solid", "gauss2"="solid", "signal"="solid" )) +
          theme(legend.direction = "vertical") # legend.key = element_blank()) +
        # guides(linetype=guide_legend(keywidth = 3, keyheight = 1),colour=guide_legend(keywidth = 3, keyheight = 1))

      }

      # name of individual plot file
      plot.file <- paste("Sample_",colnames(signal)[j], ".png", sep = "") # plot file names

      # set file path
      fin.path <- file.path(plot.folder, plot.file)

      # save plot
      ggplot2::ggsave(filename = fin.path)

    }

 }

  cat("Save output file...")
  if(is.na(output)){
    write.table(discretized.keep, file=file.path(output.folder, "dge.txt"), row.names=TRUE, col.names=TRUE)
  }else{
    write.table(discretized.keep, file=output, row.names=TRUE, col.names=TRUE)
  }

}

main <- function() {

  # Command line arguments

  args <- commandArgs(trailingOnly = TRUE)

  option_list <- list(
    make_option(c("-f", "--file"),
                type="character",
                dest="filename",
                default=NULL,
                help="dataset file name",
                metavar="character"),
    make_option(c("-p", "--plot"),
                dest = "figflag" ,
                default = FALSE,
                action = "store_true",
                help="Should the program plot the resulting densities of each sample? [default %default]",
                metavar="character"),
    make_option(c("-b", "--binary"),
                dest = "binaryflag" ,
                default = FALSE,
                action = "store_true",
                help="Should the program's output be a binary classification? [default %default]",
                metavar="character"),
    make_option(c("-o", "--out"),
                type="character",
                dest = "output",
                default=NA,
                help="output file name [default= %default]",
                metavar="character")
  )


  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);

  # input file name mandatory
  if (is.null(opt$file)){
    print_help(opt_parser)
    stop("The argument 'input file' needs to be supplied", call.=FALSE)
  }

  discretize(opt$filename, opt$figflag, opt$binaryflag, opt$output)

}

main()









