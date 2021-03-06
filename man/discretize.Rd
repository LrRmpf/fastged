% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discretization.R
\name{discretize}
\alias{discretize}
\title{Discretize function}
\usage{
discretize(file, figflag, binaryflag, output)
}
\arguments{
\item{file}{(Path) input .txt file containing fpkm values (Columns=Samples, Rows=Genes)}

\item{figflag}{Binary parameter for plotting (TRUE = plot output)}

\item{binaryflag}{Binary parameter for binary discretization (TRUE = 0,1 classes)}

\item{output}{(Path) name output file}
}
\value{
output file containing a discrete value for every fpkm value of the input file (0,1 values for binaryflag = TRUE, else -1,0,1)
}
\description{
'discretize' takes an .txt input file containing a matrix of fpkm values and returns an output file with
the discretized values (either classes 0,1 or -1,0,1)
}
