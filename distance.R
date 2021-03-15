#' @title distance function
#'
#' @description compute the Euclidian distance between two points
#'
#' @param df1 a \code{data.frame} with x, y, and z components in this order
#' @param df2 a \code{data.frame} with x, y, and z components in this order
#'
#' @return a \code{vector} with the distance
#' @import tidyverse
#' @import purrr
#' @import stringr
#' @import dplyr
#'
#'
#' @export
#'
#'

distance <- function(df1,df2) {
  sqrt((df2[,1]-df1[,1])^2+(df2[,2]-df1[,2])^2+(df2[,3]-df1[,3])^2)
}
