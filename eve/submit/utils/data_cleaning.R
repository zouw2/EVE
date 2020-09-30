
# impute missing values in the full data
imputeWithSummaryStat <- function(dsin, FUN=median, flag.var = '_F'){
  stopifnot(is.data.frame(dsin))
  m1 <- sapply(names(dsin), function(x) sum(is.na(dsin[[x]])))
  m1 <- m1[m1> 0]
  
  print(paste('imputing', length(m1), 'variables from', deparse(substitute(dsin)), 'with per variable', deparse(substitute(FUN))))
  
  for (v in names(m1)) {
    stopifnot(v %in% colnames(dsin))
    sel <- is.na(dsin[[v]])
    
    if(is.numeric(dsin[[v]])){
      s1 <- FUN(dsin[!sel, v])
      if(is.integer(dsin[[v]])) {
        s1 <- as.integer(round(s1))
      }
      stopifnot(!is.na(s1))
      print(paste('imputing', sum(sel),'missing', v,'with', round( s1, 2),', estimated from', sum(!sel),'observed values'))
      dsin[sel, v] <- s1
    }else{
      if(is.factor(dsin[[v]])) dsin[[v]] <- as.character(dsin[[v]])
      if(is.character(dsin[[v]])) {
        print(paste('set', sum(sel),'missing', v,'as UNK'))
        dsin[sel, v] <- 'UNK'
      }else{
        stop(paste('do not know how to handle variable', v,'as a', class(dsin[[v]])))
      }
    }
    
    if(!is.na(flag.var) && nchar(flag.var) > 0 ){ # create a flag variable
      f1 <- paste(v, flag.var, sep='')
      assertthat::assert_that(!f1 %in% colnames(dsin), msg = paste(f1, 'a flag variable to indicate missingness is already in the data matrix'))
      
      dsin[[f1]] <- as.integer(sel)
    }
    
  }
  dsin
}



#' Clustering the missing pattern in the input data and remove clusters of features. 
#'
#' @param ds1 input data. It must have row names and column names
#' @param dis_measure distance measure; since we will cluster only binary variables (missing or not), the default is jaccard dissimilarity measurements for both rows and columns
#' @param linkage the agglomeration method fed to hclust
#' @param numCluster The number of clusters for rows and columns. An assumption of this function (to remove high missingness) is that the majority of rows and columns have non-missing data; rows/columns with high missingness will be clustered togther and form small subgroups. Only the subgroup (of rows/columns) with the highest frequency will be returnes; other records are discorded  
#'
#' @return an input data with rows in the most frequent row subgroup and columns in the most frequent column subgroup
#' @export
#'
#' @examples c1 <- clusterMiss(ds1, dis_measure = list('row'='jaccard', 'col'='jaccard'), numCluster = c('row'=2, 'col'=4  ))
clusterMiss <- function(ds1, dis_measure =  list('row' = 'jaccard', 'col'='jaccard'), linkage = "complete", numCluster = c('row'=2, 'col'=2  )) {
  
  
  stopifnot(!is.null(row.names(ds1)))
  stopifnot(!is.null(colnames(ds1)))
  
  cat('raw input dimension:', dim(ds1), '\n')
  
  ds1 <- ds1[!apply(ds1, 1, function(x) all(is.na(x))) , !apply(ds1, 2, function(x) all(is.na(x)))]
  
  cat('dimension after removing full missing cols/rows:', dim(ds1), '\n')
  
  ds2 <- !is.na(ds1) 
  
  
  ds2b <- ds2[apply(ds2, 1, function(x) length( unique(x) ) > 1), 
              apply(ds2, 2, function(x) length( unique(x) ) > 1)] - 0.5 # remove rows and columns without any missing
  
  
  cat('matrix used for clustering:', dim(ds2b),'\n')
  
  cat('using', paste(dis_measure, collapse = ',') ,'distance' )
  
  dis_measure_name <-  names(dis_measure)
  dis_measure <- lapply(dis_measure_name, function(d) {
    d_text <- dis_measure[[d]]
    if(d_text == 'jaccard') {
      if(d == 'row') return( as.dist( philentropy::distance(ds2b > 0, method=d_text) ) )
      if(d == 'col') return( as.dist( philentropy::distance(t(ds2b) > 0, method=d_text) ) )
    }else{
      d_text
    }
  })
  
  names(dis_measure) <- dis_measure_name
  
  pheatmap(ds2b, color=c('yellow', 'black'),show_rownames = F, show_colnames = F, scale='none',clustering_distance_rows = dis_measure[['row']], clustering_distance_cols = dis_measure[['col']], clustering_method = linkage, cutree_rows=numCluster['row'] , cutree_cols=numCluster['col'] , main= paste('hclust of missingness (>0 : not missing)')  )
  
  pt_groups <- cutree(hclust(dist(ds2b, method= dis_measure[['row']]), method=linkage), k=numCluster['row']  )
  print('row groups:')
  print(table(pt_groups))
  
  ft_groups <- cutree(hclust(dist(t(ds2b, method= dis_measure[['col']])),method=linkage ), k=numCluster['col']  )
  print('column groups:')
  print(table(ft_groups))
  
  # find the label for the smaller group
  
  # support functions defined internally
  summary1 <- function(x, n=5){
    if(length(x) <= n) { print(x); return()}
    
    cat('summary for', length(x),'items\n')
    print( summary(x) )
    
  }
  
  smaller_group <- function(v){
    f <- table(v)
    head(names(f)[order(f)], length(f) -1 )
  }
  
  pt2exclude <- row.names(ds2b)[pt_groups %in% smaller_group(pt_groups)]
  ft2exclude <- colnames(ds2b)[ft_groups %in% smaller_group(ft_groups)]
  
  # summary for excluded
  print('summary of percent non-missing values among excluded items in the original matrix')
  summary1( apply(ds2[pt2exclude, ], 1, mean) )
  
  summary1( apply(ds2[, ft2exclude], 2, mean) )
  
  # summary for not excluded
  print('summary of percent non-missing values after removing patients and samples with high missing')
  
  filtered <- ds1[setdiff(row.names(ds1), pt2exclude), setdiff(colnames(ds1), ft2exclude)]
  
  ds2c <- !is.na(filtered)
  
  summary1( apply(ds2c, 1, mean) )
  
  summary1( apply(ds2c, 2, mean) )
  
  invisible(filtered )
}


#' Add a small positive number to zeros to enable division
#'
#' @param esIn expressionSet, with one or more elements from assayData()
#' @param target default to the 'exprs' element of assayData. This function will look for 0 in the varibles from this matrix (it will not check 0s from the other elements from assayData()). A positive value (half of smallest positive value for a variable) will be added to the corresponding variable in each element of assayData()
#' @param fList A subset of features to perform the addition.
#' @param smallV a positive value less than it will be considered zero
#' @param alterPositiveFeaturesOnly default to T: addition will be only applied to variables with all values positive
#'
#' @return a modified expressionSet where there is no zero in the \code{target} element of assayData() 
#' @export
#'
#' @examples
add2zero <- function(esIn, target='exprs', fList=c(), smallV= 1e-6, alterPositiveFeaturesOnly = T) {
  stopifnot(class(esIn) == 'ExpressionSet')
  
  a <- assayData(esIn)
  stopifnot(is.list(a))
  
  stopifnot(all(target %in% names(a)))
  
  if(length(fList)==0) fList <- featureNames(esIn)
  
  m1 <- sapply(fList, function(f){
    values <- do.call(c, lapply(target, function(x)a[[x]][f, ]))
    m <- 0
    if(alterPositiveFeaturesOnly && any(values < 0, na.rm=T)) return(m)
    sel <- abs(values) <= smallV
    if (any(sel, na.rm=T) && any(!sel, na.rm=T))  {
      m <- min(values[!sel], na.rm = T)/2
      if(!is.na(m)) print(paste('adding', m, 'to all values of', f))
    }
    m
  } )
  
  m1[is.na(m1) | is.infinite(m1)] <- 0
  
  stopifnot(all(names(m1) == featureNames(esIn)))
  
  a <- lapply(a, function(x){
    sweep(x, 1, m1, '+')
  })     
  
  assayData(esIn) <- a
  
  esIn
}