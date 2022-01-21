# Function that optimizes hyperparameters by spectral clustering.

specc_rec <- function(
                      pred_results, hyperparameter_names,
                      # Function parameters (default values)
                      eval_col="sum_smape_vc",
                      norm_fun = c("none","minmax","zscore"),
                      cl_sel_fun = c("min_error","centrality"),
                      rec_fun = c("mode","min_error","centrality"),
                      centr_fun = NULL, centers=2, 
                      # Spectral clustering specific parameters
                      kernel = c("rbfdot","polydot","vanilladot","tanhdot","laplacedot","besseldot","anovadot","splinedot","stringdot"),
                      kpar = "automatic",
                      nystrom.red = FALSE, 
                      nystrom.sample = dim(x)[1]/6, iterations = 200, 
                      mod.sample = 0.75, 
                      na.action = na.omit) {
  
  library(kernlab)
  
  set.seed(1234)
  
  #===== Extracting hyperparameters from the dataset =====
  
  hyperparameters <- select(pred_results, all_of(c(hyperparameter_names, eval_col)))
  
  #===== Normalizing table columns =====
  norm_fun_str <- match.arg(norm_fun)
  if(norm_fun_str == "none") norm_fun <- NULL
  else if(norm_fun_str == "minmax") norm_fun <- function(x) if(max(x)==min(x)) 0 else (x - min(x)) / (max(x) - min(x))
  else if(norm_fun_str == "zscore") norm_fun <- base::scale
  
  norm <- function(data){
    do.call(cbind,lapply(data,
                         function(col) {
                           n_col <- col
                           if (is.character(col)) n_col <- as.numeric(as.factor(col))
                           if(!is.null(norm_fun)) n_col <- norm_fun(n_col)
                           n_col[is.nan(n_col)] <- 0
                           n_col
                         }
                  )
            )
  }
  n_data <- norm(hyperparameters)
  
  
  #===== Spectral clustering =====
  kernel <- match.arg(kernel)
  
  specc_f <- kernlab::specc(as.matrix(n_data), 
                            centers=centers, kernel = kernel, kpar = kpar, 
                            nystrom.red = nystrom.red, nystrom.sample = nystrom.sample, 
                            iterations = iterations, mod.sample = mod.sample, na.action = na.action)
  
  kernel_fun <- specc_f@kernelf
  
  pred_results$cluster <- specc_f@.Data
  hyperparameters$cluster <- specc_f@.Data
  
  
  #===== Generating the graph =====
  get_graph <- function(data,kernel){
    adjmatrix <- kernlab::kernelMatrix(kernel,as.matrix(data)) + 1
    
    #print(adjmatrix)
    
    g <- igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    igraph::E(g)$weight <- igraph::E(g)$weight - 1
    
    return(g)
  }
  
  g <- get_graph(n_data,kernel_fun)
  
  
  #===== Selection of the candidate cluster of hyperparameters =====
  
  # Selects the group with the smallest mean value of forecast error
  min_error_cl <-  function(pred_results){
    eval_cl_errors <- aggregate(pred_results[eval_col], list(pred_results$cluster), base::mean)
    names(eval_cl_errors) <- c("cluster","mean_error")
    which.min(eval_cl_errors$mean_error)
  }
  # Selects the most central combinations from each cluster, and chooses the cluster
  # whose most central combination has the smallest prediction error.
  central_cl <-  function(pred_results){
    #browser()
    centrals <- sapply(1:centers, function(clstr) {
          #browser()
          clstr_data <- pred_results[pred_results$cluster == clstr,]
          n_clstr_data <- norm(clstr_data)
          
          clstr_g <- get_graph(n_clstr_data,kernel_fun)

          centr_v <- which.max(centr_fun(clstr_g))
          #print(clstr_data[centr_v, eval_col])
          clstr_data[centr_v,]
        }
    )
    #print(centrals)
    #n_centrals <- norm(centrals)
    #centrals_g <- get_graph(n_centrals)
    #centr_cl <- which.max(1/centr_fun(centrals_g))
    centr_cl <- which.min(centrals[eval_col, ])
    #print(centr_cl)
    centr_cl
  }
  
  cl_sel_fun_str <- match.arg(cl_sel_fun)
  if(cl_sel_fun_str == "min_error") cl_sel_fun <- min_error_cl
  else if(cl_sel_fun_str == "centrality") cl_sel_fun <- central_cl
  
  candidate_cluster <- hyperparameters[hyperparameters$cluster == cl_sel_fun(hyperparameters),]
  cluster <- candidate_cluster$cluster[1]
  pred_results_candidate <- pred_results[pred_results$cluster == cluster,]
  
  #===== Intra-cluster recommendation of a combination of hyperparameters =====
  
  # Calculates the mode of each hyperparameter in the candidate cluster, and
  # selects the combination with the values of each mode (not yet used)
  var_mode_rec <- function(candidate_cluster){
    mode <- function(v){
      uniqv <- unique(v)
      mod <- which.max(tabulate(match(v, uniqv)))
      return(uniqv[mod])
    }
    mode_v <- apply(candidate_cluster,2,function(col) mode(col))
    return(mode_v)
  }
  # Selects the combination with the smallest prediction error in the candidate cluster
  min_error_rec <- function(candidate_cluster){
    min_error <- which.min(candidate_cluster[,eval_col])
    return(min_error)
  }
  # Selects the most central combination in the candidate cluster
  central_rec <-  function(candidate_cluster){
    n_cand_clstr <- norm(candidate_cluster)
    
    cand_clstr_g <- get_graph(n_cand_clstr,kernel_fun)
    
    centr_v <- which.max(centr_fun(cand_clstr_g))
    return(centr_v)
  }
  
  rec_fun_str <- match.arg(rec_fun)
  if(rec_fun_str == "mode") rec_fun <- var_mode_rec
  else if(rec_fun_str == "min_error") rec_fun <- min_error_rec
  else if(rec_fun_str == "centrality") rec_fun <- central_rec
  
  # Returns the recommended (optimized) hyperparameter combination
  
  recommendation_index = rec_fun(candidate_cluster)
  recommendation = pred_results_candidate[recommendation_index,]

  #print(recommendation)
  #recommendation <- as.data.frame(recommendation)
  return(as.list(recommendation))
  #return(list(rec=recommendation,graph=g,specc=specc_f,data=pred_results,sel_cluster=candidate_cluster))
}
