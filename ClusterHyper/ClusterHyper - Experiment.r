library(factoextra)
library(dplyr)
library(reshape2)
library(grid)
library(gridExtra)
library(scales)

#Kernel options:
# rbfdot Radial Basis kernel function "Gaussian"
# polydot Polynomial kernel function
# vanilladot Linear kernel function
# tanhdot Hyperbolic tangent kernel function
# laplacedot Laplacian kernel function
# besseldot Bessel kernel function
# anovadot ANOVA RBF kernel function
# splinedot Spline kernel
# stringdot String kernel
kernels <- c("rbfdot")
#Number of clusters options: 1:10
centers <- 1:10
#Number of algorithm iterations: {100;200;500;1000;5000;10000}
iterations <- c(100,200,500,1000,5000,10000)
#Normalization options for the vector variables: {minmax;zscore}
norms <- c("zscore")
#Candidate cluster selection method options:
# minimum mean of prediction errors
# centrality between the most central vertices of each cluster
cl_sel <- c("min_error","centrality")
#Prediction error evaluation measure:
# mean ('mse_test')
# mean+var ('mse_sd_sum_test' ='mse_test'+'mse_test_sd')
evals <- c("sum_smape_vc")
#Intra-cluster recommendation method options:
# variable mode
# minimum prediction error (minimum 'mse_test')
# centrality (intra-cluster)
cl_rec <- c("min_error","centrality")

#Creating the experiment options table
options <- function(kernels, norms, cl_sel, evals, cl_rec) {
  opts <- as.data.frame(expand.grid(kernels, norms, cl_sel, evals, cl_rec))
  colnames(opts) <- c("kernels", "norms", "cl_sel", "evals", "cl_rec")
  rownames(opts) <- LETTERS[1:nrow(opts)]
  return(opts)
}

#Creating the experiment options table
options_kmeans <- function(norms, cl_sel, evals, cl_rec) {
  opts <- as.data.frame(expand.grid(norms, cl_sel, evals, cl_rec))
  colnames(opts) <- c("norms", "cl_sel", "evals", "cl_rec")
  rownames(opts) <- LETTERS[1:nrow(opts)]
  return(opts)
}

opts <- options(kernels, norms, cl_sel, evals, cl_rec)

# Fertilizer list
fertilizers <- c("Brazil_K2O",
                 "Brazil_P2O5",  
                 "Brazil_N",  
                 "Canada_K2O", 
                 "Canada_P2O5", 
                 "Canada_N",
                 "China_K2O", 
                 "China_P2O5", 
                 "China_N",
                 "France_K2O", 
                 "France_P2O5", 
                 "France_N", 
                 "India_K2O",
                 "India_P2O5", 
                 "India_N", 
                 "Indonesia_K2O", 
                 "Indonesia_P2O5", 
                 "Indonesia_N", 
                 "Pakistan_K2O", 
                 "Pakistan_P2O5", 
                 "Pakistan_N",
                 "Russia_K2O", 
                 "Russia_P2O5", 
                 "Russia_N", 
                 "Turkey_K2O", 
                 "Turkey_P2O5", 
                 "Turkey_N",
                 "USA_K2O", 
                 "USA_P2O5", 
                 "USA_N")

# To select all fertilizer files in the defined range
# for the number of cross-validated partitions

minParticoesVC <- 2
maxParticoesVC <- 4
fert <- vector()
for(i in 1:length(fertilizers)){
  for(j in minParticoesVC:maxParticoesVC){
    str <- sprintf("%s_%d", fertilizers[i], j)
    fert <- c(fert, str)
  }
}

fertilizers <- fert
rm(fert)
rm(minParticoesVC)
rm(maxParticoesVC)

################################################################################

# Spectral clustering optimization

# Number of fertilizers to be optimized
val_size <- length(fertilizers)

# Which hyperparameters will be optimized
hyperparameter_names <- c("input_size", "size", "decay")

fertilizers_analysis_list <- list()

for (i in 1:val_size){
  
  # Table to store optimization output recommendations
  fertilizers_analysis <- data.frame(
                                      fertilizer = character(),
                                      config = character(),
                                      ml_method = character(),
                                      preprocessing = character(),
                                      sw_size = numeric(),
                                      input_size = numeric(),
                                      size = numeric(),
                                      decay = numeric(),
                                      offset = numeric(),
                                      steps_ahead = numeric(),
                                      pred_size = numeric(),
                                      smape_train = numeric(),
                                      smape_test = numeric(),
                                      sum_smape_vc = numeric(),
                                      grupo = numeric()
                                    )
  
  # Creation of the column "smape_sum_vc" which will be the sum of the forecast errors
  # in each partition of the cross validation
  
  fertilizer = c("fertilizer" = fertilizers[i])
  path <- sprintf("agrupamento_espectral/data/%s.RData", fertilizer)
  load(file = path)
  s <- grep(pattern = "smape_ho", x = names(results), value = TRUE)
  results$sum_smape_vc <- rowSums(results[, s])
  
  # Hyperparameters that were not chosen for optimization must be
  # fixed
  
  results <- filter(results, method == "mlp")
  results <- filter(results, preprocessing == "global")
  results <- filter(results, sw_size == 10)
  results <- filter(results, offset == 0)
  results <- filter(results, pred_size == 8)
  results <- filter(results, steps_ahead == 1)

  # For each experiment setup
  
  for(j in 1:nrow(opts)){
    
    # Extracts the configuration identifier and all its components
    
    conf = c("config" = rownames(opts)[j])
    eval_col = as.character(opts$evals[j])
    norm_fun = as.character(opts$norms[j])
    kernel = as.character(opts$kernels[j])
    cl_sel_fun = as.character(opts$cl_sel[j])
    rec_fun = as.character(opts$cl_rec[j])
    
    # Pass settings to hyperparameter optimization function
    
    specc_fertilizer <- specc_rec(results,
                                  hyperparameter_names,
                                  eval_col=eval_col,
                                  centers=3, 
                                  kernel=kernel,
                                  norm_fun = norm_fun,
                                  cl_sel_fun = cl_sel_fun,
                                  rec_fun = rec_fun,
                                  centr_fun = function(g) igraph::eigen_centrality(g)$`vector`)
    
    # Transforming the output of "specc_rec" for easier parsing
    
    specc_fertilizer <- as.data.frame(specc_fertilizer)
    subtab1 <- specc_fertilizer[, c(3:11)]
    subtab2 <- select(specc_fertilizer, smape_train, smape_test, sum_smape_vc, cluster)
    subtab1 <- cbind(subtab1, subtab2)
    specc_fertilizer <- subtab1
    fertilizers_analysis <- rbind(fertilizers_analysis, 
                                  as.list(c(fertilizer, conf, specc_fertilizer)))
    
    # The set of recommendations obtained for a fertilizer is stored
    # in a list
    fertilizers_analysis_list[[fertilizers[i]]] <- fertilizers_analysis
  }
}

# Saves spectral clustering optimization results to a file
path <- sprintf("agrupamento_espectral/fertilizers_analysis_specc.RData")
save(fertilizers_analysis_list, file = path)

################################################################################

# K-means optimization
  
# Number of fertilizers to be optimized
val_size <- length(fertilizers)

# Which hyperparameters will be optimized
hyperparameter_names <- c("input_size", "size", "decay")

fertilizers_analysis_list <- list()
  
for (i in 1:val_size){
  
  # Table to store optimization output recommendations
  fertilizers_analysis <- data.frame(
                                      fertilizer = character(),
                                      config = character(),
                                      ml_method = character(),
                                      preprocessing = character(),
                                      sw_size = numeric(),
                                      input_size = numeric(),
                                      size = numeric(),
                                      decay = numeric(),
                                      offset = numeric(),
                                      steps_ahead = numeric(),
                                      pred_size = numeric(),
                                      smape_train = numeric(),
                                      smape_test = numeric(),
                                      sum_smape_vc = numeric(),
                                      grupo = numeric()
                                    )
    
  # Creation of the column "smape_sum_vc" which will be the sum of the forecast errors
  # in each partition of the cross validation
  fertilizer = c("fertilizer" = fertilizers[i])
  path <- sprintf("agrupamento_espectral/data/%s.RData", fertilizers[i])
  load(file = path)
  s <- grep(pattern = "smape_ho", x = names(results), value = TRUE)
  results$sum_smape_vc <- rowSums(results[, s])
  
  # Hyperparameters that were not chosen for optimization must be
  # fixed
  
  results <- filter(results, method == "mlp")
  results <- filter(results, preprocessing == "global")
  results <- filter(results, sw_size == 10)
  results <- filter(results, offset == 0)
  results <- filter(results, pred_size == 8)
  results <- filter(results, steps_ahead == 1)
  
  # For each experiment setup
    
  for(j in 1:nrow(opts)){
    
    # Extracts the configuration identifier and all its components
    
    conf = c("config" = rownames(opts)[j])
    eval_col = as.character(opts$evals[j])
    norm_fun = as.character(opts$norms[j])
    cl_sel_fun = as.character(opts$cl_sel[j])
    rec_fun = as.character(opts$cl_rec[j])
    
    # Pass settings to hyperparameter optimization function
    
    kmeans_fertilizer <- kmeans_rec(results,
                                    hyperparameter_names,
                                    eval_col=eval_col,
                                    centers=3, 
                                    norm_fun = norm_fun,
                                    cl_sel_fun = cl_sel_fun,
                                    rec_fun = rec_fun,
                                    centr_fun = function(g) igraph::eigen_centrality(g)$`vector`)
    
    # Transforming the output of "specc_rec" for easier parsing
    
    kmeans_fertilizer <- as.data.frame(kmeans_fertilizer)
    subtab1 <- kmeans_fertilizer[, c(3:11)]
    subtab2 <- select(kmeans_fertilizer, smape_train, smape_test, sum_smape_vc, cluster)
    subtab1 <- cbind(subtab1, subtab2)
    kmeans_fertilizer <- subtab1
    fertilizers_analysis <- rbind(fertilizers_analysis, 
                                  as.list(c(fertilizer, conf, kmeans_fertilizer)))
    
    # The set of recommendations obtained for a fertilizer is stored in
    # a list
    fertilizers_analysis_list[[fertilizers[i]]] <- fertilizers_analysis
  }
}

# Save k-means optimization results to a file
path <- sprintf("agrupamento_espectral/fertilizers_analysis_kmeans.RData")
save(fertilizers_analysis_list, file = path)

################################################################################

# Compara????o entre as configura????es. Para cada fertilizante, avalia-se qual foi
# a configura????o usada para encontrar a recomenda????o com menor erro de previs??o
# no treinamento. Conta-se os resultados de cada fertilizante

clust_method = c("specc", "kmeans")

list_best_config = list()
list_clust_stats <- list()

for(i in 1:length(clust_method)){
  
  # Carrega a lista de conjuntos de recomenda????es de cada fertilizante
  path <- sprintf("agrupamento_espectral/fertilizers_analysis_%s.RData", clust_method[i])
  load(file = path)
  
  # Fun????o para descobrir qual configura????o tem menor erro de previs??o de 
  # treinamento
  
  getMinTrain <- function(data){
    index <- which.min(data$smape_train)
    x <- (data[index,])$config
    return(x)
  }
  
  # Caso alguma configura????o tenha sido a melhor em nenhum fertilizante, essa
  # configura????o estar?? ausente na tabela de compara????o entre configura????es.
  # Para resolver o problema, a fun????o abaixo insere a configura????o na tabela
  # com valor 0.
  
  completaTabelaConfiguracoes <- function(data, vConfigs){
    configsAusentes <- setdiff(vConfigs, data$config)
    for(i in 1:length(configsAusentes)){
      row <- c(config=configsAusentes[i], counts=0)
      data <- rbind(data, row)
    }
    return(data)  
  }
  
  # Descobre a melhor configura????o de cada recomenda????o
  l <- sapply(fertilizers_analysis_list, getMinTrain)
  # Conta quantas vezes cada recomenda????o foi a melhor
  l <- as.data.frame(table(l), stringsAsFactors = FALSE)
  colnames(l) <- c("config", "counts")
  m <- completaTabelaConfiguracoes(l, rownames(opts))
  m$counts <- as.numeric(m$counts)
  # Transforma a coluna de quantidades em porcentagens
  m <- mutate(m, perc = counts/sum(counts))
  
  end <- sprintf("comparacaoConfigs_%s.csv", i)
  write.csv(m, end)
  
  # Imprime uma figura comparando as configura????es
  
  path <- sprintf("agrupamento_espectral/images/best_train_%s.png", clust_method[i])
  png(file = path)
  g <- ggplot(m, aes(y=perc, x=config, fill=config)) +
    labs(x = "", y = "Porcentagem", fill = "Configura????o") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_y_continuous(labels=percent) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    geom_bar(position="dodge", color = "black", stat = "identity")
  print(g)
  dev.off()
  
  # Armazena as melhores configura????es de cada m??todo de otimiza????o de
  # hiperpar??metros
  
  list_best_config[[clust_method[i]]] <- l[which.max(l$counts),]$config
  list_fert_stats <- list()
  
  # Para cada fertilizante
  
  for(k in 1:val_size){
    
    # Seleciona a melhor recomenda????o de hiperpar??metros para o fertilizante
    
    fertilizer = c("fertilizer" = fertilizers[k])
    best_config = list_best_config[[clust_method[i]]]
    path <- sprintf("agrupamento_espectral/data/%s.RData", fertilizers[k])
    load(file = path)
    candidate_vector <- fertilizers_analysis_list[[fertilizer]]
    candidate_vector <- candidate_vector[candidate_vector$config == best_config,]
    # Seleciona a combina????o de hiperpar??metros com menor erro de previs??o
    # global no treinamento (Busca em Grade)
    min_train_vector <- results[which.min(results$smape_train),]
    # Seleciona a combina????o de hiperpar??metros com menor erro de previs??o
    # global no teste
    min_test_vector <- results[which.min(results$smape_test),]
    
    # Organiza as recomenda????es em uma tabela, e guarda em uma lista
    list_fert_stats[[fertilizer]] <- data.frame(candidate_test = candidate_vector$smape_test, 
                                                min_train = min_train_vector$smape_test, 
                                                min_test = min_test_vector$smape_test)
  }
  # Salva os resultados do agrupamento espectral e do k-means em uma lista
  list_clust_stats[[clust_method[i]]] <- list_fert_stats
}

# Salva as combina????es de hiperpar??metros recomendadas em cada fertilizante
save(list_clust_stats, file = "agrupamento_espectral/cluster_statistics.RData")

################################################################################

# Pr??xima etapa: descobrir qual ?? o melhor m??todo de otimiza????o de
# hiperpar??metros (agrupamento espectral, k-means, Busca em Grade). As
# recomenda????es obtidas por cada m??todo ser??o comparadas com base no valor de 
# erro de previs??o na etapa de testes

# vetor que guarda os erros de previs??o de cada recomenda????o feita pelo
# agrupamento espectral
specc_vector = c()
# vetor que guarda os erros de previs??o de cada recomenda????o feita pelo k-means
kmeans_vector = c()
# vetor que guarda os erros de previs??o de cada recomenda????o feita pela Busca
# em grade
gridsearch_vector = c()
# vetor que guarda os menores erros de previs??o globais em cada fertilizante
min_test_vector = c()

# Para cada m??todo de otimiza????o (espectral, k-means)
for(i in 1:length(clust_method)){
  
  # Carrega as combina????es de hiperpar??metros recomendadas em cada fertilizante
  load("agrupamento_espectral/cluster_statistics.RData")
  
  fertilizer_list <- list_clust_stats[[clust_method[i]]]
  # Para cada fertilizante
  for(j in 1:length(fertilizer_list)){
    # Extrai os valores de erro de previs??o de teste do agrupamento espectral,
    # Busca em Grade e os erros m??nimos globais
    if(clust_method[i] == "specc"){
      specc_vector = c(specc_vector, fertilizer_list[[j]][['candidate_test']])
      gridsearch_vector = c(gridsearch_vector, fertilizer_list[[j]][['min_train']])
      min_test_vector = c(min_test_vector, fertilizer_list[[j]][['min_test']])
    }
    # Extrai os valores de erro de previs??o de teste do k-means
    if(clust_method[i] == "kmeans"){
      kmeans_vector = c(kmeans_vector, fertilizer_list[[j]][['candidate_test']])
    }
  }
}

# Organiza cada distribui????o de erros de previs??o em uma tabela

tabelaDistribuicoes <- data.frame(v_specc = specc_vector, 
                                  v_kmeans = kmeans_vector, 
                                  v_gridsearch = gridsearch_vector,
                                  v_min_test = min_test_vector)
rownames(tabelaDistribuicoes) <- fertilizers[1:length(rownames(tabelaDistribuicoes))]

gridSearchMelhor <- filter(tabelaDistribuicoes, v_gridsearch < v_specc |
                             v_gridsearch < v_kmeans)

# Extrai as linhas em que o erro de previs??o do espectral ou do k-means s??o 
# menores que o da Busca em Grade
result <- filter(tabelaDistribuicoes, (v_specc < v_gridsearch) & (v_kmeans < v_gridsearch))
# Conta em quantos fertilizantes o espectral atingiu o menor erro de previs??o
n_specc <- nrow(filter(result, v_specc < v_kmeans))
# Conta em quantos fertilizantes o k-means atingiu o menor erro de previs??o
n_kmeans <- nrow(filter(result, v_specc > v_kmeans))
# Conta em quantos fertilizantes o espectral e o k-means atingiram o menor erro 
# de previs??o
n_equal <- nrow(filter(result, v_specc == v_kmeans))
# Conta em quantos fertilizantes a Busca em Grade atingiu o menor erro de 
# previs??o
n_gridsearch <- nrow(filter(gridSearchMelhor, v_gridsearch < v_specc & v_gridsearch < v_kmeans))
n_gridsearch2 <- nrow(filter(gridSearchMelhor, v_gridsearch > v_specc))
n_gridsearch3 <- nrow(filter(gridSearchMelhor, v_gridsearch > v_kmeans))

# Compara em quantos fertilizantes cada m??todo de otimiza????o de hiperpar??metros
# obteve o menor erro de previs??o
compTabela = data.frame(x=c("Spectral", "K-means", "Spectral and K-means", "Grid Search"),
                        y=c(n_specc+n_gridsearch2, n_kmeans+n_gridsearch3, n_equal, n_gridsearch))
compTabela <- mutate(compTabela, perc = y/sum(y))
png(file = "agrupamento_espectral/images/compMetodosOtimizacao.png")
p <- ggplot(data = compTabela, aes(x = x, y = perc, fill = x)) +
  labs(x = "", y = "Percentage", fill = "Optimization Method") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_continuous(labels=percent) +
  scale_fill_manual(values=c("green", "red", "yellow", "blue")) +
  theme_bw() +
  geom_bar(position="dodge", color = "black", stat = "identity")
print(p)
dev.off()
write.csv(compTabela, "compTabela.csv")

# Extrai o fertilizante no qual o agrupamento espectral obteve os melhores
# resultados (comparando com os outros m??todos)
result <- result[which.min(result$v_specc),]
result <- melt(result)
# Apresenta uma figura comparando os valores de erro de previs??o para cada m??todo
# no fertilizante selecionado
png(file = "agrupamento_espectral/images/resultFertilizer.png")
p <- ggplot(data = result, aes(x = variable, y = value, fill = variable)) +
  labs(x = "", y = "Erro de previs??o no teste", fill = "Vetor de hiperpar??metros") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_fill_manual(values=c("red", "blue", "green", "purple")) +
  theme_bw() +
  geom_bar(position="dodge", color = "black", stat = "identity")
print(p)
dev.off()

# Apresenta as distribui????es de smape_test para espectral, k-means, Busca em 
# Grade e m??nimos globais

tabelaDistribuicoes_melt <- melt(tabelaDistribuicoes)
png(file = "agrupamento_espectral/images/distTestError.png")
p <- ggplot(data = tabelaDistribuicoes_melt, aes(x = variable, y = value)) +
  geom_boxplot() + 
  theme_minimal() + 
  labs(x = "", y = "smape_test")
print(p)
dev.off()

#Boxplot dos 60% dos casos em que o Grid Search padr??o foi o melhor m??todo
# de otimiza????o

gridSearchMelhor2 <- filter(tabelaDistribuicoes, v_gridsearch < v_specc &
                              v_gridsearch < v_kmeans)
gridSearchMelhor2 <- gridSearchMelhor2[1:3]
p <- ggplot(data = melt(gridSearchMelhor2), aes(x = variable, y = value)) +
  geom_boxplot() + 
  theme_minimal() + 
  labs(x = "", y = "smape_test")
print(p)

# Prepara os testes de hip??tese de Wilcoxon, selecionando cada par de
# distribui????es de erros de previs??o

tabDist1 <- melt(tabelaDistribuicoes[, c('v_specc', 'v_kmeans')])
tabDist2 <- melt(tabelaDistribuicoes[, c('v_specc', 'v_gridsearch')])
tabDist3 <- melt(tabelaDistribuicoes[, c('v_kmeans', 'v_gridsearch')])

# Faz os testes de hip??tese para cada par de distribui????es, e extrai os 
# p-valores resultantes

p1_1 <- (wilcox.test(value ~ variable, data = tabDist1, alternative='two.sided', exact=FALSE))[['p.value']]
p1_2 <- (wilcox.test(value ~ variable, data = tabDist2, alternative='two.sided', exact=FALSE))[['p.value']]
p1_3 <- (wilcox.test(value ~ variable, data = tabDist3, alternative='two.sided', exact=FALSE))[['p.value']]
p2_1 <- (wilcox.test(value ~ variable, data = tabDist1, alternative='less', exact=FALSE))[['p.value']]
p2_2 <- (wilcox.test(value ~ variable, data = tabDist2, alternative='less', exact=FALSE))[['p.value']]
p2_3 <- (wilcox.test(value ~ variable, data = tabDist3, alternative='less', exact=FALSE))[['p.value']]

# Organiza os p-valores em duas tabelas, uma para a hip??tese alternativa != 0 e
# outra para a hip??tese alternativa < 0
pValoresTeste1 <- data.frame(v_specc=c(1, p1_1, p1_2), 
                             v_kmeans=c(p1_1, 1, p1_3),
                             v_gridsearch=c(p1_2, p1_3, 1))
pValoresTeste1 <- as.data.frame(sapply(pValoresTeste1, round, 7))
rownames(pValoresTeste1) <- colnames(pValoresTeste1)

pValoresTeste2 <- data.frame(v_specc=c(1, p2_1, p2_2), 
                             v_kmeans=c(p2_1, 1, p2_3),
                             v_gridsearch=c(p2_2, p2_3, 1))
pValoresTeste2 <- as.data.frame(sapply(pValoresTeste2, round, 7))
rownames(pValoresTeste2) <- colnames(pValoresTeste2)

# Gera os vetores de diferen??a, a partir das distribui????es de erro de previs??o
# de cada m??todo de otimiza????o de hiperpar??metros
geraVetoresDiferencas <- function(v){
  diff_vector <- v - min_test_vector
  return(diff_vector)
}

tabelaDiferencas <- as.data.frame(sapply(tabelaDistribuicoes[1:3], geraVetoresDiferencas))
colnames(tabelaDiferencas) <- c("v_diff_specc", 
                                "v_diff_kmeans", 
                                "v_diff_gridsearch")

tabelaDiferencas_melt <- melt(tabelaDiferencas)

# Apresenta os vetores de diferen??a como boxplots
png(file = "agrupamento_espectral/images/distDiffs.png")
p <- ggplot(data = tabelaDiferencas_melt, aes(x = variable, y = value)) +
  geom_boxplot() + 
  theme_minimal() + 
  labs(x = "", y = "Erro de previs??o no teste")
print(p)
dev.off()

# Prepara os testes de hip??tese de Wilcoxon, selecionando cada par de
# distribui????es de erros de previs??o

tabDiff1 <- melt(tabelaDiferencas[, c('v_diff_specc', 'v_diff_kmeans')])
tabDiff2 <- melt(tabelaDiferencas[, c('v_diff_specc', 'v_diff_gridsearch')])
tabDiff3 <- melt(tabelaDiferencas[, c('v_diff_kmeans', 'v_diff_gridsearch')])

# Faz os testes de hip??tese para cada par de distribui????es, e extrai os 
# p-valores resultantes

p1_1 <- (wilcox.test(value ~ variable, data = tabDiff1, alternative='two.sided', exact=FALSE))[['p.value']]
p1_2 <- (wilcox.test(value ~ variable, data = tabDiff2, alternative='two.sided', exact=FALSE))[['p.value']]
p1_3 <- (wilcox.test(value ~ variable, data = tabDiff3, alternative='two.sided', exact=FALSE))[['p.value']]
p2_1 <- (wilcox.test(value ~ variable, data = tabDiff1, alternative='less', exact=FALSE))[['p.value']]
p2_2 <- (wilcox.test(value ~ variable, data = tabDiff2, alternative='less', exact=FALSE))[['p.value']]
p2_3 <- (wilcox.test(value ~ variable, data = tabDiff3, alternative='less', exact=FALSE))[['p.value']]

# Organiza os p-valores em duas tabelas, uma para a hip??tese alternativa != 0 e
# outra para a hip??tese alternativa < 0
pValoresDiff1 <- data.frame(v_diff_specc=c(1, p1_1, p1_2), 
                             v_diff_kmeans=c(p1_1, 1, p1_3),
                             v_diff_gridsearch=c(p1_2, p1_3, 1))
pValoresDiff1 <- as.data.frame(sapply(pValoresDiff1, round, 7))
rownames(pValoresDiff1) <- colnames(pValoresDiff1)

pValoresDiff2 <- data.frame(v_diff_specc=c(1, p2_1, p2_2), 
                             v_diff_kmeans=c(p2_1, 1, p2_3),
                             v_diff_gridsearch=c(p2_2, p2_3, 1))
pValoresDiff2 <- as.data.frame(sapply(pValoresDiff2, round, 7))
rownames(pValoresDiff2) <- colnames(pValoresDiff2)

# Fun????o para exibir uma tabela como figura
geraGrid <- function(x){
  myTable <- tableGrob(
    x, 
    theme = ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)),
                           rowhead=list(fg_params=list(hjust=0, x=0)))
  )
  grid.draw(myTable)
}

# Exibe todas as tabelas de p-valores
geraGrid(pValoresTeste1)
geraGrid(pValoresTeste2)
geraGrid(pValoresDiff1)
geraGrid(pValoresDiff2)

################################################################################
# Gera os vetores de diferen??a, a partir das distribui????es de erro de previs??o
# de cada m??todo de otimiza????o de hiperpar??metros
geraVetoresDiferencasGrid <- function(v, v_gridsearch){
  p <- v - v_gridsearch
  return(p)
}

tabelaMelhores <- pmin(tabelaDistribuicoes[1], tabelaDistribuicoes[2])
tabelaDiferencas2 <- as.data.frame(sapply(tabelaDistribuicoes, geraVetoresDiferencasGrid, tabelaDistribuicoes[3]))
colnames(tabelaDiferencas2) <- c("v_diff")

p <- ggplot(data = tabelaDiferencas2, aes(x = variable, y = value)) +
  geom_boxplot() + 
  theme_minimal() + 
  labs(x = "", y = "smape_test")
print(p)
