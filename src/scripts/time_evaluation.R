utils_path = file.path("src", "utils", "utils.R")
source(utils_path)


datasets = rownames(read.csv(DATASET_TABLE_PATH, row.names = 1))

if ("nutrimouse" %in% datasets) {
  pos <- which(datasets == "nutrimouse")
  datasets <- c(datasets[1:(pos-1)], "nutrimouse_genotype", "nutrimouse_diet", datasets[(pos+1):length(datasets)])
}

algorithms = c("IntNMF", "COCA", "jNMF", "NEMO")

results = read.csv(TIME_RESULTS_PATH, stringsAsFactors= FALSE)
convert_to_logical <- function(column) {
  if (is.character(column) && all(column %in% c("True", "False"))) {
    return(as.logical(ifelse(column == "True", TRUE, FALSE)))
  } else {
    return(column)
  }
}
results[] <- lapply(results, convert_to_logical)

for (dataset_name in datasets) {
  names <- strsplit(dataset_name, "_")[[1]]
  if ("simulated" %in% names) {
    names <- paste(names, collapse = "_")
  }
  x_name <- names[1]
  y_name <- ifelse(length(names) > 1, names[2], "X0")

  unfinished_results <- results[results[["dataset"]] == dataset_name,]
  unfinished_results <- unfinished_results[!unfinished_results[["finished"]],]

  if(nrow(unfinished_results) == 0){
        next
    }
  
  data = load_dataset(dataset_name=x_name, return_y=T, shuffle=TRUE, seed= RANDOM_STATE)
  Xs <- data[["X"]]
  y <- data[["y"]]
  y <- y[,y_name]
  n_clusters <- length(unique(y))

  for (i in rownames(unfinished_results)) {
    row <- results[i,]
    if(row[["finished"]]){
        next
    }
    alg_name <- row[["algorithm"]]
    if(!(alg_name %in% algorithms)){
        next
    }
    cat(paste("\n", dataset_name, alg_name, Sys.time(), sep = "\t"),
        file=TIME_LOGS_PATH,append=TRUE)

    start_time <- Sys.time()
    clustering <- tryCatch({
      if (alg_name == "IntNMF") {
        train_Xs = lapply(Xs, function(x){x[,apply(x, 2, var, na.rm=TRUE) != 0]})
        train_Xs = lapply(train_Xs, normalize)
        train_Xs = lapply(train_Xs, as.matrix)
        set.seed(RANDOM_STATE)
        clusters <- IntNMF::nmf.mnnals(dat=train_Xs, k=n_clusters, seed=RANDOM_STATE)$clusters
      } else if (alg_name == "COCA") {
        train_Xs = lapply(Xs, function(x){x[,apply(x, 2, var, na.rm=TRUE) != 0]})
        train_Xs = lapply(train_Xs, scale)
        set.seed(RANDOM_STATE)
        clusters <- coca::buildMOC(train_Xs, M = length(train_Xs), K = n_clusters)$moc
        clusters <- coca::coca(clusters, K = n_clusters)$clusterLabels
      } else if (alg_name == "jNMF") {
        train_Xs = lapply(Xs, function(x){x[,apply(x, 2, var, na.rm=TRUE) != 0]})
        train_Xs = lapply(train_Xs, normalize)
        train_Xs = lapply(train_Xs, as.matrix)
        set.seed(RANDOM_STATE)
        clusters <- nnTensor::jNMF(X= train_Xs, J = n_clusters)$W
        set.seed(RANDOM_STATE)
        clusters <- kmeans(clusters, n_clusters)$cluster
      } else if (alg_name == "NEMO") {
        train_Xs = lapply(Xs, function(x){x[,apply(x, 2, var, na.rm=TRUE) != 0]})
        train_Xs = lapply(train_Xs, scale)
        train_Xs = lapply(train_Xs, t)
        set.seed(RANDOM_STATE)
        clusters <- nemo.clustering(train_Xs, num.clusters= n_clusters)
      }
      stopifnot(length(y) == length(clusters))
      err <- "{}"
      completed <- TRUE
      clustering <- list(err, completed)
    }, error = function(e) {
      cat(paste("\n", dataset_name, alg_name, e, Sys.time(), sep = "\t"),
          file=TIME_ERRORS_PATH,append=TRUE)
      err <- e
      completed <- FALSE
      return(list(err, completed))
    })
    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    results[i, "finished"] = TRUE
    results[i, "completed"] = clustering[[2]]
    results[i, "time"] = elapsed_time
    results[i, "comments"] = as.character(clustering[[1]])

    write.csv(results, file = TIME_RESULTS_PATH, row.names = F)
  }
}
cat("Completed successfully!")
cat("Completed successfully", file=TIME_LOGS_PATH,append=TRUE)
