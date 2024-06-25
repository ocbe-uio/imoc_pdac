utils_path = path.file("src", "utils", "utils.R")
source(utils_path)


datasets = c(
  "simulated_gm",
  "simulated_InterSIM",
  "simulated_netMUG",
  "nutrimouse_genotype",
  "nutrimouse_diet",
  "bbcsport",
  "buaa",
  "metabric",
  "digits",
  "bdgp",
  "tcga",
  "caltech101",
  "nuswide"
)

algorithms = c("IntNMF", "COCA", "jNMF", "NEMO")

results = read.csv(TIME_RESULTS_PATH, row.names = 1)


for (dataset_name in datasets) {
  names <- strsplit(dataset_name, "_")[[1]]
  if ("simulated" %in% names) {
    names <- paste(names, collapse = "_")
  }
  x_name <- names[1]
  y_name <- ifelse(length(names) > 1, names[2], "X0")
  
  
  data = load_dataset(dataset_name=x_name, return_y=T, shuffle=TRUE, seed= RANDOM_STATE)
  Xs <- data[["X"]]
  y <- data[["y"]]
  y <- y[,y_name]
  n_clusters <- length(unique(y))
  
  for (alg_name in algorithms) {
    if (!alg_name %in% row.names(results)) {
      results[alg_name,] <- 0
    }
    time_execution <- results[alg_name, dataset_name]
    if (time_execution > 0 | is.na(time_execution)) {
      next
    }
    cat(paste("\n", dataset_name, alg_name, Sys.time(), sep = "\t"),
        file=TIME_LOGS_PATH,append=TRUE)

    elapsed_time <- tryCatch({
      start_time <- Sys.time()
      if (alg_name == "IntNMF") {
        train_Xs = lapply(Xs, normalize)
        train_Xs = lapply(train_Xs, as.matrix)
        set.seed(RANDOM_STATE)
        clusters <- IntNMF::nmf.mnnals(dat=train_Xs, k=n_clusters, seed=RANDOM_STATE)$clusters
      } else if (alg_name == "COCA") {
        train_Xs = lapply(Xs, scale)
        set.seed(RANDOM_STATE)
        clusters <- coca::buildMOC(train_Xs, M = length(train_Xs), K = n_clusters)$moc
        clusters <- coca::coca(clusters, K = n_clusters)$clusterLabels
      } else if (alg_name == "jNMF") {
        train_Xs = lapply(Xs, normalize)
        train_Xs = lapply(train_Xs, as.matrix)
        set.seed(RANDOM_STATE)
        clusters <- nnTensor::jNMF(X= train_Xs, J = n_clusters)$W
        set.seed(RANDOM_STATE)
        clusters <- kmeans(clusters, n_clusters)$cluster
      } else if (alg_name == "NEMO") {
        train_Xs = lapply(Xs, t)
        train_Xs = lapply(train_Xs, scale)
        set.seed(RANDOM_STATE)
        clusters <- nemo.clustering(train_Xs, num.clusters= n_clusters)
      }
      stopifnot(length(y) == length(clusters))
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    }, error = function(e) {
      cat(paste("\n", dataset_name, alg_name, paste(class(e), e), Sys.time(), sep = "\t"),
          file=TIME_ERRORS_PATH,append=TRUE)
      return(NA)
    })

    results[alg_name, dataset_name] <- elapsed_time
    write.csv(results, file = TIME_RESULTS_PATH)
  }
}
cat("Completed successfully!")
cat("Completed successfully", file=TIME_LOGS_PATH,append=TRUE)
