isoForest = function(pseq){
  
  # Isolation forest
  library(h2o)
  require(phyloseq)
  
  datos = as.data.frame(t(otu_table(pseq)@.Data))
  
  # Se inicializa el cluster H2O
  print('Initializing Isolation Forest...')
  h2o.init(ip = "localhost",
           # Todos los cores disponibles.
           nthreads = -1,
           # Máxima memoria disponible para el cluster.
           max_mem_size = "8g")
  
  h2o.removeAll()
  h2o.no_progress()
  
  # Carga de datos en el cluster H2O
  datos_h2o <- as.h2o(x = datos)
  
  # Modelo isolation forest
  print('Creating Isolation Forest Model...')
  isoforest <- h2o.isolationForest(
    model_id = "isoforest",
    training_frame = datos_h2o,
    x              = colnames(datos_h2o),
    max_depth      = 350, # Profundidad máxima de los árboles
    ntrees         = 100, # Número de los árboles
    sample_rate    = -1 # Ratio de observaciones empleadas en cada árbol
  )
  isoforest
  
  # Predicción
  predicciones_h2o <- h2o.predict(
    object  = isoforest,
    newdata = datos_h2o
  )
  predicciones <- as.data.frame(predicciones_h2o)
  head(predicciones)
  
  deciles <- quantile(x = predicciones$mean_length, probs = seq(0, 1, 0.1))
  deciles
  
  rmPats = which(predicciones$mean_length < deciles[2])
  Pats = setdiff(rownames(datos), rownames(datos)[rmPats])
  
  print(paste0('Removing ', length(rmPats), ' outliers patients'))
  print(rmPats)
  
  # res = prune_samples(Pats, pseq)
  
  print('Isolation forest analysis done!')
  
  return(Pats)
}

