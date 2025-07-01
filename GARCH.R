library(tseries)
library(rugarch)
library(PerformanceAnalytics)

#----------------------------------
# Funcion para ajustar modelo GARCH y calcular VaR y ES
#----------------------------------
calcular_garch_var_es <- function(data_returns, niveles_confianza = c(0.10, 0.05, 0.01)) {
  
  # Limpiar NA automáticamente
  data_returns <- na.omit(data_returns)

  resultados <- list()
  
  # Especificar el modelo GARCH(1,1)
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm"
  )
  
  # Ajustar el modelo
  fit <- ugarchfit(spec, data = data_returns)
  
  # Extraer las volatilidades condicionales
  sigma <- sigma(fit)
  media <- fitted(fit)
  
  # Calcular VaR y ES para cada nivel de confianza
  for (nivel in niveles_confianza) {
    var <- media + qnorm(nivel) * sigma
    es <- media - sigma * dnorm(qnorm(nivel)) / nivel
    
    resultados[[paste0("VaR_", nivel * 100, "%")]] <- var
    resultados[[paste0("ES_", nivel * 100, "%")]] <- es
  }
  
  resultados[["Fit"]] <- fit
  resultados[["Sigma"]] <- sigma
  resultados[["Media"]] <- media
  
  return(resultados)
}

#----------------------------------
# Funcion para graficar resultados
#----------------------------------

graficar_histograma_var_garch <- function(retornos, resultados_garch, niveles_confianza = c(0.10, 0.05, 0.01), nombre_moneda = "") {
  data_clean <- na.omit(retornos)
  
  df_plot <- data.frame(Retornos = data_clean)
  
  p <- ggplot(df_plot, aes(x = Retornos)) +
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "lightblue") +
    geom_density(color = "black", linewidth = 0.5,  alpha = 0.7) +
    labs(title = paste("Estimación del VaR y ES mediante GARCH -", nombre_moneda),
         x = "Retornos",
         y = "Densidad") +
    theme_minimal()
  
  colores <- c("orange", "red", "purple")
  line_data <- data.frame(Nivel = character(), Valor = numeric(), Tipo = character(), Color = character(), stringsAsFactors = FALSE)
  
  for (i in seq_along(niveles_confianza)) {
    nivel <- niveles_confianza[i]
    var_key <- paste0("VaR_", nivel * 100, "%")
    es_key <- paste0("ES_", nivel * 100, "%")
    
    var_val <- resultados_garch[[var_key]][length(resultados_garch[[var_key]])]
    es_val <- resultados_garch[[es_key]][length(resultados_garch[[es_key]])]
    
    line_data <- rbind(line_data, 
                       data.frame(Nivel = paste0("VaR ", (1 - nivel) * 100, "%"), Valor = var_val, Tipo = "VaR", Color = colores[i]),
                       data.frame(Nivel = paste0("ES ", (1 - nivel) * 100, "%"), Valor = es_val, Tipo = "ES", Color = colores[i]))
  }
  
  p <- p +
    geom_vline(data = line_data, aes(xintercept = Valor, color = Nivel, linetype = Tipo), linewidth = 1) +
    scale_color_manual(values = setNames(line_data$Color, line_data$Nivel)) +
    scale_linetype_manual(values = c("VaR" = "solid", "ES" = "dashed")) +
    guides(color = guide_legend(title = "Niveles de Confianza"))
  
  return(p)
}


#----------------------------------
# Funcion para comparar VaR y ES entre niveles de confianza
#----------------------------------

comparar_var_es<- function(resultados_garch, niveles_confianza = c(0.10, 0.05, 0.01)) {
  resumen <- data.frame(Nivel = numeric(), Tipo = character(), Valor_Final = numeric(), stringsAsFactors = FALSE)
  
  for (nivel in niveles_confianza) {
    var_key <- paste0("VaR_", nivel * 100, "%")
    es_key <- paste0("ES_", nivel * 100, "%")
    
    if (!is.null(resultados_garch[[var_key]]) && length(resultados_garch[[var_key]]) > 0) {
      var_val <- tail(resultados_garch[[var_key]], 1)
      resumen <- rbind(resumen, data.frame(Nivel = (1 - nivel) * 100, Tipo = "VaR", Valor_Final = var_val))
    }
    
    if (!is.null(resultados_garch[[es_key]]) && length(resultados_garch[[es_key]]) > 0) {
      es_val <- tail(resultados_garch[[es_key]], 1)
      resumen <- rbind(resumen, data.frame(Nivel = (1 - nivel) * 100, Tipo = "ES", Valor_Final = es_val))
    }
  }
  
  return(resumen)
}


###########

calcular_garch_var_es_COLON <- function(data_returns, niveles_confianza = c(0.10, 0.05, 0.01)) {
  
  # Limpieza de NA e infinitos
  data_returns <- na.omit(data_returns)
  data_returns <- data_returns[is.finite(data_returns)]
  
  # Eliminar outliers extremos mayores a 5 desviaciones estándar
  media_ret <- mean(data_returns)
  sd_ret <- sd(data_returns)
  data_returns <- data_returns[abs(data_returns - media_ret) <= 5 * sd_ret]
  
  # Verificar que queden suficientes datos
  if (length(data_returns) < 50) {
    warning("No hay suficientes datos después de limpiar la serie. Se necesitan al menos 50 observaciones.")
    return(NULL)
  }
  
  resultados <- list()
  
  # Especificar modelo GARCH(1,1) con distribución t-Student
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"   # Distribución t-Student
  )
  
  # Ajuste con control de errores y solver robusto
  fit <- tryCatch({
    ugarchfit(spec, data = data_returns, solver = "hybrid", solver.control = list(trace = 0))
  }, error = function(e) {
    message("El modelo no convergió: ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) {
    warning("El modelo GARCH no pudo ajustarse a la serie proporcionada.")
    return(NULL)
  }
  
  sigma <- sigma(fit)
  media <- fitted(fit)
  
  # Calcular VaR y ES para cada nivel de confianza
  for (nivel in niveles_confianza) {
    var <- media + qnorm(nivel) * sigma
    es <- media - sigma * dnorm(qnorm(nivel)) / nivel
    
    resultados[[paste0("VaR_", nivel * 100, "%")]] <- var
    resultados[[paste0("ES_", nivel * 100, "%")]] <- es
  }
  
  resultados[["Fit"]] <- fit
  resultados[["Sigma"]] <- sigma
  resultados[["Media"]] <- media
  
  return(resultados)
}



