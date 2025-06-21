calcular_var_es_copula <- function(datos, niveles_confianza = c(0.95, 0.99), 
                                   especificacion_garch = "sGARCH", n_sim = 10000) {
  
  # Esta línea asegura que los datos no posean valores que impidan la ejecución de la función optim y que no se 
  # incluyan valores atípicos.
  datos <- datos[is.finite(datos$Returns), ] 
  retornos <- datos$Returns
  retornos_limpios <- retornos[abs(retornos - mean(retornos)) <= 5 * sd(retornos)]

  retornos <- retornos_limpios

  # Se prueban diferentes especificaciones GARCH si la primera falla, en particular: 
  # GARCH(1,1) simple con distribución normal,
  # GARCH(1,1) con distribución t,
  # GARCH simple sin términos ARMA
  especificaciones_garch <- list(

    list(varianza = list(model = especificacion_garch, garchOrder = c(1, 1)),
         media = list(armaOrder = c(0, 0), include.mean = TRUE),
         distribucion = "norm"),
    
    list(varianza = list(model = especificacion_garch, garchOrder = c(1, 1)),
         media = list(armaOrder = c(0, 0), include.mean = TRUE),
         distribucion = "std"),
    
    list(varianza = list(model = especificacion_garch, garchOrder = c(1, 1)),
         media = list(armaOrder = c(0, 0), include.mean = FALSE),
         distribucion = "norm")
  )
  
  ajuste_garch <- NULL
  
  for (i in seq_along(especificaciones_garch)) {
    tryCatch({
      espec <- especificaciones_garch[[i]]
      objeto_spec_garch <- ugarchspec(
        variance.model = espec$varianza,
        mean.model = espec$media,
        distribution.model = espec$distribucion
      )
      
      ajuste_garch <- ugarchfit(spec = objeto_spec_garch, data = retornos, 
                                solver = "hybrid", fit.control = list(scale = 1))
      
      if (convergence(ajuste_garch) == 0) {
        cat("Modelo GARCH ajustado exitosamente con especificación", i, "\n")
        break
      }
    }, error = function(e) {
      cat("Especificación GARCH", i, "falló:", e$message, "\n")
    })
  }
  
  if (is.null(ajuste_garch) || convergence(ajuste_garch) != 0) {
    cat("Ajuste GARCH falló. Usando enfoque empírico...\n")
    residuos_estandarizados <- scale(retornos)[,1]  # Retornos estandarizados
    usar_garch <- FALSE
  } else {
    residuos_estandarizados <- residuals(ajuste_garch, standardize = TRUE)
    usar_garch <- TRUE
  }
  
  # Se precisa transformar los valores a márgenes uniformes usando la CDF empírica. Para ello se 
  # Crea una configuración bivariada para ajustar la cópula utilizando los retornos rezagados para 
  # capturar la dependencia temporal.
  
  residuos_estandarizados <- as.numeric(residuos_estandarizados)
  residuos_estandarizados <- residuos_estandarizados[is.finite(residuos_estandarizados)]
  
  n <- length(residuos_estandarizados)

  retornos_t <- residuos_estandarizados[2:n]
  retornos_t1 <- residuos_estandarizados[1:(n-1)]
  
  u1 <- rank(retornos_t) / (length(retornos_t) + 1) # Acá se transofrman ambas series a márgenes uniformes usando CDF empírica
  u2 <- rank(retornos_t1) / (length(retornos_t1) + 1)
  U <- cbind(u1, u2)
  
  if (any(U <= 0) || any(U >= 1)) {
    U[U <= 0] <- 1e-6
    U[U >= 1] <- 1 - 1e-6
  }
  
  # Lo siguiente es ajustar varios modelos de cópula y seleccionar el mejor entre las familias 
  # "Normal", "t", "Clayton", "Gumbel" y "Frank"

  familias_copula <- c(1, 2, 3, 4, 5)  
  nombres_copula <- c("Normal", "t", "Clayton", "Gumbel", "Frank")
  
  valores_aic <- numeric(length(familias_copula))
  ajustes_copula <- list()
  
  for (i in seq_along(familias_copula)) {
    tryCatch({
      ajuste <- BiCopEst(U[,1], U[,2], family = familias_copula[i], 
                         method = "mle", max.df = 30)
      
      if (is.finite(ajuste$par) && (!is.na(ajuste$par2) && is.finite(ajuste$par2) || is.na(ajuste$par2))) {
        ajustes_copula[[i]] <- ajuste
        valores_aic[i] <- ajuste$AIC # se verifican si los parámetros son válidos
      } else {
        ajustes_copula[[i]] <- NULL
        valores_aic[i] <- Inf
      }
    }, error = function(e) {
      cat("Falló el ajuste de cópula", nombres_copula[i], ":", e$message, "\n")
      ajustes_copula[[i]] <- NULL
      valores_aic[i] <- Inf
    })
  }
  
  if (all(is.infinite(valores_aic))) { # Si todos los AIC dan infinito, entonces se recurre a cópula de independencia
    cat("Todos los modelos de cópula fallaron. Usando supuesto de independencia...\n")
    mejor_copula <- list(family = 0, par = 0, par2 = NA, AIC = 0)
    nombre_mejor_familia <- "Independencia"
    usar_copula <- FALSE
  } else {
    # En caso contrario, se selecciona la mejor cópula basada en AIC
    mejor_idx <- which.min(valores_aic)
    mejor_copula <- ajustes_copula[[mejor_idx]]
    nombre_mejor_familia <- nombres_copula[mejor_idx]
    usar_copula <- TRUE
    cat("Mejor modelo de cópula:", nombre_mejor_familia, "con AIC:", valores_aic[mejor_idx], "\n")
  }
  
  # Una vez que se tiene la cópula ajustada, se generan simulaciones usando la misma o mediante bootstrap

  if (usar_copula) {
    tryCatch({
      simulacion_copula <- BiCopSim(n_sim, mejor_copula$family, mejor_copula$par, mejor_copula$par2)
    }, error = function(e) {
      cat("Simulación desde cópula falló. Usando enfoque bootstrap...\n")
      usar_copula <<- FALSE
    })
  }
  
  if (!usar_copula) { # Este enfoque es el del bootstrap como alternativa
    indices_bootstrap <- sample(1:length(retornos_t), n_sim, replace = TRUE)
    simulacion_copula <- cbind(
      rank(retornos_t[indices_bootstrap]) / (length(retornos_t) + 1),
      rank(retornos_t1[indices_bootstrap]) / (length(retornos_t1) + 1)
    )
  }
  
  # Finalmente se transforman de vuelta a escala original y se utilizan los cuantiles, asegurándose de que los valores
  # estén en rango válido.
  retornos_ordenados <- sort(retornos_t)
  n_retornos <- length(retornos_ordenados)
  
  cuantil_empirico <- function(u) {
    u <- pmax(1e-6, pmin(1-1e-6, u))
    indices <- pmax(1, pmin(n_retornos, ceiling(u * n_retornos)))
    return(retornos_ordenados[indices])
  }
  
  retornos_sim_1 <- cuantil_empirico(simulacion_copula[,1])
  retornos_sim_2 <- cuantil_empirico(simulacion_copula[,2])
  
  retornos_portafolio_sim <- (retornos_sim_1 + retornos_sim_2) / 2
  
  resultados <- list()
  
  for (nivel_conf in niveles_confianza) {
    alpha <- 1 - nivel_conf
    
    valor_var <- quantile(retornos_portafolio_sim, alpha, type = 8, names = FALSE)
    
    valor_es <- mean(retornos_portafolio_sim[retornos_portafolio_sim <= valor_var])
    
    resultados[[paste0("VaR_", nivel_conf)]] <- valor_var
    resultados[[paste0("ES_", nivel_conf)]] <- valor_es
  }
  
  # Esta lista es la de los posibles analisis a considerar, puede ser necesario eliminar algunos.
  list(
    resultados_var_es = resultados,
    modelo_copula = mejor_copula,
    familia_copula = nombre_mejor_familia,
    modelo_garch = if(usar_garch) ajuste_garch else NULL,
    retornos_simulados = retornos_portafolio_sim,
    retornos_originales = retornos,
    usar_garch = usar_garch,
    usar_copula = usar_copula,
    comparacion_aic = data.frame(
      Familia = nombres_copula,
      AIC = valores_aic
    )
  )
}

crear_graficos_diagnosticos <- function(resultados_copula) {

  p1 <- ggplot() +
    geom_density(aes(x = resultados_copula$retornos_originales, color = "Históricos"), alpha = 0.7) +
    geom_density(aes(x = resultados_copula$retornos_simulados, color = "Simulados"), alpha = 0.7) +
    labs(title = "Distribución de Retornos Históricos vs Simulados",
         x = "Retornos", y = "Densidad") +
    theme_minimal() +
    scale_color_manual(values = c("Históricos" = "blue", "Simulados" = "red")) +
    theme(legend.title = element_blank())
  
  var_95 <- resultados_copula$resultados_var_es$VaR_0.95
  var_99 <- resultados_copula$resultados_var_es$VaR_0.99
  es_95 <- resultados_copula$resultados_var_es$ES_0.95
  es_99 <- resultados_copula$resultados_var_es$ES_0.99
  
  p2 <- ggplot() +
    geom_density(aes(x = resultados_copula$retornos_simulados), fill = "lightblue", alpha = 0.7) +
    geom_vline(aes(xintercept = var_95, color = "VaR 95%"), size = 1) +
    geom_vline(aes(xintercept = var_99, color = "VaR 99%"), size = 1) +
    geom_vline(aes(xintercept = es_95, color = "ES 95%"), size = 1, linetype = "dashed") +
    geom_vline(aes(xintercept = es_99, color = "ES 99%"), size = 1, linetype = "dashed") +
    labs(title = "Visualización de VaR y Déficit Esperado",
         x = "Retornos Simulados", y = "Densidad") +
    theme_minimal() +
    scale_color_manual(values = c("VaR 95%" = "orange", "VaR 99%" = "red",
                                  "ES 95%" = "purple", "ES 99%" = "darkred")) +
    theme(legend.title = element_blank())
  
  return(list(comparacion_distribucion = p1, grafico_var_es = p2))
}
