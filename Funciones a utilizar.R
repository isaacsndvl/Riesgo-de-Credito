library(dplyr)
library(ggplot2)
library(moments)


calcular_metricas_riesgo <- function(retornos, nivel_confianza = 0.05) {
  
  # Esta función se utiliza para calcular el VaR y el ES usando los 3 métodos: 
  # Histórico, Gaussiano y Cornish-Fisher
  
  retornos <- na.omit(retornos)
  media_retorno <- mean(retornos)
  desv_retorno <- sd(retornos)
  asimetria <- skewness(retornos)
  curtosis <- kurtosis(retornos)
  
  var_hist <- quantile(retornos, nivel_confianza) # Método Histórico
  es_hist <- mean(retornos[retornos <= var_hist])
  
  var_gaussiano <- media_retorno + qnorm(nivel_confianza) * desv_retorno # Método Gaussiano
  es_gaussiano <- media_retorno - desv_retorno * dnorm(qnorm(nivel_confianza)) / nivel_confianza

  z_alpha <- qnorm(nivel_confianza) # Método Cornish-Fisher
  ajuste_cf <- (1/6) * (z_alpha^2 - 1) * asimetria + 
    (1/24) * (z_alpha^3 - 3*z_alpha) * (curtosis - 3) - 
    (1/36) * (2*z_alpha^3 - 5*z_alpha) * asimetria^2
  
  cuantil_cf <- z_alpha + ajuste_cf
  var_cf <- media_retorno + cuantil_cf * desv_retorno

  phi_cf <- dnorm(cuantil_cf)
  es_cf_base <- media_retorno - desv_retorno * phi_cf / nivel_confianza
  
  # Se agrega una corrección adicional para la asimetría y la curtosis en el ES
  correccion_es_asimetria <- (desv_retorno * asimetria / 6) * (phi_cf * (cuantil_cf^2 - 1)) / nivel_confianza
  correccion_es_curtosis <- (desv_retorno * (curtosis - 3) / 24) * (phi_cf * (cuantil_cf^3 - 3*cuantil_cf)) / nivel_confianza
  
  es_cf <- es_cf_base + correccion_es_asimetria + correccion_es_curtosis
  
  return(list(
    historico = list(VaR = var_hist, ES = es_hist),
    gaussiano = list(VaR = var_gaussiano, ES = es_gaussiano),
    cornish_fisher = list(VaR = var_cf, ES = es_cf),
    estadisticas_resumen = list(
      media = media_retorno,
      desviacion = desv_retorno,
      asimetria = asimetria,
      curtosis = curtosis,
      n_obs = length(retornos)
    ),
    detalles_cf = list(
      z_alpha = z_alpha,
      ajuste_cf = ajuste_cf,
      cuantil_cf = cuantil_cf
    )
  ))
}

analizar_riesgo_divisa <- function(datos_divisa, niveles_confianza = c(0.01, 0.05, 0.10)) {

  retornos <- datos_divisa$Returns
  df_resultados <- data.frame() # Dataframe para guardar los resultados del ES y el VaR
  
  for (nivel_conf in niveles_confianza) {
    
    metricas_riesgo <- calcular_metricas_riesgo(retornos, nivel_conf)
    
    fila <- data.frame(
      Nivel_Confianza = paste0((1-nivel_conf)*100, "%"),
      Alpha = nivel_conf,
      VaR_Historico = metricas_riesgo$historico$VaR,
      ES_Historico = metricas_riesgo$historico$ES,
      VaR_Gaussiano = metricas_riesgo$gaussiano$VaR,
      ES_Gaussiano = metricas_riesgo$gaussiano$ES,
      VaR_CornishFisher = metricas_riesgo$cornish_fisher$VaR,
      ES_CornishFisher = metricas_riesgo$cornish_fisher$ES
    )
    
    df_resultados <- rbind(df_resultados, fila)
  }

  cat("=== Resultados VaR y ES ===\n\n")
  
  for (i in 1:nrow(df_resultados)) {
    cat(sprintf("Nivel de Confianza: %s\n", df_resultados$Nivel_Confianza[i]))
    cat(sprintf("Histórico      - VaR: %8.4f, ES: %8.4f\n", 
                df_resultados$VaR_Historico[i], df_resultados$ES_Historico[i]))
    cat(sprintf("Gaussiano      - VaR: %8.4f, ES: %8.4f\n", 
                df_resultados$VaR_Gaussiano[i], df_resultados$ES_Gaussiano[i]))
    cat(sprintf("Cornish-Fisher - VaR: %8.4f, ES: %8.4f\n", 
                df_resultados$VaR_CornishFisher[i], df_resultados$ES_CornishFisher[i]))
    cat("\n")
  }
  
  return(df_resultados)
}


graficar_comparacion_var <- function(datos_divisa, nivel_confianza = 0.05) {
  
  # Esta función se utiliza para graficar los retornos junto con sus repectivos VaRs y ES's, dado un
  # nivel de confianza en particular. Este puede cambiarse entre 0.01, 0.05 y 0.1
  
  retornos <- na.omit(datos_divisa$Returns)
  metricas_riesgo <- calcular_metricas_riesgo(retornos, nivel_confianza)
  
  p <- ggplot(data.frame(retornos = retornos), aes(x = retornos)) +
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, fill = "lightblue") +
    geom_density(color = "blue", linewidth = 1) +
    
    geom_vline(xintercept = metricas_riesgo$historico$VaR, 
               color = "red", linetype = "solid", linewidth = 1) +
    geom_vline(xintercept = metricas_riesgo$gaussiano$VaR, 
               color = "green", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = metricas_riesgo$cornish_fisher$VaR, 
               color = "purple", linetype = "dotted", linewidth = 1) +
    
    annotate("text", x = metricas_riesgo$historico$VaR, 
             y = max(density(retornos)$y) * 0.9, 
             label = sprintf("Hist: %.4f", metricas_riesgo$historico$VaR), 
             color = "red", hjust = -0.1, size = 3) +
    annotate("text", x = metricas_riesgo$gaussiano$VaR, 
             y = max(density(retornos)$y) * 0.8, 
             label = sprintf("Gauss: %.4f", metricas_riesgo$gaussiano$VaR), 
             color = "green", hjust = -0.1, size = 3) +
    annotate("text", x = metricas_riesgo$cornish_fisher$VaR, 
             y = max(density(retornos)$y) * 0.7, 
             label = sprintf("CF: %.4f", metricas_riesgo$cornish_fisher$VaR), 
             color = "purple", hjust = -0.1, size = 3) +
    
    labs(title = paste("Comparación VaR al", (1-nivel_confianza)*100, "% de Nivel de Confianza"),
         x = "Retornos",
         y = "Densidad") +
    theme_minimal()
  
  return(p)
}

comparar_metodos <- function(datos_divisa, niveles_confianza = c(0.01, 0.05, 0.10)) {
  
  # Esta función compara 1 a 1 las diferentes metodologías utilizadas. Es decir, el VaR Gaussiano contra
  # el histórico, el ES de Cornish-Fisher contra el ES histórico, y así para cada una de las combinaciones.
  
  df_resultados <- analizar_riesgo_divisa(datos_divisa, niveles_confianza)
  df_comparacion <- data.frame() # Dataframe para guardar las comparaciones de los métodos.
  
  for (i in 1:nrow(df_resultados)) {
    fila <- data.frame(
      Nivel_Confianza = df_resultados$Nivel_Confianza[i],
      VaR_Gauss_vs_Hist = df_resultados$VaR_Gaussiano[i] - df_resultados$VaR_Historico[i],
      VaR_CF_vs_Hist = df_resultados$VaR_CornishFisher[i] - df_resultados$VaR_Historico[i],
      VaR_CF_vs_Gauss = df_resultados$VaR_CornishFisher[i] - df_resultados$VaR_Gaussiano[i],
      ES_Gauss_vs_Hist = df_resultados$ES_Gaussiano[i] - df_resultados$ES_Historico[i],
      ES_CF_vs_Hist = df_resultados$ES_CornishFisher[i] - df_resultados$ES_Historico[i],
      ES_CF_vs_Gauss = df_resultados$ES_CornishFisher[i] - df_resultados$ES_Gaussiano[i]
    )
    df_comparacion <- rbind(df_comparacion, fila)
  }
  
  cat("=== Comparación de Métodos (Diferencias) ===\n")
  print(df_comparacion)
  
  return(list(resultados = df_resultados, comparacion = df_comparacion))
}
