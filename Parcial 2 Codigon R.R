# === 1. LIBRERÍAS E INSTALACIONES ===

if(!require(quantmod)) install.packages("quantmod")
if(!require(xts)) install.packages("xts")
if(!require(zoo)) install.packages("zoo")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(FinTS)) install.packages("FinTS")
if(!require(tseries)) install.packages("tseries")
if(!require(fGarch)) install.packages("fGarch")   
if(!require(tidyr)) install.packages("tidyr")
if(!require(rugarch)) install.packages("rugarch")
if(!require(broom)) install.packages("broom")

library(broom)
library(rugarch)
library(quantmod)
library(xts)
library(zoo)
library(ggplot2)
library(FinTS)
library(tseries)
library(fGarch)   
library(tidyr)

options(xts_check_TZ = FALSE)
# === 2. DESCARGA Y PREPARACIÓN DE DATOS ===

start_date <- as.Date("2018-01-01")
end_date   <- as.Date("2025-05-07")
getSymbols("^IBEX", src = "yahoo", from = start_date, to = end_date)

# Limpia precios de cierre
close_prices <- na.omit(Cl(IBEX))

# Data frame limpio
df <- data.frame(
  date  = as.Date(index(close_prices)),
  close = as.numeric(close_prices)
)

# Calcular retornos logarítmicos diarios 
log_returns <- diff(log(df$close))
fechas_log_returns <- df$date[-1]   

# Si diff (que nunca da NA salvo para precios 0 o NA) devolviera algún NA, alíneas:
valid_idx <- !is.na(log_returns)
log_returns <- log_returns[valid_idx]
fechas_log_returns <- fechas_log_returns[valid_idx]

# === 3. ANÁLISIS EXPLORATORIO ===

# (a) Gráfico precios
ggplot(df, aes(x = date, y = close)) +
  geom_line(color = "red", linewidth = 1) +
  labs(title = "IBEX 35 - Precio cierre diario", y = "Precio", x = "Fecha") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# (b) Gráfico retornos
log_returns_df <- data.frame(date = fechas_log_returns, log_return = log_returns)
ggplot(log_returns_df, aes(x = date, y = log_return)) +
  geom_line(color = "blue", alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "IBEX 35 - Retorno logarítmico diario", x = "Fecha", y = "Retorno logarítmico") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# === 4. PRUEBAS DE HETEROCEDASTICIDAD ===

returns_sq <- log_returns^2
acf(returns_sq, lag.max = 20, main = "ACF - Retornos^2 de IBEX 35")
pacf(returns_sq, lag.max = 20, main = "PACF - Retornos^2 de IBEX 35")

arch_test <- ArchTest(as.numeric(log_returns), lags = 10)
print("ARCH-LM Test:"); print(arch_test)

ljung_box <- Box.test(as.numeric(log_returns)^2, lag = 10, type = "Ljung-Box")
cat("\nLjung-Box Test sobre (Retornos)^2:\n"); print(ljung_box)

# === 5. SELECCIÓN AUTOMÁTICA DEL MEJOR MODELO GARCH(p, q) ===

max_p <- 2
max_q <- 2
results <- data.frame(p=integer(), q=integer(), AIC=numeric(), BIC=numeric())

for (p in 0:max_p) {
  for (q in 0:max_q) {
    if (p == 0 & q == 0) next  # Skip (0,0)
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
      mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
      distribution.model = "norm"
    )
    fit <- tryCatch(
      ugarchfit(spec = spec, data = as.numeric(log_returns), solver = "hybrid"),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      aic <- infocriteria(fit)[1]
      bic <- infocriteria(fit)[2]
      results <- rbind(results, data.frame(p=p, q=q, AIC=aic, BIC=bic))
    }
  }
}

# Ve los resultados ordenados por AIC y selecciona el óptimo
print(results[order(results$AIC), ])

best_p <- results$p[which.min(results$AIC)]
best_q <- results$q[which.min(results$AIC)]

cat(sprintf("\nMejor modelo según AIC: GARCH(%d,%d)\n", best_p, best_q))

# Ajusta el modelo óptimo final
spec_best <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(best_p, best_q)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)
fit_opt <- ugarchfit(spec = spec_best, data = as.numeric(log_returns), solver = "hybrid")

show(fit_opt)   # Resumen básico (sin AIC/BIC)

# Mostrar AIC y BIC
ic <- infocriteria(fit_opt)
cat(sprintf("AIC: %.4f\n", ic[1]))
cat(sprintf("BIC: %.4f\n", ic[2]))

# === 6. Volatilidad Condicional Estimada ===

volatilidad_condicional <- sigma(fit_opt)
df_vol <- data.frame(Fecha = fechas_log_returns, VolCondicional = volatilidad_condicional)

ggplot(df_vol, aes(x = Fecha, y = VolCondicional)) +
  geom_line(color = "red", size = 1) +
  labs(
    title = sprintf("Volatilidad Condicional Estimada - IBEX 35 (GARCH(%d,%d))", best_p, best_q),
    x = "Fecha", y = "Desviación estándar condicional"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# === 7. Extraer automáticamente los parámetros ajustados ===

params <- coef(fit_opt)
mu    <- ifelse("mu"    %in% names(params), params["mu"],    0)
omega <- ifelse("omega" %in% names(params), params["omega"], 0)
alpha1 <- ifelse("alpha1" %in% names(params), params["alpha1"], 0)
beta1  <- ifelse("beta1"  %in% names(params), params["beta1"], 0)

ultimo_indice <- length(log_returns)
r_tm1   <- log_returns[ultimo_indice]
eps_tm1 <- r_tm1 - mu

if (best_p == 1 && best_q == 0) {
  h_t <- omega + alpha1 * (eps_tm1^2)
} else {
  h_tm1 <- sigma(fit_opt)[ultimo_indice]^2
  h_t   <- omega + alpha1 * (eps_tm1^2) + beta1 * h_tm1
}

cat(sprintf("h_t para t = %d: %.8f\n", ultimo_indice + 1, h_t))

# === 8. Extraer la volatilidad condicional estimada del ajuste ===

volatilidad_condicional <- sigma(fit_opt)   

vol_min <- min(volatilidad_condicional)
cat(sprintf("Volatilidad mínima estimada: %.4f%%\n", vol_min * 100))

# === 9. Supuestos del ejemplo y parámetros del mejor modelo ===

# Chequeo de existencia de la media en el modelo ajustado
if ("mu" %in% names(coef(fit_opt))) {
  mu_arch <- coef(fit_opt)["mu"]
  cat("El modelo incluye media ('mu').\n")
} else {
  mu_arch <- 0
  cat("El modelo NO incluye media, se asume mu = 0.\n")
}

# Última varianza condicional estimada (de la serie sigma.t)
sigma_arch <- tail(sigma(fit_opt), 1)

# VaR condicional (ARCH)
z_99 <- qnorm(0.01) # Cuantil al 1%
var_arch <- mu_arch + z_99 * sigma_arch

# VaR tradicional
mu_trad <- mean(log_returns)
sigma_trad <- sd(log_returns)
var_trad <- mu_trad + z_99 * sigma_trad

# Mostrar ambos resultados
cat("== Comparación del Value at Risk (VaR) ==\n\n")
cat("Valor del portafolio: $100,000\n")
cat("Nivel de confianza: 99%\n\n")

cat("Método Tradicional:\n")
cat(sprintf("  VaR diario (%%) : %.4f%%\n", var_trad * 100))
cat(sprintf("  VaR diario (USD): $%.2f\n\n", var_trad * 100000))

cat("Método ARCH (condicional):\n")
cat(sprintf("  VaR diario (%%) : %.4f%%\n", var_arch * 100))
cat(sprintf("  VaR diario (USD): $%.2f\n\n", var_arch * 100000))

if (abs(var_arch) < abs(var_trad)) {
  cat("✅ El VaR ARCH estima un riesgo menor porque la volatilidad actual es baja.\n")
} else {
  cat("⚠️ El VaR ARCH estima un riesgo mayor porque la volatilidad actual está elevada.\n")
}

# === 10. Comparación de Volatilidad Condicional vs Histórica ===
 
vol_cond <- sigma(fit_opt)  # volatilidad condicional diaria
vol_hist <- sd(log_returns)   # volatilidad histórica (constante)

# Para máximo y mínimo
i_max <- which.max(vol_cond)
i_min <- which.min(vol_cond)
fecha_max <- fechas_log_returns[i_max]
fecha_min <- fechas_log_returns[i_min]
valor_max <- vol_cond[i_max]
valor_min <- vol_cond[i_min]

df_vol <- data.frame(
  Fecha = fechas_log_returns,
  Volatilidad = vol_cond
)

# Graficar comparación de volatilidad condicional vs histórica
ggplot(df_vol, aes(x = Fecha, y = Volatilidad)) +
  geom_line(color = "#DC143C", size = 1, alpha = 0.9, show.legend = TRUE) +
  geom_hline(yintercept = vol_hist, linetype = "dashed", color = "gray", size = 1, alpha = 0.7) +
  # Máximo y mínimo usando annotate:
  annotate("point", x = fecha_max, y = valor_max, color = "black", size = 3) +
  annotate("point", x = fecha_min, y = valor_min, color = "blue", size = 3) +
  annotate("text", x = fecha_max, y = valor_max, 
           label = paste0("Máximo (", as.character(fecha_max), ")"), 
           vjust = -1.2, hjust = 0.2, color = "black", size = 3.5) +
  annotate("text", x = fecha_min, y = valor_min, 
           label = paste0("Mínimo (", as.character(fecha_min), ")"), 
           vjust = 1.5, hjust = 0.3, color = "blue", size = 3.5) +
  labs(
    title = 'Comparación de Volatilidad Condicional vs Histórica - IBEX 35',
    x = 'Fecha', y = 'Volatilidad diaria (desviación estándar)'
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# === 11. Comparación del Value at Risk (VaR) con Diferentes Volatilidades ===

# Valor del portafolio y nivel de confianza
valor_portafolio <- 100000
nivel_confianza <- 0.99
z_99 <- qnorm(1 - nivel_confianza)

# Media estimada (constante del modelo)
if("mu" %in% names(coef(fit_opt))) {
  mu <- coef(fit_opt)["mu"]
} else {
  mu <- 0
}

# Volatilidad histórica
sigma_hist <- sd(log_returns)

# Volatilidad mínima y máxima estimadas con ARCH(1)
volatilidad_condicional <- sigma(fit_opt)
sigma_min <- min(volatilidad_condicional)
sigma_max <- max(volatilidad_condicional)

# VaR tradicional
var_hist <- mu + z_99 * sigma_hist

# VaR con volatilidad mínima
var_min <- mu + z_99 * sigma_min

# VaR con volatilidad máxima
var_max <- mu + z_99 * sigma_max

# Mostrar resultados
cat("== Comparación del Value at Risk (VaR) con Diferentes Volatilidades ==\n\n")
cat(sprintf("Portafolio: $%s\n", format(valor_portafolio, big.mark = ",", scientific = FALSE)))
cat(sprintf("Nivel de confianza: %d%%\n\n", nivel_confianza * 100))

cat("VaR Tradicional (volatilidad histórica):\n")
cat(sprintf("  VaR diario (%%) : %.4f%%\n", var_hist * 100))
cat(sprintf("  VaR diario (USD): $%s\n\n", format(var_hist * valor_portafolio, big.mark = ",", digits = 2, nsmall = 2)))

cat("VaR ARCH con mínima volatilidad estimada:\n")
cat(sprintf("  VaR diario (%%) : %.4f%%\n", var_min * 100))
cat(sprintf("  VaR diario (USD): $%s\n\n", format(var_min * valor_portafolio, big.mark = ",", digits = 2, nsmall = 2)))

cat("VaR ARCH con máxima volatilidad estimada:\n")
cat(sprintf("  VaR diario (%%) : %.4f%%\n", var_max * 100))
cat(sprintf("  VaR diario (USD): $%s\n", format(var_max * valor_portafolio, big.mark = ",", digits = 2, nsmall = 2)))

# === 12. Ajusta un GARCH(1,1) SÓLO si el óptimo no es GARCH(1,1) ===

if (best_p != 1 || best_q != 1) {
  cat("\nEl modelo óptimo NO es GARCH(1,1), ajustando GARCH(1,1) clásico para comparación...\n\n")
  spec_garch <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm"
  )
  fit_garch <- ugarchfit(spec = spec_garch, data = as.numeric(log_returns), solver = "hybrid")
  show(fit_garch)
} else {
  cat("\nEl modelo óptimo ya ES GARCH(1,1), el ajuste adicional no es necesario.\n")
  # Si ya es GARCH(1,1), volvemos a usar fit_opt como fit_garch para comparación posterior
  fit_garch <- fit_opt
}

# === 13. Comparacion entre ARCH y un GARCH ===

# Ajustar ARCH(1)
spec_arch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 0)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)
fit_arch <- ugarchfit(spec = spec_arch, data = as.numeric(log_returns), solver = "hybrid")
# Extraer volatilidad condicional estimada para ARCH(1)
vol_arch <- sigma(fit_arch)

# Extraer volatilidad del GARCH(1,1) ajustado antes (fit_garch)
vol_garch <- sigma(fit_garch)   # fit_garch fue ajustado arriba según corresponda

# Crear dataframe para comparación
vol_df <- data.frame(
  Fecha = fechas_log_returns,
  ARCH_1 = vol_arch,
  GARCH_1_1 = vol_garch
)

# Convertir a formato largo 
vol_df_long <- pivot_longer(
  vol_df,
  cols = c("ARCH_1", "GARCH_1_1"),
  names_to = "Modelo",
  values_to = "Volatilidad"
)
vol_df_long$Modelo <- as.factor(vol_df_long$Modelo)

# Graficar volatilidad condicional de ambos modelos
ggplot(vol_df_long, aes(x = Fecha, y = Volatilidad, color = Modelo)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("ARCH_1" = "tomato", "GARCH_1_1" = "navy")) +
  labs(
    title = "Comparación de Volatilidad Condicional: ARCH(1) vs GARCH(1,1)",
    x = "Fecha",
    y = "Desviación estándar condicional"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = "Modelo"))

# === 14. Comparación entre ARCH y GARCH en VaR ===

# --- Parámetros generales ---
valor_portafolio <- 100000
nivel_confianza <- 0.99
z_99 <- qnorm(1 - nivel_confianza)  # cuantil al 1%

# --- Media estimada ---
mu_arch <- if("mu" %in% names(coef(fit_arch))) coef(fit_arch)["mu"] else 0
mu_garch <- if("mu" %in% names(coef(fit_garch))) coef(fit_garch)["mu"] else 0

# --- Volatilidad condicional del último día para cada modelo ---
# (Recuerda: sigma.t es la desviación estándar diaria estimada por el modelo)
sigma_arch_last <- tail(sigma(fit_arch), 1)
sigma_garch_last <- tail(sigma(fit_garch), 1)

# --- Cálculo del VaR diario ---
var_arch <- mu_arch + z_99 * sigma_arch_last
var_garch <- mu_garch + z_99 * sigma_garch_last

# --- Mostrar resultados ---
cat("== Comparación del Value at Risk (VaR) - ARCH vs GARCH ==\n\n")
cat(sprintf("Portafolio: $%s\n", format(valor_portafolio, big.mark = ",", scientific = FALSE)))
cat(sprintf("Nivel de confianza: %d%%\n\n", nivel_confianza * 100))

cat("Modelo ARCH(1):\n")
cat(sprintf("  Volatilidad estimada: %.4f%%\n", sigma_arch_last * 100))
cat(sprintf("  VaR diario (%%): %.4f%%\n", var_arch * 100))
cat(sprintf("  VaR diario (USD): $%s\n\n", 
            format(var_arch * valor_portafolio, big.mark = ",", nsmall = 2)))

cat("Modelo GARCH(1,1):\n")
cat(sprintf("  Volatilidad estimada: %.4f%%\n", sigma_garch_last * 100))
cat(sprintf("  VaR diario (%%): %.4f%%\n", var_garch * 100))
cat(sprintf("  VaR diario (USD): $%s\n", 
            format(var_garch * valor_portafolio, big.mark = ",", nsmall = 2)))

# === 15. Simulación y escenarios con el modelo óptimo ===

# --- Volatilidad condicional diaria (sigma_t) ---
sigma_t <- sigma(fit_opt)
cat(sprintf("   Volatilidad condicional estimada (sigma_t), último valor: %.6f\n", tail(sigma_t, 1)))

# --- Varianza condicional estimada (h_t = sigma_t^2) ---
ht <- sigma_t^2
cat("   Varianza condicional estimada (h_t), últimos 5 valores:\n")
print(tail(ht, 5))

# --- Residuos crudos (no estandarizados) ---
et <- residuals(fit_opt)
cat("   Residuos del modelo (et), últimos 5 valores:\n")
print(tail(et, 5))

# --- Residuos Estandarizados ---
z_t <- et / sigma_t
cat("   Residuos estandarizados (zt), últimos 5 valores:\n")
print(tail(z_t, 5))

# === 16. Simulación de 50 trayectorias de retornos (modelo óptimo) ===

params <- coef(fit_opt)
mu    <- if("mu"    %in% names(params)) params["mu"]    else 0
omega <- if("omega" %in% names(params)) params["omega"] else stop("omega no está en el modelo óptimo.")
alpha <- if("alpha1" %in% names(params)) params["alpha1"] else stop("alpha1 no está en el modelo óptimo.")
beta  <- if("beta1" %in% names(params)) params["beta1"]  else 0

eps_t <- tail(residuals(fit_opt), 1)                # Último residuo (r_t - mu)
ht    <- tail(sigma(fit_opt), 1)^2                # Última varianza condicional estimada (sigma^2)

set.seed(42)
n_simulaciones <- 1000
horizonte <- 10

simulaciones <- matrix(NA, nrow = n_simulaciones, ncol = horizonte)

for (i in 1:n_simulaciones) {
  h_sim <- numeric(horizonte)
  r_sim <- numeric(horizonte)
  z <- rnorm(horizonte)
  
  # Día 1
  h_sim[1] <- omega + alpha * eps_t^2 + beta * ht
  r_sim[1] <- mu + sqrt(h_sim[1]) * z[1]
  
  # Días 2 al horizonte
  for (t in 2:horizonte) {
    h_sim[t] <- omega + alpha * (r_sim[t-1] - mu)^2 + beta * h_sim[t-1]
    r_sim[t] <- mu + sqrt(h_sim[t]) * z[t]
  }
  
  simulaciones[i, ] <- r_sim
}

matplot(
  t(simulaciones[1:50, ]), type = "l", lty = 1, col = rgb(0.5,0.5,0.5,0.3),
  xlab = "Día futuro",
  ylab = "Retorno simulado",
  main = "Simulación de 50 trayectorias de retornos (modelo óptimo)"
)
grid()

# === 17. Retornos Simulados con modelo óptimo  ===

# Calcular media y percentiles por día:
medias <- apply(simulaciones, 2, mean)
p5   <- apply(simulaciones, 2, quantile, probs = 0.05)
p25  <- apply(simulaciones, 2, quantile, probs = 0.25)
p75  <- apply(simulaciones, 2, quantile, probs = 0.75)
p95  <- apply(simulaciones, 2, quantile, probs = 0.95)

dias <- 1:horizonte

# FAN CHART 
plot(
  dias, medias, type = "l", lwd = 2, col = "black",
  ylim = range(p5, p95),
  xlab = "Día futuro", ylab = "Retorno",
  main = "Fan Chart - Retornos Simulados con modelo óptimo"
)

# Rellenar bandas de incertidumbre
polygon(c(dias, rev(dias)), c(p25, rev(p75)), col = rgb(0, 0, 1, 0.2), border = NA)
polygon(c(dias, rev(dias)), c(p5, rev(p95)), col = rgb(0, 0, 1, 0.1), border = NA)

# Re-dibujar la media encima de las bandas
lines(dias, medias, lwd = 2, col = "black")

legend(
  "topright",
  legend = c("Media simulada", "50% IC", "90% IC"),
  col = c("black", rgb(0,0,1,0.2), rgb(0,0,1,0.1)),
  lwd = c(2, 10, 10), pch = NA, pt.cex = 2, bty = "n"
)
grid()

# === 18. Distribución del Retorno Acumulado en 10 Días (1000 simulaciones)  ===

# === Histograma del retorno acumulado al día 10 ===
retorno_acumulado <- rowSums(simulaciones)  # suma de cada simulación a lo largo de los 10 días

# Retorno acumulado en escala porcentual
retorno_acumulado <- rowSums(simulaciones) * 100

# Calcular percentiles
perc_1 <- quantile(retorno_acumulado, 0.01)
perc_99 <- quantile(retorno_acumulado, 0.99)

# Dibujar histograma
hist(
  retorno_acumulado,
  breaks = 40,
  col = "slateblue",
  border = "black",
  main = "Distribución del Retorno Acumulado en 10 Días (1000 simulaciones)",
  xlab = "Retorno acumulado (%)",
  ylab = "Frecuencia"
)
abline(v = perc_1, col = "red", lty = 2, lwd = 2)    # Percentil 1%
abline(v = perc_99, col = "green", lty = 2, lwd = 2) # Percentil 99%

legend(
  "topright",
  legend = c("Percentil 1%", "Percentil 99%"),
  col = c("red", "green"),
  lwd = 2,
  lty = 2,
  bty = "n"
)
grid()

# === 19. Pronósticos ===

# Crear dataframe de precios de cierre
df_precios <- data.frame(
  Fecha = index(close_prices),
  close = as.numeric(close_prices)
)

# -----------------------------------------------
# Pronóstico de 1 día adelante con rugarch
# -----------------------------------------------
# Usamos ugarchforecast y extraemos los resultados con fitted() y sigma().
# -----------------------------------------------

pred_1d <- ugarchforecast(fit_opt, n.ahead = 1)      # Pronóstico para 1 día adelante

# Extraer el retorno esperado y volatilidad pronosticada
predicted_return <- as.numeric(fitted(pred_1d))
predicted_volatility <- as.numeric(sigma(pred_1d))

cat(sprintf("Pronóstico para el 08/05/2025:\nRetorno esperado: %.6f\nVolatilidad esperada: %.6f\n", 
            predicted_return, predicted_volatility))

# Calcular los puntos predichos usando el último punto
last_price <- tail(df_precios$close, 1)
predicted_price <- last_price * exp(predicted_return)

cat(sprintf("Punto predicho: %.2f\n", predicted_price))
grid()

# -----------------------------------------------
# Pronóstico para los próximos 10 días con el modelo óptimo
# -----------------------------------------------
# Usamos ugarchforecast con n.ahead = 10, y extraemos los retornos y volatilidad esperados.
# -----------------------------------------------

pred_multi <- ugarchforecast(fit_opt, n.ahead = 10)

# Extraer retornos y volatilidades pronosticados (vectores de longitud 10)
predicted_returns <- as.numeric(fitted(pred_multi))
predicted_volatilities <- as.numeric(sigma(pred_multi))

# Calcular puntos proyectados usando el último punto observado
last_price <- tail(df_precios$close, 1)
predicted_prices <- last_price * exp(cumsum(predicted_returns))

# Crear fechas futuras (ajusta según la última fecha real observada)
fecha_inicio <- as.Date("2025-05-07")
future_dates <- seq(fecha_inicio, by = "days", length.out = 10)

# DataFrame con pronóstico para los próximos 10 días
df_forecast <- data.frame(
  Fecha = future_dates,
  Retorno = predicted_returns,
  Punto_Pronosticado = predicted_prices,
  Volatilidad = predicted_volatilities
)

print(df_forecast)

# Graficar puntos pronosticados
plot(
  df_forecast$Fecha, df_forecast$Precio_Pronosticado, type = "l",
  col = "blue", lwd = 2, xlab = "Fecha", ylab = "Punto proyectado",
  main = "Pronóstico de Puntos para los Próximos 10 Días"
)
points(df_forecast$Fecha, df_forecast$Precio_Pronosticado, pch = 16, col = "darkblue")
grid()

# -----------------------------------------------
# Graficar pronóstico para los próximos 10 días + un hsitorico de un año
# -----------------------------------------------

# Definir fechas finales y de corte
fecha_final <- as.Date("2025-05-07")
fecha_inicio <- fecha_final - 365

# Filtrar el histórico sólo al último año
df_precios_filtrado <- df_precios[df_precios$Fecha >= fecha_inicio, ]

# Preparar el dataframe combinado (histórico + pronóstico)
df_plot <- data.frame(
  Fecha = c(df_precios_filtrado$Fecha, df_forecast$Fecha),
  Precio = c(df_precios_filtrado$close, df_forecast$Precio_Pronosticado),
  Tipo = c(rep("Histórico", nrow(df_precios_filtrado)), rep("Pronóstico", nrow(df_forecast)))
)

# Intervalos de confianza para los 10 días futuros
alpha <- 1.96
df_plot$IC_sup <- df_plot$Precio
df_plot$IC_inf <- df_plot$Precio
pron <- df_plot$Tipo == "Pronóstico"
df_plot$IC_sup[pron] <- df_forecast$Precio_Pronosticado * exp(alpha * df_forecast$Volatilidad)
df_plot$IC_inf[pron] <- df_forecast$Precio_Pronosticado * exp(-alpha * df_forecast$Volatilidad)

# Graficar
ggplot(df_plot, aes(x = Fecha, y = Precio)) +
  geom_line(data = subset(df_plot, Tipo == "Histórico"), color = "black") +
  geom_line(data = subset(df_plot, Tipo == "Pronóstico"), color = "blue") +
  geom_ribbon(
    data = subset(df_plot, Tipo == "Pronóstico"),
    aes(ymin = IC_inf, ymax = IC_sup),
    fill = "blue", alpha = 0.2
  ) +
  labs(
    title = "Pronóstico de puntos con intervalo de confianza (último año + 10 días)",
    x = "Fecha", y = "Puntos"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# Definir fechas finales y de corte para el histórico (último año)
fecha_final <- as.Date("2025-05-07")    # Ajusta según tus datos
fecha_inicio <- fecha_final - 365

# Filtrar retornos históricos sólo del último año
df_retornos_hist_un_ano <- data.frame(
  Fecha = fechas_log_returns,
  Retorno = log_returns
)
df_retornos_hist_un_ano <- subset(df_retornos_hist_un_ano, Fecha >= fecha_inicio)
df_retornos_hist_un_ano$Tipo <- "Histórico"

# Data frame con el pronóstico para 10 días adelante (ya tienes)
df_retornos_pred <- data.frame(
  Fecha = future_dates,
  Retorno = predicted_returns,
  Volatilidad = predicted_volatilities,
  Tipo = "Pronóstico"
)

alpha <- 1.96  # nivel confianza 95%

# Calcular intervalos de confianza para el pronóstico
df_retornos_pred$IC_sup <- df_retornos_pred$Retorno + alpha * df_retornos_pred$Volatilidad
df_retornos_pred$IC_inf <- df_retornos_pred$Retorno - alpha * df_retornos_pred$Volatilidad

# Agregar columnas vacías y columna Volatilidad al histórico para igualar estructuras
df_retornos_hist_un_ano$IC_sup <- NA
df_retornos_hist_un_ano$IC_inf <- NA
df_retornos_hist_un_ano$Volatilidad <- NA

# Reordenar columnas para que coincidan
df_retornos_hist_un_ano <- df_retornos_hist_un_ano[, names(df_retornos_pred)]

# Combinar histórico de 1 año + pronóstico
df_retornos_plot <- rbind(df_retornos_hist_un_ano, df_retornos_pred)

# Graficar con ggplot2
ggplot(df_retornos_plot, aes(x = Fecha, y = Retorno)) +
  geom_line(data = subset(df_retornos_plot, Tipo == "Histórico"), color = "blue") +
  geom_line(data = subset(df_retornos_plot, Tipo == "Pronóstico"), color = "red") +
  geom_ribbon(data = subset(df_retornos_plot, Tipo == "Pronóstico"),
              aes(ymin = IC_inf, ymax = IC_sup),
              fill = "red", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "IBEX 35 - Retornos logarítmicos: Último año + Pronóstico 10 días",
    x = "Fecha",
    y = "Retorno logarítmico diario"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))