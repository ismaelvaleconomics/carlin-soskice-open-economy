# =============================================================================
#  MODELO DE 3 ECUACIONES EN ECONOMÍA ABIERTA — Carlin & Soskice (2023)
#  Archivo: modelo_3eq_simulaciones.R
#  Descripción: Ejemplos de simulación y visualización de los 4 tipos de choque.
#  Instrucciones: ejecutar DESPUÉS de cargar modelo_3eq_funciones.R
# =============================================================================
#
#  USO:
#    source("modelo_3eq_funciones.R")   # Cargar el núcleo del modelo
#    source("modelo_3eq_simulaciones.R") # Ejecutar las simulaciones
#
#  O ejecutar sección por sección con Ctrl+Enter en RStudio
# =============================================================================

source("modelo_3eq_funciones.R")

# =============================================================================
#  CONFIGURACIÓN GLOBAL DE PARÁMETROS
#  Puedes modificar estos valores y volver a correr para experimentar
# =============================================================================

params <- list(
  alpha      = 0.5,    # Pendiente de la PC
  beta       = 2.0,    # Peso de la inflación en la pérdida del BC
  a          = 0.5,    # Sensibilidad demanda a tasa de interés
  b          = 0.3,    # Sensibilidad demanda a tipo de cambio real
  y_e        = 0.0,    # Producto de equilibrio (normalizado)
  pi_T       = 2.0,    # Meta de inflación (%)
  r_star     = 2.0,    # Tasa de interés real mundial (%)
  q_bar      = 0.0,    # log tipo de cambio real de equilibrio
  pi_0       = 2.0,    # Inflación inicial
  r_0        = 2.0,    # Tasa inicial
  q_0        = 0.0,    # q inicial
  y_0        = 0.0,    # Producto inicial
  n_periodos = 20      # Períodos a simular
)

# Mostrar parámetros derivados clave
lambda   <- calcular_lambda(params$alpha, params$beta)
coef_rx  <- pendiente_RX(params$a, params$b, lambda)
cat("\n📐 Parámetros derivados del modelo:\n")
cat(sprintf("   λ (lambda)  = %.4f\n", lambda))
cat(sprintf("   Coef. IS    = %.4f   (solo canal de tasa de interés)\n", params$a))
cat(sprintf("   Coef. RX    = %.4f   (IS + canal cambiario)\n", coef_rx))
cat(sprintf("   La RX es %.1fx más potente que la IS pura.\n\n", coef_rx / params$a))

# =============================================================================
#  CHOQUE 1: CHOQUE DE INFLACIÓN (inflation shock)
#  ─────────────────────────────────────────────────
#  Un evento externo (ej. alza de precios de alimentos) empuja la inflación
#  2 pp por encima de la meta durante 1 período.
#  El BC y el forex reaccionan simultáneamente.
# =============================================================================

cat("═══════════════════════════════════════════════════════\n")
cat("  CHOQUE 1: Choque de inflación (+2 pp, temporal)     \n")
cat("═══════════════════════════════════════════════════════\n")

sim_inflacion <- simular_modelo(
  params       = params,
  tipo_choque  = "inflacion",
  magnitud     = 2.0,      # +2 puntos porcentuales sobre la meta
  periodo_choque = 1
)

imprimir_tabla(sim_inflacion, n_filas = 8)

# Gráficas de respuesta al impulso
irf_inflacion <- graficar_IRF(sim_inflacion,
                               titulo = "CHOQUE 1: Choque de inflación (+2 pp)")
print(irf_inflacion)

# Diagrama PC-MR
pc_mr_inflacion <- graficar_PC_MR(sim_inflacion,
                                   "Choque de inflación temporal")
print(pc_mr_inflacion)

# Diagrama IS-RX
rx_inflacion <- graficar_RX_diagram(sim_inflacion,
                                     "Choque de inflación temporal")
print(rx_inflacion)

# Diagrama AD-ERU
aderu_inflacion <- graficar_AD_ERU(sim_inflacion,
                                    "Choque de inflación temporal")
print(aderu_inflacion)

# =============================================================================
#  CHOQUE 2: CHOQUE DE DEMANDA NEGATIVO PERMANENTE (negative demand shock)
#  ─────────────────────────────────────────────────────────────────────────
#  Una caída permanente en la inversión autónoma (ej. crisis de confianza)
#  reduce la demanda agregada. El nuevo equilibrio requiere un tipo de cambio
#  real más depreciado para compensar con exportaciones netas.
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CHOQUE 2: Caída permanente de demanda (−1 unidad en A)      \n")
cat("═══════════════════════════════════════════════════════════════\n")

sim_demanda_neg <- simular_modelo(
  params       = params,
  tipo_choque  = "demanda_neg",
  magnitud     = 1.0,      # Reducción de 1 unidad en la demanda autónoma A
  periodo_choque = 1
)

imprimir_tabla(sim_demanda_neg, n_filas = 8)

irf_demanda_neg <- graficar_IRF(sim_demanda_neg,
                                 titulo = "CHOQUE 2: Caída permanente de demanda")
print(irf_demanda_neg)

pc_mr_demanda_neg <- graficar_PC_MR(sim_demanda_neg,
                                     "Choque de demanda negativo permanente")
print(pc_mr_demanda_neg)

aderu_demanda_neg <- graficar_AD_ERU(sim_demanda_neg,
                                      "Choque de demanda negativo permanente")
print(aderu_demanda_neg)

# =============================================================================
#  CHOQUE 3: CHOQUE DE OFERTA POSITIVO (positive supply shock)
#  ─────────────────────────────────────────────────────────────
#  Un avance tecnológico eleva la productividad, desplazando la ERU a la derecha.
#  El nuevo equilibrio tiene mayor y_e y un tipo de cambio real más depreciado
#  (se necesita más demanda para absorber la mayor producción potencial).
# =============================================================================

cat("════════════════════════════════════════════════════════════\n")
cat("  CHOQUE 3: Choque de oferta positivo (+0.5 en y_e)        \n")
cat("════════════════════════════════════════════════════════════\n")

sim_oferta <- simular_modelo(
  params       = params,
  tipo_choque  = "oferta",
  magnitud     = 0.5,      # y_e aumenta en 0.5 (ej. mejora tecnológica)
  periodo_choque = 1
)

imprimir_tabla(sim_oferta, n_filas = 8)

irf_oferta <- graficar_IRF(sim_oferta,
                            titulo = "CHOQUE 3: Choque de oferta positivo (↑ y_e)")
print(irf_oferta)

aderu_oferta <- graficar_AD_ERU(sim_oferta,
                                 "Choque de oferta positivo (ERU → derecha)")
print(aderu_oferta)

# =============================================================================
#  CHOQUE 4: ALZA DE LA TASA DE INTERÉS MUNDIAL r* (external rate shock)
#  ─────────────────────────────────────────────────────────────────────────
#  La Fed u otro banco central sube tasas globalmente, elevando r*.
#  Para una economía pequeña abierta, esto implica que la curva RX se desplaza
#  y el BC doméstico debe ajustar su tasa. El tipo de cambio también reacciona.
# =============================================================================

cat("════════════════════════════════════════════════════════════════\n")
cat("  CHOQUE 4: Alza de tasa de interés mundial (+1 pp en r*)     \n")
cat("════════════════════════════════════════════════════════════════\n")

sim_rstar <- simular_modelo(
  params       = params,
  tipo_choque  = "r_star",
  magnitud     = 1.0,      # r* sube 1 punto porcentual
  periodo_choque = 1
)

imprimir_tabla(sim_rstar, n_filas = 8)

irf_rstar <- graficar_IRF(sim_rstar,
                           titulo = "CHOQUE 4: Alza de tasa de interés mundial (↑ r*)")
print(irf_rstar)

aderu_rstar <- graficar_AD_ERU(sim_rstar,
                                "Alza de tasa de interés mundial (+1 pp en r*)")
print(aderu_rstar)

# =============================================================================
#  COMPARACIÓN: ECONOMÍA CERRADA vs. ABIERTA (ante el mismo choque de inflación)
#  ─────────────────────────────────────────────────────────────────────────────
#  En economía cerrada, el BC usa la IS como curva de implementación.
#  En economía abierta, usa la RX (más plana). El BC no necesita subir tanto r.
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  COMPARACIÓN: Economía Cerrada vs. Abierta               \n")
cat("  (Choque de inflación +2 pp)                             \n")
cat("═══════════════════════════════════════════════════════════\n")

# Simulación economía cerrada: b = 0 (no hay canal cambiario)
params_cerrada <- params
params_cerrada$b <- 0  # Sin canal cambiario → RX = IS

sim_cerrada <- simular_modelo(
  params       = params_cerrada,
  tipo_choque  = "inflacion",
  magnitud     = 2.0,
  periodo_choque = 1
)

# Combinar resultados
sim_inflacion$economia <- "Abierta (con canal cambiario)"
sim_cerrada$economia   <- "Cerrada (sin canal cambiario)"
df_comp <- rbind(sim_inflacion, sim_cerrada)

# Comparación de tasa de interés
p_r_comp <- ggplot(df_comp, aes(x = periodo, y = r, color = economia)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = params$r_star, linetype = "dashed",
             color = GRIS, linewidth = 0.6) +
  scale_color_manual(values = c("Abierta (con canal cambiario)" = AZUL,
                                "Cerrada (sin canal cambiario)" = ROJO)) +
  labs(
    title    = "Tasa de interés real (r) — Comparación",
    subtitle = "El BC sube MENOS en economía abierta (canal cambiario lo ayuda)",
    x = "Período", y = "Tasa de interés real (%)"
  ) +
  tema_modelo()

# Comparación de inflación (debe ser idéntica o muy similar)
p_pi_comp <- ggplot(df_comp, aes(x = periodo, y = pi, color = economia)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = params$pi_T, linetype = "dashed",
             color = GRIS, linewidth = 0.6) +
  scale_color_manual(values = c("Abierta (con canal cambiario)" = AZUL,
                                "Cerrada (sin canal cambiario)" = ROJO)) +
  labs(
    title    = "Inflación (π) — Comparación",
    subtitle = "El camino de desinflación es similar en ambas economías",
    x = "Período", y = "Inflación (%)"
  ) +
  tema_modelo()

# Tipo de cambio (solo economía abierta)
p_q_comp <- ggplot(sim_inflacion, aes(x = periodo, y = q)) +
  geom_hline(yintercept = params$q_bar, linetype = "dashed",
             color = GRIS, linewidth = 0.6) +
  geom_line(color = NARANJA, linewidth = 1.2) +
  geom_point(color = NARANJA, size = 2.5) +
  annotate("text", x = 2, y = min(sim_inflacion$q) - 0.01,
           label = "Apreciación inicial\n(overshooting)", color = NARANJA,
           size = 3.2, fontface = "italic") +
  labs(
    title    = "Tipo de cambio real (q) — Solo economía abierta",
    subtitle = "El tipo de cambio se aprecia inicialmente (overshooting) y luego deprecia",
    x = "Período", y = "log(Q) — ↑ depreciación"
  ) +
  tema_modelo()

comparacion_final <- (p_r_comp + p_pi_comp) / p_q_comp +
  plot_annotation(
    title    = "COMPARACIÓN: Economía Cerrada vs. Abierta",
    subtitle = "Choque de inflación +2 pp | Parámetros: α=0.5, β=2, a=0.5, b=0.3",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, color = "#1F4E79"),
      plot.subtitle = element_text(size = 11, color = "#404040")
    )
  )

print(comparacion_final)

cat("\n✅ Simulaciones completadas.\n")
cat("   Puedes modificar los parámetros al inicio del archivo y volver a ejecutar.\n")
cat("   Para la app interactiva: ejecuta app.R\n\n")

# =============================================================================
#  FIN DE modelo_3eq_simulaciones.R
# =============================================================================
