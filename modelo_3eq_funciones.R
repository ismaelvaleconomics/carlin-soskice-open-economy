# =============================================================================
#  MODELO DE 3 ECUACIONES EN ECONOMÍA ABIERTA — Carlin & Soskice (2023), Cap. 11
#  Archivo: modelo_3eq_funciones.R
#  Descripción: Núcleo del modelo. Define parámetros, ecuaciones y simulador.
#  Nivel: Licenciatura (últimos años) — conocimiento de ggplot2 y funciones básicas
# =============================================================================
#
#  ESTRUCTURA DEL MODELO
#  ─────────────────────
#  1. Curva de Phillips (PC):   π_t  = π_{t-1} + α·(y_t - y_e)
#  2. Regla Monetaria (MR):     y_t - y_e = -α·β·(π_t - π^T)
#  3. IS de economía abierta:   y_t = A - a·r_{t-1} + b·q_{t-1}
#  4. Condición UIP (real):     r_t - r* = q^E_{t+1} - q_t
#  5. Curva RX:                 y_t - y_e = -(a + b/(1-λ))·(r_{t-1} - r*)
#                               donde λ = 1/(1 + α²·β)
#
#  VARIABLES ENDÓGENAS (por período):
#   y  = producto (output)
#   π  = inflación
#   r  = tasa de interés real doméstica
#   q  = log del tipo de cambio real (↑q = depreciación)
#
#  PARÁMETROS:
#   y_e = producto de equilibrio (normalizado a 0 = brecha)
#   π_T = meta de inflación
#   r*  = tasa de interés real mundial (r_star)
#   q_bar = tipo de cambio real de equilibrio
#   α   = pendiente de la PC (sensibilidad inflación-brecha de producto)
#   β   = peso de la inflación en la función de pérdida del BC
#   a   = sensibilidad de la demanda a la tasa de interés
#   b   = sensibilidad de la demanda al tipo de cambio real
#   A   = demanda autónoma (desplaza la IS)
# =============================================================================

# ── Paquetes requeridos ───────────────────────────────────────────────────────
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",    quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("patchwork",quietly = TRUE)) install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
#  1. PARÁMETROS DEFAULT DEL MODELO
# =============================================================================

parametros_default <- list(
  # Parámetros estructurales
  alpha   = 0.5,    # Pendiente de la PC (α)
  beta    = 2.0,    # Peso inflación en pérdida del BC (β)
  a       = 0.5,    # Sensibilidad demanda a tasa de interés (a)
  b       = 0.3,    # Sensibilidad demanda a tipo de cambio real (b)

  # Equilibrio de mediano plazo
  y_e     = 0.0,    # Producto de equilibrio (normalizado)
  pi_T    = 2.0,    # Meta de inflación, % (π^T)
  r_star  = 2.0,    # Tasa de interés real mundial, % (r*)
  q_bar   = 0.0,    # log del tipo de cambio real de equilibrio (q̄)

  # Estado inicial (en equilibrio)
  pi_0    = 2.0,    # Inflación inicial = meta
  r_0     = 2.0,    # Tasa inicial = r*
  q_0     = 0.0,    # q inicial = q_bar
  y_0     = 0.0,    # Producto inicial = y_e

  # Horizonte de simulación
  n_periodos = 20
)

# =============================================================================
#  2. FUNCIONES DEL MODELO
# =============================================================================

#' Calcula λ (lambda) del modelo
#' λ = 1 / (1 + α²·β)
#' Aparece en la curva RX y mide cuánto de la estabilización recae en el BC.
calcular_lambda <- function(alpha, beta) {
  1 / (1 + alpha^2 * beta)
}

#' Calcula la pendiente efectiva de la curva RX
#' Coeficiente = a + b/(1-λ)
#' Siempre mayor que 'a' (la IS pura): el canal cambiario amplifica el efecto.
pendiente_RX <- function(a, b, lambda) {
  a + b / (1 - lambda)
}

#' Regla Monetaria (MR): dado π_t y la PC esperada,
#' el BC elige el output gap óptimo en t+1.
#' Devuelve el output gap deseado: (y_{t+1} - y_e)
MR_output_gap <- function(pi_t, pi_T, alpha, beta) {
  -alpha * beta * (pi_t - pi_T)
}

#' Curva RX: dado el output gap deseado, ¿qué r debe fijar el BC?
#' Despeja r_{t} de: (y_{t+1} - y_e) = -(a + b/(1-λ)) · (r_t - r*)
RX_tasa_interes <- function(gap_deseado, r_star, a, b, lambda) {
  coef <- pendiente_RX(a, b, lambda)
  r_star - gap_deseado / coef
}

#' Condición UIP real: dado el diferencial (r - r*) y el tipo de cambio
#' esperado, el tipo de cambio salta para igualar retornos.
#' q_t = q^E_{t+1} - (r_t - r*)
UIP_tipo_cambio <- function(r_t, r_star, q_esperado) {
  q_esperado - (r_t - r_star)
}

#' IS de economía abierta: calcula el producto en t dado r y q del período anterior
IS_producto <- function(A, a, r_prev, b, q_prev) {
  A - a * r_prev + b * q_prev
}

#' Curva de Phillips: inflación en t dado inflación y brecha en t
PC_inflacion <- function(pi_prev, alpha, gap) {
  pi_prev + alpha * gap
}

# =============================================================================
#  3. FUNCIÓN PRINCIPAL: SIMULADOR DEL MODELO
# =============================================================================
#
#  Recibe parámetros + descripción de un choque y simula T períodos.
#  Devuelve un data.frame con las series temporales de y, π, r, q.
#
#  TIPOS DE CHOQUE (tipo_choque):
#   "inflacion"   — choque temporal en la PC (dura 1 período)
#   "demanda_neg" — caída permanente en demanda autónoma (A baja)
#   "demanda_pos" — aumento permanente en demanda autónoma
#   "oferta"      — cambio permanente en y_e (ERU se desplaza)
#   "r_star"      — cambio en la tasa de interés mundial
# =============================================================================

simular_modelo <- function(
    params       = parametros_default,
    tipo_choque  = "inflacion",    # Tipo de choque
    magnitud     = 2.0,            # Magnitud del choque (en puntos porcentuales o unidades)
    periodo_choque = 1             # En qué período ocurre el choque (default: período 1)
) {

  # ── Desempacar parámetros ─────────────────────────────────────────────────
  alpha      <- params$alpha
  beta       <- params$beta
  a          <- params$a
  b          <- params$b
  y_e        <- params$y_e
  pi_T       <- params$pi_T
  r_star     <- params$r_star
  q_bar      <- params$q_bar
  n          <- params$n_periodos

  # ── Calcular parámetros derivados ─────────────────────────────────────────
  lambda     <- calcular_lambda(alpha, beta)
  coef_RX    <- pendiente_RX(a, b, lambda)

  # ── Inicializar vectores de resultados ────────────────────────────────────
  y     <- numeric(n + 1)   # Producto
  pi    <- numeric(n + 1)   # Inflación
  r     <- numeric(n + 1)   # Tasa de interés real doméstica
  q     <- numeric(n + 1)   # log tipo de cambio real
  A_vec <- numeric(n + 1)   # Demanda autónoma (puede cambiar con el choque)
  ye_vec<- numeric(n + 1)   # Producto de equilibrio (puede cambiar con oferta)
  rs_vec<- numeric(n + 1)   # r* (puede cambiar con choque de r*)

  # ── Condición inicial (período 0 = equilibrio) ────────────────────────────
  y[1]      <- params$y_0
  pi[1]     <- params$pi_0
  r[1]      <- params$r_0
  q[1]      <- params$q_0
  A_vec[1]  <- y_e + a * r_star - b * q_bar  # A consistente con el equilibrio
  ye_vec[1] <- y_e
  rs_vec[1] <- r_star

  # ── Aplicar choque inicial (si ocurre en período 1) ──────────────────────
  # El choque de inflación se maneja dentro del loop

  # ── Simulación período a período ─────────────────────────────────────────
  for (t in 2:(n + 1)) {

    # Período relativo al choque
    t_rel <- t - 1  # t_rel = 1 es el primer período "activo"

    # ── Actualizar parámetros según el choque ──────────────────────────────

    # Demanda autónoma A
    if (tipo_choque %in% c("demanda_neg", "demanda_pos") && t_rel >= periodo_choque) {
      signo  <- ifelse(tipo_choque == "demanda_neg", -1, 1)
      A_vec[t] <- A_vec[1] + signo * magnitud
    } else {
      A_vec[t] <- A_vec[1]
    }

    # Producto de equilibrio y_e (choque de oferta)
    if (tipo_choque == "oferta" && t_rel >= periodo_choque) {
      ye_vec[t] <- y_e + magnitud  # magnitud positiva = expansión de oferta
    } else {
      ye_vec[t] <- y_e
    }

    # Tasa mundial r* (choque externo)
    if (tipo_choque == "r_star" && t_rel >= periodo_choque) {
      rs_vec[t] <- r_star + magnitud
    } else {
      rs_vec[t] <- r_star
    }

    # ── Usar valores actualizados ───────────────────────────────────────────
    ye_t   <- ye_vec[t]
    rs_t   <- rs_vec[t]
    A_t    <- A_vec[t]

    # ── Paso 1: Calcular inflación efectiva en período t ────────────────────
    # (La IS y RX del período anterior ya determinaron y[t-1] y r[t-1])
    # El output de ESTE período es resultado de las decisiones del ANTERIOR.
    y_t_actual <- IS_producto(A_t, a, r[t-1], b, q[t-1])

    # Choque de inflación: eleva la PC temporalmente en período 'periodo_choque'
    shock_pc <- 0
    if (tipo_choque == "inflacion" && t_rel == periodo_choque) {
      shock_pc <- magnitud
    }

    pi_t_actual <- PC_inflacion(pi[t-1], alpha, y_t_actual - ye_vec[t-1]) + shock_pc

    # Guardar valores realizados en t
    y[t]  <- y_t_actual
    pi[t] <- pi_t_actual

    # ── Paso 2: BC observa π_t y elige el output gap deseado en t+1 ─────────
    # Para ello, proyecta la PC del próximo período: π_{t+1} = π_t + α·gap_{t+1}
    # La MR dice: gap_{t+1} = -α·β·(π_{t+1} - π^T)
    # Sustituyendo la PC en la MR (algebra del Capítulo 3):
    # El BC elige gap_{t+1} tal que está en la MR.

    # Proyección: ¿qué PC enfrentará en t+1?
    # π_{t+1}^E = π_t (expectativas adaptativas)
    # La MR resuelve: gap_{t+1}* = -α·β·(π_t - π^T) / (1 + α²·β)
    # (Derivación: ver Apéndice 11.5 de Carlin & Soskice)
    gap_deseado_tp1 <- -alpha * beta * (pi_t_actual - pi_T) / (1 + alpha^2 * beta)

    # ── Paso 3: BC fija r usando la curva RX ─────────────────────────────────
    # Despeja r_t de: gap_deseado = -(a + b/(1-λ))·(r_t - r*)
    r_t_nueva <- rs_t - gap_deseado_tp1 / coef_RX

    # ── Paso 4: UIP — el tipo de cambio salta ────────────────────────────────
    # El mercado forex anticipa que r > r* por varios períodos.
    # El tipo de cambio esperado en el futuro converge a q_bar del nuevo equilibrio.
    # Para simplicidad: q^E = q_bar del nuevo equilibrio (mediano plazo).
    # El nuevo q_bar depende del tipo de choque:
    q_bar_nuevo <- calcular_q_bar_nuevo(
      tipo_choque, magnitud, periodo_choque, t_rel,
      y_e, a, b, r_star, q_bar, A_vec[1], ye_vec[1], alpha, beta
    )

    q_t_nueva <- UIP_tipo_cambio(r_t_nueva, rs_t, q_bar_nuevo)

    # Guardar decisiones de política para el período siguiente
    r[t] <- r_t_nueva
    q[t] <- q_t_nueva
  }

  # ── Construir data.frame de resultados ────────────────────────────────────
  resultado <- data.frame(
    periodo   = 0:n,
    y         = y,
    pi        = pi,
    r         = r,
    q         = q,
    gap       = y - ye_vec,
    ye        = ye_vec,
    r_star    = rs_vec,
    A         = A_vec
  )

  # Añadir metadata
  attr(resultado, "params")      <- params
  attr(resultado, "tipo_choque") <- tipo_choque
  attr(resultado, "magnitud")    <- magnitud
  attr(resultado, "lambda")      <- lambda
  attr(resultado, "coef_RX")     <- coef_RX

  return(resultado)
}

# =============================================================================
#  4. FUNCIÓN AUXILIAR: CALCULA EL NUEVO q̄ DE EQUILIBRIO
# =============================================================================
#  En economía abierta, el tipo de cambio real de equilibrio cambia con
#  choques permanentes de demanda, oferta o r*.
#  Esta función calcula q_bar_nuevo dado el tipo de choque.
# =============================================================================

calcular_q_bar_nuevo <- function(
    tipo_choque, magnitud, periodo_choque, t_rel,
    y_e, a, b, r_star, q_bar, A_base, ye_base, alpha, beta
) {

  if (t_rel < periodo_choque) return(q_bar)  # Antes del choque, q_bar no cambia

  if (tipo_choque == "inflacion") {
    # Choque de inflación: el equilibrio de mediano plazo NO cambia
    # q_bar permanece igual
    return(q_bar)

  } else if (tipo_choque == "demanda_neg") {
    # Caída permanente de A: el nuevo MRE requiere mayor depreciación real
    # AD nueva: y = (A - Δ) - a·r* + b·q_bar_nuevo = y_e
    # => q_bar_nuevo = (y_e - (A_base - magnitud) + a·r*) / b
    A_nuevo <- A_base - magnitud
    return((y_e - A_nuevo + a * r_star) / b)

  } else if (tipo_choque == "demanda_pos") {
    # Aumento permanente de A: nuevo MRE requiere apreciación real
    A_nuevo <- A_base + magnitud
    return((y_e - A_nuevo + a * r_star) / b)

  } else if (tipo_choque == "oferta") {
    # Cambio en y_e: nueva ERU. La AD se ajusta con depreciación/apreciación
    # AD: y_e_nuevo = A_base - a·r* + b·q_bar_nuevo
    # => q_bar_nuevo = (y_e + magnitud - A_base + a·r*) / b
    ye_nuevo <- y_e + magnitud
    return((ye_nuevo - A_base + a * r_star) / b)

  } else if (tipo_choque == "r_star") {
    # Cambio en r*: la AD se desplaza con r*
    # y_e = A_base - a·(r* + Δ) + b·q_bar_nuevo
    # => q_bar_nuevo = (y_e - A_base + a·(r* + magnitud)) / b
    rs_nuevo <- r_star + magnitud
    return((y_e - A_base + a * rs_nuevo) / b)

  } else {
    return(q_bar)
  }
}

# =============================================================================
#  5. FUNCIONES DE VISUALIZACIÓN
# =============================================================================

# Tema visual consistente para todas las gráficas
tema_modelo <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", color = "#1F4E79", size = 13),
      plot.subtitle   = element_text(color = "#404040", size = 11),
      axis.title      = element_text(color = "#404040", size = 11),
      panel.grid.minor= element_blank(),
      panel.grid.major= element_line(color = "#E8E8E8"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background= element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.title    = element_blank()
    )
}

# Paleta de colores
AZUL    <- "#2E75B6"
ROJO    <- "#C0392B"
VERDE   <- "#27AE60"
NARANJA <- "#E67E22"
MORADO  <- "#8E44AD"
GRIS    <- "#7F8C8D"

#' Gráfica de funciones de respuesta al impulso (IRF)
#' Muestra la evolución de y, π, r y q en el tiempo.
graficar_IRF <- function(sim, titulo = "Respuesta al impulso") {

  params     <- attr(sim, "params")
  tipo_choque<- attr(sim, "tipo_choque")

  # Etiquetas según tipo de choque
  labels_choque <- c(
    "inflacion"   = "Choque de inflación (temporal)",
    "demanda_neg" = "Choque de demanda negativo (permanente)",
    "demanda_pos" = "Choque de demanda positivo (permanente)",
    "oferta"      = "Choque de oferta (permanente)",
    "r_star"      = "Cambio en tasa de interés mundial (r*)"
  )
  subtitulo <- labels_choque[tipo_choque]

  # ── Gráfica 1: Producto (brecha del producto) ──────────────────────────────
  p1 <- ggplot(sim, aes(x = periodo, y = gap)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = GRIS, linewidth = 0.6) +
    geom_line(color = AZUL, linewidth = 1.2) +
    geom_point(color = AZUL, size = 2.5) +
    labs(title = "Brecha del producto (y − yₑ)",
         x = "Período", y = "Brecha (pp)") +
    tema_modelo()

  # ── Gráfica 2: Inflación ────────────────────────────────────────────────────
  p2 <- ggplot(sim, aes(x = periodo, y = pi)) +
    geom_hline(aes(yintercept = params$pi_T),
               linetype = "dashed", color = ROJO, linewidth = 0.7) +
    geom_line(color = ROJO, linewidth = 1.2) +
    geom_point(color = ROJO, size = 2.5) +
    annotate("text", x = max(sim$periodo) * 0.85, y = params$pi_T + 0.1,
             label = "Meta π^T", color = ROJO, size = 3.5, fontface = "italic") +
    labs(title = "Inflación (π)",
         x = "Período", y = "Inflación (%)") +
    tema_modelo()

  # ── Gráfica 3: Tasa de interés real ────────────────────────────────────────
  p3 <- ggplot(sim, aes(x = periodo, y = r)) +
    geom_line(aes(y = r_star), linetype = "dashed", color = GRIS, linewidth = 0.7) +
    geom_line(color = VERDE, linewidth = 1.2) +
    geom_point(color = VERDE, size = 2.5) +
    annotate("text", x = max(sim$periodo) * 0.85,
             y = params$r_star + 0.08,
             label = "r*", color = GRIS, size = 3.5, fontface = "italic") +
    labs(title = "Tasa de interés real (r)",
         x = "Período", y = "Tasa (%)") +
    tema_modelo()

  # ── Gráfica 4: Tipo de cambio real ─────────────────────────────────────────
  p4 <- ggplot(sim, aes(x = periodo, y = q)) +
    geom_hline(yintercept = params$q_bar, linetype = "dashed",
               color = GRIS, linewidth = 0.7) +
    geom_line(color = NARANJA, linewidth = 1.2) +
    geom_point(color = NARANJA, size = 2.5) +
    annotate("text", x = max(sim$periodo) * 0.85,
             y = params$q_bar + 0.02,
             label = "q̄ inicial", color = GRIS, size = 3.5, fontface = "italic") +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    labs(title = "Tipo de cambio real (log q) — ↑ = depreciación",
         x = "Período", y = "log(Q)") +
    tema_modelo()

  # ── Combinar gráficas ───────────────────────────────────────────────────────
  combinado <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title    = titulo,
      subtitle = subtitulo,
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 15, color = "#1F4E79"),
        plot.subtitle = element_text(size = 12, color = "#404040")
      )
    )

  return(combinado)
}

#' Gráfica del diagrama IS-RX-MR (espacio r-y)
#' Muestra la curva RX, la IS inicial, y el camino de ajuste.
graficar_RX_diagram <- function(sim, tipo_choque_label = "") {

  params  <- attr(sim, "params")
  alpha   <- params$alpha
  beta    <- params$beta
  a       <- params$a
  b       <- params$b
  y_e     <- params$y_e
  r_star  <- params$r_star
  pi_T    <- params$pi_T
  lambda  <- attr(sim, "lambda")
  coef_RX <- attr(sim, "coef_RX")

  # Rango para los ejes
  r_rng <- range(sim$r, na.rm = TRUE)
  r_min <- min(r_rng[1] - 0.5, r_star - 1)
  r_max <- max(r_rng[2] + 0.5, r_star + 1)
  y_rng <- range(sim$y, na.rm = TRUE)
  y_min <- min(y_rng[1] - 0.3, y_e - 0.5)
  y_max <- max(y_rng[2] + 0.3, y_e + 0.5)

  # Grilla de r para dibujar las curvas
  r_seq <- seq(r_min, r_max, length.out = 200)

  # Curva RX: y = y_e - coef_RX * (r - r*)
  df_RX <- data.frame(
    r = r_seq,
    y = y_e - coef_RX * (r_seq - r_star)
  )

  # MR: y - y_e = -α·β·(π - π^T) → es vertical en y = y_e cuando π = π^T
  # Representamos la MR como función de π implícita; aquí la mostramos
  # como línea vertical en y_e para simplificar.

  # Camino de ajuste (r, y) período a período
  df_path <- data.frame(
    r = sim$r[-nrow(sim)],  # r fijado en t, produce y en t+1
    y = sim$y[-1],
    periodo = sim$periodo[-1]
  )

  ggplot() +
    # Curva RX
    geom_line(data = df_RX, aes(x = y, y = r, color = "Curva RX"),
              linewidth = 1.2, linetype = "solid") +
    # Línea r*
    geom_hline(yintercept = r_star, linetype = "dashed",
               color = GRIS, linewidth = 0.7) +
    # Línea y_e
    geom_vline(xintercept = y_e, linetype = "dashed",
               color = GRIS, linewidth = 0.7) +
    # Camino de ajuste
    geom_path(data = df_path, aes(x = y, y = r),
              arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
              color = AZUL, linewidth = 0.8, alpha = 0.7) +
    geom_point(data = df_path, aes(x = y, y = r, color = "Camino de ajuste"),
               size = 2.5, alpha = 0.8) +
    # Punto de equilibrio
    geom_point(aes(x = y_e, y = r_star),
               size = 4, color = "#1F4E79", shape = 18) +
    # Etiquetas de ejes
    annotate("text", x = y_e + 0.02, y = r_min + 0.1,
             label = "yₑ", color = GRIS, size = 3.5, fontface = "italic") +
    annotate("text", x = y_min + 0.05, y = r_star + 0.08,
             label = "r*", color = GRIS, size = 3.5, fontface = "italic") +
    scale_color_manual(values = c("Curva RX" = VERDE, "Camino de ajuste" = AZUL)) +
    labs(
      title    = "Diagrama IS-RX: Tasa de interés y Producto",
      subtitle = tipo_choque_label,
      x        = "Producto (y)",
      y        = "Tasa de interés real (r)"
    ) +
    coord_cartesian(xlim = c(y_min, y_max), ylim = c(r_min, r_max)) +
    tema_modelo()
}

#' Gráfica del diagrama AD-ERU (espacio q-y)
#' Muestra el equilibrio de mediano plazo y el camino de ajuste del tipo de cambio.
graficar_AD_ERU <- function(sim, tipo_choque_label = "") {

  params  <- attr(sim, "params")
  a       <- params$a
  b       <- params$b
  y_e     <- params$y_e
  r_star  <- params$r_star
  q_bar   <- params$q_bar
  tipo_ch <- attr(sim, "tipo_choque")
  mag     <- attr(sim, "magnitud")

  # Calcular nuevo q_bar de equilibrio
  A_base  <- y_e + a * r_star - b * q_bar
  q_bar_new <- calcular_q_bar_nuevo(
    tipo_ch, mag, 1, 1,
    y_e, a, b, r_star, q_bar, A_base, y_e,
    params$alpha, params$beta
  )
  ye_new <- ifelse(tipo_ch == "oferta", y_e + mag, y_e)

  # Rango
  q_rng <- range(c(sim$q, q_bar, q_bar_new), na.rm = TRUE)
  q_min <- q_rng[1] - 0.2
  q_max <- q_rng[2] + 0.2
  y_rng <- range(sim$y, na.rm = TRUE)
  y_min <- min(y_rng[1] - 0.3, y_e - 0.5)
  y_max <- max(y_rng[2] + 0.3, ye_new + 0.5)

  # AD inicial: y = A - a·r* + b·q
  q_seq <- seq(q_min, q_max, length.out = 200)
  df_AD <- data.frame(q = q_seq, y = A_base - a * r_star + b * q_seq)

  # AD nueva (si hay choque de demanda)
  A_new <- A_base + ifelse(tipo_ch == "demanda_pos", mag,
                    ifelse(tipo_ch == "demanda_neg", -mag, 0))
  rs_new <- r_star + ifelse(tipo_ch == "r_star", mag, 0)
  df_AD_new <- data.frame(q = q_seq, y = A_new - a * rs_new + b * q_seq)

  p <- ggplot() +
    # ERU inicial (vertical en y_e)
    geom_vline(xintercept = y_e, color = ROJO, linewidth = 1.1, linetype = "solid") +
    annotate("text", x = y_e + 0.03, y = q_max - 0.05,
             label = "ERU", color = ROJO, size = 3.5, fontface = "bold") +
    # AD inicial
    geom_line(data = df_AD, aes(x = y, y = q, color = "AD inicial"),
              linewidth = 1.1) +
    # Equilibrio inicial
    geom_point(aes(x = y_e, y = q_bar),
               size = 4, color = "#1F4E79", shape = 18) +
    annotate("text", x = y_e + 0.03, y = q_bar - 0.03,
             label = "A (inicial)", color = "#1F4E79", size = 3)

  # Si hay choque permanente, mostrar nueva AD / nueva ERU
  if (tipo_ch %in% c("demanda_neg", "demanda_pos", "r_star")) {
    p <- p +
      geom_line(data = df_AD_new, aes(x = y, y = q, color = "AD nueva"),
                linewidth = 1.1, linetype = "dashed") +
      geom_point(aes(x = y_e, y = q_bar_new),
                 size = 4, color = NARANJA, shape = 18) +
      annotate("text", x = y_e + 0.03, y = q_bar_new + 0.03,
               label = "Z (nuevo MRE)", color = NARANJA, size = 3)
  }

  if (tipo_ch == "oferta") {
    p <- p +
      geom_vline(xintercept = ye_new, color = ROJO,
                 linewidth = 1.1, linetype = "dashed") +
      annotate("text", x = ye_new + 0.03, y = q_max - 0.05,
               label = "ERU nueva", color = ROJO, size = 3.5, fontface = "bold") +
      geom_point(aes(x = ye_new, y = q_bar_new),
                 size = 4, color = NARANJA, shape = 18)
  }

  # Camino de ajuste en el espacio (y, q)
  df_path_aderu <- data.frame(y = sim$y, q = sim$q, periodo = sim$periodo)
  p <- p +
    geom_path(data = df_path_aderu, aes(x = y, y = q),
              arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
              color = AZUL, linewidth = 0.8, alpha = 0.6) +
    geom_point(data = df_path_aderu, aes(x = y, y = q),
               color = AZUL, size = 2, alpha = 0.7)

  p <- p +
    scale_color_manual(values = c("AD inicial" = VERDE, "AD nueva" = NARANJA)) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    labs(
      title    = "Diagrama AD–ERU: Tipo de cambio real y Producto",
      subtitle = tipo_choque_label,
      x        = "Producto (y)",
      y        = "Tipo de cambio real (log q) — ↑ depreciación"
    ) +
    coord_cartesian(xlim = c(y_min, y_max), ylim = c(q_min, q_max)) +
    tema_modelo()

  return(p)
}

#' Gráfica del diagrama PC-MR (espacio π-y)
graficar_PC_MR <- function(sim, tipo_choque_label = "") {

  params <- attr(sim, "params")
  alpha  <- params$alpha
  beta   <- params$beta
  y_e    <- params$y_e
  pi_T   <- params$pi_T

  # Rango
  y_rng  <- range(sim$y, na.rm = TRUE)
  pi_rng <- range(sim$pi, na.rm = TRUE)
  y_min  <- min(y_rng[1] - 0.2, y_e - 0.5)
  y_max  <- max(y_rng[2] + 0.2, y_e + 0.5)
  pi_min <- min(pi_rng[1] - 0.3, pi_T - 0.5)
  pi_max <- max(pi_rng[2] + 0.3, pi_T + 0.5)

  y_seq <- seq(y_min, y_max, length.out = 200)

  # MR: gap = -α·β·(π - π^T) => π = π^T - gap/(α·β)
  df_MR <- data.frame(
    y  = y_seq,
    pi = pi_T - (y_seq - y_e) / (alpha * beta)
  )

  # PCs observadas (una por período)
  # PC_t: π = π_{t-1} + α·(y - y_e)  (evalúa sobre un rango de y)
  pc_list <- lapply(seq_len(min(nrow(sim)-1, 6)), function(t) {
    data.frame(
      y      = y_seq,
      pi     = sim$pi[t] + alpha * (y_seq - y_e),
      periodo= paste0("PC período ", t - 1)
    )
  })
  df_PC <- do.call(rbind, pc_list)

  # Paleta para PCs
  n_pc <- length(unique(df_PC$periodo))
  cols_pc <- colorRampPalette(c("#D6E4F3", "#1F4E79"))(n_pc)

  # Camino de ajuste observado
  df_path <- data.frame(y = sim$y, pi = sim$pi, periodo = sim$periodo)

  ggplot() +
    # PCs
    geom_line(data = df_PC, aes(x = y, y = pi, color = periodo),
              linewidth = 0.8, alpha = 0.7) +
    scale_color_manual(values = setNames(cols_pc, unique(df_PC$periodo))) +
    # MR
    geom_line(data = df_MR, aes(x = y, y = pi),
              color = ROJO, linewidth = 1.2, linetype = "solid") +
    annotate("text", x = y_max - 0.1, y = pi_T - 0.1,
             label = "MR", color = ROJO, size = 4, fontface = "bold") +
    # Meta inflación
    geom_hline(yintercept = pi_T, linetype = "dotted",
               color = GRIS, linewidth = 0.6) +
    # Línea y_e
    geom_vline(xintercept = y_e, linetype = "dotted",
               color = GRIS, linewidth = 0.6) +
    # Camino de ajuste
    geom_path(data = df_path, aes(x = y, y = pi),
              arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
              color = AZUL, linewidth = 0.9, alpha = 0.8) +
    geom_point(data = df_path, aes(x = y, y = pi),
               color = AZUL, size = 2.5, alpha = 0.8) +
    # Punto inicial y final
    geom_point(data = df_path[1, ], aes(x = y, y = pi),
               color = "#1F4E79", size = 4, shape = 18) +
    annotate("text", x = df_path$y[1] + 0.03, y = df_path$pi[1] + 0.03,
             label = "A", color = "#1F4E79", size = 3.5, fontface = "bold") +
    coord_cartesian(xlim = c(y_min, y_max), ylim = c(pi_min, pi_max)) +
    labs(
      title    = "Diagrama PC–MR: Inflación y Producto",
      subtitle = tipo_choque_label,
      x        = "Producto (y)",
      y        = "Inflación (π, %)"
    ) +
    tema_modelo() +
    guides(color = guide_legend(nrow = 2))
}

# =============================================================================
#  6. FUNCIÓN PARA IMPRIMIR TABLA DE RESULTADOS
# =============================================================================

imprimir_tabla <- function(sim, n_filas = 10) {
  params     <- attr(sim, "params")
  tipo_choque<- attr(sim, "tipo_choque")
  magnitud   <- attr(sim, "magnitud")
  lambda     <- attr(sim, "lambda")
  coef_RX    <- attr(sim, "coef_RX")

  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  RESULTADOS DE SIMULACIÓN — Modelo 3 Ecuaciones Abierto   \n")
  cat("  Carlin & Soskice (2023), Cap. 11                         \n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("  Choque:    %s (magnitud = %.2f)\n", tipo_choque, magnitud))
  cat(sprintf("  α = %.2f   β = %.2f   a = %.2f   b = %.2f\n",
              params$alpha, params$beta, params$a, params$b))
  cat(sprintf("  λ = %.4f  |  Coef. RX = %.4f (IS pura = %.2f)\n",
              lambda, coef_RX, params$a))
  cat(sprintf("  r* = %.2f%%  |  π^T = %.2f%%  |  y_e = %.2f\n",
              params$r_star, params$pi_T, params$y_e))
  cat("───────────────────────────────────────────────────────────\n")
  cat(sprintf("  %-7s  %-8s  %-10s  %-8s  %-10s\n",
              "Período", "y (brecha)", "π (%)", "r (%)", "q (log)"))
  cat("───────────────────────────────────────────────────────────\n")

  for (i in seq_len(min(n_filas + 1, nrow(sim)))) {
    cat(sprintf("  %-7d  %+8.3f  %10.3f  %8.3f  %+10.4f\n",
                sim$periodo[i], sim$gap[i], sim$pi[i], sim$r[i], sim$q[i]))
  }
  if (nrow(sim) > n_filas + 1) {
    cat(sprintf("  ... (%d períodos más)\n", nrow(sim) - n_filas - 1))
  }
  cat("═══════════════════════════════════════════════════════════\n\n")
}

# =============================================================================
#  FIN DE modelo_3eq_funciones.R
#  Carga este archivo con: source("modelo_3eq_funciones.R")
# =============================================================================
