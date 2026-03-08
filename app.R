# =============================================================================
#  MODELO DE 3 ECUACIONES EN ECONOMÍA ABIERTA — Carlin & Soskice (2023)
#  Archivo: app.R
#  Descripción: Shiny app interactiva. Simula el modelo en tiempo real
#               con sliders para parámetros y choques.
#
#  USO:
#    1. Abre RStudio
#    2. Asegúrate de que modelo_3eq_funciones.R esté en el mismo directorio
#    3. Abre este archivo y haz clic en "Run App"
#
#  Paquetes necesarios: shiny, ggplot2, dplyr, tidyr, patchwork, bslib
# =============================================================================

# ── Paquetes ─────────────────────────────────────────────────────────────────
for (pkg in c("shiny", "ggplot2", "dplyr", "tidyr", "patchwork", "bslib")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(bslib)

# Cargar el núcleo del modelo
source("modelo_3eq_funciones.R")

# =============================================================================
#  INTERFAZ DE USUARIO (UI)
# =============================================================================

ui <- page_sidebar(
  title = "Modelo de 3 Ecuaciones en Economía Abierta",

  # ── Tema visual ─────────────────────────────────────────────────────────────
  theme = bs_theme(
    version     = 5,
    bg          = "#FFFFFF",
    fg          = "#404040",
    primary     = "#2E75B6",
    secondary   = "#6C757D",
    base_font   = font_google("Inter"),
    heading_font= font_google("Inter"),
    bootswatch  = "flatly"
  ),

  # ── Panel lateral (controles) ───────────────────────────────────────────────
  sidebar = sidebar(
    width = 320,
    bg    = "#F8FAFD",

    # ── Tipo de choque ────────────────────────────────────────────────────────
    h5("🌩️ Tipo de Choque", style = "color:#1F4E79; font-weight:bold;"),
    selectInput(
      "tipo_choque", label = NULL,
      choices = c(
        "Choque de inflación (temporal)"           = "inflacion",
        "Caída de demanda (permanente)"            = "demanda_neg",
        "Aumento de demanda (permanente)"          = "demanda_pos",
        "Choque de oferta / ↑ productividad"       = "oferta",
        "↑ Tasa de interés mundial (r*)"           = "r_star"
      ),
      selected = "inflacion"
    ),
    sliderInput("magnitud", "Magnitud del choque:",
                min = 0.1, max = 5, value = 2, step = 0.1),
    sliderInput("periodo_choque", "Período en que ocurre el choque:",
                min = 1, max = 5, value = 1, step = 1),

    hr(style = "border-color:#D6E4F3;"),

    # ── Parámetros del BC y la economía ───────────────────────────────────────
    h5("🔧 Parámetros del Banco Central", style = "color:#1F4E79; font-weight:bold;"),
    sliderInput("alpha", withMathJax("Pendiente de la PC (α):"),
                min = 0.1, max = 1.5, value = 0.5, step = 0.05),
    sliderInput("beta", withMathJax("Peso inflación en pérdida del BC (β):"),
                min = 0.5, max = 5, value = 2.0, step = 0.1),

    hr(style = "border-color:#D6E4F3;"),

    h5("📈 Parámetros de Demanda", style = "color:#1F4E79; font-weight:bold;"),
    sliderInput("a", "Sensibilidad demanda a r (a):",
                min = 0.1, max = 2, value = 0.5, step = 0.05),
    sliderInput("b", "Sensibilidad demanda a q (b):",
                min = 0.0, max = 1.5, value = 0.3, step = 0.05),

    hr(style = "border-color:#D6E4F3;"),

    h5("🌍 Equilibrio de Mediano Plazo", style = "color:#1F4E79; font-weight:bold;"),
    sliderInput("pi_T", "Meta de inflación π^T (%):",
                min = 0, max = 5, value = 2, step = 0.5),
    sliderInput("r_star", "Tasa de interés mundial r* (%):",
                min = 0, max = 6, value = 2, step = 0.25),
    sliderInput("n_periodos", "Períodos a simular:",
                min = 10, max = 40, value = 20, step = 2),

    hr(style = "border-color:#D6E4F3;"),

    # ── Botón de reset ────────────────────────────────────────────────────────
    actionButton("reset", "↺ Restaurar valores default",
                 class = "btn btn-outline-primary btn-sm w-100")
  ),

  # ── Panel principal ─────────────────────────────────────────────────────────
  # Tres filas:
  #   (1) Tarjetas de métricas clave
  #   (2) Gráficas de respuesta al impulso (4 paneles)
  #   (3) Diagramas IS-RX, AD-ERU, PC-MR

  # Métricas clave
  layout_columns(
    col_widths = c(3, 3, 3, 3),

    value_box(
      title     = "λ (lambda)",
      value     = textOutput("val_lambda"),
      showcase  = bsicons::bs_icon("calculator"),
      theme     = "primary",
      p("Fracción ajustada por BC")
    ),
    value_box(
      title     = "Coeficiente RX",
      value     = textOutput("val_rx"),
      showcase  = bsicons::bs_icon("arrow-left-right"),
      theme     = "success",
      p("vs. IS pura (= a)")
    ),
    value_box(
      title     = "Tasa pico (r máx.)",
      value     = textOutput("val_r_pico"),
      showcase  = bsicons::bs_icon("graph-up"),
      theme     = "warning",
      p("En el período de mayor ajuste")
    ),
    value_box(
      title     = "Overshooting q",
      value     = textOutput("val_overshoot"),
      showcase  = bsicons::bs_icon("currency-exchange"),
      theme     = "danger",
      p("Desviación máx. del tipo de cambio")
    )
  ),

  # Pestañas de gráficas
  navset_card_underline(
    title = NULL,

    # ── Pestaña 1: IRF ────────────────────────────────────────────────────────
    nav_panel(
      title = "📊 Respuesta al Impulso (IRF)",
      plotOutput("plot_irf", height = "550px")
    ),

    # ── Pestaña 2: PC-MR ──────────────────────────────────────────────────────
    nav_panel(
      title = "📉 Diagrama PC–MR",
      plotOutput("plot_pcmr", height = "480px")
    ),

    # ── Pestaña 3: IS-RX ──────────────────────────────────────────────────────
    nav_panel(
      title = "📈 Diagrama IS–RX",
      plotOutput("plot_rx", height = "480px")
    ),

    # ── Pestaña 4: AD-ERU ─────────────────────────────────────────────────────
    nav_panel(
      title = "🌐 Diagrama AD–ERU",
      plotOutput("plot_aderu", height = "480px")
    ),

    # ── Pestaña 5: Comparación cerrada vs. abierta ────────────────────────────
    nav_panel(
      title = "⚖️ Cerrada vs. Abierta",
      plotOutput("plot_comp", height = "550px")
    ),

    # ── Pestaña 6: Tabla de datos ─────────────────────────────────────────────
    nav_panel(
      title = "📋 Tabla de Resultados",
      br(),
      div(
        style = "background:#F8FAFD; padding:12px; border-radius:6px; font-size:0.85em;",
        h6("ℹ️ Interpretación de columnas:", style = "color:#1F4E79;"),
        tags$ul(
          tags$li(strong("gap"), " = y − yₑ: brecha del producto (positivo = expansión)"),
          tags$li(strong("q"), " = log(Q): ↑ depreciación del tipo de cambio real"),
          tags$li(strong("r"), " = tasa de interés real doméstica fijada por el BC"),
          tags$li(strong("π"), " = inflación del período")
        )
      ),
      br(),
      tableOutput("tabla_resultados")
    )
  )
)

# =============================================================================
#  SERVIDOR (SERVER)
# =============================================================================

server <- function(input, output, session) {

  # ── Reset a valores default ────────────────────────────────────────────────
  observeEvent(input$reset, {
    updateSliderInput(session, "alpha",        value = 0.5)
    updateSliderInput(session, "beta",         value = 2.0)
    updateSliderInput(session, "a",            value = 0.5)
    updateSliderInput(session, "b",            value = 0.3)
    updateSliderInput(session, "pi_T",         value = 2.0)
    updateSliderInput(session, "r_star",       value = 2.0)
    updateSliderInput(session, "magnitud",     value = 2.0)
    updateSliderInput(session, "periodo_choque", value = 1)
    updateSliderInput(session, "n_periodos",   value = 20)
    updateSelectInput(session, "tipo_choque",  selected = "inflacion")
  })

  # ── Función reactiva: construir parámetros desde los inputs ───────────────
  params_reactivos <- reactive({
    list(
      alpha      = input$alpha,
      beta       = input$beta,
      a          = input$a,
      b          = input$b,
      y_e        = 0.0,
      pi_T       = input$pi_T,
      r_star     = input$r_star,
      q_bar      = 0.0,
      pi_0       = input$pi_T,   # Empieza en meta
      r_0        = input$r_star, # Empieza en r*
      q_0        = 0.0,
      y_0        = 0.0,
      n_periodos = input$n_periodos
    )
  })

  # ── Simulación principal ──────────────────────────────────────────────────
  sim_principal <- reactive({
    simular_modelo(
      params         = params_reactivos(),
      tipo_choque    = input$tipo_choque,
      magnitud       = input$magnitud,
      periodo_choque = input$periodo_choque
    )
  })

  # ── Simulación economía cerrada (para comparación) ────────────────────────
  sim_cerrada <- reactive({
    params_c <- params_reactivos()
    params_c$b <- 0   # Sin canal cambiario
    simular_modelo(
      params         = params_c,
      tipo_choque    = input$tipo_choque,
      magnitud       = input$magnitud,
      periodo_choque = input$periodo_choque
    )
  })

  # ── Etiqueta del choque ───────────────────────────────────────────────────
  label_choque <- reactive({
    labels <- c(
      "inflacion"   = "Choque de inflación (temporal)",
      "demanda_neg" = "Caída permanente de demanda",
      "demanda_pos" = "Aumento permanente de demanda",
      "oferta"      = "Choque de oferta positivo",
      "r_star"      = "Alza de tasa mundial (r*)"
    )
    labels[input$tipo_choque]
  })

  # ── Métricas clave ────────────────────────────────────────────────────────
  output$val_lambda <- renderText({
    l <- calcular_lambda(input$alpha, input$beta)
    sprintf("%.4f", l)
  })

  output$val_rx <- renderText({
    l   <- calcular_lambda(input$alpha, input$beta)
    rx  <- pendiente_RX(input$a, input$b, l)
    sprintf("%.3f  (%.1fx IS)", rx, rx / input$a)
  })

  output$val_r_pico <- renderText({
    sim <- sim_principal()
    r_pico <- max(abs(sim$r - input$r_star))
    r_max  <- sim$r[which.max(abs(sim$r - input$r_star))]
    sprintf("%.2f%%", r_max)
  })

  output$val_overshoot <- renderText({
    sim <- sim_principal()
    q_bar_final <- tail(sim$q, 3) |> mean()
    q_max_dev   <- sim$q[which.max(abs(sim$q - q_bar_final))]
    sprintf("Δq = %+.4f", q_max_dev - q_bar_final)
  })

  # ── Gráficas ──────────────────────────────────────────────────────────────
  output$plot_irf <- renderPlot({
    sim    <- sim_principal()
    titulo <- sprintf("Respuesta al Impulso | %s | Magnitud = %.1f",
                      label_choque(), input$magnitud)
    graficar_IRF(sim, titulo = titulo)
  }, res = 110)

  output$plot_pcmr <- renderPlot({
    graficar_PC_MR(sim_principal(), label_choque())
  }, res = 110)

  output$plot_rx <- renderPlot({
    graficar_RX_diagram(sim_principal(), label_choque())
  }, res = 110)

  output$plot_aderu <- renderPlot({
    graficar_AD_ERU(sim_principal(), label_choque())
  }, res = 110)

  # Comparación cerrada vs. abierta
  output$plot_comp <- renderPlot({
    sim_ab <- sim_principal()
    sim_ce <- sim_cerrada()

    sim_ab$economia <- "Abierta (con canal cambiario)"
    sim_ce$economia <- "Cerrada (sin canal cambiario)"
    df_comp <- rbind(sim_ab, sim_ce)

    p1 <- ggplot(df_comp, aes(x = periodo, y = r, color = economia)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
      geom_hline(yintercept = input$r_star, linetype="dashed", color=GRIS) +
      scale_color_manual(values = c("Abierta (con canal cambiario)" = AZUL,
                                    "Cerrada (sin canal cambiario)" = ROJO)) +
      labs(title = "Tasa de interés real (r)",
           subtitle = "El BC sube MENOS en economía abierta",
           x = "Período", y = "r (%)") + tema_modelo()

    p2 <- ggplot(df_comp, aes(x = periodo, y = pi, color = economia)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
      geom_hline(yintercept = input$pi_T, linetype="dashed", color=GRIS) +
      scale_color_manual(values = c("Abierta (con canal cambiario)" = AZUL,
                                    "Cerrada (sin canal cambiario)" = ROJO)) +
      labs(title = "Inflación (π)",
           subtitle = "El camino de desinflación es similar",
           x = "Período", y = "π (%)") + tema_modelo()

    p3 <- ggplot(sim_ab, aes(x = periodo, y = q)) +
      geom_hline(yintercept = 0, linetype="dashed", color=GRIS) +
      geom_line(color = NARANJA, linewidth = 1.2) + geom_point(color=NARANJA, size=2.5) +
      labs(title = "Tipo de cambio real (q) — Solo economía abierta",
           subtitle = "Apreciación inicial (overshooting) → depreciación gradual",
           x = "Período", y = "log(Q) — ↑ depreciación") + tema_modelo()

    (p1 + p2) / p3 +
      plot_annotation(
        title    = "Economía Cerrada vs. Abierta",
        subtitle = sprintf("%s | Magnitud = %.1f | α=%.2f β=%.2f a=%.2f b=%.2f",
                           label_choque(), input$magnitud,
                           input$alpha, input$beta, input$a, input$b),
        theme = theme(
          plot.title    = element_text(face="bold", size=14, color="#1F4E79"),
          plot.subtitle = element_text(size=10, color="#404040")
        )
      )
  }, res = 110)

  # ── Tabla de resultados ───────────────────────────────────────────────────
  output$tabla_resultados <- renderTable({
    sim <- sim_principal()
    sim |>
      select(periodo, gap, pi, r, q) |>
      mutate(
        gap = round(gap, 4),
        pi  = round(pi,  3),
        r   = round(r,   3),
        q   = round(q,   4)
      ) |>
      rename(
        "Período"      = periodo,
        "Brecha (y−yₑ)" = gap,
        "Inflación π (%)" = pi,
        "Tasa r (%)"   = r,
        "Tipo cambio q" = q
      )
  },
  striped   = TRUE,
  hover     = TRUE,
  bordered  = TRUE,
  width     = "100%",
  align     = "c"
  )
}

# =============================================================================
#  LANZAR LA APP
# =============================================================================
shinyApp(ui = ui, server = server)

# =============================================================================
#  FIN DE app.R
# =============================================================================
