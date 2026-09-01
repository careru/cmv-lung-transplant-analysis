# =============================================================================
# ESTUDIO 1. Genotipado de las glucoproteinas gB y gH del citomegalovirus
# =============================================================================
#
# Tesis doctoral "Aspectos virologicos de la infeccion por citomegalovirus en
# receptores de trasplante pulmonar". Universidad Autonoma de Madrid, 2026.
#
# Cohorte historica cerrada de receptores de trasplante pulmonar del Hospital
# Universitario 12 de Octubre trasplantados entre 2009 y 2021, con seguimiento
# hasta el 31 de diciembre de 2022. Genotipado de gB (UL55) y gH (UL75) por PCR
# en tiempo real especifica de genotipo en muestras con carga viral de CMV
# >= 1.000 UI/mL. Analisis multinivel en tres unidades de observacion: la
# muestra (n = 483), el episodio de DNAemia (n = 234) y el paciente (n = 146).
#
# -----------------------------------------------------------------------------
# CORRESPONDENCIA CON LA TESIS
# -----------------------------------------------------------------------------
#
#  Codigo  Apartado de la tesis                                    Salida
#  ------  -----------------------------------------------------   ------------
#  [1]     5.1  Proceso de seleccion de la cohorte                  Figura 15
#  [2]     5.2  Rendimiento del genotipado                          Figura 16
#  [3]     5.3  Caracteristicas de la cohorte de estudio            Tabla 7
#  [4a]    5.4.1 Distribucion temporal de las muestras              Figuras 17-18
#  [4b]    5.4.2 Frecuencias genotipicas a nivel de muestra         Figuras 19-21
#  [4c]    5.4.2.4 Asociacion gB y gH                               Figura 22
#  [4d]    5.4.3 Relacion de los genotipos con la carga viral       Figuras 23-24
#  [4e]    5.4.4 Relacion de los genotipos con el Ct                Figura 25
#  [5a]    5.5.1 Caracteristicas de los episodios                   Tabla 8
#  [5b]    5.5.2 Frecuencias genotipicas por episodio               Tablas 9-11
#  [5c]    5.5.2 Genotipos mixtos por perfil serologico             Figura 26
#  [5d]    5.5.3 Cinetica viral segun genotipo individual           Tabla 12, Fig. 27
#  [5e]    5.5.3 Modelo multivariante en genotipos unicos de gB     Tabla 13
#  [5f]    5.5.4 Cinetica viral segun genotipos mixtos              Tabla 14, Fig. 28
#  [5g]    5.5.4 Modelo multivariante en genotipos mixtos           Tabla 15
#  [5h]    5.5.4 Efecto de gB3 en contexto de genotipos mixtos      Figura 29
#  [6a]    5.6.1 Variabilidad intraepisodio                         Tabla 16
#  [6b]    5.6.1 Distancia al P75 del Ct                            Figura 31
#  [6c]    5.6.2 Variabilidad interepisodio                         Tabla 17, Fig. 32
#  [7a]    5.7.1 Caracteristicas basales segun genotipos mixtos     Tabla 19
#  [7b]    5.7.2 Enfermedad sintomatica por CMV                     Tabla 20
#  [7c]    5.7.2 Tiempo hasta el primer episodio de DNAemia         Figura 33
#  [7d]    5.7.2 Numero de episodios por paciente                   Figura 34
#  [7e]    5.7.2 Tasa de episodios y exposicion viral acumulada     Figuras 35-36
#  [7f]    5.7.2 Recurrencia de DNAemia                             Tabla 21, Fig. 37
#  [7g]    5.7.3 Infecciones oportunistas no relacionadas con CMV   Figura 38
#  [7h]    5.7.3 Infecciones por virus respiratorios                Figura 39
#  [7i]    5.7.4.1 Rechazo agudo                                    Figura 40
#  [7j]    5.7.4.2 Disfuncion cronica del injerto pulmonar          (texto)
#  [7k]    5.7.5 Neoplasias y tumores solidos no cutaneos           Tabla 22, Fig. 41-42
#  [7l]    5.7.6 Mortalidad                                         Figura 43
#
# No se incluye el codigo de las figuras de elaboracion propia (esquemas,
# diagramas conceptuales e imagenes de la introduccion y de material y metodos),
# ni el de las Tablas 1 a 6, 18, 23 y 26, elaboradas manualmente.
#
# -----------------------------------------------------------------------------
# METODOS ESTADISTICOS
# -----------------------------------------------------------------------------
#
#  Variables continuas ....... mediana y rango intercuartilico; comparacion
#                              mediante U de Mann-Whitney o Kruskal-Wallis, con
#                              post hoc de Dunn y correccion de Bonferroni
#  Variables categoricas ..... n (%); comparacion mediante test exacto de Fisher
#  Asociacion gB-gH .......... chi-cuadrado de Pearson, V de Cramer y prueba
#                              binomial exacta por genotipo
#  Tiempo de aclaramiento .... regresion lineal multiple; diagnostico de
#                              colinealidad mediante factor de inflacion de la
#                              varianza (VIF)
#  Tasa de episodios ......... regresion binomial negativa con offset de
#                              persona-anos (IRR)
#  Desenlaces de tiempo a
#  evento .................... Kaplan-Meier con log-rank y modelos de Cox
#                              ajustados por edad al trasplante y perfil de
#                              riesgo serologico
#
#  Nivel de significacion: p < 0,05 bilateral.
#
# -----------------------------------------------------------------------------
# DATOS
# -----------------------------------------------------------------------------
#
# Este script no incluye datos de pacientes. Las rutas se han sustituido por
# rutas relativas del repositorio. La estructura esperada de cada archivo de
# entrada se describe en data/README.md.
#
# Las secciones son plegables en RStudio (## ---- [x] ---- ). Cada seccion es
# autocontenida una vez ejecutado el bloque [0] de carga y depuracion.
# =============================================================================


## ---- [0] Carga de librerias ----


## ---- [0a] Librerias ----
# Paquetes utilizados en todo el script.

library(readxl)
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(scales)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(purrr)
library(flextable)
library(survival)
library(survminer)
library(ggstatsplot)
library(showtext)
library(ggpubr)
library(rstatix)


## ---- [0b] Importacion de las muestras genotipadas ----
# Apartado 5.2. Lectura de las 618 muestras procesadas.

muestras <- read_excel("data/genotipos_muestras.xls",
                       col_types = c(
                         rep("numeric", 2), "date", "text", "numeric", "text",
                         rep("numeric", 14), "text", rep("skip", 3)))

length(unique(muestras$Nreg)) #692 muestras de 177 pacientes
length(unique(muestras$NHC))

muestras <- muestras %>% #quito 74 que no se encontraron / volumen insuf. (NE=.)
  filter(NE != ".")

head(muestras)
length(unique(muestras$Nreg)) #618 muestras de 169 pacientes
length(unique(muestras$NHC))

#Borrar entradas vacias  
muestras <- muestras %>% 
  filter(!is.na(Freg))


## ---- [0c] Clasificacion del resultado de genotipado ----
# Apartado 5.2. Genotipado completo, parcial o negativo para gB y gH.

library(dplyr)
muestras %>%
  mutate(
    gb_sum = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh_sum = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    categoria = case_when(
      gb_sum > 0 & gh_sum > 0 ~ "Genotipado completo",
      gb_sum > 0 & gh_sum == 0 ~ "Parcial - Solo gB",
      gh_sum > 0 & gb_sum == 0 ~ "Parcial - solo gH",
      gb_sum == 0 & gh_sum == 0 ~ "Negativo gB y gH" )) %>%
  count(categoria) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1))


## ---- [0d] Seleccion de las muestras con genotipado completo ----
# Apartado 5.2. Cohorte final de 483 muestras y 146 pacientes.

muestras_sel <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    )
  ) %>%
  filter(gt_cat == "Genotipado completo")

#Número de pacientes
n_distinct(muestras_sel$NHC)

# Número de muestras por paciente mínimo y maximo
muestras_por_paciente <- muestras_sel %>%
  count(NHC, name = "n_muestras")

min(muestras_por_paciente$n_muestras)
max(muestras_por_paciente$n_muestras)

# Mediana + IQR
quantile(muestras_por_paciente$n_muestras,
         probs = c(0.25, 0.5, 0.75),
         na.rm = TRUE)


## ---- [1] Figura 15. Categoria de DNAemia segun el perfil serologico D/R ----
# Apartado 5.1. Incluye las comparaciones pareadas entre perfiles serologicos.

library(ggplot2)
library(dplyr)

datos_ser <- data.frame(
  seroestatus = rep(c("D+/R+", "D-/R+", "D+/R-", "D-/R-"), each = 3),
  grupo = rep(c("Siempre CV negativa",
                "CV positiva <1000 UI/mL",
                "CV positiva ≥1000 UI/mL"), 4),
  n = c(
    # D+/R+: neg, <1000, ≥1000
    84, 46, 135,
    # D-/R+: neg, <1000, ≥1000
    27,  5,  18,
    # D+/R-: neg, <1000, ≥1000
    4,  3,  22,
    # D-/R-: neg, <1000, ≥1000
    9,  1,   1
  )
) %>%
  group_by(seroestatus) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100
  ) %>%
  ungroup() %>%
  mutate(
    seroestatus = factor(seroestatus, levels = c("D+/R+", "D-/R+", "D+/R-", "D-/R-")),
    grupo = factor(grupo, levels = c(
      "Siempre CV negativa",
      "CV positiva <1000 UI/mL",
      "CV positiva ≥1000 UI/mL"))
  )

# Comparaciones pareadas: proporción de CV ≥1000 UI/mL entre seroestatus
sero_high <- datos_ser %>%
  filter(grupo == "CV positiva ≥1000 UI/mL") %>%
  mutate(no = total - n) %>%
  select(seroestatus, si = n, no, total)

calc_p <- function(g1, g2) {
  d1 <- sero_high[sero_high$seroestatus == g1, ]
  d2 <- sero_high[sero_high$seroestatus == g2, ]
  m <- matrix(c(d1$si, d1$no, d2$si, d2$no), nrow = 2)
  fisher.test(m)$p.value
}

pares <- list(
  c("D+/R+", "D+/R-"),
  c("D-/R+", "D+/R-"),
  c("D+/R-", "D-/R-")
)

pvals <- sapply(pares, function(p) calc_p(p[1], p[2]))
print(data.frame(
  comparacion = sapply(pares, function(p) paste(p[1], "vs", p[2])),
  p_valor = round(pvals, 4)
))

stars_fn <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else NA_character_
}

# Filtrar solo las parejas con diferencias significativas
sig <- !is.na(sapply(pvals, stars_fn))
sig_pares <- pares[sig]
sig_stars <- sapply(pvals[sig], stars_fn)

# Mapeo seroestatus -> posición exacta de la barra ≥1000 UI/mL en el eje X
# Con position_dodge(width = 0.8) y 3 categorías, la 3ª barra (≥1000) queda
# desplazada +0.8/3 ≈ 0.267 respecto al centro del grupo
dodge_offset <- 0.8 / 3
x_map <- setNames((1:4) + dodge_offset, c("D+/R+", "D-/R+", "D+/R-", "D-/R-"))
xmin_v <- sapply(sig_pares, function(p) x_map[p[1]])
xmax_v <- sapply(sig_pares, function(p) x_map[p[2]])

# Apilar los corchetes a alturas crecientes (los más cortos abajo)
orden <- order(xmax_v - xmin_v, xmin_v)
xmin_v <- xmin_v[orden]
xmax_v <- xmax_v[orden]
sig_stars <- sig_stars[orden]

y_base <- 90
y_step <- 7
y_pos <- y_base + (seq_along(sig_stars) - 1) * y_step

# Etiquetas de riesgo bajo cada seroestatus
risk_labels <- c(
  "D+/R+" = "D+/R+\n(riesgo medio D+)\nn=265",
  "D-/R+" = "D-/R+\n(riesgo medio D-)\nn=50",
  "D+/R-" = "D+/R-\n(riesgo alto)\nn=29",
  "D-/R-" = "D-/R-\n(riesgo bajo)\nn=11"
)

# Construcción de los corchetes con annotate (sin dependencias extra
brackets <- data.frame(
  x    = xmin_v,
  xend = xmax_v,
  y    = y_pos,
  star = sig_stars
)

ggplot(datos_ser, aes(x = seroestatus, y = pct, fill = grupo)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = paste0(round(pct), "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.4, size = 4, fontface = "bold", color = "grey20"
  ) +
  # Línea horizontal del corchete
  geom_segment(
    data = brackets,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  # Patas verticales del corchete
  geom_segment(
    data = brackets,
    aes(x = x, xend = x, y = y - 1.2, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  geom_segment(
    data = brackets,
    aes(x = xend, xend = xend, y = y - 1.2, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  # Asteriscos
  geom_text(
    data = brackets,
    aes(x = (x + xend) / 2, y = y + 1.5, label = star),
    inherit.aes = FALSE,
    size = 5.5, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Siempre CV negativa"     = "#95A5A6",
      "CV positiva <1000 UI/mL" = "#5DADE2",
      "CV positiva ≥1000 UI/mL" = "#C0392B"
    ),
    name = ""
  ) +
  scale_x_discrete(labels = risk_labels) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 115),
    breaks = seq(0, 100, 25)
  ) +
  labs(x = "", y = "% pacientes") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(12, "pt"),
    legend.spacing.x = unit(8, "pt"),
    axis.text.x      = element_text(size = 12, lineheight = 1),
    axis.text.y      = element_text(size = 13),
    axis.title       = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    plot.margin      = margin(10, 20, 10, 10)
  )


## ---- [2] Figura 16. Distribucion anual del genotipado por tecnica de cuantificacion ----
# Apartado 5.2. Proporcion de genotipado completo por ano y plataforma.

library(patchwork)

df <- muestras %>%
  mutate(
    anio   = factor(format(as.Date(Freg), "%Y")),
    gb     = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh     = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    ),
    gt_cat = factor(
      gt_cat,
      levels = c("Genotipado completo",
                 "Genotipado parcial",
                 "Genotipado negativo")))

totales <- df %>%
  count(anio, gt_cat, name = "n") %>%
  group_by(anio) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  filter(gt_cat == "Genotipado completo") %>%
  mutate(
    pct     = n / total * 100,
    pct_fmt = paste0(formatC(pct, format = "f", digits = 0, decimal.mark = ","), "%")
  ) %>%
  mutate(anio = factor(anio, levels = levels(df$anio)))

# ── Gráfico principal ──────────────────────────────────────────────────────────
graph_temporal <- ggplot(df, aes(x = anio, fill = gt_cat)) +
  geom_bar(width = 0.85) +
  geom_text(
    data = totales,
    aes(x = anio, y = total, label = pct_fmt),
    inherit.aes = FALSE,
    vjust = -0.35, size = 4.5, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#4F94CD",
    "Genotipado parcial"  = "tan2",
    "Genotipado negativo" = "#CDC9C9"
  ), name = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "", y = "Número de muestras") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title.y     = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(11, "pt"),
    plot.margin      = margin(t = 10, r = 10, b = 0, l = 10))

# ── Franja de plataformas ──────────────────────────────────────────────────────
plataformas <- tibble::tribble(
  ~label,            ~xmin, ~xmax, ~color,
  "ARGENE\nSANGRE",   0.55,   4.45,  "#FF6B6B",
  "ALTONA\nSANGRE",   4.55,   6.45,  "#FF6B6B",
  "VERIS\nPLASMA",    6.55,   9.45,  "#FFD700",
  "ALTONA\nPLASMA",   9.55,  12.45,  "#FFD700"
) %>% mutate(xmid = (xmin + xmax) / 2)

franja <- ggplot(plataformas) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = color),
    alpha = 0.35, color = plataformas$color, linewidth = 0.5
  ) +
  geom_text(
    aes(x = xmid, y = 0.5, label = label),
    size = 3.8, fontface = "bold", lineheight = 0.9, color = "black"
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.5, 12.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),      expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 10, b = 5, l = 10))

# ── Combinar ───────────────────────────────────────────────────────────────────
graph_temporal / franja + plot_layout(heights = c(10, 1.2))

tabla_genotipos <- muestras %>%
  mutate(
    anio = factor(format(as.Date(Freg), "%Y")),
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipo completo",
      (gb >= 1 | gh >= 1) & !(gb >= 1 & gh >= 1) ~ "Genotipo parcial",
      TRUE ~ "Genotipo negativo"
    )
  ) %>%
  group_by(anio) %>%
  summarise(
    total      = n(),
    negativos  = sum(gt_cat == "Genotipo negativo"),
    parciales  = sum(gt_cat == "Genotipo parcial"),
    completos  = sum(gt_cat == "Genotipo completo"),
    perc_neg   = round(100 * negativos / total, 1),
    perc_par   = round(100 * parciales / total, 1),
    perc_comp  = round(100 * completos / total, 1),
    .groups    = "drop"
  )


tabla_genotipos


## ---- [3] Tabla 7. Caracteristicas basales de la cohorte ----
# Apartado 5.3. Comparacion por perfil serologico D/R.

library(readxl)
library(dplyr)
library(flextable)

pacientes <- read_excel("data/base_clinica.xlsx",
                        sheet = "pacientes_sel")

pacientes <- pacientes %>%
  mutate(
    edad_receptor           = as.numeric(difftime(Ftx, Fnac, units = "days")) / 365.25,
    tiempo_seguimiento_anos = if_else(
      exitus == 1,
      as.numeric(difftime(Fexitus, Ftx, units = "days")),
      as.numeric(difftime(as.Date("2022-12-31"), Ftx, units = "days"))
    ) / 365.25,
    serDR = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    tac   = grepl("TACRO", IS_base),
    cicl  = grepl("CICLO", IS_base),
    mmf   = grepl("MICO",  IS_base),
    aza   = grepl("AZA",   IS_base)
  )

gl <- pacientes
dm <- pacientes %>% filter(serDR == "D-/R+")
pp <- pacientes %>% filter(serDR == "D+/R+")
pm <- pacientes %>% filter(serDR == "D+/R-")

# ── Helpers ────────────────────────────────────────────────────────────────────
med <- function(x) sprintf("%.1f ± %.1f", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
mdn <- function(x) sprintf("%.1f (%.1f–%.1f)",
                           median(x, na.rm=TRUE),
                           quantile(x, 0.25, na.rm=TRUE),
                           quantile(x, 0.75, na.rm=TRUE))
pct <- function(x, tot=length(x)) sprintf("%d (%.1f%%)",
                                          sum(x, na.rm=TRUE),
                                          100*sum(x, na.rm=TRUE)/tot)
fp  <- function(p) if(is.na(p)) "-" else if(p<0.001) "<0,001" else gsub("\\.", ",", sprintf("%.3f", p))

eb  <- function(g, cat) {
  n <- sum(g$ebase == cat, na.rm=TRUE)
  sprintf("%d (%.1f%%)", n, 100*n/nrow(g))
}

# ── p-valores ──────────────────────────────────────────────────────────────────
p_edad_rx  <- fp(summary(aov(edad_receptor           ~ serDR, data=pacientes))[[1]][["Pr(>F)"]][1])
p_edad_don <- fp(summary(aov(edad_don                ~ serDR, data=pacientes))[[1]][["Pr(>F)"]][1])
p_sex_rx   <- fp(chisq.test(table(pacientes$SexF,     pacientes$serDR))$p.value)
p_sex_don  <- fp(chisq.test(table(pacientes$SexF_don, pacientes$serDR))$p.value)
p_tipotx   <- fp(chisq.test(table(pacientes$tipotx,   pacientes$serDR))$p.value)
p_ebase    <- fp(fisher.test(table(pacientes$ebase,   pacientes$serDR),
                             simulate.p.value=TRUE, B=10000)$p.value)
p_seg      <- fp(kruskal.test(tiempo_seguimiento_anos ~ serDR, data=pacientes)$p.value)

p_tac  <- tryCatch(fp(chisq.test(table(pacientes$tac,  pacientes$serDR), simulate.p.value=TRUE)$p.value), error=function(e) "-")
p_cicl <- tryCatch(fp(fisher.test(table(pacientes$cicl, pacientes$serDR))$p.value), error=function(e) "-")
p_mmf  <- tryCatch(fp(chisq.test(table(pacientes$mmf,  pacientes$serDR), simulate.p.value=TRUE)$p.value), error=function(e) "-")
p_aza  <- tryCatch(fp(fisher.test(table(pacientes$aza,  pacientes$serDR))$p.value), error=function(e) "-")

# ── Tabla ──────────────────────────────────────────────────────────────────────
df <- data.frame(
  Caracteristica = c(
    "Edad del receptor (años), media ± DE",
    "Edad del donante (años), media ± DE",
    "Género masculino del receptor, n (%)",
    "Género masculino del donante, n (%)",
    "Tipo de trasplante, n (%)",
    "   Unilateral",
    "   Bilateral",
    "Enfermedad basal, n (%)",
    "   EPOC / enfisema",
    "   EPID",
    "   Fibrosis quística",
    "   Hipertensión pulmonar",
    "   Otros",
    "Inmunosupresión inicial, n (%)",
    "   Esteroides",
    "   Tacrolimus",
    "   Ciclosporina",
    "   Micofenolato mofetilo / ácido micofenólico",
    "   Azatioprina",
    "Tiempo de seguimiento (años), mediana (RIC)"
  ),
  Global = c(
    med(gl$edad_receptor), med(gl$edad_don),
    pct(gl$SexF==0), pct(gl$SexF_don==0),
    "", pct(gl$tipotx=="UNI"), pct(gl$tipotx=="BI"),
    "", eb(gl,"EPOC"), eb(gl,"EPID"), eb(gl,"FQ"), eb(gl,"HAP"), eb(gl,"OTROS"),
    "", pct(rep(TRUE,nrow(gl))), pct(gl$tac), pct(gl$cicl), pct(gl$mmf), pct(gl$aza),
    mdn(gl$tiempo_seguimiento_anos)
  ),
  DmRp = c(
    med(dm$edad_receptor), med(dm$edad_don),
    pct(dm$SexF==0, nrow(dm)), pct(dm$SexF_don==0, nrow(dm)),
    "", pct(dm$tipotx=="UNI", nrow(dm)), pct(dm$tipotx=="BI", nrow(dm)),
    "", eb(dm,"EPOC"), eb(dm,"EPID"), eb(dm,"FQ"), eb(dm,"HAP"), eb(dm,"OTROS"),
    "", pct(rep(TRUE,nrow(dm)),nrow(dm)), pct(dm$tac,nrow(dm)),
    pct(dm$cicl,nrow(dm)), pct(dm$mmf,nrow(dm)), pct(dm$aza,nrow(dm)),
    mdn(dm$tiempo_seguimiento_anos)
  ),
  DpRp = c(
    med(pp$edad_receptor), med(pp$edad_don),
    pct(pp$SexF==0, nrow(pp)), pct(pp$SexF_don==0, nrow(pp)),
    "", pct(pp$tipotx=="UNI", nrow(pp)), pct(pp$tipotx=="BI", nrow(pp)),
    "", eb(pp,"EPOC"), eb(pp,"EPID"), eb(pp,"FQ"), eb(pp,"HAP"), eb(pp,"OTROS"),
    "", pct(rep(TRUE,nrow(pp)),nrow(pp)), pct(pp$tac,nrow(pp)),
    pct(pp$cicl,nrow(pp)), pct(pp$mmf,nrow(pp)), pct(pp$aza,nrow(pp)),
    mdn(pp$tiempo_seguimiento_anos)
  ),
  DpRm = c(
    med(pm$edad_receptor), med(pm$edad_don),
    pct(pm$SexF==0, nrow(pm)), pct(pm$SexF_don==0, nrow(pm)),
    "", pct(pm$tipotx=="UNI", nrow(pm)), pct(pm$tipotx=="BI", nrow(pm)),
    "", eb(pm,"EPOC"), eb(pm,"EPID"), eb(pm,"FQ"), eb(pm,"HAP"), eb(pm,"OTROS"),
    "", pct(rep(TRUE,nrow(pm)),nrow(pm)), pct(pm$tac,nrow(pm)),
    pct(pm$cicl,nrow(pm)), pct(pm$mmf,nrow(pm)), pct(pm$aza,nrow(pm)),
    mdn(pm$tiempo_seguimiento_anos)
  ),
  p = c(
    p_edad_rx, p_edad_don, p_sex_rx, p_sex_don,
    p_tipotx, "", "",
    p_ebase, "", "", "", "", "",
    "", "-", p_tac, p_cicl, p_mmf, p_aza,
    p_seg
  ),
  stringsAsFactors = FALSE
)

colnames(df) <- c("Característica", "Global\n(n=146)", "D-/R+\n(n=15)",
                  "D+/R+\n(n=110)", "D+/R-\n(n=21)", "p")

# ── Flextable ──────────────────────────────────────────────────────────────────
apartados <- c("Tipo de trasplante, n (%)", "Enfermedad basal, n (%)",
               "Inmunosupresión inicial, n (%)")

flextable(df) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = ~ Característica %in% apartados, part = "body") %>%
  italic(i = ~ grepl("^   ", Característica), j = "Característica", part = "body") %>%
  align(j = "Característica", align = "left",   part = "all") %>%
  align(j = -1,               align = "center", part = "all") %>%
  fontsize(part = "all", size = 10) %>%
  padding(part = "all", padding = 3) %>%
  add_footer_lines("*p < 0,05. EPID: enfermedad pulmonar intersticial difusa; EPOC: enfermedad pulmonar obstructiva crónica. Variables categóricas comparadas mediante test exacto de Fisher con simulación de Monte Carlo (B=10.000). Cada fármaco inmunosupresor se analizó como variable dicotómica independiente.") %>%
  fontsize(part = "footer", size = 9) %>%
  autofit()


## ---- [4a1] Figura 17. Distribucion anual de muestras y numero acumulado de trasplantes ----
# Apartado 5.4.1.

pacientes <- read_excel("data/base_clinica.xlsx",
                        sheet = "pacientes_sel")
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

## 1) Muestras completamente genotipadas por año
df_comp <- muestras %>%
  mutate(
    anio = year(as.Date(Freg)),
    gb   = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh   = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    )
  ) %>%
  filter(gt_cat == "Genotipado completo") %>%
  count(anio, name = "n_muestras")

## 2) Trasplantes por año y acumulados a partir de `pacientes$Ftx`
tx_por_anio <- pacientes %>%
  mutate(anio = year(as_date(Ftx))) %>%
  count(anio, name = "n_tx") %>%
  arrange(anio) %>%
  mutate(tx_acum = cumsum(n_tx))

## 3) Unir tablas, asegurar años 2008–2022 y asignar metodología
anios_levels <- 2008:2022

metodologia_levels <- c(
  "ARGENE SANGRE (2008-2013)",
  "ALTONA SANGRE (2014-2016)",
  "VERIS PLASMA (2017-2019)",
  "ALTONA PLASMA (2020-2022)"
)

plot_df <- full_join(df_comp, tx_por_anio, by = "anio") %>%
  complete(anio = anios_levels) %>%
  arrange(anio) %>%
  mutate(
    anio_factor = factor(anio, levels = anios_levels),
    metodologia = case_when(
      anio %in% 2008:2013 ~ metodologia_levels[1],
      anio %in% 2014:2016 ~ metodologia_levels[2],
      anio %in% 2017:2019 ~ metodologia_levels[3],
      anio %in% 2020:2022 ~ metodologia_levels[4]
    ),
    metodologia = factor(metodologia, levels = metodologia_levels)
  )

## 4) Paleta: sangre en rojos, plasma en naranjas y amarillos
colores_metodologia <- c(
  "ARGENE SANGRE (2008-2013)" = "#8B3A3A",  # rojo granate
  "ALTONA SANGRE (2014-2016)" = "lightsalmon",  # rojo medio
  
  "VERIS PLASMA (2017-2019)"  = "darkorange3",  # naranja intenso
  "ALTONA PLASMA (2020-2022)" = "#FDD049"   # amarillo dorado
)

## 5) Ratio para el eje secundario
ratio <- max(plot_df$n_muestras, na.rm = TRUE) /
  max(plot_df$tx_acum, na.rm = TRUE)

## 6) Gráfico final
graph_acumulados <- ggplot(plot_df, aes(x = anio_factor)) +
  geom_col(aes(y = n_muestras, fill = metodologia),
           width = 0.8, na.rm = TRUE) +
  geom_text(aes(y = 0, label = n_muestras),
            vjust = 1.4, size = 4, fontface = "bold", na.rm = TRUE) +
  geom_line(aes(y = tx_acum * ratio, group = 1),
            linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(y = tx_acum * ratio),
             size = 2, na.rm = TRUE) +
  geom_text(aes(y = tx_acum * ratio, label = tx_acum),
            vjust = -1, size = 4, na.rm = TRUE) +
  scale_fill_manual(name = "Metodología:", values = colores_metodologia) +
  scale_y_continuous(
    name = "Número de muestras",
    expand = expansion(mult = c(0.05, 0.12)),
    sec.axis = sec_axis(~ . / ratio,
                        name = "Número de RTP (acumulado)")
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )
graph_acumulados


## ---- [4a2] Figura 18. Muestras con DNAemia cuantificable segun categoria de carga viral ----
# Apartado 5.4.1.

# Preparar datos para el gráfico (en %)
df_cv_cat <- DTA_CV_TXP_cv2 %>%
  group_by(anio, cat_cv) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(anio) %>%
  mutate(
    total_anio = sum(n),
    pct = round(100 * n / total_anio, 0)
  ) %>%
  ungroup()

pos_dodge <- position_dodge(width = 0.85)

graph_cv_cat <- ggplot(df_cv_cat,
                       aes(x = anio, y = pct, fill = cat_cv)) +
  geom_col(position = pos_dodge, width = 0.8) +
  
  # Etiquetas SOLO para CV ≥ 1.000 UI/mL (las rojas oscuras, a la derecha)
  geom_text(
    aes(
      label = ifelse(cat_cv == "CV ≥ 1.000 UI/mL",
                     paste0(round(pct), "%"),
                     "")
    ),
    position = pos_dodge,
    vjust = -0.5,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "CV < 1.000 UI/mL" = "#FFAEB9",
      "CV ≥ 1.000 UI/mL" = "#8B0000"
    ),
    name = "",
    labels = c(
      "CV < 1.000 UI/mL",
      "CV ≥ 1.000 UI/mL"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    limits = c(0, 100),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "",
    y = "% muestras CV detectable"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.text     = element_text(size = 12),
    axis.text       = element_text(size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
    axis.title      = element_text(size = 14),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_cat


## ---- [4b1] Frecuencias genotipicas a nivel de muestra: tablas resumen ----
# Apartado 5.4.2. Recuentos de gB, gH y de la combinacion gB/gH.

df_gb <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~replace_na(., 0))),  # suma de GB1-GB4
    combo  = apply(select(., GB1:GB4) == 1, 1, function(r)
      paste(colnames(select(., GB1:GB4))[which(r)], collapse = "+"))
  ) %>%
  filter(gb_sum %in% 1:3)

tabla_resumen <- df_gb %>%
  count(gb_sum, name = "n") %>%
  mutate(
    tipo = factor(paste(gb_sum, "genotipo(s)"),
                  levels = c("1 genotipo(s)","2 genotipo(s)","3 genotipo(s)")),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(tipo) %>%
  select(tipo, n, porcentaje)

tabla_resumen

desglose <- df_gb %>%
  count(gb_sum, combo, name = "n") %>%
  group_by(gb_sum) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  arrange(gb_sum, desc(n)) %>%
  ungroup()

desglose_1 <- desglose %>% filter(gb_sum == 1) %>% select(combo, n, porcentaje)
desglose_2 <- desglose %>% filter(gb_sum == 2) %>% select(combo, n, porcentaje)
desglose_3 <- desglose %>% filter(gb_sum == 3) %>% select(combo, n, porcentaje)

desglose_1  # combinaciones cuando hay 1 genotipo
desglose_2  # combinaciones cuando hay 2 genotipos
desglose_3  # combinaciones cuando hay 3 genotipos


df_gh <- muestras_sel %>%
  mutate(
    gh_sum = rowSums(across(GH1:GH2, ~replace_na(., 0))),  # suma GH1–GH2
    combo  = apply(select(., GH1:GH2) == 1, 1, function(r)
      paste(colnames(select(., GH1:GH2))[which(r)], collapse = "+"))
  ) %>%
  filter(gh_sum %in% 1:2)

tabla_resumen_gh <- df_gh %>%
  count(gh_sum, name = "n") %>%
  mutate(
    tipo = factor(paste(gh_sum, "genotipo(s)"),
                  levels = c("1 genotipo(s)", "2 genotipo(s)")),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(tipo) %>%
  select(tipo, n, porcentaje)

tabla_resumen_gh

desglose_gh <- df_gh %>%
  count(gh_sum, combo, name = "n") %>%
  group_by(gh_sum) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  arrange(gh_sum, desc(n)) %>%
  ungroup()

desglose_1_gh <- desglose_gh %>%
  filter(gh_sum == 1) %>%
  select(combo, n, porcentaje)

desglose_2_gh <- desglose_gh %>%
  filter(gh_sum == 2) %>%
  select(combo, n, porcentaje)

desglose_1_gh   # combinaciones cuando hay 1 genotipo gH
desglose_2_gh   # combinaciones cuando hay 2 genotipos gH


library(dplyr)
library(tidyr)
library(stringr)

# Partimos de `muestras_sel` (solo muestras con genotipado completo)

df_gbgh <- muestras_sel %>%
  mutate(
    # nº de genotipos presentes
    gb_sum = rowSums(across(GB1:GB4, ~ tidyr::replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~ tidyr::replace_na(., 0))),
    
    # combinaciones gB (formato gB1, gB1 + gB3, etc.)
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    
    # combinaciones gH (formato gH1, gH1 + gH2)
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    
    # tipo único/mixto por proteína
    gB_tipo = if_else(gb_sum == 1, "gBunico", "gBmixto"),
    gH_tipo = if_else(gh_sum == 1, "gHunico", "gHmixto"),
    
    # grupo combinado
    grupo = paste(gB_tipo, gH_tipo, sep = "_"),
    
    # combinación conjunta gB/gH
    combo_total = paste(gB_combo, gH_combo, sep = " / ")
  ) %>%
  # por seguridad, nos quedamos solo con filas con gB y gH definidos
  filter(!is.na(gB_combo), !is.na(gH_combo))

## 1) Tabla resumen de grupos (gBunico_gHunico, etc.)

niveles_grupo <- c(
  "gBunico_gHunico",
  "gBmixto_gHunico",
  "gBunico_gHmixto",
  "gBmixto_gHmixto"
)

tabla_grupos <- df_gbgh %>%
  count(grupo, name = "n") %>%
  mutate(
    grupo = factor(grupo, levels = niveles_grupo),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(grupo)

tabla_grupos

## Subdesglose gBunico_gHunico
desglose_gBunico_gHunico <- df_gbgh %>%
  filter(grupo == "gBunico_gHunico") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBunico_gHunico

## Subdesglose gBmixto_gHunico
desglose_gBmixto_gHunico <- df_gbgh %>%
  filter(grupo == "gBmixto_gHunico") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBmixto_gHunico

## Subdesglose gBunico_gHmixto
desglose_gBunico_gHmixto <- df_gbgh %>%
  filter(grupo == "gBunico_gHmixto") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBunico_gHmixto

## 5) Subdesglose gBmixto_gHmixto
desglose_gBmixto_gHmixto <- df_gbgh %>%
  filter(grupo == "gBmixto_gHmixto") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBmixto_gHmixto


## ---- [4b2] Figura 19. Numero de genotipos detectados por muestra ----
# Apartado 5.4.2. Diagramas de anillo para gB, gH y la combinacion.

library(dplyr)
library(tidyr)


## 1) Filtrar muestras con genotipado completo gB + gH
df_completo <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(c(GB1, GB2, GB3, GB4), ~ replace_na(., 0))),
    gh_sum = rowSums(across(c(GH1, GH2),           ~ replace_na(., 0)))
  ) %>%
  filter(gb_sum > 0, gh_sum > 0)

## 1b) Porcentajes de genotipos mixtos, calculados desde los datos
pct_mixtos <- df_completo %>%
  summarise(
    gb   = 100 * mean(gb_sum > 1),
    gh   = 100 * mean(gh_sum > 1),
    gbgh = 100 * mean(gb_sum > 1 | gh_sum > 1)
  )

fmt_mix <- function(x) formatC(x, format = "f", digits = 1, decimal.mark = ",")

cat(sprintf("gB mixtos: %s%% | gH mixtos: %s%% | gB/gH mixtos: %s%%\n",
            fmt_mix(pct_mixtos$gb), fmt_mix(pct_mixtos$gh), fmt_mix(pct_mixtos$gbgh)))

## 2) Resumen gB
tabla_gb <- df_completo %>%
  mutate(n_gen = rowSums(across(GB1:GB4, ~ replace_na(., 0)))) %>%
  filter(n_gen %in% 1:3) %>%
  count(n_gen, name = "n") %>%
  mutate(
    tipo = case_when(
      n_gen == 1 ~ "1 genotipo",
      n_gen == 2 ~ "2 genotipos",
      n_gen == 3 ~ "3 genotipos"
    ),
    porcentaje = round(100 * n / sum(n), 1),
    etiqueta = paste0(tipo, "\n", n, " (", gsub("\\.", ",", porcentaje), "%)")
  )
Prop_gb <- tabla_gb$n
cols_gb <- c("#63B8FF", "dodgerblue3", "#27408B")[seq_along(Prop_gb)]
n_total_gb <- sum(Prop_gb)

## 3) Resumen gH
tabla_gh <- df_completo %>%
  mutate(n_gen = rowSums(across(GH1:GH2, ~ replace_na(., 0)))) %>%
  filter(n_gen %in% 1:2) %>%
  count(n_gen, name = "n") %>%
  mutate(
    tipo = ifelse(n_gen == 1, "1 genotipo", "2 genotipos"),
    porcentaje = round(100 * n / sum(n), 1),
    etiqueta = paste0(tipo, "\n", n, " (", gsub("\\.", ",", porcentaje), "%)")
  )
Prop_gh <- tabla_gh$n
cols_gh <- c("#DDA0DD", "#CD69C9")[seq_along(Prop_gh)]
n_total_gh <- sum(Prop_gh)

## 4) Resumen combinado gB + gH
tabla_gbgh <- tabla_grupos %>%
  mutate(
    grupo_label = case_when(
      grupo == "gBunico_gHunico" ~ "Genotipos únicos gB / gH",
      grupo == "gBmixto_gHunico" ~ "Genotipos mixtos gB / únicos gH",
      grupo == "gBunico_gHmixto" ~ "Genotipos únicos gB / mixtos gH",
      grupo == "gBmixto_gHmixto" ~ "Genotipos mixtos gB / gH"
    ),
    color = case_when(
      grupo == "gBunico_gHunico" ~ "#3CB371",
      grupo == "gBmixto_gHunico" ~ "#F08080",
      grupo == "gBunico_gHmixto" ~ "#CD5C5C",
      grupo == "gBmixto_gHmixto" ~ "#8B0000"
    ),
    etiqueta = paste0(grupo_label, "\n", n, " (", gsub("\\.", ",", porcentaje), "%)")
  )

## 5) Funciones de dibujo con más aire
plot_donut <- function(props, labels, cols, center_text) {
  pie(
    props,
    labels = labels,
    col = cols,
    border = "white",
    cex = 1.4,
    radius = 0.9,
    init.angle = 45,
    main = ""
  )
  symbols(0, 0,
          circles = 0.5,
          inches = FALSE,
          add = TRUE,
          bg = "white",
          fg = "white")
  text(0, 0, center_text, cex = 1.5, font = 2)
}

plot_donut_gbgh <- function(props, labels, cols, center_text) {
  pie(
    props,
    labels = labels,
    col = cols,
    border = "white",
    cex = 1.4,
    radius = 0.9,
    init.angle = 45,
    main = ""
  )
  symbols(0, 0,
          circles = 0.5,
          inches = FALSE,
          add = TRUE,
          bg = "white",
          fg = "white")
  text(0, 0, center_text, cex = 1.5, font = 2)
}

## 5b) Etiqueta con recuadro de color bajo el donut
box_label <- function(txt, col) {
  x0 <- 0; y0 <- -1.28
  w  <- strwidth(txt,  cex = 1.6, font = 2) + 0.15
  h  <- strheight(txt, cex = 1.6, font = 2) + 0.16
  rect(x0 - w/2, y0 - h/2, x0 + w/2, y0 + h/2,
       col = col, border = NA, xpd = NA)
  text(x0, y0, txt, col = "white", font = 2, cex = 1.5, xpd = NA)
}

## 6) Layout compuesto
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(1, 1))
par(mar = c(5, 2, 3, 2))

# (A) Donut gB
plot_donut(Prop_gb, tabla_gb$etiqueta, cols_gb,
           paste0("gB\nn = ", n_total_gb))
title(main = "A", adj = 0, cex.main = 3)
box_label(paste0("gB mixtos: ", fmt_mix(pct_mixtos$gb), "%"), "#27408B")

# (B) Donut gH
plot_donut(Prop_gh, tabla_gh$etiqueta, cols_gh,
           paste0("gH\nn = ", n_total_gh))
title(main = "B", adj = 0, cex.main = 3)
box_label(paste0("gH mixtos: ", fmt_mix(pct_mixtos$gh), "%"), "#CD69C9")

# (C) Donut combinado gB + gH
plot_donut_gbgh(tabla_gbgh$n, tabla_gbgh$etiqueta, tabla_gbgh$color,
                paste0("gB + gH\nn = ", sum(tabla_gbgh$n)))
title(main = "C", adj = 0.02, cex.main = 3)
box_label(paste0("gB/gH mixtos: ", fmt_mix(pct_mixtos$gbgh), "%"), "#8B0000")

## 8) Guardar para media pagina A4 vertical
png("output/Figura_19_donuts_gB_gH.png",
    width = 16.5, height = 12, units = "cm", res = 300, bg = "white",
    pointsize = 10)

layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(1, 1.05))
par(mar = c(3.2, 1, 1.8, 1))

plot_donut(Prop_gb, tabla_gb$etiqueta, cols_gb,
           paste0("gB\nn = ", n_total_gb))
title(main = "A", adj = 0, cex.main = 1.6)
box_label(paste0("gB mixtos: ", fmt_mix(pct_mixtos$gb), "%"), "#27408B")

plot_donut(Prop_gh, tabla_gh$etiqueta, cols_gh,
           paste0("gH\nn = ", n_total_gh))
title(main = "B", adj = 0, cex.main = 1.6)
box_label(paste0("gH mixtos: ", fmt_mix(pct_mixtos$gh), "%"), "#CD69C9")

plot_donut_gbgh(tabla_gbgh$n, tabla_gbgh$etiqueta, tabla_gbgh$color,
                paste0("gB + gH\nn = ", sum(tabla_gbgh$n)))
title(main = "C", adj = 0.02, cex.main = 1.6)
box_label(paste0("gB/gH mixtos: ", fmt_mix(pct_mixtos$gbgh), "%"), "#8B0000")

dev.off()


## ---- [4b3] Figura 20. Frecuencia de los genotipos de gB y gH ----
# Apartado 5.4.2. Genotipos unicos de gB y gH y combinaciones mixtas de gB.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## Helper: porcentaje con coma decimal
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")

## 0) Bases comunes
df_gb <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    combo  = apply(select(., GB1:GB4) == 1, 1, function(r)
      paste0("gB", which(r), collapse = " + "))
  ) %>%
  filter(gb_sum %in% 1:3)

df_gh <- muestras_sel %>%
  mutate(
    gh_sum = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    combo  = apply(select(., GH1:GH2) == 1, 1, function(r)
      paste0("gH", which(r), collapse = " + "))
  ) %>%
  filter(gh_sum %in% 1:2)

## 1) Únicos gB (gB1-gB4)
tabla_uni_gb <- df_gb %>%
  filter(gb_sum == 1) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni_gb <- sum(tabla_uni_gb$n)

cols_uni_gb <- c("gB1" = "#C6DBEF", "gB2" = "#6BAED6",
                 "gB3" = "#3182BD", "gB4" = "#08519C")

## 2) Únicos gH (gH1-gH2)
tabla_uni_gh <- df_gh %>%
  filter(gh_sum == 1) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni_gh <- sum(tabla_uni_gh$n)

cols_uni_gh <- c("gH1" = "#DDA0DD", "gH2" = "#CD69C9")

## 3) Mixtos gB (todas las combinaciones, incl. n=1)
tabla_mix_gb <- df_gb %>%
  filter(gb_sum %in% 2:3) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_mix_gb <- sum(tabla_mix_gb$n)

cols_mix_gb <- grDevices::hcl.colors(nrow(tabla_mix_gb), "Dark 3")
names(cols_mix_gb) <- as.character(tabla_mix_gb$combo)

## 4) Panel únicos gB
p_uni_gb <- ggplot(tabla_uni_gb, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_uni_gb, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 1.30))) +
  labs(x = "", y = "",
       title = paste0("(A) Genotipos únicos gB (n = ", n_uni_gb, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 5) Panel únicos gH
p_uni_gh <- ggplot(tabla_uni_gh, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_uni_gh, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 150, 300),
    limits = c(0, 500),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "", y = "",
       title = paste0("(B) Genotipos únicos gH (n = ", n_uni_gh, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 6) Panel mixtos gB
p_mix_gb <- ggplot(tabla_mix_gb, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_mix_gb, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Número de muestras", y = "",
       title = paste0("(C) Genotipos mixtos gB (n = ", n_mix_gb, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 13,
                                margin = margin(b = 9, l = -35)),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 7) Composición: fila superior (gB | gH) y debajo mixtos gB a todo el ancho
panel_gb <- (p_uni_gb | p_uni_gh) / p_mix_gb +
  plot_layout(heights = c(max(nrow(tabla_uni_gb), nrow(tabla_uni_gh)),
                          nrow(tabla_mix_gb)))
panel_gb


## ---- [4b4] Figura 21. Combinaciones genotipicas de gB/gH ----
# Apartado 5.4.2.3. Combinaciones unicas y mixtas.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## Helper: porcentaje con coma decimal
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")

## 0) Base gB+gH
df_gbgh <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    # única solo si gB y gH son ambos únicos
    es_unica = (gb_sum == 1 & gh_sum == 1)
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

## 1) Combinaciones ÚNICAS (gB único / gH único)
tabla_uni <- df_gbgh %>%
  filter(es_unica) %>%
  count(combo_total, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo_total = factor(combo_total, levels = combo_total),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni <- sum(tabla_uni$n)

cols_uni <- colorRampPalette(c("#C7E9C0", "#00441B"))(nrow(tabla_uni))

## 2) Combinaciones MIXTAS (mezcla en gB, en gH o en ambas)
tabla_mix <- df_gbgh %>%
  filter(!es_unica) %>%
  count(combo_total, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo_total = factor(combo_total, levels = combo_total),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_mix <- sum(tabla_mix$n)

cols_mix <- colorRampPalette(c("#FCBBA1", "#67000D"))(nrow(tabla_mix))

## 3) Panel A: combinaciones únicas
p_uni <- ggplot(tabla_uni, aes(x = n, y = combo_total, fill = combo_total)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4) +
  scale_fill_manual(values = cols_uni, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = "", y = "",
       title = paste0("(A) Combinaciones únicas gB/gH (n = ", n_uni, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 12,
                                margin = margin(b = 7, l = -30)),  
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 4) Panel B: combinaciones mixtas
p_mix <- ggplot(tabla_mix, aes(x = n, y = combo_total, fill = combo_total)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4) +
  scale_fill_manual(values = cols_mix, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.50))) +
  labs(x = "Número de muestras", y = "",
       title = paste0("(B) Combinaciones mixtas gB/gH (n = ", n_mix, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 12,
                                margin = margin(b = 7, l = -30)),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 5) Composición: únicas arriba, mixtas abajo (alturas según nº de barras)
panel_gbgh <- p_uni / p_mix +
  plot_layout(heights = c(nrow(tabla_uni), nrow(tabla_mix)))
panel_gbgh


## ---- [4c] Figura 22. Asociacion entre los genotipos de gB y gH ----
# Apartado 5.4.2.4. Heatmap, chi-cuadrado y V de Cramer.

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

# Preparación de datos
df_asoc <- df_gbgh %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

# Tabla de contingencia (frecuencia de combinaciones gB ↔ gH)
tabla_asoc <- df_asoc %>%
  count(gB_combo, gH_combo, name = "n") %>%
  pivot_wider(names_from = gH_combo, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("gB_combo")

tabla_asoc

# Test estadístico de independencia (χ² y Cramér’s V)
m_asoc <- as.matrix(tabla_asoc)
test_chi <- suppressWarnings(chisq.test(m_asoc))

n <- sum(m_asoc)
phi2 <- test_chi$statistic / n
k <- min(nrow(m_asoc), ncol(m_asoc))
cramers_v <- sqrt(phi2 / (k - 1))

cat("=== Asociación gB ↔ gH ===\n")
cat("Chi-cuadrado:", round(test_chi$statistic, 3), "\n")
cat("Grados de libertad:", test_chi$parameter, "\n")
cat("p-valor:", round(test_chi$p.value, 4), "\n")
cat("Cramér’s V:", round(cramers_v, 3), "\n\n")

# Mapa de calor ordenado
tabla_larga <- df_asoc %>%
  count(gB_combo, gH_combo, name = "Frecuencia")

# Definir orden manual de los ejes
orden_gB <- c("gB1", "gB2", "gB3", "gB4",
              "gB1 + gB2", "gB1 + gB3", "gB1 + gB4",
              "gB2 + gB3", "gB2 + gB4", "gB3 + gB4",
              "gB1 + gB2 + gB3", "gB1 + gB2 + gB4", "gB1 + gB3 + gB4",
              "gB2 + gB3 + gB4")
orden_gH <- c("gH1", "gH2", "gH1 + gH2")

tabla_larga <- tabla_larga %>%
  mutate(
    gB_combo = factor(gB_combo, levels = orden_gB),
    gH_combo = factor(gH_combo, levels = orden_gH)
  )

ggplot(tabla_larga, aes(x = gH_combo, y = gB_combo, fill = Frecuencia)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Frecuencia), color = "black", size = 3.5) +
  scale_fill_gradient(low = "#E0F3DB", high = "#0868AC") +
  scale_y_discrete(limits = rev(levels(tabla_larga$gB_combo))) +  # ← invierte el eje Y
  labs(x = "", y = "", title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 11)
  )

# Interpretación automática
if (test_chi$p.value > 0.05) {
  cat("No se detectan asociaciones estadísticamente significativas entre los genotipos gB y gH (p > 0.05).\n")
} else {
  cat("Se detecta asociación significativa entre gB y gH (p ≤ 0.05).\n")
}

cat("Fuerza de asociación (Cramér’s V): ",
    ifelse(cramers_v < 0.1, "muy débil",
           ifelse(cramers_v < 0.3, "moderada", "fuerte")), "\n")


## ---- [4d1] Figura 23. Carga viral segun el resultado del genotipado ----
# Apartado 5.4.3. Genotipado completo, parcial y negativo.

df <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    ),
    gt_cat = factor(gt_cat, levels = c(
      "Genotipado completo",
      "Genotipado parcial",
      "Genotipado negativo"
    )),
    logCV = log10(CV)
  ) %>%
  drop_na(gt_cat, logCV)

comparaciones <- list(
  c("Genotipado completo", "Genotipado parcial"),
  c("Genotipado completo", "Genotipado negativo"),
  c("Genotipado parcial", "Genotipado negativo")
)

stats <- df %>%
  pairwise_wilcox_test(logCV ~ gt_cat, comparisons = comparaciones, p.adjust.method = "none") %>%
  add_significance("p") %>%
  add_y_position(fun = "max", step.increase = 0.15)

y_top <- max(stats$y.position, na.rm = TRUE) * 1.05

p <- ggplot(df, aes(x = gt_cat, y = logCV, fill = gt_cat)) +
  geom_boxplot(
    width = 0.50, color = "grey25", linewidth = 0.45,
    outlier.shape = 21, outlier.size = 1.8,
    outlier.fill = "white", outlier.color = "grey50",
    outlier.stroke = 0.4, outlier.alpha = 0.7
  ) +
  geom_jitter(
    aes(color = gt_cat),
    width = 0.12, size = 1.0, alpha = 0.35, show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#63B8FF",
    "Genotipado parcial"  = "#FFA54F",
    "Genotipado negativo" = "#BDBDBD"
  )) +
  scale_color_manual(values = c(
    "Genotipado completo" = "#1A6EBF",
    "Genotipado parcial"  = "#CC7A1A",
    "Genotipado negativo" = "#757575"
  )) +
  labs(x = NULL, y = "Carga viral CMV Log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.12)),
    limits = c(NA, y_top)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
    axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 11, color = "grey20", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    data = stats,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.012,
    bracket.size = 0.4,
    size = 4.0,
    color = "grey30",
    inherit.aes = FALSE
  )

p


## ---- [4d2] Figura 24. Carga viral segun genotipos unicos y mixtos ----
# Apartado 5.4.3. Cuatro paneles: gB, gH y combinacion.

## Base (solo genotipado completo gB y gH presentes)
dfc <- muestras_sel %>%
  mutate(
    gb_sum  = rowSums(across(GB1:GB4, ~ tidyr::replace_na(., 0))),
    gh_sum  = rowSums(across(GH1:GH2, ~ tidyr::replace_na(., 0))),
    CVlog10 = log10(CV)
  ) %>%
  filter(gb_sum >= 1, gh_sum >= 1) %>%
  drop_na(CV, CVlog10)

## Secciones
gb_ind <- dfc %>%
  filter(gb_sum == 1) %>%
  pivot_longer(GB1:GB4, names_to = "categoria_raw", values_to = "v") %>%
  filter(v == 1) %>%
  transmute(
    seccion   = "gB (individuales únicos)",
    categoria = paste0("g", substr(categoria_raw, 2, nchar(categoria_raw))),
    CV, CVlog10
  )

gh_ind <- dfc %>%
  filter(gh_sum == 1) %>%
  pivot_longer(GH1:GH2, names_to = "categoria_raw", values_to = "v") %>%
  filter(v == 1) %>%
  transmute(
    seccion   = "gH (individuales únicos)",
    categoria = paste0("g", substr(categoria_raw, 2, nchar(categoria_raw))),
    CV, CVlog10
  )

gb_cat <- dfc %>%
  transmute(
    seccion   = "gB (único vs mixto)",
    categoria = ifelse(gb_sum == 1, "gB único", "gB mixto"),
    CV, CVlog10
  )

gh_cat <- dfc %>%
  transmute(
    seccion   = "gH (único vs mixto)",
    categoria = ifelse(gh_sum == 1, "gH único", "gH mixto"),
    CV, CVlog10
  )

gbgh_cat <- dfc %>%
  transmute(
    seccion   = "gB+gH (único vs mixto)",
    categoria = ifelse(gb_sum == 1 & gh_sum == 1, "gB+gH único", "gB+gH mixto"),
    CV, CVlog10
  )

## Niveles
sec_levels <- c(
  "gB (individuales únicos)",
  "gH (individuales únicos)",
  "gB (único vs mixto)",
  "gH (único vs mixto)",
  "gB+gH (único vs mixto)"
)

cat_levels <- c(
  "gB1", "gB2", "gB3", "gB4",
  "gH1", "gH2",
  "gB único", "gB mixto",
  "gH único", "gH mixto",
  "gB+gH único", "gB+gH mixto"
)

df_plot <- bind_rows(gb_ind, gh_ind, gb_cat, gh_cat, gbgh_cat) %>%
  mutate(
    seccion   = factor(seccion, levels = sec_levels),
    categoria = factor(categoria, levels = cat_levels)
  )

## Paleta refinada
cols <- c(
  "gB1"         = "#5B9BD5",
  "gB2"         = "#5B9BD5",
  "gB3"         = "#5B9BD5",
  "gB4"         = "#5B9BD5",
  "gH1"         = "#D98FD6",
  "gH2"         = "#D98FD6",
  "gB único"    = "#4CAF7D",
  "gB mixto"    = "#E07070",
  "gH único"    = "#4CAF7D",
  "gH mixto"    = "#E07070",
  "gB+gH único" = "#4CAF7D",
  "gB+gH mixto" = "#E07070"
)

## Estadísticas
max_log_by_seccion <- function(data, sec_label) {
  data %>% filter(seccion == sec_label) %>% pull(CVlog10) %>% max(na.rm = TRUE)
}

st_gb_ind <- gb_ind %>%
  pairwise_wilcox_test(CV ~ categoria, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(seccion = sec_levels[1]) %>%
  arrange(group1, group2)

max_gb_log <- max_log_by_seccion(df_plot, sec_levels[1])
st_gb_ind$y.position <- max_gb_log + seq(0.12, by = 0.12, length.out = nrow(st_gb_ind))

st_gh_ind <- gh_ind %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[2], y.position = max_log_by_seccion(df_plot, sec_levels[2]) + 0.12)

st_gb <- gb_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[3], y.position = max_log_by_seccion(df_plot, sec_levels[3]) + 0.12)

st_gh <- gh_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[4], y.position = max_log_by_seccion(df_plot, sec_levels[4]) + 0.12)

st_gbgh <- gbgh_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[5], y.position = max_log_by_seccion(df_plot, sec_levels[5]) + 0.12)

stats_all <- bind_rows(st_gb_ind, st_gh_ind, st_gb, st_gh, st_gbgh) %>%
  mutate(
    p.signif = dplyr::coalesce(p.signif, p.adj.signif),
    seccion  = factor(seccion, levels = sec_levels)
  )

stats_sig <- stats_all %>% filter(p.signif != "ns")

## Plot
p <- ggplot(df_plot, aes(categoria, CVlog10, fill = categoria)) +
  geom_boxplot(
    width          = 0.55,
    color          = "grey20",
    linewidth      = 0.42,
    outlier.shape  = 21,
    outlier.size   = 1.7,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.35,
    outlier.alpha  = 0.75
  ) +
  geom_jitter(
    aes(color = categoria),
    width = 0.10, size = 0.9, alpha = 0.30, show.legend = FALSE
  ) +
  facet_grid(
    ~ seccion,
    scales = "free_x",
    space  = "free_x",
    switch = "x",
    drop   = FALSE
  ) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  scale_color_manual(
    values = c(
      "gB1" = "#1A5FA8", "gB2" = "#1A5FA8", "gB3" = "#1A5FA8", "gB4" = "#1A5FA8",
      "gH1" = "#A04EA0", "gH2" = "#A04EA0",
      "gB único" = "#2E7D52", "gB mixto" = "#B04040",
      "gH único" = "#2E7D52", "gH mixto" = "#B04040",
      "gB+gH único" = "#2E7D52", "gB+gH mixto" = "#B04040"
    )
  ) +
  labs(x = NULL, y = "Carga viral CMV Log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.12)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    # Grid: solo horizontales suaves, sin verticales
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    # Cuadrante con borde fino y fondo muy levemente tintado
    panel.border       = element_rect(color = "grey55", fill = NA, linewidth = 0.45),
    panel.background   = element_rect(fill = "#FAFAFA", color = NA),
    # Strip labels debajo
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "grey90", color = "grey70", linewidth = 0.35),
    strip.text         = element_text(size = 10.5, color = "grey20", face = "bold",
                                      margin = margin(t = 4, b = 4)),
    # Ejes
    axis.line          = element_blank(),   # el borde del panel ya hace de eje
    axis.ticks         = element_line(color = "grey50", linewidth = 0.35),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 10, color = "grey20", margin = margin(t = 3)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    # Separación entre paneles
    panel.spacing.x    = grid::unit(10, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    stats_sig,
    label        = "p.signif",
    xmin         = "group1",
    xmax         = "group2",
    y.position   = "y.position",
    tip.length   = 0.012,
    bracket.size = 0.42,
    size         = 3.8,
    color        = "grey25",
    inherit.aes  = FALSE
  )

p

## Resumen numérico
resumen_CV_gt <- df_plot %>%
  group_by(seccion, categoria) %>%
  summarise(
    n   = dplyr::n(),
    med = median(CV, na.rm = TRUE),
    p25 = quantile(CV, .25, na.rm = TRUE),
    p75 = quantile(CV, .75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(seccion, categoria) %>%
  mutate(
    med_fmt = scales::number(med, accuracy = 1, big.mark = ".", decimal.mark = ","),
    p25_fmt = scales::number(p25, accuracy = 1, big.mark = ".", decimal.mark = ","),
    p75_fmt = scales::number(p75, accuracy = 1, big.mark = ".", decimal.mark = ","),
    resumen_texto = paste0(
      "mediana de ", med_fmt,
      " UI/mL (IQR: ", p25_fmt, "–", p75_fmt,
      "; n = ", n, ")"
    )
  )

resumen_CV_gt

## Tabla de p-valores
tabla_pvalores <- stats_all %>%
  select(seccion, group1, group2, p, p.adj, p.signif) %>%
  mutate(
    p_lab       = ifelse(is.na(p), NA_character_,
                         formatC(p, format = "f", digits = 3, decimal.mark = ",")),
    p_adj_lab   = ifelse(is.na(p.adj), NA_character_,
                         formatC(p.adj, format = "f", digits = 3, decimal.mark = ",")),
    comparacion = paste0(group1, " vs ", group2)
  ) %>%
  select(seccion, comparacion, p_lab, p_adj_lab, p.signif)

tabla_pvalores


## ---- [4e] Figura 25. Valores de Ct segun el genotipo ----
# Apartado 5.4.4. Dos paneles, gB y gH.

## ---- [a] *R14. Boxplot valores Ct — dos paneles gB / gH ====

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Formato largo
df_ct_long <- muestras_sel %>%
  select(GB1CT, GB2CT, GB3CT, GB4CT, GH1CT, GH2CT) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "Genotipo",
    values_to = "CT"
  ) %>%
  mutate(
    Genotipo = factor(
      Genotipo,
      levels = c("GB1CT", "GB2CT", "GB3CT", "GB4CT", "GH1CT", "GH2CT")
    ),
    Etiqueta = recode(
      Genotipo,
      "GB1CT" = "gB1", "GB2CT" = "gB2", "GB3CT" = "gB3", "GB4CT" = "gB4",
      "GH1CT" = "gH1", "GH2CT" = "gH2"
    ),
    Glucoproteina = if_else(grepl("^gB", Etiqueta), "", "")
  ) %>%
  filter(!is.na(CT))

# Paletas
cols_fill <- c(
  "GB1CT" = "#5B9BD5",
  "GB2CT" = "#5B9BD5",
  "GB3CT" = "#5B9BD5",
  "GB4CT" = "#5B9BD5",
  "GH1CT" = "#D98FD6",
  "GH2CT" = "#D98FD6"
)

cols_jitter <- c(
  "GB1CT" = "#1A5FA8",
  "GB2CT" = "#1A5FA8",
  "GB3CT" = "#1A5FA8",
  "GB4CT" = "#1A5FA8",
  "GH1CT" = "#A04EA0",
  "GH2CT" = "#A04EA0"
)

y_max_plot <- max(df_ct_long$CT, na.rm = TRUE) + 1.5

# Plot
p_ct <- ggplot(df_ct_long, aes(x = Etiqueta, y = CT, fill = Genotipo)) +
  geom_boxplot(
    width          = 0.50,
    color          = "grey25",
    linewidth      = 0.45,
    outlier.shape  = 21,
    outlier.size   = 1.8,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.4,
    outlier.alpha  = 0.7
  ) +
  geom_jitter(
    aes(color = Genotipo),
    width = 0.12, size = 1.0, alpha = 0.35, show.legend = FALSE
  ) +
  facet_grid(
    ~ Glucoproteina,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values  = cols_fill)  +
  scale_color_manual(values = cols_jitter) +
  scale_y_continuous(
    name   = "Valor Ct",
    limits = c(24, y_max_plot),
    breaks = seq(25, 40, 5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    # Línea de eje solo en el panel izquierdo para Y, y en ambos para X
    axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
    axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 11, color = "grey20", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    # Strip
    strip.background   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold", color = "grey20",
                                      margin = margin(b = 6)),
    panel.spacing.x    = unit(20, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off")

p_ct


## ---- [5a] Tabla 8. Caracteristicas de los episodios de DNAemia ----
# Apartado 5.5.1. Comparacion por perfil serologico D/R.

library(dplyr)
library(tibble)
library(flextable)

## ── Preparación de datos ──────────────────────────────────────────────────

dnaemias_t <- dnaemias %>%
  mutate(
    Fini      = as.Date(Fini),
    Ffin      = as.Date(Ffin),
    FultCVpos = as.Date(FultCVpos),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    no_tratado  = case_when(tto == 0 ~ 1L, tto == 1 ~ 0L, TRUE ~ NA_integer_),
    sintomatico = case_when(clinica == "VA" ~ 0L, !is.na(clinica) ~ 1L, TRUE ~ NA_integer_),
    alt_pp_bin  = if_else(!is.na(alt_pp) & alt_pp == 1, 1L, 0L)
  )

grupos <- c("Global", "D-/R+", "D+/R+", "D+/R-")

## ── Funciones auxiliares ──────────────────────────────────────────────────

fmt_iqr <- function(data, var, grp) {
  d <- if (grp == "Global") data else filter(data, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  if (nrow(d) == 0) return("\u2013")
  med <- median(d[[var]])
  q1  <- quantile(d[[var]], 0.25)
  q3  <- quantile(d[[var]], 0.75)
  sprintf("%.1f (%.1f\u2013%.1f)", med, q1, q3)
}

fmt_npct <- function(data, var, grp) {
  d <- if (grp == "Global") data else filter(data, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  if (nrow(d) == 0) return("\u2013")
  n   <- sum(d[[var]] == 1)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

fmt_p <- function(p) {
  if (is.na(p) || is.null(p)) return("")
  if (p < 0.001) return("< 0,001")
  sub("\\.", ",", sprintf("%.3f", p))
}

## Episodios por paciente (nivel paciente)
base_pac <- dnaemias_t %>% select(NHC, serDR, nepi_gtc) %>% distinct()

fmt_iqr_pac <- function(var, grp) {
  d <- if (grp == "Global") base_pac else filter(base_pac, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  med <- median(d[[var]])
  q1  <- quantile(d[[var]], 0.25)
  q3  <- quantile(d[[var]], 0.75)
  sprintf("%.1f (%.1f\u2013%.1f)", med, q1, q3)
}

## Episodios con >1 muestra genotipada
fmt_mas1 <- function(grp) {
  d   <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n   <- sum(d$PCR_gt > 1, na.rm = TRUE)
  tot <- nrow(d)
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## alt_pp: NA → 0, excluir D-/R+
fmt_altpp <- function(grp) {
  if (grp == "D-/R+") return("\u2013")
  d   <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n   <- sum(d$alt_pp_bin, na.rm = TRUE)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## Resistencia: mostrar – si n = 0
fmt_resist <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n   <- sum(d$CMV_R == 1, na.rm = TRUE)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## No tratado: mostrar – si n = 0
fmt_notrat <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  d <- filter(d, !is.na(no_tratado))
  n   <- sum(d$no_tratado == 1)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## ── p-values ──────────────────────────────────────────────────────────────

calc_p_cat <- function(data, var) {
  tab <- data %>%
    filter(!is.na(serDR), !is.na(.data[[var]]), serDR != "Global") %>%
    group_by(serDR) %>%
    summarise(ev = sum(.data[[var]] == 1), no = sum(.data[[var]] == 0), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(chisq.test(tab)$p.value, error = function(e) fisher.test(tab)$p.value)
}

calc_p_cont <- function(data, var) {
  d <- filter(data, !is.na(.data[[var]]), !is.na(serDR))
  tryCatch(
    kruskal.test(as.formula(paste(var, "~ serDR")), data = d)$p.value,
    error = function(e) NA_real_
  )
}

p_epi_pac  <- kruskal.test(nepi_gtc ~ serDR, data = base_pac)$p.value
p_pcr_gt   <- calc_p_cont(dnaemias_t, "PCR_gt")
p_mas1_pcr <- calc_p_cat(
  dnaemias_t %>% mutate(mas1 = as.integer(PCR_gt > 1)), "mas1"
)
p_primer   <- calc_p_cat(dnaemias_t, "primer_mil")

p_altpp <- {
  d <- filter(dnaemias_t, serDR %in% c("D+/R+", "D+/R-"))
  tab <- d %>%
    group_by(serDR) %>%
    summarise(ev = sum(alt_pp_bin), no = n() - sum(alt_pp_bin), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

p_sint   <- calc_p_cat(dnaemias_t, "sintomatico")
p_ingr   <- calc_p_cat(dnaemias_t, "ingr")
p_ingcmv <- calc_p_cat(dnaemias_t, "ing_cmv")

p_notrat <- {
  d <- filter(dnaemias_t, !is.na(no_tratado), serDR %in% c("D-/R+", "D+/R+"))
  tab <- d %>%
    group_by(serDR) %>%
    summarise(ev = sum(no_tratado), no = n() - sum(no_tratado), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

p_cvini  <- calc_p_cont(dnaemias_t, "cvini")
p_cvpico <- calc_p_cont(dnaemias_t, "cvpico")
p_aclar  <- calc_p_cont(filter(dnaemias_t, !is.na(aclaramiento_dias)), "aclaramiento_dias")
p_auc    <- calc_p_cont(filter(dnaemias_t, !is.na(AUCCV)), "AUCCV")

## ── Construcción de filas ─────────────────────────────────────────────────

make_row <- function(variable, vals, p_val, es_header = FALSE) {
  tibble(
    Variable = variable,
    Global   = vals[1],
    `D-/R+`  = vals[2],
    `D+/R+`  = vals[3],
    `D+/R-`  = vals[4],
    p        = if (es_header) "" else fmt_p(p_val),
    header   = es_header,
    p_raw    = if (es_header) NA_real_ else p_val
  )
}

tabla_r2 <- bind_rows(
  
  ## ── Sección 1: Distribución ──
  make_row("Distribución de episodios y muestras", rep("", 4), NA, es_header = TRUE),
  make_row(
    "Episodios por paciente, mediana (RIC)",
    sapply(grupos, fmt_iqr_pac, var = "nepi_gtc"),
    p_epi_pac
  ),
  make_row(
    "Muestras genotipadas por episodio, mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(dnaemias_t, "PCR_gt", g)),
    p_pcr_gt
  ),
  make_row(
    "Episodios con >1 muestra genotipada, n (%)",
    sapply(grupos, fmt_mas1),
    p_mas1_pcr
  ),
  
  ## ── Sección 2: Características clínicas ──
  make_row("Características clínicas", rep("", 4), NA, es_header = TRUE),
  make_row(
    "Primer episodio de CV \u2265 1.000 UI/mL, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "primer_mil", g)),
    p_primer
  ),
  make_row(
    "Suspensión precoz de la profilaxis primaria, n (%)",
    sapply(grupos, fmt_altpp),
    p_altpp
  ),
  make_row(
    "Episodios sintomáticos, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "sintomatico", g)),
    p_sint
  ),
  make_row(
    "Hospitalización durante el episodio, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "ingr", g)),
    p_ingr
  ),
  make_row(
    "Hospitalización atribuible a CMV, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "ing_cmv", g)),
    p_ingcmv
  ),
  make_row(
    "Resistencia antiviral demostrada, n (%)",
    sapply(grupos, fmt_resist),
    NA_real_
  ),
  make_row(
    "Episodio no tratado con antiviral, n (%)",
    sapply(grupos, fmt_notrat),
    p_notrat
  ),
  
  ## ── Sección 3: Cinética viral ──
  make_row("Características de cinética viral", rep("", 4), NA, es_header = TRUE),
  make_row(
    "CV inicial (UI/mL), mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(dnaemias_t, "cvini", g)),
    p_cvini
  ),
  make_row(
    "CV pico (UI/mL), mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(dnaemias_t, "cvpico", g)),
    p_cvpico
  ),
  make_row(
    "Tiempo a aclaramiento (días), mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(filter(dnaemias_t, !is.na(aclaramiento_dias)),
                                       "aclaramiento_dias", g)),
    p_aclar
  ),
  make_row(
    "AUC de exposición viral, mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(filter(dnaemias_t, !is.na(AUCCV)), "AUCCV", g)),
    p_auc
  )
)

## ── Flextable ─────────────────────────────────────────────────────────────

idx_headers <- which(tabla_r2$header)
idx_sig     <- which(!is.na(tabla_r2$p_raw) & tabla_r2$p_raw < 0.05)

ft_r2 <- tabla_r2 %>%
  select(Variable, Global, `D-/R+`, `D+/R+`, `D+/R-`, p) %>%
  flextable() %>%
  theme_booktabs() %>%
  set_header_labels(
    Variable = "Característica",
    Global   = "Global\nn = 234",
    `D-/R+`  = "D\u2212/R+\nn = 21",
    `D+/R+`  = "D+/R+\nn = 176",
    `D+/R-`  = "D+/R\u2212\nn = 37",
    p        = "p"
  ) %>%
  align(j = "Variable", align = "left",   part = "all") %>%
  align(j = c("Global","D-/R+","D+/R+","D+/R-","p"), align = "center", part = "all") %>%
  bold(i = idx_headers, j = "Variable", part = "body") %>%
  bg(i = idx_headers, bg = "#F2F2F2", part = "body") %>%
  bold(i = idx_sig, j = "p", part = "body") %>%
  padding(
    i = setdiff(seq_len(nrow(tabla_r2)), idx_headers),
    j = "Variable", padding.left = 14, part = "body"
  ) %>%
  bold(part = "header") %>%
  bg(bg = "#F2F2F2", part = "header") %>%
  add_footer_lines(values = c(
    "Los datos se presentan como mediana (RIC) o n (%), según corresponda.",
    "CV: carga viral; AUC: área bajo la curva de exposición viral acumulada.",
    "Comparaciones entre grupos: prueba de Kruskal-Wallis para variables continuas; ji-cuadrado o test exacto de Fisher para variables categóricas.",
    "– : ningún caso registrado en ese grupo.",
    "* p < 0,05 (valor destacado en negrita)."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left") %>%
  fontsize(part = "body",   size = 10) %>%
  fontsize(part = "header", size = 10) %>%
  autofit()

ft_r2


## ---- [5b] Tablas 9, 10 y 11. Frecuencias genotipicas en muestras y episodios ----
# Apartado 5.5.2. gB, gH y la combinacion gB/gH.

library(dplyr)
library(tidyr)
library(flextable)

## ── gB muestras ──────────────────────────────────────────────────────────

df_gb_m <- muestras_sel %>%
  mutate(
    gb_sum  = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(gb_sum >= 1, !is.na(gB_combo))

n_m_total   <- nrow(df_gb_m)                          # 483
n_m_gb_uni  <- sum(df_gb_m$gb_sum == 1)               # únicos
n_m_gb_mix  <- sum(df_gb_m$gb_sum > 1)                # mixtos

freq_gb_uni_m <- df_gb_m %>%
  filter(gb_sum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_uni, 1))

freq_gb_mix_m <- df_gb_m %>%
  filter(gb_sum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_mix, 1)) %>%
  arrange(desc(n))

## ── gH muestras ──────────────────────────────────────────────────────────

df_gh_m <- muestras_sel %>%
  mutate(
    gh_sum   = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    })
  ) %>%
  filter(gh_sum >= 1, !is.na(gH_combo))

n_m_gh_uni <- sum(df_gh_m$gh_sum == 1)
n_m_gh_mix <- sum(df_gh_m$gh_sum > 1)

freq_gh_uni_m <- df_gh_m %>%
  filter(gh_sum == 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gh_uni, 1))

freq_gh_mix_m <- df_gh_m %>%
  filter(gh_sum > 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gh_mix, 1))

## ── gB+gH muestras ───────────────────────────────────────────────────────

df_gbgh_m <- muestras_sel %>%
  mutate(
    gb_sum   = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gh_sum   = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    mixto = gb_sum > 1 | gh_sum > 1
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

n_m_gbgh_uni <- sum(!df_gbgh_m$mixto)
n_m_gbgh_mix <- sum(df_gbgh_m$mixto)

freq_gbgh_uni_m <- df_gbgh_m %>%
  filter(!mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gbgh_uni, 1)) %>%
  arrange(desc(n))

freq_gbgh_mix_m <- df_gbgh_m %>%
  filter(mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gbgh_mix, 1)) %>%
  arrange(desc(n))


## ── gB episodios ─────────────────────────────────────────────────────────

df_gb_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gB_combo))

n_e_total  <- nrow(df_gb_e)                           # 234
n_e_gb_uni <- sum(df_gb_e$gbsum == 1)
n_e_gb_mix <- sum(df_gb_e$gbsum > 1)

freq_gb_uni_e <- df_gb_e %>%
  filter(gbsum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gb_uni, 1))

freq_gb_mix_e <- df_gb_e %>%
  filter(gbsum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gb_mix, 1)) %>%
  arrange(desc(n))

## ── gH episodios ─────────────────────────────────────────────────────────

df_gh_e <- dnaemias %>%
  mutate(
    gH_combo = apply(select(., gh1:gh2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gH_combo))

n_e_gh_uni <- sum(df_gh_e$ghsum == 1)
n_e_gh_mix <- sum(df_gh_e$ghsum > 1)

freq_gh_uni_e <- df_gh_e %>%
  filter(ghsum == 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gh_uni, 1))

freq_gh_mix_e <- df_gh_e %>%
  filter(ghsum > 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gh_mix, 1))

## ── gB+gH episodios ──────────────────────────────────────────────────────

df_gbgh_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., gh1:gh2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    mixto = gbsum > 1 | ghsum > 1
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

n_e_gbgh_uni <- sum(!df_gbgh_e$mixto)
n_e_gbgh_mix <- sum(df_gbgh_e$mixto)

freq_gbgh_uni_e <- df_gbgh_e %>%
  filter(!mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gbgh_uni, 1)) %>%
  arrange(desc(n))

freq_gbgh_mix_e <- df_gbgh_e %>%
  filter(mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gbgh_mix, 1)) %>%
  arrange(desc(n))


fmt_np <- function(n, pct) sprintf("%d (%s)", n, format(pct, nsmall = 1))

hacer_ft <- function(df, titulo, nota) {
  
  idx_headers <- which(df$header)
  idx_indent  <- which(!df$header)
  
  ft <- df %>%
    select(Variable, Muestras, Episodios) %>%
    flextable() %>%
    theme_booktabs() %>%
    set_header_labels(
      Variable  = titulo,
      Muestras  = paste0("Muestras\n(n = ", n_m_total, ")"),
      Episodios = paste0("Episodios\n(n = ", n_e_total, ")")
    ) %>%
    align(j = "Variable",              align = "left",   part = "all") %>%
    align(j = c("Muestras","Episodios"), align = "center", part = "all") %>%
    bold(i = idx_headers, j = "Variable", part = "body") %>%
    bg(i = idx_headers, bg = "#F2F2F2", part = "body") %>%
    padding(i = idx_indent, j = "Variable", padding.left = 14, part = "body") %>%
    bold(part = "header") %>%
    bg(bg = "#F2F2F2", part = "header") %>%
    add_footer_lines(values = nota) %>%
    fontsize(part = "footer", size = 8) %>%
    align(part = "footer", align = "left") %>%
    fontsize(part = "body",   size = 10) %>%
    fontsize(part = "header", size = 10) %>%
    autofit()
  
  ft
}


combos_gb <- sort(union(freq_gb_uni_m$gB_combo, freq_gb_uni_e$gB_combo))
combos_gb_mix <- union(freq_gb_mix_m$gB_combo, freq_gb_mix_e$gB_combo)

get_m <- function(df, key) {
  r <- df[df[[1]] == key, ]
  if (nrow(r) == 0) return("\u2013")
  fmt_np(r$n, r$pct)
}

filas_r3 <- bind_rows(
  tibble(Variable = "Genotipos únicos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_uni, round(100*n_m_gb_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gb_uni, round(100*n_e_gb_uni/n_e_total,1)),
         header = TRUE),
  
  purrr::map_dfr(combos_gb, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gb_uni_m, k),
           Episodios = get_m(freq_gb_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_mix, round(100*n_m_gb_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gb_mix, round(100*n_e_gb_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gb_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(rename(freq_gb_mix_m, gB_combo=gB_combo), k),
           Episodios = get_m(rename(freq_gb_mix_e, gB_combo=gB_combo), k),
           header    = FALSE)
  })
)

ft_r3 <- hacer_ft(
  filas_r3,
  "Genotipo gB",
  c("Los porcentajes de gB1\u2013gB4 se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r3


filas_r4 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gH, n (%)",
         Muestras  = fmt_np(n_m_gh_uni, round(100*n_m_gh_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gh_uni, round(100*n_e_gh_uni/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(sort(union(freq_gh_uni_m$gH_combo, freq_gh_uni_e$gH_combo)), function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gh_uni_m, k),
           Episodios = get_m(freq_gh_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gH, n (%)",
         Muestras  = fmt_np(n_m_gh_mix, round(100*n_m_gh_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gh_mix, round(100*n_e_gh_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(union(freq_gh_mix_m$gH_combo, freq_gh_mix_e$gH_combo), function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gh_mix_m, k),
           Episodios = get_m(freq_gh_mix_e, k),
           header    = FALSE)
  })
)

ft_r4 <- hacer_ft(
  filas_r4,
  "Genotipo gH",
  c("Los porcentajes de gH1 y gH2 se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r4


combos_gbgh_uni  <- union(freq_gbgh_uni_m$combo_total,  freq_gbgh_uni_e$combo_total)
combos_gbgh_mix  <- union(freq_gbgh_mix_m$combo_total,  freq_gbgh_mix_e$combo_total)

## Ordenar por frecuencia en muestras (desc)
combos_gbgh_uni <- freq_gbgh_uni_m %>%
  full_join(freq_gbgh_uni_e, by = "combo_total", suffix = c("_m","_e")) %>%
  arrange(desc(replace_na(n_m, 0))) %>%
  pull(combo_total)

combos_gbgh_mix <- freq_gbgh_mix_m %>%
  full_join(freq_gbgh_mix_e, by = "combo_total", suffix = c("_m","_e")) %>%
  arrange(desc(replace_na(n_m, 0))) %>%
  pull(combo_total)

filas_r5 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gB / gH, n (%)",
         Muestras  = fmt_np(n_m_gbgh_uni, round(100*n_m_gbgh_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gbgh_uni, round(100*n_e_gbgh_uni/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gbgh_uni, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gbgh_uni_m, k),
           Episodios = get_m(freq_gbgh_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB / gH, n (%)",
         Muestras  = fmt_np(n_m_gbgh_mix, round(100*n_m_gbgh_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gbgh_mix, round(100*n_e_gbgh_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gbgh_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gbgh_mix_m, k),
           Episodios = get_m(freq_gbgh_mix_e, k),
           header    = FALSE)
  })
)

ft_r5 <- hacer_ft(
  filas_r5,
  "Combinación gB / gH",
  c("Los porcentajes de genotipos únicos se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r3
ft_r4
ft_r5


## ---- [5c] Figura 26. Genotipos mixtos por episodio segun perfil serologico ----
# Apartado 5.5.2.

library(ggplot2)
library(dplyr)
library(tidyr)

# ── 1. Etiquetas del eje x con n= ─────────────────────────────────────────────
n_total <- nrow(df)
n_DR    <- sum(df$serDR == "D+/R-",  na.rm = TRUE)
n_DpR   <- sum(df$serDR == "D+/R+",  na.rm = TRUE)
n_DmR   <- sum(df$serDR == "D-/R+",  na.rm = TRUE)

etiq_x <- c(
  "Global" = paste0("Global\n(n = ", n_total, ")"),
  "D+/R-"  = paste0("D+/R-\n(n = ", n_DR,    ")"),
  "D+/R+"  = paste0("D+/R+\n(n = ", n_DpR,   ")"),
  "D-/R+"  = paste0("D-/R+\n(n = ", n_DmR,   ")")
)

# ── 2a. Tabla por seroestatus ─────────────────────────────────────────────────
tab_ser <- df %>%
  group_by(serDR) %>%
  summarise(
    episodios     = n(),
    `gB mixto`    = 100 * sum(gb_mixto,   na.rm = TRUE) / episodios,
    `gH mixto`    = if (all(is.na(gh_mixto)))    NA_real_ else
      100 * sum(gh_mixto,   na.rm = TRUE) / episodios,
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ else
      100 * sum(gbgh_mixto, na.rm = TRUE) / episodios,
    .groups = "drop"
  )

# ── 2b. Fila global ───────────────────────────────────────────────────────────
tab_global <- df %>%
  summarise(
    serDR         = "Global",
    episodios     = n(),
    `gB mixto`    = 100 * sum(gb_mixto,   na.rm = TRUE) / episodios,
    `gH mixto`    = if (all(is.na(gh_mixto)))    NA_real_ else
      100 * sum(gh_mixto,   na.rm = TRUE) / episodios,
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ else
      100 * sum(gbgh_mixto, na.rm = TRUE) / episodios
  )

# ── 2c. Combinar y pivotar ────────────────────────────────────────────────────
tab_fig <- bind_rows(tab_global, tab_ser) %>%
  pivot_longer(
    cols      = c(`gB mixto`, `gH mixto`, `gB/gH mixto`),
    names_to  = "tipo_genotipo",
    values_to = "porcentaje"
  ) %>%
  filter(!is.na(porcentaje)) %>%
  mutate(
    serDR = factor(serDR, levels = c("Global", "D+/R-", "D+/R+", "D-/R+")),
    tipo_genotipo = factor(tipo_genotipo,
                           levels = c("gB mixto", "gH mixto", "gB/gH mixto"))
  )

# ── 3. Figura ─────────────────────────────────────────────────────────────────
figura_R15 <- ggplot(tab_fig,
                     aes(x = serDR, y = porcentaje, fill = tipo_genotipo)) +
  geom_col(
    position  = position_dodge(width = 0.8),
    width     = 0.7,
    color     = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 1.5,           # separador visual Global | por seroestatus
    color      = "grey75",
    linewidth  = 0.4,
    linetype   = "dashed"
  ) +
  geom_text(
    aes(label = paste0(formatC(porcentaje, digits = 0, format = "f",
                               decimal.mark = ","), "%")),
    position  = position_dodge(width = 0.75),
    vjust     = -0.5,
    size      = 4,
    color     = "grey20",
    na.rm     = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "gB mixto"    = "#7EC0EE",
      "gH mixto"    = "#EEA9B8",
      "gB/gH mixto" = "tomato3"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = etiq_x) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.08)),
    breaks = seq(0, 50, 10),
    labels = scales::label_number(suffix = "%", decimal.mark = ",")
  ) +
  labs(x = "", y = "Frecuencia (%)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_line(color = "grey60", linewidth = 0.3),
    axis.text          = element_text(color = "grey30", size = 11),
    axis.title         = element_text(color = "grey20", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.45, "cm"),
    legend.spacing.x   = unit(0.4, "cm"),
    plot.margin        = margin(t = 8, r = 10, b = 4, l = 4)
  )

figura_R15


## ---- [5d1] Tabla 12. Cinetica viral segun genotipo individual ----
# Apartado 5.5.3. Episodios con genotipo unico.

library(dplyr)
library(purrr)
library(flextable)

# ── 1. Función genérica ───────────────────────────────────────────────────────
resumen_bloque <- function(data, geno_var, geno_levels, geno_labels,
                           vars, var_labels, tipo_comparacion) {
  map_dfr(seq_along(vars), function(i) {
    v   <- vars[i]
    lab <- var_labels[i]
    
    celdas <- map_chr(geno_levels, function(g) {
      x <- data %>% filter(.data[[geno_var]] == g) %>%
        pull(.data[[v]]) %>% na.omit()
      if (length(x) == 0) return("–")
      sprintf("%s (%s–%s)\nn = %d",
              formatC(median(x),          format = "f", digits = 0),
              formatC(quantile(x, 0.25),  format = "f", digits = 0),
              formatC(quantile(x, 0.75),  format = "f", digits = 0),
              length(x))
    })
    
    grupos <- map(geno_levels, function(g)
      data %>% filter(.data[[geno_var]] == g) %>%
        pull(.data[[v]]) %>% na.omit())
    
    p_val <- tryCatch({
      if (tipo_comparacion == "kruskal") kruskal.test(grupos)$p.value
      else wilcox.test(grupos[[1]], grupos[[2]], exact = FALSE)$p.value
    }, error = function(e) NA_real_)
    
    p_fmt <- if (is.na(p_val)) "" else
      if (p_val < 0.001) "<0,001" else
        sub("\\.", ",", sprintf("%.3f", p_val))
    
    row        <- as_tibble_row(setNames(celdas, geno_labels))
    row$Parámetro <- lab
    row$p         <- p_fmt
    row$p_raw     <- p_val
    row
  })
}

# ── 2. Preparar datos ─────────────────────────────────────────────────────────
dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )

vars     <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
var_labs <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
              "Aclaramiento (días)", "AUC viral")

# ── 3. Bloque gB ──────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  filter(gbsum == 1) %>%
  mutate(gb_tipo = case_when(
    gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
    gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
  ))

bloque_gb <- resumen_bloque(
  data             = df_gb,
  geno_var         = "gb_tipo",
  geno_levels      = c("gB1","gB2","gB3","gB4"),
  geno_labels      = c("gB1","gB2","gB3","gB4"),
  vars             = vars,
  var_labels       = var_labs,
  tipo_comparacion = "kruskal"
)

# ── 4. Bloque gH ──────────────────────────────────────────────────────────────
df_gh <- dnaemias %>%
  filter(ghsum == 1) %>%
  mutate(gh_tipo = case_when(
    gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"
  ))

bloque_gh_raw <- resumen_bloque(
  data             = df_gh,
  geno_var         = "gh_tipo",
  geno_levels      = c("gH1","gH2"),
  geno_labels      = c("gH1","gH2"),
  vars             = vars,
  var_labels       = var_labs,
  tipo_comparacion = "wilcox"
)

# Añadir columnas gB3/gB4 vacías y renombrar gH1/gH2 → gB1/gB2
# (el flextable usará los nombres de columna, las cabeceras se personalizan)
bloque_gh <- bloque_gh_raw %>%
  rename(gB1 = gH1, gB2 = gH2) %>%
  mutate(gB3 = NA_character_, gB4 = NA_character_) %>%
  select(Parámetro, gB1, gB2, gB3, gB4, p, p_raw)

# ── 5. Filas cabecera ─────────────────────────────────────────────────────────
cab <- function(txt) tibble(
  Parámetro = txt,
  gB1 = NA_character_, gB2 = NA_character_,
  gB3 = NA_character_, gB4 = NA_character_,
  p   = NA_character_,  p_raw = NA_real_
)

# ── 6. Tabla completa ─────────────────────────────────────────────────────────
tabla_final <- bind_rows(
  cab("Glucoproteína gB"),
  bloque_gb,
  cab("Glucoproteína gH"),
  bloque_gh
)

# ── 7. Índices para formato ───────────────────────────────────────────────────
idx_cab <- which(tabla_final$Parámetro %in%
                   c("Glucoproteína gB", "Glucoproteína gH"))
idx_sig <- which(!is.na(tabla_final$p_raw) & tabla_final$p_raw < 0.05)
idx_gh  <- which(tabla_final$Parámetro == "Glucoproteína gH") +
  seq_len(length(vars))   # filas de datos de gH

# ── 8. Flextable ──────────────────────────────────────────────────────────────
ft <- tabla_final %>%
  select(-p_raw) %>%
  flextable() %>%
  # cabecera con nombres correctos por sección
  set_header_labels(
    Parámetro = "Parámetro",
    gB1 = "gB1 / gH1",
    gB2 = "gB2 / gH2",
    gB3 = "gB3",
    gB4 = "gB4",
    p   = "p"
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  # cabeceras de sección
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab,   bg = "#F2F2F2", part = "body") %>%
  # p significativas en negrita
  bold(i = idx_sig, j = "p", part = "body") %>%
  # alineación
  align(j = "Parámetro", align = "left",   part = "all") %>%
  align(j = c("gB1","gB2","gB3","gB4","p"),
        align = "center", part = "all") %>%
  # sangría en filas de datos
  padding(i = setdiff(seq_len(nrow(tabla_final)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  # gris para celdas vacías (gB3/gB4 en filas de gH)
  bg(i = idx_gh, j = c("gB3","gB4"), bg = "#F8F8F8", part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Datos presentados como mediana (RIC). CV: carga viral; AUC: área bajo la curva viral.",
    "Solo se incluyen episodios con genotipo único de cada glucoproteína.",
    "Glucoproteína gB: comparación entre gB1, gB2, gB3 y gB4 (Kruskal-Wallis).",
    "Glucoproteína gH: comparación entre gH1 y gH2 (U de Mann-Whitney).",
    "Las columnas gB3 y gB4 no aplican a gH (–). Los valores de p < 0,05 en negrita."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft


## ---- [5d2] Figura 27. Cinetica viral entre los genotipos unicos de gB ----
# Apartado 5.5.3. Incluye post hoc de Dunn con correccion de Bonferroni.

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyr)
library(patchwork)

# ── 1. Datos ──────────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    cvini_log         = log10(cvini),
    cvpico_log        = log10(cvpico)
  ) %>%
  filter(gbsum == 1) %>%
  mutate(
    gb_tipo = factor(case_when(
      gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
    ), levels = c("gB1","gB2","gB3","gB4"))
  )

# ── 2. Comparaciones post-hoc significativas ──────────────────────────────────
comp_aclar <- tibble(
  group1       = c("gB1", "gB2"),
  group2       = c("gB3", "gB3"),
  p.adj.signif = c("**", "**"),
  y.position   = c(160, 175)
)

comp_auc <- tibble(
  group1       = "gB1",
  group2       = "gB3",
  p.adj.signif = "*",
  y.position   = 390
)

comp_vacio <- tibble(
  group1       = character(0),
  group2       = character(0),
  p.adj.signif = character(0),
  y.position   = numeric(0)
)

# ── 3. Colores: paleta de azules equilibrada (gB3 no tan oscuro) ───────────────
cols_gb <- c(
  "gB1" = "#D6EAF8",
  "gB2" = "#85C1E9",
  "gB3" = "#3498DB",
  "gB4" = "#A9CCE3"
)

cols_jitter <- c(
  "gB1" = "#2E86C1",
  "gB2" = "#2471A3",
  "gB3" = "#1B4F72",
  "gB4" = "#5499C7"
)

# ── 4. Función para construir cada panel ──────────────────────────────────────
make_panel <- function(data, var, ylab, comparaciones,
                       cap_vis = NULL, ylim_max = NULL) {
  
  d <- data %>%
    select(gb_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  if (!is.null(cap_vis)) d <- d %>% filter(valor <= cap_vis)
  
  p <- ggplot(d, aes(x = gb_tipo, y = valor, fill = gb_tipo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    # Mediana resaltada en blanco por encima del relleno
    stat_summary(
      fun = median, geom = "crossbar",
      width = 0.55, fatten = 0,
      color = "white", linewidth = 0.7
    ) +
    geom_jitter(
      aes(color = gb_tipo),
      width = 0.12, size = 0.9, alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_gb) +
    scale_color_manual(values = cols_jitter) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 15) +
    theme(
      legend.position    = "none",
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black", linewidth = 0.4),
      axis.text.x        = element_text(size = 14, color = "black"),
      axis.text.y        = element_text(size = 12, color = "black"),
      axis.title.y       = element_text(size = 14, color = "black", face = "bold"),
      plot.margin        = margin(8, 12, 6, 6)
    )
  
  if (!is.null(ylim_max)) {
    p <- p + scale_y_continuous(
      name   = ylab,
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    )
  } else {
    p <- p + scale_y_continuous(
      name   = ylab,
      expand = expansion(mult = c(0.02, 0.15))
    )
  }
  
  if (nrow(comparaciones) > 0) {
    p <- p + stat_pvalue_manual(
      comparaciones,
      label        = "p.adj.signif",
      xmin         = "group1",
      xmax         = "group2",
      y.position   = "y.position",
      tip.length   = 0.01,
      bracket.size = 0.45,
      size         = 6,
      inherit.aes  = FALSE
    )
  }
  
  p
}

# ── 5. Construir los 4 paneles ────────────────────────────────────────────────
p_cvini <- make_panel(
  data          = df_gb,
  var           = "cvini_log",
  ylab          = "log10 CV inicial (UI/mL)",
  comparaciones = comp_vacio
)

p_cvpico <- make_panel(
  data          = df_gb,
  var           = "cvpico_log",
  ylab          = "log10 CV pico (UI/mL)",
  comparaciones = comp_vacio
)

p_aclar <- make_panel(
  data          = df_gb,
  var           = "aclaramiento_dias",
  ylab          = "Aclaramiento (días)",
  comparaciones = comp_aclar,
  cap_vis       = 200,
  ylim_max      = 200
)

p_auc <- make_panel(
  data          = df_gb,
  var           = "AUCCV",
  ylab          = "AUC viral (UI/mL\u00b7d\u00eda)",
  comparaciones = comp_auc,
  cap_vis       = 500,
  ylim_max      = 500
)

# ── 6. Combinar 2x2 ──────────────────────────────────────────────────────────
figura_gb_cinetica <- (p_cvini | p_cvpico) / (p_aclar | p_auc)

figura_gb_cinetica


library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyr)
library(patchwork)

# ── 1. Datos ──────────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    cvini_log         = log10(cvini),
    cvpico_log        = log10(cvpico)
  ) %>%
  filter(gbsum >= 1) %>%
  mutate(
    gb_tipo = factor(case_when(
      gbsum >  1 ~ "gB mixto",
      gb1   == 1 ~ "gB1",
      gb2   == 1 ~ "gB2",
      gb3   == 1 ~ "gB3",
      gb4   == 1 ~ "gB4"
    ), levels = c("gB1","gB2","gB3","gB4","gB mixto"))
  )

# ── 2. Colores ────────────────────────────────────────────────────────────────
cols_gb <- c(
  "gB1"      = "#AED6F1",
  "gB2"      = "#5DADE2",
  "gB3"      = "#1A5276",
  "gB4"      = "#85929E",
  "gB mixto" = "#D7BDE2"
)

cols_jitter <- c(
  "gB1"      = "#2980B9",
  "gB2"      = "#1A6FA8",
  "gB3"      = "#0D2B45",
  "gB4"      = "#555F66",
  "gB mixto" = "#7D3C98"
)

col_mixto <- "#7D3C98"   # color para asteriscos/corchetes con gB mixto
col_otros <- "grey20"    # color para el resto

# ── 3. Tests Kruskal-Wallis (informativos) ────────────────────────────────────
cat("\n--- Kruskal-Wallis ---\n")
print(df_gb %>% kruskal_test(cvini_log         ~ gb_tipo))
print(df_gb %>% kruskal_test(cvpico_log        ~ gb_tipo))
print(df_gb %>% filter(aclaramiento_dias <= 200) %>%
        kruskal_test(aclaramiento_dias   ~ gb_tipo))
print(df_gb %>% filter(AUCCV <= 500) %>%
        kruskal_test(AUCCV               ~ gb_tipo))

# ── 4. Función para calcular comparaciones post-hoc ──────────────────────────
calcular_comp <- function(data, var, y_base, y_step) {
  d <- data %>%
    select(gb_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  res <- d %>%
    dunn_test(valor ~ gb_tipo, p.adjust.method = "BH") %>%
    filter(p.adj < 0.05) %>%
    mutate(
      p.adj.signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      ),
      involucra_mixto = group1 == "gB mixto" | group2 == "gB mixto"
    ) %>%
    arrange(involucra_mixto, p.adj) %>%       # primero las no-mixto, luego mixto
    mutate(y.position = y_base + (row_number() - 1) * y_step) %>%
    select(group1, group2, p.adj.signif, y.position, involucra_mixto)
  
  res
}

# Calcular comparaciones para cada variable
comp_cvini  <- calcular_comp(df_gb, "cvini_log",
                             y_base =   6.8, y_step = 0.35)
comp_cvpico <- calcular_comp(df_gb, "cvpico_log",
                             y_base =   7.0, y_step = 0.35)
comp_aclar  <- calcular_comp(df_gb %>% filter(aclaramiento_dias <= 200),
                             "aclaramiento_dias",
                             y_base = 155, y_step = 13)
comp_auc    <- calcular_comp(df_gb %>% filter(AUCCV <= 500),
                             "AUCCV",
                             y_base = 380, y_step = 32)

cat("\n--- Comparaciones significativas Dunn (BH) ---\n")
print(list(cvini = comp_cvini, cvpico = comp_cvpico,
           aclar = comp_aclar,    auc    = comp_auc))

# ── 5. Función para construir cada panel ──────────────────────────────────────
make_panel <- function(data, var, ylab, comparaciones,
                       cap_vis = NULL, ylim_max = NULL) {
  
  d <- data %>%
    select(gb_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  if (!is.null(cap_vis)) d <- d %>% filter(valor <= cap_vis)
  
  p <- ggplot(d, aes(x = gb_tipo, y = valor, fill = gb_tipo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    geom_jitter(
      aes(color = gb_tipo),
      width = 0.12, size = 0.9, alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_gb) +
    scale_color_manual(values = cols_jitter) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 10, color = "grey25"),
      axis.text.y        = element_text(size = 10, color = "grey30"),
      axis.title.y       = element_text(size = 11, color = "grey20"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(8, 12, 6, 6)
    )
  
  if (!is.null(ylim_max)) {
    p <- p + scale_y_continuous(
      name   = ylab,
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    )
  } else {
    p <- p + scale_y_continuous(
      name   = ylab,
      expand = expansion(mult = c(0.02, 0.20))
    )
  }
  
  # Separar comparaciones por color y añadir cada subconjunto en una capa
  if (nrow(comparaciones) > 0) {
    
    comp_otros <- comparaciones %>% filter(!involucra_mixto)
    comp_mix   <- comparaciones %>% filter( involucra_mixto)
    
    if (nrow(comp_otros) > 0) {
      p <- p + stat_pvalue_manual(
        comp_otros,
        label        = "p.adj.signif",
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        tip.length   = 0.01,
        bracket.size = 0.45,
        size         = 4.5,
        color        = col_otros,
        inherit.aes  = FALSE
      )
    }
    
    if (nrow(comp_mix) > 0) {
      p <- p + stat_pvalue_manual(
        comp_mix,
        label        = "p.adj.signif",
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        tip.length   = 0.01,
        bracket.size = 0.45,
        size         = 4.5,
        color        = col_mixto,
        inherit.aes  = FALSE
      )
    }
  }
  
  p
}

# ── 6. Construir los 4 paneles ────────────────────────────────────────────────
p_cvini <- make_panel(
  data          = df_gb,
  var           = "cvini_log",
  ylab          = "log10 CV inicial",
  comparaciones = comp_cvini
)

p_cvpico <- make_panel(
  data          = df_gb,
  var           = "cvpico_log",
  ylab          = "log10 CV pico",
  comparaciones = comp_cvpico
)

p_aclar <- make_panel(
  data          = df_gb,
  var           = "aclaramiento_dias",
  ylab          = "Aclaramiento (días)",
  comparaciones = comp_aclar,
  cap_vis       = 200,
  ylim_max      = 215
)

p_auc <- make_panel(
  data          = df_gb,
  var           = "AUCCV",
  ylab          = "AUC viral (UI/mL\u00b7d\u00eda)",
  comparaciones = comp_auc,
  cap_vis       = 500,
  ylim_max      = 520
)

# ── 7. Combinar 2x2 ───────────────────────────────────────────────────────────
figura_gb_cinetica <- (p_cvini | p_cvpico) / (p_aclar | p_auc)

figura_gb_cinetica
library(dplyr)
library(dunn.test)

# ── Preparar datos ────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  ) %>%
  filter(gbsum == 1) %>%
  mutate(gb_tipo = case_when(
    gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
    gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
  ))

# ── Post-hoc Dunn — Aclaramiento ──────────────────────────────────────────────
dunn.test(
  x      = df_gb$aclaramiento_dias,
  g      = df_gb$gb_tipo,
  method = "bonferroni",
  kw     = TRUE,
  label  = TRUE
)

# ── Post-hoc Dunn — AUC ───────────────────────────────────────────────────────
dunn.test(
  x      = df_gb$AUCCV,
  g      = df_gb$gb_tipo,
  method = "bonferroni",
  kw     = TRUE,
  label  = TRUE
)


## ---- [5e] Tabla 13. Modelo multivariante del aclaramiento en genotipos unicos de gB ----
# Apartado 5.5.3. Regresion lineal multiple y diagnostico de colinealidad.

# ── Modelo clínico gB3 (variables elegidas a priori) ──────────────────────────
df_clin_gb3 <- df_gb3_uni %>%
  select(NHC, aclaramiento_dias,
         gb3_f, risk_f, cvini, cvpico, profilaxis_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N modelo clínico: %d episodios / %d pacientes ──\n",
            nrow(df_clin_gb3), n_distinct(df_clin_gb3$NHC)))

modelo_clinico_gb3 <- lm(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f,
  data = df_clin_gb3
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO CLÍNICO gB3 (variables a priori)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_clinico_gb3))

cat("\n── IC 95% ──\n")
print(confint(modelo_clinico_gb3))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_clinico_gb3))

# Efecto ajustado de gB3 en una línea
coef_gb3 <- coef(modelo_clinico_gb3)["gb3_fgB3"]
ci_gb3   <- confint(modelo_clinico_gb3)["gb3_fgB3", ]
p_gb3    <- summary(modelo_clinico_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

cat(sprintf("\n── Efecto ajustado gB3: %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f ──\n",
            coef_gb3, ci_gb3[1], ci_gb3[2], p_gb3))


library(lmerTest)

modelo_mixto_clinico <- lmer(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + (1 | NHC),
  data = df_clin_gb3
)

print(summary(modelo_mixto_clinico))
print(confint(modelo_mixto_clinico, method = "Wald")["gb3_fgB3", ])

vc  <- as.data.frame(VarCorr(modelo_mixto_clinico))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\nICC: %.3f\n", icc))


## ---- [5f1] Tabla 14. Cinetica viral segun presencia de genotipos mixtos ----
# Apartado 5.5.4.

library(dplyr)
library(flextable)
library(tibble)

# ── 1. Calcular aclaramiento si no existe ─────────────────────────────────────
dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )

# ── 2. Funciones auxiliares ───────────────────────────────────────────────────
fmt_med <- function(x) {
  n   <- sum(!is.na(x))
  med <- median(x, na.rm = TRUE)
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  sprintf("%s (%s\u2013%s)\n(n=%d)",
          formatC(med, format="f", digits=0),
          formatC(q1,  format="f", digits=0),
          formatC(q3,  format="f", digits=0), n)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0,001*")
  txt <- sprintf("%.3f", p)
  txt <- sub("\\.", ",", txt)
  if (p < 0.05) paste0(txt, "*") else txt
}

wilcox_p <- function(data, var, grupo) {
  d <- data %>% filter(!is.na(.data[[var]]), !is.na(.data[[grupo]]))
  if (nrow(d) < 4) return(NA_real_)
  tryCatch(
    wilcox.test(as.formula(paste(var, "~", grupo)), data = d)$p.value,
    error = function(e) NA_real_
  )
}

# ── 3. Función que construye 4 filas para un grupo genotípico ─────────────────
bloque <- function(data, var_mix) {
  d0 <- data %>% filter(.data[[var_mix]] == 0)
  d1 <- data %>% filter(.data[[var_mix]] == 1)
  
  vars  <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
  labs  <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
             "Aclaramiento (días)", "AUC viral")
  
  tibble(
    Parámetro = labs,
    `Genotipo único` = sapply(vars, function(v) fmt_med(d0[[v]])),
    `Genotipo mixto`  = sapply(vars, function(v) fmt_med(d1[[v]])),
    p = sapply(vars, function(v) fmt_p(wilcox_p(data, v, var_mix)))
  )
}

# ── 4. Construir tabla completa ───────────────────────────────────────────────
cabecera <- function(txt) tibble(
  Parámetro        = txt,
  `Genotipo único` = NA_character_,
  `Genotipo mixto` = NA_character_,
  p                = NA_character_
)

tabla_cin <- bind_rows(
  cabecera("Glucoproteína gB"),
  bloque(dnaemias, "gbmix"),
  cabecera("Glucoproteína gH"),
  bloque(dnaemias, "ghmix"),
  cabecera("Glucoproteínas gB y gH"),
  bloque(dnaemias, "gbghmix")
)

# ── 5. Flextable ──────────────────────────────────────────────────────────────
# Índices de filas cabecera
idx_cab <- which(tabla_cin$Parámetro %in%
                   c("Glucoproteína gB", "Glucoproteína gH",
                     "Glucoproteínas gB y gH"))

ft_cin <- flextable(tabla_cin) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  align(j = "Parámetro", align = "left",   part = "all") %>%
  align(j = c("Genotipo único", "Genotipo mixto", "p"),
        align = "center", part = "all") %>%
  padding(i = setdiff(seq_len(nrow(tabla_cin)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Los datos se presentan como mediana (RIC).",
    "CV: carga viral; AUC: área bajo la curva.",
    "Comparaciones mediante la prueba U de Mann-Whitney.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft_cin


## ---- [5f2] Figura 28. Cinetica viral de los genotipos mixtos por perfil serologico ----
# Apartado 5.5.4.

library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggtext)

# ── Variables de control ──────────────────────────────────────────────────────
id_patient <- "NHC"
cv_y_max   <- 7

df_sub <- dnaemias

vars_to_pivot <- c("diasdur_vir", "AUCCV", "cvpico", "cvini")

# ── Formato largo ─────────────────────────────────────────────────────────────
df_long <- df_sub %>%
  mutate(
    gB   = ifelse(gbmix   == 1, "gB mixto",    "gB único"),
    gH   = ifelse(ghmix   == 1, "gH mixto",    "gH único"),
    gBGH = ifelse(gbghmix == 1, "gB+gH mixto", "gB+gH único")
  ) %>%
  pivot_longer(cols = c(gB, gH, gBGH),
               names_to  = "tipo_gen",
               values_to = "genotipo") %>%
  select(serDR, tipo_gen, genotipo,
         !!sym(id_patient), all_of(vars_to_pivot)) %>%
  pivot_longer(cols      = all_of(vars_to_pivot),
               names_to  = "variable",
               values_to = "valor_raw") %>%
  mutate(
    valor = case_when(
      variable %in% c("cvpico", "cvini") & !is.na(valor_raw) ~
        log10(as.numeric(valor_raw) + 1),
      TRUE ~ as.numeric(valor_raw)
    ),
    variable = factor(recode(
      variable,
      "diasdur_vir" = "Aclaramiento (días)",
      "AUCCV"       = "AUC (UI/mL·día)",
      "cvpico"      = "log10 CV pico (UI/mL)",
      "cvini"       = "log10 CV inicial (UI/mL)"
    ), levels = c(
      "log10 CV inicial (UI/mL)",
      "log10 CV pico (UI/mL)",
      "Aclaramiento (días)",
      "AUC (UI/mL·día)"
    )),
    serDR    = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    genotipo = factor(genotipo, levels = c(
      "gB único",    "gB mixto",
      "gH único",    "gH mixto",
      "gB+gH único", "gB+gH mixto"
    ))
  )

# ── Posiciones X ──────────────────────────────────────────────────────────────
x_map <- c(
  "gB único"    = 1,   "gB mixto"    = 2,
  "gH único"    = 3.3, "gH mixto"    = 4.3,
  "gB+gH único" = 5.6, "gB+gH mixto" = 6.6
)
df_plot <- df_long %>%
  mutate(x_pos = x_map[as.character(genotipo)])

# ── Límites visuales por panel ────────────────────────────────────────────────
caps_dyn <- tibble(
  variable = c("log10 CV inicial (UI/mL)", "log10 CV pico (UI/mL)",
               "Aclaramiento (días)", "AUC (UI/mL·día)"),
  cap_vis  = c(cv_y_max, cv_y_max, 100, 700)
)

df_plot_vis <- df_plot %>%
  left_join(caps_dyn, by = "variable") %>%
  filter(is.na(valor) | valor <= cap_vis)

# ── Comparaciones Wilcoxon sobre datos SIN filtrar (solo significativas) ──────
comparaciones <- df_plot %>%
  group_by(variable, serDR, tipo_gen) %>%
  wilcox_test(valor ~ genotipo) %>%
  add_significance("p") %>%
  ungroup() %>%
  filter(p.signif != "ns") %>%
  mutate(
    group1 = case_when(
      tipo_gen == "gB"   ~ "gB único",
      tipo_gen == "gH"   ~ "gH único",
      tipo_gen == "gBGH" ~ "gB+gH único"
    ),
    group2 = case_when(
      tipo_gen == "gB"   ~ "gB mixto",
      tipo_gen == "gH"   ~ "gH mixto",
      tipo_gen == "gBGH" ~ "gB+gH mixto"
    ),
    xmin = x_map[group1],
    xmax = x_map[group2]
  )

base_positions <- df_plot_vis %>%
  group_by(variable, serDR) %>%
  summarise(maxval  = max(valor,   na.rm = TRUE),
            cap_vis = first(cap_vis), .groups = "drop") %>%
  mutate(base = pmin(maxval * 1.12, cap_vis * 0.90))

comp_pos <- comparaciones %>%
  left_join(base_positions, by = c("variable", "serDR")) %>%
  group_by(variable, serDR) %>%
  arrange(p) %>%
  mutate(idx = row_number(),
         y   = base + (idx - 1) * 0.05 * cap_vis) %>%
  ungroup()

# ── Etiquetas serológicas con n ───────────────────────────────────────────────
counts_serDR <- df_sub %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    pacientes = n_distinct(.data[[id_patient]]),
    .groups   = "drop"
  ) %>%
  mutate(label = paste0(
    "<b>", serDR, "</b><br>",
    "<span style='font-size:10pt'>",
    pacientes, " pac. / ", episodios, " ep.</span>"
  ))

serdr_labels        <- counts_serDR$label
names(serdr_labels) <- counts_serDR$serDR

# ── Colores ───────────────────────────────────────────────────────────────────
colores <- c(
  "gB único"    = "#63B8FF",
  "gB mixto"    = "dodgerblue3",
  "gH único"    = "#DDA0DD",
  "gH mixto"    = "#CD69C9",
  "gB+gH único" = "#3CB371",
  "gB+gH mixto" = "#F08080"
)

# ── Figura ────────────────────────────────────────────────────────────────────
p_global <- ggplot(df_plot_vis,
                   aes(x = x_pos, y = valor, fill = genotipo)) +
  geom_boxplot(
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4,
    outlier.shape = 21,
    outlier.size  = 1.4,
    outlier.fill  = "white",
    outlier.color = "grey50",
    outlier.alpha = 0.6
  ) +
  geom_vline(xintercept = c(2.65, 4.95),
             color = "grey80", linewidth = 0.35, linetype = "dashed") +
  facet_grid(variable ~ serDR, scales = "free_y",
             labeller = labeller(serDR = serdr_labels)) +
  geom_text(data = comp_pos,
            aes(x = (xmin + xmax) / 2, y = y, label = p.signif),
            inherit.aes = FALSE, size = 3.8, vjust = 0, color = "grey20") +
  geom_segment(data = comp_pos,
               aes(x    = xmin, xend = xmax,
                   y    = y - 0.012 * cap_vis,
                   yend = y - 0.012 * cap_vis),
               inherit.aes = FALSE,
               linewidth = 0.45, color = "grey30") +
  scale_fill_manual(values = colores) +
  scale_x_continuous(
    breaks = unname(x_map),
    labels = names(x_map),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "none",
    panel.border       = element_rect(color = "grey70", fill = NA,
                                      linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_blank(),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 10, color = "grey25",
                                      angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 12, color = "grey30"),
    strip.background   = element_rect(fill = "grey95", color = "grey70",
                                      linewidth = 0.4),
    strip.text.x       = element_markdown(size = 13, face = "bold",
                                          color = "grey20"),
    strip.text.y       = element_text(size = 12, face = "bold",
                                      color = "grey20"),
    panel.spacing.x    = unit(6,  "pt"),
    panel.spacing.y    = unit(10, "pt"),
    plot.margin        = margin(8, 10, 6, 6)
  )

ggsave("figura_cinetica_global.pdf",
       plot   = p_global,
       width  = 17,
       height = 22,
       units  = "cm",
       device = "pdf")

p_global


## ---- [5f3] Analisis univariante de factores de confusion ----
# Apartado 5.5.4. Seleccion de covariables para el modelo multivariante.

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# ── 1. Preparar datos ─────────────────────────────────────────────────────────
df_uni <- dnaemias %>%
  left_join(pacientes %>% select(NHC, edad_tx, tipotx, ebase), by = "NHC") %>%
  mutate(
    risk_f = factor(
      ifelse(serDR == "D+/R-", "Alto riesgo", "Riesgo intermedio"),
      levels = c("Riesgo intermedio", "Alto riesgo")
    ),
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    Tdesdetx_meses    = as.numeric(Fini - as.Date(Ftx)) / 30.44,
    sintomatico       = case_when(
      clinica == "VA"  ~ 0L,
      !is.na(clinica)  ~ 1L,
      TRUE             ~ NA_integer_
    ),
    gbghmix_f    = factor(gbghmix,  levels = c(0,1), labels = c("Único","Mixto")),
    ingr_f       = factor(ingr,     levels = c(0,1), labels = c("No","Sí")),
    ing_cmv_f    = factor(ing_cmv,  levels = c(0,1), labels = c("No","Sí")),
    sintom_f     = factor(sintomatico, levels = c(0,1), labels = c("No","Sí")),
    profilaxis_f = factor(antiviral, levels = c(0,1), labels = c("No","Sí")),
    tto_f        = factor(tto,      levels = c(0,1), labels = c("No","Sí")),
    cmvr_f       = factor(CMV_R,    levels = c(0,1), labels = c("No","Sí")),
    tipotx_f     = factor(tipotx,   levels = c("UNI","BI")),
    ebase_f      = factor(ebase)
  ) %>%
  filter(!is.na(aclaramiento_dias))

# ── 2. Correlaciones de Spearman (variables continuas) ────────────────────────
cor_vars <- list(
  "Edad receptor (años)"        = "edad_tx",
  "CV inicial (UI/mL)"          = "cvini",
  "CV pico (UI/mL)"             = "cvpico",
  "Tiempo desde TXP (meses)"     = "Tdesdetx_meses"
)

res_cor <- map_dfr(names(cor_vars), function(lab) {
  v <- cor_vars[[lab]]
  x <- as.numeric(df_uni[[v]])
  if (v %in% c("cvini","cvpico")) x <- log10(x + 1)
  idx  <- !is.na(x) & !is.na(df_uni$aclaramiento_dias)
  x_ok <- x[idx]
  y_ok <- df_uni$aclaramiento_dias[idx]
  if (length(x_ok) < 5) return(tibble(Variable = lab, tipo = "continua",
                                      estimador = NA_real_, p_val = NA_real_,
                                      label_est = "n insuf."))
  test <- tryCatch(
    cor.test(x_ok, y_ok, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) return(tibble(Variable = lab, tipo = "continua",
                                   estimador = NA_real_, p_val = NA_real_,
                                   label_est = "error"))
  tibble(
    Variable  = lab,
    tipo      = "continua",
    estimador = round(test$estimate, 3),
    p_val     = test$p.value,
    label_est = sprintf("rho = %+.3f", test$estimate)
  )
})

# ── 3. Comparaciones de grupos (variables categóricas 2 grupos) ───────────────
cat_vars <- list(
  "Genotipos mixtos gB/gH"            = "gbghmix_f",
  "Perfil de riesgo serológico D/R"   = "risk_f",      
  "Hospitalización cualquier causa"   = "ingr_f",
  "Hospitalización por CMV"           = "ing_cmv_f",
  "Enfermedad sintomática"            = "sintom_f",
  "Profilaxis al inicio del episodio" = "profilaxis_f",
  "Tratamiento antiviral"             = "tto_f",
  "Tipo de trasplante (BI vs UNI)"    = "tipotx_f"
)

res_cat <- map_dfr(names(cat_vars), function(lab) {
  v  <- cat_vars[[lab]]
  g  <- df_uni[[v]]
  ac <- df_uni$aclaramiento_dias
  d0 <- ac[g == levels(g)[1] & !is.na(g)]
  d1 <- ac[g == levels(g)[2] & !is.na(g)]
  test <- wilcox.test(d1, d0, exact = FALSE)
  dif  <- median(d1, na.rm = TRUE) - median(d0, na.rm = TRUE)
  tibble(
    Variable  = lab,
    tipo      = "categorica",
    estimador = round(dif, 1),
    p_val     = test$p.value,
    label_est = sprintf("\u0394 = %+.1f d\u00edas", dif)
  )
})

# ── 4. Kruskal-Wallis (variables categóricas >2 grupos) ───────────────────────
kw_eb  <- kruskal.test(aclaramiento_dias ~ ebase_f, data = df_uni)

res_eb <- tibble(
  Variable  = "Enfermedad de base",
  tipo      = "categorica",
  estimador = NA_real_,
  p_val     = kw_eb$p.value,
  label_est = "Kruskal-Wallis"
)

res_eb <- tibble(
  Variable  = "Enfermedad de base",
  tipo      = "categorica",
  estimador = NA_real_,
  p_val     = kw_eb$p.value,
  label_est = "Kruskal-Wallis"
)

# ── 5. Combinar y preparar para figura ────────────────────────────────────────
orden_vars <- c(
  "Edad receptor (años)",
  "Tipo de trasplante (BI vs UNI)",
  "Enfermedad de base",
  "Perfil de riesgo serológico D/R",    # <-- corregido
  "Tiempo desde TXP (meses)",
  "CV inicial (UI/mL)",
  "CV pico (UI/mL)",
  "Genotipos mixtos gB/gH",
  "Enfermedad sintomática",
  "Profilaxis al inicio del episodio",
  "Tratamiento antiviral",
  "Hospitalización cualquier causa",
  "Hospitalización por CMV"
)

grupos_map <- c(
  "Edad receptor (años)"              = "Demográficas",
  "Tipo de trasplante (BI vs UNI)"    = "Demográficas",
  "Enfermedad de base"                = "Demográficas",
  "Perfil de riesgo serológico D/R"   = "Demográficas",   
  "Tiempo desde TXP (meses)"          = "Cinéticas",
  "CV inicial (UI/mL)"                = "Cinéticas",
  "CV pico (UI/mL)"                   = "Cinéticas",
  "Genotipos mixtos gB/gH"            = "Cinéticas",
  "Enfermedad sintomática"            = "Clínicas",
  "Profilaxis al inicio del episodio" = "Clínicas",
  "Tratamiento antiviral"             = "Clínicas",
  "Hospitalización cualquier causa"   = "Clínicas",
  "Hospitalización por CMV"           = "Clínicas"
)

res_all <- bind_rows(res_cor, res_cat, res_eb) %>%
  mutate(
    sig     = p_val < 0.05,
    p_fmt   = case_when(
      p_val < 0.001 ~ "<0,001",
      TRUE          ~ sub("\\.", ",", sprintf("%.3f", p_val))
    ),
    p_label = paste0("p = ", p_fmt),
    grupo   = factor(grupos_map[as.character(Variable)],
                     levels = c("Demográficas", "Cinéticas", "Clínicas")),
    Variable = factor(Variable, levels = rev(orden_vars))
  )

# ── 6. Figura ─────────────────────────────────────────────────────────────────
figura_univariante <- ggplot(res_all, aes(x = p_val, y = Variable)) +
  geom_vline(
    xintercept = 0.05,
    color = "firebrick3",
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0.20, color = "grey65",
             linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(color = sig, size = sig)) +
  geom_text(aes(label = p_label, color = sig),
            hjust = -0.15, size = 4) +
  geom_text(aes(x = 1.1, label = label_est),
            hjust = 1, size = 4, color = "grey40") +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    name   = "valor p",
    limits = c(0, 1.15),
    breaks = c(0.05, 0.20, 0.50, 1.00),
    labels = c("0,05", "0,20", "0,50", "1,00")
  ) +
  scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey50")) +
  scale_size_manual(values  = c("TRUE" = 3.5,          "FALSE" = 2.5)) +
  labs(y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.y        = element_text(size = 12, color = "grey20"),
    axis.text.x        = element_text(size = 12, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20"),
    strip.text.y       = element_text(angle = -90, size = 12,
                                      face = "bold", color = "grey20",
                                      margin = margin(l = 6, r = 6)),
    strip.background   = element_rect(fill = "grey93", color = NA),
    panel.spacing.y    = unit(0.6, "lines"),
    plot.margin        = margin(8, 60, 8, 8)
  )

figura_univariante


## ---- [5g] Tabla 15. Modelo multivariante del aclaramiento en genotipos mixtos ----
# Apartado 5.5.4. Regresion lineal multiple y VIF.

# Variables candidatas (p < 0,20 en univariante)
df_multi <- df_uni %>%
  filter(
    !is.na(aclaramiento_dias),
    !is.na(gbghmix_f),
    !is.na(risk_f),
    !is.na(cvini),
    !is.na(cvpico),
    !is.na(Tdesdetx_meses),
    !is.na(sintom_f),
    !is.na(profilaxis_f),
    !is.na(ing_cmv_f)
  )

modelo_completo <- lm(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    Tdesdetx_meses + sintom_f + profilaxis_f + ing_cmv_f,
  data = df_multi
)

summary(modelo_completo)


# ── Selección backward ────────────────────────────────────────────────────────
modelo_backward <- step(modelo_completo, direction = "backward", trace = TRUE)
summary(modelo_backward)


# ── Intervalos de confianza y VIF ─────────────────────────────────────────────
library(car)

confint(modelo_backward, level = 0.95)
vif(modelo_backward)


confint(modelo_backward, level = 0.95)


## ---- [5h] Figura 29. Efecto de gB3 en el contexto de genotipos mixtos ----
# Apartado 5.5.4. Tiempo de aclaramiento estratificado.

## ---- [i]  Análisis gB3 en contexto de genotipos mixtos ----------------------

# ── 1. Crear variables de clasificación ───────────────────────────────────────
df_gb3 <- df_uni %>%
  mutate(
    # gB3 presente en el episodio (independientemente de si es único o mixto)
    gb3_presente = as.integer(gb3 == 1),
    
    # Cuatro grupos: único/mixto x gB3 sí/no
    grupo_gb3 = case_when(
      gbmix == 0 & gb3 == 0 ~ "Único sin gB3",
      gbmix == 0 & gb3 == 1 ~ "Único con gB3",
      gbmix == 1 & gb3 == 0 ~ "Mixto sin gB3",
      gbmix == 1 & gb3 == 1 ~ "Mixto con gB3",
      TRUE ~ NA_character_
    ),
    grupo_gb3 = factor(grupo_gb3,
                       levels = c("Único sin gB3", "Único con gB3",
                                  "Mixto sin gB3", "Mixto con gB3"))
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(grupo_gb3))

# ── 2. Análisis A: entre episodios mixtos, gB3 vs no gB3 ─────────────────────
df_mixtos <- df_gb3 %>% filter(gbmix == 1)

cat("\n── Análisis A: episodios mixtos, gB3 vs no gB3 ──\n")
df_mixtos %>%
  group_by(gb3_presente) %>%
  summarise(
    n        = n(),
    mediana  = median(aclaramiento_dias, na.rm = TRUE),
    q1       = quantile(aclaramiento_dias, 0.25, na.rm = TRUE),
    q3       = quantile(aclaramiento_dias, 0.75, na.rm = TRUE)
  ) %>% print()

wilcox.test(aclaramiento_dias ~ gb3_presente, data = df_mixtos, exact = FALSE)

# ── 3. Análisis B: cuatro grupos ──────────────────────────────────────────────
cat("\n── Análisis B: cuatro grupos ──\n")
df_gb3 %>%
  group_by(grupo_gb3) %>%
  summarise(
    n        = n(),
    mediana  = median(aclaramiento_dias, na.rm = TRUE),
    q1       = quantile(aclaramiento_dias, 0.25, na.rm = TRUE),
    q3       = quantile(aclaramiento_dias, 0.75, na.rm = TRUE)
  ) %>% print()

kruskal.test(aclaramiento_dias ~ grupo_gb3, data = df_gb3)

# Post-hoc si Kruskal significativo
library(dunn.test)
dunn.test(df_gb3$aclaramiento_dias, df_gb3$grupo_gb3,
          method = "bonferroni", kw = TRUE)

etiquetas_n <- c(
  "Único sin gB3" = paste0("Único sin gB3\n(n=", sum(df_gb3$grupo_gb3 == "Único sin gB3", na.rm=TRUE), ")"),
  "Único con gB3" = paste0("Único con gB3\n(n=", sum(df_gb3$grupo_gb3 == "Único con gB3", na.rm=TRUE), ")"),
  "Mixto sin gB3" = paste0("Mixto sin gB3\n(n=", sum(df_gb3$grupo_gb3 == "Mixto sin gB3", na.rm=TRUE), ")"),
  "Mixto con gB3" = paste0("Mixto con gB3\n(n=", sum(df_gb3$grupo_gb3 == "Mixto con gB3", na.rm=TRUE), ")")
)

figura_gb3 <- ggplot(
  df_gb3,
  aes(
    x = grupo_gb3,
    y = aclaramiento_dias,
    fill = ifelse(grepl("con gB3", grupo_gb3), "con gB3", "sin gB3")
  )
) +
  geom_boxplot(
    alpha         = 0.75,
    outlier.shape = NA,      # oculta outliers nativos del boxplot
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4
  ) +
  # Puntos individuales solo hasta 150
  geom_jitter(
    data        = df_gb3 %>% filter(aclaramiento_dias <= 150),
    aes(x = grupo_gb3, y = aclaramiento_dias,
        fill = ifelse(grepl("con gB3", grupo_gb3), "con gB3", "sin gB3")),
    shape       = 21,
    size        = 1.5,
    color       = "grey40",
    alpha       = 0.6,
    width       = 0.15,
    height      = 0,
    inherit.aes = FALSE
  ) + 
  geom_vline(
    xintercept = 2.5,
    linetype   = "dashed",
    color      = "grey70",
    linewidth  = 0.4
  ) +
  # Barras de significación
  annotate("segment", x = 1,   xend = 2, y = 170, yend = 170, linewidth = 0.4) +
  annotate("segment", x = 1,   xend = 1, y = 165, yend = 170, linewidth = 0.4) +
  annotate("segment", x = 2,   xend = 2, y = 165, yend = 170, linewidth = 0.4) +
  annotate("text",    x = 1.5, y = 174,  label = "***",  size = 5) +
  annotate("segment", x = 1,   xend = 4, y = 185, yend = 185, linewidth = 0.4) +
  annotate("segment", x = 1,   xend = 1, y = 180, yend = 185, linewidth = 0.4) +
  annotate("segment", x = 4,   xend = 4, y = 180, yend = 185, linewidth = 0.4) +
  annotate("text",    x = 2.5, y = 189,  label = "****", size = 5) +
  # Etiquetas de bloque
  annotate("text", x = 1.5, y = 208, label = "Genotipo único",
           size = 4, color = "grey30", fontface = "bold") +
  annotate("text", x = 3.5, y = 208, label = "Genotipo mixto",
           size = 4, color = "grey30", fontface = "bold") +
  scale_fill_manual(
    values = c(
      "sin gB3" = "#CDC0B0",
      "con gB3" = "#4F94CD"
    )
  ) +
  scale_x_discrete(labels = etiquetas_n) +
  scale_y_continuous(
    name   = "Tiempo hasta aclaramiento viral (días)",
    breaks = seq(0, 150, by = 50)
  ) +
  coord_cartesian(ylim = c(0, 215), clip = "off") +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.x        = element_text(size = 11, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey40"),
    axis.title.y       = element_text(size = 11, color = "grey20"),
    plot.margin        = margin(8, 8, 8, 8)
  )

figura_gb3


## ---- [6a1] Variabilidad intraepisodio: preparacion de los datos ----
# Apartado 5.6.1. Clasificacion de episodios estables y variables.

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(flextable)
library(officer)

# ── 1. Cargar muestras ────────────────────────────────────────────────────────
muestras <- read_excel("data/base_clinica.xlsx",
                       sheet = "genotipado",
                       col_types = c("numeric", "numeric", "date", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "text"))

# ── 1b. Filtrar solo muestras con genotipado completo ─────────────────────────
muestras <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE)
  ) %>%
  filter(gb >= 1 & gh >= 1)

cat(sprintf("Muestras con genotipado completo: %d\n", nrow(muestras)))
cat(sprintf("Pacientes: %d\n", n_distinct(muestras$NHC)))

# ── 2. Asignar cada muestra a su episodio ─────────────────────────────────────
muestras_ep <- muestras %>%
  inner_join(
    dnaemias %>%
      select(NHC, Fini, Ffin, gbghmix) %>%
      mutate(Fini = as.Date(Fini), Ffin = as.Date(Ffin)),
    by = "NHC",
    relationship = "many-to-many"
  ) %>%
  filter(
    Freg >= Fini &
      (is.na(Ffin) | Freg <= Ffin)
  )

cat(sprintf("Muestras asignadas: %d\n", nrow(muestras_ep)))
cat(sprintf("Episodios representados: %d\n",
            n_distinct(muestras_ep %>% select(NHC, Fini))))

# ── 3. Episodios con ≥2 muestras ──────────────────────────────────────────────
ep_n_muestras <- muestras_ep %>%
  group_by(NHC, Fini, Ffin, gbghmix) %>%
  summarise(n_muestras = n(), .groups = "drop")

n_total_ep <- 234
n_intra_ep <- ep_n_muestras %>% filter(n_muestras >= 2) %>% nrow()

cat(sprintf("Episodios con ≥2 muestras: %d (%.1f%%)\n",
            n_intra_ep, 100 * n_intra_ep / n_total_ep))

ep_intra <- ep_n_muestras %>% filter(n_muestras >= 2)

cat(sprintf("Mediana muestras por episodio: %.0f (RIQ %.0f-%.0f)\n",
            median(ep_intra$n_muestras),
            quantile(ep_intra$n_muestras, 0.25),
            quantile(ep_intra$n_muestras, 0.75)))

# ── 4. Clasificar episodios únicos vs mixtos ──────────────────────────────────
cat("\n── Clasificación episodios con ≥2 muestras ──\n")
ep_intra %>%
  mutate(clase = ifelse(gbghmix == 0, "Único", "Mixto")) %>%
  count(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 5. Variabilidad intraepisodio ─────────────────────────────────────────────
muestras_intra <- muestras_ep %>%
  semi_join(ep_intra, by = c("NHC", "Fini")) %>%
  mutate(patron = paste(GB1, GB2, GB3, GB4, GH1, GH2, sep = "-"))

variabilidad_ep <- muestras_intra %>%
  group_by(NHC, Fini, gbghmix) %>%
  summarise(
    n_muestras = n(),
    n_patrones = n_distinct(patron),
    .groups    = "drop"
  ) %>%
  mutate(
    clase    = ifelse(gbghmix == 0, "Único", "Mixto"),
    variable = ifelse(n_patrones > 1, "Variable", "Estable")
  )

cat("\n── Variabilidad por clase ──\n")
variabilidad_ep %>%
  group_by(clase, variable) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 6. Episodios mixtos variables ─────────────────────────────────────────────
mixtos_var <- variabilidad_ep %>%
  filter(clase == "Mixto", variable == "Variable")

cat(sprintf("\n── Episodios mixtos variables: %d episodios, %d pacientes ──\n",
            nrow(mixtos_var), n_distinct(mixtos_var$NHC)))

mixtos_var %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  distinct(NHC, serDR) %>%
  count(serDR) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()


## ---- [6a2] Tabla 16. Cambios intraepisodio en episodios con genotipos mixtos ----
# Apartado 5.6.1.

## ---- [b] Análisis detallado episodios mixtos variables ----------------------
muestras_mv <- muestras_ep %>%
  semi_join(mixtos_var, by = c("NHC", "Fini")) %>%
  mutate(patron = paste(GB1, GB2, GB3, GB4, GH1, GH2, sep = "-")) %>%
  arrange(NHC, Fini, Freg)

# Número de cambios por episodio
cambios_ep <- muestras_mv %>%
  group_by(NHC, Fini) %>%
  arrange(Freg, .by_group = TRUE) %>%
  mutate(cambio = patron != lag(patron)) %>%
  summarise(
    n_muestras = n(),
    n_cambios  = sum(cambio, na.rm = TRUE),
    .groups    = "drop"
  )

cat("── Número de cambios por episodio ──\n")
cat(sprintf("Rango: %d – %d\n", min(cambios_ep$n_cambios), max(cambios_ep$n_cambios)))
cat(sprintf("Mediana: %.0f (RIQ %.0f–%.0f)\n",
            median(cambios_ep$n_cambios),
            quantile(cambios_ep$n_cambios, 0.25),
            quantile(cambios_ep$n_cambios, 0.75)))

# Estabilidad por genotipo
genotipos <- c("GB1","GB2","GB3","GB4","GH1","GH2")
nombres   <- c("gB1","gB2","gB3","gB4","gH1","gH2")

estabilidad_gt <- muestras_mv %>%
  group_by(NHC, Fini) %>%
  summarise(
    across(all_of(genotipos), ~ {
      presente <- sum(. == 1, na.rm = TRUE)
      n_tot    <- n()
      if (presente == 0) "ausente"
      else if (presente == n_tot) "estable"
      else "variable"
    }),
    .groups = "drop"
  )

cat("\n── Frecuencia genotipos estables ──\n")
for (i in seq_along(genotipos)) {
  g <- genotipos[i]
  n_est <- sum(estabilidad_gt[[g]] == "estable")
  cat(sprintf("%s estable: %d/35 (%.1f%%)\n", nombres[i], n_est, 100*n_est/35))
}

cat("\n── Frecuencia genotipos variables ──\n")
for (i in seq_along(genotipos)) {
  g <- genotipos[i]
  n_var <- sum(estabilidad_gt[[g]] == "variable")
  cat(sprintf("%s variable: %d/35 (%.1f%%)\n", nombres[i], n_var, 100*n_var/35))
}

estabilidad_gt <- estabilidad_gt %>%
  mutate(tiene_estable = rowSums(across(all_of(genotipos), ~ . == "estable")) > 0)

cat(sprintf("\n── Episodios con ≥1 genotipo estable: %d/35 (%.1f%%)\n",
            sum(estabilidad_gt$tiene_estable),
            100 * mean(estabilidad_gt$tiene_estable)))

cat(sprintf("── Episodios con cambio completo: %d\n",
            sum(!estabilidad_gt$tiene_estable)))

tabla_mv <- estabilidad_gt %>%
  select(-any_of("n_muestras")) %>%
  left_join(
    cambios_ep %>% select(NHC, Fini, n_cambios, n_muestras),
    by = c("NHC", "Fini")
  ) %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  mutate(
    estables = pmap_chr(
      list(GB1, GB2, GB3, GB4, GH1, GH2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "estable") "gB1",
          if (gb2 == "estable") "gB2",
          if (gb3 == "estable") "gB3",
          if (gb4 == "estable") "gB4",
          if (gh1 == "estable") "gH1",
          if (gh2 == "estable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    ),
    variables = pmap_chr(
      list(GB1, GB2, GB3, GB4, GH1, GH2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "variable") "gB1",
          if (gb2 == "variable") "gB2",
          if (gb3 == "variable") "gB3",
          if (gb4 == "variable") "gB4",
          if (gh1 == "variable") "gH1",
          if (gh2 == "variable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    )
  ) %>%
  arrange(NHC, Fini) %>%
  mutate(
    pac_num   = as.integer(factor(NHC, levels = unique(NHC))),
    ep_num    = row_number(),
    pac_label = paste0(pac_num, " ", serDR)
  ) %>%
  select(ep_num, pac_label, n_muestras, estables, variables, n_cambios)

tabla_ft <- tabla_mv %>%
  rename(
    "Episodio"              = ep_num,
    "Paciente"              = pac_label,
    "N muestras"            = n_muestras,
    "Patrón dominante"      = estables,
    "Genotipos fluctuantes" = variables,
    "N cambios"             = n_cambios
  )

ft <- flextable(tabla_ft) %>%
  width(j = "Episodio",              width = 0.7) %>%
  width(j = "Paciente",              width = 1.6) %>%
  width(j = "N muestras",            width = 0.9) %>%
  width(j = "Patrón dominante",      width = 2.0) %>%
  width(j = "Genotipos fluctuantes", width = 2.0) %>%
  width(j = "N cambios",             width = 0.9) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Patrón dominante", "Genotipos fluctuantes"),
        align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = 1, part = "body") %>%
  color(i = 1, color = "firebrick3", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(
    i = which(diff(as.integer(factor(tabla_ft$Paciente))) != 0),
    border = fp_border(color = "grey50", width = 0.5),
    part = "body"
  ) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "Patrón dominante: genotipos presentes en todas las muestras del episodio.",
    "Genotipos fluctuantes: genotipos con aparición o desaparición entre muestras consecutivas.",
    "El episodio en rojo corresponde al único caso con cambio completo del patrón dominante.",
    "N cambios: número de transiciones de patrón genotípico entre muestras consecutivas."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(caption = "Tabla R8. Episodios con variabilidad intraepisodio en genotipos mixtos gB/gH.") %>%
  set_table_properties(layout = "autofit")

ft


## ---- [6b] Figura 31. Distancia al P75 del Ct en genotipos estables y fluctuantes ----
# Apartado 5.6.1.

## ---- [d] *Figura 21: Distancia al P75 estables vs fluctuantes ----------------
ct_cols  <- c("GB1CT","GB2CT","GB3CT","GB4CT","GH1CT","GH2CT")
gt_cols  <- c("GB1","GB2","GB3","GB4","GH1","GH2")
gt_names <- c("gB1","gB2","gB3","gB4","gH1","gH2")

gt_largo <- muestras_mv %>%
  select(NHC, Fini, Freg, all_of(gt_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "presente") %>%
  filter(presente == 1) %>%
  mutate(genotipo = gt_names[match(genotipo_var, gt_cols)])

ct_largo <- muestras_mv %>%
  select(NHC, Fini, Freg, all_of(ct_cols)) %>%
  pivot_longer(cols = all_of(ct_cols),
               names_to = "ct_var", values_to = "Ct") %>%
  mutate(genotipo_var = gt_cols[match(ct_var, ct_cols)])

det_largo <- gt_largo %>%
  left_join(ct_largo %>% select(NHC, Fini, Freg, genotipo_var, Ct),
            by = c("NHC", "Fini", "Freg", "genotipo_var")) %>%
  filter(!is.na(Ct))

estabilidad_largo <- estabilidad_gt %>%
  select(NHC, Fini, all_of(gt_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "estabilidad") %>%
  filter(estabilidad != "ausente")

det_largo <- det_largo %>%
  left_join(estabilidad_largo, by = c("NHC", "Fini", "genotipo_var")) %>%
  filter(!is.na(estabilidad))

# P75 por genotipo sobre todas las muestras
todas_ct <- muestras %>%
  select(NHC, Freg, all_of(gt_cols), all_of(ct_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "presente") %>%
  filter(presente == 1) %>%
  left_join(
    muestras %>%
      select(NHC, Freg, all_of(ct_cols)) %>%
      pivot_longer(cols = all_of(ct_cols),
                   names_to = "ct_var", values_to = "Ct") %>%
      mutate(genotipo_var = gt_cols[match(ct_var, ct_cols)]) %>%
      select(NHC, Freg, genotipo_var, Ct),
    by = c("NHC", "Freg", "genotipo_var")
  ) %>%
  filter(!is.na(Ct))

p75_genotipo <- todas_ct %>%
  group_by(genotipo_var) %>%
  summarise(P75 = quantile(Ct, 0.75, na.rm = TRUE), .groups = "drop")

det_largo <- det_largo %>%
  left_join(p75_genotipo, by = "genotipo_var") %>%
  mutate(dist_p75 = P75 - Ct)

cat("\n── Distancia al P75: estable vs fluctuante ──\n")
det_largo %>%
  group_by(estabilidad) %>%
  summarise(
    n             = n(),
    mediana       = median(dist_p75, na.rm = TRUE),
    q1            = quantile(dist_p75, 0.25, na.rm = TRUE),
    q3            = quantile(dist_p75, 0.75, na.rm = TRUE),
    pct_sobre_p75 = 100 * mean(dist_p75 < 0, na.rm = TRUE)
  ) %>% print()

wilcox.test(dist_p75 ~ estabilidad, data = det_largo, exact = FALSE)

etiq <- det_largo %>%
  group_by(estabilidad) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(ifelse(estabilidad == "estable",
                               "Genotipos estables", "Genotipos fluctuantes"),
                        "\n(n=", n, ")"))

det_fig <- det_largo %>%
  mutate(
    grupo = ifelse(estabilidad == "estable",
                   "Patrón estable", "Genotipos fluctuantes"),
    grupo = factor(grupo, levels = c("Patrón estable", "Genotipos fluctuantes"))
  )

etiq_x <- etiq %>%
  mutate(grupo = ifelse(estabilidad == "estable",
                        "Patrón estable", "Genotipos fluctuantes"),
         grupo = factor(grupo, levels = c("Patrón estable", "Genotipos fluctuantes")))

figura_R22 <- ggplot(det_fig, aes(x = grupo, y = dist_p75, fill = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "firebrick3", linewidth = 0.6) +
  annotate("text", x = 2.42, y = 0.8, label = "P75",
           size = 5, color = "firebrick3", fontface = "bold") +
  geom_jitter(aes(color = grupo), width = 0.15, size = 0.9,
              alpha = 0.35, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.5,
               color = "grey15", linewidth = 0.7) +
  # Mediana resaltada en blanco
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, fatten = 0,
               color = "white", linewidth = 0.7) +
  annotate("segment", x = 1, xend = 2, y = 9,   yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = 8.5, yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 8.5, yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = 9.7, label = "***",
           size = 7, color = "black", fontface = "bold") +
  scale_x_discrete(labels = setNames(etiq_x$label, etiq_x$grupo)) +
  scale_fill_manual(
    values = c("Patrón estable"        = "#2E7D32",
               "Genotipos fluctuantes" = "#E65100")
  ) +
  scale_color_manual(
    values = c("Patrón estable"        = "#2E7D32",
               "Genotipos fluctuantes" = "#E65100")
  ) +
  scale_y_continuous(
    name   = "Distancia al P75 del genotipo (ciclos Ct)",
    breaks = seq(-6, 10, by = 2)
  ) +
  coord_cartesian(ylim = c(-7, 11)) +
  labs(x = NULL) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.line       = element_line(color = "black", linewidth = 0.6),
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    axis.text.x     = element_text(size = 14, color = "black"),
    axis.text.y     = element_text(size = 13, color = "black"),
    axis.title.y    = element_text(size = 14, color = "black", face = "bold"),
    plot.margin     = margin(8, 8, 8, 8)
  )

figura_R22


## ---- [6c1] Variabilidad interepisodio: clasificacion de los patrones ----
# Apartado 5.6.2.

library(dplyr)
library(tidyr)
library(purrr)
library(flextable)
library(officer)

# ── 1. Pacientes con >1 episodio genotipado ───────────────────────────────────
ep_por_paciente <- dnaemias %>%
  filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  group_by(NHC) %>%
  summarise(n_episodios = n(), .groups = "drop") %>%
  filter(n_episodios > 1)

cat(sprintf("Pacientes con >1 episodio: %d\n", nrow(ep_por_paciente)))
cat(sprintf("Total episodios: %d\n", sum(ep_por_paciente$n_episodios)))

cat("\n── Distribución por número de episodios ──\n")
ep_por_paciente %>%
  count(n_episodios) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 2. Episodios de pacientes con >1 episodio ─────────────────────────────────
ep_multi <- dnaemias %>%
  filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  semi_join(ep_por_paciente, by = "NHC")

cat("\n── Clasificación único/mixto ──\n")
ep_multi %>%
  mutate(clase = ifelse(gbghmix == 0, "Único", "Mixto")) %>%
  count(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 3. Clasificar pacientes según comportamiento interepisodio ────────────────
patron_paciente <- ep_multi %>%
  group_by(NHC) %>%
  summarise(
    n_episodios = n(),
    n_unicos    = sum(gbghmix == 0),
    n_mixtos    = sum(gbghmix == 1),
    .groups     = "drop"
  ) %>%
  mutate(
    patron = case_when(
      n_unicos > 0 & n_mixtos == 0 ~ "Siempre único",
      n_unicos == 0 & n_mixtos > 0 ~ "Siempre mixto",
      TRUE                          ~ "Mixtos y únicos"
    )
  )

cat("\n── Clasificación pacientes ──\n")
patron_paciente %>%
  count(patron) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

cat("\n── Episodios por grupo ──\n")
patron_paciente %>%
  group_by(patron) %>%
  summarise(n_episodios = sum(n_episodios), .groups = "drop") %>%
  mutate(pct = 100 * n_episodios / sum(n_episodios)) %>%
  print()

# ── 4. Perfil serológico por grupo ────────────────────────────────────────────
cat("\n── Perfil serológico por grupo ──\n")
patron_paciente %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  group_by(patron, serDR) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patron) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 5. Tiempo entre episodios por grupo ───────────────────────────────────────
cat("\n── Tiempo entre episodios (Ffin → Fini siguiente, días) ──\n")
ep_multi %>%
  mutate(Ffin = as.Date(Ffin)) %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>%
  mutate(tiempo_entre = as.numeric(Fini - lag(Ffin))) %>%
  filter(!is.na(tiempo_entre), tiempo_entre >= 0) %>%
  left_join(patron_paciente %>% select(NHC, patron), by = "NHC") %>%
  group_by(patron) %>%
  summarise(
    mediana = median(tiempo_entre, na.rm = TRUE),
    q1      = quantile(tiempo_entre, 0.25, na.rm = TRUE),
    q3      = quantile(tiempo_entre, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ── 6. Análisis genotipos estables/variables interepisodio ────────────────────
pacientes_cambio <- patron_paciente %>%
  filter(patron == "Mixtos y únicos") %>%
  pull(NHC)

ep_cambio <- ep_multi %>%
  filter(NHC %in% pacientes_cambio) %>%
  arrange(NHC, Fini)

genotipos <- c("gb1","gb2","gb3","gb4","gh1","gh2")
nombres   <- c("gB1","gB2","gB3","gB4","gH1","gH2")

estab_inter <- ep_cambio %>%
  group_by(NHC) %>%
  summarise(
    n_ep = n(),
    across(all_of(genotipos), ~ {
      presente <- sum(. == 1, na.rm = TRUE)
      n_tot    <- n()
      if (presente == 0) "ausente"
      else if (presente == n_tot) "estable"
      else "variable"
    }),
    .groups = "drop"
  ) %>%
  mutate(
    tiene_estable = rowSums(across(all_of(genotipos), ~ . == "estable")) > 0
  )

cat(sprintf("\n── Pacientes con ≥1 genotipo estable: %d/21 (%.1f%%)\n",
            sum(estab_inter$tiene_estable),
            100 * mean(estab_inter$tiene_estable)))

cat(sprintf("── Pacientes sin ningún genotipo estable: %d\n",
            sum(!estab_inter$tiene_estable)))

cat("\n── Frecuencia genotipos estables ──\n")
for (i in seq_along(genotipos)) {
  n_est <- sum(estab_inter[[genotipos[i]]] == "estable")
  cat(sprintf("%s estable: %d/21 (%.1f%%)\n", nombres[i], n_est, 100*n_est/21))
}

cat("\n── Frecuencia genotipos variables ──\n")
for (i in seq_along(genotipos)) {
  n_var <- sum(estab_inter[[genotipos[i]]] == "variable")
  cat(sprintf("%s variable: %d/21 (%.1f%%)\n", nombres[i], n_var, 100*n_var/21))
}

# ── 7. Número de cambios entre episodios consecutivos ─────────────────────────
cambios_inter <- ep_cambio %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>%
  mutate(
    patron_ep = paste(gb1, gb2, gb3, gb4, gh1, gh2, sep = "-"),
    cambio    = patron_ep != lag(patron_ep)
  ) %>%
  summarise(
    n_cambios = sum(cambio, na.rm = TRUE),
    .groups   = "drop"
  )


## ---- [6c2] Tabla 17. Pacientes con cambios genotipicos interepisodio ----
# Apartado 5.6.2.

# ── 8. Construir tabla ────────────────────────────────────────────────────────
tabla_inter <- estab_inter %>%
  left_join(cambios_inter, by = "NHC") %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  mutate(
    estables = pmap_chr(
      list(gb1, gb2, gb3, gb4, gh1, gh2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "estable") "gB1",
          if (gb2 == "estable") "gB2",
          if (gb3 == "estable") "gB3",
          if (gb4 == "estable") "gB4",
          if (gh1 == "estable") "gH1",
          if (gh2 == "estable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    ),
    variables = pmap_chr(
      list(gb1, gb2, gb3, gb4, gh1, gh2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "variable") "gB1",
          if (gb2 == "variable") "gB2",
          if (gb3 == "variable") "gB3",
          if (gb4 == "variable") "gB4",
          if (gh1 == "variable") "gH1",
          if (gh2 == "variable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    )
  ) %>%
  arrange(desc(n_ep), desc(n_cambios)) %>%
  mutate(pac_num = row_number()) %>%
  select(pac_num, serDR, n_ep, estables, variables, n_cambios)

# ── 9. Flextable ──────────────────────────────────────────────────────────────
tabla_inter_ft <- tabla_inter %>%
  rename(
    "Paciente"              = pac_num,
    "Serología D/R"         = serDR,
    "N episodios"           = n_ep,
    "Patrón dominante"      = estables,
    "Genotipos fluctuantes" = variables,
    "N cambios"             = n_cambios
  )

fila_sin_estable <- which(tabla_inter_ft$`Patrón dominante` == "—")

ft_inter <- flextable(tabla_inter_ft) %>%
  width(j = "Paciente",              width = 0.7) %>%
  width(j = "Serología D/R",         width = 1.2) %>%
  width(j = "N episodios",           width = 0.9) %>%
  width(j = "Patrón dominante",      width = 2.0) %>%
  width(j = "Genotipos fluctuantes", width = 2.2) %>%
  width(j = "N cambios",             width = 0.9) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Patrón dominante", "Genotipos fluctuantes"),
        align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = fila_sin_estable, part = "body") %>%
  color(i = fila_sin_estable, color = "firebrick3", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "Patrón dominante: genotipos presentes en todos los episodios del paciente.",
    "Genotipos fluctuantes: genotipos con aparición o desaparición entre episodios.",
    "El paciente en rojo corresponde al único caso sin ningún genotipo estable entre episodios.",
    "N cambios: número de transiciones de patrón genotípico entre episodios consecutivos."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R9. Patrón de cambio de genotipos gB/gH en pacientes con infecciones múltiples que presentaron variabilidad interepisodio."
  ) %>%
  set_table_properties(layout = "autofit")

ft_inter


## ---- [6c3] Figura 32. Variabilidad genotipica interepisodio ----
# Apartado 5.6.2. Diagrama de flujo de la clasificacion.

library(DiagrammeR)

grViz("
digraph flujo {

  graph [layout = dot, rankdir = TB,
         fontname = 'Arial', splines = ortho,
         nodesep = 0.6, ranksep = 0.8]

  node [shape = rectangle, style = filled,
        fontname = 'Arial', fontsize = 19,
        width = 2.0, height = 0.7, margin = '0.18,0.12']

  edge [color = '#999999', arrowsize = 0.5,
        arrowhead = vee]

subgraph cluster_A {
  label = < <B>A. Variabilidad intraepisodio</B> >
  fontsize = 25
  labeljust = 'c'
  margin = 40
  color = white

    A1 [label = '234 episodios\ngenotipado gB/gH',     fillcolor = '#607D8B', fontcolor = white]
    A2 [label = '131 (56,0%)\n1 sola muestra',         fillcolor = '#90A4AE', fontcolor = white]
    A3 [label = '103 (44,0%)\n≥2 muestras',            fillcolor = '#607D8B', fontcolor = white]
    A4 [label = '61 (59,2%)\ngenotipo único',           fillcolor = '#1565C0', fontcolor = white]
    A5 [label = '42 (40,8%)\ngenotipo mixto',           fillcolor = '#BF360C', fontcolor = white]
    A6 [label = '61 (100%)\npatrón estable',            fillcolor = '#2E7D32', fontcolor = white]
    A7 [label = '7 (16,7%)\nestable',                   fillcolor = '#2E7D32', fontcolor = white]
    A8 [label = '35 (83,3%)\nfluctuante',               fillcolor = '#E65100', fontcolor = white]
    A9 [label = '34 (97,1%)\n≥1 dominante',             fillcolor = '#2E7D32', fontcolor = white]
    A10[label = '1 (2,9%)\ncambio completo',            fillcolor = '#7B1FA2', fontcolor = white]

    A1 -> A2
    A1 -> A3
    A3 -> A4
    A3 -> A5
    A4 -> A6
    A5 -> A7
    A5 -> A8
    A8 -> A9
    A8 -> A10

    {rank = same; A2; A3}
    {rank = same; A4; A5}
    {rank = same; A6; A7; A8}
    {rank = same; A9; A10}
  }

  subgraph cluster_B {
    label = < <B>B. Variabilidad interepisodio</B> >
    fontsize = 25
    labeljust = 'c'
      margin = 40
    color = white

    B1 [label = '146 pacientes\ngenotipado gB/gH',     fillcolor = '#607D8B', fontcolor = white]
    B2 [label = '90 (61,6%)\n1 solo episodio',         fillcolor = '#90A4AE', fontcolor = white]
    B3 [label = '56 (38,4%)\n>1 episodio',             fillcolor = '#607D8B', fontcolor = white]
    B4 [label = '35 (62,5%)\nmantienen patrón',        fillcolor = '#1565C0', fontcolor = white]
    B5 [label = '21 (37,5%)\ncambian patrón',          fillcolor = '#BF360C', fontcolor = white]
    B6 [label = '26 (74,3%)\nsiempre único',           fillcolor = '#1565C0', fontcolor = white]
    B7 [label = '9 (25,7%)\nsiempre mixto',            fillcolor = '#BF360C', fontcolor = white]
    B8 [label = '20 (95,2%)\n≥1 dominante',            fillcolor = '#2E7D32', fontcolor = white]
    B9 [label = '1 (4,8%)\ncambio completo',           fillcolor = '#7B1FA2', fontcolor = white]

    B1 -> B2
    B1 -> B3
    B3 -> B4
    B3 -> B5
    B4 -> B6
    B4 -> B7
    B5 -> B8
    B5 -> B9

    {rank = same; B2; B3}
    {rank = same; B4; B5}
    {rank = same; B6; B7; B8; B9}
  }
}
")


## ---- [7a] Tabla 19. Caracteristicas basales segun genotipos mixtos gB/gH ----
# Apartado 5.7.1.

library(dplyr)
library(lubridate)
library(readxl)
library(flextable)
library(officer)

pacientes <- read_excel("data/base_clinica.xlsx",
                        sheet = "pacientes_sel")

# ── 1. Preparar datos ─────────────────────────────────────────────────────────
pacientes_tab <- pacientes %>%
  mutate(
    Ftx         = as.Date(Ftx),
    Ffinseg     = as.Date(Ffinseg),
    t_seg_anios = as.numeric(Ffinseg - Ftx) / 365.25,
    grupo       = factor(gbghm, levels = c(0,1),
                         labels = c("Genotipos únicos", "Genotipos mixtos"))
  )

g0  <- pacientes_tab %>% filter(grupo == "Genotipos únicos")
g1  <- pacientes_tab %>% filter(grupo == "Genotipos mixtos")
tot <- pacientes_tab
n0  <- nrow(g0); n1 <- nrow(g1); nt <- nrow(tot)

# ── 2. Funciones auxiliares ───────────────────────────────────────────────────
med_ric <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f (%.1f\u2013%.1f)",
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}

n_pct <- function(x, n_total) {
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE),
          100 * sum(x, na.rm = TRUE) / n_total)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0,001",
         sub("\\.", ",", sprintf("%.3f", p)))
}

p_cont <- function(x, g) {
  fmt_p(wilcox.test(x ~ g, exact = FALSE)$p.value)
}

p_cat <- function(x, g) {
  tab <- table(x, g)
  test <- tryCatch(
    fisher.test(tab, simulate.p.value = TRUE, B = 10000),
    error = function(e) chisq.test(tab)
  )
  fmt_p(test$p.value)
}

# ── 3. Construir filas ────────────────────────────────────────────────────────
filas <- bind_rows(
  
  tibble(variable = "n",
         total = as.character(nt),
         unicos = as.character(n0),
         mixtos = as.character(n1),
         p = ""),
  
  # ── Demográficas ──────────────────────────────────────────────────────────
  tibble(variable = "CARACTERÍSTICAS DEMOGRÁFICAS",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "Edad receptor (años), mediana (RIC)",
         total  = med_ric(tot$edad_tx),
         unicos = med_ric(g0$edad_tx),
         mixtos = med_ric(g1$edad_tx),
         p = p_cont(pacientes_tab$edad_tx, pacientes_tab$grupo)),
  
  tibble(variable = "Edad donante (años), mediana (RIC)",
         total  = med_ric(tot$edad_don),
         unicos = med_ric(g0$edad_don),
         mixtos = med_ric(g1$edad_don),
         p = p_cont(pacientes_tab$edad_don, pacientes_tab$grupo)),
  
  tibble(variable = "Sexo femenino receptor, n (%)",
         total  = n_pct(tot$SexF == 1, nt),
         unicos = n_pct(g0$SexF  == 1, n0),
         mixtos = n_pct(g1$SexF  == 1, n1),
         p = p_cat(pacientes_tab$SexF, pacientes_tab$grupo)),
  
  tibble(variable = "Sexo femenino donante, n (%)",
         total  = n_pct(tot$SexF_don == 1, nt),
         unicos = n_pct(g0$SexF_don  == 1, n0),
         mixtos = n_pct(g1$SexF_don  == 1, n1),
         p = p_cat(pacientes_tab$SexF_don, pacientes_tab$grupo)),
  
  # ── Seguimiento ──────────────────────────────────────────────────────────
  tibble(variable = "CARACTERÍSTICAS DE SEGUIMIENTO",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "Tiempo de seguimiento (años), mediana (RIC)",
         total  = med_ric(tot$t_seg_anios),
         unicos = med_ric(g0$t_seg_anios),
         mixtos = med_ric(g1$t_seg_anios),
         p = p_cont(pacientes_tab$t_seg_anios, pacientes_tab$grupo)),
  
  tibble(variable = "Episodios CMV totales, mediana (RIC)",
         total  = med_ric(tot$nepi_tot),
         unicos = med_ric(g0$nepi_tot),
         mixtos = med_ric(g1$nepi_tot),
         p = p_cont(pacientes_tab$nepi_tot, pacientes_tab$grupo)),
  
  tibble(variable = "Episodios CMV genotipados, mediana (RIC)",
         total  = med_ric(tot$nepi_gtc),
         unicos = med_ric(g0$nepi_gtc),
         mixtos = med_ric(g1$nepi_gtc),
         p = p_cont(pacientes_tab$nepi_gtc, pacientes_tab$grupo)),
  
  tibble(variable = "Determinaciones CV totales, mediana (RIC)",
         total  = med_ric(tot$Npcrs),
         unicos = med_ric(g0$Npcrs),
         mixtos = med_ric(g1$Npcrs),
         p = p_cont(pacientes_tab$Npcrs, pacientes_tab$grupo)),
  
  tibble(variable = "Antigenemias totales, mediana (RIC)",
         total  = med_ric(tot$Nag),
         unicos = med_ric(g0$Nag),
         mixtos = med_ric(g1$Nag),
         p = p_cont(pacientes_tab$Nag, pacientes_tab$grupo)),
  
  # ── Tipo de trasplante ────────────────────────────────────────────────────
  tibble(variable = "Tipo de trasplante, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$tipotx, pacientes_tab$grupo)),
  
  tibble(variable = "   Unipulmonar",
         total  = n_pct(tot$tipotx == "UNI", nt),
         unicos = n_pct(g0$tipotx  == "UNI", n0),
         mixtos = n_pct(g1$tipotx  == "UNI", n1), p = ""),
  
  tibble(variable = "   Bipulmonar",
         total  = n_pct(tot$tipotx == "BI", nt),
         unicos = n_pct(g0$tipotx  == "BI", n0),
         mixtos = n_pct(g1$tipotx  == "BI", n1), p = ""),
  
  # ── Enfermedad de base ────────────────────────────────────────────────────
  tibble(variable = "Enfermedad de base, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$ebase, pacientes_tab$grupo)),
  
  tibble(variable = "   EPID",
         total  = n_pct(tot$ebase == "EPID", nt),
         unicos = n_pct(g0$ebase  == "EPID", n0),
         mixtos = n_pct(g1$ebase  == "EPID", n1), p = ""),
  
  tibble(variable = "   EPOC",
         total  = n_pct(tot$ebase == "EPOC", nt),
         unicos = n_pct(g0$ebase  == "EPOC", n0),
         mixtos = n_pct(g1$ebase  == "EPOC", n1), p = ""),
  
  tibble(variable = "   FPI",
         total  = n_pct(tot$ebase == "FPI", nt),
         unicos = n_pct(g0$ebase  == "FPI", n0),
         mixtos = n_pct(g1$ebase  == "FPI", n1), p = ""),
  
  tibble(variable = "   FQ",
         total  = n_pct(tot$ebase == "FQ", nt),
         unicos = n_pct(g0$ebase  == "FQ", n0),
         mixtos = n_pct(g1$ebase  == "FQ", n1), p = ""),
  
  tibble(variable = "   HAP",
         total  = n_pct(tot$ebase == "HAP", nt),
         unicos = n_pct(g0$ebase  == "HAP", n0),
         mixtos = n_pct(g1$ebase  == "HAP", n1), p = ""),
  
  tibble(variable = "   Otros",
         total  = n_pct(!tot$ebase %in% c("EPID","EPOC","FPI","FQ","HAP"), nt),
         unicos = n_pct(!g0$ebase  %in% c("EPID","EPOC","FPI","FQ","HAP"), n0),
         mixtos = n_pct(!g1$ebase  %in% c("EPID","EPOC","FPI","FQ","HAP"), n1),
         p = ""),
  
  # ── Perfil serológico ─────────────────────────────────────────────────────
  tibble(variable = "Perfil serológico D/R, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$serDR, pacientes_tab$grupo)),
  
  tibble(variable = "   D-/R+",
         total  = n_pct(tot$serDR == "D-/R+", nt),
         unicos = n_pct(g0$serDR  == "D-/R+", n0),
         mixtos = n_pct(g1$serDR  == "D-/R+", n1), p = ""),
  
  tibble(variable = "   D+/R+",
         total  = n_pct(tot$serDR == "D+/R+", nt),
         unicos = n_pct(g0$serDR  == "D+/R+", n0),
         mixtos = n_pct(g1$serDR  == "D+/R+", n1), p = ""),
  
  tibble(variable = "   D+/R-",
         total  = n_pct(tot$serDR == "D+/R-", nt),
         unicos = n_pct(g0$serDR  == "D+/R-", n0),
         mixtos = n_pct(g1$serDR  == "D+/R-", n1), p = ""),
  
  # ── Inmunosupresión ───────────────────────────────────────────────────────
  tibble(variable = "INMUNOSUPRESIÓN BASAL, n (%)",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "   Tacrolimus",
         total  = n_pct(grepl("TACRO", tot$IS_base), nt),
         unicos = n_pct(grepl("TACRO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("TACRO", g1$IS_base),  n1),
         p = p_cat(grepl("TACRO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Ciclosporina",
         total  = n_pct(grepl("CICLO", tot$IS_base), nt),
         unicos = n_pct(grepl("CICLO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("CICLO", g1$IS_base),  n1),
         p = p_cat(grepl("CICLO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Micofenolato",
         total  = n_pct(grepl("MICO", tot$IS_base), nt),
         unicos = n_pct(grepl("MICO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("MICO", g1$IS_base),  n1),
         p = p_cat(grepl("MICO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Azatioprina",
         total  = n_pct(grepl("AZA", tot$IS_base), nt),
         unicos = n_pct(grepl("AZA", g0$IS_base),  n0),
         mixtos = n_pct(grepl("AZA", g1$IS_base),  n1),
         p = p_cat(grepl("AZA", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Esteroides",
         total  = "146 (100%)",
         unicos = "95 (100%)",
         mixtos = "51 (100%)",
         p = "\u2014")
)

# ── 4. Flextable ──────────────────────────────────────────────────────────────
filas_seccion <- which(filas$variable %in% c(
  "CARACTERÍSTICAS DEMOGRÁFICAS",
  "CARACTERÍSTICAS DE SEGUIMIENTO",
  "INMUNOSUPRESIÓN BASAL, n (%)"
))

filas_padre <- which(filas$variable %in% c(
  "Tipo de trasplante, n (%)",
  "Enfermedad de base, n (%)",
  "Perfil serológico D/R, n (%)"
))

colnames(filas) <- c("Variable",
                     "Total\n(n = 146)",
                     "Genotipos únicos\n(n = 95)",
                     "Genotipos mixtos\n(n = 51)",
                     "p")

# p significativas para resaltar
p_sig <- c("<0,001")

ft <- flextable(filas) %>%
  width(j = "Variable",                   width = 4.5) %>%
  width(j = "Total\n(n = 146)",           width = 2.0) %>%
  width(j = "Genotipos únicos\n(n = 95)", width = 2.0) %>%
  width(j = "Genotipos mixtos\n(n = 51)", width = 2.0) %>%
  width(j = "p",                          width = 0.9) %>%
  align(align = "left",   part = "all") %>%
  align(j = c(2,3,4,5), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = filas_seccion, part = "body") %>%
  bg(i = filas_seccion, bg = "#F0F0F0", part = "body") %>%
  italic(i = filas_padre, part = "body") %>%
  color(i = ~ p %in% p_sig, j = "p",
        color = "firebrick3", part = "body") %>%
  bold(i  = ~ p %in% p_sig, j = "p", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(i = filas_seccion - 1,
        border = fp_border(color = "grey70", width = 0.4),
        part = "body") %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "RIC: rango intercuartílico; EPID: enfermedad pulmonar intersticial difusa; EPOC: enfermedad pulmonar obstructiva crónica; FPI: fibrosis pulmonar idiopática; FQ: fibrosis quística; HAP: hipertensión arterial pulmonar.",
    "Variables continuas: test de Wilcoxon. Variables categóricas: test exacto de Fisher o chi-cuadrado.",
    "Valores de p en rojo indican diferencia estadísticamente significativa (p < 0,05)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R11. Características basales de los pacientes según presencia de genotipos mixtos gB/gH."
  ) %>%
  set_table_properties(layout = "autofit")

ft


## ---- [7b] Tabla 20. Enfermedad sintomatica por CMV y genotipos mixtos ----
# Apartado 5.7.2. Odds ratio cruda y ajustada.


library(dplyr)
library(lme4)

# ── 1. Preparar datos a nivel paciente ────────────────────────────────────────
pacientes_sintom <- pacientes %>%
  mutate(
    sintomatico = as.integer(vira == 0),
    gbghm_f     = factor(gbghm, levels = c(0,1), labels = c("Único","Mixto")),
    gbm_f       = factor(gbm,   levels = c(0,1), labels = c("Único","Mixto")),
    ghm_f       = factor(ghm,   levels = c(0,1), labels = c("Único","Mixto")),
    serDR_f     = factor(serDR)
  )

cat("── Distribución sintomático vs grupo ──\n")
pacientes_sintom %>%
  count(sintomatico, gbghm_f) %>%
  group_by(gbghm_f) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 2. Funciones ──────────────────────────────────────────────────────────────
calc_or <- function(var_mix, data, label) {
  m_crudo <- glm(sintomatico ~ get(var_mix),
                 data = data, family = binomial)
  or_c <- exp(coef(m_crudo)[2])
  ci_c <- exp(confint.default(m_crudo)[2,])
  p_c  <- summary(m_crudo)$coefficients[2,4]
  
  m_ajus <- glm(sintomatico ~ get(var_mix) + serDR_f,
                data = data, family = binomial)
  or_a <- exp(coef(m_ajus)[2])
  ci_a <- exp(confint.default(m_ajus)[2,])
  p_a  <- summary(m_ajus)$coefficients[2,4]
  
  cat(sprintf("\n── %s ──\n", label))
  cat(sprintf("OR cruda:    %.2f (%.2f–%.2f)  p = %.3f\n",
              or_c, ci_c[1], ci_c[2], p_c))
  cat(sprintf("OR ajustada: %.2f (%.2f–%.2f)  p = %.3f\n",
              or_a, ci_a[1], ci_a[2], p_a))
  
  tibble(
    genotipo    = label,
    or_cruda    = sprintf("%.2f (%.2f\u2013%.2f)", or_c, ci_c[1], ci_c[2]),
    p_cruda     = ifelse(p_c < 0.001, "<0,001",
                         sub("\\.", ",", sprintf("%.3f", p_c))),
    or_ajustada = sprintf("%.2f (%.2f\u2013%.2f)", or_a, ci_a[1], ci_a[2]),
    p_ajustada  = ifelse(p_a < 0.001, "<0,001",
                         sub("\\.", ",", sprintf("%.3f", p_a)))
  )
}

calc_or_mixto <- function(var_mix, data, label) {
  m <- tryCatch(
    glmer(sintomatico ~ get(var_mix) + (1|NHC),
          data = data, family = binomial,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(m)) {
    cat(sprintf("\n── %s: error en modelo mixto ──\n", label))
    return(tibble(genotipo = label, or_mixto = "error", p_mixto = "error"))
  }
  coefs <- summary(m)$coefficients
  or    <- exp(coefs[2, 1])
  se    <- coefs[2, 2]
  ci_lo <- exp(coefs[2, 1] - 1.96 * se)
  ci_hi <- exp(coefs[2, 1] + 1.96 * se)
  p     <- coefs[2, 4]
  
  cat(sprintf("\n── %s (modelo mixto) ──\n", label))
  cat(sprintf("OR: %.2f (%.2f–%.2f)  p = %.3f\n", or, ci_lo, ci_hi, p))
  
  tibble(
    genotipo = label,
    or_mixto = sprintf("%.2f (%.2f\u2013%.2f)", or, ci_lo, ci_hi),
    p_mixto  = ifelse(p < 0.001, "<0,001",
                      sub("\\.", ",", sprintf("%.3f", p)))
  )
}

freq_gt <- function(var, label, data) {
  tab    <- table(data[[var]], data$sintomatico)
  test   <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
  n0     <- sum(data$sintomatico == 0)
  n1     <- sum(data$sintomatico == 1)
  n_var0 <- sum(data[[var]] == 1 & data$sintomatico == 0, na.rm = TRUE)
  n_var1 <- sum(data[[var]] == 1 & data$sintomatico == 1, na.rm = TRUE)
  tibble(
    genotipo = label,
    asintom  = sprintf("%d (%.1f%%)", n_var0, 100*n_var0/n0),
    sintom   = sprintf("%d (%.1f%%)", n_var1, 100*n_var1/n1),
    p        = ifelse(test$p.value < 0.001, "<0,001",
                      sub("\\.", ",", sprintf("%.3f", test$p.value)))
  )
}

# ── 3. OR cruda y ajustada a nivel paciente ───────────────────────────────────
res_or <- bind_rows(
  calc_or("gbm_f",   pacientes_sintom, "gB mixtos"),
  calc_or("ghm_f",   pacientes_sintom, "gH mixtos"),
  calc_or("gbghm_f", pacientes_sintom, "gB+gH mixtos")
)

cat("\n── Tabla OR resumen ──\n")
print(res_or)

# ── 4. Preparar datos a nivel episodio ────────────────────────────────────────
dnaemias_sintom <- dnaemias %>%
  mutate(
    sintomatico = as.integer(!is.na(clinica) & clinica != "VA"),
    gbghm_f     = factor(gbghmix, levels = c(0,1), labels = c("Único","Mixto")),
    gbm_f       = factor(gbmix,   levels = c(0,1), labels = c("Único","Mixto")),
    ghm_f       = factor(ghmix,   levels = c(0,1), labels = c("Único","Mixto"))
  ) %>%
  left_join(
    pacientes %>% select(NHC, serDR_pac = serDR),
    by = "NHC"
  ) %>%
  mutate(serDR_f = factor(serDR_pac))

cat(sprintf("\nTotal episodios: %d\n", nrow(dnaemias_sintom)))
cat(sprintf("Sintomáticos: %d\n",     sum(dnaemias_sintom$sintomatico)))

# ── 5. Modelo mixto por episodio ──────────────────────────────────────────────
res_mixto <- bind_rows(
  calc_or_mixto("gbm_f",   dnaemias_sintom, "gB mixtos"),
  calc_or_mixto("ghm_f",   dnaemias_sintom, "gH mixtos"),
  calc_or_mixto("gbghm_f", dnaemias_sintom, "gB+gH mixtos")
)

cat("\n── Tabla OR mixto resumen ──\n")
print(res_mixto)

# ── 6. Frecuencias genotípicas por episodio ───────────────────────────────────
res_freq <- bind_rows(
  freq_gt("gbmix",   "gB mixtos",    dnaemias_sintom),
  freq_gt("ghmix",   "gH mixtos",    dnaemias_sintom),
  freq_gt("gbghmix", "gB+gH mixtos", dnaemias_sintom)
)

cat("\n── Frecuencias genotípicas ──\n")
print(res_freq)


library(flextable)
library(officer)

# ── Combinar resultados ───────────────────────────────────────────────────────
tabla_sintom <- res_or %>%
  left_join(res_mixto, by = "genotipo") %>%
  rename(
    "Genotipos mixtos"          = genotipo,
    "OR cruda\n(IC 95%)"        = or_cruda,
    "p"                         = p_cruda,
    "OR ajustada*\n(IC 95%)"    = or_ajustada,
    "p "                        = p_ajustada,
    "OR por episodio**\n(IC 95%)"= or_mixto,
    "p  "                       = p_mixto
  )

ft_sintom <- flextable(tabla_sintom) %>%
  width(j = "Genotipos mixtos",           width = 1.6) %>%
  width(j = "OR cruda\n(IC 95%)",         width = 2.0) %>%
  width(j = "p",                          width = 0.8) %>%
  width(j = "OR ajustada*\n(IC 95%)",     width = 2.0) %>%
  width(j = "p ",                         width = 0.8) %>%
  width(j = "OR por episodio**\n(IC 95%)",width = 2.0) %>%
  width(j = "p  ",                        width = 0.8) %>%
  align(align = "center", part = "all") %>%
  align(j = "Genotipos mixtos", align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bg(i = seq(2, nrow(tabla_sintom), 2), bg = "#F7F7F7", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  padding(padding.top = 3, padding.bottom = 3,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "*Ajustada por perfil serológico D/R.",
    "**Modelo de regresión logística mixta con efecto aleatorio por paciente."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R12. Asociación entre genotipos mixtos gB/gH y enfermedad sintomática por CMV."
  ) %>%
  set_table_properties(layout = "autofit")

ft_sintom


## ---- [7c] Figura 33. Tiempo desde el trasplante hasta el primer episodio ----
# Apartado 5.7.2. Doble panel: nivel paciente y nivel episodio.

## Panel A: Genotipos mixtos a lo largo del seguimiento (nivel paciente)
##          -> variables gbm/ghm/gbghm desde tabla `pacientes`
## Panel B: Genotipos mixtos en el primer episodio (nivel episodio)
##          -> variables gbmix/ghmix/gbghmix desde tabla `dnaemias`
## Filtro: primer_ev = 1
## Orden de paneles: D+/R- -> D+/R+ -> D-/R+

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

DODGE <- 0.7

#  Paletas (compartidas)
cols_fill <- c(
  "gB_\u00danico"     = "#AED6F1",  "gB_Mixto"     = "#4682B4",
  "gH_\u00danico"     = "pink",     "gH_Mixto"     = "hotpink3",
  "gB+gH_\u00danico"  = "#4CAF7D",  "gB+gH_Mixto"  = "#E07070"
)

cols_jit <- c(
  "gB_\u00danico"     = "#2980B9",  "gB_Mixto"     = "#0D2B45",
  "gH_\u00danico"     = "pink3",    "gH_Mixto"     = "palevioletred4",
  "gB+gH_\u00danico"  = "#2E7D52",  "gB+gH_Mixto"  = "#B04040"
)

df_prof <- tibble::tibble(
  serDR = factor(c("D+/R-", "D+/R+", "D-/R+"),
                 levels = c("D+/R-", "D+/R+", "D-/R+")),
  prof_meses = c(12, 6, 6)
)

#  FUNCION CONSTRUCTORA DEL PANEL
construir_panel <- function(df_long, etiquetas_riesgo, titulo_panel) {
  
  # Tests significativos del panel
  stat_test <- df_long %>%
    dplyr::group_by(serDR, glucoproteina) %>%
    rstatix::wilcox_test(meses_tx_inicio ~ genotipo) %>%
    rstatix::add_significance("p") %>%
    dplyr::ungroup() %>%
    dplyr::filter(p.signif != "ns") %>%
    dplyr::mutate(
      x_pos      = as.numeric(glucoproteina),
      xmin       = x_pos - DODGE / 4,
      xmax       = x_pos + DODGE / 4,
      y.position = 22
    )
  
  ggplot(df_long,
         aes(x = glucoproteina, y = meses_tx_inicio,
             fill = cat_color,
             group = interaction(glucoproteina, genotipo))) +
    geom_hline(data = df_prof,
               aes(yintercept = prof_meses),
               color = "grey55", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(
      width          = 0.55,
      position       = position_dodge(width = DODGE),
      color          = "grey25",
      linewidth      = 0.4,
      outlier.shape  = 21,
      outlier.size   = 1.6,
      outlier.fill   = "white",
      outlier.color  = "grey50",
      outlier.stroke = 0.3,
      outlier.alpha  = 0.7
    ) +
    geom_jitter(
      aes(color = cat_color),
      position    = position_jitterdodge(jitter.width = 0.10,
                                         dodge.width  = DODGE),
      size        = 1.0,
      alpha       = 0.4,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cols_fill, guide = "none") +
    scale_color_manual(values = cols_jit, guide = "none") +
    facet_wrap(~ serDR, labeller = labeller(serDR = etiquetas_riesgo)) +
    {
      if (nrow(stat_test) > 0)
        ggpubr::stat_pvalue_manual(
          stat_test,
          label        = "p.signif",
          xmin         = "xmin",
          xmax         = "xmax",
          y.position   = "y.position",
          tip.length   = 0.01,
          bracket.size = 0.5,
          size         = 6,
          color        = "grey20",
          inherit.aes  = FALSE
        )
    } +
    labs(
      x = NULL,
      y = "Meses desde TX a 1er episodio",
      title = titulo_panel
    ) +
    scale_y_continuous(
      breaks = seq(0, 24, 6),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    coord_cartesian(ylim = c(0, 24)) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      panel.spacing.x    = unit(15, "pt"),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 12, color = "grey20"),
      axis.text.y        = element_text(size = 11, color = "grey30"),
      axis.title.y       = element_text(size = 12, color = "grey15",
                                        margin = margin(r = 8)),
      plot.title         = element_text(size = 13, color = "grey15",
                                        face = "bold",
                                        margin = margin(b = 8)),
      strip.background   = element_rect(fill = "grey94",
                                        color = "grey70", linewidth = 0.4),
      strip.text         = element_text(size = 12, color = "grey15",
                                        face = "bold",
                                        margin = margin(t = 5, b = 5)),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(10, 12, 6, 10)
    )
}

#  PANEL A - Nivel PACIENTE (gbm, ghm, gbghm desde tabla pacientes)
df_fig_A <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44
  ) %>%
  dplyr::select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  dplyr::left_join(pacientes %>% dplyr::select(NHC, gbm, ghm, gbghm), by = "NHC") %>%
  dplyr::mutate(
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long_A <- df_fig_A %>%
  dplyr::select(NHC, serDR, meses_tx_inicio, gbm, ghm, gbghm) %>%
  tidyr::pivot_longer(
    cols      = c(gbm, ghm, gbghm),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "\u00danico"),
                           levels = c("\u00danico", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_\u00danico",   "gB_Mixto",
                                  "gH_\u00danico",   "gH_Mixto",
                                  "gB+gH_\u00danico","gB+gH_Mixto"))
  )

counts_A <- df_fig_A %>%
  dplyr::group_by(serDR) %>%
  dplyr::summarise(n = dplyr::n_distinct(NHC), .groups = "drop")

etiquetas_A <- setNames(
  sprintf("%s (n = %d)", counts_A$serDR, counts_A$n),
  as.character(counts_A$serDR)
)

panel_A <- construir_panel(
  df_long          = df_long_A,
  etiquetas_riesgo = etiquetas_A,
  titulo_panel     = "A. Genotipos únicos/mixtos a lo largo del seguimiento"
)

#  PANEL B - Nivel EPISODIO (gbmix, ghmix, gbghmix desde tabla dnaemias)
df_fig_B <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44,
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long_B <- df_fig_B %>%
  dplyr::select(NHC, serDR, meses_tx_inicio, gbmix, ghmix, gbghmix) %>%
  tidyr::pivot_longer(
    cols      = c(gbmix, ghmix, gbghmix),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbmix", "ghmix", "gbghmix"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "\u00danico"),
                           levels = c("\u00danico", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_\u00danico",   "gB_Mixto",
                                  "gH_\u00danico",   "gH_Mixto",
                                  "gB+gH_\u00danico","gB+gH_Mixto"))
  )

counts_B <- df_fig_B %>%
  dplyr::group_by(serDR) %>%
  dplyr::summarise(n = dplyr::n_distinct(NHC), .groups = "drop")

etiquetas_B <- setNames(
  sprintf("%s (n = %d)", counts_B$serDR, counts_B$n),
  as.character(counts_B$serDR)
)

panel_B <- construir_panel(
  df_long          = df_long_B,
  etiquetas_riesgo = etiquetas_B,
  titulo_panel     = "B. Genotipos únicos/mixtos en el primer episodio"
)

#  COMBINAR PANELES A y B
figura_R20_doble <- panel_A / panel_B +
  patchwork::plot_layout(heights = c(1, 1))

figura_R20_doble

ggsave("output/figura_R20_doble_paciente_vs_episodio.png",
       figura_R20_doble,
       width = 10, height = 9, dpi = 320, bg = "white")


## ---- [7d] Figura 34. Numero de episodios por paciente segun genotipos mixtos ----
# Apartado 5.7.2.

## ---- [X] *Figura R35bis. Distribución del nº de episodios por paciente

library(dplyr)
library(ggplot2)
library(scales)

# ── 1. Preparar datos ────────────────────────────────────────────────────────
df_dist <- pacientes %>%
  filter(!is.na(nepi_tot), !is.na(gbghm)) %>%
  mutate(
    genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                      levels = c("Único", "Mixto")),
    cat_ep = factor(case_when(
      nepi_tot == 1 ~ "1",
      nepi_tot == 2 ~ "2",
      nepi_tot == 3 ~ "3",
      nepi_tot == 4 ~ "4",
      nepi_tot >= 5 ~ "\u22655"
    ), levels = c("1", "2", "3", "4", "\u22655"))
  )

# Conteos por grupo y categoría
df_pir <- df_dist %>%
  count(genotipo, cat_ep) %>%
  group_by(genotipo) %>%
  mutate(
    pct       = 100 * n / sum(n),
    pct_signo = ifelse(genotipo == "Único", -pct, pct),
    etiqueta  = sprintf("%d (%.1f%%)", n, pct) %>%
      sub("\\.", ",", .)
  ) %>%
  ungroup()

# Mostrar conteos por consola para confirmar
cat("--- Distribución por número de episodios ---\n")
print(df_pir)

# ── 2. Paleta coherente con el resto del documento ──────────────────────────
cols_pir <- c("Único" = "#4CAF7D", "Mixto" = "#E07070")

n_unico <- sum(df_dist$genotipo == "Único")
n_mixto <- sum(df_dist$genotipo == "Mixto")

# ── 3. Figura pirámide ──────────────────────────────────────────────────────
fig_pir <- ggplot(df_pir, aes(x = pct_signo, y = cat_ep, fill = genotipo)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_text(
    aes(label = etiqueta,
        hjust = ifelse(genotipo == "Único", 1.05, -0.05)),
    size = 3.4, color = "grey15"
  ) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_x_continuous(
    name   = "Porcentaje de pacientes (%)",
    breaks = seq(-50, 50, 10),
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-55, 55),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = cols_pir,
                    labels = c(sprintf("Genotipos \u00fanicos (n = %d)", n_unico),
                               sprintf("Genotipos mixtos (n = %d)", n_mixto))) +
  labs(
    y    = "N episodios por paciente",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(12, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.y        = element_text(size = 12, color = "grey15", face = "bold"),
    axis.text.x        = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20",
                                      margin = margin(t = 6)),
    axis.title.y       = element_text(size = 11, color = "grey20",
                                      margin = margin(r = 8)),
    plot.margin        = margin(10, 14, 8, 8)
  )

fig_pir

ggsave("output/figura_R35bis_distribucion_episodios.png",
       fig_pir, width = 6, height = 3, dpi = 320, bg = "white")


library(dplyr)

df_test <- pacientes %>%
  filter(!is.na(nepi_tot), !is.na(gbghm)) %>%
  mutate(genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                           levels = c("Único", "Mixto")))

# Mann-Whitney
test_mw <- wilcox.test(nepi_tot ~ genotipo, data = df_test, exact = FALSE)
cat("Mann-Whitney p =", format.pval(test_mw$p.value, digits = 3), "\n")

# Chi-cuadrado de la distribución categórica (≥4 vs <4 episodios)
df_test <- df_test %>% mutate(cat_4 = ifelse(nepi_tot >= 4, "≥4", "<4"))
test_chi <- chisq.test(table(df_test$genotipo, df_test$cat_4))
cat("Chi-cuadrado p =", format.pval(test_chi$p.value, digits = 3), "\n")
print(table(df_test$genotipo, df_test$cat_4))

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

# ── 1. AUCCV acumulada por paciente ──────────────────────────────────────────
# Limpieza previa de columnas si existen de ejecuciones anteriores
pacientes <- pacientes %>%
  select(-any_of(c("AUCCV_AC", "n_episodios_AUC",
                   "AUCCV_AC.x", "AUCCV_AC.y",
                   "n_episodios_AUC.x", "n_episodios_AUC.y")))

auc_paciente <- dnaemias %>%
  group_by(NHC) %>%
  summarise(
    AUCCV_AC        = sum(AUCCV, na.rm = TRUE),
    n_episodios_AUC = sum(!is.na(AUCCV)),
    .groups         = "drop"
  )

pacientes <- pacientes %>%
  left_join(auc_paciente, by = "NHC")

# ── 2. Datos para comparaciones ──────────────────────────────────────────────
df_comp <- pacientes %>%
  filter(!is.na(n_episodios_AUC), n_episodios_AUC > 0) %>%
  mutate(
    gbm_lab   = factor(ifelse(gbm   == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto")),
    ghm_lab   = factor(ifelse(ghm   == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto")),
    gbghm_lab = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto"))
  )

# ── 3. Paletas por glucoproteína ─────────────────────────────────────────────
# gB: azul claro / azul oscuro
cols_fill_gb <- c("Único" = "#AED6F1", "Mixto" = "#1A5276")
cols_jit_gb  <- c("Único" = "#2980B9", "Mixto" = "#0D2B45")

# gH: rosa / rosa oscuro
cols_fill_gh <- c("Único" = "pink",    "Mixto" = "hotpink3")
cols_jit_gh  <- c("Único" = "pink3",   "Mixto" = "palevioletred4")

# gB+gH: verde / rojo
cols_fill_gbgh <- c("Único" = "#4CAF7D", "Mixto" = "#E07070")
cols_jit_gbgh  <- c("Único" = "#2E7D52", "Mixto" = "#B04040")

# ── 4. Función para construir cada panel ─────────────────────────────────────
# Límite Y fijo en 2000; asteriscos en vez de p
make_panel_paciente <- function(data, var_grupo, titulo,
                                cols_fill, cols_jit,
                                ylim_max = 1000) {
  
  d <- data %>%
    select(grupo = all_of(var_grupo), AUCCV_AC) %>%
    filter(!is.na(AUCCV_AC), !is.na(grupo))
  
  # Test sobre datos completos (sin filtrar por cap visual)
  stat <- d %>%
    rstatix::wilcox_test(AUCCV_AC ~ grupo) %>%
    rstatix::add_significance("p") %>%
    mutate(y.position = ylim_max * 0.95)
  
  # Datos con cap visual aplicado
  d_vis <- d %>% filter(AUCCV_AC <= ylim_max)
  
  ggplot(d_vis, aes(x = grupo, y = AUCCV_AC, fill = grupo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    geom_jitter(
      aes(color = grupo),
      width = 0.12, size = 1.1, alpha = 0.45,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_fill) +
    scale_color_manual(values = cols_jit) +
    stat_pvalue_manual(
      stat,
      label        = "p.signif",
      tip.length   = 0.01,
      bracket.size = 0.45,
      size         = 5,
      inherit.aes  = FALSE
    ) +
    labs(x = NULL,
         y = "AUCCV acumulada (UI/mL\u00b7d\u00eda)",
         title = titulo) +
    scale_y_continuous(
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 11, color = "grey25"),
      axis.text.y        = element_text(size = 10, color = "grey30"),
      axis.title.y       = element_text(size = 10, color = "grey20"),
      plot.title         = element_text(size = 13, color = "grey15",
                                        face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(8, 10, 6, 6)
    )
}

# ── 5. Construir los 3 paneles ───────────────────────────────────────────────
p_gb <- make_panel_paciente(
  df_comp, "gbm_lab", "gB",
  cols_fill = cols_fill_gb, cols_jit = cols_jit_gb
)

p_gh <- make_panel_paciente(
  df_comp, "ghm_lab", "gH",
  cols_fill = cols_fill_gh, cols_jit = cols_jit_gh
)

p_gbgh <- make_panel_paciente(
  df_comp, "gbghm_lab", "gB/gH",
  cols_fill = cols_fill_gbgh, cols_jit = cols_jit_gbgh
)

# ── 6. Combinar y mostrar ────────────────────────────────────────────────────
figura_AUCCV_AC <- p_gb | p_gh | p_gbgh

figura_AUCCV_AC

# ── 7. Guardar ───────────────────────────────────────────────────────────────
# ggsave("output/figura_R22_AUCCV_AC_mixto.png",
#        figura_AUCCV_AC,
#        width = 10, height = 4, dpi = 320, bg = "white")


## ---- [7e1] Tasa de episodios de DNAemia (IRR) ----
# Apartado 5.7.2. Regresion binomial negativa con offset de persona-anos.

library(dplyr)
library(MASS)
library(ggplot2)
library(survival)
library(survminer)

# Resolver conflictos
select <- dplyr::select
rm(offset)

# ── 1. Preparar datos IRR ─────────────────────────────────────────────────────
df_irr <- pacientes %>%
  dplyr::mutate(
    Ftx_dt       = as.Date(Ftx),
    Ffinseg_dt   = as.Date(Ffinseg),
    persona_anos = as.numeric(Ffinseg_dt - Ftx_dt) / 365.25,
    risk_f       = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo  = factor(gbghm, levels = c(0,1),
                          labels = c("Genotipos únicos",
                                     "Genotipos mixtos"))
  ) %>%
  dplyr::filter(persona_anos > 0)

# ── 2. Estadísticos descriptivos ──────────────────────────────────────────────
cat("── Total episodios CMV ──\n")
cat(sprintf("n = %d\n", sum(pacientes$nepi_tot, na.rm = TRUE)))

cat("\n── Tasas por grupo ──\n")
df_irr %>%
  dplyr::group_by(gbghm_grupo) %>%
  dplyr::summarise(
    n               = n(),
    total_episodios = sum(nepi_tot, na.rm = TRUE),
    total_pa        = round(sum(persona_anos, na.rm = TRUE), 1),
    tasa_pa         = round(sum(nepi_tot) / sum(persona_anos), 3),
    media_ep        = round(mean(nepi_tot, na.rm = TRUE), 2),
    mediana_ep      = median(nepi_tot, na.rm = TRUE),
    .groups         = "drop"
  ) %>% print()

# ── 3. Sobredispersión ────────────────────────────────────────────────────────
cat("\n── Sobredispersión ──\n")
df_irr %>%
  dplyr::group_by(gbghm_grupo) %>%
  dplyr::summarise(
    media    = round(mean(nepi_tot), 2),
    varianza = round(var(nepi_tot), 2),
    ratio    = round(var(nepi_tot) / mean(nepi_tot), 2),
    .groups  = "drop"
  ) %>% print()

m_pois <- glm(nepi_tot ~ gbghm_grupo + offset(log(persona_anos)),
              family = poisson, data = df_irr)
cat(sprintf("Poisson deviance/df = %.2f\n",
            m_pois$deviance / m_pois$df.residual))

# ── 4. Modelos binomial negativa ──────────────────────────────────────────────
cat("\n── IRR crudo (binomial negativa) ──\n")
m_nb <- glm.nb(nepi_tot ~ gbghm_grupo + offset(log(persona_anos)),
               data = df_irr)
print(round(exp(cbind(IRR = coef(m_nb), confint(m_nb))), 3))
cat("\np-valor gbghm:\n")
print(round(summary(m_nb)$coefficients, 3))

cat("\n── IRR ajustado (binomial negativa) ──\n")
m_nb_adj <- glm.nb(nepi_tot ~ gbghm_grupo + edad_tx + risk_f +
                     offset(log(persona_anos)), data = df_irr)
print(round(exp(cbind(IRR = coef(m_nb_adj), confint(m_nb_adj))), 3))
cat("\nCoeficientes con p-valores:\n")
print(round(summary(m_nb_adj)$coefficients, 3))


## ---- [7e2] Figura 35. Comparacion de la tasa de episodios por paciente ----
# Apartado 5.7.2.

library(dplyr)
library(tidyr)
library(MASS)
library(ggplot2)

select <- dplyr::select

# ── 1. Preparar datos base (a nivel paciente) ────────────────────────────────
df_irr3 <- pacientes %>%
  dplyr::mutate(
    Ftx_dt       = as.Date(Ftx),
    Ffinseg_dt   = as.Date(Ffinseg),
    persona_anos = as.numeric(Ffinseg_dt - Ftx_dt) / 365.25,
    risk_f       = factor(risk, levels = c("MEDIO", "ALTO"))
  ) %>%
  dplyr::filter(persona_anos > 0)

# ── 2. Funcion: IRR ajustado para una variable de grupo ──────────────────────
calc_irr <- function(data, var_grupo) {
  data$grupo <- data[[var_grupo]]
  m <- glm.nb(nepi_tot ~ grupo + edad_tx + risk_f + offset(log(persona_anos)),
              data = data)
  irr <- round(exp(coef(m)[2]), 2)
  ci  <- round(exp(confint.default(m)[2, ]), 2)   # Wald (rapido y estable)
  p   <- summary(m)$coefficients[2, 4]
  list(irr = irr, ci_lo = ci[1], ci_hi = ci[2], p = p)
}

# ── 3. Calcular IRR para las 3 glucoproteinas ────────────────────────────────
vars   <- c(gbm = "gbm", ghm = "ghm", gbghm = "gbghm")
labs_g <- c(gbm = "gB", ghm = "gH", gbghm = "gB/gH")

res_irr <- lapply(names(vars), function(v) {
  d <- df_irr3 %>%
    dplyr::mutate(grupo = factor(.data[[v]], levels = c(0, 1),
                                 labels = c("Únicos", "Mixtos"))) %>%
    dplyr::filter(!is.na(grupo))
  r <- calc_irr(d %>% dplyr::mutate(!!v := grupo), v)
  data.frame(
    glucoproteina = labs_g[[v]],
    irr = r$irr, ci_lo = r$ci_lo, ci_hi = r$ci_hi, p = r$p
  )
}) %>% bind_rows() %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina, levels = c("gB", "gH", "gB/gH")),
    label_irr = sprintf("IRR = %.2f (%.2f\u2013%.2f)\np = %s",
                        irr, ci_lo, ci_hi,
                        ifelse(p < 0.001, "<0,001",
                               sub("\\.", ",", sprintf("%.3f", p))))
  )

cat("── IRR ajustados ──\n"); print(res_irr)

# ── 4. Datos en formato largo: paciente x glucoproteina ──────────────────────
df_long3 <- df_irr3 %>%
  dplyr::select(NHC, nepi_tot, persona_anos, gbm, ghm, gbghm) %>%
  tidyr::pivot_longer(c(gbm, ghm, gbghm),
                      names_to = "gen", values_to = "mixto") %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(gen, levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB/gH")),
    grupo     = factor(ifelse(mixto == 1, "Mixtos", "Únicos"),
                       levels = c("Únicos", "Mixtos")),
    cat_color = paste(glucoproteina,
                      ifelse(mixto == 1, "Mixto", "Unico"), sep = "_")
  )

# ── 5. Tasas (ep/paciente-año) por panel y grupo ─────────────────────────────
tasas3 <- df_long3 %>%
  dplyr::group_by(glucoproteina, grupo) %>%
  dplyr::summarise(
    tasa = round(sum(nepi_tot) / sum(persona_anos), 2),
    .groups = "drop"
  )

# ── 6. n por panel y grupo (para etiquetas del eje) ──────────────────────────
n_grupo <- df_long3 %>%
  dplyr::count(glucoproteina, grupo)

# ── 7. Paletas (gB azul | gH rosa | gB/gH verde-rojo) ────────────────────────
cols_fill <- c(
  "gB_Unico"    = "#AED6F1", "gB_Mixto"    = "#4682B4",
  "gH_Unico"    = "pink",    "gH_Mixto"    = "hotpink3",
  "gB/gH_Unico" = "darkseagreen3", "gB/gH_Mixto" = "#c0392b"
)
cols_jit <- c(
  "gB_Unico"    = "#2980B9",        "gB_Mixto"    = "#0D2B45",
  "gH_Unico"    = "palevioletred3", "gH_Mixto"    = "palevioletred4",
  "gB/gH_Unico" = "#698B69",        "gB/gH_Mixto" = "#922b21"
)

# ── 8. Figura ────────────────────────────────────────────────────────────────
figura_R35_triple <- ggplot(df_long3,
                            aes(x = grupo, y = nepi_tot)) +
  geom_boxplot(aes(fill = cat_color),
               width = 0.45, outlier.shape = NA, alpha = 0.8,
               color = "black", linewidth = 0.6) +
  geom_jitter(aes(color = cat_color), width = 0.15, size = 1.8,
              alpha = 0.6, shape = 16, show.legend = FALSE) +
  # Etiqueta de tasa sobre cada caja
  geom_label(data = tasas3,
             aes(x = grupo, y = 10,
                 label = paste0(gsub("\\.", ",", tasa))),
             fill = "white", color = "black", size = 4,
             label.size = 0.25, fontface = "bold") +
  # Barra de significacion + IRR por panel
  
  geom_text(data = res_irr,
            aes(x = 1.5, y = 13, label = label_irr),
            size = 4, color = "black", fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~ glucoproteina) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  scale_y_continuous(
    name   = "Número de episodios de DNAemia",
    breaks = seq(0, 12, 2),
    limits = c(0, 15),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.spacing.x    = unit(14, "pt"),
    axis.title.y       = element_text(size = 13, color = "black",
                                      margin = margin(r = 10)),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 14, color = "black",
                                      face = "bold",
                                      margin = margin(t = 2, b = 6)),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 16, 10, 12)
  )

print(figura_R35_triple)

# ════════════════════════════════════════════════════════════════════════════
#  DATOS NUMERICOS — Figura R35 triple (episodios CMV por gB, gH, gB/gH)
# ════════════════════════════════════════════════════════════════════════════

# ── 1. n, tasa, media, mediana (RIC) por panel y grupo ───────────────────────
cat("\n========== DESCRIPTIVO POR GRUPO ==========\n")
tabla_desc <- df_long3 %>%
  dplyr::group_by(glucoproteina, grupo) %>%
  dplyr::summarise(
    n          = dplyr::n(),
    episodios  = sum(nepi_tot),
    persona_a  = round(sum(persona_anos), 1),
    tasa_pa    = round(sum(nepi_tot) / sum(persona_anos), 2),
    media      = round(mean(nepi_tot), 2),
    mediana    = median(nepi_tot),
    q1         = quantile(nepi_tot, 0.25),
    q3         = quantile(nepi_tot, 0.75),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    `mediana (RIC)` = sprintf("%g (%g\u2013%g)", mediana, q1, q3)
  )
print(tabla_desc, n = Inf)

# ── 2. IRR ajustados con IC y p (formato tesis) ──────────────────────────────
cat("\n========== IRR AJUSTADOS (NB: edad + risk) ==========\n")
res_irr %>%
  dplyr::mutate(
    `IRR (IC95%)` = sprintf("%.2f (%.2f\u2013%.2f)", irr, ci_lo, ci_hi),
    `p`           = ifelse(p < 0.001, "<0,001",
                           sub("\\.", ",", sprintf("%.3f", p)))
  ) %>%
  dplyr::select(glucoproteina, `IRR (IC95%)`, `p`) %>%
  print()

# ── 3. p crudas de comparacion de distribucion (Wilcoxon) por panel ──────────
#     (util si quieres contrastar tasa vs distribucion de episodios)
cat("\n========== WILCOXON nepi_tot ~ grupo (crudo) ==========\n")
df_long3 %>%
  dplyr::group_by(glucoproteina) %>%
  dplyr::summarise(
    p_raw = wilcox.test(nepi_tot ~ grupo, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_fmt = ifelse(p_raw < 0.001, "<0,001",
                   sub("\\.", ",", sprintf("%.3f", p_raw)))
  ) %>%
  print()


## ---- [7e3] Figura 36. Exposicion viral acumulada (AUC) por paciente ----
# Apartado 5.7.2.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ── 1. Episodios: todos los de dnaemias, excluyendo los sin aclaramiento ──────
n_inicial <- nrow(dnaemias)

df_ep <- dnaemias %>%
  filter(!is.na(aclaramiento_dias), !is.na(AUCCV))

n_excluidos <- n_inicial - nrow(df_ep)

cat(sprintf("Episodios totales en dnaemias: %d\n", n_inicial))
cat(sprintf("Excluidos (sin aclaramiento/AUC): %d\n", n_excluidos))
cat(sprintf("Episodios analizados: %d\n", nrow(df_ep)))

# ── 2. AUC acumulada por paciente (suma de AUCCV de sus episodios) ────────────
auc_pac <- df_ep %>%
  group_by(NHC) %>%
  summarise(AUCa = sum(AUCCV, na.rm = TRUE), .groups = "drop") %>%
  left_join(pacientes %>% select(NHC, gbm, ghm, gbghm), by = "NHC")

cat(sprintf("Pacientes con AUC acumulada: %d\n", nrow(auc_pac)))

# ── 3. Formato largo: una fila por paciente x glucoproteina ──────────────────
df_long <- auc_pac %>%
  pivot_longer(c(gbm, ghm, gbghm),
               names_to = "gen", values_to = "mixto") %>%
  filter(!is.na(mixto)) %>%
  mutate(
    glucoproteina = factor(gen,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB/gH")),
    genotipo  = factor(ifelse(mixto == 1, "Mixtos", "\u00danicos"),
                       levels = c("\u00danicos", "Mixtos")),
    # clave de color independiente del texto del eje (usa 0/1, no la etiqueta)
    cat_color = paste(glucoproteina, ifelse(mixto == 1, "Mixto", "\u00danico"),
                      sep = "_")
  )

# ── 4. Medianas, RIC y test de Wilcoxon por glucoproteina ────────────────────
resumen <- df_long %>%
  group_by(glucoproteina, genotipo) %>%
  summarise(
    n       = n(),
    mediana = round(median(AUCa)),
    q1      = round(quantile(AUCa, 0.25)),
    q3      = round(quantile(AUCa, 0.75)),
    .groups = "drop"
  ) %>%
  mutate(resumen = sprintf("%d (RIC %d\u2013%d)", mediana, q1, q3))

cat("\n── Medianas (RIC) por grupo ──\n"); print(resumen)

stat_test <- df_long %>%
  group_by(glucoproteina) %>%
  rstatix::wilcox_test(AUCa ~ genotipo) %>%
  rstatix::add_significance("p")

cat("\n── Test de Wilcoxon ──\n"); print(stat_test)

# ── 5. Paletas (gB azul | gH rosa | gB/gH verde-rojo) ────────────────────────
cols_fill <- c(
  "gB_\u00danico"    = "#AED6F1", "gB_Mixto"    = "#4682B4",
  "gH_\u00danico"    = "pink",    "gH_Mixto"    = "hotpink3",
  "gB/gH_\u00danico" = "#4CAF7D", "gB/gH_Mixto" = "#E07070"
)
cols_jit <- c(
  "gB_\u00danico"    = "#2980B9",        "gB_Mixto"    = "#0D2B45",
  "gH_\u00danico"    = "palevioletred3", "gH_Mixto"    = "palevioletred4",
  "gB/gH_\u00danico" = "#2E7D52",        "gB/gH_Mixto" = "#B04040"
)

# ── 6. Posicion de los corchetes de significacion (dentro del rango 0-1200) ──
stat_test <- stat_test %>% mutate(y.position = 1050)

# ── 7. Figura ────────────────────────────────────────────────────────────────
fig_auc <- ggplot(df_long, aes(x = genotipo, y = AUCa)) +
  geom_boxplot(aes(fill = cat_color),
               width = 0.45, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color = cat_color),
              width = 0.15, size = 1.8, alpha = 0.5,
              shape = 16, show.legend = FALSE) +
  facet_wrap(~ glucoproteina) +
  stat_pvalue_manual(
    stat_test,
    label        = "p.signif",
    tip.length   = 0.01,
    bracket.size = 0.6,
    size         = 5,
    color        = "black"
  ) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  scale_y_continuous(
    name   = "AUC acumulada (UI/mL\u00b7d\u00eda)",
    breaks = seq(0, 1200, 250),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 1200)) +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "none",
    panel.spacing.x    = unit(14, "pt"),
    axis.title.y       = element_text(size = 13, color = "black",
                                      margin = margin(r = 10)),
    axis.text.x        = element_text(size = 13, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 14, color = "black",
                                      face = "bold",
                                      margin = margin(t = 2, b = 6)),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 16, 10, 12)
  )

print(fig_auc)


## ---- [7f1] Recurrencia de DNAemia: construccion de la base y modelos de Cox ----
# Apartado 5.7.2. Tres ventanas temporales.

# ── 1. Identificar primer episodio por paciente ─────────────────────────────
# Definición: el episodio con Fini más cercano a Ftx (mínimo Fini por NHC)

primer_episodio <- dnaemias %>%
  mutate(
    Fini = as.Date(Fini),
    Ffin = as.Date(Ffin)
  ) %>%
  group_by(NHC) %>%
  arrange(Fini, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Ftx              = as.Date(Ftx),
    dias_tx_a_inicio = as.numeric(Fini - Ftx)
  )

# ── 2. Excluir pacientes sin Ffin documentado ────────────────────────────────
# (4 fallecidos antes de aclaramiento + 3 trasladados = 7 excluidos)
primer_ep_evaluables <- primer_episodio %>%
  filter(!is.na(Ffin))

cat(sprintf("── Pacientes evaluables: %d ──\n", nrow(primer_ep_evaluables)))

# ── 3. Identificar segundo episodio por paciente ────────────────────────────
segundo_episodio <- dnaemias %>%
  mutate(
    Fini = as.Date(Fini),
    Ffin = as.Date(Ffin)
  ) %>%
  inner_join(
    primer_ep_evaluables %>% select(NHC, Ffin_primer = Ffin),
    by = "NHC"
  ) %>%
  filter(Fini > Ffin_primer) %>%
  group_by(NHC) %>%
  arrange(Fini, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(NHC,
         Fini_segundo    = Fini,
         Ffin_segundo    = Ffin,
         gbghmix_segundo = gbghmix,
         gbmix_segundo   = gbmix,
         ghmix_segundo   = ghmix,
         cvini_segundo   = cvini,
         clinica_segundo = clinica)

# ── 4. Construir base df_recur ──────────────────────────────────────────────
df_recur <- primer_ep_evaluables %>%
  select(NHC,
         Ftx,
         Fini_primer    = Fini,
         Ffin_primer    = Ffin,
         gbghmix_primer = gbghmix,
         gbmix_primer   = gbmix,
         ghmix_primer   = ghmix,
         cvini_primer   = cvini,
         cvpico_primer  = cvpico,
         dias_tx_a_inicio) %>%
  left_join(segundo_episodio, by = "NHC") %>%
  left_join(
    pacientes %>%
      mutate(
        Ffinseg = as.Date(Ffinseg),
        Fexitus = as.Date(Fexitus)
      ) %>%
      select(NHC, Ffinseg, exitus, Fexitus, serDR),
    by = "NHC"
  )

# ── 5. Calcular tiempos y eventos para análisis ──────────────────────────────
df_recur <- df_recur %>%
  mutate(
    fecha_fin_seg = pmin(Ffinseg, Fexitus, na.rm = TRUE),
    dias_a_recur  = as.numeric(Fini_segundo - Ffin_primer),
    dias_seguim   = as.numeric(fecha_fin_seg - Ffin_primer),
    
    # Recurrencia a 6 meses
    recur_180 = case_when(
      !is.na(dias_a_recur) & dias_a_recur <= 180 ~ 1L,
      TRUE                                       ~ 0L
    ),
    tiempo_km = case_when(
      recur_180 == 1     ~ dias_a_recur,
      dias_seguim < 180  ~ dias_seguim,
      TRUE               ~ 180
    ),
    evento_km = recur_180,
    
    # Recurrencia a 12 meses
    recur_365 = case_when(
      !is.na(dias_a_recur) & dias_a_recur <= 365 ~ 1L,
      TRUE                                       ~ 0L
    ),
    tiempo_km_365 = case_when(
      recur_365 == 1     ~ dias_a_recur,
      dias_seguim < 365  ~ dias_seguim,
      TRUE               ~ 365
    ),
    evento_km_365 = recur_365,
    
    # Recurrencia en seguimiento completo
    recur_global = as.integer(!is.na(Fini_segundo)),
    tiempo_completo = case_when(
      recur_global == 1 ~ dias_a_recur,
      TRUE              ~ dias_seguim
    ),
    
    # Factores de exposición
    exp_grupo = factor(gbghmix_primer, levels = c(0, 1),
                       labels = c("Único", "Mixto")),
    gb_grupo  = factor(gbmix_primer,   levels = c(0, 1),
                       labels = c("gB único", "gB mixto")),
    gh_grupo  = factor(ghmix_primer,   levels = c(0, 1),
                       labels = c("gH único", "gH mixto")),
    serDR_f   = factor(serDR)
  )


# Objetos Surv
surv_180      <- Surv(time = df_recur$tiempo_km,     event = df_recur$evento_km)
surv_365      <- Surv(time = df_recur$tiempo_km_365, event = df_recur$evento_km_365)
surv_completo <- Surv(time = df_recur$tiempo_completo, event = df_recur$recur_global)

# ── 6 meses ─────────────────────────────────────────────────────────────────
cox_gbgh_6m     <- coxph(surv_180 ~ exp_grupo, data = df_recur)
cox_gb_6m       <- coxph(surv_180 ~ gb_grupo,  data = df_recur)
cox_gh_6m       <- coxph(surv_180 ~ gh_grupo,  data = df_recur)
cox_gbgh_6m_adj <- coxph(surv_180 ~ exp_grupo + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gb_6m_adj   <- coxph(surv_180 ~ gb_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gh_6m_adj   <- coxph(surv_180 ~ gh_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)

# ── 12 meses ────────────────────────────────────────────────────────────────
cox_gbgh_12m     <- coxph(surv_365 ~ exp_grupo, data = df_recur)
cox_gb_12m       <- coxph(surv_365 ~ gb_grupo,  data = df_recur)
cox_gh_12m       <- coxph(surv_365 ~ gh_grupo,  data = df_recur)
cox_gbgh_12m_adj <- coxph(surv_365 ~ exp_grupo + serDR_f + cvpico_primer,
                          data = df_recur)
cox_gb_12m_adj   <- coxph(surv_365 ~ gb_grupo  + serDR_f + cvpico_primer,
                          data = df_recur)
cox_gh_12m_adj   <- coxph(surv_365 ~ gh_grupo  + serDR_f + cvpico_primer,
                          data = df_recur)

# ── Seguimiento completo ────────────────────────────────────────────────────
cox_gbgh_sc     <- coxph(surv_completo ~ exp_grupo, data = df_recur)
cox_gb_sc       <- coxph(surv_completo ~ gb_grupo,  data = df_recur)
cox_gh_sc       <- coxph(surv_completo ~ gh_grupo,  data = df_recur)
cox_gbgh_sc_adj <- coxph(surv_completo ~ exp_grupo + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gb_sc_adj   <- coxph(surv_completo ~ gb_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gh_sc_adj   <- coxph(surv_completo ~ gh_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)


## ---- [7f2] Tabla 21. Recurrencia de episodios segun el patron genotipico ----
# Apartado 5.7.2.

# ── Funciones auxiliares ────────────────────────────────────────────────────
extraer_hr_fmt <- function(cox) {
  hr <- exp(coef(cox)[1])
  ci <- exp(confint(cox)[1, ])
  p  <- summary(cox)$coefficients[1, 5]
  list(
    hr_ci = sprintf("%.2f (%.2f\u2013%.2f)", hr, ci[1], ci[2]),
    p     = ifelse(p < 0.001, "<0,001",
                   sub("\\.", ",", sprintf("%.3f", p)))
  )
}

fmt_n_pct <- function(n, N) {
  sub("\\.", ",", sprintf("%d (%.1f%%)", n, 100 * n / N))
}

# ── Construir filas ─────────────────────────────────────────────────────────
filas_R13 <- bind_rows(
  
  # Sección 1: 6 meses
  tibble(variable = "RECURRENCIA A 6 MESES (n = 33)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(21, 99),
         mixtos   = fmt_n_pct(12, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_6m_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(25, 108),
         mixtos   = fmt_n_pct(8, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_6m_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(24, 112),
         mixtos   = fmt_n_pct(9, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_6m_adj)$p),
  
  # Sección 2: 12 meses
  tibble(variable = "RECURRENCIA A 12 MESES (n = 42)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(25, 99),
         mixtos   = fmt_n_pct(17, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_12m_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(30, 108),
         mixtos   = fmt_n_pct(12, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_12m_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(29, 112),
         mixtos   = fmt_n_pct(13, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_12m_adj)$p),
  
  # Sección 3: Seguimiento completo
  tibble(variable = "RECURRENCIA EN SEGUIMIENTO COMPLETO (n = 56)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(35, 99),
         mixtos   = fmt_n_pct(21, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_sc_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(42, 108),
         mixtos   = fmt_n_pct(14, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_sc_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(40, 112),
         mixtos   = fmt_n_pct(16, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_sc_adj)$p)
)

colnames(filas_R13) <- c("Genotipos",
                         "Únicos\nn (%)",
                         "Mixtos\nn (%)",
                         "HR cruda\n(IC 95%)",
                         "p",
                         "HR ajustada*\n(IC 95%)",
                         "p ")

# ── Filas especiales por índice ─────────────────────────────────────────────
filas_seccion_R13 <- c(1, 5, 9)

es_p_sig <- function(x) {
  num <- suppressWarnings(as.numeric(sub(",", ".", x)))
  !is.na(num) & num < 0.05
}

filas_p_sig_crudo <- which(es_p_sig(filas_R13[[5]]))
filas_p_sig_adj   <- which(es_p_sig(filas_R13[[7]]))

# ── Construir flextable ─────────────────────────────────────────────────────
ft_R13 <- flextable(filas_R13) %>%
  width(j = 1, width = 2.0) %>%
  width(j = 2, width = 1.4) %>%
  width(j = 3, width = 1.4) %>%
  width(j = 4, width = 1.8) %>%
  width(j = 5, width = 0.9) %>%
  width(j = 6, width = 1.8) %>%
  width(j = 7, width = 0.9) %>%
  align(align = "left", part = "all") %>%
  align(j = c(2, 3, 4, 5, 6, 7), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = filas_seccion_R13, part = "body") %>%
  bg(i = filas_seccion_R13, bg = "#F0F0F0", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(i = c(4, 8),
        border = fp_border(color = "grey70", width = 0.4),
        part = "body") %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "*Modelo de regresión de Cox ajustado por perfil serológico D/R y carga viral pico del primer episodio.",
    "Las HR comparan el grupo mixto frente al grupo único, considerado como referencia.",
    "Valores de p en rojo indican diferencia estadísticamente significativa (p < 0,05)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R13. Recurrencia de DNAemia por CMV según el patrón genotípico del primer episodio en tres ventanas temporales."
  )

if (length(filas_p_sig_crudo) > 0) {
  ft_R13 <- ft_R13 %>%
    color(i = filas_p_sig_crudo, j = 5, color = "firebrick3", part = "body") %>%
    bold( i = filas_p_sig_crudo, j = 5, part = "body")
}
if (length(filas_p_sig_adj) > 0) {
  ft_R13 <- ft_R13 %>%
    color(i = filas_p_sig_adj, j = 7, color = "firebrick3", part = "body") %>%
    bold( i = filas_p_sig_adj, j = 7, part = "body")
}

ft_R13 <- ft_R13 %>%
  set_table_properties(layout = "autofit")

ft_R13


## ---- [7f3] Figura 37. Probabilidad libre de recurrencia a 12 meses ----
# Apartado 5.7.2. Kaplan-Meier para gB, gH y gB/gH.

library(survival)
library(survminer)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── 1. Ajustes KM a 12 meses (objeto surv_365 y grupos ya definidos) ─────────
km_gb   <- survfit(surv_365 ~ gb_grupo,  data = df_recur)
km_gh   <- survfit(surv_365 ~ gh_grupo,  data = df_recur)
km_gbgh <- survfit(surv_365 ~ exp_grupo, data = df_recur)

# ── 2. p de log-rank ─────────────────────────────────────────────────────────
p_lr <- function(formula) {
  ld <- survdiff(formula, data = df_recur)
  1 - pchisq(ld$chisq, df = length(ld$n) - 1)
}
p_gb   <- p_lr(surv_365 ~ gb_grupo)
p_gh   <- p_lr(surv_365 ~ gh_grupo)
p_gbgh <- p_lr(surv_365 ~ exp_grupo)

fmt_p <- function(p) sub("\\.", ",", formatC(p, format = "f", digits = 3))

# ── Paleta: VERDE = genotipo unico | ROJO = genotipo mixto ───────────────────
pal_km <- c("#2E8B57", "firebrick3")   # verde mar / rojo

# ── 3. Eje X en meses (marcas a 0, 3, 6, 9, 12) ──────────────────────────────
meses_breaks <- c(0, 3, 6, 9, 12)
dias_breaks  <- meses_breaks * 30.44

# ── 4. Funcion para construir un panel limpio ────────────────────────────────
panel_12m <- function(km_fit, p_val, titulo, mostrar_y) {
  gs <- ggsurvplot(
    km_fit,
    pval           = FALSE,
    conf.int       = TRUE,
    conf.int.alpha = 0.08,                 # ← menos opacidad en el sombreado
    conf.int.style = "ribbon",
    risk.table     = TRUE,
    risk.table.height = 0.26,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 3.5,
    palette        = pal_km,
    legend.title   = "",
    legend.labs    = c("Genotipos \u00fanicos", "Genotipos mixtos"),
    xlab           = "Meses desde aclaramiento",
    ylab           = if (mostrar_y) "Libre de recurrencia" else "",
    title          = titulo,
    risk.table.title = "Pacientes en riesgo",
    ggtheme        = theme_classic(base_size = 12),
    font.legend    = c(11, "plain", "black"),
    font.x         = c(11, "plain", "black"),
    font.y         = c(11, "plain", "black"),
    font.tickslab  = c(10, "plain", "black"),
    xlim           = c(0, 365),
    ylim           = c(0, 1)
  )
  
  # Curva: sin grid, p anotada, eje en meses
  gs$plot <- gs$plot +
    scale_x_continuous(breaks = dias_breaks, labels = meses_breaks,
                       limits = c(0, 365),
                       expand = expansion(mult = c(0.05, 0.02))) +
    annotate("text", x = 20, y = 0.08,
             label = paste0("italic('p')*' = '*'", fmt_p(p_val), "'"),
             parse = TRUE, hjust = 0, size = 3.6, color = "black") +
    labs(x = "Meses desde aclaramiento") +
    theme(
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      legend.position  = "none"
    )
  
  # Tabla de riesgo: sin grid, eje en meses
  gs$table <- gs$table +
    scale_x_continuous(breaks = dias_breaks, labels = meses_breaks,
                       limits = c(0, 365),
                       expand = expansion(mult = c(0.09, 0))) +
    labs(x = NULL, title = "Pacientes en riesgo") +
    theme(
      plot.title       = element_text(size = 10, face = "plain", color = "black"),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank()
    )
  
  gs
}

# ── 5. Construir los tres paneles ────────────────────────────────────────────
g_gb   <- panel_12m(km_gb,   p_gb,   "gB",    mostrar_y = TRUE)
g_gh   <- panel_12m(km_gh,   p_gh,   "gH",    mostrar_y = TRUE)
g_gbgh <- panel_12m(km_gbgh, p_gbgh, "gB/gH", mostrar_y = TRUE)

# ── 6. Combinar curva + su tabla de riesgo en un bloque vertical ─────────────
bloque <- function(gs) {
  (gs$plot / gs$table) + plot_layout(heights = c(4, 1.1))
}
b_gb   <- bloque(g_gb)
b_gh   <- bloque(g_gh)
b_gbgh <- bloque(g_gbgh)

# ── 7. Leyenda compartida en DOS FILAS (una por grupo) ───────────────────────
leyenda_plot <- ggsurvplot(
  km_gb, palette = pal_km,
  legend.title = "", legend.labs = c("Genotipos \u00fanicos", "Genotipos mixtos"),
  ggtheme = theme_classic(base_size = 13),
  font.legend = c(11, "plain", "black")
)$plot +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.position   = "bottom",
        legend.key.height = unit(14, "pt"),
        legend.key.width  = unit(20, "pt"),
        legend.spacing.y  = unit(2, "pt"),
        legend.text       = element_text(size = 11, margin = margin(l = 4)),
        legend.margin     = margin(t = 0, r = 0, b = 60, l = 0))  # ← b sube la leyenda

leyenda <- cowplot::get_legend(leyenda_plot)
celda_leyenda <- patchwork::wrap_elements(full = leyenda)


leyenda <- cowplot::get_legend(leyenda_plot)
celda_leyenda <- patchwork::wrap_elements(full = leyenda)

# ── 8. Composicion 2x2: gB | gH  /  gB/gH | leyenda ──────────────────────────
figura_km_12m <-
  (b_gb   | b_gh) /
  (b_gbgh | celda_leyenda) &
  theme(plot.background = element_rect(fill = "white", color = NA))

print(figura_km_12m)

ggsave("output/Figura_KM_recurrencia_12m_2x2.png",
       figura_km_12m,
       width = 24, height = 20, units = "cm", dpi = 300, bg = "white")


## ---- [7g1] Infecciones oportunistas no relacionadas con CMV: analisis ----
# Apartado 5.7.3.

library(survival)
library(survminer)

# ── 6. Preparar datos ─────────────────────────────────────────────────────────
df_opo <- pacientes %>%
  mutate(
    Ftx_dt           = as.Date(Ftx),
    F_opo_dt         = as.Date(F_opo),
    Ffinseg_dt       = as.Date(Ffinseg),
    fecha_fin_real   = if_else(inf_opo == 1, F_opo_dt, Ffinseg_dt),
    t_analisis_anios = as.numeric(fecha_fin_real - Ftx_dt) / 365.25,
    evento_opo       = if_else(inf_opo == 1, 1L, 0L),
    risk_f           = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo      = factor(gbghm, levels = c(0,1),
                              labels = c("Genotipos únicos",
                                         "Genotipos mixtos"))
  ) %>%
  filter(t_analisis_anios > 0)

# ── 7. Descriptivo ────────────────────────────────────────────────────────────
cat("\n── Infecciones oportunistas por grupo ──\n")
df_opo %>%
  group_by(gbghm_grupo) %>%
  summarise(
    n     = n(),
    n_opo = sum(evento_opo),
    pct   = round(n_opo / n * 100, 1),
    .groups = "drop"
  ) %>% print()

cat("\n── Total pacientes con IO ──\n")
cat(sprintf("n = %d (%.1f%%)\n",
            sum(pacientes$inf_opo == 1, na.rm = TRUE),
            100 * mean(pacientes$inf_opo == 1, na.rm = TRUE)))

cat("\n── Número de episodios por paciente ──\n")
table(pacientes$n_opo[pacientes$inf_opo == 1], useNA = "always") %>% print()

cat("\n── Etiología IO ──\n")
table(pacientes$det_opo, useNA = "always") %>%
  sort(decreasing = TRUE) %>% print()

cat("\n── Test chi-cuadrado IO vs grupo ──\n")
print(chisq.test(table(df_opo$gbghm_grupo, df_opo$evento_opo)))

# ── 8. KM + log-rank ─────────────────────────────────────────────────────────
surv_opo <- Surv(df_opo$t_analisis_anios, df_opo$evento_opo)
km_opo   <- survfit(surv_opo ~ gbghm_grupo, data = df_opo)
lr_opo   <- survdiff(surv_opo ~ gbghm_grupo, data = df_opo)
p_lr     <- 1 - pchisq(lr_opo$chisq, df = 1)
cat(sprintf("\n── Log-rank IO: p = %.3f ──\n", p_lr))

cat("\n── Supervivencia libre IO a 5 años ──\n")
print(summary(km_opo, times = 5))

# ── 9. Cox univariante y ajustado ─────────────────────────────────────────────
cox_opo_crudo <- coxph(surv_opo ~ gbghm_grupo, data = df_opo)
cox_opo_ajust <- coxph(surv_opo ~ gbghm_grupo + edad_tx + risk_f,
                       data = df_opo)

cat("\n── Cox IO crudo ──\n")
hr_oc <- exp(coef(cox_opo_crudo)[1])
ci_oc <- exp(confint(cox_opo_crudo)[1,])
p_oc  <- summary(cox_opo_crudo)$coefficients[1,5]
cat(sprintf("HR: %.2f (%.2f–%.2f)  p = %.3f\n",
            hr_oc, ci_oc[1], ci_oc[2], p_oc))

cat("\n── Cox IO ajustado ──\n")
data.frame(
  HR    = round(exp(coef(cox_opo_ajust)), 2),
  IC_lo = round(exp(confint(cox_opo_ajust))[,1], 2),
  IC_hi = round(exp(confint(cox_opo_ajust))[,2], 2),
  p     = round(summary(cox_opo_ajust)$coefficients[,5], 3)
) %>% print()

# ── IO segun genotipo mixto: gB, gH y gB/gH ──────────────────────────────────
library(dplyr)

# Funcion: descriptivo + test para una variable de mezcla
io_por_genotipo <- function(data, var_mix, etiqueta) {
  d <- data %>%
    mutate(grupo = factor(.data[[var_mix]], levels = c(0, 1),
                          labels = c("\u00danicos", "Mixtos"))) %>%
    filter(!is.na(grupo))
  
  tab <- table(d$grupo, d$evento_opo)
  
  resumen <- d %>%
    group_by(grupo) %>%
    summarise(
      n     = n(),
      n_opo = sum(evento_opo),
      pct   = round(100 * n_opo / n, 1),
      .groups = "drop"
    )
  
  # Fisher (robusto con n pequenos) y chi-cuadrado de Yates
  p_fisher <- fisher.test(tab)$p.value
  p_chisq  <- suppressWarnings(chisq.test(tab)$p.value)
  
  cat(sprintf("\n══════ %s ══════\n", etiqueta))
  print(resumen)
  cat(sprintf("p (Fisher)        = %s\n",
              sub("\\.", ",", sprintf("%.3f", p_fisher))))
  cat(sprintf("p (chi-cuadrado)  = %s\n",
              sub("\\.", ",", sprintf("%.3f", p_chisq))))
  
  invisible(resumen)
}

# Aplicar a las tres clasificaciones (usa df_opo, que ya tiene evento_opo)
io_por_genotipo(df_opo, "gbm",   "gB mixtos")
io_por_genotipo(df_opo, "ghm",   "gH mixtos")
io_por_genotipo(df_opo, "gbghm", "gB/gH mixtos")


## ---- [7g2] Figura 38. Distribucion etiologica de las infecciones oportunistas ----
# Apartado 5.7.3.

# ── Etiquetas personalizadas (edita aquí a tu gusto) ─────────────────────────
etiquetas_custom <- c(
  "VVZ"                    = "Herpes Zóster",
  "Candidiasis"            = "Candidiasis invasiva",
  "IFI"                    = "Infección fúngica invasiva",
  "VEB/PTLD"               = "PTLD asociado a VEB",
  "TBC"                    = "Tuberculosis diseminada",
  "Cryptosporidium"        = "Criptosporidiosis",
  "VHH8/Sarcoma de Kaposi" = "Sarcoma de Kaposi (VHH-8)",
  "Virus JC"               = "LMP por virus JC"
)

figura_opo_barras <- ggplot(df_counts %>% arrange(n),
                            aes(x = reorder(categoria, n), y = n,
                                fill = categoria)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = paste0("n=", n, " (",
                       gsub("\\.", ",", as.character(pct)), "%)")),
    hjust = -0.1, size = 3.8, color = "grey20"
  )  +
  scale_fill_manual(values = colores, labels = etiquetas_custom) +
  scale_x_discrete(labels = etiquetas_custom) +  # ← etiquetas eje Y
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.45))
  ) +
  coord_flip() +
  labs(x = NULL, y = "Número de episodios") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    axis.text.y        = element_text(size = 11, color = "grey20"),
    axis.text.x        = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20"),
    axis.line          = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    plot.margin        = margin(10, 20, 10, 10)
  )

figura_opo_barras


## ---- [7h1] Infecciones por virus respiratorios comunitarios: analisis ----
# Apartado 5.7.3.

library(readxl)
VR <- read_excel("data/base_clinica.xlsx", 
                 sheet = "VR")
head(VR)

library(tidyverse)
library(readxl)


VR <- read_excel("data/base_clinica.xlsx", sheet = "VR") %>%
  select(-gbghm, -genotipo) %>%
  left_join(pacientes %>% select(NHC, gbghm), by = "NHC") %>%
  mutate(año = as.integer(año),
         mes = as.integer(mes),
         genotipo_label = ifelse(gbghm == 0, "Genotipos únicos", "Genotipos mixtos"),
         virus = factor(virus, levels = c("SARS-CoV-2", "Influenza", "VRS",
                                          "Rhinovirus", "Metapneumovirus",
                                          "Parainfluenza", "Adenovirus",
                                          "Coronavirus estacional")))

        
cat("=== PREVALENCIA INFECCIÓN VR ===\n")
cat("Total pacientes con inf_VR:", sum(pacientes$inf_VR, na.rm = TRUE), "\n")
cat("Prevalencia global:", round(mean(pacientes$inf_VR, na.rm = TRUE)*100, 1), "%\n\n")

prev_geno <- pacientes %>%
  group_by(gbghm) %>%
  summarise(n     = n(),
            n_inf = sum(inf_VR, na.rm = TRUE),
            prev  = round(n_inf / n * 100, 1), .groups = "drop")
print(prev_geno)

cat("\nTest Fisher prevalencia VR por genotipo:\n")
print(fisher.test(table(pacientes$gbghm, pacientes$inf_VR)))

# 2. NÚMERO DE EPISODIOS
cat("\n=== EPISODIOS POR PACIENTE ===\n")

n_ep_pac <- VR %>% count(NHC, gbghm, genotipo_label, name = "n_ep")

n_ep_all <- pacientes %>%
  select(NHC, gbghm) %>%
  left_join(n_ep_pac %>% select(NHC, n_ep), by = "NHC") %>%
  replace_na(list(n_ep = 0)) %>%
  mutate(genotipo_label = ifelse(gbghm == 0, "Genotipos únicos", "Genotipos mixtos"))

cat("Total episodios - únicos:", sum(n_ep_all$n_ep[n_ep_all$gbghm == 0]), "\n")
cat("Total episodios - mixtos:", sum(n_ep_all$n_ep[n_ep_all$gbghm == 1]), "\n")
cat("Total episodios global:",   sum(n_ep_all$n_ep), "\n\n")

cat("Mediana episodios (pacientes infectados):\n")
n_ep_pac %>%
  group_by(genotipo_label) %>%
  summarise(mediana = median(n_ep),
            p25     = quantile(n_ep, 0.25),
            p75     = quantile(n_ep, 0.75),
            .groups = "drop") %>%
  print()

cat("\nTest Wilcoxon episodios (todos los pacientes, incluyendo 0):\n")
print(wilcox.test(n_ep ~ gbghm, data = n_ep_all))

# 3. DISTRIBUCIÓN POR VIRUS
cat("\n=== DISTRIBUCIÓN POR VIRUS ===\n")
dist_virus <- VR %>%
  count(genotipo_label, virus, .drop = FALSE) %>%
  group_by(genotipo_label) %>%
  mutate(total = sum(n),
         pct   = round(n / total * 100, 1)) %>%
  ungroup()
print(dist_virus, n = 30)

# 4. ANÁLISIS TEMPORAL
cat("\n=== ANÁLISIS TEMPORAL ===\n")
VR_con_año <- VR %>% filter(!is.na(año))
n_post  <- sum(VR_con_año$año >= 2020)
n_total <- nrow(VR_con_año)
cat("Episodios 2020-2022:", n_post, "/", n_total,
    "(", round(n_post / n_total * 100, 1), "%)\n")

# 5. COVID-19
cat("\n=== COVID-19 ===\n")
covid_nhc <- VR %>% filter(virus == "SARS-CoV-2") %>% distinct(NHC)
pacientes <- pacientes %>%
  mutate(covid = ifelse(NHC %in% covid_nhc$NHC, 1, 0))

pacientes %>%
  group_by(gbghm) %>%
  summarise(n       = n(),
            n_covid = sum(covid),
            pct     = round(n_covid / n * 100, 1), .groups = "drop") %>%
  print()

cat("\nTest Fisher COVID por genotipo:\n")
print(fisher.test(table(pacientes$gbghm, pacientes$covid)))

cat("\nNeumonía COVID por genotipo:\n")
pacientes %>%
  group_by(gbghm) %>%
  summarise(n          = n(),
            n_neumonia = sum(neumoniaCOVID, na.rm = TRUE),
            pct        = round(n_neumonia / n * 100, 1),
            .groups    = "drop") %>%
  print()

cat("\n=== VR Y RECHAZO AGUDO ===\n")

cat("Mediana episodios VR según RA:\n")
n_ep_all %>%
  left_join(pacientes %>% select(NHC, RA), by = "NHC") %>%
  group_by(RA) %>%
  summarise(n       = n(),
            mediana = median(n_ep),
            p25     = quantile(n_ep, 0.25),
            p75     = quantile(n_ep, 0.75),
            .groups = "drop") %>%
  print()

print(wilcox.test(n_ep ~ RA,
                  data = n_ep_all %>%
                    left_join(pacientes %>% select(NHC, RA), by = "NHC")))

cat("\nOR infección VR ~ RA en genotipos únicos:\n")
pac_unicos <- pacientes %>% filter(gbghm == 0)
print(fisher.test(table(pac_unicos$inf_VR, pac_unicos$RA)))

cat("\nOR infección VR ~ RA en genotipos mixtos:\n")
pac_mixtos <- pacientes %>% filter(gbghm == 1)
print(fisher.test(table(pac_mixtos$inf_VR, pac_mixtos$RA)))

# 7. VR + RA → DCIP
cat("\n=== VR + RA → DCIP ===\n")
pacientes %>%
  mutate(grupo = case_when(
    inf_VR == 1 & RA == 1 ~ "VR + RA",
    inf_VR == 0 & RA == 0 ~ "Sin VR ni RA",
    inf_VR == 1 & RA == 0 ~ "Solo VR",
    inf_VR == 0 & RA == 1 ~ "Solo RA"
  )) %>%
  group_by(grupo) %>%
  summarise(n      = n(),
            n_DCIP = sum(DCIP, na.rm = TRUE),
            pct    = round(n_DCIP / n * 100, 1),
            .groups = "drop") %>%
  print()


## ---- [7h2] Figura 39. Distribucion anual de los episodios por virus respiratorios ----
# Apartado 5.7.3.

virus_colors <- c(
  "SARS-CoV-2"             = "#EE6363",
  "Influenza"              = "dodgerblue2",
  "VRS"                    = "seagreen3",
  "Rhinovirus"             = "#9B5DE5",
  "Metapneumovirus"        = "#F77F00",
  "Parainfluenza"          = "#00E5EE",
  "Adenovirus"             = "#FF69B4",
  "Coronavirus estacional" = "#ADB5BD"
)

años_rng <- 2010:2022

df_todos <- VR_con_año %>%
  mutate(grupo_label = "Todos los pacientes (n = 146)") %>%
  count(año, virus, grupo_label)

df_geno <- VR_con_año %>%
  mutate(grupo_label = ifelse(gbghm == 0,
                              "Genotipos gB+gH únicos (n = 95)",
                              "Genotipos gB+gH mixtos (n = 51)")) %>%
  count(año, virus, grupo_label)

df_fig <- bind_rows(df_todos, df_geno) %>%
  mutate(grupo_label = factor(grupo_label,
                              levels = c("Todos los pacientes (n = 146)",
                                         "Genotipos gB+gH únicos (n = 95)",
                                         "Genotipos gB+gH mixtos (n = 51)")))

df_completo <- expand.grid(
  año         = años_rng,
  virus       = factor(names(virus_colors), levels = names(virus_colors)),
  grupo_label = factor(levels(df_fig$grupo_label), levels = levels(df_fig$grupo_label))
) %>%
  left_join(df_fig, by = c("año", "virus", "grupo_label")) %>%
  replace_na(list(n = 0))

fig_VR <- ggplot(df_completo,
                 aes(x = año, y = n, fill = virus)) +
  geom_col(width = 0.75, position = "stack") +
  geom_vline(xintercept = 2019.5, linetype = "dashed",
             color = "gray40", linewidth = 0.45) +
  annotate("text", x = 2019.65, y = Inf,
           label = "Pandemia\nCOVID-19",
           hjust = 0, vjust = 1.3,
           size = 3.8, color = "gray40") +
  facet_wrap(~ grupo_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = virus_colors, breaks = names(virus_colors)) +
  scale_x_continuous(breaks = años_rng, labels = as.character(años_rng)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Año",
       y = "Núm. infecciones",
       fill = NULL) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13),
    axis.text.x      = element_text(size = 11),
    axis.text.y      = element_text(size = 11),
    axis.title       = element_text(size = 12),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    panel.grid       = element_blank()
  ) +
  guides(fill = guide_legend(
    nrow = 2,
    byrow = TRUE,
    keywidth  = 0.8,
    keyheight = 0.8
  ))

print(fig_VR)

ggsave("output/Figura_R37.png", fig_VR,
       width = 27, height = 16, units = "cm",
       dpi = 300)


# Para cada paciente con RA y VR, ver qué ocurrió antes
pacientes %>%
  filter(RA == 1, inf_VR == 1) %>%
  select(NHC, F_RA1, gbghm) %>%
  left_join(VR %>% select(NHC, mes, año) %>%
              mutate(F_VR_primera = as.Date(paste(año, mes, "01", sep="-"))) %>%
              group_by(NHC) %>%
              summarise(F_VR_primera = min(F_VR_primera, na.rm = TRUE)),
            by = "NHC") %>%
  mutate(F_RA1 = as.Date(F_RA1),
         orden = case_when(
           F_RA1 < F_VR_primera ~ "RA antes de VR",
           F_VR_primera < F_RA1 ~ "VR antes de RA",
           TRUE ~ "simultáneo"
         )) %>%
  count(orden)


## ---- [7i] Figura 40. Probabilidad libre de rechazo agudo ----
# Apartado 5.7.4.1. Kaplan-Meier y modelos de Cox para gB, gH y gB/gH.

library(survival)
library(survminer)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── 1. Anadir grupos gB y gH a nivel paciente sobre df_ra ────────────────────
df_ra <- df_ra %>%
  mutate(
    gbm_grupo   = factor(gbm,   levels = c(0,1),
                         labels = c("Genotipos \u00fanicos", "Genotipos mixtos")),
    ghm_grupo   = factor(ghm,   levels = c(0,1),
                         labels = c("Genotipos \u00fanicos", "Genotipos mixtos"))
    # gbghm_grupo ya existe
  )

surv_ra <- Surv(time = df_ra$t_analisis_anios, event = df_ra$evento_ra)

km_gb   <- survfit(surv_ra ~ gbm_grupo,   data = df_ra)
km_gh   <- survfit(surv_ra ~ ghm_grupo,   data = df_ra)
km_gbgh <- survfit(surv_ra ~ gbghm_grupo, data = df_ra)

# ── 2. p de log-rank ─────────────────────────────────────────────────────────
p_lr <- function(formula) {
  ld <- survdiff(formula, data = df_ra)
  1 - pchisq(ld$chisq, df = length(ld$n) - 1)
}
p_gb   <- p_lr(surv_ra ~ gbm_grupo)
p_gh   <- p_lr(surv_ra ~ ghm_grupo)
p_gbgh <- p_lr(surv_ra ~ gbghm_grupo)

fmt_p <- function(p) sub("\\.", ",", formatC(p, format = "f", digits = 3))

# ── Paleta y eje X ───────────────────────────────────────────────────────────
pal_km  <- c("#2E8B57", "firebrick3")   # verde unico / rojo mixto
xmax    <- max(df_ra$t_analisis_anios, na.rm = TRUE)
xbreaks <- seq(0, floor(xmax), 2)

# ── 3. Funcion de panel limpio ───────────────────────────────────────────────
panel_ra <- function(km_fit, p_val, titulo, mostrar_y) {
  gs <- ggsurvplot(
    km_fit,
    data              = df_ra,
    pval              = FALSE,
    conf.int          = TRUE,
    conf.int.alpha    = 0.08,
    conf.int.style    = "ribbon",
    risk.table        = TRUE,
    risk.table.height = 0.22,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 3.3,
    palette           = pal_km,
    legend.title      = "",
    legend.labs       = c("Genotipos \u00fanicos", "Genotipos mixtos"),
    xlab              = "Tiempo desde trasplante (a\u00f1os)",
    ylab              = if (mostrar_y) "Libre de rechazo agudo" else "",
    title             = titulo,
    risk.table.title  = "Pacientes en riesgo",
    ggtheme           = theme_classic(base_size = 12),
    font.x            = c(11, "plain", "black"),
    font.y            = c(11, "plain", "black"),
    font.tickslab     = c(10, "plain", "black"),
    xlim              = c(0, xmax),
    ylim              = c(0.7, 1)
  )
  
  gs$plot <- gs$plot +
    scale_x_continuous(breaks = xbreaks,
                       expand = expansion(mult = c(0.06, 0.02))) +
    annotate("text", x = 0.1, y = 0.735,
             label = paste0("italic('p')*' = '*'", fmt_p(p_val), "'"),
             parse = TRUE, hjust = 0, size = 3.6, color = "black") +
    theme(
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      legend.position  = "none"
    )
  
  gs$table <- gs$table +
    scale_x_continuous(breaks = xbreaks,
                       expand = expansion(mult = c(0.07, 0.02))) +
    labs(x = NULL, title = "Pacientes en riesgo") +
    theme(
      plot.title       = element_text(size = 10, face = "plain", color = "black"),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank()
    )
  
  gs
}

# ── 4. Construir los tres paneles ────────────────────────────────────────────
g_gb   <- panel_ra(km_gb,   p_gb,   "gB",    mostrar_y = TRUE)
g_gh   <- panel_ra(km_gh,   p_gh,   "gH",    mostrar_y = TRUE)
g_gbgh <- panel_ra(km_gbgh, p_gbgh, "gB/gH", mostrar_y = TRUE)

# ── 5. Bloque curva + tabla ──────────────────────────────────────────────────
bloque <- function(gs) (gs$plot / gs$table) + plot_layout(heights = c(4, 1))
b_gb   <- bloque(g_gb)
b_gh   <- bloque(g_gh)
b_gbgh <- bloque(g_gbgh)

# ── 6. Leyenda compartida en DOS FILAS ───────────────────────────────────────
leyenda_plot <- ggsurvplot(
  km_gbgh, data = df_ra, palette = pal_km,
  legend.title = "", legend.labs = c("Genotipos \u00fanicos", "Genotipos mixtos"),
  ggtheme = theme_classic(base_size = 12),
  font.legend = c(11, "plain", "black")
)$plot +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.position   = "bottom",
        legend.key.height = unit(14, "pt"),
        legend.key.width  = unit(20, "pt"),
        legend.spacing.y  = unit(2, "pt"),
        legend.text       = element_text(size = 11, margin = margin(l = 4)))

leyenda <- cowplot::get_legend(leyenda_plot)
celda_leyenda <- patchwork::wrap_elements(full = leyenda)

# ── 7. Composicion 2x2: gB | gH  /  gB/gH | leyenda ──────────────────────────
figura_R31_2x2 <-
  (b_gb   | b_gh) /
  (b_gbgh | celda_leyenda) &
  theme(plot.background = element_rect(fill = "white", color = NA))

print(figura_R31_2x2)

ggsave("output/Figura_R31_KM_rechazo_2x2.png",
       figura_R31_2x2,
       width = 24, height = 20, units = "cm", dpi = 300, bg = "white")


# ════════════════════════════════════════════════════════════════════════════
#  DATOS NUMERICOS — KM rechazo agudo por gB, gH y gB/gH
# ════════════════════════════════════════════════════════════════════════════

library(survival)
library(dplyr)

fmt_p <- function(p) ifelse(p < 0.001, "<0,001",
                            sub("\\.", ",", sprintf("%.3f", p)))

# Funcion que resume un analisis completo para una variable de mezcla
resumen_ra <- function(var_grupo, etiqueta) {
  
  d <- df_ra %>%
    mutate(grupo = factor(.data[[var_grupo]], levels = c(0, 1),
                          labels = c("\u00danicos", "Mixtos"))) %>%
    filter(!is.na(grupo))
  
  surv <- Surv(d$t_analisis_anios, d$evento_ra)
  
  cat(sprintf("\n══════════════ %s ══════════════\n", etiqueta))
  
  # ── Descriptivo por grupo ──
  cat("\n-- Rechazo agudo por grupo --\n")
  d %>%
    group_by(grupo) %>%
    summarise(
      n        = n(),
      n_evento = sum(evento_ra),
      pct      = round(100 * mean(evento_ra), 1),
      .groups  = "drop"
    ) %>% print()
  
  # ── Log-rank ──
  ld <- survdiff(surv ~ grupo, data = d)
  p_lr <- 1 - pchisq(ld$chisq, df = length(ld$n) - 1)
  cat(sprintf("\nLog-rank: chi2 = %.2f  p = %s\n", ld$chisq, fmt_p(p_lr)))
  
  # ── Supervivencia libre de rechazo a tiempos clave ──
  cat("\n-- Libre de rechazo (1, 2, 3, 5 anios) --\n")
  km <- survfit(surv ~ grupo, data = d)
  print(summary(km, times = c(1, 2, 3, 5)))
  
  # ── Cox crudo ──
  cox_c <- coxph(surv ~ grupo, data = d)
  hr_c  <- exp(coef(cox_c)[1])
  ci_c  <- exp(confint(cox_c)[1, ])
  p_c   <- summary(cox_c)$coefficients[1, 5]
  cat(sprintf("\nCox crudo:    HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              hr_c, ci_c[1], ci_c[2], fmt_p(p_c)))
  
  # ── Cox ajustado (edad + risk) ──
  cox_a <- coxph(surv ~ grupo + edad_tx + risk, data = d)
  hr_a  <- exp(coef(cox_a)[1])
  ci_a  <- exp(confint(cox_a)[1, ])
  p_a   <- summary(cox_a)$coefficients[1, 5]
  cat(sprintf("Cox ajustado: HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              hr_a, ci_a[1], ci_a[2], fmt_p(p_a)))
  
  invisible(NULL)
}

# ── Ejecutar para las tres clasificaciones ───────────────────────────────────
resumen_ra("gbm",   "gB mixtos")
resumen_ra("ghm",   "gH mixtos")
resumen_ra("gbghm", "gB/gH mixtos")


## ---- [7j] Disfuncion cronica del injerto pulmonar: incidencia y modelos de Cox ----
# Apartado 5.7.4.2. Resultados descritos en el texto, sin figura asociada.

# ════════════════════════════════════════════════════════════════════════════
#  DATOS DESCRIPTIVOS GENERALES — DCIP
# ════════════════════════════════════════════════════════════════════════════
library(dplyr)

# ── Total DCIP ──
cat("── Total DCIP ──\n")
cat(sprintf("n = %d (%.1f%%)\n",
            sum(df_dcip$DCIP == 1, na.rm = TRUE),
            100 * mean(df_dcip$DCIP == 1, na.rm = TRUE)))

# ── Fenotipo DCIP ──
cat("\n── Fenotipo DCIP ──\n")
table(pacientes$tipoDCI, useNA = "always") %>% print()

# ── Tiempo hasta DCIP (anios) ──
cat("\n── Tiempo hasta DCIP (anios) ──\n")
df_dcip %>%
  filter(DCIP == 1) %>%
  summarise(
    mediana = round(median(t_hasta_dcip_anios, na.rm = TRUE), 1),
    q1      = round(quantile(t_hasta_dcip_anios, 0.25, na.rm = TRUE), 1),
    q3      = round(quantile(t_hasta_dcip_anios, 0.75, na.rm = TRUE), 1)
  ) %>% print()

# ── Seguimiento global (anios) ──
cat("\n── Seguimiento global (anios) ──\n")
df_dcip %>%
  summarise(
    mediana = round(median(t_analisis_anios, na.rm = TRUE), 1),
    q1      = round(quantile(t_analisis_anios, 0.25, na.rm = TRUE), 1),
    q3      = round(quantile(t_analisis_anios, 0.75, na.rm = TRUE), 1),
    max     = round(max(t_analisis_anios, na.rm = TRUE), 1)
  ) %>% print()

# ── DCIP por grupo (las 3 clasificaciones) ──
for (v in c("gbm", "ghm", "gbghm")) {
  cat(sprintf("\n── DCIP por grupo (%s) ──\n", v))
  df_dcip %>%
    mutate(grupo = factor(.data[[v]], levels = c(0,1),
                          labels = c("Unicos", "Mixtos"))) %>%
    filter(!is.na(grupo)) %>%
    group_by(grupo) %>%
    summarise(
      n      = n(),
      n_dcip = sum(DCIP == 1, na.rm = TRUE),
      pct    = round(100 * mean(DCIP == 1, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>% print()
}

# ── Temporalidad: DNAemia mixta antes de DCIP ──
cat("\n── Temporalidad: DNAemia mixta antes de DCIP ──\n")
primera_mixta <- dnaemias %>%
  filter(gbghmix == 1) %>%
  mutate(Fini = as.Date(Fini)) %>%
  group_by(NHC) %>%
  summarise(primera_mixta = min(Fini), .groups = "drop")

df_dcip %>%
  filter(DCIP == 1) %>%
  left_join(primera_mixta, by = "NHC") %>%
  mutate(
    orden = case_when(
      is.na(primera_mixta)   ~ "Sin DNAemia mixta",
      primera_mixta < F_DCIP ~ "DNAemia antes de DCIP",
      TRUE                   ~ "DNAemia despues de DCIP"
    ),
    t_entre = as.numeric(F_DCIP - primera_mixta) / 365.25
  ) %>%
  group_by(gbghm_grupo, orden) %>%
  summarise(
    n       = n(),
    mediana = round(median(t_entre, na.rm = TRUE), 1),
    q1      = round(quantile(t_entre, 0.25, na.rm = TRUE), 1),
    q3      = round(quantile(t_entre, 0.75, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>% print()


## ---- [7k1] Neoplasias: estadisticos descriptivos ----
# Apartado 5.7.5.1.

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(patchwork)

# ── 1. Preparar datos ─────────────────────────────────────────────────────────
df_onco <- pacientes %>%
  mutate(
    Ftx     = as.Date(Ftx),
    F_onco  = as.Date(F_onco),
    Ffinseg = as.Date(Ffinseg),
    fecha_fin_real   = if_else(onco == 1, F_onco, Ffinseg),
    t_analisis_dias  = as.numeric(fecha_fin_real - Ftx),
    t_analisis_anios = t_analisis_dias / 365.25,
    t_hasta_onco_anios = as.numeric(F_onco - Ftx) / 365.25,
    evento_onco      = if_else(onco == 1, 1L, 0L),
    # Tumor sólido no cutáneo
    # Tumor sólido no cutáneo (incluye SOL y SOL-P = 20 casos)
    solido         = as.integer(onco == 1 & onco_tipo == "SOL"),
    fecha_fin_sol  = if_else(solido == 1, F_onco, Ffinseg),
    t_analisis_sol = as.numeric(fecha_fin_sol - Ftx) / 365.25,
    evento_solido  = if_else(solido == 1, 1L, 0L),
    risk_f           = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo      = factor(gbghm, levels = c(0,1),
                              labels = c("Genotipos únicos",
                                         "Genotipos mixtos"))
  )

# ── 2. Estadísticos generales ─────────────────────────────────────────────────
cat("── Total neoplasias ──\n")
cat(sprintf("n = %d (%.1f%%)\n",
            sum(df_onco$onco == 1, na.rm = TRUE),
            100 * mean(df_onco$onco == 1, na.rm = TRUE)))

cat("\n── Tipo de neoplasia ──\n")
table(pacientes$onco_tipo, useNA = "always") %>% print()

cat("\n── Detalle onco_tipo2 ──\n")
table(pacientes$onco_tipo2, useNA = "always") %>% print()

cat("\n── Tiempo hasta neoplasia (años) ──\n")
df_onco %>%
  filter(onco == 1) %>%
  summarise(
    mediana = median(t_hasta_onco_anios, na.rm = TRUE),
    q1      = quantile(t_hasta_onco_anios, 0.25, na.rm = TRUE),
    q3      = quantile(t_hasta_onco_anios, 0.75, na.rm = TRUE)
  ) %>% print()

cat("\n── Seguimiento global (años) ──\n")
df_onco %>%
  summarise(
    mediana = median(t_analisis_anios, na.rm = TRUE),
    q1      = quantile(t_analisis_anios, 0.25, na.rm = TRUE),
    q3      = quantile(t_analisis_anios, 0.75, na.rm = TRUE),
    max     = max(t_analisis_anios, na.rm = TRUE)
  ) %>% print()

cat("\n── Neoplasia por grupo genotipo ──\n")
df_onco %>%
  group_by(gbghm_grupo) %>%
  summarise(
    n       = n(),
    n_onco  = sum(onco == 1, na.rm = TRUE),
    pct     = 100 * mean(onco == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

cat("\n── Tumores sólidos por grupo ──\n")
df_onco %>%
  group_by(gbghm_grupo) %>%
  summarise(
    n        = n(),
    n_solido = sum(solido == 1, na.rm = TRUE),
    pct      = 100 * mean(solido == 1, na.rm = TRUE),
    .groups  = "drop"
  ) %>% print()

# ── 3. Tabaquismo como confusor ───────────────────────────────────────────────
cat("\n── Tabaquismo ──\n")
cat(sprintf("Exfumadores: %d (%.1f%%)\n",
            sum(pacientes$exfum == 1, na.rm = TRUE),
            100 * mean(pacientes$exfum == 1, na.rm = TRUE)))

cat("\n── Tabaquismo vs genotipos mixtos ──\n")
cat(sprintf("p: %.3f\n",
            fisher.test(table(df_onco$exfum,
                              df_onco$gbghm))$p.value))

cat("\n── Cox tabaco vs neoplasia ──\n")
surv_tab <- Surv(df_onco$t_analisis_anios, df_onco$evento_onco)
cox_tab  <- coxph(surv_tab ~ exfum, data = df_onco)
hr_t  <- exp(coef(cox_tab)[1])
ci_t  <- exp(confint(cox_tab)[1,])
p_t   <- summary(cox_tab)$coefficients[1,5]
cat(sprintf("HR tabaco vs neoplasia: %.2f (%.2f–%.2f)  p = %.3f\n",
            hr_t, ci_t[1], ci_t[2], p_t))

cat("\n── Cox tabaco vs tumor sólido ──\n")
surv_sol <- Surv(df_onco$t_analisis_sol, df_onco$evento_solido)
cox_tab_sol <- coxph(surv_sol ~ exfum, data = df_onco)
hr_ts <- exp(coef(cox_tab_sol)[1])
ci_ts <- exp(confint(cox_tab_sol)[1,])
p_ts  <- summary(cox_tab_sol)$coefficients[1,5]
cat(sprintf("HR tabaco vs tumor sólido: %.2f (%.2f–%.2f)  p = %.3f\n",
            hr_ts, ci_ts[1], ci_ts[2], p_ts))

# ── 4. KM + log-rank neoplasia global ────────────────────────────────────────
surv_onco    <- Surv(df_onco$t_analisis_anios, df_onco$evento_onco)
km_onco      <- survfit(surv_onco ~ gbghm_grupo, data = df_onco)
lr_onco      <- survdiff(surv_onco ~ gbghm_grupo, data = df_onco)
p_onco       <- 1 - pchisq(lr_onco$chisq, df = length(lr_onco$n) - 1)
cat(sprintf("\n── Log-rank neoplasia global: p = %.3f ──\n", p_onco))

# ── 5. Cox neoplasia global ───────────────────────────────────────────────────
cox_onco_crudo <- coxph(surv_onco ~ gbghm_grupo, data = df_onco)
cox_onco_ajust <- coxph(surv_onco ~ gbghm_grupo + edad_tx + risk_f,
                        data = df_onco)

cat("\n── Cox neoplasia global crudo ──\n")
hr_oc <- exp(coef(cox_onco_crudo)[1])
ci_oc <- exp(confint(cox_onco_crudo)[1,])
p_oc  <- summary(cox_onco_crudo)$coefficients[1,5]
cat(sprintf("HR: %.2f (%.2f–%.2f)  p = %.3f\n", hr_oc, ci_oc[1], ci_oc[2], p_oc))

cat("\n── Cox neoplasia global ajustado ──\n")
data.frame(
  HR    = round(exp(coef(cox_onco_ajust)), 2),
  IC_lo = round(exp(confint(cox_onco_ajust))[,1], 2),
  IC_hi = round(exp(confint(cox_onco_ajust))[,2], 2),
  p     = round(summary(cox_onco_ajust)$coefficients[,5], 3)
) %>% print()

# ── 6. KM + Cox tumor sólido ─────────────────────────────────────────────────
km_sol      <- survfit(surv_sol ~ gbghm_grupo, data = df_onco)
lr_sol      <- survdiff(surv_sol ~ gbghm_grupo, data = df_onco)
p_sol       <- 1 - pchisq(lr_sol$chisq, df = length(lr_sol$n) - 1)
cat(sprintf("\n── Log-rank tumor sólido: p = %.3f ──\n", p_sol))

cox_sol_ajust <- coxph(surv_sol ~ gbghm_grupo + edad_tx + risk_f,
                       data = df_onco)

cat("\n── Cox tumor sólido ajustado ──\n")
data.frame(
  HR    = round(exp(coef(cox_sol_ajust)), 2),
  IC_lo = round(exp(confint(cox_sol_ajust))[,1], 2),
  IC_hi = round(exp(confint(cox_sol_ajust))[,2], 2),
  p     = round(summary(cox_sol_ajust)$coefficients[,5], 3)
) %>% print()

# ── 7. Supervivencia libre tumor sólido a 7 años ─────────────────────────────
cat("\n── Supervivencia libre tumor sólido a 7 años ──\n")
print(summary(km_sol, times = 7))

# ── 8. Temporalidad: DNAemia mixta antes de neoplasia ────────────────────────
cat("\n── Temporalidad: genotipos mixtos vs neoplasia ──\n")

primera_mixta_onco <- dnaemias %>%
  filter(gbghmix == 1) %>%
  mutate(Fini = as.Date(Fini)) %>%
  group_by(NHC) %>%
  summarise(primera_mixta = min(Fini), .groups = "drop")

df_onco %>%
  filter(onco == 1) %>%
  left_join(primera_mixta_onco, by = "NHC") %>%
  mutate(
    orden = case_when(
      is.na(primera_mixta)    ~ "Sin DNAemia mixta",
      primera_mixta < F_onco  ~ "DNAemia antes de neoplasia",
      TRUE                    ~ "DNAemia después de neoplasia"
    ),
    t_entre = as.numeric(F_onco - primera_mixta) / 365.25
  ) %>%
  group_by(gbghm_grupo, orden) %>%
  summarise(
    n       = n(),
    mediana = median(t_entre, na.rm = TRUE),
    q1      = quantile(t_entre, 0.25, na.rm = TRUE),
    q3      = quantile(t_entre, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

cat("\n── Temporalidad: genotipos mixtos vs tumor sólido ──\n")
df_onco %>%
  filter(solido == 1) %>%
  left_join(primera_mixta_onco, by = "NHC") %>%
  mutate(
    orden = case_when(
      is.na(primera_mixta)    ~ "Sin DNAemia mixta",
      primera_mixta < F_onco  ~ "DNAemia antes de tumor sólido",
      TRUE                    ~ "DNAemia después de tumor sólido"
    ),
    t_entre = as.numeric(F_onco - primera_mixta) / 365.25
  ) %>%
  group_by(gbghm_grupo, orden) %>%
  summarise(
    n       = n(),
    mediana = median(t_entre, na.rm = TRUE),
    q1      = quantile(t_entre, 0.25, na.rm = TRUE),
    q3      = quantile(t_entre, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()


## ---- [7k2] Figuras 41 y 42. Supervivencia libre de neoplasia y de tumor solido ----
# Apartado 5.7.5.2 y 5.7.5.3. Kaplan-Meier para gB, gH y gB/gH.

library(survival)
library(survminer)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── 1. Anadir grupos gB y gH a nivel paciente sobre df_onco ──────────────────
df_onco <- df_onco %>%
  mutate(
    gbm_grupo = factor(gbm, levels = c(0,1),
                       labels = c("Genotipos \u00fanicos", "Genotipos mixtos")),
    ghm_grupo = factor(ghm, levels = c(0,1),
                       labels = c("Genotipos \u00fanicos", "Genotipos mixtos"))
    # gbghm_grupo ya existe
  )

surv_onco <- Surv(df_onco$t_analisis_anios, df_onco$evento_onco)

km_gb   <- survfit(surv_onco ~ gbm_grupo,   data = df_onco)
km_gh   <- survfit(surv_onco ~ ghm_grupo,   data = df_onco)
km_gbgh <- survfit(surv_onco ~ gbghm_grupo, data = df_onco)

# ── 2. p de log-rank ─────────────────────────────────────────────────────────
p_lr <- function(formula) {
  ld <- survdiff(formula, data = df_onco)
  1 - pchisq(ld$chisq, df = length(ld$n) - 1)
}
p_gb   <- p_lr(surv_onco ~ gbm_grupo)
p_gh   <- p_lr(surv_onco ~ ghm_grupo)
p_gbgh <- p_lr(surv_onco ~ gbghm_grupo)

fmt_p <- function(p) sub("\\.", ",", formatC(p, format = "f", digits = 3))

# ── Paleta y eje X ───────────────────────────────────────────────────────────
pal_km  <- c("#2E8B57", "firebrick3")
xmax    <- max(df_onco$t_analisis_anios, na.rm = TRUE)
xbreaks <- seq(0, floor(xmax), 2)

# ── 3. Funcion de panel limpio ───────────────────────────────────────────────
panel_onco <- function(km_fit, p_val, titulo, mostrar_y) {
  gs <- ggsurvplot(
    km_fit,
    data              = df_onco,
    pval              = FALSE,
    conf.int          = TRUE,
    conf.int.alpha    = 0.08,
    conf.int.style    = "ribbon",
    risk.table        = TRUE,
    risk.table.height = 0.22,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 3.3,
    palette           = pal_km,
    legend.title      = "",
    legend.labs       = c("Genotipos \u00fanicos", "Genotipos mixtos"),
    xlab              = "Tiempo desde trasplante (a\u00f1os)",
    ylab              = if (mostrar_y) "Libre de neoplasia" else "",
    title             = titulo,
    risk.table.title  = "Pacientes en riesgo",
    ggtheme           = theme_classic(base_size = 12),
    font.x            = c(11, "plain", "black"),
    font.y            = c(11, "plain", "black"),
    font.tickslab     = c(10, "plain", "black"),
    xlim              = c(0, xmax),
    ylim              = c(0, 1)
  )
  
  gs$plot <- gs$plot +
    scale_x_continuous(breaks = xbreaks,
                       expand = expansion(mult = c(0.06, 0.02))) +
    annotate("text", x = 0.1, y = 0.08,
             label = paste0("italic('p')*' = '*'", fmt_p(p_val), "'"),
             parse = TRUE, hjust = 0, size = 3.6, color = "black") +
    theme(
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      legend.position  = "none"
    )
  
  gs$table <- gs$table +
    scale_x_continuous(breaks = xbreaks,
                       expand = expansion(mult = c(0.06, 0.02))) +
    labs(x = NULL, title = "Pacientes en riesgo") +
    theme(
      plot.title       = element_text(size = 10, face = "plain", color = "black"),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank()
    )
  
  gs
}

# ── 4. Construir los tres paneles ────────────────────────────────────────────
g_gb   <- panel_onco(km_gb,   p_gb,   "gB",    mostrar_y = TRUE)
g_gh   <- panel_onco(km_gh,   p_gh,   "gH",    mostrar_y = TRUE)
g_gbgh <- panel_onco(km_gbgh, p_gbgh, "gB/gH", mostrar_y = TRUE)

# ── 5. Bloque curva + tabla ──────────────────────────────────────────────────
bloque <- function(gs) (gs$plot / gs$table) + plot_layout(heights = c(4, 1))
b_gb   <- bloque(g_gb)
b_gh   <- bloque(g_gh)
b_gbgh <- bloque(g_gbgh)

# ── 6. Leyenda compartida en DOS FILAS ───────────────────────────────────────
leyenda_plot <- ggsurvplot(
  km_gbgh, data = df_onco, palette = pal_km,
  legend.title = "", legend.labs = c("Genotipos \u00fanicos", "Genotipos mixtos"),
  ggtheme = theme_classic(base_size = 12),
  font.legend = c(11, "plain", "black")
)$plot +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.position   = "bottom",
        legend.key.height = unit(14, "pt"),
        legend.key.width  = unit(20, "pt"),
        legend.spacing.y  = unit(2, "pt"),
        legend.text       = element_text(size = 11, margin = margin(l = 4)))

leyenda <- cowplot::get_legend(leyenda_plot)
celda_leyenda <- patchwork::wrap_elements(full = leyenda)

# ── 7. Composicion 2x2 ───────────────────────────────────────────────────────
figura_R34_2x2 <-
  (b_gb   | b_gh) /
  (b_gbgh | celda_leyenda) &
  theme(plot.background = element_rect(fill = "white", color = NA))

print(figura_R34_2x2)

ggsave("output/Figura_R34_KM_neoplasia_2x2.png",
       figura_R34_2x2,
       width = 24, height = 20, units = "cm", dpi = 300, bg = "white")


## ---- [7k3] Tabla 22. Riesgo de neoplasia y de tumor solido no cutaneo ----
# Apartado 5.7.5.2. Hazard ratios crudas y ajustadas.

# ════════════════════════════════════════════════════════════════════════════
#  DATOS NUMERICOS — Neoplasia global y tumor solido por gB, gH y gB/gH
# ════════════════════════════════════════════════════════════════════════════
library(survival)
library(dplyr)

fmt_p2 <- function(p) ifelse(p < 0.001, "<0,001",
                             sub("\\.", ",", sprintf("%.3f", p)))

# Funcion generica: desenlace = nombre de la var evento; t_var = tiempo
resumen_onco <- function(var_grupo, evento, t_var, etiqueta) {
  d <- df_onco %>%
    mutate(grupo = factor(.data[[var_grupo]], levels = c(0,1),
                          labels = c("\u00danicos", "Mixtos"))) %>%
    filter(!is.na(grupo))
  surv <- Surv(d[[t_var]], d[[evento]])
  
  cat(sprintf("\n══════ %s ══════\n", etiqueta))
  cat("\n-- Eventos por grupo --\n")
  d %>% group_by(grupo) %>%
    summarise(n = n(), n_evento = sum(.data[[evento]]),
              pct = round(100 * mean(.data[[evento]]), 1), .groups = "drop") %>%
    print()
  
  ld   <- survdiff(surv ~ grupo, data = d)
  p_lr <- 1 - pchisq(ld$chisq, df = length(ld$n) - 1)
  cat(sprintf("\nLog-rank: p = %s\n", fmt_p2(p_lr)))
  
  cox_c <- coxph(surv ~ grupo, data = d)
  hr_c  <- exp(coef(cox_c)[1]); ci_c <- exp(confint(cox_c)[1, ])
  p_c   <- summary(cox_c)$coefficients[1, 5]
  cat(sprintf("Cox crudo:    HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              hr_c, ci_c[1], ci_c[2], fmt_p2(p_c)))
  
  cox_a <- coxph(surv ~ grupo + edad_tx + risk_f, data = d)
  hr_a  <- exp(coef(cox_a)[1]); ci_a <- exp(confint(cox_a)[1, ])
  p_a   <- summary(cox_a)$coefficients[1, 5]
  cat(sprintf("Cox ajustado: HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              hr_a, ci_a[1], ci_a[2], fmt_p2(p_a)))
}

cat("\n########## NEOPLASIA GLOBAL ##########\n")
resumen_onco("gbm",   "evento_onco", "t_analisis_anios", "gB mixtos")
resumen_onco("ghm",   "evento_onco", "t_analisis_anios", "gH mixtos")
resumen_onco("gbghm", "evento_onco", "t_analisis_anios", "gB/gH mixtos")

cat("\n########## TUMOR SOLIDO ##########\n")
resumen_onco("gbm",   "evento_solido", "t_analisis_sol", "gB mixtos")
resumen_onco("ghm",   "evento_solido", "t_analisis_sol", "gH mixtos")
resumen_onco("gbghm", "evento_solido", "t_analisis_sol", "gB/gH mixtos")


## ---- [7l1] Mortalidad: preparacion de los datos ----
# Apartado 5.7.6.

library(dplyr); library(survival); library(survminer)
library(ggplot2); library(patchwork); library(cmprsk)

df_mort <- pacientes %>%
  mutate(
    Ftx_dt = as.Date(Ftx), Fexitus_dt = as.Date(Fexitus),
    Ffinseg_dt = as.Date(Ffinseg), Fnac_dt = as.Date(Fnac),
    fecha_fin_real   = if_else(exitus == 1, Fexitus_dt, Ffinseg_dt),
    t_analisis_anios = as.numeric(fecha_fin_real - Ftx_dt) / 365.25,
    evento_mort      = if_else(exitus == 1, 1L, 0L),
    edad_tx          = as.numeric(Ftx_dt - Fnac_dt) / 365.25,
    risk_f           = factor(risk, levels = c("MEDIO", "ALTO")),
    gbm_grupo   = factor(gbm,   levels = c(0,1), labels = c("Genotipos \u00fanicos","Genotipos mixtos")),
    ghm_grupo   = factor(ghm,   levels = c(0,1), labels = c("Genotipos \u00fanicos","Genotipos mixtos")),
    gbghm_grupo = factor(gbghm, levels = c(0,1), labels = c("Genotipos \u00fanicos","Genotipos mixtos"))
  ) %>%
  filter(!is.na(t_analisis_anios), t_analisis_anios > 0)


## ---- [7l2] Figura 43. Supervivencia global segun genotipos mixtos ----
# Apartado 5.7.6. Kaplan-Meier para gB, gH y gB/gH.

surv_mort <- Surv(df_mort$t_analisis_anios, df_mort$evento_mort)
km_gb   <- survfit(surv_mort ~ gbm_grupo,   data = df_mort)
km_gh   <- survfit(surv_mort ~ ghm_grupo,   data = df_mort)
km_gbgh <- survfit(surv_mort ~ gbghm_grupo, data = df_mort)

p_lr <- function(f) { ld <- survdiff(f, data = df_mort); 1 - pchisq(ld$chisq, df = length(ld$n) - 1) }
p_gb   <- p_lr(surv_mort ~ gbm_grupo)
p_gh   <- p_lr(surv_mort ~ ghm_grupo)
p_gbgh <- p_lr(surv_mort ~ gbghm_grupo)

fmt_p <- function(p) sub("\\.", ",", formatC(p, format = "f", digits = 3))
pal_km  <- c("#2E8B57", "firebrick3")
xmax    <- max(df_mort$t_analisis_anios, na.rm = TRUE)
xbreaks <- seq(0, floor(xmax), 2)

panel_mort <- function(km_fit, p_val, titulo, mostrar_y) {
  gs <- ggsurvplot(
    km_fit, data = df_mort, pval = FALSE,
    conf.int = TRUE, conf.int.alpha = 0.08, conf.int.style = "ribbon",
    risk.table = TRUE, risk.table.height = 0.22, risk.table.y.text = FALSE,
    risk.table.fontsize = 3.3, palette = pal_km,
    legend.title = "", legend.labs = c("Genotipos \u00fanicos","Genotipos mixtos"),
    xlab = "Tiempo desde trasplante (a\u00f1os)",
    ylab = if (mostrar_y) "Supervivencia global" else "",
    title = titulo, risk.table.title = "Pacientes en riesgo",
    ggtheme = theme_classic(base_size = 12),
    font.x = c(11,"plain","black"), font.y = c(11,"plain","black"),
    font.tickslab = c(10,"plain","black"),
    xlim = c(0, xmax), ylim = c(0, 1)
  )
  gs$plot <- gs$plot +
    scale_x_continuous(breaks = xbreaks, expand = expansion(mult = c(0.06, 0.02))) +
    annotate("text", x = 0.1, y = 0.08,
             label = paste0("italic('p')*' = '*'", fmt_p(p_val), "'"),
             parse = TRUE, hjust = 0, size = 3.6, color = "black") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.line = element_line(color = "black", linewidth = 0.4),
          axis.ticks = element_line(color = "black", linewidth = 0.3),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          legend.position = "none")
  gs$table <- gs$table +
    scale_x_continuous(breaks = xbreaks, expand = expansion(mult = c(0.06, 0.02))) +
    labs(x = NULL, title = "Pacientes en riesgo") +
    theme(plot.title = element_text(size = 10, color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank())
  gs
}

g_gb   <- panel_mort(km_gb,   p_gb,   "gB",    TRUE)
g_gh   <- panel_mort(km_gh,   p_gh,   "gH",    TRUE)
g_gbgh <- panel_mort(km_gbgh, p_gbgh, "gB/gH", TRUE)

bloque <- function(gs) (gs$plot / gs$table) + plot_layout(heights = c(4, 1))
b_gb <- bloque(g_gb); b_gh <- bloque(g_gh); b_gbgh <- bloque(g_gbgh)

leyenda <- cowplot::get_legend(
  ggsurvplot(km_gbgh, data = df_mort, palette = pal_km,
             legend.title = "", legend.labs = c("Genotipos \u00fanicos","Genotipos mixtos"),
             ggtheme = theme_classic(base_size = 12),
             font.legend = c(11,"plain","black"))$plot +
    guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.key.height = unit(14, "pt"), legend.key.width = unit(20, "pt"),
          legend.spacing.y = unit(2, "pt"),
          legend.text = element_text(size = 11, margin = margin(l = 4)))
)
celda_leyenda <- patchwork::wrap_elements(full = leyenda)

figura_mort_2x2 <- (b_gb | b_gh) / (b_gbgh | celda_leyenda) &
  theme(plot.background = element_rect(fill = "white", color = NA))
print(figura_mort_2x2)

ggsave("output/Figura_KM_mortalidad_2x2.png",
       figura_mort_2x2, width = 24, height = 20, units = "cm", dpi = 300, bg = "white")


## ---- [7l3] Mortalidad: hazard ratios y analisis de sensibilidad ----
# Apartado 5.7.6. Modelos crudo, ajustado y de riesgos competitivos.

fmt_p2 <- function(p) ifelse(p < 0.001, "<0,001", sub("\\.", ",", sprintf("%.3f", p)))

resumen_mort <- function(var_grupo, etiqueta) {
  d <- df_mort %>%
    mutate(grupo = factor(.data[[var_grupo]], levels = c(0,1),
                          labels = c("\u00danicos","Mixtos"))) %>%
    filter(!is.na(grupo))
  surv <- Surv(d$t_analisis_anios, d$evento_mort)
  
  cat(sprintf("\n══════════════ %s ══════════════\n", etiqueta))
  cat("\n-- Fallecidos por grupo --\n")
  d %>% group_by(grupo) %>%
    summarise(n = n(), n_exitus = sum(evento_mort),
              pct = round(100*mean(evento_mort),1), .groups="drop") %>% print()
  
  ld <- survdiff(surv ~ grupo, data = d)
  cat(sprintf("\nLog-rank: p = %s\n", fmt_p2(1 - pchisq(ld$chisq, df = length(ld$n)-1))))
  
  cat("\n-- Supervivencia (1,3,5,7,10 anios) --\n")
  print(summary(survfit(surv ~ grupo, data = d), times = c(1,3,5,7,10)))
  
  cc <- coxph(surv ~ grupo, data = d)
  cat(sprintf("\nCox crudo:       HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              exp(coef(cc)[1]), exp(confint(cc)[1,1]), exp(confint(cc)[1,2]),
              fmt_p2(summary(cc)$coefficients[1,5])))
  
  ca <- coxph(surv ~ grupo + edad_tx + serDR, data = d)
  cat(sprintf("Cox ajustado:    HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              exp(coef(ca)[1]), exp(confint(ca)[1,1]), exp(confint(ca)[1,2]),
              fmt_p2(summary(ca)$coefficients[1,5])))
  
  ds <- d %>% mutate(ev = if_else(evento_mort == 1 &
                                    grepl("COVID|covid", causa_exitus, ignore.case = TRUE),
                                  0L, evento_mort))
  cs <- coxph(Surv(ds$t_analisis_anios, ds$ev) ~ grupo + edad_tx + serDR, data = ds)
  cat(sprintf("Cox sens. COVID: HR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              exp(coef(cs)[1]), exp(confint(cs)[1,1]), exp(confint(cs)[1,2]),
              fmt_p2(summary(cs)$coefficients[1,5])))
  
  dfg <- d %>% mutate(fstatus = case_when(
    evento_mort == 0 ~ 0L,
    grepl("DCIP|cardiovascular|cardio", causa_exitus, ignore.case = TRUE) ~ 2L,
    TRUE ~ 1L))
  fg <- crr(ftime = dfg$t_analisis_anios, fstatus = dfg$fstatus,
            cov1 = model.matrix(~ grupo + edad_tx, data = dfg)[,-1])
  sfg <- summary(fg)
  cat(sprintf("Fine & Gray:     SHR = %.2f (%.2f\u2013%.2f)  p = %s\n",
              sfg$conf.int[1,1], sfg$conf.int[1,3], sfg$conf.int[1,4],
              fmt_p2(sfg$coef[1,5])))
}

resumen_mort("gbm",   "gB mixtos")
resumen_mort("ghm",   "gH mixtos")
resumen_mort("gbghm", "gB/gH mixtos")
