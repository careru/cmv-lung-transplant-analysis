# =============================================================================
# ESTUDIO 2. Caracterizacion de variantes de resistencia antiviral mediante NGS
# =============================================================================
#
# Tesis doctoral "Aspectos virologicos de la infeccion por citomegalovirus en
# receptores de trasplante pulmonar". Universidad Autonoma de Madrid, 2026.
#
# Cohorte: 11 receptores de trasplante pulmonar (RTP) con DNAemia por CMV
# >= 10.000 UI/mL entre 2022 y 2024. Genes UL97, UL54 y UL56 secuenciados con
# DeepChek CMV v2.0. Clasificacion de variantes segun CHARMD (Tilloy 2024) en
# mutaciones de resistencia antiviral (MRA), polimorfismos naturales (PN),
# variantes de significado incierto (VSI) y variantes no descritas (VND).
#
# -----------------------------------------------------------------------------
# CORRESPONDENCIA CON LA TESIS
# -----------------------------------------------------------------------------
#
#  Codigo   Apartado de la tesis                              Salida
#  ------   -----------------------------------------------   -------------------
#  [0]      6.1  Seleccion y caracteristicas de la cohorte     (depuracion)
#  [1]      6.2  Analisis global de variantes                  Recuentos del texto
#  [2]      6.2.1 Variantes asociadas a resistencia antiviral  Tabla 25
#  [3]      6.2.2 Polimorfismos naturales                      Anexo II
#  [4]      6.2.4 Variantes no descritas previamente           Figuras 44 a 46
#  [5]      6.2.4 Variantes no descritas previamente           Anexo III
#  [6]      6.2.5 Influencia del umbral de deteccion           Figura 47
#  [7]      6.2.6 Dinamica temporal de la resistencia en RTP-2 Figuras 48 y 49
#
# El apartado 6.2.3 (variantes de significado incierto) no genera salida
# grafica ni tabular: su unico contenido, la reclasificacion de D605E, procede
# de la clasificacion CHARMD del bloque [0e] y de revision bibliografica.
#
# -----------------------------------------------------------------------------
# METODOS ESTADISTICOS
# -----------------------------------------------------------------------------
#
# Estudio descriptivo. Las variantes se resumen como media +/- desviacion
# estandar de la frecuencia de variante, y como mediana con rango intercuartilico
# para la cobertura de secuenciacion y la carga viral. No se aplicaron
# contrastes de hipotesis.
#
# Filtros analiticos aplicados a todas las variantes:
#   - superar el filtro del sistema experto de DeepChek
#   - frecuencia de variante >= 5%
#   - cobertura total de la posicion >= 100 lecturas
#
# -----------------------------------------------------------------------------
# DATOS
# -----------------------------------------------------------------------------
#
# Este script no incluye datos de pacientes. Los numeros de historia clinica se
# han sustituido por constantes simbolicas (NHC_RTPxx) y las rutas locales por
# rutas relativas del repositorio. La estructura esperada de cada archivo de
# entrada se describe en data/README.md.
#
# Las secciones son plegables en RStudio (## ---- [x] ---- ). Cada seccion es
# autocontenida una vez ejecutado el bloque [0] de preparacion de datos.
# =============================================================================


## ---- [0] Preparacion de los datos ----
# Apartado 6.1 de la tesis: seleccion y caracteristicas de la cohorte.

# Identificadores de los 11 RTP incluidos. En el analisis original son numeros
# de historia clinica; aqui se sustituyen por constantes simbolicas.
NHC_RTP01 <- "RTP-01"; NHC_RTP02 <- "RTP-02"; NHC_RTP03 <- "RTP-03"
NHC_RTP04 <- "RTP-04"; NHC_RTP05 <- "RTP-05"; NHC_RTP06 <- "RTP-06"
NHC_RTP07 <- "RTP-07"; NHC_RTP08 <- "RTP-08"; NHC_RTP09 <- "RTP-09"
NHC_RTP10 <- "RTP-10"; NHC_RTP11 <- "RTP-11"



## ---- [0a] Union de los archivos CSV exportados de DeepChek ----
# Apartado 6.1. Lectura y fusion de los informes de secuenciacion.

library(tidyverse)
library(openxlsx)

# ── Configuración ─────────────────────────────────────────────────────────────
path        <- "data/ngs_csv"
output_file <- "output/mutaciones.xlsx"
extra_mra   <- c("L545F")

# ── Carga y fusión de archivos ────────────────────────────────────────────────
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
cat("Archivos encontrados:", length(files), "\n")

all_data <- map(files, ~ {
  df <- read_delim(.x, delim = ";", col_types = cols(.default = "c")) %>%
    rename(Npet = 1) %>%
    mutate(Npet = trimws(Npet))
  
  names(df) <- str_replace(
    names(df),
    "Is mutation interest \\(CMV variant database .*\\)",
    "Resistencia"
  )
  names(df) <- str_replace(names(df), "Passed Expert System", "Filtro")
  
  return(df)
}) %>%
  bind_rows()

cat("Muestras únicas (pre-filtro):", n_distinct(all_data$Npet), "\n")

# ── Selección de columnas ─────────────────────────────────────────────────────
all_data <- all_data %>%
  dplyr::select(
    Npet, Protein, Position, Mutation, Prevalence, Resistencia, Filtro,
    matches("Mutation count FW|Mutation\\.count\\.FW"),
    matches("Mutation count RV|Mutation\\.count\\.RV"),
    matches("Position coverage FW|Position\\.coverage\\.FW"),
    matches("Position coverage RV|Position\\.coverage\\.RV"),
    matches("^Q-Score$|^Q\\.Score$")
  )

# ── Renombrar a nombres cortos ────────────────────────────────────────────────
old_names <- names(all_data)
names(all_data)[str_detect(old_names, "Mutation.*count.*FW|Mutation\\.count\\.FW")]       <- "mut_fw"
names(all_data)[str_detect(old_names, "Mutation.*count.*RV|Mutation\\.count\\.RV")]       <- "mut_rv"
names(all_data)[str_detect(old_names, "Position.*coverage.*FW|Position\\.coverage\\.FW")] <- "cov_fw"
names(all_data)[str_detect(old_names, "Position.*coverage.*RV|Position\\.coverage\\.RV")] <- "cov_rv"
names(all_data)[str_detect(old_names, "Q.Score|Q-Score")]                                 <- "q_score"

# ── Convertir a numérico y calcular métricas derivadas ───────────────────────
all_data <- all_data %>%
  mutate(
    across(c(mut_fw, mut_rv, cov_fw, cov_rv, q_score), as.numeric),
    cobertura_total = cov_fw + cov_rv,
    balance_ok      = mut_fw > 0 & mut_rv > 0
  )

cat("Columnas de calidad disponibles:",
    paste(c("mut_fw","mut_rv","cov_fw","cov_rv","q_score",
            "cobertura_total","balance_ok"), collapse = ", "), "\n")

# ── Filtros ───────────────────────────────────────────────────────────────────
n_dash <- sum(all_data$Mutation == "-", na.rm = TRUE)
cat("Filas con Mutation == '-' eliminadas:", n_dash, "\n")

all_data <- all_data %>%
  filter(Mutation != "-") %>%
  filter(Filtro == "Yes")

cat("Muestras únicas (post-filtro Expert System):", n_distinct(all_data$Npet), "\n")
cat("Variantes totales (post-filtro Expert System):", nrow(all_data), "\n")

# ── Corrección de Resistencia ─────────────────────────────────────────────────
all_data <- all_data %>%
  mutate(Resistencia = if_else(str_trim(Mutation) %in% extra_mra, "Yes", Resistencia))
cat("Variantes marcadas como MRA adicional:",
    sum(str_trim(all_data$Mutation) %in% extra_mra), "\n")

nas <- sum(is.na(all_data$Resistencia))
cat("NAs en Resistencia:", nas, "\n")
if (nas > 0) {
  all_data <- all_data %>% mutate(Resistencia = replace_na(Resistencia, "No"))
  cat("NAs corregidos a 'No'. Comprobación:", sum(is.na(all_data$Resistencia)), "\n")
}

# ── Exportar ──────────────────────────────────────────────────────────────────
write.xlsx(all_data, file = output_file)
cat("Archivo exportado:", output_file, "\n")



## ---- [0b] Fusion con los datos clinicos ----
# Apartado 6.1. Asignacion de cada muestra a su paciente, carga viral y fecha.

library(haven)
library(data.table)

# ── Exclusiones ───────────────────────────────────────────────────────────────
npet_excluidos <- readRDS("data/npet_excluidos.rds")
# Identificadores de muestra excluidos por control de calidad.
# No se distribuyen: son identificadores internos de laboratorio.
cat("Exclusiones únicas definidas:", length(npet_excluidos), "\n")

# ── Importar datos clínicos ───────────────────────────────────────────────────
pacientes <- read_dta("data/muestras_clinicas.dta")
setDT(pacientes)
pacientes <- pacientes %>% mutate(Npet = as.character(Npet))

# ── Aplicar exclusiones antes del join ───────────────────────────────────────
n_antes        <- n_distinct(all_data$Npet)
all_data_clean <- all_data %>% filter(!Npet %in% npet_excluidos)
n_despues      <- n_distinct(all_data_clean$Npet)
cat("Muestras antes de exclusiones:", n_antes,             "\n")
cat("Muestras excluidas:            ", n_antes - n_despues, "\n")
cat("Muestras tras exclusiones:     ", n_despues,           "\n")

# ── Fusión con datos clínicos ─────────────────────────────────────────────────
all_data_nhc <- all_data_clean %>%
  left_join(pacientes %>% dplyr::select(Npet, NHC, CV, Freg), by = "Npet") %>%
  filter(!is.na(NHC))

# ── Diagnóstico ───────────────────────────────────────────────────────────────
cat("\n── Pacientes (allsamples.dta) ──────────────────\n")
cat("Npet únicos:", n_distinct(pacientes$Npet), "\n")
cat("NHC únicos: ", n_distinct(pacientes$NHC),  "\n")
cat("\n── all_data_nhc (tras exclusiones + join) ───────\n")
cat("Npet únicos:", n_distinct(all_data_nhc$Npet), "\n")
cat("NHC únicos: ", n_distinct(all_data_nhc$NHC),  "\n")

# ── Comprobación: no deben quedar muestras sin match ─────────────────────────
sin_match_restantes <- setdiff(unique(all_data_clean$Npet), unique(pacientes$Npet))
if (length(sin_match_restantes) > 0) {
  cat("\nAVISO: Muestras aún sin match en allsamples.dta:",
      paste(sin_match_restantes, collapse = ", "), "\n")
} else {
  cat("\nTodas las muestras restantes tienen match en allsamples.dta.\n")
}



## ---- [0c] Filtro analitico: sistema experto y frecuencia de variante >= 5% ----
# Apartado 6.1. Criterios de calidad aplicados a todas las variantes.

all_data_nhc <- all_data_nhc %>%
  mutate(Prevalence = as.numeric(Prevalence)) %>%
  filter(Filtro == "Yes", Prevalence >= 5, cobertura_total >= 100)

cat("── Tras filtro Prevalence ≥5% ──────────────────\n")
cat("Variantes totales:", nrow(all_data_nhc),          "\n")
cat("Npet únicos:      ", n_distinct(all_data_nhc$Npet), "\n")
cat("NHC únicos:       ", n_distinct(all_data_nhc$NHC),  "\n")

cat("\nDistribución por gen:\n")
all_data_nhc %>%
  count(Protein, name = "n_variantes") %>%
  arrange(Protein) %>%
  print()



## ---- [0d] Seleccion de la subcohorte de receptores de trasplante pulmonar ----
# Apartado 6.1. Restriccion a los 11 RTP incluidos en el estudio.

pacientes_txp <- c(
  NHC_RTP02, NHC_RTP11, NHC_RTP03, NHC_RTP04, NHC_RTP05,
  NHC_RTP06, NHC_RTP07, NHC_RTP08, NHC_RTP09, NHC_RTP10, NHC_RTP01
)

all_data_txp <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp)

cat("── Subcohorte TXP ──────────────────────────────\n")
cat("Pacientes TXP definidos:", length(pacientes_txp),         "\n")
cat("Pacientes TXP con datos:", n_distinct(all_data_txp$NHC),  "\n")
cat("Npet únicos TXP:        ", n_distinct(all_data_txp$Npet), "\n")
cat("Variantes totales TXP:  ", nrow(all_data_txp),            "\n")

cat("\nDistribución por gen (TXP):\n")
all_data_txp %>%
  count(Protein, name = "n_variantes") %>%
  arrange(Protein) %>%
  print()

sin_datos <- setdiff(pacientes_txp, unique(all_data_txp$NHC))
if (length(sin_datos) > 0) {
  cat("\nAVISO: Pacientes TXP sin variantes tras filtro:",
      paste(sin_datos, collapse = ", "), "\n")
} else {
  cat("\nTodos los pacientes TXP tienen datos tras el filtro.\n")
}



## ---- [0e] Clasificacion CHARMD: MRA, PN y VND ----
# Apartados 6.2.1 a 6.2.4. Asignacion de clase a cada variante segun CHARMD.

# ══════════════════════════════════════════════════════════════════════════════
# VECTORES CHARMD — extraídos de CHARMD_UL54/UL56/UL97.xlsx (Tilloy et al. 2024)
# ══════════════════════════════════════════════════════════════════════════════

# ── UL54 ──────────────────────────────────────────────────────────────────────
UL54_MRA <- c(
  "D301N", "E303G", "E303D", "N408D", "N408H", "N408K", "N408S", "N410K",
  "F412C", "F412L", "F412S", "F412V", "D413A", "D413E", "D413N", "D413del",
  "D413Y", "P488R", "N495K", "K500N", "L501F", "L501I", "T503A", "T503I",
  "A505V", "K513E", "K513N", "K513Q", "K513R", "D515E", "D515Y", "L516R",
  "L516W", "I521T", "P522A", "P522S", "P522T", "C524del", "V526L", "C539G",
  "C539R", "D542E", "A543V", "L545F", "L545S", "L545W", "T552S", "T552N",
  "Q578H", "Q578L", "S585A", "D588E", "D588N", "C590F", "F595I", "T700A",
  "V715A", "V715M", "I726T", "E756D", "E756K", "E756Q", "L773V", "L776M",
  "V781I", "V787A", "V787E", "V787L", "L802M", "K805Q", "A809V", "V812L",
  "T813S", "T821I", "P829S", "A834P", "T838A", "G841A", "G841S", "M844T",
  "M844V", "A928T", "V946L", "L957F", "D981del", "981-982del", "A987G"
)
UL54_VSI <- c(
  "Q229K", "E235G", "D247N", "D262N", "D271A", "D284E", "D288N", "N345S",
  "Y380C", "F396L", "L424V", "F460L", "V476G", "V482G", "A492D", "A543S",
  "R581H", "V654G", "S660G", "G667N", "A692G", "F718L", "F718S", "I726V",
  "E793V", "Q795P", "Q795R", "L802V", "G822D", "M827I", "M828V", "P859A",
  "L862F", "V902G", "K947E", "M959T"
)
UL54_POL <- c(
  "N4H", "N4Y", "S8N", "S8T", "A15S", "A17P", "R21C", "S24L", "S32P",
  "V95E", "H142Y", "G143S", "D163H", "T242K", "C256S", "C256L", "N267Y",
  "A269V", "A336T", "I341T", "G347D", "T351A", "V355A", "F357I", "P375L",
  "L394F", "K426R", "G441S", "V450G", "F460S", "S464F", "H465Y", "N467S",
  "A469T", "A473V", "V483A", "R512C", "D515G", "Y518C", "P522L", "Q541P",
  "G604S", "S612N", "P617S", "I622L", "A626V", "P628A", "P628L", "G629S",
  "A631G", "S633F", "V634A", "M640R", "A647S", "S651E", "G653S", "S655L",
  "S660N", "S663N", "F669L", "S676G", "G678S", "V681del", "N685S", "H686Y",
  "A688V", "G689R", "T691S", "T691A", "A692V", "A692S", "A692F", "A693T",
  "S695T", "Q697H", "I722V", "L737M", "V759M", "A786P", "R792S", "R800C",
  "R847H", "N855D", "H863R", "Q868R", "D870H", "V873L", "G874R", "in884T",
  "A885S", "A885T", "in885S", "in885T", "P887S", "L890F", "L890S", "T892I",
  "S897L", "N898D", "E899K", "E903G", "G920S", "V927M", "S932N", "V953A",
  "R954H", "S956Y", "A972V", "E989Q", "G993C", "T995K", "N998D", "S1000L",
  "D1005N", "R1006C", "R1006G", "A1012V", "L1020I", "R1052C", "D1055L",
  "R1070G", "Q1071R", "L1074M", "V1079I", "T1108A", "N1116H", "A1122T",
  "P1129Q", "G1133S", "A1138T", "S1146N", "N1147S", "R1149T", "G1151del",
  "P1153S", "L1156del", "S1162L", "S1235T", "C1241G", "A1108T", "T1122A",
  "L897S", "T885A"
)

# ── UL56 ──────────────────────────────────────────────────────────────────────
UL56_MRA <- c(
  "C25F", "V231A", "V231L", "N232Y", "V236A", "V236L", "V236M", "E237D",
  "L241P", "T244K", "L254F", "L257F", "L257I", "K258E", "F261C", "F261L",
  "Y321C", "C325F", "C325R", "C325Y", "C325W", "M329T", "A365S", "N368D",
  "R369G", "R369M", "R369S", "R369T", "R396S"
)
UL56_VSI <- c(
  "S229F", "L328V", "F345L"
)
UL56_POL <- c(
  "E2K", "V12L", "V33A", "A36V", "E75G", "L122P", "R129H", "A221V",
  "D242G", "A327V", "E393K", "A425V", "G436E", "M442T", "A444G", "S445N",
  "N446del", "N446S", "449-451del", "450-453del", "T452I", "S454N", "G460V",
  "A464V", "A464T", "V471A", "V476A", "E485G", "R507C", "H509N", "Q577R",
  "V582M", "D586N", "A648R", "A648S", "G649D", "L658F", "H698Q", "C721Y",
  "S749N", "V778A", "A779T", "A779V", "A788T", "V793A", "V793P", "P800L",
  "P803A", "Y806T", "T811P", "A812E", "A812G", "G813A", "A820Q", "L831M",
  "G840A", "S848N", "V425A", "N586D"
)

# ── UL97 ──────────────────────────────────────────────────────────────────────
UL97_MRA <- c(
  "L337M", "F342S", "F342Y", "G343A", "V353A", "K355del", "V356G", "K359E",
  "K359Q", "L397R", "L405P", "T409M", "H411L", "H411N", "H411Y", "D456N",
  "M460I", "M460T", "M460V", "V466G", "C480F", "C480R", "C518Y", "H520Q",
  "P521L", "590-600del", "591-594del", "591-607del", "C592F", "C592G",
  "A594E", "A594G", "A594P", "A594S", "A594T", "A594V", "L595F", "L595S",
  "L595W", "L595del", "595-603del", "E596G", "E596Y", "E596Q", "599-603del",
  "K599T", "601-602del", "601-603del", "C603W", "C603Y", "C607Y", "I610T",
  "A613V", "Y617del"
)
UL97_VSI <- c(
  "A427V", "V466M", "H469Y", "R476C", "A478V", "N510S", "M550I", "A582V",
  "A588V", "590-603del", "A591D", "A591V", "E596del", "N597del", "597-598del",
  "597-599del", "G598S", "K599E", "K599R", "L600I", "L600del", "T601M",
  "C603R", "C603S", "D605E", "C607F", "Y617H", "G623S", "E655K", "T659I",
  "V665I", "A674T"
)
UL97_POL <- c(
  "S7F", "L14F", "Q19E", "Q54E", "N68D", "N108S", "S108N", "R112C", "R112H",
  "D118N", "L126Q", "P132L", "S133P", "V171I", "D184G", "G215A", "A238T",
  "I244V", "P247S", "S249C", "F294C", "K316T", "E317A", "D329H", "V356W",
  "V356L", "K356V", "T429I", "A440V", "A442G", "Q449K", "Q449R", "M460L",
  "G561A", "H587Y", "M615V", "A639T", "A669V", "T75A", "Q126L", "V244I",
  "D68N", "S249N"
)

cat("\n── Vectores CHARMD (Tilloy et al. 2024) ────────\n")
cat(sprintf("UL54 → MRA: %d | VSI: %d | POL: %d\n",
            length(UL54_MRA), length(UL54_VSI), length(UL54_POL)))
cat(sprintf("UL56 → MRA: %d | VSI: %d | POL: %d\n",
            length(UL56_MRA), length(UL56_VSI), length(UL56_POL)))
cat(sprintf("UL97 → MRA: %d | VSI: %d | POL: %d\n",
            length(UL97_MRA), length(UL97_VSI), length(UL97_POL)))

# ── Asignación de clase (prioridad: MRA > VSI > POL > NEW) ───────────────────
all_data_nhc <- all_data_nhc %>%
  mutate(class = case_when(
    Protein == "UL54" & Mutation %in% UL54_MRA ~ "MRA",
    Protein == "UL54" & Mutation %in% UL54_VSI ~ "VSI",
    Protein == "UL54" & Mutation %in% UL54_POL ~ "POL",
    Protein == "UL56" & Mutation %in% UL56_MRA ~ "MRA",
    Protein == "UL56" & Mutation %in% UL56_VSI ~ "VSI",
    Protein == "UL56" & Mutation %in% UL56_POL ~ "POL",
    Protein == "UL97" & Mutation %in% UL97_MRA ~ "MRA",
    Protein == "UL97" & Mutation %in% UL97_VSI ~ "VSI",
    Protein == "UL97" & Mutation %in% UL97_POL ~ "POL",
    TRUE ~ "NEW"
  ))

# Regenerar all_data_txp con la columna class ya asignada
all_data_txp <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp)

cat("\nall_data_txp actualizado con class:",
    n_distinct(all_data_txp$NHC), "pacientes |",
    nrow(all_data_txp), "variantes\n")

cat("\n── Distribución de class por gen (todos) ───────\n")
all_data_nhc %>%
  count(Protein, class) %>%
  arrange(Protein, class) %>%
  print()

cat("\n── Variantes NEW por gen (no en CHARMD) ────────\n")
all_data_nhc %>%
  filter(class == "NEW") %>%
  count(Protein, Mutation, sort = TRUE) %>%
  print(n = Inf)



## ---- [0f] Etiquetas de paciente ----
# Apartado 6.1. Codificacion RTP-1 a RTP-11 usada en tablas y figuras.

rtp_map <- tibble::tribble(
  ~NHC,      ~RTP,
  NHC_RTP03,   "RTP-1",
  NHC_RTP02,  "RTP-2",
  NHC_RTP09,    "RTP-3",
  NHC_RTP05,    "RTP-4",
  NHC_RTP08,   "RTP-5",
  NHC_RTP07,  "RTP-6",
  NHC_RTP01,   "RTP-7",
  NHC_RTP10,  "RTP-8",
  NHC_RTP04,  "RTP-9",
  NHC_RTP06,  "RTP-10",
  NHC_RTP11,  "RTP-11"
)

all_data_nhc <- all_data_nhc %>% left_join(rtp_map, by = "NHC")
all_data_txp <- all_data_txp %>% left_join(rtp_map, by = "NHC")



## ---- [0g] Objetos de trabajo ----
# Apartado 6.1.

UL97_data <- all_data_nhc %>% filter(Protein == "UL97")
UL54_data <- all_data_nhc %>% filter(Protein == "UL54")
UL56_data <- all_data_nhc %>% filter(Protein == "UL56")



## ---- [0h] Muestras con resultado valido por gen ----
# Apartado 6.2. Denominadores de UL97 (26), UL54 (24) y UL56 (11).

npet_txp_total <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet) %>%
  nrow()
cat("Total peticiones únicas (Npet) en subcohorte TXP:", npet_txp_total, "\n\n")

cat("── Npet con resultado válido por gen ────────────────────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  group_by(Protein) %>%
  summarise(
    n_npet      = n_distinct(Npet),
    pct         = round(n_distinct(Npet) / npet_txp_total * 100, 1),
    n_pacientes = n_distinct(NHC),
    .groups     = "drop"
  ) %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  arrange(Protein) %>%
  print()

npet_ul97 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL97") %>% distinct(Npet)
npet_ul54 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL54") %>% distinct(Npet)
npet_ul56 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL56") %>% distinct(Npet)

cat("\n── Npet en UL97 pero no en UL54 (cobertura insuficiente u omitidas) ─────\n")
npet_solo_ul97 <- setdiff(npet_ul97$Npet, npet_ul54$Npet)
cat("N =", length(npet_solo_ul97), "\n")
cat(npet_solo_ul97, sep = "\n")

cat("\n── Npet en UL54 pero no en UL97 ─────────────────────────────────────────\n")
npet_solo_ul54 <- setdiff(npet_ul54$Npet, npet_ul97$Npet)
cat("N =", length(npet_solo_ul54), "\n")
cat(npet_solo_ul54, sep = "\n")

cat("\n── Npet con UL56 analizado ───────────────────────────────────────────────\n")
cat("N =", nrow(npet_ul56), "\n")
cat("Porcentaje sobre total TXP:",
    round(nrow(npet_ul56) / npet_txp_total * 100, 1), "%\n")

cat("\n── Desglose por Npet: genes analizados (TRUE/FALSE) ─────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet, NHC) %>%
  left_join(npet_ul97 %>% mutate(UL97 = TRUE), by = "Npet") %>%
  left_join(npet_ul54 %>% mutate(UL54 = TRUE), by = "Npet") %>%
  left_join(npet_ul56 %>% mutate(UL56 = TRUE), by = "Npet") %>%
  mutate(across(c(UL97, UL54, UL56), ~replace_na(., FALSE))) %>%
  left_join(
    all_data_nhc %>% filter(NHC %in% pacientes_txp) %>% distinct(Npet, RTP),
    by = "Npet"
  ) %>%
  arrange(RTP, Npet) %>%
  print(n = Inf)

# ── Diagnóstico calidad: resumen de cobertura y Q-Score en subcohorte TXP ────
cat("\n── Resumen de métricas de calidad (subcohorte TXP, Prevalence ≥5%) ──────\n")
all_data_txp %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  group_by(Protein) %>%
  summarise(
    cob_media   = round(mean(cobertura_total, na.rm = TRUE), 0),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_mediana   = round(median(q_score, na.rm = TRUE), 1),
    pct_balance = round(mean(balance_ok, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  rename(
    "Gen"             = Protein,
    "Cob. media (×)"  = cob_media,
    "Cob. mediana (×)"= cob_mediana,
    "Cob. mínima (×)" = cob_min,
    "Q-Score mediana" = q_mediana,
    "% balance FW/RV" = pct_balance
  ) %>%
  print()



## ---- [1] Recuento de variantes por umbral, gen y clase ----
# Apartado 6.2. Cifras globales citadas en el texto.

library(flextable)
library(officer)
library(dplyr)
library(purrr)

thresholds <- c(5, 10, 15, 20)
genes      <- c("UL97", "UL54", "UL56")

calcular_celda <- function(data, gene, clase, thr) {
  d <- data %>%
    filter(Protein == gene, class == clase, Prevalence >= thr)
  n_var <- n_distinct(d$Mutation)
  n_pac <- n_distinct(d$NHC)
  if (n_var == 0) return("0 / 0")
  paste0(n_var, " / ", n_pac)
}

# ── Construir dataframe ───────────────────────────────────────────────────────
filas <- map_dfr(genes, function(gene) {
  total_pac <- n_distinct(filter(all_data_txp, Protein == gene)$NHC)
  map_dfr(thresholds, function(thr) {
    tibble(
      Gen    = if (thr == thresholds[1]) paste0(gene, "  (n = ", total_pac, " pac.)") else "",
      Umbral = paste0("\u2265 ", thr, "%"),
      MRA    = calcular_celda(all_data_txp, gene, "MRA",   thr),
      POL    = calcular_celda(all_data_txp, gene, "POL",   thr),
      VSI    = calcular_celda(all_data_txp, gene, "VSI",   thr),
      Nuevas = calcular_celda(all_data_txp, gene, "NEW",   thr)
    )
  })
})

# ── Identificar filas con MRA > 0 para colorear ───────────────────────────────
filas_mra <- which(filas$MRA != "0 / 0")

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft <- flextable::flextable(filas) %>%
  
  # Cabeceras
  flextable::set_header_labels(
    Gen = "Gen", Umbral = "Umbral",
    MRA = "MRA", POL = "POL", VSI = "VSI", Nuevas = "Nuevas"
  ) %>%
  flextable::add_header_row(
    values    = c("", "", "Variantes \u00fanicas / Pacientes afectados", "", "", ""),
    colwidths = c(1, 1, 1, 1, 1, 1),
    top       = TRUE
  ) %>%
  flextable::merge_at(i = 1, j = 3:6, part = "header") %>%
  
  # Fusión Gen
  flextable::merge_at(i = 1:4,  j = 1, part = "body") %>%
  flextable::merge_at(i = 5:8,  j = 1, part = "body") %>%
  flextable::merge_at(i = 9:12, j = 1, part = "body") %>%
  
  # Alineación
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  
  # Tipografía base
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  
  # Cabecera: fondo azul oscuro, texto blanco
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  
  # Cuerpo: sombreado alternado por bloque gen
  flextable::bg(i = 1:4,  bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 5:8,  bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = 9:12, bg = "#EEF4FB", part = "body") %>%
  
  # Columna Umbral: fondo gris neutro para guía visual
  flextable::bg(j = 2, bg = "#F5F5F5", part = "body") %>%
  flextable::color(j = 2, color = "#555555", part = "body") %>%
  flextable::italic(j = 2, part = "body") %>%
  
  # MRA con color si > 0
  flextable::bold(i = filas_mra, j = 3, part = "body") %>%
  flextable::color(i = filas_mra, j = 3, color = "#B22222", part = "body") %>%
  
  # Bordes
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = c(4, 8), border = b_sep, part = "body") %>%
  flextable::hline(i = 1, border = b_sep, part = "header") %>%
  
  # Anchos
  flextable::width(j = 1,   width = 4.0, unit = "cm") %>%
  flextable::width(j = 2,   width = 1.8, unit = "cm") %>%
  flextable::width(j = 3:6, width = 2.8, unit = "cm") %>%
  
  # Padding
  flextable::padding(padding = 6, part = "all") %>%
  
  # Nota al pie
  flextable::add_footer_lines(
    "Valores expresados como: variantes \u00fanicas detectadas / pacientes con al menos una variante de esa clase al umbral indicado. MRA: mutaci\u00f3n de resistencia a antivirales; POL: polimorfismo natural; VSI: variante de significado incierto; Nuevas: variantes no clasificadas en CHARMD (Tilloy et al., 2024). An\u00e1lisis restringido a la subcohorte de trasplante pulmonar (n\u00famero de pacientes con datos por gen indicado entre par\u00e9ntesis)."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft



## ---- [2] Tabla 25. Mutaciones de resistencia antiviral ----
# Apartado 6.2.1. MRA detectadas en RTP-2 con frecuencia, cobertura y carga viral.

library(dplyr)
library(flextable)
library(officer)
library(purrr)
library(tibble)

nhc_rtp2 <- NHC_RTP02

# ── Datos fijos (CHARMD) ──────────────────────────────────────────────────────
charmd <- tibble(
  Protein    = c("UL97", "UL97", "UL97", "UL54"),
  Mutation   = c("A594V", "T409M", "H411Y", "T503I"),
  Resistencia_fenotipica = c(
    "Ganciclovir",
    "Maribavir",
    "Maribavir",
    "Ganciclovir y Cidofovir"
  )
)

# ── Calcular metricas desde all_data_txp ──────────────────────────────────────
mra_raw <- all_data_txp %>%
  filter(class == "MRA", Prevalence >= 5)

mra_stats <- mra_raw %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras   = n_distinct(Npet),
    prev_media   = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de      = round(sd(Prevalence,   na.rm = TRUE), 1),
    cv_mediana   = round(median(as.numeric(CV), na.rm = TRUE), 0),
    cv_p25       = round(quantile(as.numeric(CV), 0.25, na.rm = TRUE), 0),
    cv_p75       = round(quantile(as.numeric(CV), 0.75, na.rm = TRUE), 0),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(charmd, by = c("Protein", "Mutation"))

# ── Umbral minimo de deteccion ────────────────────────────────────────────────
umbral_min_mra <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "MRA", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(umbral_min = paste0("\u2265", min(umbral), "%"), .groups = "drop")

mra_stats <- mra_stats %>%
  left_join(umbral_min_mra, by = c("Protein", "Mutation"))

# ── Total muestras validas por gen en RTP-2 (calculado) ───────────────────────
n_muestras_gen <- all_data_nhc %>%
  filter(NHC == nhc_rtp2, Protein %in% c("UL97", "UL54")) %>%
  group_by(Protein) %>%
  summarise(n_total = n_distinct(Npet), .groups = "drop")

cat("── Denominadores de la Tabla 25 ──\n")
print(n_muestras_gen)

n_ul97 <- n_muestras_gen$n_total[n_muestras_gen$Protein == "UL97"]
n_ul54 <- n_muestras_gen$n_total[n_muestras_gen$Protein == "UL54"]

# ── Helper de formato numerico sin espacios sobrantes ─────────────────────────
fnum <- function(x) format(x, big.mark = ".", trim = TRUE)

mra_stats <- mra_stats %>%
  left_join(n_muestras_gen, by = "Protein") %>%
  mutate(
    pct   = round(n_muestras / n_total * 100, 1),
    col_n = paste0(n_muestras, " (",
                   format(pct, decimal.mark = ",", trim = TRUE), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ",", trim = TRUE), " (n = 1)"),
      paste0(format(prev_media, decimal.mark = ",", trim = TRUE),
             " \u00b1 ", format(prev_de, decimal.mark = ",", trim = TRUE))
    ),
    col_cv = if_else(
      is.na(cv_mediana),
      "\u2014",
      paste0(fnum(cv_mediana), " (", fnum(cv_p25), "\u2013", fnum(cv_p75), ")")
    ),
    col_cob = if_else(
      n_muestras == 1,
      paste0(fnum(cob_mediana), "\u00d7"),
      paste0(fnum(cob_mediana), "\u00d7 (", fnum(cob_p25), "\u2013",
             fnum(cob_p75), "; m\u00edn. ", fnum(cob_min), "\u00d7)")
    ),
    col_q = if_else(
      n_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "\u2013", q_p75, ")")
    ),
    col_balance = paste0(n_balance_ok, "/", n_muestras)
  )

# ── Orden explicito ───────────────────────────────────────────────────────────
orden_mra <- tibble(
  Protein  = c("UL97", "UL97", "UL97", "UL54"),
  Mutation = c("A594V", "T409M", "H411Y", "T503I")
)

filas_mra <- orden_mra %>%
  left_join(mra_stats, by = c("Protein", "Mutation")) %>%
  mutate(
    Gen = case_when(
      Mutation == "A594V" ~ paste0("UL97\n(n = ", n_ul97, " muestras)"),
      Mutation == "T503I" ~ paste0("UL54\n(n = ", n_ul54, " muestras)"),
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(
    Gen, Mutation, col_n, col_prev, col_cv,
    col_cob, col_q, col_balance, umbral_min, Resistencia_fenotipica
  )

names(filas_mra) <- c(
  "Gen", "MRA", "Muestras, n (%)",
  "Prevalencia (media \u00b1 DE, %)",
  "CV en muestras con MRA\nmediana (RIC), UI/mL",
  "Cobertura\n(mediana, RIC; m\u00edn.)",
  "Q-score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)",
  "Umbral m\u00ednimo de\ndetecci\u00f3n",
  "Resistencia fenot\u00edpica*"
)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_mra <- flextable::flextable(filas_mra) %>%
  flextable::merge_at(i = 1:3, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = 1:3, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 4,   bg = "#FFFFFF", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = 3, border = b_sep, part = "body") %>%
  flextable::width(j = 1,  width = 2.4, unit = "cm") %>%
  flextable::width(j = 2,  width = 1.8, unit = "cm") %>%
  flextable::width(j = 3,  width = 2.0, unit = "cm") %>%
  flextable::width(j = 4,  width = 2.5, unit = "cm") %>%
  flextable::width(j = 5,  width = 2.8, unit = "cm") %>%
  flextable::width(j = 6,  width = 3.0, unit = "cm") %>%
  flextable::width(j = 7,  width = 1.8, unit = "cm") %>%
  flextable::width(j = 8,  width = 2.1, unit = "cm") %>%
  flextable::width(j = 9,  width = 2.0, unit = "cm") %>%
  flextable::width(j = 10, width = 3.0, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(
    "*Resistencia fenot\u00edpica seg\u00fan la base de datos CHARMD (Tilloy et al., 2024). CV: carga viral; DE: desviaci\u00f3n est\u00e1ndar; RIC: rango intercuart\u00edlico. El umbral m\u00ednimo de detecci\u00f3n indica el umbral de frecuencia de variante m\u00e1s bajo al que la variante es detectable. Cobertura: n\u00famero total de lecturas forward + reverse que cubren la posici\u00f3n, expresado como mediana, rango intercuart\u00edlico y valor m\u00ednimo. Q-score: puntuaci\u00f3n de calidad de llamada. Balance FW/RV: n\u00famero de muestras con soporte en ambas cadenas sobre el total de muestras con la variante."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_mra



## ---- [3a] Polimorfismos naturales: preparacion de los datos ----
# Apartado 6.2.2.

all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein) %>%
  summarise(
    n_variantes = n_distinct(Mutation),
    n_pacientes = n_distinct(NHC),
    .groups = "drop"
  )

all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    NHC_list = paste(sort(unique(NHC)), collapse = ", "),
    n_pac    = n_distinct(NHC),
    .groups  = "drop"
  ) %>%
  arrange(Protein, desc(n_pac)) %>%
  print(n = Inf)



## ---- [3b] Anexo II. Tabla de polimorfismos naturales ----
# Apartado 6.2.2. PN identificados en UL97, UL54 y UL56.

# ── Datos ─────────────────────────────────────────────────────────────────────
pol_stats <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras   = n_distinct(Npet),
    n_pacientes  = n_distinct(NHC),
    prev_media   = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de      = round(sd(Prevalence, na.rm = TRUE), 1),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  )

# ── Umbral máximo de detección ────────────────────────────────────────────────
umbral_max_pol <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "POL", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    umbral_max = paste0("\u2265", max(umbral), "%"),
    .groups = "drop"
  )

# ── Total muestras y pacientes por gen ────────────────────────────────────────
n_total_gen <- tibble(
  Protein         = c("UL97", "UL54", "UL56"),
  n_tot_muestras  = c(26, 24, 11),
  n_tot_pacientes = c(11, 11, 3)
)

# ── Unir todo ─────────────────────────────────────────────────────────────────
pol_stats <- pol_stats %>%
  left_join(umbral_max_pol, by = c("Protein", "Mutation")) %>%
  left_join(n_total_gen,    by = "Protein") %>%
  mutate(
    pct_muestras  = round(n_muestras  / n_tot_muestras  * 100, 1),
    pct_pacientes = round(n_pacientes / n_tot_pacientes * 100, 1),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ","), " (n=1)"),
      paste0(format(prev_media, decimal.mark = ","),
             " \u00b1 ",
             format(prev_de, decimal.mark = ","))
    ),
    col_muestras  = paste0(n_muestras, " (",
                           format(pct_muestras,  decimal.mark = ","), "%)"),
    col_pacientes = paste0(n_pacientes, " (",
                           format(pct_pacientes, decimal.mark = ","), "%)"),
    col_cobertura = paste0(
      format(cob_mediana, big.mark = "."), "\u00d7",
      " (", format(cob_p25, big.mark = "."),
      "\u2013", format(cob_p75, big.mark = "."), "\u00d7;",
      " m\u00edn. ", format(cob_min, big.mark = "."), "\u00d7)"
    ),
    col_qscore = paste0(
      q_med, " (", q_p25, "\u2013", q_p75, ")"
    )
  )

# ── Orden ─────────────────────────────────────────────────────────────────────
pol_stats <- pol_stats %>%
  mutate(Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))) %>%
  arrange(Protein, desc(n_muestras))

# ── Dataframe para flextable ──────────────────────────────────────────────────
pol_ft <- pol_stats %>%
  group_by(Protein) %>%
  mutate(
    Gen = if_else(row_number() == 1,
                  paste0(as.character(Protein),
                         "\n(n = ", n_tot_muestras, " muestras)"),
                  "")
  ) %>%
  ungroup() %>%
  dplyr::select(Gen, Mutation, col_muestras, col_pacientes,
                col_prev, umbral_max, col_cobertura, col_qscore)

names(pol_ft) <- c(
  "Gen", "PN",
  "Muestras, n (%)",
  "Pacientes, n (%)",
  "Prevalencia\n(media \u00b1 DE, %)",
  "Umbral m\u00e1x.\ndetecci\u00f3n",
  "Cobertura\n(mediana, RIC; m\u00edn.)",
  "Q-Score\n(mediana, RIC)"
)

# ── Filas por bloque ──────────────────────────────────────────────────────────
n_ul97 <- sum(pol_stats$Protein == "UL97")
n_ul54 <- sum(pol_stats$Protein == "UL54")
n_ul56 <- sum(pol_stats$Protein == "UL56")

filas_ul97 <- 1:n_ul97
filas_ul54 <- (n_ul97 + 1):(n_ul97 + n_ul54)
filas_ul56 <- (n_ul97 + n_ul54 + 1):(n_ul97 + n_ul54 + n_ul56)
filas_sep  <- c(n_ul97, n_ul97 + n_ul54)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_pol <- flextable::flextable(pol_ft) %>%
  flextable::merge_at(i = filas_ul97, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul54, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul56, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul54, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 3.0, unit = "cm") %>%
  flextable::width(j = 2, width = 2.0, unit = "cm") %>%
  flextable::width(j = 3, width = 2.5, unit = "cm") %>%
  flextable::width(j = 4, width = 2.5, unit = "cm") %>%
  flextable::width(j = 5, width = 3.0, unit = "cm") %>%
  flextable::width(j = 6, width = 2.0, unit = "cm") %>%
  flextable::width(j = 7, width = 3.5, unit = "cm") %>%
  flextable::width(j = 8, width = 2.5, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "PN: polimorfismo natural. DE: desviaci\u00f3n est\u00e1ndar. RIC: rango intercuart\u00edlico. ",
    "Los porcentajes de muestras y pacientes se calculan sobre el total con datos v\u00e1lidos por gen ",
    "(UL97: 26 muestras, 11 pacientes; UL54: 24 muestras, 11 pacientes; UL56: 11 muestras, 3 pacientes). ",
    "El umbral m\u00e1ximo de detecci\u00f3n indica el umbral de frecuencia de variante m\u00e1s alto al que la variante es detectable. ",
    "Cobertura: suma de lecturas forward + reverse en la posici\u00f3n (mediana y RIC entre todas las muestras con la variante; m\u00ednimo observado). ",
    "Q-Score: mediana y RIC de la puntuaci\u00f3n de calidad de llamada de variante."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_pol



## ---- [4a] Variantes no descritas previamente: preparacion de los datos ----
# Apartado 6.2.4.

# ── Resumen actualizado de variantes NEW ──────────────────────────────────────
new_vars <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("Total NEW por gen:\n")
new_vars %>% count(Protein) %>% print()

cat("\nNEW por gen y umbral:\n")
map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "NEW", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_variantes = n_distinct(Mutation),
              n_pacientes = n_distinct(NHC), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  arrange(Protein, umbral) %>%
  print(n = Inf)






new_vars_dom <- new_vars %>%
  mutate(
    dominio = case_when(
      Protein == "UL97" & Position >= 335 & Position <= 370 ~ "P-loop",
      Protein == "UL97" & Position >= 409 & Position <= 411 ~ "Gatekeeper",
      Protein == "UL97" & Position >= 440 & Position <= 480 ~ "Loop catalítico",
      Protein == "UL97" & Position >= 520 & Position <= 607 ~ "Región canónica GCV",
      Protein == "UL54" & Position >= 408 & Position <= 412 ~ "Región IV",
      Protein == "UL54" & Position >= 501 & Position <= 545 ~ "Delta-C",
      Protein == "UL54" & Position >= 700 & Position <= 715 ~ "Región II",
      Protein == "UL54" & Position >= 725 & Position <= 964 ~ "Unión ATP/dNTP",
      Protein == "UL54" & Position >= 802 & Position <= 821 ~ "Región III",
      Protein == "UL54" & Position >= 910 & Position <= 913 ~ "Región I",
      Protein == "UL54" & Position >= 962 & Position <= 987 ~ "Región V",
      Protein == "UL54" & Position >= 1153 & Position <= 1159 ~ "NLS-A",
      Protein == "UL54" & Position >= 1213 & Position <= 1242 ~ "UL44 binding",
      Protein == "UL56" & Position == 25 ~ "N-terminal",
      Protein == "UL56" & Position >= 202 & Position <= 220 ~ "LAGLIDADG",
      Protein == "UL56" & Position >= 231 & Position <= 369 ~ "Región LMV",
      Protein == "UL56" & Position >= 514 & Position <= 522 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 671 & Position <= 680 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 757 & Position <= 764 ~ "Interacción UL89",
      TRUE ~ "Fuera de dominio funcional"
    ),
    en_dominio = dominio != "Fuera de dominio funcional"
  ) %>%
  arrange(Protein, desc(en_dominio), desc(n_pacientes), desc(prev_max))

cat("── En dominio funcional ─────────────────────────\n")
new_vars_dom %>% filter(en_dominio) %>% print(n = Inf)

cat("\n── Fuera de dominio pero relevantes ─────────────\n")
new_vars_dom %>%
  filter(!en_dominio, n_pacientes >= 2 | prev_max >= 15) %>%
  print(n = Inf)



## ---- [4b] Figura 44. Distribucion de VND en UL97 ----
# Apartado 6.2.4. Posicion y frecuencia de las VND sobre el gen UL97.

library(ggplot2)
library(dplyr)
library(patchwork)

dominios_ul97 <- data.frame(
  start = c(230),
  end   = c(620),
  label = c("Dominio quinasa"),
  fill  = c("lightskyblue")
)

hotspots_ul97 <- data.frame(
  start = c(335, 397, 519, 590),
  end   = c(370, 480, 521, 607),
  label = c("", "", "", ""),
  fill  = c("#B03060", "#8E44AD", "#E67E22", "#E67E22")
)

ul97_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL97") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul97 <- new_vars_dom %>%
  filter(Protein == "UL97") %>%
  left_join(ul97_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul97 %>%
  mutate(
    y_nivel    = nivel[as.character(umbral_max)],
    en_quinasa = Position >= 230 & Position <= 620,
    y_point    = ifelse(en_quinasa, -y_nivel, y_nivel),
    color      = ifelse(en_quinasa, "#C0392B", "#2980B9"),
    label_mut  = case_when(
      en_quinasa                        ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = -(3.0 + (within_rank - 1) * 0.60)
  ) %>%
  ungroup()

ymax_plot <- max(
  max(df_lab_azul$y_lab, na.rm = TRUE),
  max(abs(df_lab_rojo$y_lab), na.rm = TRUE)
) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul97 <- ggplot() +
  
  # Líneas de referencia — empiezan después de las etiquetas, acaban en aa 707
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 707, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul97,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 4.0, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul97,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = hotspots_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = -0.20, label = label),
            size = 3, color = "white", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–707 con tapas
  geom_segment(aes(x = 1,   xend = 707, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 707, xend = 707, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 707, y = -0.28, label = "707", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_quinasa, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-100, 730),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-ymax_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots <- ggplot(data.frame(
  x     = c(1.5, 3.7, 6),
  y     = c(1, 1, 1),
  label = c("GCV/MBV-R (335-370)",
            "MBV-R (397-480)",
            "GCV-R cod\u00f3n 520 + (590-607)"),
  color = c("#B03060", "#8E44AD", "#E67E22")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 9.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul97 <- p_ul97 / leyenda_hotspots +
  plot_layout(heights = c(10, 1))

figura_ul97
ggsave(
  "output/figura_lollipop_UL97.png",
  figura_ul97,
  width = 30, height = 15, units = "cm", dpi = 400,
  bg = "white"
)



## ---- [4c] Figura 45. Distribucion de VND en UL54 ----
# Apartado 6.2.4. Posicion y frecuencia de las VND sobre el gen UL54.

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Dominios UL54 ─────────────────────────────────────────────────────────────
dominios_ul54 <- data.frame(
  start = c(290, 696, 1150),
  end   = c(550, 990, 1243),
  label = c("Exonucleasa", "Polimerasa", "NLS + unión ppUL44"),
  fill  = c("#82E0AA", "#7FB3D3", "#E8A0D0")
)

# ── Hotspots UL54 ─────────────────────────────────────────────────────────────
hotspots_ul54 <- data.frame(
  start = c(300, 408, 498, 985,   598, 802,   700,   750, 781, 945),
  end   = c(325, 420, 545, 990,   600, 825,   715,   760, 790, 955),
  fill  = c(rep("#E67E22", 4),
            rep("#8E44AD", 2),
            "#E74C3C",
            rep("#C0392B", 3))
)

# ── Datos UL54 con umbral_max ─────────────────────────────────────────────────
ul54_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL54") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul54 <- new_vars_dom %>%
  filter(Protein == "UL54") %>%
  left_join(ul54_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

# ── Regiones de resistencia ───────────────────────────────────────────────────
en_region_resistencia <- function(pos) {
  (pos >= 300 & pos <= 325) |
    (pos >= 408 & pos <= 420) |
    (pos >= 498 & pos <= 545) |
    (pos >= 598 & pos <= 600) |
    (pos >= 700 & pos <= 715) |
    (pos >= 750 & pos <= 760) |
    (pos >= 781 & pos <= 790) |
    (pos >= 802 & pos <= 825) |
    (pos >= 945 & pos <= 955) |
    (pos >= 985 & pos <= 990)
}

# ── Niveles fijos ─────────────────────────────────────────────────────────────
nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)

df <- ul54 %>%
  mutate(
    en_resist = en_region_resistencia(as.numeric(Position)),
    y_nivel   = nivel[as.character(umbral_max)],
    y_point   = ifelse(en_resist, -y_nivel, y_nivel),
    color     = ifelse(en_resist, "#C0392B", "#2980B9"),
    label_mut = case_when(
      en_resist                         ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

# ── Etiquetas con agrupación por proximidad ───────────────────────────────────
df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_resist) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 40))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_resist) %>%
  arrange(Position) %>%
  mutate(
    escalon = case_when(
      Position < 817  ~ (row_number() - 1) %% 3,
      TRUE            ~ (row_number() - 1 - sum(Position < 817)) %% 3
    ),
    y_lab = -(3.0 + escalon * 0.60)
  )%>%
  ungroup()

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5
ymin_plot <- min(df_lab_rojo$y_lab, na.rm = TRUE) - 0.5

# ── Niveles de referencia ─────────────────────────────────────────────────────
niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

# ── Plot UL54 ─────────────────────────────────────────────────────────────────
p_ul54 <- ggplot() +
  
  # Líneas de referencia — de -55 a aa 1243
  geom_segment(data = niveles_ref,
               aes(x = -55, xend = 1243, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul54,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.85, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul54 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.23, label = label),
            size = 4.0, color = "grey15", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul54,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–1243 con tapas
  geom_segment(aes(x = 1,    xend = 1243, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,    xend = 1,    y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 1243, xend = 1243, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,    y = -0.28, label = "1",    size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 1243, y = -0.28, label = "1243", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_resist, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-130, 1380),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(ymin_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 50)
  )

# ── Leyenda hotspots ──────────────────────────────────────────────────────────
leyenda_hotspots_ul54 <- ggplot(data.frame(
  x     = c(1.0, 4.2, 7.4, 10.6),
  y     = 1,
  label = c("GCV/CDV", "GCV/CDV/FOS", "FOS > GCV", "FOS"),
  color = c("#E67E22", "#8E44AD", "#E74C3C", "#C0392B")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.18, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 13)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 50)
  )

# ── Composición final ─────────────────────────────────────────────────────────
figura_ul54 <- p_ul54 / leyenda_hotspots_ul54 +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

figura_ul54

ggsave(
  "output/figura_lollipop_UL54.png",
  figura_ul54,
  width = 34, height = 18, units = "cm", dpi = 400
)



## ---- [4d] Figura 46. Distribucion de VND en UL56 ----
# Apartado 6.2.4. Posicion y frecuencia de las VND sobre el gen UL56.

library(ggplot2)
library(dplyr)
library(patchwork)

dominios_ul56 <- data.frame(
  start = c(4,   117, 280, 514, 671, 757),
  end   = c(20,  130, 310, 522, 680, 764),
  label = c("NLS", "endonucleasa", "Dedo de zinc", "Interacción pUL89", "Estab. 1", "Estab. 2"),
  fill  = c("cadetblue3", "cadetblue3", "cadetblue3", "#F0A500", "#F0A500", "#F0A500")
)

hotspots_ul56 <- data.frame(
  start = c(20,  230, 392),
  end   = c(30,  370, 398),
  label = c("",  "", ""),
  fill  = c("#B03060", "#B03060", "#B03060")
)

ul56_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL56") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul56 <- new_vars_dom %>%
  filter(Protein == "UL56") %>%
  left_join(ul56_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul56 %>%
  mutate(
    y_nivel   = nivel[as.character(umbral_max)],
    en_lmv    = Position >= 231 & Position <= 412,
    y_point   = ifelse(en_lmv, -y_nivel, y_nivel),
    color     = ifelse(en_lmv, "#C0392B", "#2980B9"),
    label_mut = Mutation
  )

df_lab_azul_base <- df %>%
  filter(!is.na(label_mut), !en_lmv) %>%
  arrange(Position) %>%
  mutate(row_idx = row_number())

df_lab_azul <- df_lab_azul_base %>%
  mutate(
    y_lab_raw = 2.5 + (row_idx - 1) * paso,
    y_lab_raw = case_when(
      Position <= 30 ~ y_lab_raw + 0.8,
      TRUE           ~ y_lab_raw
    ),
    y_lab = case_when(
      Position >= 595 ~ y_lab_raw - 1.6,
      TRUE            ~ y_lab_raw
    )
  )

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_lmv) %>%
  arrange(Position) %>%
  mutate(y_lab = -(0.9 + (row_number() - 1) * paso))

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul56 <- ggplot() +
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 850, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  # Etiquetas prevalencia — más cerca del aa 1 (x = -55, límite -75)
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  # Dominios
  geom_rect(data = dominios_ul56,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul56 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 3.5, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5, lineheight = 0.85) +
  # Hotspots
  geom_rect(data = hotspots_ul56,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  # Backbone aa 1–850
  geom_segment(aes(x = 1, xend = 850, y = 0, yend = 0),
               color = "grey10", linewidth = 1.0) +
  # Tapas extremos (palitos)
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 850, xend = 850, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_lmv, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 3.5, color = "grey10", vjust = 1) +
  annotate("text", x = 850, y = -0.28, label = "850", size = 3.5, color = "grey10", vjust = 1) +
  scale_x_continuous(
    limits = c(-100, 870),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-2.5, ymax_plot)) +
  # Sin título de panel — el nombre "UL56" va en el pie de figura
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots_ul56 <- ggplot(data.frame(
  x     = 1.5,
  y     = 1,
  label = "LMV-R",
  color = "#B03060"
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 3.5, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 4)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul56 <- p_ul56 / leyenda_hotspots_ul56 +
  plot_layout(heights = c(10, 1))

figura_ul56

ggsave(
  "output/figura_lollipop_UL56.png",
  figura_ul56,
  width = 28, height = 15, units = "cm", dpi = 400,
  bg = "white"
)



## ---- [5] Anexo III. Tabla combinada de VND por dominio y hotspot ----
# Apartado 6.2.4. Anotacion funcional de las 101 VND en los tres genes.

library(dplyr)
library(flextable)
library(officer)

# ── Anotación de región de interés (prioridad: hotspot > dominio) ─────────────
annot_ul97 <- function(pos) dplyr::case_when(
  pos >= 335 & pos <= 370  ~ "Hotspot GCV/MBV-R (aa 335\u2013370)",
  pos >= 397 & pos <= 480  ~ "Hotspot MBV-R (aa 397\u2013480)",
  pos >= 519 & pos <= 521  ~ "Hotspot GCV-R (aa 519\u2013521)",
  pos >= 590 & pos <= 607  ~ "Hotspot GCV-R (aa 590\u2013607)",
  pos >= 230 & pos <= 620  ~ "Dominio quinasa",
  TRUE                     ~ ""
)
annot_ul54 <- function(pos) dplyr::case_when(
  pos >= 300 & pos <= 325  ~ "Hotspot GCV/CDV (aa 300\u2013325)",
  pos >= 408 & pos <= 420  ~ "Hotspot GCV/CDV (aa 408\u2013420)",
  pos >= 498 & pos <= 545  ~ "Hotspot GCV/CDV (aa 498\u2013545)",
  pos >= 985 & pos <= 990  ~ "Hotspot GCV/CDV (aa 985\u2013990)",
  pos >= 598 & pos <= 600  ~ "Hotspot GCV/CDV/FOS (aa 598\u2013600)",
  pos >= 802 & pos <= 825  ~ "Hotspot GCV/CDV/FOS (aa 802\u2013825)",
  pos >= 700 & pos <= 715  ~ "Hotspot FOS/GCV (aa 700\u2013715)",
  pos >= 750 & pos <= 760  ~ "Hotspot FOS (aa 750\u2013760)",
  pos >= 781 & pos <= 790  ~ "Hotspot FOS (aa 781\u2013790)",
  pos >= 945 & pos <= 955  ~ "Hotspot FOS (aa 945\u2013955)",
  pos >= 290 & pos <= 550  ~ "Dominio exonucleasa",
  pos >= 696 & pos <= 990  ~ "Dominio polimerasa",
  pos >= 1150 & pos <= 1243 ~ "Regi\u00f3n NLS / uni\u00f3n ppUL44",
  TRUE                     ~ ""
)
annot_ul56 <- function(pos) dplyr::case_when(
  pos >= 230 & pos <= 370  ~ "Hotspot LMV-R (aa 230\u2013370)",
  TRUE                     ~ ""
)

# ── Formateadores (coma decimal, punto de millar) ─────────────────────────────
pc  <- function(x) formatC(x, format = "f", digits = 1, decimal.mark = ",")
mil <- function(x) formatC(x, format = "f", digits = 0, big.mark = ".", decimal.mark = ",")

# ── Denominadores (de tu objeto n_total_gen, bloque [3b]) ─────────────────────
n_total_gen <- tibble(
  Protein         = c("UL97", "UL54", "UL56"),
  n_tot_muestras  = c(26, 24, 11),
  n_tot_pacientes = c(11, 11, 3)
)

# ── Construcción de la tabla ──────────────────────────────────────────────────
tabla_vnd <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence,   na.rm = TRUE), 1),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(n_total_gen, by = "Protein") %>%
  mutate(
    pct_muestras  = round(n_muestras  / n_tot_muestras  * 100, 1),
    pct_pacientes = round(n_pacientes / n_tot_pacientes * 100, 1),
    col_muestras  = paste0(n_muestras,  " (", pc(pct_muestras),  "%)"),
    col_pacientes = paste0(n_pacientes, " (", pc(pct_pacientes), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(pc(prev_media), " (n=1)"),
      paste0(pc(prev_media), " \u00b1 ", pc(prev_de))
    ),
    col_cob = if_else(
      n_muestras == 1,
      paste0(mil(cob_mediana), "\u00d7"),
      paste0(mil(cob_mediana), "\u00d7 (", mil(cob_p25), "\u2013", mil(cob_p75), ")")
    ),
    dom_hot = case_when(
      Protein == "UL97" ~ annot_ul97(Position),
      Protein == "UL54" ~ annot_ul54(Position),
      Protein == "UL56" ~ annot_ul56(Position),
      TRUE              ~ ""
    ),
    Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  ) %>%
  arrange(Protein, Position)

# ── Diagnóstico (verifica que son 101) ────────────────────────────────────────
cat("Total VND (todas):", nrow(tabla_vnd), "\n")
tabla_vnd %>% count(Protein) %>% print()
cat("Con Dominio/Hotspot asignado:", sum(tabla_vnd$dom_hot != ""), "\n")

# ── Dataframe para flextable ──────────────────────────────────────────────────
vnd_ft <- tabla_vnd %>%
  transmute(
    Gen = as.character(Protein), Variante = Mutation,
    col_muestras, col_pacientes, col_prev, col_cob, dom_hot
  )
names(vnd_ft) <- c(
  "Gen", "Variante", "Muestras, n (%)", "Pacientes, n (%)",
  "Frecuencia al\u00e9lica\n(media \u00b1 DE)",
  "Cobertura de posici\u00f3n\nMediana (RIC)",
  "Dominio / Hotspot"
)

# ── Bloques por gen (para sombreado y separadores) ───────────────────────────
n_ul97 <- sum(tabla_vnd$Protein == "UL97")
n_ul54 <- sum(tabla_vnd$Protein == "UL54")
n_ul56 <- sum(tabla_vnd$Protein == "UL56")
filas_ul97 <- seq_len(n_ul97)
filas_ul54 <- n_ul97 + seq_len(n_ul54)
filas_ul56 <- n_ul97 + n_ul54 + seq_len(n_ul56)
filas_sep  <- c(n_ul97, n_ul97 + n_ul54)

b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_vnd <- flextable::flextable(vnd_ft) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 7), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul54, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 1.6, unit = "cm") %>%
  flextable::width(j = 2, width = 1.9, unit = "cm") %>%
  flextable::width(j = 3, width = 2.3, unit = "cm") %>%
  flextable::width(j = 4, width = 2.3, unit = "cm") %>%
  flextable::width(j = 5, width = 2.7, unit = "cm") %>%
  flextable::width(j = 6, width = 3.3, unit = "cm") %>%
  flextable::width(j = 7, width = 4.6, unit = "cm") %>%
  flextable::padding(padding = 4, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND: variante no descrita en CHARMD (Tilloy et al., 2024), umbral \u22655%. ",
    "Porcentajes sobre muestras con secuenciaci\u00f3n v\u00e1lida por gen (UL97: 26; UL54: 24; UL56: 11) ",
    "y sobre pacientes con dato por gen (UL97: 11; UL54: 11; UL56: 3). ",
    "Frecuencia al\u00e9lica: media \u00b1 DE; en variantes con una sola muestra se indica el valor \u00fanico (n=1). ",
    "Cobertura: lecturas FW+RV en la posici\u00f3n, mediana (RIC). ",
    "Dominio / Hotspot: detalle solo en variantes localizadas en una regi\u00f3n de inter\u00e9s; en blanco si quedan fuera. ",
    "El asterisco (*) indica cod\u00f3n de parada."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_vnd



## ---- [6] Figura 47. Variantes unicas segun el umbral de frecuencia ----
# Apartado 6.2.5. Recuento por gen y clase a los umbrales 5, 10, 15 y 20%.

library(dplyr)
library(ggplot2)
library(tibble)

# ── Paleta ──────────────────────────────────────────────────────────────────
pal_genes <- c("UL97" = "#4A90C4", "UL54" = "#E8A33D", "UL56" = "#D9534F")

# ── Datos ───────────────────────────────────────────────────────────────────
datos_umbral <- tribble(
  ~Gen,   ~Umbral,  ~Tipo,  ~Variantes,
  "UL97", "≥ 5%",   "MRA",   3,
  "UL97", "≥ 10%",  "MRA",   3,
  "UL97", "≥ 15%",  "MRA",   2,
  "UL97", "≥ 20%",  "MRA",   2,
  "UL54", "≥ 5%",   "MRA",   1,
  "UL54", "≥ 10%",  "MRA",   1,
  "UL54", "≥ 15%",  "MRA",   1,
  "UL54", "≥ 20%",  "MRA",   1,
  "UL56", "≥ 5%",   "MRA",   0,
  "UL56", "≥ 10%",  "MRA",   0,
  "UL56", "≥ 15%",  "MRA",   0,
  "UL56", "≥ 20%",  "MRA",   0,
  "UL97", "≥ 5%",   "PN",    9,      # incluye D605E (5,7%)
  "UL97", "≥ 10%",  "PN",    8,
  "UL97", "≥ 15%",  "PN",    8,
  "UL97", "≥ 20%",  "PN",    8,
  "UL54", "≥ 5%",   "PN",   12,
  "UL54", "≥ 10%",  "PN",   11,
  "UL54", "≥ 15%",  "PN",   11,
  "UL54", "≥ 20%",  "PN",   11,
  "UL56", "≥ 5%",   "PN",    3,
  "UL56", "≥ 10%",  "PN",    3,
  "UL56", "≥ 15%",  "PN",    3,
  "UL56", "≥ 20%",  "PN",    3,
  "UL97", "≥ 5%",   "VND",  25,
  "UL97", "≥ 10%",  "VND",   2,
  "UL97", "≥ 15%",  "VND",   2,
  "UL97", "≥ 20%",  "VND",   1,
  "UL54", "≥ 5%",   "VND",  67,
  "UL54", "≥ 10%",  "VND",  17,
  "UL54", "≥ 15%",  "VND",  13,
  "UL54", "≥ 20%",  "VND",   6,
  "UL56", "≥ 5%",   "VND",   9,
  "UL56", "≥ 10%",  "VND",   1,
  "UL56", "≥ 15%",  "VND",   1,
  "UL56", "≥ 20%",  "VND",   1
) %>%
  mutate(
    Umbral = factor(Umbral, levels = c("≥ 20%", "≥ 15%", "≥ 10%", "≥ 5%")),
    Tipo   = factor(Tipo,   levels = c("MRA", "PN", "VND")),
    Gen    = factor(Gen,    levels = c("UL97", "UL54", "UL56"))
  )

# ── Figura ──────────────────────────────────────────────────────────────────
p_umbral <- ggplot(datos_umbral, aes(x = Umbral, y = Variantes, fill = Gen)) +
  geom_col(
    position  = position_dodge(width = 0.8),
    width     = 0.7,
    color     = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = if_else(Variantes > 0, as.character(Variantes), ""),
        color = Gen),
    position = position_dodge(width = 0.8),
    vjust    = -0.6,
    size     = 3.9,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~ Tipo, scales = "free_y", nrow = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = pal_genes) +
  scale_color_manual(values = pal_genes) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(
    x    = "Umbral de frecuencia de variante",
    y    = "Número de variantes únicas",
    fill = "Gen"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12),
    axis.title.x       = element_text(face = "bold", size = 11,
                                      margin = margin(t = 8)),
    axis.title.y       = element_text(face = "bold", size = 11,
                                      margin = margin(r = 8)),
    axis.text          = element_text(color = "black", size = 11),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.spacing      = unit(0.7, "cm"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 11),
    legend.margin      = margin(t = 0),
    plot.margin        = margin(4, 6, 4, 4)
  )

p_umbral

ggsave(
  filename = "output/figura_47_variantes_umbral_NGS.png",
  plot     = p_umbral,
  width    = 16,
  height   = 8,
  units    = "cm",
  dpi      = 600
)



## ---- [7a] Figura 48. Evolucion temporal de las MRA en RTP-2 ----
# Apartado 6.2.6. Carga viral, tratamiento antiviral y MRA a lo largo de 455 dias.

library(tidyverse)
library(grid)
library(gridExtra)
# ── Construcción de mra_rtp2 ─────────────────────────────────────────────────
fecha_txp_rtp2 <- as.Date("2022-06-17")

mra_rtp2 <- all_data_txp %>%
  filter(NHC == NHC_RTP02, class == "MRA") %>%
  mutate(
    Freg = as.Date(Freg),
    dia  = as.numeric(Freg - fecha_txp_rtp2)
  ) %>%
  filter(!is.na(dia), Prevalence >= 5) %>%
  mutate(
    Prevalence_group = factor(
      case_when(
        Prevalence >= 20 ~ "\u226520%",
        Prevalence >= 15 ~ "\u226515%",
        Prevalence >= 10 ~ "\u226510%",
        TRUE             ~ "\u22655%"
      ),
      levels = c("\u22655%", "\u226510%", "\u226515%", "\u226520%")
    ),
    Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  )

cat("mra_rtp2 construido:", nrow(mra_rtp2), "filas |",
    n_distinct(mra_rtp2$Mutation), "mutaciones únicas\n")
cat("\nDistribución por mutación y día:\n")
mra_rtp2 %>%
  select(dia, Protein, Mutation, Prevalence, Prevalence_group) %>%
  arrange(Protein, Mutation, dia) %>%
  print(n = Inf)


# ══════════════════════════════════════════════════════════════════════════════
#  COLORES
# ══════════════════════════════════════════════════════════════════════════════
col_GCV   <- "#3D7A5A"
col_VGCV  <- "#2E86AB"
col_MBV   <- "#B8860B"
col_LEF   <- "#C0652B"
col_FOS   <- "#9B2335"
col_LMV   <- "#7B5EA7"
col_CMVIg <- "#A8C8E0"

colores_fill  <- c("UL97" = alpha("#4472C4", 0.80),
                   "UL54" = alpha("#E6A817", 0.80),
                   "UL56" = alpha("#E8735A", 0.80))
colores_color <- c("UL97" = "#4472C4",
                   "UL54" = "#C8880A",
                   "UL56" = "#C94A30")
tamaños <- c("\u22655%" = 2, "\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)

# ══════════════════════════════════════════════════════════════════════════════
#  EJE TEMPORAL — índice numérico 1..51
# ══════════════════════════════════════════════════════════════════════════════
niveles_dia <- c(
  "1","26","42","49","55","62","69","75","84","94","108","122",
  "159","181","187","201","214","227","244","271","279","291","304",
  "361","371","382","388","390","395","402","409","418","423","433",
  "445","458","466","479","497","531","571","630","668","697","709",
  "742","774","801","850","872","909"
)
n_niv   <- length(niveles_dia)   # 51

pos_map        <- tibble(dia = as.numeric(niveles_dia), pos = seq_along(niveles_dia))
interpolar_pos <- function(d) approx(pos_map$dia, pos_map$pos, xout = d, rule = 2)$y

pos_evr <- interpolar_pos(279)

# ══════════════════════════════════════════════════════════════════════════════
#  DATOS CLÍNICOS RTP-2
# ══════════════════════════════════════════════════════════════════════════════
df_rtp2 <- tibble(
  dia = c(
    1,   26,  42,  49,  55,  62,  69,  75,  84,  94,  108, 122,
    159, 181, 187, 201, 214, 227, 244, 271, 279, 291, 304,
    361, 371, 382, 388, 390, 395, 402, 409, 418, 423, 433,
    445, 458, 466, 479, 497, 531, 571, 630, 668, 697, 709,
    742, 774, 801, 850, 872, 909
  ),
  log_cv = c(
    0.000, 0.000,
    4.278, 4.842, 5.042, 4.389, 3.891, 3.461, 3.234, 2.876, 0.000, 2.241,
    2.740, 3.557, 3.389, 3.461, 3.042, 2.885, 3.346,
    2.699, 3.279, 3.653, 4.008,
    3.741, 4.276, 4.610, 4.521, 4.422, 4.157, 4.312,
    0.000, 3.828, 3.872, 3.638,
    3.752, 4.381, 4.362, 4.364, 4.173,
    4.203, 3.712, 3.328, 2.965, 2.845, 2.823, 2.823,
    3.215, 3.034, 2.645, 3.002, 0.000
  ),
  muestra = c(
    NA, NA,
    1, 2, 3, 4, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, 5,
    6, 7, 8, 9, 10, 11, NA,
    NA, NA, NA, NA,
    NA, 12, 13, 14, 15,
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA
  ),
  secuenciada = c(
    F, F,
    T, T, T, T, F, F, F, F, F, F,
    F, F, F, F, F, F, F,
    F, F, F, T,
    T, T, T, T, T, T, F,
    F, F, F, F,
    F, T, T, T, T,
    F, F, F, F, F, F, F,
    F, F, F, F, F
  )
) %>%
  mutate(
    x_pos      = interpolar_pos(dia),
    tipo_punto = ifelse(secuenciada, "Muestra secuenciada", "Muestra no secuenciada")
  )

# ══════════════════════════════════════════════════════════════════════════════
#  DATOS DE MUTACIONES
# ══════════════════════════════════════════════════════════════════════════════
dia_remap <- c(
  "41"="42","48"="49","54"="55","61"="62",
  "303"="304","360"="361","370"="371",
  "381"="382","387"="388","389"="390","394"="395",
  "457"="402","465"="466","478"="479","496"="497"
)
remap_dia <- function(d) {
  s <- as.character(d)
  ifelse(s %in% names(dia_remap), dia_remap[s], s)
}

vsi_extra <- tibble(
  Mutation = c(
    "F396L","D605E",
    rep("A497P", 12),
    "W424*","G340C","D481N","Y281H","M697I"
  ),
  dia = c(
    54, 465,
    48,54,303,360,370,381,387,389,394,465,478,496,
    61,478,478,478,478
  ),
  Prevalence = c(
    12.9, 5.7,
    5.13,5.24,6.39,5.51,6.89,5.77,7.71,6.0,5.42,5.86,5.16,NA,
    23.9,6.2,5.1,8.0,8.4
  ),
  Protein = c(
    "UL54-VSI","UL97-VSI",
    rep("UL97-VND",12),
    "UL97-VND","UL97-VND","UL97-VND","UL56-VND","UL56-VND"
  )
) %>%
  filter(!is.na(Prevalence)) %>%
  mutate(
    Prevalence_group = factor(case_when(
      Prevalence>=20~"\u226520%",Prevalence>=15~"\u226515%",
      Prevalence>=10~"\u226510%",TRUE~"\u22655%"),
      levels=c("\u22655%","\u226510%","\u226515%","\u226520%")),
    Protein  = factor(Protein,
                      levels=c("UL97","UL54","UL56",
                               "UL54-VSI","UL97-VSI",
                               "UL97-VND","UL54-VND","UL56-VND")),
    dia_str  = remap_dia(dia),
    x_pos    = interpolar_pos(as.numeric(dia_str))
  )

mra_rtp2_base <- mra_rtp2 %>%
  filter(!is.na(Mutation), !is.na(Prevalence_group)) %>%
  mutate(
    dia_str = remap_dia(dia),
    x_pos   = interpolar_pos(as.numeric(dia_str))
  )

mra_rtp2_all <- bind_rows(mra_rtp2_base, vsi_extra) %>%
  filter(as.character(Mutation) %in% c("A594V","T409M","H411Y","T503I","T700A"))

mra_rtp2_plot <- mra_rtp2_all %>%
  filter(!is.na(Mutation), !is.na(Prevalence_group)) %>%
  mutate(
    Mutation  = as.character(Mutation),
    gen_color = case_when(
      Protein %in% c("UL97","UL97-VSI","UL97-VND") ~ "UL97",
      Protein %in% c("UL54","UL54-VSI","UL54-VND") ~ "UL54",
      Protein %in% c("UL56","UL56-VND")             ~ "UL56"
    ),
    gen_color = factor(gen_color, levels=c("UL97","UL54","UL56")),
    clase_var = case_when(
      Protein %in% c("UL97","UL54","UL56")             ~ "MRA",
      Protein %in% c("UL97-VSI","UL54-VSI")            ~ "VSI",
      Protein %in% c("UL97-VND","UL54-VND","UL56-VND") ~ "VND"
    ),
    clase_var = factor(clase_var, levels=c("MRA","VSI","VND")),
    # Orden: primero por gen (UL97, UL54, UL56), luego por mutación
    gen_orden = case_when(
      gen_color=="UL97" ~ 1L,
      gen_color=="UL54" ~ 2L,
      gen_color=="UL56" ~ 3L
    )
  ) %>%
  arrange(gen_orden, Mutation) %>%
  mutate(Mutation = factor(Mutation, levels=unique(Mutation)))

# ── Posiciones Y para mutaciones (zona negativa) ─────────────────────────────
mut_levels <- levels(mra_rtp2_plot$Mutation)
y_mut_base <- -0.4
y_mut_step <- -0.55
mut_y_map  <- setNames(
  y_mut_base + (seq_along(mut_levels)-1) * y_mut_step,
  mut_levels
)
mra_rtp2_plot <- mra_rtp2_plot %>%
  mutate(y_pos = mut_y_map[as.character(Mutation)])

# Color por gen para etiquetas de mutación
mut_gen_map <- mra_rtp2_plot %>%
  select(Mutation, gen_color) %>%
  distinct() %>%
  mutate(Mutation = as.character(Mutation))

x_mut_label <- -0.3   # X nombres de mutación

# Cruces UL54 no analizado (muestras 6 y 14)
cruces_ul54_muts <- mra_rtp2_plot %>%
  filter(Protein %in% c("UL54","UL54-VSI","UL54-VND")) %>%
  pull(Mutation) %>% as.character() %>% unique()
cruces_m6_m14 <- bind_rows(
  tibble(x_pos=interpolar_pos(361), y_pos=mut_y_map[cruces_ul54_muts], Mutation=cruces_ul54_muts),
  tibble(x_pos=interpolar_pos(479), y_pos=mut_y_map[cruces_ul54_muts], Mutation=cruces_ul54_muts)
)

# Cruces UL56 para muestras 1–7 (UL56 no analizado en esas muestras)
ul56_muts <- mra_rtp2_plot %>%
  filter(gen_color=="UL56") %>%
  pull(Mutation) %>% as.character() %>% unique()
dias_m1_7 <- c(42, 49, 55, 62, 304, 361, 371)   # días remapeados de muestras 1–7
cruces_ul56_m1_7 <- bind_rows(lapply(dias_m1_7, function(d) {
  tibble(x_pos=interpolar_pos(d),
         y_pos=mut_y_map[ul56_muts],
         Mutation=ul56_muts)
}))

mra_data <- mra_rtp2_plot %>% filter(clase_var=="MRA")
vsi_data <- mra_rtp2_plot %>% filter(clase_var %in% c("VSI","VND"))

# Etiquetas eje Y mutaciones con color por gen
mut_tick_df <- tibble(
  label = mut_levels,
  y_pos = mut_y_map[mut_levels]
) %>%
  left_join(mut_gen_map, by=c("label"="Mutation")) %>%
  mutate(col_text = case_when(
    gen_color=="UL97" ~ "#4472C4",
    gen_color=="UL54" ~ "#C8880A",
    gen_color=="UL56" ~ "#C94A30"
  ))

# ══════════════════════════════════════════════════════════════════════════════
#  BARRAS DE ANTIVIRAL — más altas, más juntas, letra mayor
# ══════════════════════════════════════════════════════════════════════════════
y_inf_trat <- 6.28; y_sup_trat <- 6.58; y_top_trat <- 6.84; bar_h <- 0.24

trat_df <- tibble(
  label    = c("GCV",  "VGCV", "GCV",  "MBV",  "LEF",  "FOS",  "LMV",  "CMV-Ig"),
  xmin_d   = c( 42,     69,     271,    291,    360,    382,    409,    62),
  xmax_d   = c( 68,    270,     290,    381,    381,    408,    432,   909),
  y_row    = c(y_sup_trat,y_sup_trat,y_sup_trat,y_sup_trat,
               y_top_trat,y_sup_trat,y_sup_trat,y_inf_trat),
  fill_col = c(col_GCV,col_VGCV,col_GCV,col_MBV,col_LEF,col_FOS,col_LMV,col_CMVIg),
  txt_col  = c("white","white","white","white","white","white","white","#1A3A52")
) %>%
  mutate(
    xmin_pos = interpolar_pos(xmin_d),
    xmax_pos = interpolar_pos(xmax_d),
    x_mid    = (xmin_pos + xmax_pos) / 2,
    ancho    = xmax_pos - xmin_pos,
    txt_size = case_when(ancho>3.5~3.4, ancho>1.8~3.0, TRUE~2.6)
  )

make_rounded_bar <- function(xmin, xmax, ymid, h, fill_col, id, n_arc=40) {
  r <- min(h/2, (xmax-xmin)/2)
  theta_l <- seq(pi/2, 3*pi/2, length.out=n_arc)
  theta_r <- seq(-pi/2, pi/2,  length.out=n_arc)
  tibble(
    x        = c(xmin+r+r*cos(theta_l), xmax-r+r*cos(theta_r)),
    y        = c(ymid+r*sin(theta_l),   ymid+r*sin(theta_r)),
    group    = id,
    fill_col = fill_col
  )
}

polys_trat <- pmap_dfr(
  list(trat_df$xmin_pos, trat_df$xmax_pos, trat_df$y_row,
       rep(bar_h, nrow(trat_df)), trat_df$fill_col, seq_len(nrow(trat_df))),
  function(xmin,xmax,ymid,h,fill_col,id) make_rounded_bar(xmin,xmax,ymid,h,fill_col,id)
)

# Límites verticales de la figura
y_lim_inf <- min(mut_y_map) - 0.4
y_lim_sup <- 7.25

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURA ÚNICA
# ══════════════════════════════════════════════════════════════════════════════
p <- ggplot() +
  
  # ── Fondo blanco zona tratamientos (encima del CV) ────────────────────────
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=6.05, ymax=y_lim_sup,
           fill="white", color=NA) +
  
  # ── Grid vertical solo entre Y=0 y Y=5.8 (no llega a zona antivirales) ───
  geom_vline(xintercept=seq_along(niveles_dia),
             color="#DDDDDD", linewidth=0.25,
             ymin=0, ymax=1) +   # ymin/ymax en unidades de npc — se recorta con annotate rect
  
  # ── Rect blanco que tapa el grid vertical en zona antivirales ────────────
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=5.85, ymax=y_lim_sup,
           fill="white", color=NA) +
  
  # ── Barras de antiviral ───────────────────────────────────────────────────
  geom_polygon(data=polys_trat,
               aes(x=x, y=y, group=group, fill=I(fill_col)),
               color=NA, inherit.aes=FALSE) +
  geom_text(data=trat_df %>% mutate(r=row_number()),
            aes(x=x_mid, y=y_row, label=label),
            color=trat_df$txt_col,
            size=trat_df$txt_size, fontface="bold",
            hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
  
  # ── Carga viral ───────────────────────────────────────────────────────────
  geom_line(data=df_rtp2, aes(x=x_pos, y=log_cv, group=1),
            linewidth=0.6, color="darkred") +
  geom_point(data=df_rtp2,
             aes(x=x_pos, y=log_cv, color=tipo_punto, shape=tipo_punto),
             size=2.5) +
  geom_text(data=df_rtp2 %>% filter(!is.na(muestra)),
            aes(x=x_pos, y=log_cv, label=muestra),
            vjust=-0.8, size=3.5, color="darkred", fontface="bold") +
  
  # ── Grid horizontal solo en zona mutaciones (desde etiqueta hasta fin) ───
  geom_segment(data=mut_tick_df,
               aes(x=x_mut_label, xend=n_niv+0.5, y=y_pos, yend=y_pos),
               color="#DDDDDD", linewidth=0.25, inherit.aes=FALSE) +
  
  # ── Etiquetas Y mutaciones coloreadas por gen ─────────────────────────────
  geom_text(data=mut_tick_df,
            aes(x=x_mut_label - 0.1, y=y_pos, label=label),
            color=mut_tick_df$col_text,
            hjust=1, size=3.8, fontface="bold",
            inherit.aes=FALSE) +
  
  # ── Etiqueta LogCV rotada, pegada al gráfico ──────────────────────────────
  annotate("text", x=-1.2, y=3.0,
           label="Log CV (UI/mL)", angle=90, size=4.0,
           color="grey30", fontface="bold", hjust=0.5, vjust=0) +
  
  # ── Etiquetas numéricas eje Y (Log CV) ──────────────────────────────────────
  annotate("text",
           x     = rep(0.55, 5),
           y     = 1:5,
           label = as.character(1:5),
           hjust = 1, size = 3.5,
           color = "grey35") +
  
  # ── Marcas de tick opcionales (muy discretas) ────────────────────────────────
  
  annotate("segment",
           x    = rep(0.75, 5), xend = rep(0.97, 5),
           y    = 1:5,          yend = 1:5,
           color = "grey70", linewidth = 0.25) +
  
  # ── Puntos mutaciones MRA (relleno sólido) ────────────────────────────────
  geom_point(data=mra_data,
             aes(x=x_pos, y=y_pos, size=Prevalence_group,
                 color=gen_color, fill=gen_color),
             shape=21, stroke=1, alpha=0.80) +
  # ── Puntos VSI + VND (relleno blanco) ────────────────────────────────────
  geom_point(data=vsi_data,
             aes(x=x_pos, y=y_pos, size=Prevalence_group, color=gen_color),
             shape=21, fill="white", stroke=1) +
  geom_point(data=cruces_m6_m14,
             aes(x=x_pos, y=y_pos),
             shape=4, size=2.4, stroke=0.8,
             color="grey35", inherit.aes=FALSE) +
  geom_point(data=cruces_ul56_m1_7,
             aes(x=x_pos, y=y_pos),
             shape=4, size=2.4, stroke=0.8,
             color="grey35", inherit.aes=FALSE) +
  
  # ── Escalas ───────────────────────────────────────────────────────────────
  scale_color_manual(
    values = c(colores_color,
               "Muestra secuenciada"    = "darkred",
               "Muestra no secuenciada" = "black"),
    breaks = c("UL97","UL54",
               "Muestra secuenciada","Muestra no secuenciada"),
    labels = c("UL97","UL54          ",
               "Muestra secuenciada","Muestra no secuenciada"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = colores_fill,
    breaks = c("UL97","UL54","UL56"),
    name   = NULL, guide="none"
  ) +
  scale_shape_manual(
    values = c("Muestra secuenciada"=16, "Muestra no secuenciada"=16),
    guide  = "none"
  ) +
  scale_size_manual(values=tamaños, name="Umbral") +
  
  scale_x_continuous(
    limits = c(-3.5, n_niv + 0.5),
    breaks = seq_along(niveles_dia),
    labels = niveles_dia,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(y_lim_inf, y_lim_sup),
    breaks = seq(0, 6, by=1),
    labels = NULL,
    expand = c(0,0)
  ) +
  
  guides(
    size = guide_legend(
      order=1, nrow=1,
      override.aes=list(shape=21, color=rep("grey50",4),
                        fill=rep("grey70",4), stroke=1)
    ),
    color = guide_legend(
      order=2, nrow=1,
      override.aes=list(
        shape  = c(21, 21, 16, 16),
        size   = c( 5,  5,  4,  4),
        alpha  = c(0.8,0.8, 1,  1),
        fill   = c(alpha("#4472C4",0.8), alpha("#E6A817",0.8),
                   "darkred", "black"),
        stroke = c(1, 1, 0, 0),
        color  = c("#4472C4","#C8880A","darkred","black")
      ),
      breaks=c("UL97","UL54",
               "Muestra secuenciada","Muestra no secuenciada")
    )
  ) +
  
  labs(x="D\u00edas postrasplante", y=NULL) +
  coord_cartesian(clip="off") +
  
  theme_minimal(base_size=13) +
  theme(
    axis.text.x       = element_text(angle=45, hjust=1, size=9),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid        = element_blank(),
    axis.title.x      = element_text(face="bold", margin=margin(t=12)),
    axis.title.y      = element_blank(),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.title      = element_text(size=11, face="bold"),
    legend.text       = element_text(size=11),
    legend.spacing.x  = unit(0.25,"cm"),
    plot.margin       = margin(5, 15, 5, 10),
    plot.background   = element_rect(fill="white", color=NA),
    panel.background  = element_rect(fill="white", color=NA)
  )

# Vista previa
grid::grid.newpage()
print(p)

# Guardado
png(filename="output/figura_RTP2_evolucion_MRA5.png",
    width=32, height=15, units="cm", res=300, bg="white")
print(p)
dev.off()



## ---- [7b] Figura 49. Evolucion temporal de MRA y VND en RTP-2 ----
# Apartado 6.2.6. Misma serie temporal incluyendo las variantes no descritas.

library(tidyverse)
library(grid)
library(gridExtra)

# ══════════════════════════════════════════════════════════════════════════════
#  COLORES
# ══════════════════════════════════════════════════════════════════════════════
col_GCV   <- "#3D7A5A"
col_VGCV  <- "#2E86AB"
col_MBV   <- "#B8860B"
col_LEF   <- "#C0652B"
col_FOS   <- "#9B2335"
col_LMV   <- "#7B5EA7"
col_CMVIg <- "#A8C8E0"

colores_fill  <- c("UL97" = alpha("#4472C4", 0.80),
                   "UL54" = alpha("#E6A817", 0.80),
                   "UL56" = alpha("#E8735A", 0.80))
colores_color <- c("UL97" = "#4472C4",
                   "UL54" = "#C8880A",
                   "UL56" = "#C94A30")
tamaños <- c("\u22655%" = 2, "\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)

# ══════════════════════════════════════════════════════════════════════════════
#  EJE TEMPORAL — índice numérico 1..51
# ══════════════════════════════════════════════════════════════════════════════
niveles_dia <- c(
  "1","26","42","49","55","62","69","75","84","94","108","122",
  "159","181","187","201","214","227","244","271","279","291","304",
  "361","371","382","388","390","395","402","409","418","423","433",
  "445","458","466","479","497","531","571","630","668","697","709",
  "742","774","801","850","872","909"
)
n_niv   <- length(niveles_dia)   # 51

pos_map        <- tibble(dia = as.numeric(niveles_dia), pos = seq_along(niveles_dia))
interpolar_pos <- function(d) approx(pos_map$dia, pos_map$pos, xout = d, rule = 2)$y

pos_evr <- interpolar_pos(279)


# ══════════════════════════════════════════════════════════════════════════════
#  DATOS CLÍNICOS RTP-2
# ══════════════════════════════════════════════════════════════════════════════
df_rtp2 <- tibble(
  dia = c(
    1,   26,  42,  49,  55,  62,  69,  75,  84,  94,  108, 122,
    159, 181, 187, 201, 214, 227, 244, 271, 279, 291, 304,
    361, 371, 382, 388, 390, 395, 402, 409, 418, 423, 433,
    445, 458, 466, 479, 497, 531, 571, 630, 668, 697, 709,
    742, 774, 801, 850, 872, 909
  ),
  log_cv = c(
    1.000, 1.000,
    4.278, 4.842, 5.042, 4.389, 3.891, 3.461, 3.234, 2.876, 1.000, 2.241,
    2.740, 3.557, 3.389, 3.461, 3.042, 2.885, 3.346,
    2.699, 3.279, 3.653, 4.008,
    3.741, 4.276, 4.610, 4.521, 4.422, 4.157, 4.312,
    1.000, 3.828, 3.872, 3.638,
    3.752, 4.381, 4.362, 4.364, 4.173,
    4.203, 3.712, 3.328, 2.965, 2.845, 2.823, 2.823,
    3.215, 3.034, 2.645, 3.002, 1.000
  ),
  muestra = c(
    NA, NA,
    1, 2, 3, 4, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, 5,
    6, 7, 8, 9, 10, 11, NA,
    NA, NA, NA, NA,
    NA, 12, 13, 14, 15,
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA
  ),
  secuenciada = c(
    F, F,
    T, T, T, T, F, F, F, F, F, F,
    F, F, F, F, F, F, F,
    F, F, F, T,
    T, T, T, T, T, T, F,
    F, F, F, F,
    F, T, T, T, T,
    F, F, F, F, F, F, F,
    F, F, F, F, F
  )
) %>%
  mutate(
    x_pos      = interpolar_pos(dia),
    tipo_punto = ifelse(secuenciada, "Muestra secuenciada", "Muestra no secuenciada")
  )

# ══════════════════════════════════════════════════════════════════════════════
#  DATOS DE MUTACIONES
# ══════════════════════════════════════════════════════════════════════════════
dia_remap <- c(
  "41"="42","48"="49","54"="55","61"="62",
  "303"="304","360"="361","370"="371",
  "381"="382","387"="388","389"="390","394"="395",
  "457"="402","465"="466","478"="479","496"="497"
)
remap_dia <- function(d) {
  s <- as.character(d)
  ifelse(s %in% names(dia_remap), dia_remap[s], s)
}

vsi_extra <- tibble(
  Mutation = c(
    # VSI
    "F396L","D605E",
    # VND UL97 — A497P
    rep("A497P", 12),
    # VND UL97 — otras
    "W424*","G340C","D481N",
    # VND UL56
    "Y281H","M697I","I508V",
    # VND UL54 — muestra 4 (dia 61→62)
    "F300L","Q323R","M491T","P704S","D711Y","P712H","D854G","S913G",
    # VND UL54 — muestra 5 (dia 303→304)
    "M491T","F702L","F817L"
  ),
  dia = c(
    54, 465,
    48,54,303,360,370,381,387,389,394,465,478,496,
    61,478,478,
    478,478,387,
    61,61,61,61,61,61,61,61,
    303,303,303
  ),
  Prevalence = c(
    12.9, 5.7,
    5.13,5.24,6.39,5.51,6.89,5.77,7.71,6.0,5.42,5.86,5.16,NA,
    23.9,6.2,5.1,
    8.0,8.4,5.2,
    7.9,25.0,10.2,7.3,8.6,8.2,19.6,15.1,
    6.9,5.2,5.2
  ),
  Protein = c(
    "UL54-VSI","UL97-VSI",
    rep("UL97-VND",12),
    "UL97-VND","UL97-VND","UL97-VND",
    "UL56-VND","UL56-VND","UL56-VND",
    rep("UL54-VND",8),
    rep("UL54-VND",3)
  )
) %>%
  filter(!is.na(Prevalence)) %>%
  mutate(
    Prevalence_group = factor(case_when(
      Prevalence>=20~"\u226520%",Prevalence>=15~"\u226515%",
      Prevalence>=10~"\u226510%",TRUE~"\u22655%"),
      levels=c("\u22655%","\u226510%","\u226515%","\u226520%")),
    Protein  = factor(Protein,
                      levels=c("UL97","UL54","UL56",
                               "UL54-VSI","UL97-VSI",
                               "UL97-VND","UL54-VND","UL56-VND")),
    dia_str  = remap_dia(dia),
    x_pos    = interpolar_pos(as.numeric(dia_str))
  )

mra_rtp2_base <- mra_rtp2 %>%
  filter(!is.na(Mutation), !is.na(Prevalence_group)) %>%
  mutate(
    dia_str = remap_dia(dia),
    x_pos   = interpolar_pos(as.numeric(dia_str))
  )

mra_rtp2_all <- bind_rows(mra_rtp2_base, vsi_extra)

mra_rtp2_plot <- mra_rtp2_all %>%
  filter(!is.na(Mutation), !is.na(Prevalence_group)) %>%
  mutate(
    Mutation  = as.character(Mutation),
    gen_color = case_when(
      Protein %in% c("UL97","UL97-VSI","UL97-VND") ~ "UL97",
      Protein %in% c("UL54","UL54-VSI","UL54-VND") ~ "UL54",
      Protein %in% c("UL56","UL56-VND")             ~ "UL56"
    ),
    gen_color = factor(gen_color, levels=c("UL97","UL54","UL56")),
    clase_var = case_when(
      Protein %in% c("UL97","UL54","UL56")             ~ "MRA",
      Protein %in% c("UL97-VSI","UL54-VSI")            ~ "VSI",
      Protein %in% c("UL97-VND","UL54-VND","UL56-VND") ~ "VND"
    ),
    clase_var = factor(clase_var, levels=c("MRA","VSI","VND")),
    # Orden: primero por gen (UL97, UL54, UL56), luego por mutación
    gen_orden = case_when(
      gen_color=="UL97" ~ 1L,
      gen_color=="UL54" ~ 2L,
      gen_color=="UL56" ~ 3L
    )
  ) %>%
  arrange(gen_orden, Mutation) %>%
  mutate(Mutation = factor(Mutation, levels=unique(Mutation)))

# ── Posiciones Y para mutaciones (zona negativa) ─────────────────────────────
mut_levels <- levels(mra_rtp2_plot$Mutation)
y_mut_base <- -0.1
y_mut_step <- -0.55
mut_y_map  <- setNames(
  y_mut_base + (seq_along(mut_levels)-1) * y_mut_step,
  mut_levels
)
mra_rtp2_plot <- mra_rtp2_plot %>%
  mutate(y_pos = mut_y_map[as.character(Mutation)])

# Color por gen para etiquetas de mutación
mut_gen_map <- mra_rtp2_plot %>%
  select(Mutation, gen_color) %>%
  distinct() %>%
  mutate(Mutation = as.character(Mutation))

x_mut_label <- -0.3   # X nombres de mutación

# Cruces UL54 no analizado (muestras 6 y 14)
cruces_ul54_muts <- mra_rtp2_plot %>%
  filter(Protein %in% c("UL54","UL54-VSI","UL54-VND")) %>%
  pull(Mutation) %>% as.character() %>% unique()
cruces_m6_m14 <- bind_rows(
  tibble(x_pos=interpolar_pos(361), y_pos=mut_y_map[cruces_ul54_muts], Mutation=cruces_ul54_muts),
  tibble(x_pos=interpolar_pos(479), y_pos=mut_y_map[cruces_ul54_muts], Mutation=cruces_ul54_muts)
)

# Cruces UL56 para muestras 1–7 (UL56 no analizado en esas muestras)
ul56_muts <- mra_rtp2_plot %>%
  filter(gen_color=="UL56") %>%
  pull(Mutation) %>% as.character() %>% unique()
dias_m1_7 <- c(42, 49, 55, 62, 304, 361, 371)   # días remapeados de muestras 1–7
cruces_ul56_m1_7 <- bind_rows(lapply(dias_m1_7, function(d) {
  tibble(x_pos=interpolar_pos(d),
         y_pos=mut_y_map[ul56_muts],
         Mutation=ul56_muts)
}))

mra_data <- mra_rtp2_plot %>% filter(clase_var=="MRA")
vsi_data <- mra_rtp2_plot %>% filter(clase_var %in% c("VSI","VND"))

# Etiquetas eje Y mutaciones con color por gen
mut_tick_df <- tibble(
  label = mut_levels,
  y_pos = mut_y_map[mut_levels]
) %>%
  left_join(mut_gen_map, by=c("label"="Mutation")) %>%
  mutate(col_text = case_when(
    gen_color=="UL97" ~ "#4472C4",
    gen_color=="UL54" ~ "#C8880A",
    gen_color=="UL56" ~ "#C94A30"
  ))

# ══════════════════════════════════════════════════════════════════════════════
#  BARRAS DE ANTIVIRAL — más altas, más juntas, letra mayor
# ══════════════════════════════════════════════════════════════════════════════
y_inf_trat <- 5.60; y_sup_trat <- 6.00; y_top_trat <- 6.38; bar_h <- 0.35

trat_df <- tibble(
  label    = c("GCV",  "VGCV", "GCV",  "MBV",  "LEF",  "FOS",  "LMV",  "CMV-Ig"),
  xmin_d   = c( 42,     69,     271,    291,    361,    382,    409,    62),
  xmax_d   = c( 68,    270,     290,    381,    381,    408,    432,   909),
  y_row    = c(y_sup_trat,y_sup_trat,y_sup_trat,y_sup_trat,
               y_top_trat,y_sup_trat,y_sup_trat,y_inf_trat),
  fill_col = c(col_GCV,col_VGCV,col_GCV,col_MBV,col_LEF,col_FOS,col_LMV,col_CMVIg),
  txt_col  = c("white","white","white","white","white","white","white","#1A3A52")
) %>%
  mutate(
    xmin_pos = interpolar_pos(xmin_d),
    xmax_pos = interpolar_pos(xmax_d),
    x_mid    = (xmin_pos + xmax_pos) / 2,
    ancho    = xmax_pos - xmin_pos,
    txt_size = case_when(ancho>3.5~3.4, ancho>1.8~3.0, TRUE~2.6)
  )

make_rounded_bar <- function(xmin, xmax, ymid, h, fill_col, id, n_arc=40) {
  r <- min(h/2, (xmax-xmin)/2)
  theta_l <- seq(pi/2, 3*pi/2, length.out=n_arc)
  theta_r <- seq(-pi/2, pi/2,  length.out=n_arc)
  tibble(
    x        = c(xmin+r+r*cos(theta_l), xmax-r+r*cos(theta_r)),
    y        = c(ymid+r*sin(theta_l),   ymid+r*sin(theta_r)),
    group    = id,
    fill_col = fill_col
  )
}

polys_trat <- pmap_dfr(
  list(trat_df$xmin_pos, trat_df$xmax_pos, trat_df$y_row,
       rep(bar_h, nrow(trat_df)), trat_df$fill_col, seq_len(nrow(trat_df))),
  function(xmin,xmax,ymid,h,fill_col,id) make_rounded_bar(xmin,xmax,ymid,h,fill_col,id)
)

# Límites verticales de la figura
y_lim_inf <- min(mut_y_map) - 0.1
y_lim_sup <- 7

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURA ÚNICA
# ══════════════════════════════════════════════════════════════════════════════
p <- ggplot() +
  
  # ── Fondo blanco zona tratamientos (encima del CV) ────────────────────────
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=5.50, ymax=y_lim_sup,
           fill="white", color=NA) +
  
  # ── Grid vertical solo entre Y=0 y Y=5.8 (no llega a zona antivirales) ───
  geom_vline(xintercept=seq_along(niveles_dia),
             color="#DDDDDD", linewidth=0.25,
             ymin=0, ymax=1) +
  
  # ── Rect blanco que tapa el grid vertical en zona antivirales ────────────
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=5.35, ymax=y_lim_sup,
           fill="white", color=NA) +
  
  # ── Barras de antiviral ───────────────────────────────────────────────────
  geom_polygon(data=polys_trat,
               aes(x=x, y=y, group=group, fill=I(fill_col)),
               color=NA, inherit.aes=FALSE) +
  geom_text(data=trat_df %>% mutate(r=row_number()),
            aes(x=x_mid, y=y_row, label=label),
            color=trat_df$txt_col,
            size=trat_df$txt_size, fontface="bold",
            hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
  
  # ── Etiquetas numéricas eje Y (Log CV) ───────────────────────────────────
  annotate("text",
           x     = rep(0.55, 5),
           y     = 1:5,
           label = as.character(1:5),
           hjust = 1, size = 3.5,
           color = "grey35") +
  
  # ── Marcas de tick ────────────────────────────────────────────────────────
  annotate("segment",
           x    = rep(0.75, 5), xend = rep(0.97, 5),
           y    = 1:5,          yend = 1:5,
           color = "grey70", linewidth = 0.25) +
  
  # ── Carga viral ───────────────────────────────────────────────────────────
  geom_line(data=df_rtp2, aes(x=x_pos, y=log_cv, group=1),
            linewidth=0.6, color="darkred") +
  geom_point(data=df_rtp2,
             aes(x=x_pos, y=log_cv, color=tipo_punto, shape=tipo_punto),
             size=2.5) +
  geom_text(data=df_rtp2 %>% filter(!is.na(muestra)),
            aes(x=x_pos, y=log_cv, label=muestra),
            vjust=-0.8, size=3.5, color="darkred", fontface="bold") +
  
  # ── Grid horizontal solo en zona mutaciones (desde etiqueta hasta fin) ───
  geom_segment(data=mut_tick_df,
               aes(x=x_mut_label, xend=n_niv+0.5, y=y_pos, yend=y_pos),
               color="#DDDDDD", linewidth=0.25, inherit.aes=FALSE) +
  
  # ── Etiquetas Y mutaciones coloreadas por gen ─────────────────────────────
  geom_text(data=mut_tick_df,
            aes(x=x_mut_label - 0.1, y=y_pos, label=label),
            color=mut_tick_df$col_text,
            hjust=1, size=3.8, fontface="bold",
            inherit.aes=FALSE) +
  
  # ── Etiqueta LogCV rotada, pegada al gráfico ──────────────────────────────
  annotate("text", x=-1.2, y=3.0,
           label="Log CV (UI/mL)", angle=90, size=4.0,
           color="grey30", fontface="bold", hjust=0.5, vjust=0) +
  
  
  
  # ── Puntos mutaciones MRA (relleno sólido) ────────────────────────────────
  geom_point(data=mra_data,
             aes(x=x_pos, y=y_pos, size=Prevalence_group,
                 color=gen_color, fill=gen_color),
             shape=21, stroke=1, alpha=0.80) +
  # ── Puntos VSI + VND (relleno blanco) ────────────────────────────────────
  geom_point(data=vsi_data,
             aes(x=x_pos, y=y_pos, size=Prevalence_group, color=gen_color),
             shape=21, fill="white", stroke=1) +
  geom_point(data=cruces_m6_m14,
             aes(x=x_pos, y=y_pos),
             shape=4, size=2.4, stroke=0.8,
             color="grey35", inherit.aes=FALSE) +
  geom_point(data=cruces_ul56_m1_7,
             aes(x=x_pos, y=y_pos),
             shape=4, size=2.4, stroke=0.8,
             color="grey35", inherit.aes=FALSE) +
  
  # ── Escalas ───────────────────────────────────────────────────────────────
  scale_color_manual(
    values = c(colores_color,
               "Muestra secuenciada"    = "darkred",
               "Muestra no secuenciada" = "black"),
    breaks = c("UL97","UL54","UL56",
               "Muestra secuenciada","Muestra no secuenciada"),
    labels = c("UL97","UL54","UL56          ",
               "Muestra secuenciada","Muestra no secuenciada"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = colores_fill,
    breaks = c("UL97","UL54","UL56"),
    name   = NULL, guide="none"
  ) +
  scale_shape_manual(
    values = c("Muestra secuenciada"=16, "Muestra no secuenciada"=16),
    guide  = "none"
  ) +
  scale_size_manual(values=tamaños, name="Umbral") +
  
  scale_x_continuous(
    limits = c(-3.5, n_niv + 0.5),
    breaks = seq_along(niveles_dia),
    labels = niveles_dia,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(y_lim_inf, y_lim_sup),
    breaks = seq(0, 6, by=1),
    labels = NULL,
    expand = c(0,0)
  ) +
  
  guides(
    size = guide_legend(
      order=1, nrow=1,
      override.aes=list(shape=21, color=rep("grey50",4),
                        fill=rep("grey70",4), stroke=1)
    ),
    color = guide_legend(
      order=2, nrow=1,
      override.aes=list(
        shape  = c(21, 21, 21, 16, 16),
        size   = c( 5,  5,  5,  4,  4),
        alpha  = c(0.8,0.8,0.8, 1,  1),
        fill   = c(alpha("#4472C4",0.8), alpha("#E6A817",0.8),
                   alpha("#E8735A",0.8), "darkred", "black"),
        stroke = c(1, 1, 1, 0, 0),
        color  = c("#4472C4","#C8880A","#C94A30","darkred","black")
      ),
      breaks=c("UL97","UL54","UL56",
               "Muestra secuenciada","Muestra no secuenciada")
    )
  ) +
  
  labs(x="D\u00edas postrasplante", y=NULL) +
  coord_cartesian(clip="off") +
  
  theme_minimal(base_size=13) +
  theme(
    axis.text.x       = element_text(angle=45, hjust=1, size=9),
    axis.ticks.x      = element_line(color="grey90", linewidth=0.3),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid        = element_blank(),
    axis.title.x      = element_text(face="bold", margin=margin(t=12)),
    axis.title.y      = element_blank(),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.title      = element_text(size=11, face="bold"),
    legend.text       = element_text(size=11),
    legend.spacing.x  = unit(0.25,"cm"),
    plot.margin       = margin(5, 15, 5, 10),
    plot.background   = element_rect(fill="white", color=NA),
    panel.background  = element_rect(fill="white", color=NA)
  )

# Vista previa
grid::grid.newpage()
print(p)

# Guardado
png(filename="output/figura_RTP2_evolucion_VSI.png",
    width=34, height=20, units="cm", res=300, bg="white")
print(p)
dev.off()
