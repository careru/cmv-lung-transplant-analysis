# CMV en receptores de trasplante pulmonar: código de análisis

Código en R utilizado para el análisis estadístico y la generación de figuras
y tablas de la tesis doctoral *Aspectos virológicos de la infección por
citomegalovirus en receptores de trasplante pulmonar* (Universidad Autónoma de
Madrid, 2026).

Cada sección del código está etiquetada con el apartado de la tesis al que
corresponde y con la figura o tabla que genera, de modo que cualquier resultado
publicado pueda trazarse hasta el código que lo produjo.

## Contenido

- `estudio1-genotipos/estudio1_genotipos.R` — Cohorte histórica cerrada de
  receptores de trasplante pulmonar trasplantados entre 2009 y 2021.
  Genotipado de las glucoproteínas gB (UL55) y gH (UL75) mediante PCR en
  tiempo real específica de genotipo. Análisis multinivel en muestra (n = 483),
  episodio de DNAemia (n = 234) y paciente (n = 146).

- `estudio2-ngs-resistencia/estudio2_ngs_resistencia.R` — Caracterización de
  variantes de resistencia antiviral en UL97, UL54 y UL56 mediante
  secuenciación masiva (DeepChek CMV v2.0) en 11 receptores con DNAemia
  ≥ 10.000 UI/mL entre 2022 y 2024. Clasificación de variantes según el marco
  CHARMD.

La correspondencia completa entre secciones de código, apartados de la tesis y
figuras o tablas se detalla en la cabecera de cada script.

## Disponibilidad de los datos

**Este repositorio contiene únicamente código.** No incluye datos de pacientes.

Los identificadores de historia clínica se han sustituido por constantes
simbólicas y las rutas locales por rutas relativas. El acceso a los datos puede
solicitarse a la autora, sujeto a la aprobación institucional y del comité de
ética correspondiente.

### Estructura esperada de los archivos de entrada

Los scripts esperan encontrar los datos en un directorio `data/` que no forma
parte de este repositorio.

**Estudio 1**

`data/genotipos_muestras.xls` — una fila por muestra procesada (n = 618), con
identificadores de muestra y paciente, fecha, carga viral, detección de cada
genotipo de gB y gH (0/1), valores de Ct, marca de repetición de la PCR y
observaciones.

`data/base_clinica.xlsx` — dos hojas. `pacientes_sel`, una fila por paciente
(n = 146) con fechas de trasplante, seguimiento y fallecimiento, perfil
serológico D/R, inmunosupresión basal y variables de desenlace clínico.
`dnaemias`, una fila por episodio (n = 234) con fechas de inicio y
aclaramiento, parámetros de cinética viral, sintomatología y genotipos
detectados.

`data/cargas_virales.dta` — todas las determinaciones de carga viral de la
cohorte (n = 12.022), con paciente, fecha y valor.

**Estudio 2**

`data/ngs_csv/` — informes CSV exportados de DeepChek CMV v2.0, uno por
muestra, con gen, posición, variante, frecuencia, cobertura por cadena y
calidad de llamada.

`data/muestras_clinicas.dta` — correspondencia entre muestra y paciente, con
carga viral y fecha.

`data/npet_excluidos.rds` — identificadores de muestra excluidos por control de
calidad.

## Requisitos

R (≥ 4.3) y los paquetes `tidyverse`, `survival`, `survminer`, `flextable`,
`officer`, `ggplot2`, `patchwork`, `ggpubr`, `rstatix`, `MASS`, `car`,
`readxl`, `haven` y `data.table`.

## Ejecución

Los scripts se organizan en secciones plegables de RStudio marcadas como
`## ---- [x] ----`. Cada sección es autocontenida una vez ejecutado el bloque
`[0]` de carga y depuración de datos.

## Cita

Reyes Ruiz CA. *Aspectos virológicos de la infección por citomegalovirus en
receptores de trasplante pulmonar* [tesis doctoral]. Madrid: Universidad
Autónoma de Madrid; 2026.

## Autoría

Carmen Alhena Reyes Ruiz — ORCID [0000-0003-0778-4971](https://orcid.org/0000-0003-0778-4971)

Tesis dirigida por la Dra. María Dolores Folgueira López, Servicio de
Microbiología, Hospital Universitario 12 de Octubre, Madrid.

## Licencia

Publicado bajo licencia MIT.
