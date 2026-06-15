#====================================================
#PACOTES
#===================================================

library(plm)
library(lmtest)
library(sandwich)
library(stargazer)
library(tidyverse)
library(ipeaplot)
library(strucchange)
library(broom)
library(patchwork)
library(purrr)
library(fixest)
library(plotly)
library(mediation)   # Mediação causal (Imai et al., 2010)

# Diagrama Causal (DAG)

library(dagitty)
library(ggdag)

dag <- dagitty("
dag {

desemprego -> homicidios
gini -> homicidios
homens_jovens -> homicidios
alcool -> homicidios
drogas -> homicidios
pib -> homicidios
densidade_pop. -> homicidios
urbanizacao -> homicidios
crescimento_pop. -> homicidios

armas -> homicidios

lei -> interacao
seguranca -> interacao
interacao -> armas
interacao -> homicidios
seguranca -> homicidios
seguranca -> armas

}
")
# Converter para tidy
dag_tidy <- tidy_dagitty(dag) |>
  mutate(
    label = case_when(
      name == "desemprego" ~ "Desemprego",
      name == "gini" ~ "Desigualdade",
      name == "homens_jovens" ~ "Homens jovens",
      name == "alcool" ~ "Álcool",
      name == "drogas" ~ "Drogas",
      name == "pib" ~ "PIB per capita",
      name == "homicidios" ~ "Mortes violentas",
      name == "armas" ~ "Armas de fogo",
      name == "seguranca" ~ "Gasto em segurança",
      name == "lei" ~ "Estatuto do Desarmamento",
      name == "interacao" ~ "Segurança × Lei",
      TRUE ~ name
    )
  )

tidy_dagitty(dag) |>
  ggdag(text = FALSE) +   # desliga texto padrão
  geom_dag_point(
    color = "red",
    fill = "lightgray",
    size = 27,
    shape = 21
  ) +
  geom_dag_text(
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  geom_dag_edges(
    edge_color = "blue"
  ) +
  theme_dag() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(color = "gray15", face = "bold"),
    plot.subtitle = element_text(color = "gray30")
  ) +
  labs(
    title = "DAG do mecanismo causal",
    subtitle = "Segurança pública, Estatuto do Desarmamento e homicídios"
  ) 

#========================================================================
# definição da variável preditiva principal e exploração básica
#=======================================================================


dados_completos <- dados_completos |>
  arrange(Estado, ano) |>
  group_by(Estado) |>
  mutate(
    despesa_seguranca = seg_pub_estad_imputado) |>
  mutate(
    lag_seg = lag(despesa_seguranca, 1),
    lag_l_seg = log(lag_seg + 1),
    lag_armas = lag(armas_de_fogo_cook,1),
    log_homicidios = log(tx_homicid_ajrand),
    log_mpaf = log(mortes_PAF),
  ) |>
  ungroup()

#Eploração básica
summary(plm(despesa_seguranca ~ plm::lag(tx_homicid_ajrand, 0:2) + 
              PIB_per_capita +
              lag_seg, 
            data = dados_completos, 
            model = "within", 
            effect = "twoway", 
            index = c("Estado", "ano")))

summary(plm(armas_de_fogo_cook ~ plm::lag(tx_homicid_ajrand, 0:2) +
              lag_armas, 
            data = dados_completos, 
            model = "within", 
            effect = "twoway", 
            index = c("Estado", "ano")))


summary(plm(tx_homicid_ajrand ~ plm::lag(tx_homicid_ajrand, 1:3) , 
            data = dados_completos, 
            model = "within", 
            effect = "twoway", 
            index = c("Estado", "ano")))

summary(plm(tx_homicid_ajrand ~ plm::lag(Gini, 0:1) +
              plm::lag(desemprego_imputado,0:1), 
            data = dados_completos, 
            model = "within", 
            effect = "twoway", 
            index = c("Estado", "ano")))

# ==============================================================================
# 4. CRIAR VARIÁVEIS PADRONIZADAS (Z-SCORE) - CÓDIGO CORRIGIDO
# ==============================================================================

dados_definidos <- dados_completos |> filter(ano >= 1985) |>
  mutate(
    s_Gini = as.numeric(scale(Gini, center = TRUE, scale = TRUE)),
    s_edu_estad = as.numeric(scale(gasto_edu_per_capita, center = TRUE, scale = TRUE)),
    s_desemprego = as.numeric(scale(desemprego_imputado, center = TRUE, scale = TRUE)),
    s_homens_jovens = as.numeric(scale(homens_jovens, center = TRUE, scale = TRUE)),
    s_uso_drogas = as.numeric(scale(uso_drogas, center = TRUE, scale = TRUE)),
    s_uso_alcool = as.numeric(scale(uso_alcool, center = TRUE, scale = TRUE)),
    s_urbano = as.numeric(scale(urbanizacao, center = TRUE, scale = TRUE)),
    s_armas_de_fogo = as.numeric(scale(armas_de_fogo_cook, center = TRUE, scale = TRUE)),
    s_lag_armas = as.numeric(scale(lag_armas, center = TRUE, scale = TRUE)),
    s_densidade = as.numeric(scale(densidade_populacional, center = TRUE, scale = TRUE)),
    log_PIB = as.numeric(log(PIB_per_capita)),
    log_homens_jovens = as.numeric(log(homens_jovens)),
    log_drogas = as.numeric(log1p(uso_drogas)),
    log_armas = as.numeric(log1p(armas_de_fogo_cook)),
    log_alcool = as.numeric(log1p(uso_alcool)),
    log_desemprego = as.numeric(log(desemprego_imputado)),
    log_Gini = as.numeric(log(Gini)),
    log_urbano = as.numeric(log(urbanizacao)),
    s_PIB_per_capita = as.numeric(scale(PIB_per_capita, center = TRUE, scale = TRUE)),
    s_despesa_seguranca = as.numeric(scale(despesa_seguranca, center = TRUE, scale = TRUE)),
    s_saude_estad = as.numeric(scale(gasto_saúde_per_capita, center = TRUE, scale = TRUE)),
    s_industria_PIB = as.numeric(scale(industria_PIB, center = TRUE, scale = TRUE)),
    s_cresc_PIB = as.numeric(scale(var_PIB, center = TRUE, scale = TRUE)),
    s_ano = as.numeric(scale(ano, center = TRUE, scale = TRUE)),
    lag_s_seg = as.numeric(scale(lag_seg, center = TRUE, scale = TRUE)),
    s_cresc_pop = as.numeric(scale(cresc_populacional, center = TRUE, scale = TRUE)),
    s_setor_publico = as.numeric(scale(setor_publico_per_capita, center = TRUE, scale = TRUE)),
    i_seg_desarm = as.numeric(scale(gasto_seguranca, center = T, scale = T) * Lei_desarmamento_10826_2003)
  )



# ==============================================================================
# Código unificado e padronizado para análise do Estatuto do Desarmamento
# Impacto do gasto em segurança sobre homicídios e mortes por armas de fogo
# ==============================================================================

# Diretórios de saída (ajuste conforme necessário)
dir_homicidios <- "/home/matheus/Documentos/IPEA/modelos/estatuto do desarmamento como estratégia de segurança pública/Homicídios"
dir_mortes     <- "/home/matheus/Documentos/IPEA/modelos/estatuto do desarmamento como estratégia de segurança pública/Mortes por armas de fogo"


#===============================================================================================
# Análise de interação entre gasto em segurança e Estatuto do Desarmamento
#===============================================================================================

# definir controles e preditores básicos
preditores_base <- "s_desemprego+s_Gini+s_homens_jovens+s_uso_drogas+s_uso_alcool+s_PIB_per_capita+s_urbano+s_densidade+s_cresc_pop"

#=================================================================================================
#Testes de Cointegração
#=================================================================================================

library(dcce)
# ---------------------------------------------------------
# Teste de Pedroni (Heterogeneidade total de coeficientes)
# ---------------------------------------------------------
teste_pedroni_hmcd <- panel_coint_test(
  data = dados_definidos,
  unit_index = "Estado",  # Identificador do painel (N)
  time_index = "ano",       # Identificador do tempo (T)
  formula = log_homicidios ~ s_desemprego+s_Gini+s_homens_jovens+s_uso_drogas+s_uso_alcool+s_PIB_per_capita+s_urbano+s_densidade+s_cresc_pop +
    s_armas_de_fogo +  
     + i_seg_desarm,
  test = "pedroni",
  lags = 1L                   # Defasagens para a regressão ADF nos resíduos
)

print(teste_pedroni_hmcd)
# Interpretação: Valores t e rho altamente negativos 
# rejeitam a hipótese nula de que NÃO há cointegração.

teste_pedroni_mpaf <- panel_coint_test(
  data = dados_definidos,
  unit_index = "Estado",  # Identificador do painel (N)
  time_index = "ano",       # Identificador do tempo (T)
  formula = log_mpaf ~ s_desemprego+s_Gini+s_homens_jovens+s_uso_drogas+s_uso_alcool+s_PIB_per_capita+s_urbano+s_densidade+s_cresc_pop +
    s_armas_de_fogo +  
    + i_seg_desarm,
  test = "pedroni",
  lags = 1L                   # Defasagens para a regressão ADF nos resíduos
)

print(teste_pedroni_mpaf)
# Interpretação: Valores t e rho altamente negativos 
# rejeitam a hipótese nula de que NÃO há cointegração.

# ---------------------------------------------------------
# Teste de Kao (Assume inclinações homogêneas no painel)
# ---------------------------------------------------------
teste_kao_hmcd <- panel_coint_test(
  data = dados_definidos,
  unit_index = "Estado",
  time_index = "ano",
  formula = log_homicidios ~ s_desemprego+s_Gini+s_homens_jovens+s_uso_drogas+s_uso_alcool+s_PIB_per_capita+s_urbano+s_densidade+s_cresc_pop +
    s_armas_de_fogo +  
    + i_seg_desarm,
  test = "kao",
  lags = 1L
)

print(teste_kao_hmcd)

# Interpretação: Valores t e rho altamente negativos 
# rejeitam a hipótese nula de que NÃO há cointegração.

teste_kao_mpaf <- panel_coint_test(
  data = dados_definidos,
  unit_index = "Estado",
  time_index = "ano",
  formula = log_mpaf ~ s_desemprego+s_Gini+s_homens_jovens+s_uso_drogas+s_uso_alcool+s_PIB_per_capita+s_urbano+s_densidade+s_cresc_pop +
    s_armas_de_fogo +  
    + i_seg_desarm,
  test = "kao",
  lags = 1L
)

print(teste_kao_mpaf)

# Interpretação: Valores t e rho altamente negativos 
# rejeitam a hipótese nula de que NÃO há cointegração.



# Efeitos marginais

# ----------------------------------------
# 1. Funções auxiliares
# ----------------------------------------

# 1.1 Estimar modelos de efeitos fixos para um dado desfecho
estimar_modelos <- function(dados, outcome) {
  # outcome: string com o nome da variável dependente (ex: "log(tx_homicid_ajrand)")
  
  # Extrair o nome da variável base para o lag (remove "log()" se existir)
  var_lag <- gsub("log\\(|\\)", "", outcome) 
  
  # ----------------------------------------
  # A. DEFINIÇÃO DA FÓRMULA BÁSICA ÚNICA
  # ----------------------------------------
  # Agrupamos todas as covariáveis que se repetem em todos os modelos.
  
  # Construímos as fórmulas concatenando a base com as variáveis específicas
  formula_base          <- as.formula(paste(outcome, "~", preditores_base, "+ s_armas_de_fogo + s_despesa_seguranca"))
  formula_base_sAF      <- as.formula(paste(outcome, "~", preditores_base, "+ s_despesa_seguranca"))

  formula_desarm        <- as.formula(paste(outcome, "~", preditores_base, "+ s_armas_de_fogo + s_despesa_seguranca * Lei_desarmamento_10826_2003"))
  formula_desarm_sAF    <- as.formula(paste(outcome, "~", preditores_base, "+ s_despesa_seguranca * Lei_desarmamento_10826_2003"))

  # ----------------------------------------
  # B. ESTIMATIVA DOS MODELOS
  # ----------------------------------------
  modelos <- list()
  modelos$base          <- plm(formula_base,          data = dados, model = "within", effect = "twoway", index = c("Estado", "ano"))
  modelos$base_sAF      <- plm(formula_base_sAF,      data = dados, model = "within", effect = "twoway", index = c("Estado", "ano"))
  modelos$desarm        <- plm(formula_desarm,        data = dados, model = "within", effect = "twoway", index = c("Estado", "ano"))
  modelos$desarm_sAF    <- plm(formula_desarm_sAF,    data = dados, model = "within", effect = "twoway", index = c("Estado", "ano"))

  # ----------------------------------------
  # C. CÁLCULO DOS ERROS (Driscoll-Kraay)
  # ----------------------------------------
  # *CORREÇÃO APLICADA: Referências corretas para os modelos _ar*
  erros <- list()
  erros$base          <- vcovSCC(modelos$base,          type = "sss", cluster = "group")
  erros$base_sAF      <- vcovSCC(modelos$base_sAF,      type = "sss", cluster = "group")
  erros$desarm        <- vcovSCC(modelos$desarm,        type = "sss", cluster = "group")
  erros$desarm_sAF    <- vcovSCC(modelos$desarm_sAF,    type = "sss", cluster = "group")

  list(modelos = modelos, erros = erros)
}

# ----------------------------------------
# Análise para HOMICÍDIOS
# ----------------------------------------

# 3.1 Modelos de efeitos fixos
modelos_hom <- estimar_modelos(dados_definidos, "log(tx_homicid_ajrand)")

# Tabela com stargazer
stargazer(modelos_hom$modelos, se = modelos_hom$erros,
          title = "Homicídios - Efeitos Fixos com Erros D-K",
          out = file.path(dir_homicidios, "homicidios_seguranca.html"),
                    type = "text", style = "jpam",
          column.labels = c("base", "s/AF", "ED", "EDs/AF"),
          model.numbers = FALSE, keep.stat = c("n", "rsq"),
          star.cutoffs = c(0.1, 0.05, 0.01))


# ----------------------------------------
# Análise para mortes por armas de fogo
# ----------------------------------------

# 4.1 Modelos de efeitos fixos
modelos_mort <- estimar_modelos(dados_definidos, "log(mortes_PAF)")

# Tabela
stargazer(modelos_mort$modelos, se = modelos_mort$erros,
          title = "Mortes por Armas - Efeitos Fixos com Erros D-K",
          type = "text", style = "jpam",
          out = file.path(dir_mortes, "mortes_seguranca.html"),
          column.labels = c("base", "s/AF", "ED", "EDs/AF"),
          model.numbers = FALSE, keep.stat = c("n", "rsq"),
          star.cutoffs = c(0.1, 0.05, 0.01))


#==============================================================================
# Efeitos Marginais
#==============================================================================

if(!require(marginaleffects)) install.packages("marginaleffects")
library(marginaleffects)

calc_efeito_marginal_dk <- function(coefs, vcov_dk, var_gasto, var_interacao) {
  efeito_0 <- coefs[var_gasto]
  se_0 <- sqrt(vcov_dk[var_gasto, var_gasto])
  efeito_1 <- coefs[var_gasto] + coefs[var_interacao]
  var_1 <- vcov_dk[var_gasto, var_gasto] + vcov_dk[var_interacao, var_interacao] + 2*vcov_dk[var_gasto, var_interacao]
  se_1 <- sqrt(var_1)
  data.frame(Periodo = c("Pré-Estatuto", "Pós-Estatuto"),
             estimate = c(efeito_0, efeito_1),
             conf.low = c(efeito_0 - 1.96*se_0, efeito_1 - 1.96*se_1),
             conf.high = c(efeito_0 + 1.96*se_0, efeito_1 + 1.96*se_1))
}

efeitos_hom_corrigido <- calc_efeito_marginal_dk(coef(modelos_hom$modelos$desarm),
                                                 modelos_hom$erros$desarm,
                                                 "s_despesa_seguranca",
                                                 "s_despesa_seguranca:Lei_desarmamento_10826_2003")

# Aplicar para homicídios
efeitos_hom <- calc_efeito_marginal_dk(
  coefs = coef(modelos_hom$modelos$desarm),
  vcov_dk = modelos_hom$erros$desarm,
  var_gasto = "s_despesa_seguranca",
  var_interacao = "s_despesa_seguranca:Lei_desarmamento_10826_2003"
)

efeitos_hom

# Plotar
marginais_homicidio <- ggplot(efeitos_hom, aes(x = Periodo, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(size = 1) +
  geom_text(aes(label = sprintf("%.3f", estimate)), vjust = -1.5, fontface = "bold") +
  labs(title = "Homicídios: Efeito marginal do gasto (com IC 95%)",
       x = "Período", y = "Efeito Estimado") +
  theme_ipea() +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

print(marginais_homicidio)

save_ipeaplot(marginais_homicidio, path = dir_homicidios,
              file.name = "efeitos marginais homicídios", 
              format = c("eps", "png"))

#Aplicar a mortes PAF

efeitos_mpaf <- calc_efeito_marginal_dk(
  coefs = coef(modelos_mort$modelos$desarm),
  vcov_dk = modelos_mort$erros$desarm,
  var_gasto = "s_despesa_seguranca",
  var_interacao = "s_despesa_seguranca:Lei_desarmamento_10826_2003"
)


efeitos_mpaf

marginais_mpaf <- ggplot(efeitos_mpaf, aes(x = Periodo, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(size = 1) +
  geom_text(aes(label = sprintf("%.3f", estimate)), vjust = -1.5, fontface = "bold") +
  labs(title = "Mortes por armas de fogo: Efeito marginal do gasto (com IC 95%)",
       x = "Período", y = "Efeito Estimado") +
  theme_ipea() +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

print(marginais_mpaf)

save_ipeaplot(marginais_mpaf, path = dir_homicidios,
              file.name = "efeitos marginais homicídios", 
              format = c("eps", "png"))


# =============================================================================
# ANÁLISE DE MEDIAÇÃO VIA EFEITOS MARGINAIS (USANDO MODELOS COM INTERAÇÃO)
# =============================================================================
# Baseia-se nos modelos de desfecho já estimados (desarm) e estima modelos
# para o mediador (armas) com a mesma estrutura de interação.
# -----------------------------------------------------------------------------

# ----------------------------------------
# 1. Função auxiliar para mediação a partir de modelos com interação
# ----------------------------------------

#' Mediação usando efeitos marginais de modelos com interação
#' @param modelo_m Modelo plm para o mediador (com interação tratamento * lei)
#' @param modelo_y Modelo plm para o desfecho (com interação tratamento * lei)
#' @param treat Nome da variável tratamento (ex: "s_despesa_seguranca")
#' @param mediator Nome da variável mediadora (ex: "s_armas_de_fogo")
#' @param lei Nome da variável indicadora da lei (ex: "Lei_desarmamento_10826_2003")
#' @param vcov_m Matriz de covariância para o modelo do mediador (opcional)
#' @param vcov_y Matriz de covariância para o modelo do desfecho (opcional)
#' @return Data.frame com ACME, ADE e Total para cada período (pré e pós)
#' 



mediacao_interacao <- function(modelo_m, modelo_y,
                               treat = "s_despesa_seguranca",
                               mediator = "s_armas_de_fogo",
                               lei = "Lei_desarmamento_10826_2003",
                               vcov_m = NULL, vcov_y = NULL) {
  
  # Extrair coeficientes e matrizes de variância-covariância
  coef_m <- coef(modelo_m)
  coef_y <- coef(modelo_y)
  
  if (is.null(vcov_m)) vcov_m <- vcov(modelo_m)
  if (is.null(vcov_y)) vcov_y <- vcov(modelo_y)
  
  # Nomes dos coeficientes de interação (assumindo formato "treat:lei")
  interacao <- paste0(treat, ":", lei)
  
  # -------------------------------------------------
  # Efeitos do tratamento no mediador (alpha)
  # -------------------------------------------------
  # alpha quando lei = 0
  alpha_0 <- coef_m[treat]
  # alpha quando lei = 1
  alpha_1 <- coef_m[treat] + coef_m[interacao]
  
  # Variâncias e covariâncias para alpha
  var_alpha_0 <- vcov_m[treat, treat]
  var_alpha_1 <- vcov_m[treat, treat] + vcov_m[interacao, interacao] + 2 * vcov_m[treat, interacao]
  se_alpha_0 <- sqrt(var_alpha_0)
  se_alpha_1 <- sqrt(var_alpha_1)
  
  # -------------------------------------------------
  # Efeito do mediador no desfecho (beta) - assumido constante
  # -------------------------------------------------
  beta <- coef_y[mediator]
  se_beta <- sqrt(vcov_y[mediator, mediator])
  
  # -------------------------------------------------
  # Efeito direto do tratamento no desfecho (ADE)
  # -------------------------------------------------
  ade_0 <- coef_y[treat]
  ade_1 <- coef_y[treat] + coef_y[interacao]
  
  var_ade_0 <- vcov_y[treat, treat]
  var_ade_1 <- vcov_y[treat, treat] + vcov_y[interacao, interacao] + 2 * vcov_y[treat, interacao]
  se_ade_0 <- sqrt(var_ade_0)
  se_ade_1 <- sqrt(var_ade_1)
  
  # -------------------------------------------------
  # Efeito indireto (ACME) = alpha * beta
  # -------------------------------------------------
  acme_0 <- alpha_0 * beta
  acme_1 <- alpha_1 * beta
  
  # Variância do ACME via delta (assumindo independência entre alpha e beta)
  var_acme_0 <- (beta^2 * var_alpha_0) + (alpha_0^2 * vcov_y[mediator, mediator])
  var_acme_1 <- (beta^2 * var_alpha_1) + (alpha_1^2 * vcov_y[mediator, mediator])
  se_acme_0 <- sqrt(var_acme_0)
  se_acme_1 <- sqrt(var_acme_1)
  
  # -------------------------------------------------
  # Efeito total = direto + indireto
  # -------------------------------------------------
  total_0 <- ade_0 + acme_0
  total_1 <- ade_1 + acme_1
  
  # Variância do total (assumindo independência entre ADE e ACME)
  var_total_0 <- var_ade_0 + var_acme_0
  var_total_1 <- var_ade_1 + var_acme_1
  se_total_0 <- sqrt(var_total_0)
  se_total_1 <- sqrt(var_total_1)
  
  # Intervalos de confiança (95%)
  ic_95 <- function(est, se) c(est - 1.96 * se, est + 1.96 * se)
  
  # Montar tabela de resultados
  resultados <- data.frame(
    Periodo = rep(c("Pré-Estatuto (Lei=0)", "Pós-Estatuto (Lei=1)"), each = 3),
    Efeito = rep(c("Indireto (ACME)", "Direto (ADE)", "Total"), times = 2),
    Estimativa = c(acme_0, ade_0, total_0, acme_1, ade_1, total_1),
    Erro_Padrao = c(se_acme_0, se_ade_0, se_total_0, se_acme_1, se_ade_1, se_total_1),
    IC_low = c(ic_95(acme_0, se_acme_0)[1], ic_95(ade_0, se_ade_0)[1], ic_95(total_0, se_total_0)[1],
               ic_95(acme_1, se_acme_1)[1], ic_95(ade_1, se_ade_1)[1], ic_95(total_1, se_total_1)[1]),
    IC_high = c(ic_95(acme_0, se_acme_0)[2], ic_95(ade_0, se_ade_0)[2], ic_95(total_0, se_total_0)[2],
                ic_95(acme_1, se_acme_1)[2], ic_95(ade_1, se_ade_1)[2], ic_95(total_1, se_total_1)[2])
  )
  
  return(resultados)
}

# ----------------------------------------
# 2. Estimar modelos do mediador (armas) com interação
# ----------------------------------------

# Fórmula para o mediador (mesmos controles, com interação)
formula_mediador <- as.formula(
  "s_armas_de_fogo ~ lag(s_desemprego, 1) + lag(s_Gini,1) + lag(s_homens_jovens,1) +
   lag(s_uso_drogas, 1) + lag(s_uso_alcool,1) + lag(s_PIB_per_capita,1) +
   s_despesa_seguranca * Lei_desarmamento_10826_2003"
)

# Modelo para o mediador (armas) com efeitos fixos twoway
modelo_m_armas <- plm(
  formula_mediador,
  data = dados_definidos,
  model = "within",
  effect = "twoway",
  index = c("Estado", "ano")
)

# Matriz de covariância Driscoll-Kraay para o mediador
vcov_m_armas <- vcovSCC(modelo_m_armas, type = "sss", cluster = "time")

# ----------------------------------------
# 3. Aplicar mediação para HOMICÍDIOS
# ----------------------------------------
cat("\n========== MEDIAÇÃO COM INTERAÇÃO - HOMICÍDIOS ==========\n")

resultados_hom_int <- mediacao_interacao(   # <--- corrigido
  modelo_m = modelo_m_armas,
  modelo_y = modelos_hom$modelos$desarm,
  treat = "s_despesa_seguranca",
  mediator = "s_armas_de_fogo",
  lei = "Lei_desarmamento_10826_2003",
  vcov_m = vcov_m_armas,
  vcov_y = modelos_hom$erros$desarm
)

print(resultados_hom_int)

# ----------------------------------------
# 4. Aplicar mediação para MORTES POR ARMAS DE FOGO
# ----------------------------------------
cat("\n========== MEDIAÇÃO COM INTERAÇÃO - MORTES PAF ==========\n")

resultados_mort_int <- mediacao_interacao(   # <--- corrigido
  modelo_m = modelo_m_armas,
  modelo_y = modelos_mort$modelos$desarm,
  treat = "s_despesa_seguranca",
  mediator = "s_armas_de_fogo",
  lei = "Lei_desarmamento_10826_2003",
  vcov_m = vcov_m_armas,
  vcov_y = modelos_mort$erros$desarm
)

print(resultados_mort_int)

# ----------------------------------------
# 5. Salvar resultados
# ----------------------------------------
write.csv(resultados_hom_int, file.path(dir_homicidios, "mediacao_interacao_homicidios.csv"), row.names = FALSE)
write.csv(resultados_mort_int, file.path(dir_mortes, "mediacao_interacao_mortes.csv"), row.names = FALSE)

#Visualização dos efeitos marginais e mediação

# ----------------------------------------
# 1. Preparar dados para plotagem
# ----------------------------------------

# Adicionar coluna de desfecho
resultados_hom_int$Desfecho <- "Homicídios"
resultados_mort_int$Desfecho <- "Mortes por Arma de Fogo"

# Combinar os dois data frames
resultados_mediacao <- bind_rows(resultados_hom_int, resultados_mort_int)

# Criar uma variável de período mais curta para os rótulos
resultados_mediacao$Periodo_short <- ifelse(
  grepl("Pré", resultados_mediacao$Periodo),
  "Pré-Estatuto",
  "Pós-Estatuto"
)

# ----------------------------------------
# 2. Gráfico de decomposição para cada desfecho
# ----------------------------------------

# Função para gerar o gráfico de um desfecho específico
plot_decomposicao <- function(dados, titulo) {
  ggplot(dados, aes(x = Efeito, y = Estimativa, fill = Periodo_short)) +
    geom_col(position = position_dodge(0.9), width = 0.7) +
    geom_errorbar(
      aes(ymin = IC_low, ymax = IC_high),
      position = position_dodge(0.9),
      width = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = titulo,
      x = "Tipo de Efeito",
      y = "Magnitude do Efeito",
      fill = "Período"
    ) +
    theme_ipea() +
    theme(legend.position = "bottom") +
    scale_fill_ipea(palette = "Red-Blue-White")
}


# Gráfico para homicídios
p_hom <- plot_decomposicao(
  filter(resultados_mediacao, Desfecho == "Homicídios"),
  "Homicídios: Decomposição dos Efeitos do Gasto em Segurança"
)

# Gráfico para mortes por arma de fogo
p_mort <- plot_decomposicao(
  filter(resultados_mediacao, Desfecho == "Mortes por Arma de Fogo"),
  "Mortes por Armas: Decomposição dos Efeitos do Gasto em Segurança"
)

# Exibir
print(p_hom)
print(p_mort)

plot_decomposicao_junto <- p_hom/p_mort

print(plot_decomposicao_junto)

# Salvar
save_ipeaplot(p_hom, path =  dir_homicidios,
              file.name = "decomposicao_homicidios", format = c("png", "eps"))

save_ipeaplot(p_mort, path = dir_homicidios,
              file.name = "decomposicao_mortes", format = c("png", "eps"))

save_ipeaplot(plot_decomposicao_junto, 
              file.name = "decomposicao_juntos", format = c("png", "eps"))


# ----------------------------------------
# 3. Gráfico comparativo: ACME e ADE entre desfechos
# ----------------------------------------

# Filtrar apenas ACME e ADE (excluir Total)
efeitos_principais <- resultados_mediacao %>%
  filter(Efeito %in% c("Indireto (ACME)", "Direto (ADE)", "Total"))

# Gráfico de pontos com intervalos
p_comp <- ggplot(efeitos_principais, 
                 aes(x = Desfecho, y = Estimativa, color = Periodo_short, shape = Efeito)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbar(aes(ymin = IC_low, ymax = IC_high),
                position = position_dodge(0.5), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Efeitos Direto e Indireto do Gasto em Segurança",
       x = "Desfecho",
       y = "Estimativa",
       color = "Período",
       shape = "Efeito") +
  theme_ipea() +  theme(legend.position = "bottom") +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

print(p_comp)

save_ipeaplot(p_comp, path = dir_homicidios, file.name = "comparacao_efeitos", format = c("eps", "png"))

# ----------------------------------------
# 4. Gráfico de barras empilhadas mostrando a proporção do efeito total que é mediada
# ----------------------------------------

# Calcular proporção mediada (ACME / Total) para cada período e desfecho
proporcao <- resultados_mediacao |>
  dplyr::filter(Efeito %in% c("Indireto (ACME)", "Total")) |>
  dplyr::select(Desfecho, Periodo_short, Efeito, Estimativa) |>
  tidyr::pivot_wider(names_from = Efeito, values_from = Estimativa) |>
  dplyr::mutate(Proporcao = `Indireto (ACME)` / Total)

# Gráfico de barras da proporção
p_prop <- ggplot(proporcao, aes(x = Desfecho, y = Proporcao, fill = Periodo_short)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Proporção do Efeito Total Mediada pelo Estoque de Armas",
       x = "Desfecho",
       y = "Proporção Mediada (ACME / Total)",
       fill = "Período") +
  theme_ipea() +
  theme(legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

print(p_prop)
save_ipeaplot(p_prop, path = dir_homicidios, file.name =  "proporcao_mediada", format = c("png", "eps"))

# ----------------------------------------
# 5. (Opcional) Visualização dos efeitos marginais do gasto (já calculados)
# ----------------------------------------
# Reaproveitando os data frames efeitos_hom e efeitos_mpaf

# Combinar em um único data frame
efeitos_marginais <- bind_rows(
  mutate(efeitos_hom, Desfecho = "Homicídios"),
  mutate(efeitos_mpaf, Desfecho = "Mortes por Arma de Fogo")
)

# Gráfico
p_marg <- ggplot(efeitos_marginais, aes(x = Periodo, y = estimate, color = Desfecho)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(0.3), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Efeito Marginal do Gasto em Segurança sobre os Desfechos",
       x = "Período",
       y = "Efeito Estimado",
       color = "Desfecho") +
  theme_ipea(legend.position = "bottom") +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

print(p_marg)
save_ipeaplot(p_marg, path = dir_homicidios, file.name = "efeitos_marginais_comparado", format = c("png", "eps"))

cat("\nVisualizações salvas nos diretórios especificados.\n")

#################################################################
# Estudo de Evento do modelo de diferenças em diferenças contínuo
# Tendências paralelas prévias e exposição dinâmica
#===============================================

#==================================================================
# ESTUDO DE EVENTO DINÂMICO (ANO A ANO)

# Criando grupos  anos
dados_definidos$quinquenio <- floor(dados_definidos$ano / 5) * 5

dados_definidos$quadrienio <- floor(dados_definidos$ano / 4) * 4

dados_definidos$trienio <- floor(dados_definidos$ano / 3) * 3

dados_definidos$bienio <- floor(dados_definidos$ano / 2) * 2



#linha de base: prevalência de armas de fogo em 2003
dados_definidos <- dados_definidos %>%
  group_by(Estado) %>%
  mutate(base_armas = mean(s_armas_de_fogo[ano >= 1998 & ano <= 2002], na.rm = TRUE)) %>%
  ungroup()

# -----------------------------------------------------------------------------------------------
# DEFINIÇÃO GLOBAL DAS VARIÁVEIS DE CONTROLE (Sem armas de fogo)
# -----------------------------------------------------------------------------------------------

# ===============================================================================================
# 1. MODELO HOMICÍDIOS (Estudo de Evento com DiD Contínuo)
# ===============================================================================================

# Montando a fórmula dinamicamente (incluindo s_cresc_pop que é exclusiva deste modelo)
formula_hmcd <- as.formula(
  paste("log(tx_homicid_ajrand) ~ i(bienio, s_despesa_seguranca, ref = '2002') +", 
        preditores_base, "| Estado + ano")
)

ev_hmcd_cont <- feols(
  formula_hmcd,
  data     = dados_definidos,
  vcov     = "DK",
  panel.id = ~Estado + ano
)

# Testes de Wald
wald_pre_h      <- wald(ev_hmcd_cont, "bienio::(199[0-9]|200[0-2]):s_despesa_seguranca", vcov = "hetero", print = FALSE)
wald_post_h     <- wald(ev_hmcd_cont, "bienio::20(04|06|08|10|12|14|16|18|20|22):s_despesa_seguranca", vcov = "hetero", print = FALSE)
wald_post_ini_h <- wald(ev_hmcd_cont, "bienio::20(04|06|08|10|12|14):s_despesa_seguranca", vcov = "hetero", print = FALSE)
wald_post_rec_h <- wald(ev_hmcd_cont, "bienio::20(16|18|20|22):s_despesa_seguranca", vcov = "hetero", print = FALSE)

# Gráfico de Event Study
iplot_hmcd <- iplot(ev_hmcd_cont,
      main = "Efeito marginal do gasto em segurança por ano sobre homicídios intencionais",
      xlab = "Ano",
      ylab = "Coeficiente do gasto (impacto de 1 d.p. extra)")
legend(
  "bottomright",
  legend = c(
    sprintf("Wald: Pré 2004: p = %.3f", wald_pre_h["p"]),
    sprintf("Wald: Pós 2004: p = %.3f", wald_post_h["p"]),
    sprintf("Wald: 2005-2014: p = %.3f", wald_post_ini_h["p"]),
    sprintf("Wald: Pós 2015-2022: p = %.3f", wald_post_rec_h["p"])
  ),
  bty = "n"
)

summary(ev_hmcd_cont)

coefs_h <- broom::tidy(ev_hmcd_cont)
coefs_h |> dplyr::filter(grepl("199|200[0-2]", term))


# ===============================================================================================
# 2. MODELO MORTES POR ARMAS DE FOGO (Estudo de Evento Contínuo)
# ===============================================================================================

# Montando a fórmula dinamicamente
formula_mpaf <- as.formula(
  paste("log(mortes_PAF) ~ i(bienio, s_despesa_seguranca, ref = '2002') + ", 
        preditores_base, 
        "| Estado + ano")
)

ev_mpaf_cont <- feols(
  formula_mpaf,
  data     = dados_definidos,
  vcov     = "DK",
  panel.id = ~Estado + ano
)

# Testes de Wald
wald_pre_p      <- wald(ev_mpaf_cont, "bienio::(199[0-9]|200[0-2]):s_despesa_seguranca", vcov = "hetero", print = FALSE)
wald_post_p     <- wald(ev_mpaf_cont, "bienio::20(04|06|08|10|12|14|16|18|20|22):s_despesa_seguranca", vcov = "hetero", print = FALSE)
wald_post_ini_p <- wald(ev_mpaf_cont, "bienio::20(04|06|08|10|12|14):s_despesa_seguranca", vcov = "hetero", print = FALSE)
wald_post_rec_p <- wald(ev_mpaf_cont, "bienio::20(16|18|20|22):s_despesa_seguranca", vcov = "hetero", print = FALSE)

# Gráfico de Event Study
iplot_mpaf <- iplot(ev_mpaf_cont,
      main = "Efeito marginal do gasto em segurança por ano sobre mortes por armas de fogo",
      xlab = "Ano",
      ylab = "Coeficiente do gasto (impacto de 1 d.p. extra)")
legend(
  "bottomright",
  legend = c(
    sprintf("Wald: Pré 2004: p = %.3f", wald_pre_p["p"]),
    sprintf("Wald: Pós 2004: p = %.3f", wald_post_p["p"]),
    sprintf("Wald: 2005-2014: p = %.3f", wald_post_ini_p["p"]),
    sprintf("Wald: Pós 2015-2022: p = %.3f", wald_post_rec_p["p"])
  ),
  bty = "n"
)

summary(ev_mpaf_cont)

coefs_p <- broom::tidy(ev_mpaf_cont)
coefs_p |> dplyr::filter(grepl("199|200[0-2]", term))



# ===============================================================================================
# 3. TABELAS DE RESULTADOS (modelsummary)
# ===============================================================================================
library(modelsummary)
library(dplyr)

# Criar tabela simples do teste de Wald
wald_table <- data.frame(
  Teste      = c("Pré-tendência", "Pós (total)", "Pós inicial", "Pós recente"),
  Homicídios = c(wald_pre_h["p"], wald_post_h["p"], wald_post_ini_h["p"], wald_post_rec_h["p"]),
  Mortes_PAF = c(wald_pre_p["p"], wald_post_p["p"], wald_post_ini_p["p"], wald_post_rec_p["p"])
)

wald_table <- wald_table |>
  mutate(across(-Teste, ~sprintf("%.3f", .)))

# Tabela conjunta com modelsummary + p-values do teste de Wald abaixo
modelsummary(
  list(
    "Homicídios" = ev_hmcd_cont,
    "Mortes PAF" = ev_mpaf_cont
  ),
  statistic = "std.error",
  stars     = TRUE,
  gof_omit  = "IC|Log|Adj",
  add_rows  = data.frame(
    term         = c("Pré-tendência (p)", "Pós total (p)", "Pós inicial (p)", "Pós recente (p)"),
    Homicídios   = c(sprintf("%.3f", wald_pre_h["p"]), sprintf("%.3f", wald_post_h["p"]), sprintf("%.3f", wald_post_ini_h["p"]), sprintf("%.3f", wald_post_rec_h["p"])),
    `Mortes PAF` = c(sprintf("%.3f", wald_pre_p["p"]), sprintf("%.3f", wald_post_p["p"]), sprintf("%.3f", wald_post_ini_p["p"]), sprintf("%.3f", wald_post_rec_p["p"]))
  )
)

# =========================================
# Gráfico único de resultados
# =========================================

library(ggplot2)
library(dplyr)
library(fixest)


# 1. Extrair os dados dos modelos
df_h <- broom::tidy(ev_hmcd_cont, conf.int = TRUE) %>% 
  filter(grepl("bienio", term)) %>%
  mutate(modelo = "Homicídios")

df_p <- broom::tidy(ev_mpaf_cont, conf.int = TRUE) %>% 
  filter(grepl("bienio", term)) %>%
  mutate(modelo = "Mortes por PAF")

# 2. Unir os dados
plot_data <- bind_rows(df_h, df_p) %>%
  # Limpeza do termo para extrair o ano (ex: "bienio::2004:s_despesa_seguranca" -> 2004)
  mutate(ano = as.numeric(gsub("bienio::(\\d{4}):.*", "\\1", term)))

# 3. Plotar com ggplot2
# 1. Preparação do texto dinâmico
# Extraímos o p-valor e formatamos (ex: < 0.001 ou com 3 casas decimais)
format_p <- function(x) ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))

texto_wald <- paste0(
  "Wald Test (p-values):\n",
  "Homicídios:\n",
  "  Pré: ", format_p(wald_pre_h["p"]), " | Pós: ", format_p(wald_post_h["p"]), "\n",
  "  2004-2014: ", format_p(wald_post_ini_h["p"]), " | 2015-2024: ", format_p(wald_post_rec_h["p"]), "\n",
  "Mortes PAF:\n",
  "  Pré: ", format_p(wald_pre_p["p"]), " | Pós: ", format_p(wald_post_p["p"]), "\n",
  "  2004-2014: ", format_p(wald_post_ini_p["p"]), " | 2015-2024: ", format_p(wald_post_rec_p["p"])
)

# 2. Plotagem com o texto automático
evento_plot <- ggplot(plot_data, aes(x = ano, y = estimate, color = modelo, shape = modelo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(2002, 2003), linetype = "dotted", color = "black", alpha = 0.6) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.5)) +
  # Inserindo o texto dinâmico (ajuste x e y conforme a escala do seu gráfico)
  annotate("text", 
           x = max(plot_data$ano), 
           y = min(plot_data$conf.low), 
           label = texto_wald, 
           hjust = 1, 
           vjust = 0, 
           size = 3, 
           color = "black", family = "mono") +
  theme_ipea() +
  labs(
    title = "Efeito do Gasto em Segurança sobre mortes por homicídios e por arma de fogo",
    x = "ano",
    y = "coeficiente (Impacto de 1 d.p. extra)",
    color = "variável:",
    shape = "variável:"
  ) +
  theme(legend.position = "top") +
  scale_color_ipea(palette = "Red-Blue-White")

print(evento_plot)
save_ipeaplot(evento_plot,
              file.name = "evento_plot",
              format = c("png", "eps"))

#==========================================
# Teste placebo flexível (desfecho e baseline)
#==========================================

anos_placebo <- c(1992, 1994, 1996, 1998, 2000)

# Função que estima o placebo para uma combinação específica
rodar_placebo <- function(desfecho,   # string: "tx_homicid_ajrand" ou "mortes_PAF"
                          baseline,  # string
                          label_desfecho,
                          label_baseline) {
  
  # Construir a fórmula dinamicamente
  formula_str <- paste0(
    "log(", desfecho, ") ~ ", baseline, " * Lei_placebo + ",
    preditores_base," | Estado + ano"
  )
  formula_placebo <- as.formula(formula_str)
  
  map_df(anos_placebo, function(a) {
    
    dados_tmp <- dados_definidos %>%
      filter(ano < 2004) %>%
      mutate(Lei_placebo = ifelse(ano >= a, 1, 0))
    
    mod <- feols(
      formula_placebo,
      data = dados_tmp,
      vcov = "DK",
      panel.id = ~Estado + ano,
      warn = FALSE, quiet = TRUE   # evitar mensagens excessivas
    )
    
    # Nome do coeficiente da interação (ex: "armas_2003:Lei_placebo")
    coef_name <- paste0(baseline, ":Lei_placebo")
    coef_val <- coef(mod)[coef_name]
    se_val   <- sqrt(vcov(mod)[coef_name, coef_name])
    
    data.frame(
      ano_placebo = a,
      coef  = coef_val,
      se    = se_val,
      desfecho = label_desfecho,
      baseline = label_baseline
    )
  })
}

# Rodar para as 4 combinações
placebo_hom_seg <- rodar_placebo(
  desfecho = "tx_homicid_ajrand",
  baseline = "s_despesa_seguranca",   # ✅ CORRETO
  label_desfecho = "Homicídios",
  label_baseline = "gasto em segurança"
)
placebo_paf_seg   <- rodar_placebo(desfecho = "mortes_PAF", baseline = "s_despesa_seguranca",
                                   label_desfecho = "Mortes PAF", label_baseline =  "Segurança anual")

# Combinar tudo
resultados_placebo_todos <- bind_rows(
  placebo_hom_seg,
  placebo_paf_seg
) %>%
  mutate(
    t     = coef / se,
    p_value = 2 * pnorm(-abs(t)),
    ci_low  = coef - 1.96 * se,
    ci_high = coef + 1.96 * se
  )

print(resultados_placebo_todos)

ggplot(resultados_placebo_todos, aes(x = ano_placebo, y = coef, color = baseline)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3,
                position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~desfecho, scales = "free_y") +
  labs(title    = "Testes de placebo (datas fictícias antes de 2004)",
       x        = "Ano do placebo",
       y        = "Coeficiente da interação Baseline (1998-2002) × Lei Placebo",
       color    = "Baseline") +
  theme_ipea(base_size = 13) +
  theme(legend.position = "bottom") +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

# ==============================================================================
# ANÁLISE DO ESTATUTO DO DESARMAMENTO
# Substituições metodológicas:
#   (A) Mediação causal via Imai et al. (2010) com análise de sensibilidade
# ==============================================================================

# ==============================================================================
# PRESSUPOSTOS: dados_definidos já está carregado conforme script anterior,
# com todas as variáveis padronizadas (s_...) e a dummy Lei_desarmamento_10826_2003
# ==============================================================================

# ==============================================================================
# BLOCO A — MEDIAÇÃO CAUSAL (Imai et al., 2010)
# ==============================================================================
# Estratégia: aplicar mediate() sobre dados com efeitos fixos removidos (demeaning)
# O demeaning por Estado e ano equivale a condicionar nos efeitos fixos bidirecionais,
# tornando a mediação causal compatível com a estrutura de painel.
#
# Hipótese causal:
#   Gasto em segurança → Estoque de armas (mediador) → Homicídios / Mortes PAF
#   O Estatuto modifica o efeito do gasto sobre o estoque de armas (caminho α)
#
# Suposição de identificação (ignorabilidade sequencial):
#   Condicional ao tratamento (gasto) e aos controles, o mediador (armas)
#   é como se fosse aleatoriamente atribuído — i.e., não há confundidores
#   não observados do par mediador→desfecho que sejam independentes do tratamento.
# ==============================================================================

# -----------------------------------------------------------------------
# A.1 — Demeaning bidirecional manual (remoção de efeitos fixos Estado + ano)
# -----------------------------------------------------------------------
# O estimador within bidirecional é algebricamente equivalente a:
#   y_demeaned = y - ȳ_i - ȳ_t + ȳ
# onde ȳ_i é a média da unidade i, ȳ_t é a média do período t, e ȳ é a
# média global. Essa transformação é aplicada a cada variável antes de
# passar ao mediate(), contornando a limitação do Within() do plm com
# pseries matriciais.


# -----------------------------------------------------------------------
# A.2 — Análise de mediação separada por período (pré e pós-Estatuto)
# -----------------------------------------------------------------------
# Separamos os subconjuntos porque queremos comparar o ACME e o ADE
# antes e depois da lei, que é o objeto central de análise


# ==============================================================================
# BLOCO A — MEDIAÇÃO CAUSAL (Imai et al., 2010)  [VERSÃO REVISADA]
# Principais mudanças em relação à versão anterior:
#   1. Bootstrap não-paramétrico (boot = TRUE, boot.ci.type = "perc") no lugar
#      do estimador quasi-Bayesiano — mais robusto em amostras finitas e não
#      depende de suposições distributivas sobre os coeficientes.
#   2. Orientações detalhadas de interpretação inline para cada saída:
#      (a) summary(mediate): ACME, ADE, Efeito Total e Proporção Mediada
#      (b) summary(medsens): tabela de ρ, rho crítico, R²_M~·R²_Y~
# ==============================================================================
#
# CADEIA CAUSAL ESTIMADA:
#
#   Gasto em Segurança ──α──▶ Estoque de Armas ──β──▶ { Homicídios / Mortes PAF }
#                  └────────────── δ (ADE) ───────────────────────┘
#
#   ACME  = α × β  (efeito indireto que transita pelo estoque de armas)
#   ADE   = δ      (efeito direto: todos os demais canais, ex. dissuasão policial)
#   Total = ACME + ADE
#
# SUPOSIÇÃO DE IDENTIFICAÇÃO (ignorabilidade sequencial, Imai et al., 2010):
#   Condicionado ao gasto em segurança (tratamento) e aos controles X, o
#   estoque de armas (mediador) é como se fosse aleatório — i.e., não existe
#   confundidor não observado que afete SIMULTANEAMENTE o mediador e o desfecho
#   independentemente do tratamento. A análise de sensibilidade (medsens) avalia
#   a robustez dos resultados a violações desta suposição.
#
# POR QUE BOOTSTRAP NÃO-PARAMÉTRICO?
#   O estimador quasi-Bayesiano (boot = FALSE) amostra diretamente da
#   distribuição assintótica dos coeficientes. É computacionalmente mais rápido,
#   mas:
#     - Assume normalidade assintótica dos coeficientes
#     - Pode ser impreciso em amostras moderadas (n ~ 500-700 aqui)
#     - Produz ICs simétricos por construção
#   O bootstrap não-paramétrico (boot = TRUE) reamostrar as OBSERVAÇÕES,
#   recalcula os modelos M e Y em cada amostra e extrai o ACME de cada réplica.
#   Produz ICs assimétricos, sem suposições distributivas, e é preferido pela
#   literatura de simulação (Imai et al., 2010; MacKinnon et al., 2004).
#   Desvantagem: custo computacional ~10-20x maior.
#   Recomendação: sims = 2000 para versão final de artigo; sims = 500 para
#   exploração inicial.
# ==============================================================================


# Criar variáveis de log antes do demeaning

dados_definidos$log_homicidios <- log(dados_definidos$tx_homicid_ajrand)
dados_definidos$log_mortes_PAF <- log(dados_definidos$mortes_PAF)


# Extração automática das variáveis individuais (remove espaços e quebra pelo "+")
vars_base_limpas <- trimws(strsplit(preditores_base, "\\+")[[1]])

# Criação automática da string com o prefixo "dm_" para as fórmulas de mediação
preditores_base_dm <- paste0("dm_", vars_base_limpas, collapse = " + ")

# ==============================================================================
# ANÁLISE DO ESTATUTO DO DESARMAMENTO — VERSÃO REVISADA
# Substituições metodológicas:
#   (A) Mediação causal via Imai et al. (2010) com análise de sensibilidade
# ==============================================================================

# -----------------------------------------------------------------------
# A.1 — Demeaning bidirecional manual (remoção de efeitos fixos Estado + ano)
# -----------------------------------------------------------------------
demean_twoway <- function(x, unit, time) {
  media_unit   <- ave(x, unit, FUN = function(z) mean(z, na.rm = TRUE))
  media_time   <- ave(x, time, FUN = function(z) mean(z, na.rm = TRUE))
  media_global <- mean(x, na.rm = TRUE)
  x - media_unit - media_time + media_global
}


# Montagem dinâmica do vetor de variáveis a demeanar
vars_demean <- c(
  "s_despesa_seguranca",
  "s_armas_de_fogo",
  vars_base_limpas, # Suas variáveis base injetadas automaticamente aqui
  "Lei_desarmamento_10826_2003",
  "log_homicidios", 
  "log_mortes_PAF"
)

# Criar dataset demeaned em um único passo
dados_dm <- dados_definidos

for (v in vars_demean) {
  dados_dm[[paste0("dm_", v)]] <- demean_twoway(
    dados_definidos[[v]],
    dados_definidos$Estado,
    dados_definidos$ano
  )
}

# Verificar: médias das variáveis demeaned devem ser ~0
cat("Verificação do demeaning (médias devem ser ~0):\n")
vars_check <- paste0("dm_", vars_demean)
print(round(sapply(dados_dm[vars_check], mean, na.rm = TRUE), 8))

# -----------------------------------------------------------------------
# A.2 — Análise de mediação separada por período (pré e pós-Estatuto)
# -----------------------------------------------------------------------
dados_dm_pre  <- filter(dados_dm, Lei_desarmamento_10826_2003 == 0)
dados_dm_pos  <- filter(dados_dm, Lei_desarmamento_10826_2003 == 1)

# ==============================================================================
# A.3 — Função de estimação com bootstrap não-paramétrico e Fórmulas Dinâmicas
# ==============================================================================
estimar_mediacao_causal <- function(dados_periodo,
                                    outcome_dm,
                                    sims    = 5000,
                                    seed    = 42,
                                    ci.type = "perc") {
  
  # 1. Fórmula do mediador construída dinamicamente usando a variáveis_base_dm
  formula_m <- as.formula(paste(
    "dm_s_armas_de_fogo ~ dm_s_despesa_seguranca +", 
    preditores_base_dm
  ))
  
  modelo_m <- lm(formula_m, data = dados_periodo)
  modelo_m$call$formula <- formula_m   # Embala o objeto no call para evitar erro de escopo no bootstrap
  
  # 2. Fórmula do desfecho construída dinamicamente
  formula_y <- as.formula(paste(
    outcome_dm, "~ dm_s_despesa_seguranca + dm_s_armas_de_fogo +", 
    preditores_base_dm
  ))
  
  modelo_y <- lm(formula_y, data = dados_periodo)
  modelo_y$call$formula <- formula_y   # Embala o objeto no call para evitar erro de escopo no bootstrap
  
  # Bootstrap não-paramétrico
  set.seed(seed)
  med_result <- mediate(
    modelo_m, modelo_y,
    treat        = "dm_s_despesa_seguranca",
    mediator     = "dm_s_armas_de_fogo",
    sims         = sims,
    boot         = TRUE,          
    boot.ci.type = ci.type        
  )
  
  return(list(modelo_m = modelo_m, modelo_y = modelo_y, mediacao = med_result))
}

# ==============================================================================
# A.4 — Estimação: HOMICÍDIOS (pré e pós-Estatuto)
# ==============================================================================
cat("\n=== MEDIAÇÃO CAUSAL BOOTSTRAP: HOMICÍDIOS (PRÉ-ESTATUTO) ===\n")
med_hom_pre <- estimar_mediacao_causal(dados_dm_pre, "dm_log_homicidios", sims = 5000)
summary(med_hom_pre$mediacao)

cat("\n=== MEDIAÇÃO CAUSAL BOOTSTRAP: HOMICÍDIOS (PÓS-ESTATUTO) ===\n")
med_hom_pos <- estimar_mediacao_causal(dados_dm_pos, "dm_log_homicidios", sims = 2000)
summary(med_hom_pos$mediacao)

# ==============================================================================
# A.5 — Estimação: MORTES POR ARMAS DE FOGO (pré e pós-Estatuto)
# ==============================================================================
cat("\n=== MEDIAÇÃO CAUSAL BOOTSTRAP: MORTES PAF (PRÉ-ESTATUTO) ===\n")
med_paf_pre <- estimar_mediacao_causal(dados_dm_pre, "dm_log_mortes_PAF", sims = 2000)
summary(med_paf_pre$mediacao)

cat("\n=== MEDIAÇÃO CAUSAL BOOTSTRAP: MORTES PAF (PÓS-ESTATUTO) ===\n")
med_paf_pos <- estimar_mediacao_causal(dados_dm_pos, "dm_log_mortes_PAF", sims = 2000)
summary(med_paf_pos$mediacao)

# ==============================================================================
# A.6 — Análise de sensibilidade (Imai & Yamamoto, 2010)
# ==============================================================================
extrair_sens <- function(sens_obj) {
  rho  <- sens_obj$rho
  acme <- sens_obj$d0  
  rho_critico <- NA_real_
  r2_critico  <- NA_real_
  r2p_critico <- NA_real_
  
  mudanca <- which(diff(sign(acme)) != 0)
  
  if (length(mudanca) == 0) {
    cat("  ACME não cruza zero no intervalo testado.\n")
  } else {
    i <- mudanca[1]
    prop <- abs(acme[i]) / (abs(acme[i]) + abs(acme[i+1]))
    rho_critico <- rho[i] + prop * (rho[i+1] - rho[i])
    
    if (!is.null(sens_obj$R2star.d.thresh)) {
      r2_critico  <- sens_obj$R2star.d.thresh
    } else {
      r2_critico <- sens_obj$R2star.prod[i] + prop * (sens_obj$R2star.prod[i+1] - sens_obj$R2star.prod[i])
    }
    
    if (!is.null(sens_obj$R2tilde.d.thresh)) {
      r2p_critico <- sens_obj$R2tilde.d.thresh
    } else {
      r2p_critico <- sens_obj$R2tilde.prod[i] + prop * (sens_obj$R2tilde.prod[i+1] - sens_obj$R2tilde.prod[i])
    }
  }
  
  list(rho_critico = rho_critico, R2_bruto = r2_critico, R2_parcial = r2p_critico)
}

interpretar_sens <- function(sens_obj, rotulo, rotulo_x) {
  s <- extrair_sens(sens_obj)
  cat("\n--- INTERPRETAÇÃO DA SENSIBILIDADE (", rotulo, ") ---\n", sep = "")
  
  if (is.na(s$rho_critico)) {
    cat("ρ crítico: não determinado no intervalo testado.\n")
  } else {
    r2p <- if (!is.na(s$R2_parcial)) s$R2_parcial else s$R2_bruto
    robustez <- dplyr::case_when(
      abs(s$rho_critico) < 0.20 ~ "FRAGIL     — qualquer confundidor leve anularia o ACME",
      abs(s$rho_critico) < 0.40 ~ "MODERADA   — exigiria confundidor de magnitude substancial",
      abs(s$rho_critico) < 0.60 ~ "ROBUSTA    — confundidor forte necessário para anular",
      TRUE                        ~ "MUITO ROBUSTA — confundidor extremo seria necessário"
    )
    cat(sprintf(
      "  ρ crítico         : %.3f\n  R²_M~·R²_Y~      : %.4f (%.1f%% da variância residual)\n  Robustez          : %s\n\n  Interpretação:\n  O ACME só seria zero se a correlação entre os resíduos dos modelos\n  de Armas e de %s fosse %.3f.\n  Isso requer um confundidor que explique ao menos %.1f%% da variância\n  residual de CADA modelo (além dos controles já incluídos).\n",
      s$rho_critico, r2p, r2p * 100, robustez, rotulo, s$rho_critico, r2p * 100
    ))
  }
  
  plot(sens_obj, main = paste0("Sensibilidade do ACME — ", rotulo), xlab = rotulo_x, ylab = "ACME estimado (IC 95%)")
  abline(h = 0, lty = 2, col = "red", lwd = 1.5)
  
  if (!is.na(s$rho_critico)) {
    abline(v = s$rho_critico, lty = 3, col = "blue", lwd = 1.2)
    r2p_label <- if (!is.na(s$R2_parcial)) s$R2_parcial else s$R2_bruto
    mtext(paste0("ρ crítico ≈ ", round(s$rho_critico, 3), " | R²_M~R²_Y~ ≈ ", round(r2p_label, 4)), side = 3, line = 0.3, cex = 0.82, col = "blue")
  }
  return(invisible(s))
}

cat("\n=== SENSIBILIDADE BOOTSTRAP: HOMICÍDIOS PÓS-ESTATUTO ===\n")
sens_hom_pos <- medsens(med_hom_pos$mediacao, rho.by = 0.05, effect.type = "indirect")
plot(sens_hom_pos)
summary(sens_hom_pos)
s_hom_pos <- interpretar_sens(sens_hom_pos, rotulo = "Homicídios Pós-Estatuto", rotulo_x = "ρ (correlação entre resíduos de Armas e Homicídios)")

cat("\n=== SENSIBILIDADE BOOTSTRAP: MORTES PAF PÓS-ESTATUTO ===\n")
sens_paf_pos <- medsens(med_paf_pos$mediacao, rho.by = 0.05, effect.type = "indirect")
summary(sens_paf_pos)
s_paf_pos <- interpretar_sens(sens_paf_pos, rotulo = "Mortes PAF Pós-Estatuto", rotulo_x = "ρ (correlação entre resíduos de Armas e Mortes PAF)")

cat("\n=== SENSIBILIDADE: HOMICÍDIOS PRÉ-ESTATUTO (baseline) ===\n")
sens_hom_pre <- medsens(med_hom_pre$mediacao, rho.by = 0.05, effect.type = "indirect")
summary(sens_hom_pre)
s_hom_pre <- interpretar_sens(sens_hom_pre, rotulo = "Homicídios Pré-Estatuto", rotulo_x = "ρ (correlação entre resíduos de Armas e Homicídios)")

cat("\n=== SENSIBILIDADE: MORTES PAF PRÉ-ESTATUTO ===\n")
sens_paf_pre <- medsens(med_paf_pre$mediacao, rho.by = 0.05, effect.type = "indirect")
summary(sens_paf_pre)
s_paf_pre <- interpretar_sens(sens_paf_pre, rotulo = "Mortes PAF Pré-Estatuto", rotulo_x = "ρ (correlação entre resíduos de Armas e Mortes PAF)")

# ==============================================================================
# A.7 — Tabela-resumo de sensibilidade comparativa
# ==============================================================================
tabela_sensibilidade <- data.frame(
  Cenario = c("Homicídios — Pré-Estatuto", "Homicídios — Pós-Estatuto", "Mortes PAF — Pré-Estatuto", "Mortes PAF — Pós-Estatuto"),
  rho_critico = c(s_hom_pre$rho_critico, s_hom_pos$rho_critico, s_paf_pre$rho_critico, s_paf_pos$rho_critico),
  R2_parcial = c(s_hom_pre$R2_parcial, s_hom_pos$R2_parcial, s_paf_pre$R2_parcial, s_paf_pos$R2_parcial)
) %>%
  mutate(
    robustez = dplyr::case_when(
      is.na(rho_critico)          ~ "Não determinado",
      abs(rho_critico) < 0.20     ~ "Frágil (|ρ| < 0.20)",
      abs(rho_critico) < 0.40     ~ "Moderada (0.20 ≤ |ρ| < 0.40)",
      abs(rho_critico) < 0.60     ~ "Robusta (0.40 ≤ |ρ| < 0.60)",
      TRUE                        ~ "Muito robusta (|ρ| ≥ 0.60)"
    ),
    rho_critico = round(rho_critico, 3),
    R2_parcial  = round(R2_parcial,  4)
  )

cat("\n=== TABELA COMPARATIVA DE SENSIBILIDADE ===\n")
print(tabela_sensibilidade)
write.csv(tabela_sensibilidade, file = "tabela_sensibilidade.csv", row.names = FALSE)

# ==============================================================================
# A.8 — Visualização comparativa dos efeitos causais (ACME e ADE)
# ==============================================================================
extrair_med <- function(med_obj, periodo, desfecho) {
  s <- summary(med_obj$mediacao)
  data.frame(
    Desfecho   = desfecho,
    Periodo    = periodo,
    Efeito     = c("ACME (Indireto)", "ADE (Direto)", "Total"),
    Estimativa = c(s$d.avg,    s$z.avg,    s$tau.coef),
    IC_low     = c(s$d.avg.ci[1], s$z.avg.ci[1], s$tau.ci[1]),
    IC_high    = c(s$d.avg.ci[2], s$z.avg.ci[2], s$tau.ci[2]),
    p_value    = c(s$d.avg.p,  s$z.avg.p,  s$tau.p)
  )
}

resultados_mediacao_causal <- bind_rows(
  extrair_med(med_hom_pre, "Pré-Estatuto", "Homicídios"),
  extrair_med(med_hom_pos, "Pós-Estatuto", "Homicídios"),
  extrair_med(med_paf_pre, "Pré-Estatuto", "Mortes PAF"),
  extrair_med(med_paf_pos, "Pós-Estatuto", "Mortes PAF")
)

cat("\n=== RESULTADOS CONSOLIDADOS DA MEDIAÇÃO CAUSAL ===\n")
print(resultados_mediacao_causal)
write.csv(resultados_mediacao_causal, file = "resultados_mediacao_causal.csv", row.names = FALSE)

# Gráfico de efeitos causais com ggplot2
p_mediacao_causal <- ggplot(
  resultados_mediacao_causal,
  aes(x = Efeito, y = Estimativa, color = Periodo, shape = Periodo)
) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = IC_low, ymax = IC_high), position = position_dodge(0.4), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Desfecho) +
  labs(
    title    = "Mediação Causal: Efeitos do Gasto em Segurança via Estoque de Armas",
    subtitle = "ACME = Efeito Causal Médio de Mediação; ADE = Efeito Direto Médio",
    x        = "Tipo de Efeito",
    y        = "Estimativa (escala log)",
    color    = "Período", shape = "Período",
    caption  = "IC 95% via bootstrap não-paramétrico (percentil), sims = 2000.\nDados demeaned por Estado e ano (efeitos fixos bidirecionais).\nImai, Keele, Tingley & Yamamoto (2010)."
  ) +
  theme_ipea() +
  theme(legend.position = "bottom") +
  scale_color_ipea(palette = "Red-Blue-White") +
  scale_fill_ipea(palette = "Red-Blue-White")

print(p_mediacao_causal)

save_ipeaplot(p_mediacao_causal, file.name = "mediacao_causal_bootstrap", format = c("png", "eps"))

# ==============================================================================
# A.9 — Salvar objetos
# ==============================================================================
save(
  med_hom_pre, med_hom_pos,
  med_paf_pre, med_paf_pos,
  sens_hom_pre, sens_hom_pos, sens_paf_pos,
  resultados_mediacao_causal,
  tabela_sensibilidade,
  file = file.path(dir_homicidios, "../resultados_mediacao_bootstrap.RData"))