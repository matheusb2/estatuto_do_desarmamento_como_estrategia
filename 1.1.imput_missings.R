library(mice)
library(miceadds)
library(corrplot)
library(tidyverse)
library(glmnet)


# ==========================================================================
# Preparar dados para a imputação múltipla
# ==========================================================================

# Garantir que 'Estado' seja fator e 'ano' numérico (já no início)
dados_BR_1980_2024 <- painel_final %>%
  mutate(Estado = as.factor(Estado),
         ano = as.numeric(ano),
         seg_pub_estad_incompleto = as.numeric(gasto_seg_per_capita),
         gasto_social = as.numeric(gasto_edu_per_capita + gasto_saúde_per_capita))


# Criar versões log das variáveis a imputar
dados_BR_1980_2024 <- dados_BR_1980_2024 %>%
  mutate(sq_ano = ano ^2)


# ---------------------------
# MATRIZ DE CORRELAÇÃO (inalterada)
# ---------------------------
variaveis <- c("ano", "desemprego_incompleto", "homens_jovens", "jovens_15_29", "Gini",
               "cresc_populacional", "densidade_populacional", "urbanizacao", 
               "PIB_per_capita", "var_PIB", "energia_per_capita",
               "setor_publico_per_capita",
               "industria_PIB", "seg_pub_estad_incompleto", "receita_per_capita",
               "despesa_per_capita", "uso_drogas", "uso_alcool", "armas_de_fogo_cook",
               "gasto_edu_per_capita",
               "gasto_prevassist_per_capita")

dados_cor <- dados_BR_1980_2024[, variaveis]
matriz_cor <- cor(dados_cor, use = "pairwise.complete.obs")
corrplot(matriz_cor, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         number.cex = 0.7, diag = FALSE)

print(matriz_cor)
write.csv(matriz_cor, "/home/matheus/Documentos/IPEA/modelos/imputar_missing/matriz_correlacao.csv")

#===========================================================================
# Seleção de variáveis por LASSO
#===========================================================================

selecao_cor_lasso <- function(df, var_dep, candidatas, limite_cor = 0.10, 
                              sempre_incluir = c("ano"), alpha_lasso = 1, 
                              nfolds = 10, lambda_regra = "min") {
  
  cat("\n=== Seleção para:", var_dep, "===\n\n")
  
  # Remove as variáveis 'sempre_incluir' da lista de candidatas para evitar duplicatas
  candidatas <- setdiff(candidatas, sempre_incluir)
  
  # ------------------------------------------------------------------
  # ESTÁGIO 1 — Filtro por correlação
  # ------------------------------------------------------------------
  # Garantir que o df seja um tibble ou data.frame padrão
  df_limpo <- df %>% dplyr::select(all_of(c(var_dep, candidatas, sempre_incluir))) %>% na.omit()
  
  cors_abs <- sapply(candidatas, function(v) {
    abs(cor(df_limpo[[v]], df_limpo[[var_dep]], use = "pairwise.complete.obs"))
  })
  
  cat("-- Estágio 1: filtro por |r| >=", limite_cor, "--\n")
  print(round(sort(cors_abs, decreasing = TRUE), 3))
  
  aprovadas_cor <- names(cors_abs[cors_abs >= limite_cor])
  rejeitadas_cor <- setdiff(candidatas, aprovadas_cor)
  
  entrada_lasso <- c(aprovadas_cor, sempre_incluir)
  
  # ------------------------------------------------------------------
  # ESTÁGIO 2 — LASSO
  # ------------------------------------------------------------------
  cat("-- Estágio 2: LASSO (alpha =", alpha_lasso, ") --\n")
  
  X <- as.matrix(df_limpo[, entrada_lasso])
  y <- df_limpo[[var_dep]]
  
  # Penalty factor: 0 para sempre_incluir (não penaliza)
  penalty_factor <- ifelse(colnames(X) %in% sempre_incluir, 0, 1)
  
  set.seed(42)
  cv_fit <- cv.glmnet(X, y, alpha = alpha_lasso, nfolds = nfolds, penalty.factor = penalty_factor)
  
  lambda_escolhido <- if (lambda_regra == "min") cv_fit$lambda.min else cv_fit$lambda.1se
  
  coefs <- coef(cv_fit, s = lambda_escolhido)
  # Extrair apenas nomes das variáveis com coeficientes != 0
  selecionadas_lasso <- rownames(coefs)[which(coefs[, 1] != 0)]
  selecionadas_lasso <- setdiff(selecionadas_lasso, "(Intercept)")
  
  rejeitadas_lasso <- setdiff(entrada_lasso, selecionadas_lasso)
  
  cat(sprintf("Lambda (%s): %.4f\n", lambda_regra, lambda_escolhido))
  cat("Selecionadas:", paste(selecionadas_lasso, collapse = ", "), "\n\n")
  
  return(list(
    preditores = selecionadas_lasso,
    formula    = as.formula(paste(var_dep, "~", paste(selecionadas_lasso, collapse = " + "))),
    cv_lasso   = cv_fit
  ))
}

# ============================================================
# APLICAÇÃO: DESEMPREGO (desemprego_incompleto)
# ============================================================

# Base de treino: anos >= 1991 e desemprego_incompleto não NA
treino_desemprego <- dados_BR_1980_2024 %>%
  filter(ano >= 1991, !is.na(desemprego_incompleto))

# Lista de candidatas quantitativas
candidatas_desemprego <- c("homens_jovens", "Gini", "homens_jovens", "uso_drogas",
                           "uso_alcool", "armas_de_fogo_cook", "txhomicid", "urbanizacao",
                           "PIB_per_capita", "cresc_populacional", "densidade_populacional",
                           "var_PIB",  "industria_PIB",
                           "receita_per_capita", "gasto_edu_per_capita")

res_lasso_desemp <- selecao_cor_lasso(
  df             = treino_desemprego,
  var_dep        = "desemprego_incompleto",
  candidatas     = candidatas_desemprego,
  limite_cor     = 0.15,
  sempre_incluir = "ano",
  alpha_lasso    = 1,       # LASSO puro; use 0.5 para Elastic Net
  lambda_regra   = "min"
)

# Salvar fórmula final para uso posterior
preditores_desemprego <- res_lasso_desemp$preditores

# =============================================================================
# IMPUTAÇÃO MÚLTIPLA – DESEMPREGO (1980–1990)
# =============================================================================

# Preditores selecionados (baseados nas correlações e plausibilidade)
pred_sel <- c(preditores_desemprego)

# Subconjunto para imputação
dados_mi <- dados_BR_1980_2024 %>%
  dplyr::select(Estado, ano, desemprego_incompleto, all_of(pred_sel))

# Configurar métodos
meth <- make.method(dados_mi)
meth["desemprego_incompleto"] <- "pmm"   # regressão linear

# Matriz de preditores: removemos 'Estado' da predição de desemprego_incompleto
pred <- make.predictorMatrix(dados_mi)

# Executar imputação múltipla (5 conjuntos)
imp <- mice(dados_mi, method = meth, predictorMatrix = pred, 
            m = 25, maxit = 5, seed = 123, print = FALSE)

plot(imp)           # traços das médias e desvios por iteração
densityplot(imp)    # distribuição imputada vs observada
stripplot(imp, desemprego_incompleto ~ .imp)
summary(imp)

# Identificar linhas de 1980-1990 (onde desemprego é NA)
linhas_pred <- which(dados_BR_1980_2024$ano >= 1979 & dados_BR_1980_2024$ano <= 1990)

# Extrair todas as imputações completas
imp_list <- complete(imp, action = "all")

# Para cada linha, calcular média e desvio padrão na escala log
imp_values <- sapply(imp_list, function(df) df[linhas_pred, "desemprego_incompleto"])
imp_mean <- rowMeans(imp_values, na.rm = TRUE)
imp_sd   <- apply(imp_values, 1, sd, na.rm = TRUE)

# Criar data frame 'predicao' com as estimativas
predicao <- dados_BR_1980_2024[linhas_pred, c("Estado", "ano")]
predicao <- predicao %>%
  mutate(desemprego_imputado_mi   = imp_mean,
         desemprego_imputado_lwr  = imp_mean - 1.96 * imp_sd,
         desemprego_imputado_upr  = imp_mean + 1.96 * imp_sd)

# Visualizar resumo
summary(predicao$desemprego_imputado_mi)

# =============================================================================
# IMPUTAÇÃO MÚLTIPLA – SEGURANÇA PÚBLICA (1980–1984)
# =============================================================================


# Base de treino: anos >= 1985 e seg_pub_estad_incompleto não NA
treino_seg <- dados_BR_1980_2024 %>%
  filter(ano >= 1985, !is.na(seg_pub_estad_incompleto))

# Preditores candidatos (pode usar a mesma lista inicial)
candidatas_seg <- c( "ano",  "homens_jovens", "Gini", "homens_jovens", "uso_drogas",
                     "uso_alcool", "armas_de_fogo_cook", "txhomicid", "urbanizacao",
                     "PIB_per_capita", "receita_per_capita")

resultado_seg <- selecao_cor_lasso(
  df             = treino_seg,
  var_dep        = "seg_pub_estad_incompleto",
  candidatas     = candidatas_seg,
  limite_cor     = 0.15,
  sempre_incluir = "ano")

preditores_seg <- resultado_seg$preditores
pred_sel_seg <- c("ano", preditores_seg)

dados_mi_seg <- dados_BR_1980_2024 %>%
  dplyr::select(Estado, ano, seg_pub_estad_incompleto, all_of(pred_sel_seg))

meth_seg <- make.method(dados_mi_seg)
meth_seg["seg_pub_estad_incompleto"] <- "pmm"

pred_seg <- make.predictorMatrix(dados_mi_seg)

imp_seg <- mice(dados_mi_seg, method = meth_seg, predictorMatrix = pred_seg,
                m = 20, maxit = 5, seed = 456, print = FALSE)

plot(imp_seg)           # traços das médias e desvios por iteração
densityplot(imp_seg)    # distribuição imputada vs observada
stripplot(imp_seg, seg_pub_estad_incompleto ~ .imp)

linhas_pred_seg <- which(dados_BR_1980_2024$ano >= 1979 & dados_BR_1980_2024$ano <= 1984)
imp_list_seg <- complete(imp_seg, action = "all")
log_imp_seg <- sapply(imp_list_seg, function(df) df[linhas_pred_seg, "seg_pub_estad_incompleto"])
imp_mean_seg <- rowMeans(log_imp_seg, na.rm = TRUE)

predicao_seg <- dados_BR_1980_2024[linhas_pred_seg, c("Estado", "ano")]
predicao_seg <- predicao_seg %>%
  mutate(seg_pub_estad_imputado_mi = imp_mean_seg)

# Visualizar resumo
summary(predicao_seg$seg_pub_estad_imputado_mi)


# =============================================================================
# JUNTAR TUDO NO DATASET ORIGINAL
# =============================================================================

dados_BR_1980_2024 <- dados_BR_1980_2024 %>%
  left_join(predicao %>% dplyr::select(Estado, ano, 
                                desemprego_imputado_mi, 
                                desemprego_imputado_lwr, 
                                desemprego_imputado_upr), 
            by = c("Estado", "ano")) %>%
  left_join(predicao_seg %>% dplyr::select(Estado, ano, seg_pub_estad_imputado_mi),
            by = c("Estado", "ano"))

# Criar variáveis completas (observado + imputado)
dados_BR_1980_2024 <- dados_BR_1980_2024 %>%
  mutate(
    desemprego_imputado = if_else(
      ano >= 1991 & !is.na(desemprego_incompleto),
      desemprego_incompleto,
      desemprego_imputado_mi
    ),
    seg_pub_estad_imputado = if_else(
      ano >= 1985 & !is.na(seg_pub_estad_incompleto),
      seg_pub_estad_incompleto,
      seg_pub_estad_imputado_mi
    )
  )

# Salvar resultado final
write.csv(dados_BR_1980_2024, 
          "/home/matheus/Documentos/IPEA/modelos/imputar_missing/dados_BR_1980_2024_com_imputacao_mi.csv", 
          row.names = FALSE)

# =================================================================================
# Integrar com os dados de mortalidade corrigidos
# =================================================================================

#==============================================================================
#Carregar nova base de dados
#==============================================================================

# Endereço dos dados imputados

dados_imputados <- "~/Documentos/IPEA/modelos/imputar_missing/dados_BR_1980_2024_com_imputacao_mi.csv"

dados_BR_1980_2024_com_imputacao_mi <- read.csv(dados_imputados)

#Estatísticas descritivas das variáveis
summary(dados_BR_1980_2024)

# Verificar estrutura dos dados
str(dados_BR_1980_2024_com_imputacao_mi)

# Nomes das colunas
print(names(dados_BR_1980_2024_com_imputacao_mi))

# ==============================================================================
# 3. FAZER JOIN ENTRE OS DATASETS (CÓDIGO CORRIGIDO)
# ==============================================================================

# Primeiro, vamos remover colunas duplicadas antes do join
# Manter as colunas de dados_BR_1980_2024, exceto aquelas que já existem em dados_mortes
colunas_para_manter_BR <- setdiff(names(dados_BR_1980_2024_com_imputacao_mi), 
                                  c("Estado", "ano", "Macrorregiao", "populacao"))

# Criar um subset de dados_BR_1980_2024 apenas com colunas únicas
dados_BR_para_join <- dados_BR_1980_2024_com_imputacao_mi %>%
  dplyr::select(Estado, ano, all_of(colunas_para_manter_BR))

# Fazer o join
dados_completos <- dados_mortes %>%
  left_join(dados_BR_para_join, by = c("Estado", "ano"))

# Verificar resultado do join
cat("\nResultado do join:\n")
cat("Linhas em dados_mortes:", nrow(dados_mortes), "\n")
cat("Linhas após join:", nrow(dados_completos), "\n")
cat("Colunas em dados_mortes:", ncol(dados_mortes), "\n")
cat("Colunas em dados_completos:", ncol(dados_completos), "\n")

# Verificar algumas linhas para confirmar o join
cat("\nPrimeiras linhas do dataset combinado (apenas algumas colunas):\n")
print(head(dados_completos %>% 
             dplyr::select(Estado, ano, agressoes_e_confrontos, desemprego_imputado, Gini), 10))

#salvar dados

arquivo__dados_completos <- "/home/matheus/Documentos/IPEA/modelos/dados_completos/dados_completos_BR_80_22.csv"

write.csv(dados_completos, arquivo__dados_completos,
          row.names = F)