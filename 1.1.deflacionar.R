library(tidyverse)

# Carregar dados

dados_BR_1980_2024 <- read.csv("~/Documentos/IPEA/modelos/interação desigualdade homicídios/dados_BR_1980_2024.csv", 
                               dec=",")
str(dados_BR_1980_2024)

# ==============================================================================
# DEFLACIONAMENTO COM ENCADEAMENTO TÉCNICO (CORRIGIDO)
# ==============================================================================

library(ipeadatar)
library(dplyr)
library(lubridate)

# 1. Coletar os índices (Mantendo o nome 'value' original para facilitar o join)
ipca_m <- ipeadata("PRECOS12_IPCA12") 
igp_m  <- ipeadata("IGP12_IGPDI12") 

# 2. Procedimento de Encadeamento (Linking)
ponto_uniao <- as.Date("1980-01-01")

# Extrair valores para o fator de escala
val_ipca_80 <- ipca_m$value[ipca_m$date == ponto_uniao]
val_igp_80  <- igp_m$value[igp_m$date == ponto_uniao]

fator_escala <- val_ipca_80 / val_igp_80

# 3. Criar a série contínua imputada
indices_ajustados <- igp_m %>%
  filter(date < ponto_uniao & date >= as.Date("1979-01-01")) %>%
  # Ajustamos o nível do IGP-DI para a escala do IPCA antes de unir
  mutate(value = value * fator_escala) %>%
  bind_rows(ipca_m %>% filter(date >= ponto_uniao)) %>%
  mutate(ano = year(date))

# 4. Calcular o Deflator Anual (Média Anual para variáveis de fluxo)
deflator_anual <- indices_ajustados %>%
  group_by(ano) %>%
  summarise(indice_medio = mean(value, na.rm = TRUE), .groups = 'drop')

# 5. Definir Ano de Referência e Fatores
ano_ref <- 2010
valor_base <- deflator_anual$indice_medio[deflator_anual$ano == ano_ref]

deflator_final <- deflator_anual %>%
  mutate(fator_deflator = valor_base / indice_medio)

# 6. Aplicação no seu Painel de Dados
# Usando across() para aplicar a mesma operação em várias colunas de uma vez
painel_final <- dados_BR_1980_2024 %>%
  left_join(deflator_final, by = "ano") %>%
  mutate(
    across(
      .cols = c(gasto_edu_nominal, gasto_saúde_nominal, 
                gasto_seg_nominal, gasto_prevassist_nominal, receita_nominal,
                despesa_nominal),
      .fns = ~ (.x * fator_deflator) / populacao,
      .names = "{.col}_per_capita"
    )
  ) %>%
  # Limpeza opcional: remover o termo "_nominal" dos novos nomes
  rename_with(~gsub("_nominal_per_capita", "_per_capita", .x), contains("_per_capita"))

# 7. Verificação da transição (1979-1981)
print(deflator_final %>% filter(ano %in% 1978:1982))

# Visualizar resultado das novas colunas
# Forma segura de visualizar sem erro de "argumentos não utilizados"
painel_final %>% 
  dplyr::select(ano, Estado, ends_with("_per_capita")) %>% 
  head()

