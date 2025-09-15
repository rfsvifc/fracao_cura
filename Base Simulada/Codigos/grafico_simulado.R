library(dplyr)
library(survival)
library(ggplot2)
library(plotly)
library(CASdatasets)
library(viridis)
library(readr)
library(ggpubr)
library(rstan)
library(loo)
library(purrr)

# Definindo diretório de trabalho na pasta raiz do script #
dirname(rstudioapi::getSourceEditorContext()$path) |> 
  setwd()

data <- readRDS("../Base/dados_simulados.rds")

################################################################################################################
# Cenário 1 - Sem Cura
ajuste <- readRDS("../Ajustes/ajuste_modelagem1_sim.rds")
df_model_all <- readRDS("../DadosGraficos/df_model_ln_all_modelagem1_sim.rds")
df_km <- readRDS("../DadosGraficos/df_km_sim.rds")

p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.6, linewidth = 1) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = surv_lower, ymax = surv_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = survival, color = grupo), linewidth = 1) +
  theme_bw() +
  labs(x = "Tempo", y = "S(t)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.6, linewidth = 1) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = haz_lower, ymax = haz_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = hazard, color = grupo), linewidth = 1) +
  theme_bw() +
  labs(x = "Tempo", y = "h(t)") +
  scale_y_continuous(limits = c(0, 0.05)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

# Coeficientes
print(ajuste, pars = c("alpha", "sigma"))
traceplot(ajuste, pars = c("alpha", "sigma"))

# Comparação modelo

# LOO-CV
loo_1 <- readRDS("../Criterios/loo_modelagem1_sim.rds")
loo_1

# WAIC
waic_1 <- readRDS("../Criterios/waic_modelagem1_sim.rds")
waic_1

################################################################################################################
# Cenário 2 - Intercepto na Cura
ajuste <- readRDS("../Ajustes/ajuste_modelagem2_sim.rds")
df_model_all <- readRDS("../DadosGraficos/df_model_all_modelagem2_sim.rds")
df_km <- readRDS("../DadosGraficos/df_km_sim.rds")

p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.6, linewidth = 1) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = surv_lower, ymax = surv_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = survival, color = grupo), linewidth = 1) +
  theme_bw() +
  labs(x = "Tempo", y = "S(t)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.6, linewidth = 1) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = haz_lower, ymax = haz_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = hazard, color = grupo), linewidth = 1) +
  theme_bw() +
  labs(x = "Tempo", y = "h(t)") +
  scale_y_continuous(limits = c(0, 0.05)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

# Coeficientes
print(ajuste, pars = c("beta", "alpha", "sigma"))
traceplot(ajuste, pars = c("beta", "alpha", "sigma"))

# Probabilidade de cura
p_amostras <- rstan::extract(ajuste, pars = "p")$p

# Como só há intercepto, todos os valores em cada linha são iguais
# Vamos pegar a coluna 1, pois representa toda a linha
p_iter <- p_amostras[, 1]  # vetor com valores de p por iteração

# Estatísticas resumidas
resumo_IC <- data.frame(
  grupo = "Global",
  p_medio = mean(p_iter),
  p_sd = sd(p_iter),
  p_inf = quantile(p_iter, 0.025),
  p_sup = quantile(p_iter, 0.975)
)

print(resumo_IC)

# Comparação modelo

# LOO-CV
loo_2 <- readRDS("../Criterios/loo_modelagem2_sim.rds")
loo_2

# WAIC
waic_2 <- readRDS("../Criterios/waic_modelagem2_sim.rds")
waic_2

################################################################################################################
# Cenário 3 - categoria_cliente na Cura
ajuste <- readRDS("../Ajustes/ajuste_modelagem3_sim.rds")
df_model_all <- readRDS("../DadosGraficos/df_model_all_modelagem3_sim.rds")
df_km <- readRDS("../DadosGraficos/df_km_sim.rds")

p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.6, linewidth = 1) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = surv_lower, ymax = surv_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = survival, color = grupo), linewidth = 1) +
  theme_bw() +
  labs(x = "Tempo", y = "S(t)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.6, linewidth = 1) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = haz_lower, ymax = haz_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = hazard, color = grupo), linewidth = 1) +
  theme_bw() +
  labs(x = "Tempo", y = "h(t)") +
  scale_y_continuous(limits = c(0, 0.05)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"))

ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

# Coeficientes
print(ajuste, pars = c("beta", "alpha", "sigma"))
traceplot(ajuste, pars = c("beta", "alpha", "sigma"))

# Probabilidade de Cura - Gráfico
p_amostras <- rstan::extract(ajuste, pars = "p")$p
p_medio <- apply(p_amostras, 2, mean)
p_ic <- apply(p_amostras, 2, quantile, probs = c(0.025, 0.975))

# Número de iterações e indivíduos
n_iter <- nrow(p_amostras)
n_ind <- ncol(p_amostras)

data$grupo <- data$forma_cobranca

# Criar um data frame "long" com o grupo de cada indivíduo
df_grupo <- data.frame(id = 1:n_ind, grupo = data$grupo)

# Para cada iteração, calcular média do p por grupo
# Inicializa uma lista para guardar médias por grupo e iteracao
medias_iter_grupo <- vector("list", n_iter)

for(i in 1:n_iter){
  p_i <- p_amostras[i, ]  # valores p da iteração i
  df_tmp <- data.frame(id = 1:n_ind, p = p_i)
  df_tmp <- left_join(df_tmp, df_grupo, by = "id")
  
  medias_iter_grupo[[i]] <- df_tmp %>%
    group_by(grupo) %>%
    summarise(p_medio = mean(p))
}

# Combina tudo num data frame
df_medias <- bind_rows(medias_iter_grupo, .id = "iter")
df_medias$iter <- as.integer(df_medias$iter)

resumo_IC <- df_medias %>%
  group_by(grupo) %>%
  summarise(p_vec = list(p_medio)) %>%
  mutate(
    p_medio = map_dbl(p_vec, mean),
    p_sd    = map_dbl(p_vec, sd),
    p_inf = map_dbl(p_vec, ~ quantile(.x, 0.025)),
    p_sup = map_dbl(p_vec, ~ quantile(.x, 0.975))
  ) %>%
  select(-p_vec)

print(resumo_IC)

# Comparação modelo

# LOO-CV
loo_3 <- readRDS("../Criterios/loo_modelagem3_sim.rds")
loo_3

# WAIC
waic_3 <- readRDS("../Criterios/waic_modelagem3_sim.rds")
waic_3