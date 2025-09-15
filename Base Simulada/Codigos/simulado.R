library(rstan)
library(dplyr)
library(survival)
library(ggsurvfit)
library(survminer)
library(ggplot2)
library(loo)
library(plotly)
library(purrr)

dirname(rstudioapi::getSourceEditorContext()$path) |> 
  setwd()

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

set.seed(123)

N <- 5000

# Parâmetros (iguais ao seu original)
true_beta <- c(-1.5, 1.5)  # Efeito de categoria_cliente na cura
true_alpha <- c(3.2, 1.3)   # Efeito de forma_cobranca no tempo
true_sigma <- 0.6

# Gerar a variável que realmente causa a cura.
categoria_cliente <- factor(sample(c("Cliente A", "Cliente B"), size = N, replace = TRUE))

# Gerar a forma_cobranca com uma probabilidade que DEPENDE da categoria_cliente.
# Vamos criar uma forte preferência:
# - Cliente A (baixa cura) vai preferir usar "Cartão".
# - Cliente B (alta cura) vai preferir usar "Dinheiro".

# Define a probabilidade de um cliente usar "Dinheiro"
# Se for Cliente B, a chance é de 80%. Se for Cliente A, a chance é de 20%.
prob_dinheiro <- ifelse(categoria_cliente == "Cliente B", 0.80, 0.20)

# Amostra a forma de cobrança para cada cliente com base na sua probabilidade
forma_cobranca <- factor(
  sapply(prob_dinheiro, function(p) {
    sample(c("Dinheiro", "Cartão"), size = 1, prob = c(p, 1 - p))
  }),
  levels = c("Cartão", "Dinheiro") # FAZ CARTÃO SER O INTERCEPTO
)

# Criando o dataframe
dados_simulados <- data.frame(
  id = 1:N,
  categoria_cliente = categoria_cliente,
  forma_cobranca = forma_cobranca
)

# Matrizes para o modelo (com intercepto)
Z <- model.matrix(~ categoria_cliente, data = dados_simulados)  # para fração de cura
X <- model.matrix(~ forma_cobranca, data = dados_simulados)    # para tempo (mu)

# Calculando a probabilidade de cura (p_i)
prob_cura <- boot::inv.logit(Z %*% true_beta)  

# Gerando status de cura:
# status_cura == 0 => curado
# status_cura == 1 => suscetível (pode ter o evento)
status_cura <- rbinom(N, 1, prob = 1 - prob_cura)

dados_simulados$status_real <- factor(ifelse(status_cura == 1, "Suscetível", "Curado"))

# Gerando tempo de evento para suscetíveis
mu_individual <- as.vector(X %*% true_alpha)
tempo_evento <- rep(Inf, N)  # curados não terão evento
idx_suscetiveis <- which(status_cura == 1)

tempo_evento[idx_suscetiveis] <- rlnorm(
  n = length(idx_suscetiveis),
  meanlog = mu_individual[idx_suscetiveis],
  sdlog = true_sigma
)

# Tempo de censura
tempo_censura <- sample(1:150, size = N, replace = TRUE)

# Tempo observado é o mínimo entre tempo de evento e censura
tempo_obs_continuo <- pmin(tempo_evento, tempo_censura)

# Arredondar tempo observado para inteiro (meses, por exemplo)
dados_simulados$tempo <- ceiling(tempo_obs_continuo)

# Indicador delta: 1 se evento observado, 0 se censurado (incluindo curados)
dados_simulados$delta <- as.numeric(tempo_evento <= tempo_censura)

# Visualização para conferir
print("Tabela de contingência entre categoria_cliente e forma_cobranca:")
print(table(dados_simulados$categoria_cliente, dados_simulados$forma_cobranca))

print("Amostra dos dados finais:")
print(head(dados_simulados))

# Salvar os dados simulados
saveRDS(dados_simulados, file = "../Base/dados_simulados.rds")

######################################################################################################
# Análise Exploratória
data <- readRDS("../Base/dados_simulados.rds")

# forma_cobranca
km_fit <- survfit2(Surv(tempo, delta) ~ forma_cobranca, data = data)
ggsurvplot(km_fit, data = data)

# categoria_cliente
km_fit <- survfit2(Surv(tempo, delta) ~ categoria_cliente, data = data)
ggsurvplot(km_fit, data = data)

#######################################################################################################
# Modelo Stan
modelo_lognormal_cura <- "
data {
  int<lower=1> N;
  vector<lower=0>[N] t;
  int<lower=0,upper=1> delta[N];
  
  int<lower=1> Q;
  int<lower=1> P;
  
  matrix[N, Q] Z;
  matrix[N, P] X;
}
parameters {
  vector[Q] beta;
  vector[P] alpha;
  real<lower=0> sigma2;
}
transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  vector[N] p = inv_logit(Z * beta);
  vector[N] mu = X * alpha;
}
model {
  beta ~ multi_normal(rep_vector(0, Q), diag_matrix(rep_vector(25, Q)));
  alpha ~ multi_normal(rep_vector(0, P), diag_matrix(rep_vector(25, P)));
  sigma2 ~ inv_gamma(0.01, 0.01);

  for (i in 1:N) {
    if (delta[i] == 1) {
      target += log1m(p[i]) + lognormal_lpdf(t[i] | mu[i], sigma);
    } else {
      target += log_mix(p[i],
                        0,
                        lognormal_lccdf(t[i] | mu[i], sigma));
    }
  }
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    if (delta[i] == 1) {
      log_lik[i] = log1m(p[i]) + lognormal_lpdf(t[i] | mu[i], sigma);
    } else {
      log_lik[i] = log_mix(p[i], 0, lognormal_lccdf(t[i] | mu[i], sigma));
    }
  }
}
"

data <- readRDS("../Base/dados_simulados.rds")

# Cenário 2: cura sendo intercepto
X <- model.matrix(~ forma_cobranca, data = data)
Z <- model.matrix(~ 1, data = data)

data_stan <- list(N=nrow(data),
                  t=data$tempo,
                  delta=data$delta, 
                  Q=ncol(Z), 
                  P=ncol(X), 
                  Z=Z,
                  X=X)

ajuste_2 <- stan(model_code=modelo_lognormal_cura, 
                 data=data_stan, 
                 iter = 2000, 
                 chains = 4, 
                 seed = 1, 
                 control = list(adapt_delta=0.99),
                 save_warmup = FALSE)

saveRDS(ajuste_2, file = "../Ajustes/ajuste_modelagem2_sim.rds")

# LOO-CV
loo_2 <- loo(ajuste_2, cores = 2)
saveRDS(loo_2, file = "../Criterios/loo_modelagem2_sim.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_2, pars = "log_lik")$log_lik
waic_2 <- waic(log_lik_array)
saveRDS(waic_2, file = "../Criterios/waic_modelagem2_sim.rds")

# Gráficos
data$grupo <- data$forma_cobranca

# Posterior
posterior <- extract(ajuste_2)
n_amostras <- length(posterior$sigma)
t_seq <- seq(1, max(data$tempo))

# Curvas modeladas com IC
df_model_all <- data.frame()
grupos_unicos <- unique(data$grupo)

for (g in grupos_unicos) {
  idx <- which(data$grupo == g)
  
  S_samples <- matrix(NA, nrow = n_amostras, ncol = length(t_seq))
  h_samples <- matrix(NA, nrow = n_amostras, ncol = length(t_seq))
  
  for (s in 1:n_amostras) {
    beta_s <- posterior$beta[s, ]
    alpha_s <- posterior$alpha[s, ]
    sigma_s <- posterior$sigma[s]
    
    Z_g <- as.matrix(Z[idx, , drop = FALSE])
    X_g <- as.matrix(X[idx, , drop = FALSE])
    
    p_i <- plogis(Z_g %*% beta_s)
    mu_i <- X_g %*% alpha_s
    
    # Sobrevivência e densidade para cada tempo
    S_matrix <- sapply(t_seq, function(t) {
      p_i + (1 - p_i) * (1 - plnorm(t, mu_i, sigma_s))
    })
    f_matrix <- sapply(t_seq, function(t) {
      (1 - p_i) * dlnorm(t, mu_i, sigma_s)
    })
    
    # Sobrevivência e risco médios do grupo
    S_i <- colMeans(S_matrix)
    f_i <- colMeans(f_matrix)
    h_i <- f_i / S_i
    
    S_samples[s, ] <- S_i
    h_samples[s, ] <- h_i
  }
  
  # Estatísticas da curva de sobrevivência
  S_median <- apply(S_samples, 2, median)
  S_lower <- apply(S_samples, 2, quantile, probs = 0.025)
  S_upper <- apply(S_samples, 2, quantile, probs = 0.975)
  
  # Estatísticas da hazard
  h_median <- apply(h_samples, 2, median)
  h_lower <- apply(h_samples, 2, quantile, probs = 0.025)
  h_upper <- apply(h_samples, 2, quantile, probs = 0.975)
  
  df_model_all <- rbind(df_model_all, data.frame(
    time = t_seq,
    survival = S_median,
    surv_lower = S_lower,
    surv_upper = S_upper,
    hazard = h_median,
    haz_lower = h_lower,
    haz_upper = h_upper,
    grupo = g
  ))
}

# Kaplan-Meier
df_km <- data.frame()

for (g in grupos_unicos) {
  subset_data <- data[data$grupo == g, ]
  surv_fit <- survfit(Surv(tempo, delta) ~ 1, data = subset_data)
  
  df_km <- rbind(df_km, data.frame(
    time = surv_fit$time,
    survival = surv_fit$surv,
    grupo = g
  ))
}

df_km <- df_km %>%
  group_by(grupo) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(haz = 1 - (survival / lag(survival)))

# Salvando
saveRDS(df_model_all, file = "../DadosGraficos/df_model_all_modelagem2_sim.rds")
saveRDS(df_km, file = "../DadosGraficos/df_km_sim.rds")

# Gráfico sobrevivência com IC
p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.6) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = surv_lower, ymax = surv_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = survival, color = grupo)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank())

# Gráfico hazard com IC
p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.6) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = haz_lower, ymax = haz_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = hazard, color = grupo)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank())

# Visualização
subplot(ggplotly(p1), plotly::style(ggplotly(p2), showlegend = FALSE))

# Sumário
print(ajuste_2, pars = c("beta", "alpha", "sigma"))
traceplot(ajuste_2, pars = c("beta", "alpha", "sigma"))

##############################################################################################################
# Cenário 3: cura sendo categoria_cliente
X <- model.matrix(~ forma_cobranca, data = data)
Z <- model.matrix(~ categoria_cliente, data = data)

data_stan <- list(N=nrow(data),
                  t=data$tempo,
                  delta=data$delta, 
                  Q=ncol(Z), 
                  P=ncol(X), 
                  Z=Z,
                  X=X)

ajuste_3 <- stan(model_code=modelo_lognormal_cura, 
                 data=data_stan, 
                 iter = 2000, 
                 chains = 4, 
                 seed = 1, 
                 control = list(adapt_delta=0.99),
                 save_warmup = FALSE)

saveRDS(ajuste_3, file = "../Ajustes/ajuste_modelagem3_sim.rds")

# LOO-CV
loo_3 <- loo(ajuste_3, cores = 2)
saveRDS(loo_3, file = "../Criterios/loo_modelagem3_sim.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_3, pars = "log_lik")$log_lik
waic_3 <- waic(log_lik_array)
saveRDS(waic_3, file = "../Criterios/waic_modelagem3_sim.rds")

# Gráficos
data$grupo <- data$forma_cobranca

# Posterior
posterior <- extract(ajuste_3)
n_amostras <- length(posterior$sigma)
t_seq <- seq(1, max(data$tempo))

# Curvas modeladas com IC
df_model_all <- data.frame()
grupos_unicos <- unique(data$grupo)

for (g in grupos_unicos) {
  idx <- which(data$grupo == g)
  
  S_samples <- matrix(NA, nrow = n_amostras, ncol = length(t_seq))
  h_samples <- matrix(NA, nrow = n_amostras, ncol = length(t_seq))
  
  for (s in 1:n_amostras) {
    beta_s <- posterior$beta[s, ]
    alpha_s <- posterior$alpha[s, ]
    sigma_s <- posterior$sigma[s]
    
    Z_g <- as.matrix(Z[idx, , drop = FALSE])
    X_g <- as.matrix(X[idx, , drop = FALSE])
    
    p_i <- plogis(Z_g %*% beta_s)
    mu_i <- X_g %*% alpha_s
    
    # Sobrevivência e densidade para cada tempo
    S_matrix <- sapply(t_seq, function(t) {
      p_i + (1 - p_i) * (1 - plnorm(t, mu_i, sigma_s))
    })
    f_matrix <- sapply(t_seq, function(t) {
      (1 - p_i) * dlnorm(t, mu_i, sigma_s)
    })
    
    # Sobrevivência e risco médios do grupo
    S_i <- colMeans(S_matrix)
    f_i <- colMeans(f_matrix)
    h_i <- f_i / S_i
    
    S_samples[s, ] <- S_i
    h_samples[s, ] <- h_i
  }
  
  # Estatísticas da curva de sobrevivência
  S_median <- apply(S_samples, 2, median)
  S_lower <- apply(S_samples, 2, quantile, probs = 0.025)
  S_upper <- apply(S_samples, 2, quantile, probs = 0.975)
  
  # Estatísticas da hazard
  h_median <- apply(h_samples, 2, median)
  h_lower <- apply(h_samples, 2, quantile, probs = 0.025)
  h_upper <- apply(h_samples, 2, quantile, probs = 0.975)
  
  df_model_all <- rbind(df_model_all, data.frame(
    time = t_seq,
    survival = S_median,
    surv_lower = S_lower,
    surv_upper = S_upper,
    hazard = h_median,
    haz_lower = h_lower,
    haz_upper = h_upper,
    grupo = g
  ))
}

# Kaplan-Meier
df_km <- data.frame()

for (g in grupos_unicos) {
  subset_data <- data[data$grupo == g, ]
  surv_fit <- survfit(Surv(tempo, delta) ~ 1, data = subset_data)
  
  df_km <- rbind(df_km, data.frame(
    time = surv_fit$time,
    survival = surv_fit$surv,
    grupo = g
  ))
}

df_km <- df_km %>%
  group_by(grupo) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(haz = 1 - (survival / lag(survival)))

# Salvando
saveRDS(df_model_all, file = "../DadosGraficos/df_model_all_modelagem3_sim.rds")

# Gráfico sobrevivência com IC
p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.6) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = surv_lower, ymax = surv_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = survival, color = grupo)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank())

# Gráfico hazard com IC
p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.6) +
  geom_ribbon(data = df_model_all, aes(x = time, ymin = haz_lower, ymax = haz_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_all, aes(x = time, y = hazard, color = grupo)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank())

# Visualização
subplot(ggplotly(p1), plotly::style(ggplotly(p2), showlegend = FALSE))

# Sumário
print(ajuste_3, pars = c("beta", "alpha", "sigma"))
traceplot(ajuste_3, pars = c("beta", "alpha", "sigma"))

#######################################################################################################
modelo_lognormal_sem_cura <- "
data {
  int<lower=1> N;
  vector<lower=0>[N] t;
  int<lower=0,upper=1> delta[N];
  
  int<lower=1> P;
  matrix[N, P] X;
}
parameters {
  vector[P] alpha;
  real<lower=0> sigma2;
}
transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  vector[N] mu = X * alpha;
}
model {
  alpha ~ multi_normal(rep_vector(0, P), diag_matrix(rep_vector(25, P)));
  sigma2 ~ inv_gamma(0.01, 0.01);
  
  for (i in 1:N) {
    if (delta[i] == 1) {
      // Evento: usa a densidade (lpdf)
      target += lognormal_lpdf(t[i] | mu[i], sigma);
    } else {
      // Censura: usa a função de sobrevivência (lccdf)
      target += lognormal_lccdf(t[i] | mu[i], sigma);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    if (delta[i] == 1) {
      log_lik[i] = lognormal_lpdf(t[i] | mu[i], sigma);
    } else {
      log_lik[i] = lognormal_lccdf(t[i] | mu[i], sigma);
    }
  }
}
"

# Cenário 1: Log-Normal Padrão (Sem Cura)
X <- model.matrix(~ forma_cobranca, data = data)

data_stan <- list(N=nrow(data),
                  t=data$tempo,
                  delta=data$delta,
                  P=ncol(X),
                  X=X)

ajuste_1 <- stan(model_code=modelo_lognormal_sem_cura, 
               data=data_stan, 
               iter = 2000, 
               chains = 4, 
               seed = 1, 
               control = list(adapt_delta=0.99),
               save_warmup = FALSE)

saveRDS(ajuste_1, file = "../Ajustes/ajuste_modelagem1_sim.rds")

# LOO-CV
loo_1 <- loo(ajuste_1, cores = 2)
saveRDS(loo_1, file = "../Criterios/loo_modelagem1_sim.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_1, pars = "log_lik")$log_lik
waic_1 <- waic(log_lik_array)
saveRDS(waic_1, file = "../Criterios/waic_modelagem1_sim.rds")

data$grupo <- data$forma_cobranca

# Posterior do modelo log-normal padrão (sem cura)
posterior_ln <- extract(ajuste_1)
n_amostras_ln <- length(posterior_ln$sigma)
t_seq <- seq(1, max(data$tempo))

# Curvas modeladas com IC
grupos_unicos <- unique(data$grupo)

# Data frame para armazenar resultados para todos os grupos
df_model_ln_all <- data.frame()

for (g in grupos_unicos) {
  idx <- which(data$grupo == g)
  
  # Matrizes para armazenar as amostras de sobrevivência e hazard
  S_samples_ln <- matrix(NA, nrow = n_amostras_ln, ncol = length(t_seq))
  h_samples_ln <- matrix(NA, nrow = n_amostras_ln, ncol = length(t_seq))
  
  for (s in 1:n_amostras_ln) {
    alpha_s <- posterior_ln$alpha[s, ]
    sigma_s <- posterior_ln$sigma[s]
    
    mu_i <- as.vector(X[idx, , drop = FALSE] %*% alpha_s)
    
    # Calcula S(t) e f(t) para cada tempo t em t_seq
    S_matrix_ln <- sapply(t_seq, function(t) {
      1 - plnorm(t, mu_i, sigma_s)
    })
    
    f_matrix_ln <- sapply(t_seq, function(t) {
      dlnorm(t, mu_i, sigma_s)
    })
    
    # Média dos indivíduos do grupo para cada t
    S_i <- colMeans(S_matrix_ln)
    f_i <- colMeans(f_matrix_ln)
    h_i <- f_i / S_i
    
    S_samples_ln[s, ] <- S_i
    h_samples_ln[s, ] <- h_i
  }
  
  # Estatísticas da curva de sobrevivência
  S_median <- apply(S_samples_ln, 2, median)
  S_lower <- apply(S_samples_ln, 2, quantile, probs = 0.025)
  S_upper <- apply(S_samples_ln, 2, quantile, probs = 0.975)
  
  # Estatísticas da hazard
  h_median <- apply(h_samples_ln, 2, median)
  h_lower <- apply(h_samples_ln, 2, quantile, probs = 0.025)
  h_upper <- apply(h_samples_ln, 2, quantile, probs = 0.975)
  
  df_model_ln_all <- rbind(df_model_ln_all, data.frame(
    time = t_seq,
    survival = S_median,
    surv_lower = S_lower,
    surv_upper = S_upper,
    hazard = h_median,
    haz_lower = h_lower,
    haz_upper = h_upper,
    grupo = g
  ))
}

# Kaplan-Meier
df_km <- data.frame()

for (g in grupos_unicos) {
  subset_data <- data[data$grupo == g, ]
  surv_fit <- survfit(Surv(tempo, delta) ~ 1, data = subset_data)
  
  df_km <- rbind(df_km, data.frame(
    time = surv_fit$time,
    survival = surv_fit$surv,
    grupo = g
  ))
}

df_km <- df_km %>%
  group_by(grupo) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(haz = 1 - (survival / lag(survival)))

# Salvando
saveRDS(df_model_ln_all, file = "../DadosGraficos/df_model_ln_all_modelagem1_sim.rds")

# Sobrevivência
p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.7) +
  geom_line(data = df_model_ln_all, aes(x = time, y = survival, color = grupo), size = 0.8) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_d() +
  labs(y = "S(t)", x = "Tempo", color = "Grupo") +
  theme(legend.title = element_text())

# Risco
p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.7) +
  geom_line(data = df_model_ln_all, aes(x = time, y = hazard, color = grupo), size = 0.8) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.1)) +
  scale_color_viridis_d() +
  labs(y = "h(t)", x = "Tempo", color = "Grupo") +
  theme(legend.title = element_text())

# Visualização
subplot(ggplotly(p1), plotly::style(ggplotly(p2), showlegend = FALSE))

print(ajuste_1, pars = c("alpha", "sigma"))
traceplot(ajuste_1, pars = c("alpha", "sigma"))