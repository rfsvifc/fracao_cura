# Pacotes
library(rstan)
library(dplyr)
library(CASdatasets)
library(survival)
library(ggplot2)
library(plotly)
library(loo)
library(purrr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Dados
data("uslapseagent")
data <- uslapseagent

data <- data %>%
  mutate(duration_modificada = floor(duration + 0.5)) %>%
  filter(duration_modificada > 0 & duration_modificada <= 60)

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

# === MODELAGEM 2: Cura Constante ===
X2 <- model.matrix(~ acc.death.rider + risk.state, data = data)
Z2 <- model.matrix(~ 1, data = data)

data_stan2 <- list(N = nrow(data), 
                   t = data$duration_modificada, 
                   delta = data$surrender, 
                   Q = ncol(Z2), 
                   P = ncol(X2), 
                   Z = Z2, 
                   X = X2)

ajuste_2 <- stan(model_code = modelo_lognormal_cura, 
                 data = data_stan2, 
                 iter = 2000, 
                 chains = 4, 
                 seed = 1, 
                 control = list(adapt_delta=0.99),
                 save_warmup = FALSE)

saveRDS(ajuste_2, file = "../Ajustes/ajuste_modelagem2.rds")

# LOO-CV
loo_2 <- loo(ajuste_2, cores = 2)
saveRDS(loo_2, file = "../Criterios/loo_modelagem2.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_2, pars = "log_lik")$log_lik
waic_2 <- waic(log_lik_array)
saveRDS(waic_2, file = "../Criterios/waic_modelagem2.rds")

# === MODELAGEM 3 ===
X3 <- model.matrix(~ acc.death.rider + risk.state, data = data)
Z3 <- model.matrix(~ premium.frequency, data = data)

data_stan3 <- list(N = nrow(data), 
                   t = data$duration_modificada, 
                   delta = data$surrender, 
                   Q = ncol(Z3), 
                   P = ncol(X3), 
                   Z = Z3, 
                   X = X3)

ajuste_3 <- stan(model_code = modelo_lognormal_cura,
                 data = data_stan3,
                 iter = 2000,
                 chains = 4,
                 seed = 1, 
                 control = list(adapt_delta=0.99),
                 save_warmup = FALSE)

saveRDS(ajuste_3, file = "../Ajustes/ajuste_modelagem3.rds")

# LOO-CV
loo_3 <- loo(ajuste_3, cores = 2)
saveRDS(loo_3, file = "../Criterios/loo_modelagem3.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_3, pars = "log_lik")$log_lik
waic_3 <- waic(log_lik_array)
saveRDS(waic_3, file = "../Criterios/waic_modelagem3.rds")

# === CENÁRIO 4 ===
X4 <- model.matrix(~ acc.death.rider + risk.state, data = data)
Z4 <- model.matrix(~ premium.frequency + underwriting.age, data = data)

data_stan4 <- list(N = nrow(data), 
                   t = data$duration_modificada,
                   delta = data$surrender,
                   Q = ncol(Z4),
                   P = ncol(X4),
                   Z = Z4,
                   X = X4)

ajuste_4 <- stan(model_code = modelo_lognormal_cura, 
                 data = data_stan4, 
                 iter = 2000, 
                 chains = 4, 
                 seed = 1, 
                 control = list(adapt_delta=0.99),
                 save_warmup = FALSE)

saveRDS(ajuste_4, file = "../Ajustes/ajuste_modelagem4.rds")

# LOO-CV
loo_4 <- loo(ajuste_4, cores = 2)
saveRDS(loo_4, file = "../Criterios/loo_modelagem4.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_4, pars = "log_lik")$log_lik
waic_4 <- waic(log_lik_array)
saveRDS(waic_4, file = "../Criterios/waic_modelagem4.rds")

######################################################################################
# Gráficos
ajuste <- readRDS("../Ajustes/ajuste_modelagem3.rds")

data$grupo <- interaction(data$acc.death.rider, data$risk.state, drop = TRUE)

X <- model.matrix(~ acc.death.rider + risk.state, data = data)
Z <- model.matrix(~ premium.frequency, data = data)

# Posterior
posterior <- extract(ajuste)
n_amostras <- length(posterior$sigma)
t_seq <- seq(1, max(data$duration_modificada))

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
  surv_fit <- survfit(Surv(duration_modificada, surrender) ~ 1, data = subset_data)
  
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
saveRDS(df_model_all, file = "../DadosGraficos/df_model_all_modelagem3.rds")

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

###########################################################################################
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

# === MODELAGEM 1: Log-Normal Padrão (Sem Cura) ===
X1 <- model.matrix(~ acc.death.rider + risk.state, data = data)

data_stan1 <- list(N = nrow(data), 
                   t = data$duration_modificada,
                   delta = data$surrender, 
                   P = ncol(X1), 
                   X = X1)

ajuste_1 <- stan(model_code = modelo_lognormal_sem_cura, 
                 data = data_stan1, 
                 iter = 2000, 
                 chains = 4, 
                 seed = 1, 
                 control = list(adapt_delta=0.99),
                 save_warmup = FALSE)

saveRDS(ajuste_1, file = "../Ajustes/ajuste_modelagem1.rds")

# LOO-CV
loo_1 <- loo(ajuste_1, cores = 2)
saveRDS(loo_1, file = "../Criterios/loo_modelagem1.rds")

# WAIC
log_lik_array <- rstan::extract(ajuste_1, pars = "log_lik")$log_lik
waic_1 <- waic(log_lik_array)
saveRDS(waic_1, file = "../Criterios/waic_modelagem1.rds")

####################################################################################################
# Comparação Final dos Modelos
loo_1 <- readRDS("../Criterios/loo_modelagem1.rds")
loo_2 <- readRDS("../Criterios/loo_modelagem2.rds")
loo_3 <- readRDS("../Criterios/loo_modelagem3.rds")
loo_4 <- readRDS("../Criterios/loo_modelagem4.rds")

comp <- loo_compare(loo_1, loo_2, loo_3, loo_4)
print(comp)

####################################################################################################
ajuste_ln <- readRDS("../Ajustes/ajuste_modelagem1.rds")

data$grupo <- interaction(data$acc.death.rider, data$risk.state, drop = TRUE)

X <- model.matrix(~ acc.death.rider + risk.state, data = data)

# Posterior do modelo log-normal padrão (sem cura)
posterior_ln <- extract(ajuste_ln)
n_amostras_ln <- length(posterior_ln$sigma)
t_seq <- seq(1, max(data$duration_modificada))

# Curvas modeladas com IC
df_model_ln_all <- data.frame()
grupos_unicos <- unique(data$grupo)

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
  surv_fit <- survfit(Surv(duration_modificada, surrender) ~ 1, data = subset_data)
  
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
saveRDS(df_model_ln_all, file = "../DadosGraficos/df_model_ln_all_modelagem1.rds")

# Gráfico sobrevivência com IC
p1 <- ggplot() +
  geom_step(data = df_km, aes(x = time, y = survival, color = grupo), alpha = 0.6) +
  geom_ribbon(data = df_model_ln_all, aes(x = time, ymin = surv_lower, ymax = surv_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_ln_all, aes(x = time, y = survival, color = grupo)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank())

# Gráfico hazard com IC
p2 <- ggplot() +
  geom_line(data = df_km, aes(x = time, y = haz, color = grupo), alpha = 0.6) +
  geom_ribbon(data = df_model_ln_all, aes(x = time, ymin = haz_lower, ymax = haz_upper, fill = grupo), alpha = 0.2) +
  geom_line(data = df_model_ln_all, aes(x = time, y = hazard, color = grupo)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank())

# Visualização
subplot(ggplotly(p1), plotly::style(ggplotly(p2), showlegend = FALSE))

print(ajuste_ln, pars = c("alpha", "sigma"))
traceplot(ajuste_ln, pars = c("alpha", "sigma"))