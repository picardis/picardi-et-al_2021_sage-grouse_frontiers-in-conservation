# # # # North Dakota Sage-Grouse Translocation # # # #
# Temporal dynamics of post-release habitat selection #
# # # # # # # Frontiers in Conservation # # # # # # # 
# # # Simona Picardi, Nathan Ranc, Brian Smith # # # 

# This is the code associated with manuscript:
# "Individual variation in temporal dynamics of post-release habitat selection"
# Picardi et al., Frontiers in Conservation Science

# Load packages ####

library(tidyverse)
library(amt)
library(MuMIn)
library(patchwork)

# Load data ####

dat <- read.csv("Picardi-et-al_Frontiers-in-Conservation_Data.csv")

nd_amt <- dat %>% 
  nest(rsteps = -hen_name)

# Model with time dependence : herbaceous ####

issa_t1h <- nd_amt %>% 
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- NA
    res <- try(fit_issf(data = x, 
                        formula = case_ ~ 
                          perennial_herb_scaled : log_dst +
                          dist_to_roads_log_scaled +
                          perennial_herb_scaled +
                          sagebrush_scaled +
                          slope_scaled +
                          sl_ + 
                          log_sl_ + 
                          cos_ta_ +
                          strata(step_id_),
                        model = TRUE))
  }))

# Model with no time dependence ####

issa_t0 <- nd_amt %>% 
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- NA
    res <- try(fit_issf(data = x, 
                        formula = case_ ~ 
                          dist_to_roads_log_scaled +
                          perennial_herb_scaled +
                          sagebrush_scaled +
                          slope_scaled +
                          sl_ + 
                          log_sl_ + 
                          cos_ta_ +
                          strata(step_id_),
                        model = TRUE))
  }))

# Model selection tables ####

mods <- issa_t0 %>% 
  dplyr::select(hen_name, t0 = issa) %>% 
  left_join(issa_t1h, by = "hen_name") %>% 
  dplyr::select(hen_name, t0, t1h = issa) %>% 
  pivot_longer(cols = t0:t1h, names_to = "model", values_to = "issa") %>% 
  nest(models = -hen_name) 

ms_ind <- lapply(mods$models, function(x) {
  mod_set <- lapply(x$issa, function(x) {
    return(x$model)
  })
  mod_sel <- as.data.frame(model.sel(mod_set))
  mod_sel$model_name <- ifelse(!is.na(mod_sel$`log_dst:perennial_herb_scaled`), "T1H", "T0")
  return(mod_sel)
})

names(ms_ind) <- mods$hen_name

# What is the top-ranked model for each individual (delta >2)?

topmod_delta2 <- data.frame(hen_name = names(ms_ind),
                            top = unlist(lapply(ms_ind, function(x) {
                              if (x$delta[2] > 2) {
                                return(x$model_name[1]) 
                              } else {return(NA)}
                            })),
                            row.names = NULL)

topmod_delta2 %>% 
  group_by(top) %>% 
  tally()

# Individuals with model t0 as top-ranked (delta > 2)

hens_t0 <- topmod_delta2 %>% 
  filter(top == "T0") %>% 
  pull(hen_name)

# Individuals with model t1h as top-ranked (delta > 2)

hens_t1h <- topmod_delta2 %>% 
  filter(top == "T1H") %>% 
  pull(hen_name)

# Correlations with status ####

status <- nd_amt %>% unnest(rsteps) %>% dplyr::select(hen_name, status) %>% unique()

topmod_delta2 <- topmod_delta2 %>% 
  left_join(status, by = "hen_name") 

# Status

counts <- topmod_delta2 %>% 
  filter(!is.na(top)) %>% 
  group_by(status, top) %>% 
  tally() %>% 
  pivot_wider(names_from = top, values_from = n) 

count_mat <- unname(as.matrix(counts))[, 2:3]
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0

fisher.test(count_mat)

# Individual parameter estimates ####

# Model T1 for individuals that had T1 as top-ranked (delta >2)
t1_t1t <- issa_t1h %>% 
  filter(hen_name %in% hens_t1h) 

# Model T0 for individuals that had T0 as top-ranked (delta >2)
t0_t0t <- issa_t0 %>% 
  filter(hen_name %in% hens_t0) 

# Model T0 for individuals that had T1 as top-ranked (delta >2)
t0_t1t <- issa_t0 %>% 
  filter(hen_name %in% hens_t1h) 

# Model T1 for individuals that had T0 as top-ranked (delta >2)
t1_t0t <- issa_t1h %>% 
  filter(hen_name %in% hens_t0)

# Bind model results for each model/group combo

t1_t1t <- bind_rows(lapply(t1_t1t$issa, function(x) {
  broom::tidy(x$model)
}), .id = "hen_name")

t0_t0t <- bind_rows(lapply(t0_t0t$issa, function(x) {
  broom::tidy(x$model)
}), .id = "hen_name")

t1_t0t <- bind_rows(lapply(t1_t0t$issa, function(x) {
  broom::tidy(x$model)
}), .id = "hen_name")

t0_t1t <- bind_rows(lapply(t0_t1t$issa, function(x) {
  broom::tidy(x$model)
}), .id = "hen_name")

# Calculate CIs at the individual level

t1_t1top <- t1_t1t %>% 
  mutate(lwr = estimate - 1.96 * std.error,
         upr = estimate + 1.96 * std.error) %>% 
  filter(!term %in% c("sl_", "cos_ta_", "log_sl_"))

t0_t0top <- t0_t0t %>% 
  mutate(lwr = estimate - 1.96 * std.error,
         upr = estimate + 1.96 * std.error) %>% 
  filter(!term %in% c("sl_", "cos_ta_", "log_sl_")) %>% 
  mutate(hen_name = as.character(
    as.numeric(hen_name) + length(unique(t1_t1top$hen_name))))

t1_t0top <- t1_t0t %>% 
  mutate(lwr = estimate - 1.96 * std.error,
         upr = estimate + 1.96 * std.error) %>% 
  filter(!term %in% c("sl_", "cos_ta_", "log_sl_")) %>% 
  mutate(hen_name = as.character(
    as.numeric(hen_name) + length(unique(t1_t1top$hen_name))))

t0_t1top <- t0_t1t %>% 
  mutate(lwr = estimate - 1.96 * std.error,
         upr = estimate + 1.96 * std.error) %>% 
  filter(!term %in% c("sl_", "cos_ta_", "log_sl_"))

# Make 4-panel plot of individual-level estimates:

pd <- position_dodge(0.75)

t1_t1top$term <- factor(t1_t1top$term)
t0_t0top$term <- factor(t0_t0top$term, levels = levels(t1_t1top$term))
t1_t0top$term <- factor(t1_t0top$term, levels = levels(t1_t1top$term))
t0_t1top$term <- factor(t0_t1top$term, levels = levels(t1_t1top$term))

colorz <- RColorBrewer::brewer.pal(n = length(levels(t1_t1top$term)), 
                                   name = "Dark2")
names(colorz) <- levels(t1_t1top$term)

p_t1_t1top <- ggplot(t1_t1top, aes(y = estimate, x = term,
                                   ymin = lwr, ymax = upr,
                                   group = factor(hen_name),
                                   color = term,
                                   fill = term
)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0.2) +
  geom_label(position = pd, size = 3, label = LETTERS[as.numeric(t1_t1top$hen_name)],
             color = "white") +
  theme_bw() +
  labs(x = " ", y = "log-RSS", color = " ", title = "Model T1, Group T1") +
  scale_color_manual(values = colorz) +
  scale_fill_manual(values = colorz) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape(guide = "none") +
  coord_cartesian(ylim = c(-4, 4))

p_t0_t0top <- ggplot(t0_t0top, aes(y = estimate, x = term,
                                   ymin = lwr, ymax = upr,
                                   group = factor(hen_name, levels = 6:12),
                                   color = term,
                                   fill = term
)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0.2) +
  geom_label(position = pd, size = 3, label = LETTERS[as.numeric(t0_t0top$hen_name)],
             color = "white") +
  theme_bw() +
  labs(x = " ", y = "log-RSS", color = " ", title = "Model T0, Group T0") +
  scale_color_manual(values = colorz) +
  scale_fill_manual(values = colorz) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape(guide = "none") +
  coord_cartesian(ylim = c(-4, 4))

p_t1_t0top <- ggplot(t1_t0top, aes(y = estimate, x = term,
                                   ymin = lwr, ymax = upr,
                                   group = factor(hen_name, levels = 6:12),
                                   color = term,
                                   fill = term
)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0.2) +
  geom_label(position = pd, size = 3, label = LETTERS[as.numeric(t1_t0top$hen_name)],
             color = "white") +
  theme_bw() +
  labs(x = " ", y = "log-RSS", color = " ", title = "Model T1, Group T0") +
  scale_color_manual(values = colorz) +
  scale_fill_manual(values = colorz) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape(guide = "none") +
  coord_cartesian(ylim = c(-4, 4))

p_t0_t1top <- ggplot(t0_t1top, aes(y = estimate, x = term,
                                   ymin = lwr, ymax = upr,
                                   group = factor(hen_name),
                                   color = term,
                                   fill = term
)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0.2) +
  geom_label(position = pd, size = 3, label = LETTERS[as.numeric(t0_t1top$hen_name)],
             color = "white") +
  theme_bw() +
  labs(x = " ", y = "log-RSS", color = " ", title = "Model T0, Group T1") +
  scale_color_manual(values = colorz) +
  scale_fill_manual(values = colorz) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape(guide = "none") +
  coord_cartesian(ylim = c(-4, 4))

((p_t0_t0top / p_t1_t0top) | (p_t1_t1top / p_t0_t1top))

# Plot model predictions through time ####

source("log-rss-boot.R")

unnested <- nd_amt %>% 
  dplyr::select(hen_name, rsteps) %>% 
  unnest(cols = rsteps)

# Compare median to 3rd quartile (corresponds to 46% vs 52% herbaceous)

x1 <- data.frame(log_dst = log(1:60),
                 perennial_herb_scaled = quantile(unnested$perennial_herb_scaled, .75), 
                 dist_to_roads_log_scaled = median(unnested$dist_to_roads_log_scaled),
                 sagebrush_scaled = median(unnested$sagebrush_scaled),
                 slope_scaled = median(unnested$slope_scaled),
                 sl_ = median(unnested$sl_),
                 log_sl_ = median(unnested$log_sl_),
                 cos_ta_ = median(unnested$cos_ta_))

x2 <- data.frame(log_dst = log(1),
                 perennial_herb_scaled = median(unnested$perennial_herb_scaled),
                 dist_to_roads_log_scaled = median(unnested$dist_to_roads_log_scaled),
                 sagebrush_scaled = median(unnested$sagebrush_scaled),
                 slope_scaled = median(unnested$slope_scaled),
                 sl_ = median(unnested$sl_),
                 log_sl_ = median(unnested$log_sl_),
                 cos_ta_ = median(unnested$cos_ta_))

# Compute log-RSS at the individual level

# Group T1
res_t1 <- data.frame()

for (i in 1:length(unique(t1_t1t$hen_name))) {
  
  who <- unique(t1_t1t$hen_name)[i]
  
  mod <- issa_t1h %>% 
    filter(hen_name %in% hens_t1h) %>% 
    slice(i) 
  
  mod <- mod$issa[[1]]$model
  
  terms <- data.frame(term = c("perennial_herb_scaled",
                               "dist_to_roads_log_scaled",
                               "sagebrush_scaled",
                               "slope_scaled",
                               "sl_",
                               "log_sl_",
                               "cos_ta_",
                               "perennial_herb_scaled:log_dst"))
  
  betas <- t1_t1t %>% 
    filter(hen_name == i) %>% 
    dplyr::select(term, estimate)
  
  betas <- terms %>% left_join(betas, by = "term")
  
  logRSS <- my_log_rss(formula = ~ 0 + 
                         perennial_herb_scaled +
                         perennial_herb_scaled : log_dst +
                         dist_to_roads_log_scaled +
                         sagebrush_scaled +
                         slope_scaled +
                         sl_ + 
                         log_sl_ + 
                         cos_ta_,
                       betas = betas$estimate,
                       x1 = x1,
                       x2 = x2,
                       model = mod,
                       ci_level = 0.90)
  
  logRSS <- x1 %>% 
    bind_cols(logRSS) %>% 
    mutate(hen_name = who)
  
  res_t1 <- bind_rows(res_t1, logRSS)
  
}

res_t1$group <- "Group T1"

# Group T0
res_t0 <- data.frame()

for (i in 1:length(unique(t1_t0t$hen_name))) {
  
  who <- unique(t1_t0t$hen_name)[i]
  
  mod <- issa_t1h %>% 
    filter(hen_name %in% hens_t0) %>% 
    slice(i) 
  
  mod <- mod$issa[[1]]$model
  
  terms <- data.frame(term = c("perennial_herb_scaled",
                               "dist_to_roads_log_scaled",
                               "sagebrush_scaled",
                               "slope_scaled",
                               "sl_",
                               "log_sl_",
                               "cos_ta_",
                               "perennial_herb_scaled:log_dst"))
  
  betas <- t1_t0t %>% 
    filter(hen_name == i) %>% 
    dplyr::select(term, estimate)
  
  betas <- terms %>% left_join(betas, by = "term")
  
  logRSS <- my_log_rss(formula = ~ 0 + 
                         perennial_herb_scaled +
                         perennial_herb_scaled : log_dst +
                         dist_to_roads_log_scaled +
                         sagebrush_scaled +
                         slope_scaled +
                         sl_ + 
                         log_sl_ + 
                         cos_ta_,
                       betas = betas$estimate,
                       x1 = x1,
                       x2 = x2,
                       model = mod,
                       ci_level = 0.90)
  
  logRSS <- x1 %>% 
    bind_cols(logRSS) %>% 
    mutate(hen_name = who)
  
  res_t0 <- bind_rows(res_t0, logRSS)
  
}

res_t0$group <- "Group T0"

# Hens with no top-ranked model (delta <=2)
res_noT <- data.frame()

t1_noT <- issa_t1h %>% 
  filter(!hen_name %in% c(hens_t1h, hens_t0))

t1_noT_tidy <- bind_rows(lapply(t1_noT$issa, function(x) {
  broom::tidy(x$model)
}), .id = "hen_name")

for (i in 1:length(unique(t1_noT$hen_name))) {
  
  who <- unique(t1_noT_tidy$hen_name)[i]
  
  mod <- issa_t1h %>% 
    filter(!hen_name %in% c(hens_t0, hens_t1h)) %>% 
    slice(i) 
  
  mod <- mod$issa[[1]]$model
  
  terms <- data.frame(term = c("perennial_herb_scaled",
                               "dist_to_roads_log_scaled",
                               "sagebrush_scaled",
                               "slope_scaled",
                               "sl_",
                               "log_sl_",
                               "cos_ta_",
                               "perennial_herb_scaled:log_dst"))
  
  betas <- t1_noT_tidy %>% 
    filter(hen_name == i) %>% 
    dplyr::select(term, estimate)
  
  betas <- terms %>% left_join(betas, by = "term")
  
  logRSS <- my_log_rss(formula = ~ 0 + 
                         perennial_herb_scaled +
                         perennial_herb_scaled : log_dst +
                         dist_to_roads_log_scaled +
                         sagebrush_scaled +
                         slope_scaled +
                         sl_ + 
                         log_sl_ + 
                         cos_ta_,
                       betas = betas$estimate,
                       x1 = x1,
                       x2 = x2,
                       model = mod,
                       ci_level = 0.90)
  
  logRSS <- x1 %>% 
    bind_cols(logRSS) %>% 
    mutate(hen_name = who)
  
  res_noT <- bind_rows(res_noT, logRSS)
  
}

res_noT$group <- "No group"

# Plot
res <- bind_rows(res_t1, res_t0, res_noT) %>% 
  mutate(group = factor(group, levels = c("Group T1", "Group T0", "No group")))

ggplot(res, aes(x = exp(log_dst), y = log_rss, group = hen_name, 
                color = hen_name)) +
  facet_wrap(~ group) +
  geom_path(size = 1) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() +
  labs(x = "Days since translocation", y = "log-RSS") +
  theme(legend.position = "none") 
