# Calculate log-RSS given model formula, betas, x1, and x2

my_log_rss <- function(model, formula, betas, x1, x2, ci_level){
  
  # Difference between x1 and x2
  diff_x <- sweep(model.matrix(formula, x1), 
                  2, 
                  model.matrix(formula, x2))
  
  # Get mean estimate
  log_rss <- diff_x %*% betas
  
  # Get variance
  var_pred <- diag(diff_x %*% vcov(model) %*% t(diff_x))
  # Get SE
  se <- unname(sqrt(var_pred))
  # Get critical value 
  p <- 1 - ((1 - ci_level)/2)
  zstar <- qnorm(p)
  # Compute bounds
  log_rss_lwr <- log_rss - zstar * se
  log_rss_upr <- log_rss + zstar * se
  
  # Bind
  out <- data.frame(log_rss,
                    log_rss_lwr, 
                    log_rss_upr)
  
  # Return
  return(out)
  
}

# Empirical bootstrapping
# To generate confidence intervals, sample from individual model betas with 
# replacement, refit IVW regression to get mean parameter estimate, and
# re-calculate log-RSS.

boot <- function(dat, issa_formula, x1, x2) {
  
  # Number of individuals
  n_indiv <- nrow(dat$data[[1]])
  
  # Resample the rows (individuals)
  rows <- sample(1:n_indiv, size = n_indiv, replace = TRUE)
  iw <- dat %>% 
    mutate(new = lapply(data, function(x){
      return(x[rows, ])
    })) %>% 
    # Fit the model 
    mutate(iw = lapply(new, function(x) {
      
      mod <- lm(estimate ~ 1, data = x, weights = 1/(std.error)^2)
      
      return(mod)
    })) %>% 
    # Predict the betas
    mutate(pred = lapply(iw, function(x) {
      
      pred <- predict(x, 
                      newdata = data.frame(dummy = NA),
                      se.fit = TRUE)
      
      est <- data.frame(mean = pred$fit,
                        lwr = pred$fit - 1.96 * pred$se.fit,
                        upr = pred$fit + 1.96 * pred$se.fit)
      
      return(est)
    })) %>% 
    # Unnest
    select(term, pred) %>% 
    unnest(cols = pred)
  
  # Calculate log-RSS
  res <- my_log_rss(formula =  issa_formula,
                    betas = iw$mean,
                    x1 = x1,
                    x2 = x2)
  
  # Return
  return(res)
}