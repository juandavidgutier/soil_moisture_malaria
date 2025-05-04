library(tidyr); library(dlnm) ; library(mixmeta) ; library(splines) ; library(tibble)
library(dplyr); library(ggplot2); library(RCurl); library(MKdescr); library(purrr); library(broom)

df_na <- read.csv("D:/top50.csv")
df <- df_na %>% drop_na()

# sd units
df$Rain <- zscore(df$Rain, na.rm = TRUE)
df$Temperature <- zscore(df$Temperature, na.rm = TRUE)
df$SST12 <- zscore(df$SST12, na.rm = TRUE)
df$TROP <- zscore(df$TROP, na.rm = TRUE)
df$HFP <- as.numeric(df$HFP)
df$Forest <- zscore(df$Forest, na.rm = TRUE)
df$deforest <- as.numeric(df$deforest)
df$mining <- as.numeric(df$mining)
df$MPI <- as.numeric(df$MPI)
df$vectors78 <- as.numeric(df$vectors78)
df$vectors48 <- as.numeric(df$vectors48)
df$consenso3 <- as.factor(df$consenso3)
df$Humidity <- as.numeric(df$Humidity)

# list of regions and municipalities
dlist <- df %>%
  split(.$region) %>%
  map(~ .x %>%
        split(.$municipality) %>%
        map(~ .x %>% 
              select(-region) %>% 
              mutate(period1 = factor(period1)) %>% 
              as.data.frame()))


# model specification
varfun <- "ns"
tlag <- 2
vardf <- 3

# model formula
fmla <- as.formula(Cases ~ cb + ns(Rain,knots=quantile(c(0.17,0.95))) + 
                        ns(SST12,knots=quantile(c(0.27,0.45,0.75))) + 
                        ns(TROP,knots=quantile(c(0.17,0.35,0.75))) +
                        ns(Temperature,knots=quantile(c(0.85,0.92))) +
                        Forest + vectors78 + log(lag1_cases+1)  + 
                        MPI + mining + HFP + deforest + 
                        ns(Year,df=2) + 
                        ns(Month,df=2))


# blank coefficients and vcov matrices
dist <- names(dlist)
b.coef <- matrix(data=NA,nrow=length(dist),ncol=vardf,dimnames=list(dist))
b.vcov <- vector("list",length(dist)); names(b.vcov) <- dist

################################
############# BLUP ############# 
coef1 <- b.coef
vcov1 <- b.vcov
for (j in 1:length(dist)) {
  municipios <- dlist[[dist[j]]]
    for (k in names(municipios)) {
    dat <- municipios[[k]]
    var1 <- dat$SoilMoi
    
    # model and prediction
    argvar <- list(fun = varfun, knots=quantile(var1,c(0.81, 0.85)))  
    arglag <- list(fun = varfun, df = 2)
    cb <- crossbasis(var1, lag=tlag, argvar=argvar, arglag=arglag)
    dat$cb <- cb 
    mod <- glm(fmla, data=dat, family=quasipoisson, maxit=100, 
               offset=log(dat$total_population))
    dlist[[dist[j]]][[k]] <- dat 
    
    # reduction
    red <- crossreduce(cb, mod, cen=median(var1)) 
    coef1[j, ] <- coef(red)
    vcov1[[j]] <- vcov(red)
    rm(dat, argvar, arglag, cb, mod, red, var1)
  }
}

m1 <- mixmeta(formula=coef1~1,S=vcov1,control=list(showiter=TRUE))
summary(m1)
blup1 <- blup(m1,vcov=TRUE) 

###########################################
############# Pooled plot ##############
# model and prediction
t1 <- unlist(lapply(dlist[dist], function(region) {
  lapply(region, function(municipality) municipality$SoilMoi)
}))

tpct <- round(quantile(t1,c(0.01,0.05,0.5,0.95,0.99)),1)
argvar1 <- list(x = t1, fun = varfun, knots = quantile(t1, c(0.81, 0.85), na.rm = TRUE))  
bvar1 <- do.call(onebasis,argvar1)
pred1 <- crosspred(bvar1,coef=m1$coefficients,vcov=m1$vcov,
                   model.link="log",by=0.1,cen=tpct["50%"])
cen1 <- as.numeric(names(which.min(pred1$allRRfit[as.character(seq(tpct["1%"],tpct["99%"],by=0.1))])))
pred1 <- crosspred(bvar1,coef=m1$coefficients,vcov=m1$vcov,
                   model.link="log",by=0.1,cen=cen1)

# Pooled plot
plot(pred1,ylab="Relative Risk",ylim=c(0,5),xlab=" Monthly soil moisture (kg."~m^{-2}*")",xaxt="n",
     main="Cases of malaria",
     lwd=4,col="red",ci.arg=list(density=30,angle=45,col="pink"))
axis(1,at=seq(floor(min(t1)),ceiling(max(t1)),by=0.5))
breaks <- c(min(t1)-1,seq(pred1$predvar[1],
                          pred1$predvar[length(pred1$predvar)],length=30),max(t1)+1)
hist <- hist(t1,breaks=breaks,plot=F)
hist$density <- hist$density/max(hist$density)*0.7
prop <- max(hist$density)/max(hist$counts)
counts <- pretty(hist$count,3)
plot(hist,ylim=c(0,max(hist$density)*3.5),axes=F,ann=F,col="lightgray",
     breaks=breaks,freq=F,add=T)
axis(4,at=counts*prop,labels=counts,cex.axis=0.7)
abline(v=cen1,col="red",lty=1) 


#################################################
############# Municipality level Forest plot ##############

# effects for each municipality
effects <- list()
for (j in 1:length(dist)) {
  municipios <- dlist[[dist[j]]]
  for (k in names(municipios)) {
    dat <- municipios[[k]]
    var1 <- dat$SoilMoi  
    argvar <- list(fun = varfun, knots=quantile(var1,c(0.81, 0.85)))
    arglag <- list(fun = varfun, df = 2)
    cb <- crossbasis(var1, lag=tlag, argvar=argvar, arglag=arglag)
    mod <- glm(fmla, data=dat, family=quasipoisson, maxit=100, 
               offset=log(dat$total_population))
    tpct <- round(quantile(var1, c(0.01, 0.99)), 1)
    argvar1 <- list(x=var1, fun=varfun, df=vardf)
    bvar1 <- do.call(onebasis, argvar1)
    red_temp <- crossreduce(cb, mod, cen=median(var1))
    pred_temp <- crosspred(bvar1, coef=coef(red_temp), vcov=vcov(red_temp),
                           model.link="log", by=0.1, cen=median(var1))
    valid_range <- seq(tpct["1%"], tpct["99%"], by=0.1)
    valid_range_str <- as.character(valid_range)
    valid_range_str <- valid_range_str[valid_range_str %in% names(pred_temp$allRRfit)]
    
    if(length(valid_range_str) > 0) {
      cen1 <- as.numeric(names(which.min(pred_temp$allRRfit[valid_range_str])))
    } else {
      cen1 <- median(var1)
    }
    red <- crossreduce(cb, mod, cen=cen1)  
    effects[[paste(dist[j], k, sep="_")]] <- list(
      coef = coef(red),
      vcov = vcov(red),
      region = dist[j],
      municipality = k,
      mrt = cen1  
    )
  }
}

effects_df <- do.call(rbind, lapply(effects, function(x) {
  data.frame(
    region = x$region,
    municipality = x$municipality,
    log_estimate = x$coef[1],
    log_se = sqrt(diag(x$vcov))[1]
  )
}))

# log scale
effects_df$log_ci_lower <- effects_df$log_estimate - 1.96 * effects_df$log_se
effects_df$log_ci_upper <- effects_df$log_estimate + 1.96 * effects_df$log_se

effects_df$estimate <- exp(effects_df$log_estimate)
effects_df$ci_lower <- exp(effects_df$log_ci_lower)
effects_df$ci_upper <- exp(effects_df$log_ci_upper)

# Pooled effect
pooled_effect <- mixmeta(effects_df$log_estimate ~ 1, S = effects_df$log_se^2, method = "ml")
pooled_log_estimate <- coef(pooled_effect)
pooled_log_se <- sqrt(vcov(pooled_effect))

# CI of pooled effect
pooled_df <- data.frame(
  region = "Pooled",
  municipality = "Pooled",
  log_estimate = pooled_log_estimate[1],
  log_se = pooled_log_se[1],
  log_ci_lower = pooled_log_estimate[1] - 1.96 * pooled_log_se[1],
  log_ci_upper = pooled_log_estimate[1] + 1.96 * pooled_log_se[1]
)

# Transform to scale
pooled_df$estimate <- exp(pooled_df$log_estimate)
pooled_df$ci_lower <- exp(pooled_df$log_ci_lower)
pooled_df$ci_upper <- exp(pooled_df$log_ci_upper)

# Whole effects
all_cols <- c("region", "municipality", "log_estimate", "log_se", "log_ci_lower", 
              "log_ci_upper", "estimate", "ci_lower", "ci_upper")
effects_df_complete <- effects_df[, all_cols]
pooled_df_complete <- pooled_df[, all_cols]
effects_df_complete <- rbind(effects_df_complete, pooled_df_complete)

# forest plot
ggplot(effects_df_complete, aes(x = estimate, y = reorder(municipality, estimate), 
                                xmin = ci_lower, xmax = ci_upper)) +
  geom_point() +
  geom_errorbarh(height = 0.2) +
  facet_grid(region ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(x = "Relative Risk", y = "Municipality", title = "Forest Plot by Municipality and Region") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  theme(strip.text.y = element_text(angle = 0)) 


#################################################
############# Sensitivity analysis ##############

# models
fmla_main <- as.formula(Cases ~ cb + ns(Rain,knots=quantile(c(0.17,0.95))) + 
                          ns(SST12,knots=quantile(c(0.27,0.45,0.75))) + 
                          ns(TROP,knots=quantile(c(0.17,0.35,0.75))) +
                          ns(Temperature,knots=quantile(c(0.85,0.92))) +
                          Forest + vectors78 + log(lag1_cases+1)  + 
                          MPI + mining + HFP + deforest + 
                          ns(Year,df=2) + ns(Month,df=2))

fmla_1 <- as.formula(Cases ~ cb + ns(Rain,knots=quantile(c(0.17,0.95))) + 
                       ns(SST12,knots=quantile(c(0.27,0.45,0.75))) + 
                       ns(TROP,knots=quantile(c(0.17,0.35,0.75))) +
                       ns(Temperature,knots=quantile(c(0.85,0.92))) +
                       Forest + vectors78 + log(lag1_cases+1)  + 
                       MPI + mining + HFP + deforest + 
                       ns(Year,df=3) + ns(Month,df=3))

fmla_2 <- as.formula(Cases ~ cb + ns(Rain,knots=quantile(c(0.17,0.95))) + 
                       ns(SST12,knots=quantile(c(0.27,0.45,0.75))) + 
                       ns(TROP,knots=quantile(c(0.17,0.35,0.75))) +
                       ns(Temperature,knots=quantile(c(0.85,0.92))) +
                       Forest + vectors78 + log(lag1_cases+1)  + 
                       MPI + mining + HFP + deforest + ns(Humidity,df=3) +
                       ns(Year,df=2) + ns(Month,df=2))

fmla_3 <- as.formula(Cases ~ cb + ns(Rain,knots=quantile(c(0.17,0.95))) + 
                       ns(SST12,knots=quantile(c(0.27,0.45,0.75))) + 
                       ns(TROP,knots=quantile(c(0.17,0.35,0.75))) +
                       ns(Temperature,knots=quantile(c(0.85,0.92))) +
                       Forest + vectors48 + log(lag1_cases+1)  + 
                       MPI + mining + HFP + deforest + 
                       ns(Year,df=2) + ns(Month,df=2))

fmla_4 <- as.formula(Cases ~ cb + ns(Rain,knots=quantile(c(0.17,0.95))) + 
                       ns(SST12,knots=quantile(c(0.27,0.45,0.75))) + 
                       ns(TROP,knots=quantile(c(0.17,0.35,0.75))) +
                       ns(Temperature,knots=quantile(c(0.85,0.92))) +
                       Forest + vectors78 + log(lag1_cases+1)  + 
                       MPI + mining + HFP + deforest + consenso3 +
                       ns(Year,df=2) + ns(Month,df=2))

formula_list <- list(
  "Main" = fmla_main,
  "With Year and Month df=3" = fmla_1,
  "With Humidity df=3" = fmla_2,
  "With 4 to 8 vectors" = fmla_3, # same as main model
  "With consensus3" = fmla_4,
  "With lag=3" = fmla_main 
)

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")


#function
process_models <- function(formula_list, dlist) {
  resultados <<- lapply(seq_along(formula_list), function(i) {
    fname <- names(formula_list)[i]
    fmla <- formula_list[[i]]
    
    b.coef <- matrix(NA, nrow=length(dist), ncol=vardf, dimnames=list(dist))
    b.vcov <- vector("list", length(dist))
    names(b.vcov) <- dist
    
    for (j in seq_along(dist)) {
      municipios <- dlist[[dist[j]]]
      coef_temp <- matrix(NA, nrow=length(municipios), ncol=vardf)
      vcov_temp <- vector("list", length(municipios))
      
      for (k in seq_along(municipios)) {
        dat <- municipios[[k]]
        var1 <- dat$SoilMoi
        
        current_lag <- if(i == 6) 3 else tlag
        
        argvar <- list(fun = varfun, knots = quantile(var1, c(0.81, 0.85)))  
        arglag <- list(fun = varfun, df = 2)
        
        cb <- crossbasis(var1, lag = current_lag, argvar = argvar, arglag = arglag)
        dat$cb <- cb
        
        mod <- tryCatch({
          glm(fmla, data = dat, family = quasipoisson, 
              offset = log(dat$total_population), maxit = 100)
        }, error = function(e) NULL)
        
        if (!is.null(mod)) {
          red <- crossreduce(cb, mod, cen = median(var1))
          coef_temp[k, ] <- coef(red)
          vcov_temp[[k]] <- vcov(red)
        }
      }
      
      b.coef[j, ] <- colMeans(coef_temp, na.rm = TRUE)
      b.vcov[[j]] <- Reduce("+", vcov_temp)/length(vcov_temp)
    }
    
    m <- mixmeta(b.coef ~ 1, S = b.vcov)
    blup <- blup(m, vcov = TRUE)
    
    argvar1 <- list(x = t1, fun = varfun, knots = quantile(t1, c(0.81, 0.85), na.rm = TRUE))
    bvar1 <- do.call(onebasis, argvar1)
    pred_temp <- crosspred(bvar1, coef = m$coefficients, vcov = m$vcov,
                           model.link = "log", by = 0.1, cen = median(t1))
    
    
    tpct <- round(quantile(t1, c(0.01, 0.99)), 1)
    valid_range <- seq(tpct[1], tpct[2], by = 0.1)
    valid_range_str <- as.character(valid_range)
    valid_range_str <- valid_range_str[valid_range_str %in% names(pred_temp$allRRfit)]
    
    if(length(valid_range_str) > 0) {
      model_cen <- as.numeric(names(which.min(pred_temp$allRRfit[valid_range_str])))
    } else {
      model_cen <- median(t1)
    }
    
    print(paste("Model:", fname, "- reference value:", model_cen))
    
    
    pred <- crosspred(bvar1, coef = m$coefficients, vcov = m$vcov,
                      model.link = "log", by = 0.1, cen = model_cen,
                      at = seq(min(t1, na.rm=TRUE), max(t1, na.rm=TRUE), by=0.1))
    
    data.frame(
      SoilMoi = pred$predvar,
      RR = pred$allRRfit,
      Low = pred$allRRlow,
      High = pred$allRRhigh,
      Model = fname,
      MRE = model_cen  
    )
  })
  
  combined_data <- do.call(rbind, resultados)
  return(combined_data)
}

combined_data <- process_models(formula_list, dlist)

model_mre_values <- lapply(resultados, function(model_data) {
  data.frame(
    Model = unique(model_data$Model),
    MRE = unique(model_data$MRE)
  )
})
model_mre_df <- do.call(rbind, model_mre_values)

# plot
fig <- ggplot(combined_data, aes(x = SoilMoi, y = RR, color = Model)) +
  geom_line(linewidth = 0.8) +
  ylim(0.5, 2) +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression("Soil moisture (kg."~m^{-2}*")"), 
    y = "RR",
    title = "") +
  geom_vline(xintercept = 42.4, linetype = "dashed", color = "red") 

print(fig)


#negative control
for (j in seq_along(dist)) {
  municipios <- dlist[[dist[j]]]
  
  for (k in names(municipios)) {
    dat <- municipios[[k]]
    var1 <- dat$SoilMoi
    
    # crossbasis
    argvar <- list(fun = varfun, knots = quantile(var1, c(0.81, 0.85), na.rm = TRUE))
    arglag <- list(fun = varfun, df = 2)
    cb <- crossbasis(var1, lag = tlag, argvar = argvar, arglag = arglag)
    dat$cb <- cb
    
    mod <- glm(fmla, data = dat, family = quasipoisson, 
               offset = log(total_population), maxit = 100)
    
    n <- nrow(dat)
    pred_temp <- se_temp <- rep(NA, n)
    
    if(n > tlag) {
      pred_values <- predict(mod, type = "response")
      se_values <- predict(mod, se.fit = TRUE)$se.fit
      pred_temp[(tlag + 1):n] <- pred_values
      se_temp[(tlag + 1):n] <- se_values
    }
    
    dat$predicted_cases <- pred_temp
    dat$predicted_se <- se_temp
    
    dlist[[dist[j]]][[k]] <- dat %>% select(-cb)  
    
    red <- crossreduce(cb, mod, cen = median(var1, na.rm = TRUE))
    b.coef[j, ] <- coef(red)
    b.vcov[[j]] <- vcov(red)
    
    rm(dat, cb, mod, red)
  }
}

#dataframe with predictions
df_pred <- map_df(dlist, ~bind_rows(.x, .id = "municipality"), .id = "region") %>%
  select(region, municipality, Year, Month, period1,
         predicted_cases, predicted_se, SoilMoi, Cases) %>%
  group_by(region, municipality) %>%
  mutate(obs_id = row_number())  

#join data
df_final <- df %>%
  group_by(region, municipality) %>%
  mutate(obs_id = row_number()) %>% 
  left_join(df_pred, by = c("region", "municipality" = "municipality", 
                            "Year", "Month", "period1", 
                            "SoilMoi", "Cases", "obs_id")) %>%
  arrange(region, municipality, obs_id) %>%
  select(-obs_id)


fmla_exten <- update(fmla, . ~ . + lead1_SoilMoi)

cb <- crossbasis(df$SoilMoi, lag = tlag,  
           argvar = list(fun = "ns", knots = quantile(df$SoilMoi, c(0.81, 0.85))),
           arglag = list(fun = "ns", df = 2))

mod_exten <- glm(fmla_exten, data = df_final,
                     family = quasipoisson,
                     offset = log(total_population))

#main model
mod_main <- glm(fmla, data=df, family=quasipoisson, maxit=100, 
           offset=log(df$total_population))

#compare
results <- list(
  Original = tidy(mod_main, conf.int = TRUE),
  Exten = tidy(mod_exten, conf.int = TRUE)
)

#Extract results
confusion_residual <- results$Exten %>%
  filter(term == "lead1_SoilMoi") %>%
  select(term, estimate, p.value, conf.low, conf.high)

print(confusion_residual) 

