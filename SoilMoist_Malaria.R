library(dlnm) ; library(mvmeta) ; library(splines) ; library(dplyr); library(ggplot2); library(RCurl); library(MKdescr)

# LOAD THE DATASET
url_path = "https://raw.githubusercontent.com/juandavidgutier/soil_moisture_malaria/master/top50.csv"
muni_Col <- read.csv(url_path)
dim(muni_Col)
head(muni_Col)

#sd units
muni_Col$Rain <- zscore(muni_Col$Rain, na.rm = TRUE)
muni_Col$Runoff <- zscore(muni_Col$Runoff, na.rm = TRUE)
muni_Col$SST12 <- zscore(muni_Col$SST12, na.rm = TRUE)
muni_Col$SST3 <- zscore(muni_Col$SST3, na.rm = TRUE)
muni_Col$SST34 <- zscore(muni_Col$SST34, na.rm = TRUE)
muni_Col$SST4 <- zscore(muni_Col$SST4, na.rm = TRUE)
muni_Col$NATL <- zscore(muni_Col$NATL, na.rm = TRUE)
muni_Col$SATL <- zscore(muni_Col$SATL, na.rm = TRUE)
muni_Col$TROP <- zscore(muni_Col$TROP, na.rm = TRUE)

# REGIONS
regions <- as.character(unique(muni_Col$municipality)) 

# CREATE A LIST WITH THE REGIONAL SERIES
data <- lapply(regions,function(x) muni_Col[muni_Col$municipality==x,])
names(data) <- regions
m <- length(regions)


# SoilMoi RANGES
ranges <- t(sapply(data, function(x) range(x$SoilMoi,na.rm=T)))

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}


################################################################################
#Lag 0-1 months
lag <- c(0,1)

#FIRST STAGE
bound <- colMeans(ranges)

argvar <- list(fun="poly", degree=2, cen=0.0) 
arglag <- list(fun="ns", df=1, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=2, cen=0.0) 
arglag2 <- list(fun="poly", degree=1,intercept=FALSE)
argvar3 <- list(fun="bs", df=2, cen=0.0) 
arglag3 <- list(fun="poly", degree=1, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall3 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep="")))

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "reml" # PLEASE, CHANGE IT TO "ml" and "mm" TO DEVELOP SENSITIVITY TEST

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF SoilMoi AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:100/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF SoilMoi AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.0,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])



# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Fig 2A
SoilMoi <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

SoilMoi_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(SoilMoi, best_model, best_model_h, best_model_l))
names <- c("SoilMoi", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

f2A = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= best_model_l,  
                  ymax= best_model_h), 
              fill='grey',alpha=0.3) + 
  geom_line(aes(x=SoilMoi,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("a") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f2A)


# Q TEST AND I-SQUARE

(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

#90TH AND 95TH PERCENTILES
print(SoilMoi_tiles)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["37.8",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["38.2",]),3)


#META-REGRESSION 
# INPUT THE META-VARIABLE: occu4.8
occu_df <- muni_Col %>%
  group_by(Code.DANE) %>%
  summarise(occu4.8 = mean(occu4.8))

occu4.8 <- occu_df$occu4.8

# MULTIVARIATE META-REGRESSION
(mvallhfp <- update(mvall,.~occu4.8))
summary(mvallhfp)


# PREDICTION FROM META-REGRESSION
val <- round(quantile(occu4.8,c(10,90)/100),1)
predall <- predict(mvallhfp,data.frame(occu4.8=val),vcov=T)

cpallhfpat10 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallhfpat90 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallhfp <- qtest(mvallhfp))
round(((qallhfp$Q-qallhfp$df)/qallhfp$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallhfp,"occu4.8"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","occu4.8")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall, Sall, method = "reml") 
mvhfp <- mvmeta(yall ~ occu4.8, Sall, method = "reml", ) 

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvhfp, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

#Fig 3A 
SoilMoi <- as.data.frame(cpallhfpat90$predvar)
at90 <- as.data.frame(cpallhfpat90$allRRfit)
at90_h <- as.data.frame(cpallhfpat90$allRRhigh)
at90_l <- as.data.frame(cpallhfpat90$allRRlow)
at10 <- as.data.frame(cpallhfpat10$allRRfit)
at10_h <- as.data.frame(cpallhfpat10$allRRhigh)
at10_l <- as.data.frame(cpallhfpat10$allRRlow)
data_runoff <- as.data.frame(cbind(SoilMoi, at90, at90_h, at90_l, at10, at10_h, at10_l))
names <- c("SoilMoi", "at90", "at90_h", "at90_l", "at10", "at10_h", "at10_l")
colnames(data_runoff) <-  names

f3A = ggplot(data_runoff) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at90 - at90_l, #ymin= best_model - conf_int[1],
                  ymax= at90 + at90_h), 
              fill='red',alpha=0.3) + #percentile 90th of occu4.8
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at10 - at10_l,
                  ymax= at10 + at10_h),
              fill='blue',alpha=0.3) + #percentile 10th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at90), linewidth=1.05,
            color='red') + #percentile 90th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at10), linewidth=1.05,
            color='blue') + #percentile 10th of occu4.8
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("a") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f3A)



####################################################################################
#Lag 0-2 months
lag <- c(0,2)

#FIRST STAGE
bound <- colMeans(ranges)

argvar <- list(fun="poly", degree=2, cen=0.0) 
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.0) 
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.0) 
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
yall3 <- yall2 

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "reml" ## PLEASE, CHANGE IT TO "ml" and "mm" TO DEVELOP SENSITIVITY TEST

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF SoilMoi AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:200/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF SoilMoi AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.0,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])



# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Fig 2B
SoilMoi <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

SoilMoi_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(SoilMoi, best_model, best_model_h, best_model_l))
names <- c("SoilMoi", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

f2B = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.3) + 
  geom_line(aes(x=SoilMoi,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("b") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f2B)



# Q TEST AND I-SQUARE

(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

#90TH AND 95TH PERCENTILES
print(SoilMoi_tiles)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["37.8",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["38.2",]),3)


#META-REGRESSION 
# INPUT THE META-VARIABLE: occu4.8
occu_df <- muni_Col %>%
  group_by(Code.DANE) %>%
  summarise(occu4.8 = mean(occu4.8))

occu4.8 <- occu_df$occu4.8

# MULTIVARIATE META-REGRESSION
(mvallhfp <- update(mvall,.~occu4.8))
summary(mvallhfp)


# PREDICTION FROM META-REGRESSION
val <- round(quantile(occu4.8,c(10,90)/100),1)
predall <- predict(mvallhfp,data.frame(occu4.8=val),vcov=T)

cpallhfpat10 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallhfpat90 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallhfp <- qtest(mvallhfp))
round(((qallhfp$Q-qallhfp$df)/qallhfp$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallhfp,"occu4.8"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","occu4.8")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "reml") 
mvhfp <- mvmeta(yall2 ~ occu4.8, Sall2, method = "reml", ) 

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvhfp, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

#Fig 3B
SoilMoi <- as.data.frame(cpallhfpat90$predvar)
at90 <- as.data.frame(cpallhfpat90$allRRfit)
at90_h <- as.data.frame(cpallhfpat90$allRRhigh)
at90_l <- as.data.frame(cpallhfpat90$allRRlow)
at10 <- as.data.frame(cpallhfpat10$allRRfit)
at10_h <- as.data.frame(cpallhfpat10$allRRhigh)
at10_l <- as.data.frame(cpallhfpat10$allRRlow)
data_runoff <- as.data.frame(cbind(SoilMoi, at90, at90_h, at90_l, at10, at10_h, at10_l))
names <- c("SoilMoi", "at90", "at90_h", "at90_l", "at10", "at10_h", "at10_l")
colnames(data_runoff) <-  names

f3B = ggplot(data_runoff) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at90 - at90_l, #ymin= best_model - conf_int[1],
                  ymax= at90 + at90_h), 
              fill='red',alpha=0.3) + #percentile 90th of occu4.8
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at10 - at10_l,
                  ymax= at10 + at10_h),
              fill='blue',alpha=0.3) + #percentile 10th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at90), linewidth=1.05,
            color='red') + #percentile 90th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at10), linewidth=1.05,
            color='blue') + #percentile 10th of occu4.8
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("b") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f3B)




####################################################################################
#Lag 0-3 months
lag <- c(0,3)

#FIRST STAGE
bound <- colMeans(ranges)

argvar <- list(fun="poly", degree=2, cen=0.0) 
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.0) 
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.0) 
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 

yall3 <- yall2 

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "reml" # PLEASE, CHANGE IT TO "ml" and "mm" TO DEVELOP SENSITIVITY TEST

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF SoilMoi AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:300/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF SoilMoi AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.0,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])



# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Fig 2C
SoilMoi <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

SoilMoi_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(SoilMoi, best_model, best_model_h, best_model_l))
names <- c("SoilMoi", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

f2C = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.3) + 
  geom_line(aes(x=SoilMoi,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("c") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f2C)



# Q TEST AND I-SQUARE

(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

#90TH AND 95TH PERCENTILES
print(SoilMoi_tiles)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["37.8",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["38.2",]),3)


#META-REGRESSION 
# INPUT THE META-VARIABLE: occu4.8
occu_df <- muni_Col %>%
  group_by(Code.DANE) %>%
  summarise(occu4.8 = mean(occu4.8))

occu4.8 <- occu_df$occu4.8

# MULTIVARIATE META-REGRESSION
(mvallhfp <- update(mvall,.~occu4.8))
summary(mvallhfp)


# PREDICTION FROM META-REGRESSION
val <- round(quantile(occu4.8,c(10,90)/100),1)
predall <- predict(mvallhfp,data.frame(occu4.8=val),vcov=T)

cpallhfpat10 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallhfpat90 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallhfp <- qtest(mvallhfp))
round(((qallhfp$Q-qallhfp$df)/qallhfp$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallhfp,"occu4.8"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","occu4.8")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "reml") 
mvhfp <- mvmeta(yall2 ~ occu4.8, Sall2, method = "reml", ) 

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvhfp, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

#Fig 3C 
SoilMoi <- as.data.frame(cpallhfpat90$predvar)
at90 <- as.data.frame(cpallhfpat90$allRRfit)
at90_h <- as.data.frame(cpallhfpat90$allRRhigh)
at90_l <- as.data.frame(cpallhfpat90$allRRlow)
at10 <- as.data.frame(cpallhfpat10$allRRfit)
at10_h <- as.data.frame(cpallhfpat10$allRRhigh)
at10_l <- as.data.frame(cpallhfpat10$allRRlow)
data_runoff <- as.data.frame(cbind(SoilMoi, at90, at90_h, at90_l, at10, at10_h, at10_l))
names <- c("SoilMoi", "at90", "at90_h", "at90_l", "at10", "at10_h", "at10_l")
colnames(data_runoff) <-  names

f3C = ggplot(data_runoff) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at90 - at90_l, #ymin= best_model - conf_int[1],
                  ymax= at90 + at90_h), 
              fill='red',alpha=0.3) + #percentile 90th of occu4.8
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at10 - at10_l,
                  ymax= at10 + at10_h),
              fill='blue',alpha=0.3) + #percentile 10th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at90), linewidth=1.05,
            color='red') + #percentile 90th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at10), linewidth=1.05,
            color='blue') + #percentile 10th of occu4.8
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("c") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f3C)





################################################################################
#Lag 0-4 months
lag <- c(0,4)

#FIRST STAGE
# MAIN MODEL
argvar <- list(fun="poly", degree=2, cen=0.0) 
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.0) 
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.0) 
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 

yall3 <- yall2 


# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$SoilMoi,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + Rain + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})


# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "reml" # PLEASE, CHANGE IT TO "ml" and "mm" TO DEVELOP SENSITIVITY TEST


# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)



# BASES OF SoilMoi AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 1:400/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF SoilMoi AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.0,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])

#print(exp(cpall$coefficients))
#print(exp(confint(cpall, method="Wald")))

#RR for SoilMoi for each city
#exp(mvall[["model"]])

# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Fig 2D
SoilMoi <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

SoilMoi_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(SoilMoi, best_model, best_model_h, best_model_l))
names <- c("SoilMoi", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

f2D = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.3) + 
  geom_line(aes(x=SoilMoi,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=SoilMoi_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("d") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f2D)


# Q TEST AND I-SQUARE


(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["37.8",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["38.2",]),3)

#META-REGRESSION 
# INPUT THE META-VARIABLE: occu4.8
occu_df <- muni_Col %>%
  group_by(Code.DANE) %>%
  summarise(occu4.8 = mean(occu4.8))

occu4.8 <- occu_df$occu4.8

# MULTIVARIATE META-REGRESSION
(mvallhfp <- update(mvall,.~occu4.8))
summary(mvallhfp)


# PREDICTION FROM META-REGRESSION
val <- round(quantile(occu4.8,c(10,90)/100),1)
predall <- predict(mvallhfp,data.frame(occu4.8=val),vcov=T)

cpallhfpat10 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallhfpat90 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallhfp <- qtest(mvallhfp))
round(((qallhfp$Q-qallhfp$df)/qallhfp$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallhfp,"occu4.8"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","occu4.8")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "reml") 
mvhfp <- mvmeta(yall2 ~ occu4.8, Sall2, method = "reml", ) 

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvhfp, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)


#Fig 3D 
SoilMoi <- as.data.frame(cpallhfpat90$predvar)
at90 <- as.data.frame(cpallhfpat90$allRRfit)
at90_h <- as.data.frame(cpallhfpat90$allRRhigh)
at90_l <- as.data.frame(cpallhfpat90$allRRlow)
at10 <- as.data.frame(cpallhfpat10$allRRfit)
at10_h <- as.data.frame(cpallhfpat10$allRRhigh)
at10_l <- as.data.frame(cpallhfpat10$allRRlow)
data_runoff <- as.data.frame(cbind(SoilMoi, at90, at90_h, at90_l, at10, at10_h, at10_l))
names <- c("SoilMoi", "at90", "at90_h", "at90_l", "at10", "at10_h", "at10_l")
colnames(data_runoff) <-  names

f3D = ggplot(data_runoff) +
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at90 - at90_l, #ymin= best_model - conf_int[1],
                  ymax= at90 + at90_h), 
              fill='red',alpha=0.3) + #percentile 90th of occu4.8
  geom_ribbon(aes(x=SoilMoi,
                  ymin= at10 - at10_l,
                  ymax= at10 + at10_h),
              fill='blue',alpha=0.3) + #percentile 10th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at90), linewidth=1.05,
            color='red') + #percentile 90th of occu4.8
  geom_line(aes(x=SoilMoi,
                y=at10), linewidth=1.05,
            color='blue') + #percentile 10th of occu4.8
  geom_hline(aes(yintercept=1), color="black", linetype="dashed") +
  theme_bw() +
  labs(x = expression("Soil moisture (kg."~m^{-2}*")"), y = "RR", size = 14) +
  ggtitle("d") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))
print(f3D)


