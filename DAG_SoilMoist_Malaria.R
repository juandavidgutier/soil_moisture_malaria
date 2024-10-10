library(ggdag)
library(dagitty)
library(lavaan)
library(dplyr)
library(GGally)
library(tidyr)
library(MKdescr)



#implied Conditional Independences
url_path = "https://raw.githubusercontent.com/juandavidgutier/soil_moisture_malaria/master/top50.csv"
dataset <- read.csv(url_path)
dataset <- select(dataset, Cases, SoilMoi, Rain, Runoff, SST12, SST3, SST34, SST4, NATL, SATL, TROP)
dataset <- dataset[complete.cases(dataset), ] 

#sd units
dataset$SoilMoi <- zscore(dataset$SoilMoi, na.rm = TRUE)
dataset$Rain <- zscore(dataset$Rain, na.rm = TRUE)
dataset$Runoff <- zscore(dataset$Runoff, na.rm = TRUE)
dataset$SST12 <- zscore(dataset$SST12, na.rm = TRUE)
dataset$SST3 <- zscore(dataset$SST3, na.rm = TRUE)
dataset$SST34 <- zscore(dataset$SST34, na.rm = TRUE)
dataset$SST4 <- zscore(dataset$SST4, na.rm = TRUE)
dataset$NATL <- zscore(dataset$NATL, na.rm = TRUE)
dataset$SATL <- zscore(dataset$SATL, na.rm = TRUE)
dataset$TROP <- zscore(dataset$TROP, na.rm = TRUE)

str(dataset)

#descriptive analysis
#ggpairs(dataset)

#DAG 
dag <- dagitty('dag {
Cases [pos="0, 0.5"]
SoilMoi  [pos="-1, 0.5"]

Rain [pos="-2.4, -1.2"]
Runoff [pos="-1.0, -1.0"]

SST3 [pos="-1.8, 1.3"]
SST4 [pos="-1.9, 1.4"]
SST34 [pos="-2, 1.5"]
SST12 [pos="-2.1, 1.6"]
NATL [pos="-2.2, 1.7"]
SATL [pos="-2.3, 1.8"]
TROP [pos="-2.4, 1.9"]


SST12 -> SST3
SST12 -> SST34
SST12 -> SST4
SST12 -> NATL
SST12 -> SATL
SST12 -> TROP

SST3 -> SST34
SST3 -> SST4
SST3 -> NATL
SST3 -> SATL
SST3 -> TROP

SST34 -> SST4
SST34 -> NATL
SST34 -> SATL
SST34 -> TROP

SST4 -> NATL
SST4 -> SATL
SST4 -> TROP

NATL -> SATL
NATL -> TROP

SATL -> TROP

SST12 -> SoilMoi
SST3 -> SoilMoi
SST34 -> SoilMoi
SST4 -> SoilMoi
NATL -> SoilMoi
SATL -> SoilMoi
TROP -> SoilMoi

SST12 -> Cases
SST3 -> Cases
SST34 -> Cases
SST4 -> Cases
NATL -> Cases
SATL -> Cases
TROP -> Cases

SST12 -> Runoff
SST3 -> Runoff
SST34 -> Runoff
SST4 -> Runoff
NATL -> Runoff
SATL -> Runoff
TROP -> Runoff


SST12 -> Rain
SST3 -> Rain
SST34 -> Rain
SST4 -> Rain
NATL -> Rain
SATL -> Rain
TROP -> Rain

Rain  -> SoilMoi
Rain  -> Runoff
Rain  -> Cases

Runoff -> Cases
SoilMoi -> Runoff
SST12 -> SoilMoi
SST3  -> SoilMoi
SST34 -> SoilMoi
SST4 -> SoilMoi
NATL -> SoilMoi
SATL -> SoilMoi
TROP -> SoilMoi


SoilMoi -> Cases

}')  


plot(dag)


## check whether any correlations are perfect (i.e., collinearity)
myCov <- cov(dataset)
round(myCov, 2)

myCor <- cov2cor(myCov)
noDiag <- myCor
diag(noDiag) <- 0
any(noDiag == 1)

## if not, check for multicollinearity (i.e., is one variable a linear combination of 2+ variables?)
det(myCov) < 0
## or
any(eigen(myCov)$values < 0)


## Independencias condicionales
impliedConditionalIndependencies(dag, max.results=3)
corr <- lavCor(dataset, group=dataset$Code.DANE)

summary(corr)

# Note there are not conditional independences
localTests(dag, sample.cov=corr, sample.nobs=nrow(dataset), max.conditioning.variables=3)



#identification
simple_dag <- dagify(
  Cases ~  SoilMoi + SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP  + Rain + Runoff,
  SoilMoi ~ SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP + Rain ,
  SST12 ~ SST3 + SST34 + SST4 + NATL + SATL +  TROP,
  SST3 ~ SST34 + SST4 + NATL + SATL +  TROP,
  SST34 ~ SST4 + NATL + SATL +  TROP,
  SST4 ~ NATL + SATL +  TROP,
  NATL ~ SATL +  TROP,
  SATL ~  TROP,
  Rain ~ SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP,
  Runoff ~ SST12 + SST3 + SST34 + SST4 + NATL + SATL + TROP  + Rain + SoilMoi,
  exposure = "SoilMoi",
  outcome = "Cases",
  coords = list(x = c(SoilMoi=2, Rain=1, Cases=2, SST12=3, SST3=3.1, SST34=3.2, SST4=3.3, NATL=3.4, SATL=3.5, TROP=3.6,
                      Runoff=0.5),
                y = c(SoilMoi=2, Rain=3, Cases=1, SST12=3, SST3=3.1, SST34=3.2, SST4=3.3, NATL=3.4, SATL=3.5, TROP=3.6,
                      Runoff=1.9))
)


# theme_dag
ggdag(simple_dag) + 
  theme_dag()

ggdag_status(simple_dag) +
  theme_dag()


#adjust
adjustmentSets(simple_dag,  type = "minimal")

ggdag_adjustment_set(simple_dag, shadow = TRUE) +
  theme_dag()

