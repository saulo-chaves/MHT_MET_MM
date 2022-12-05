rm(list=ls())

# Required packages --------
library(asreml)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(gghighlight)

# Loading the data -------
data = read.table("Data/D2.txt", header = T)
head(data)

# Exploratory analysis --------

## Box-plot ----------

data %>% ggplot(aes(x = env, y = dbh))+
  geom_boxplot(color = 'black', fill = "tomato")+
  labs(y = "Diameter at breast height (cm)", x = "Environment", 
       title = "Diameter at breast height at different environments")

data %>% ggplot(aes(x = gen, y = dbh))+
  geom_boxplot(color = 'black', fill = "tomato")+
  labs(y = "Clones", x = "Environments", 
       title = "Mean diameter at breast height of the evaluated hybrids")+
  theme(axis.text.x = element_text(angle = 90))

## Performance over the sites ------------

data %>% group_by(env, gen) %>% summarise(dbh = mean(dbh, na.rm=T)) %>%
  ggplot(aes(x = env, y = dbh)) +
  stat_summary(aes(group = gen),
               fun = mean, geom = "line", na.rm = TRUE, size = 1)+
  theme(legend.position = "none", legend.title = element_text(),
        axis.text.x = element_text(angle = 90))+
  ylab("Diameter at breast height (cm)")+xlab("Sites")


## Set the factors ------
data = transform(data, env = factor(env),
                 gen = factor(gen),
                 block = factor(block))
ngen = nlevels(data$gen)
nenv = nlevels(data$env)

## Models -------

### Base-line model

m0 = asreml(fixed = dbh ~ env,
            random = ~gen + gen:env + block,
            data = data, 
            maxiter = 100)
sum.m0 = summary(m0)$varcomp
aic.m0 = summary(m0)$aic
bic.m0 = summary(m0)$bic
loglik.m0 = summary(m0)$loglik

mvd.m0 = mean(((predict.asreml(m0,classify = "gen",sed = T)$sed)^2)[
  upper.tri((predict.asreml(m0,classify = "gen",sed = T)$sed)^2,diag=F)
])
PEV.m0 = mean(diag((predict.asreml(m0,classify = "gen",vcov = T)$vcov)))
acc.m0 = sqrt(1-(PEV.m0/sum.m0["gen","component"]))
her.m0 = 1-(mvd.m0/(2*sum.m0["gen","component"]))

blup.m0 = predict.asreml(m0,classify = "gen")$pvals

gain.m0 = (mean(blup.m0[order(blup.m0$predicted.value, 
                              decreasing = T),]$predicted.value[1:30])-
             mean(data$dbh,na.rm=T))/mean(data$dbh,na.rm=T) * 100

### Modelling the residual effects ---------

#### Diagonal

m1 = asreml(fixed = dbh ~ env,
            random = ~gen + gen:env + block,
            residual = ~dsum(~id(units)|env),
            data = data, 
            maxiter = 100)
sum.m1 = summary(m1)$varcomp
aic.m1 = summary(m1)$aic
bic.m1 = summary(m1)$bic
loglik.m1 = summary(m1)$loglik

mvd.m1 = mean(((predict.asreml(m1,classify = "gen",sed = T)$sed)^2)[
  upper.tri((predict.asreml(m1,classify = "gen",sed = T)$sed)^2,diag=F)
])
PEV.m1 = mean(diag((predict.asreml(m1,classify = "gen",vcov = T)$vcov)))
acc.m1 = sqrt(1-(PEV.m1/sum.m1["gen","component"]))
her.m1 = 1-(mvd.m1/(2*sum.m1["gen","component"]))

blup.m1 = predict.asreml(m1,classify = "gen")$pvals

gain.m1 = (mean(blup.m1[order(blup.m1$predicted.value, 
                              decreasing = T),]$predicted.value[1:30])-
             mean(data$dbh, na.rm=T))/mean(data$dbh, na.rm=T) * 100

### Modelling the block effects ---------

#### Diagonal

m2 = asreml(fixed = dbh ~ env,
            random = ~gen + gen:env + at(env):block,
            residual = ~dsum(~id(units)|env),
            data = data, 
            maxiter = 100)
sum.m2 = summary(m2)$varcomp
aic.m2 = summary(m2)$aic
bic.m2 = summary(m2)$bic
loglik.m2 = summary(m2)$loglik

mvd.m2 = mean(((predict.asreml(m2,classify = "gen",sed = T)$sed)^2)[
  upper.tri((predict.asreml(m2,classify = "gen",sed = T)$sed)^2,diag=F)
])
PEV.m2 = mean(diag((predict.asreml(m2,classify = "gen",vcov = T)$vcov)))
acc.m2 = sqrt(1-(PEV.m2/sum.m2["gen","component"]))
her.m2 = 1-(mvd.m2/(2*sum.m2["gen","component"]))

blup.m2 = predict.asreml(m2,classify = "gen")$pvals

gain.m2 = (mean(blup.m2[order(blup.m2$predicted.value, 
                              decreasing = T),]$predicted.value[1:30])-
             mean(data$dbh, na.rm=T))/mean(data$dbh, na.rm=T) * 100

### Modelling the genetic effects ---------

#### Heterogeneous compound symmetry

m3 = asreml(fixed = dbh ~ env,
            random = ~corh(env):gen + at(env):block,
            residual = ~dsum(~id(units)|env),
            data = data, 
            maxiter = 100) 
sum.m3 = summary(m3)$varcomp
aic.m3 = summary(m3)$aic
bic.m3 = summary(m3)$bic
loglik.m3 = summary(m3)$loglik

acc.m3 = NULL
her.m3 = NULL
for(i in 1:nlevels(data$env)){
  pred_vcov = predict(m3, classify = "gen:env",
                      level=list(env = i), vcov = T)
  pred_sed = predict(m3, classify = "gen:env",
                     level=list(env = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m3[i] = sqrt(1-(PEV/sum.m3[grep('gen',rownames(sum.m3)),1][i+1]))
  her.m3[i] = 1-(mvd/(2*sum.m3[grep('gen',rownames(sum.m3)),1][i+1]))
}

acc.m3.mean = mean(acc.m3, na.rm = T)
her.m3.mean = mean(her.m3, na.rm = T)

blup.m3 = predict.asreml(m3,classify = "gen:env")$pvals
blup.m3.mean = blup.m3 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m3 = (mean(blup.m3.mean[order(blup.m3.mean$blup, 
                                   decreasing = T),]$blup[1:30])-
             mean(data$dbh, na.rm=T))/mean(data$dbh, na.rm=T) * 100


#### Unstructured

m4 = asreml(fixed = dbh ~ env,
            random = ~us(env):gen + at(env):block,
            residual = ~dsum(~id(units)|env),
            data = data, 
            maxiter = 100) 
sum.m4 = summary(m4)$varcomp
aic.m4 = summary(m4)$aic
bic.m4 = summary(m4)$bic
loglik.m4 = summary(m4)$loglik

Gvcov = matrix(0, nenv, nenv)
Gvcov[upper.tri(Gvcov, diag = T)] = sum.m4[grep('gen',rownames(sum.m4)),1]
Gvcov[lower.tri(Gvcov, diag = F)] = t(Gvcov)[lower.tri(t(Gvcov), diag = F)]
gencorr = cov2cor(Gvcov)

acc.m4 = NULL
her.m4 = NULL
for(i in 1:nlevels(data$env)){
  pred_vcov = predict(m4, classify = "gen:env",
                      level=list(env = i), vcov = T)
  pred_sed = predict(m4, classify = "gen:env",
                     level=list(env = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m4[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m4[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m4.mean = mean(acc.m4, na.rm = T)
her.m4.mean = mean(her.m4, na.rm = T)

blup.m4 = predict.asreml(m4,classify = "gen:env")$pvals
blup.m4.mean = blup.m4 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m4 = (mean(blup.m4.mean[order(blup.m4.mean$blup, 
                                   decreasing = T),]$blup[1:30])-
             mean(data$dbh, na.rm=T))/mean(data$dbh, na.rm=T) * 100

#### First-order factor analytic

m5 = asreml(fixed = dbh ~ env,
            random = ~fa(env,1):gen + at(env):block,
            residual = ~dsum(~id(units)|env),
            data = data, 
            maxiter = 100,
            predict = predict.asreml(classify = "gen:env")) 
sum.m5 = summary(m5)$varcomp
aic.m5 = summary(m5)$aic
bic.m5 = summary(m5)$bic
loglik.m5 = summary(m5)$loglik

lambda = sum.m5[grep('fa1',rownames(sum.m5)),1]
psi = diag(sum.m5[grep('var',rownames(sum.m5)),1])
Gvcov = lambda %*% t(lambda) + psi
expvar1 = (sum(diag(lambda %*% t(lambda)))/
             sum(diag(lambda %*% t(lambda) + psi)))*100

acc.m5 = NULL
her.m5 = NULL
for(i in 1:nlevels(data$env)){
  pred_vcov = predict(m5, classify = "gen:env",
                      level=list(env = i), vcov = T)
  pred_sed = predict(m5, classify = "gen:env",
                     level=list(env = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m5[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m5[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m5.mean = mean(acc.m5, na.rm = T)
her.m5.mean = mean(her.m5, na.rm = T)

blup.m5 = m5$predictions$pvals
blup.m5.mean = blup.m5 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m5 = (mean(blup.m5.mean[order(blup.m5.mean$blup, 
                                     decreasing = T),]$blup[1:30])-
              mean(data$dbh, na.rm=T))/mean(data$dbh, na.rm=T) * 100


#### Second-order factor analytic

m6 = asreml(fixed = dbh ~ env,
             random = ~fa(env,2):gen + at(env):block,
             residual = ~dsum(~id(units)|env),
             data = data, 
             maxiter = 100,
             predict = predict.asreml(classify = "gen:env")) 
sum.m6 = summary(m6)$varcomp
aic.m6 = summary(m6)$aic
bic.m6 = summary(m6)$bic
loglik.m6 = summary(m6)$loglik

fa1.loadings = sum.m6[grep('fa1', rownames(sum.m6)),1]
fa2.loadings = sum.m6[grep('fa2', rownames(sum.m6)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings = svd(mat.loadings)
mat.loadings.star = mat.loadings %*% svd.mat.loadings$v * - 1
colnames(mat.loadings.star) = c("fa1","fa2")
psi = diag(sum.m6[grep("var", rownames(sum.m6)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar2 = (sum(diag(lamblamb.star))/
             sum(diag(lamblamb.star + psi)))*100
Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi

acc.m6 = NULL
her.m6 = NULL
for(i in 1:nlevels(data$env)){
  pred_vcov = predict(m6, classify = "gen:env",
                      level=list(env = i), vcov = T)
  pred_sed = predict(m6, classify = "gen:env",
                     level=list(env = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m6[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m6[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m6.mean = mean(acc.m6, na.rm = T)
her.m6.mean = mean(her.m6, na.rm = T)

blup.m6 = m6$predictions$pvals
blup.m6.mean = blup.m6 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m6 = (mean(blup.m6.mean[order(blup.m6.mean$blup, 
                                     decreasing = T),]$blup[1:30])-
              mean(data$dbh,na.rm=T))/mean(data$dbh,na.rm=T) * 100

#### Third-order factor analytic

m7 = asreml(fixed = dbh ~ env,
             random = ~fa(env,3):gen + at(env):block,
             residual = ~dsum(~id(units)|env),
             data = data, 
             maxiter = 100,
             predict = predict.asreml(classify = "gen:env")) 
sum.m7 = summary(m7)$varcomp
aic.m7 = summary(m7)$aic
bic.m7 = summary(m7)$bic
loglik.m7 = summary(m7)$loglik

fa1.loadings = sum.m7[grep('fa1', rownames(sum.m7)),1]
fa2.loadings = sum.m7[grep('fa2', rownames(sum.m7)),1]
fa3.loadings = sum.m7[grep('fa3', rownames(sum.m7)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings,fa3.loadings))
svd.mat.loadings = svd(mat.loadings)
mat.loadings.star = mat.loadings %*% svd.mat.loadings$v * - 1
colnames(mat.loadings.star) = c("fa1","fa2","fa3")
psi = diag(sum.m7[grep("var", rownames(sum.m7)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar3 = (sum(diag(lamblamb.star))/
             sum(diag(lamblamb.star + psi)))*100
Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi

acc.m7 = NULL
her.m7 = NULL
for(i in 1:nlevels(data$env)){
  pred_vcov = predict(m7, classify = "gen:env",
                      level=list(env = i), vcov = T)
  pred_sed = predict(m7, classify = "gen:env",
                     level=list(env = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m7[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m7[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m7.mean = mean(acc.m7, na.rm = T)
her.m7.mean = mean(her.m7, na.rm = T)

blup.m7 = m7$predictions$pvals
blup.m7.mean = blup.m7 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m7 = (mean(blup.m7.mean[order(blup.m7.mean$blup, 
                                     decreasing = T),]$blup[1:30])-
              mean(data$dbh,na.rm=T))/mean(data$dbh,na.rm=T) * 100


## Summary: Model selection ------

### Dataframe

modsel = data.frame(
  'Model' = paste0("MET",1:8),
  "R" = c("ID","DIAG","DIAG","DIAG","DIAG","DIAG","DIAG","DIAG"),
  "B" = c("ID","ID","DIAG","DIAG","DIAG","DIAG","DIAG","DIAG"),
  "G" = c("CS", "CS", "CS","CSH","US","FA1","FA2","FA3"),
  'AIC' = c(aic.m0, aic.m1,aic.m2,aic.m3,aic.m4,aic.m5,aic.m6,aic.m7),
  'BIC' = c(bic.m0, bic.m1,bic.m2,bic.m3,bic.m4,bic.m5,bic.m6,bic.m7),
  "Accuracy" = c(acc.m0, acc.m1,acc.m2,acc.m3.mean,acc.m4.mean,acc.m5.mean,
                 acc.m6.mean,acc.m7.mean),
  "ExpeGain" = c(gain.m0,gain.m1,gain.m2,gain.m3,gain.m4,gain.m5,gain.m6,gain.m7),
  "ExpVar" = c('-','-','-','-',"-",expvar1,expvar2,expvar3)
)

### Ranking comparison (Spearman correlation)

rankings = cbind(blup.m0$predicted.value,blup.m1$predicted.value,blup.m2$predicted.value,
                 blup.m3.mean$blup,blup.m4.mean$blup,blup.m5.mean$blup,
                 blup.m6.mean$blup,blup.m7.mean$blup)
rownames(rankings) = blup.m0$gen
colnames(rankings) = c(paste0("MET",1:8))

corrank = cor(rankings, method = 'spearman')

Heatmap(corrank,col=colorRampPalette(brewer.pal(8, "YlOrBr"))(25),
        rect_gp = gpar(col = "white", lwd = 1, type = "none"),
        heatmap_legend_param = list(title="Cor",at=c(0.70,0.85,1),
                                    labels=c("0.70","0.85","1.00")),
        row_order = rownames(corrank), show_heatmap_legend = T,
        column_order = colnames(corrank),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j)  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          if(i >= j) grid.text(sprintf("%.2f", corrank[i, j]), x, y, gp = gpar(fontsize = 20))})

## Genetic correlation between sites ------------

rownames(gencorr) = colnames(gencorr) = levels(data$env)

col_fun = colorRamp2(c(-1,-.5,0,.5,1), c("red4","red1","white",
                                         "darkolivegreen2","darkolivegreen4"))

Heatmap(gencorr,col=col_fun,
        rect_gp = gpar(col = "white", lwd = 1, type = "none"),
        column_dend_height = unit(2, "cm"), 
        row_dend_width = unit(2, "cm"),
        heatmap_legend_param = list(title="Cor",at=c(0.0,0.5,1.0),
                                    labels=c("0.0","0.5","1.0")),
        row_order = rownames(gencorr),
        column_order = colnames(gencorr),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j)  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          if(i >= j) grid.text(sprintf("%.2f", gencorr[i, j]), x, y, gp = gpar(fontsize = 20))})

## Ranking comparison: Base model and chosen model -----------

base = blup.m0[,1:2]
colnames(base) = c('gen','blup.base')
chosen = blup.m4.mean
comparison = cbind(base, blup.chosen = chosen$blup)

plotly::ggplotly(
  comparison %>% arrange(-blup.base) %>% mutate(rank.base = 1:ngen) %>% 
    arrange(-blup.chosen) %>% mutate(rank.chosen = 1:ngen) %>% 
    ggplot(aes(x = rank.base, y = rank.chosen))+
    geom_point(aes(color = gen), size = 2)+
    theme(axis.line = element_line(colour = 'black'), legend.position = 'none',
          panel.background = element_blank(), text = element_text(size = 14))+
    labs(x = "Model 1 - Ranking", y = "Model 12 - Ranking", 
         color = "Genotypes",caption = paste("Ranking correlation =", 
                                             round(cor(comparison$blup.base,
                                                       comparison$blup.chosen,
                                                       method = 'spearman'),3)))+
    scale_x_reverse(breaks = c(1,seq(5,ngen,by=5))) + 
    scale_y_reverse(breaks = c(1,seq(5,ngen,by=5)))
)
