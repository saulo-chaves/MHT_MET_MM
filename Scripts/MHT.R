rm(list=ls())

# Required packages --------
library(asreml)
library(tidyverse)
library(desplot)
library(ComplexHeatmap)
library(circlize)
library(gifski)

# Loading the data -------
data = read.csv("Data/D1.csv", header = T, sep = ";")

# Exploratory analysis --------

## Box-plot ----------

data %>% ggplot(aes(x = yr, y = fy))+
  geom_boxplot(color = 'black', fill = "tomato")+
  labs(y = "Fruit yield", x = "Harvest years", 
       title = "Fruit yield over the harvest years")

data %>% ggplot(aes(x = gen, y = fy))+
  geom_boxplot(color = 'black', fill = "tomato")+
  labs(y = "Fruit yield", x = "Hybrids", 
       title = "Mean fruit yield of the evaluated hybrids")+
  theme(axis.text.x = element_text(angle = 90))

## Performance over the harvest years ------------

data %>% group_by(yr, gen) %>% summarise(fy = mean(fy, na.rm=T)) %>%
  ggplot(aes(x = yr, y = fy)) +
  stat_summary(aes(group = gen),
               fun = mean, geom = "line", na.rm = TRUE, size = 1)+
  theme(legend.position = "none", legend.title = element_text(),
        axis.text.x = element_text(angle = 90))+
  ylab("Fruit yield")+xlab("Harvest years")

## Spatial distribution --------

plotlist = list()

for(i in unique(data$yr)){
  plotlist[[i]] = ggdesplot(fy ~ col * row, text = gen, cex = 1,
                            data = data[data$yr == i,], flip = T, out1 = block,
                            shorten = "no", show.key = T,ticks = T, 
                            col.regions=colorRampPalette(brewer.pal(8, "YlOrBr"))(25),
                            ylab = "Row", xlab = "Column", main = paste("Harvest year:",i)) +
    labs(fill = "Fruit yield (kg)") +
    theme(text = element_text(size = 15))
  
}

gif_file = "Scripts/Outputs/MET/Spat_distr.gif"
save_gif(
  for(i in seq_along(plotlist)){
    plot(plotlist[[i]])
  }, gif_file,width = 1000, height = 720, delay = 3 
)


# Statistical analyses -----
## Set the factors ------
data = transform(data, yr = factor(yr),
                       gen = factor(gen),
                       block = factor(block),
                       plot = factor(plot),
                       row = factor(row),
                       col = factor(col))
ngen = nlevels(data$gen)
nharv = nlevels(data$yr)

data = data %>% arrange(yr, col, row)

## Models -------

### Base-line model

m0 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen + gen:yr + plot,
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
                             decreasing = T),]$predicted.value[1:7])-
             mean(data$fy))/mean(data$fy) * 100

### Modelling the residual effects ---------

#### Adding the spatial adjustment

m1 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen + gen:yr + plot,
            residual = ~id(yr):ar1(col):ar1(row),
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
                              decreasing = T),]$predicted.value[1:7])-
             mean(data$fy))/mean(data$fy) * 100

#### Diagonal

m2 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen + gen:yr + plot,
            residual = ~dsum(~id(yr):ar1(col):ar1(row)|yr),
            data = data, 
            maxiter = 100)
m2 = update(m2)
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
                              decreasing = T),]$predicted.value[1:7])-
             mean(data$fy))/mean(data$fy) * 100

#### First-order autoregressive

m3 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen + gen:yr + plot,
            residual = ~ar1(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100)

sum.m3 = summary(m3)$varcomp
aic.m3 = summary(m3)$aic
bic.m3 = summary(m3)$bic
loglik.m3 = summary(m3)$loglik

mvd.m3 = mean(((predict.asreml(m3,classify = "gen",sed = T)$sed)^2)[
  upper.tri((predict.asreml(m3,classify = "gen",sed = T)$sed)^2,diag=F)
])
PEV.m3 = mean(diag((predict.asreml(m3,classify = "gen",vcov = T)$vcov)))
acc.m3 = sqrt(1-(PEV.m3/sum.m3["gen","component"]))
her.m3 = 1-(mvd.m3/(2*sum.m3["gen","component"]))

blup.m3 = predict.asreml(m3,classify = "gen")$pvals

gain.m3 = (mean(blup.m3[order(blup.m3$predicted.value, 
                              decreasing = T),]$predicted.value[1:7])-
             mean(data$fy))/mean(data$fy) * 100

#### First-order heterogeneous autoregressive

m4 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen + gen:yr + plot,
            residual = ~ar1h(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100)

sum.m4 = summary(m4)$varcomp
aic.m4 = summary(m4)$aic
bic.m4 = summary(m4)$bic
loglik.m4 = summary(m4)$loglik

mvd.m4 = mean(((predict.asreml(m4,classify = "gen",sed = T)$sed)^2)[
  upper.tri((predict.asreml(m4,classify = "gen",sed = T)$sed)^2,diag=F)
])
PEV.m4 = mean(diag((predict.asreml(m4,classify = "gen",vcov = T)$vcov)))
acc.m4 = sqrt(1-(PEV.m4/sum.m4["gen","component"]))
her.m4 = 1-(mvd.m4/(2*sum.m4["gen","component"]))

blup.m4 = predict.asreml(m4,classify = "gen")$pvals

gain.m4 = (mean(blup.m4[order(blup.m4$predicted.value, 
                              decreasing = T),]$predicted.value[1:7])-
             mean(data$fy))/mean(data$fy) * 100

### Modelling the genetic effects ---------

#### Heterogeneous compound symmetry

m5 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen:corh(yr) + plot,
            residual = ~ar1h(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100) 
m5 = update(m5)
sum.m5 = summary(m5)$varcomp
aic.m5 = summary(m5)$aic
bic.m5 = summary(m5)$bic
loglik.m5 = summary(m5)$loglik

acc.m5 = NULL
her.m5 = NULL
for(i in 1:nlevels(data$yr)){
  pred_vcov = predict(m5, classify = "gen:yr",
                        level=list(yr = i), vcov = T)
  pred_sed = predict(m5, classify = "gen:yr",
                     level=list(yr = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m5[i] = sqrt(1-(PEV/sum.m5[grep('gen',rownames(sum.m5)),1][i+1]))
  her.m5[i] = 1-(mvd/(2*sum.m5[grep('gen',rownames(sum.m5)),1][i+1]))
}

acc.m5.mean = mean(acc.m5, na.rm = T)
her.m5.mean = mean(her.m5, na.rm = T)

blup.m5 = predict.asreml(m5,classify = "gen:yr")$pvals
blup.m5.mean = blup.m5 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m5 = (mean(blup.m5.mean[order(blup.m5.mean$blup, 
                              decreasing = T),]$blup[1:7])-
             mean(data$fy))/mean(data$fy) * 100


#### First-order autoregressive

m6 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen:ar1(yr) + plot,
            residual = ~ar1h(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100)
sum.m6 = summary(m6)$varcomp
aic.m6 = summary(m6)$aic
bic.m6 = summary(m6)$bic
loglik.m6 = summary(m6)$loglik

mvd.m6 = mean(((predict.asreml(m6,classify = "gen",sed = T)$sed)^2)[
  upper.tri((predict.asreml(m6,classify = "gen",sed = T)$sed)^2,diag=F)
])
PEV.m6 = mean(diag((predict.asreml(m6,classify = "gen",vcov = T)$vcov)))
acc.m6 = sqrt(1-(PEV.m6/sum.m6["gen:yr","component"]))
her.m6 = 1-(mvd.m6/(2*sum.m6["gen:yr","component"]))

blup.m6 = predict.asreml(m6,classify = "gen:yr")$pvals
blup.m6.mean = blup.m6 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m6 = (mean(blup.m6.mean[order(blup.m6.mean$blup, 
                                   decreasing = T),]$blup[1:7])-
             mean(data$fy))/mean(data$fy) * 100

#### First-order heterogeneous autoregressive

m7 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen:ar1h(yr) + plot,
            residual = ~ar1h(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100)
m7 = update(m7)
sum.m7 = summary(m7)$varcomp
aic.m7 = summary(m7)$aic
bic.m7 = summary(m7)$bic
loglik.m7 = summary(m7)$loglik

acc.m7 = NULL
her.m7 = NULL
for(i in 1:nlevels(data$yr)){
  pred_vcov = predict(m7, classify = "gen:yr",
                      level=list(yr = i), vcov = T)
  pred_sed = predict(m7, classify = "gen:yr",
                     level=list(yr = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m7[i] = sqrt(1-(PEV/sum.m7[grep('gen',rownames(sum.m7)),1][i+1]))
  her.m7[i] = 1-(mvd/(2*sum.m7[grep('gen',rownames(sum.m7)),1][i+1]))
}

acc.m7.mean = mean(acc.m7, na.rm = T)
her.m7.mean = mean(her.m7, na.rm = T)

blup.m7 = predict.asreml(m7,classify = "gen:yr")$pvals
blup.m7.mean = blup.m7 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m7 = (mean(blup.m7.mean[order(blup.m7.mean$blup, 
                                   decreasing = T),]$blup[1:7])-
             mean(data$fy))/mean(data$fy) * 100


#### Unstructured

m8 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen:us(yr) + plot,
            residual = ~ar1h(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100) # did not converge

#### First-order factor analytic

m9 = asreml(fixed = fy ~ yr + block:yr,
            random = ~gen:fa(yr,1) + plot,
            residual = ~ar1h(yr):ar1(col):ar1(row),
            data = data, 
            maxiter = 100, 
            predict = predict.asreml(classify = "gen:yr"))
m9 = update(m9)
sum.m9 = summary(m9)$varcomp
aic.m9 = summary(m9)$aic
bic.m9 = summary(m9)$bic
loglik.m9 = summary(m9)$loglik

lambda = sum.m9[grep('fa1',rownames(sum.m9)),1]
psi = diag(sum.m9[grep('var',rownames(sum.m9)),1])
Gvcov = lambda %*% t(lambda) + psi
expvar1 = (sum(diag(lambda %*% t(lambda)))/
             sum(diag(lambda %*% t(lambda) + psi)))*100

acc.m9 = NULL
her.m9 = NULL
for(i in 1:nlevels(data$yr)){
  pred_vcov = predict(m9, classify = "gen:yr",
                      level=list(yr = i), vcov = T)
  pred_sed = predict(m9, classify = "gen:yr",
                     level=list(yr = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m9[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m9[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m9.mean = mean(acc.m9, na.rm = T)
her.m9.mean = mean(her.m9, na.rm = T)

blup.m9 = m9$predictions$pvals
blup.m9.mean = blup.m9 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m9 = (mean(blup.m9.mean[order(blup.m9.mean$blup, 
                                   decreasing = T),]$blup[1:7])-
             mean(data$fy))/mean(data$fy) * 100


#### Second-order factor analytic

m10 = asreml(fixed = fy ~ yr + block:yr,
             random = ~gen:fa(yr,2) + plot,
             residual = ~ar1h(yr):ar1(col):ar1(row),
             data = data, 
             maxiter = 100, 
             predict = predict.asreml(classify = "gen:yr"))

sum.m10 = summary(m10)$varcomp
aic.m10 = summary(m10)$aic
bic.m10 = summary(m10)$bic
loglik.m10 = summary(m10)$loglik

fa1.loadings = sum.m10[grep('fa1', rownames(sum.m10)),1]
fa2.loadings = sum.m10[grep('fa2', rownames(sum.m10)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings = svd(mat.loadings)
mat.loadings.star = mat.loadings %*% svd.mat.loadings$v * - 1
colnames(mat.loadings.star) = c("fa1","fa2")
psi = diag(sum.m10[grep("var", rownames(sum.m10)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar2 = (sum(diag(lamblamb.star))/
             sum(diag(lamblamb.star + psi)))*100
Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi
gencorr = cov2cor(Gvcov)

acc.m10 = NULL
her.m10 = NULL
for(i in 1:nlevels(data$yr)){
  pred_vcov = predict(m10, classify = "gen:yr",
                      level=list(yr = i), vcov = T)
  pred_sed = predict(m10, classify = "gen:yr",
                     level=list(yr = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m10[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m10[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m10.mean = mean(acc.m10, na.rm = T)
her.m10.mean = mean(her.m10, na.rm = T)

blup.m10 = m10$predictions$pvals
blup.m10.mean = blup.m10 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m10 = (mean(blup.m10.mean[order(blup.m10.mean$blup, 
                                     decreasing = T),]$blup[1:7])-
              mean(data$fy))/mean(data$fy) * 100

#### Third-order factor analytic

m11 = asreml(fixed = fy ~ yr + block:yr,
             random = ~gen:fa(yr,3) + plot,
             residual = ~ar1h(yr):ar1(col):ar1(row),
             data = data, 
             maxiter = 100, 
             predict = predict.asreml(classify = "gen:yr"))
m11 = update(m11)
sum.m11 = summary(m11)$varcomp
aic.m11 = summary(m11)$aic
bic.m11 = summary(m11)$bic
loglik.m11 = summary(m11)$loglik

fa1.loadings = sum.m11[grep('fa1', rownames(sum.m11)),1]
fa2.loadings = sum.m11[grep('fa2', rownames(sum.m11)),1]
fa3.loadings = sum.m11[grep('fa3', rownames(sum.m11)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings,fa3.loadings))
svd.mat.loadings = svd(mat.loadings)
mat.loadings.star = mat.loadings %*% svd.mat.loadings$v * - 1
colnames(mat.loadings.star) = c("fa1","fa2","fa3")
psi = diag(sum.m11[grep("var", rownames(sum.m11)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar3 = (sum(diag(lamblamb.star))/
             sum(diag(lamblamb.star + psi)))*100
Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi
gencorr = cov2cor(Gvcov)

acc.m11 = NULL
her.m11 = NULL
for(i in 1:nlevels(data$yr)){
  pred_vcov = predict(m11, classify = "gen:yr",
                      level=list(yr = i), vcov = T)
  pred_sed = predict(m11, classify = "gen:yr",
                     level=list(yr = i), sed = T)
  pred_sed$pvals 
  PEV = mean(diag(pred_vcov$vcov))
  mvd = mean((pred_sed$sed^2)[upper.tri(pred_sed$sed^2,diag = F)])
  acc.m11[i] = sqrt(1-(PEV/Gvcov[i,i]))
  her.m11[i] = 1-(mvd/(2*Gvcov[i,i]))
}

acc.m11.mean = mean(acc.m11, na.rm = T)
her.m11.mean = mean(her.m11, na.rm = T)

blup.m11 = m11$predictions$pvals
blup.m11.mean = blup.m11 %>% group_by(gen) %>%
  summarise(blup = mean(predicted.value))

gain.m11 = (mean(blup.m11.mean[order(blup.m11.mean$blup, 
                                     decreasing = T),]$blup[1:7])-
              mean(data$fy))/mean(data$fy) * 100


## Summary: Model selection ------

### Dataframe

modsel = data.frame(
  'Model' = paste0("MHT",1:12),
  "R" = c("ID", "ID + SPT","DIAG + SPT","AR1 + SPT","AR1H + SPT","AR1H + SPT",
          "AR1H + SPT","AR1H + SPT","AR1H + SPT","AR1H + SPT","AR1H + SPT","AR1H + SPT"),
  "G" = c("CS", "CS", "CS","CS","CS","CSH","AR1","AR1H","US","FA1","FA2","FA3"),
  'AIC' = c(aic.m0, aic.m1,aic.m2,aic.m3,aic.m4,aic.m5,aic.m6,
            aic.m7,"-",aic.m9, aic.m10,aic.m11),
  'BIC' = c(bic.m0, bic.m1,bic.m2,bic.m3,bic.m4,bic.m5,bic.m6,
            bic.m7,"-",bic.m9, bic.m10,bic.m11),
  "Accuracy" = c(acc.m0, acc.m1,acc.m2,acc.m3,acc.m4,acc.m5.mean,acc.m6,
                 acc.m7.mean,"-",acc.m9.mean, acc.m10.mean,acc.m11.mean),
  "ExpeGain" = c(gain.m0,gain.m1,gain.m2,gain.m3,gain.m4,gain.m5,gain.m6,gain.m7,
                 "-",gain.m9,gain.m10,gain.m11),
  "ExpVar" = c('-','-','-','-','-','-','-','-','-',expvar1,expvar2,expvar3)
)

write.csv(modsel, file = "Scripts/Outputs/MHT/modsel_MHT.csv",row.names = F)

### Ranking comparison (Spearman correlation)

rankings = cbind(blup.m0$predicted.value,blup.m1$predicted.value,blup.m2$predicted.value,
                 blup.m3$predicted.value,blup.m4$predicted.value,blup.m5.mean$blup,
                 blup.m6.mean$blup,blup.m7.mean$blup,blup.m9.mean$blup,
                 blup.m10.mean$blup,blup.m11.mean$blup)
rownames(rankings) = blup.m0$gen
colnames(rankings) = c(paste0("MHT",1:7),paste0("MHT",9:12))

corrank = cor(rankings, method = 'spearman')

Heatmap(corrank,col=colorRampPalette(brewer.pal(8, "YlOrBr"))(25),
        rect_gp = gpar(col = "white", lwd = 1, type = "none"),
        column_dend_height = unit(2, "cm"), 
        row_dend_width = unit(2, "cm"),show_heatmap_legend = F,
        row_order = rownames(corrank),
        column_order = colnames(corrank),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j)  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          if(i >= j) grid.text(sprintf("%.2f", corrank[i, j]), x, y, gp = gpar(fontsize = 20))})

## Genetic correlation between harvest-years ------------

rownames(gencorr) = colnames(gencorr) = levels(data$yr)

col_fun = colorRamp2(c(-1,-.5,0,.5,1), c("red4","red1","white",
                                         "darkolivegreen2","darkolivegreen4"))

Heatmap(gencorr,col=col_fun,
        rect_gp = gpar(col = "white", lwd = 1, type = "none"),
        column_dend_height = unit(2, "cm"), 
        row_dend_width = unit(2, "cm"),show_heatmap_legend = T,
        row_order = rownames(gencorr),
        column_order = colnames(gencorr),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j)  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          if(i >= j) grid.text(sprintf("%.2f", gencorr[i, j]), x, y, gp = gpar(fontsize = 20))})

## Ranking comparison: Base model and chosen model -----------

base = blup.m0[,1:2]
colnames(base) = c('gen','blup.base')
chosen = blup.m11.mean
comparison = cbind(base, blup.chosen = chosen$blup)

comparison %>% arrange(-blup.base) %>% mutate(rank.base = 1:ngen) %>% 
  arrange(-blup.chosen) %>% mutate(rank.chosen = 1:ngen) %>% 
  ggplot(aes(x = rank.base, y = rank.chosen))+
  geom_point(aes(color = gen), size = 2)+
  theme(axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), text = element_text(size = 14))+
  labs(x = "Model 1 - Ranking", y = "Model 12 - Ranking", 
       color = "Genotypes",caption = paste("Ranking correlation =", 
                                           round(cor(comparison$blup.base,
                                                     comparison$blup.chosen,
                                                     method = 'spearman'),4)))+
  scale_x_reverse(breaks = c(1,seq(5,ngen,by=5))) + 
  scale_y_reverse(breaks = c(1,seq(5,ngen,by=5))) +
  gghighlight(rank.base %in% 1:7 |
                rank.chosen %in% 1:7,
              label_key = gen,keep_scales = F, use_direct_label = F)


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
                                                                      method = 'spearman'),4)))+
                   scale_x_reverse(breaks = c(1,seq(5,ngen,by=5))) + 
                   scale_y_reverse(breaks = c(1,seq(5,ngen,by=5)))
)
                  
