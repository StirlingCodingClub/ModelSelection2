### This is a simple analysis of the Takahe dataset from the appendix in Grueber et al (2011).
### The dataset (takahe.csv) is a CSV version of Table S1 in their Appendix.

### This analysis serves primarily as a demonstration of model averaging, and is not intended to exactly duplicate
### Grueber et al's analysis.

### Jeroen Minderman, March 2021.

takahe = read.csv("takahe.csv", header = T)
head(takahe)
str(takahe)

takahe$IndID = factor(takahe$IndID)
#takahe$YearID = factor(takahe$YearID)
takahe$Island = factor(takahe$Island)
str(takahe)

nrow(takahe)
nlevels(takahe$IndID)
# nlevels(takahe$YearID)
# nlevels(takahe$Island)

PrFledge = cbind(takahe$Fledge, takahe$Hatch-takahe$Fledge)

library(lme4)
library(arm)
global.model = glmer(PrFledge ~ I(Age^2) + Age + Island + YearID + F + F:Age + (1|IndID), data = takahe, 
                     family = "binomial", na.action = "na.fail")
### Note many "singular fit" warnings when using R version 3.6.3 and lme4_1.1-2.
### I am unsure of the reason for this, the estimates further down seem to be the same as those reported by Grueber et al.
stdz.model = standardize(global.model, standardize.y = FALSE)
stdz.model

### This is just to ensure that arm::standardize() behaves as expected:
### Scale and standardise inputs:
#takahe$sAge = (takahe$Age-mean(takahe$Age))/(2*sd(takahe$Age))
#takahe$sAge2 = takahe$sAge^2
#takahe$sF = (takahe$F-mean(takahe$F))/(2*sd(takahe$F))
#takahe$sYearID = (takahe$YearID-mean(takahe$YearID))/(2*sd(takahe$YearID))
#takahe$sIsland = (takahe$Island-mean(takahe$Island))/(2*sd(takahe$Island))
#man.stdz.model = glmer(PrFledge ~ sAge2 + sAge + sIsland + sYearID + sF + sF:sAge + (1|IndID), data = takahe, family = "binomial")
#summary(man.stdz.model)
#stdz.model
#man.stdz.model

library(MuMIn)
model.set = dredge(stdz.model)
nrow(model.set) # 40 models

### Top model sets
# By delta<2
top.models = get.models(model.set, subset = delta<2)
top.models # Note this is a list of model fits rather than a comparison table
class(top.models)
top.models.tab = model.sel(top.models)
top.models.tab

# By 95% confidence set:
top.model.95 = get.models(model.set, cumsum(weight) <= 0.95)
top.model.95.tab = model.sel(top.model.95)

### Averaging
### Basedon delta <2 set, as per Grueber et al.
top.models.tab
nrow(top.models.tab) # 6 models in this set
# Note in this set, all but one of the original predictors (Age^2) is included in at least one of the models.
top.models.tab
# Also note how some of the predictors occur in only a small number of models in this set; indeed yearID occurs in only one of them.

### "Model averaging" essentially just means taking the estimates from each model in the set, and calculating its
### average weighted by model weights. model.avg() does this for you:
ma = model.avg(top.models)
ma
### We can replciate these esimates by hand, by just calculating the weighted average.
### Here, for the intercept:
# The intercept estimates for each model
top.models.tab$`(Intercept)`
# The weights for each model
top.models.tab$weight
# Multiplied and summed == its weighted average:
sum(top.models.tab$`(Intercept)` * top.models.tab$weight)
# Note how this is the same as the averaged coefficient reported by model.avg():
ma$coefficients

### Difference between "subset" and "full" coefficients. Equivalent to "natural averageing" and "zero" methods. First
### averages estimates only over models each predictor occurs in. Second "assumes" zero coefficient for models where the
### predictor is missing.
### Thus, in current case, for Intercept and z.Age full and subset estimates are the same - because they occur in all models in the set:
top.models.tab
### For all others, the full and subset estimates are different for the others, which do not occur in all models in the
### set. In all cases the full estimates are closer to zero than the subset estimates. This is the result of "shrinkage"
### to zero, because in a number  of cases, "zeros" are averaged. The shrinkage to zero is most pronounced where the
### predictor occurs in fewer models. Because of this our quick and dirty calculation doesn't work:
sum(top.models.tab$z.F * top.models.tab$weight)
### This is because we are trying to multiply NA's:
top.models.tab$z.F

### To replicate these for the Full approach we have to replace the NA's with zeros:
z.F_all = top.models.tab$z.F
z.F_all[is.na(z.F_all)] = 0
### Now the averaging should work:
sum(z.F_all * top.models.tab$weight)
### To replicate the subset ones, we need to ignore the missing ones instead:
z.F_sub = top.models.tab$z.F                       # just take the estimates for z.F
z.F_sub_w = top.models.tab$weight[!is.na(z.F_sub)] # Take the weights for all models, but only for those were z.F is not NA
z.F_sub = z.F_sub[!is.na(z.F_sub)]                 # Remove the NA's from the esimtates for z.F
sum(z.F_sub * z.F_sub_w)                           # Now the weighted average works for the subset, too.
ma

### The same as above applies for all of the predictor averaging.

### To extract more comprehensively summarised estimates, including averaged SE's you can use summary() and the
### resulting $coefmat.full or $coefmat.subset for zero method and NA method averages:
summary(ma)$coefmat.subset
### This should be approximately the same as estimates and SE's in Table A1 in Grueber et al. Just sorted differently:
new_sort = c("(Intercept)","Island2","Island3","Island4","z.Age","z.F","z.YearID","z.Age:z.F")
ma_ests = summary(ma)$coefmat.subset[new_sort,c("Estimate","Adjusted SE")]
### Rounded fo easier reading
ma_ests = round(ma_ests,3)
ma_ests

### We can extract the confidence intervals for these estimates directly, using confint():
round(confint(ma),3)[new_sort,]
### And add these to our ests/SE table:
tab_a1 = as.data.frame(cbind(ma_ests,round(confint(ma),3)[new_sort,]))

### To finish "rebuilding" Table A1 in Grueber et al - variable importance:
imp = importance(ma)
imp
### Match to the table by variable/row name:
tab_a1$Importance = round(as.vector(imp[match(row.names(tab_a1), names(imp))]),2)
### Manually add the importance for Island (as this couldn't be matched due to multiple levels):
tab_a1$Importance[2] = round(imp["Island"],2)
tab_a1

### It would be interesting to compare these estimates with those calculated using the "zero" method instead:
ma_ests_full = summary(ma)$coefmat.full[new_sort,c("Estimate","Adjusted SE")]
ma_ests_full = round(as.data.frame(ma_ests_full),3)
round(ma_ests_full,3)
ma_ests[,1:2]
### Note how the "full" estimates in some cases are closer to zero, owing to shrinkage.

### Finally, compare full model and zero-method average:
round(summary(stdz.model)$coef,3)
round(summary(ma)$coefmat.full,3)
### Note the difference in effect size of z.F. This can be seen in prediction plots - weaker effect for (red) model
### averaged predictions.
source("takahe_grueber_pred.R")
