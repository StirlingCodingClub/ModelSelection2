### This is an attempt at plots of predictions of the Takahe data set, comparing the "full" model predictions
###  with the model averaged predictions

### Load the output from "takahe_grueber.R":
#load("takahe_averaged.Rdata")

### Plot the effect of z.F on fledging success, as predicted by both the full model and the model averaged models.

takahe$z.F = (takahe$F-mean(takahe$F))/(2*sd(takahe$F))
z.F_pred = seq(min(takahe$z.F),max(takahe$z.F),0.1)
takahe$psx = paste(takahe$z.F,takahe$Fledge/takahe$Hatch)
psize = as.data.frame(table(takahe$psx))
psize$Freq = log(scale(psize$Freq, center = FALSE, scale = 5))
psize$Freq = psize$Freq-(2*min(psize$Freq))

takahe$psize = psize$Freq[match(takahe$psx, psize$Var1)]
plot(takahe$z.F, takahe$Fledge/takahe$Hatch, cex = takahe$psize,
     xlim = c(-0.35,2), ylim = c(-0.05, 1.05),
     col = scales::alpha("darkgrey", 0.75))

### Predictions from the FULL model:
# full_mod_pred = 
#   predict(stdz.model, type = "response", re.form = NA,
#           newdata = data.frame(z.F = z.F_pred,
#                                Island = factor(1, levels = levels(takahe$Island)),
#                                `I(z.Age^2)` = 0,
#                                z.Age = 0,
#                                z.YearID = 0
#           )
#   )
# lines(z.F_pred, full_mod_pred, lwd = 3, col = "red")

### repeat with manual predictions:
pred_frame = data.frame(Age = 0, 
                        Island = factor("1", levels(takahe$Island)), 
                        YearID = 0, 
                        F = z.F_pred)
pred_mm = model.matrix(~I(Age^2) + Age + Island + YearID + F + F:Age, data = pred_frame)
ppred = pred_mm %*% fixef(stdz.model)
lines(z.F_pred, plogis(ppred), lwd = 3)
pred_vars = diag(pred_mm %*% tcrossprod(vcov(stdz.model), pred_mm))
pred_lo = plogis(ppred-1.96*sqrt(pred_vars))
pred_hi = plogis(ppred+1.96*sqrt(pred_vars))
lines(z.F_pred, pred_lo)
lines(z.F_pred, pred_hi)


### Predictions using the AVERAGED parameter estimates:

### Remake model matrix w/o Age^2 to match coefficient estimates from model averageing.
pred_mm_ma = model.matrix(~Age + F + F:Age + Island + YearID, data = pred_frame)
# Note the model average coefs are in a different order, so we need to reorder these:
ma_coef = ma$coefficients["full",]
ma_coef = ma_coef[c("(Intercept)", "z.Age","z.F","Island2","Island3","Island4","z.YearID","z.Age:z.F")]
ma_coef
colnames(pred_mm_ma)

ppred_ma = pred_mm_ma %*% ma_coef
#lines(z.F_pred, plogis(ppred_ma), col = "green", lwd = 3)
### At this point there is an issue because in vcov(ma) there are NaN values:
vcov(ma, full = FALSE)
### Presumably this is because of the singularity fitting issues. As a result the calculation of prediction errors can't
### be done, so I don't know how to address this using this method.

### We can repeat this instead using AICcmodavg::modavgPred():
predf_ma = data.frame(z.Age = 0, Island = factor("1", levels = levels(takahe$Island)), z.YearID = 0, z.F = z.F_pred)
preds_ma = modavgPred(top.models, newdata = predf_ma, uncond.se = "old")
lines(z.F_pred, preds_ma$mod.avg.pred, col = "red", lwd = 2)
lines(z.F_pred, preds_ma$lower.CL, col = "red", lwd = 1, lty = "dashed")
lines(z.F_pred, preds_ma$upper.CL, col = "red", lwd = 1, lty = "dashed")
### From the comparison of the manual calculation above, the latter method only plots the "FULL" averages (i.e. zero method!)
