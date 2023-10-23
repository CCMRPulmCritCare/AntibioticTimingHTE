## Project: HTE Analysis

## Author: Jennifer Cano (modified Patricia Kipnis's code)

## Description: Running causal forest on entire VA/KP cohort

library(haven)
library(grf)
library(dplyr)
library(Hmisc)
library(car)
library(aod)
library(knitr)

setwd("filepath") 

sepsis_risk<-read_sas("causalforest_ptid_r.sas7bdat") 

sepsis_risk_df<-as.data.frame(sepsis_risk)

sepsis_risk_df$surv=as.factor(sepsis_risk$survival)

#create age var
summary(sepsis_risk_df$age)
sum(is.na(sepsis_risk_df$age))

sepsis_risk_df$age_over75 = ifelse(sepsis_risk_df$age <=75,0,1)

tapply(sepsis_risk_df$age, sepsis_risk_df$age_over75, summary)

#confirm cancer vars are mutually exclusive
table(sepsis_risk_df$cancer_nonmet, sepsis_risk_df$cancer_met)

#deleted some columns that were later used in the call to the causal forest
sepsis_risk0<-subset(sepsis_risk_df,!is.na(Tx0),select=-c(aod_ind,survival,mort30_ed,Tx,septic_shock)) 

#variables I had deleted: HTN,COAG,SIRS_RR,DM_UNCOMP
df<-sepsis_risk0

#converting all columns to numerical
df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))))

#define outcome as survival
names(df)[names(df) == "surv"] <- "Y"
names(df)[names(df) == "Tx0"] <- "W"

#fit the forest
covariate_names<-c("sirs_temp","sirs_pulse","sirs_rr","sirs_wbc","aod_lung","aod_kidney","aod_liver",
                   "aod_heme","aod_lactate", "pressor_in_72hr", "htn","chf","cardic_arrhym","valvular_d2","pulm_circ",
                   "pvd","paralysis","neuro","pulm","dm_uncomp","dm_comp","hypothyroid","renal",
                   "liver","pud","lymphoma","cancer_met","cancer_nonmet","ra","coag","obesity","wtloss","fen",
                   "anemia_new","male","age_over75")


train_fraction <- 0.8  # Use train_fraction % of the dataset to train our models
n <- dim(df)[1]
train_idx <- sample.int(n, replace=F, size=floor(n*train_fraction))
df_train <- df[train_idx,]
df_test <- df[-train_idx,]


cf <- causal_forest(
  X = as.matrix(df_train[,c("sirs_temp","sirs_pulse","sirs_rr","sirs_wbc","aod_lung","aod_kidney","aod_liver",
                            "aod_heme","aod_lactate","pressor_in_72hr","htn","chf","cardic_arrhym","valvular_d2","pulm_circ",
                            "pvd","paralysis","neuro","pulm","dm_uncomp","dm_comp","hypothyroid","renal",
                            "liver","pud","lymphoma","cancer_met","cancer_nonmet","ra","coag","obesity","wtloss","fen",
                            "anemia_new","male","age_over75")]),
  Y = df_train$Y,
  W = df_train$W,
  num.trees=1000,
  seed = 218604)


## Get predictions from the causal forest
oob_pred <- predict(cf, estimate.variance=TRUE)
head(oob_pred)
oob_tauhat_cf <- oob_pred$predictions
oob_tauhat_cf_se <- sqrt(oob_pred$variance.estimates)

test_pred <- predict(cf, newdata=as.matrix(df_test[covariate_names]), estimate.variance=TRUE)
tauhat_cf_test <- test_pred$predictions

tauhat_cf_test_se <- sqrt(test_pred$variance.estimates)
hist(oob_tauhat_cf, main="Causal forests: out-of-bag CATE")

# Compute Y-star
p <- mean(df_test$W)
Y_star <- ((df_test$W - p)/(p*(1-p)))*df_test$Y

# Compute test mse for all methods
mse<-(Y_star - tauhat_cf_test)^2
mse_summary <- describe(mse)[, c('mean', 'se')]
mses<-rbind(mse,mse_summary)

var_imp <- c(variable_importance(cf))
names(var_imp) <- covariate_names
sorted_var_imp <- sort(var_imp, decreasing=TRUE)

sorted_var_imp

sorted_var_imp_df = data.frame(sorted_var_imp)

##save variable importance values
sink("CF_VarImp_KP_VA_v6.txt")
# print tables with knitr package
sorted_var_imp_df
sink()
closeAllConnections()

#import table after closing session
sorted_var_imp_df = read.delim("CF_VarImp_KP_VA_v6.txt", header=TRUE, sep="", dec=".")

#make row names into a column 
sorted_var_imp_df$Variables = rownames(sorted_var_imp_df)

#move to first column
sorted_var_imp_df = sorted_var_imp_df[, c(2, 1)]

#export table into Word
library(rtf)
rtffile = RTF("filepath/Sorted Var Imp.doc")
addTable(rtffile, sorted_var_imp_df)
done(rtffile)

#create plot
library(ggplot2)

#Y-axis labels
labs = c('Comorbid paralysis','Comorbid peptic ulcer disease',
         'Comorbid obesity','Comorbid lymphoma','Abnormal heart rate','Comorbid diabetes, complicated',
         'Comorbid peripheral vascular disorders','Comorbid renal failure','Comorbid hypothyroidism','Comorbid chronic pulmonary disease',
         'Comorbid rheumatoid arthritis/collagen vascular diseases', 'Abnormal WBC','Comorbid solid tumor, without metastasis',
         'Comorbid hypertension','Abnormal respiratory rate','Comorbid cardiac arrhythmias','Comorbid neurodegenerative disorders',
         'Comorbid fluid and electrolyte disorders','Male sex','Acute lactate elevation',
         'Comorbid weight loss', 'Comorbid valvular disease', 'Comorbid diabetes, uncomplicated','Comorbid congestive heart failure',
         'Comorbid pulmonary circulation disorders','Age >75', 'Comorbid anemia','Acute respiratory dysfunction',
         'Acute renal dysfunction','Abnormal temperature','Comorbid coagulopathy',
         'Acute liver dysfunction','Acute hematologic dysfunction','Comorbid liver disease','Shock',
         'Comorbid metastatic cancer')

svi_plot = ggplot(sorted_var_imp_df) +
  geom_col(aes(x = sorted_var_imp, y = reorder(Variables, sorted_var_imp)), fill = "#076fa2", width = 0.8) +
  xlab("Variable Importance") +
  ylab("") +
  scale_x_continuous(position = "top") +
  scale_y_discrete(labels = labs) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x.top = element_text(margin = margin(0, 0, 10, 0)),
        axis.text.x.top = element_text(vjust = 2))

#export figure into Word
rtffile = RTF("filepath/Sorted Var Imp Plot.doc")
addPlot(rtffile, plot.fun=print, width=8, height=5, res=500, svi_plot)
done(rtffile)

# Manually creating subgroups
num_tiles <- 4  # ntiles = CATE is above / below the median
df_train$cate <- oob_tauhat_cf
df_train$ntile <- factor(ntile(oob_tauhat_cf, n=num_tiles))

head(df_train)

#generate augmented inverse probability weights (AIPW)
ols_sample_ate <- lm("Y ~ ntile + ntile:W", data=df_train)
estimated_sample_ate <- coef(summary(ols_sample_ate))[(num_tiles+1):(2*num_tiles), c("Estimate", "Std. Error")]
hypothesis_sample_ate <- paste0("ntile1:W = ", paste0("ntile", seq(2, num_tiles), ":W"))
ftest_pvalue_sample_ate <- linearHypothesis(ols_sample_ate, hypothesis_sample_ate)[2,"Pr(>F)"]
aipw <- average_treatment_effect(cf)

estimated_aipw_ate <- lapply(
  seq(num_tiles), function(w) {
    ate <- average_treatment_effect(cf, subset = df_train$ntile == w)
  })
estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))

aipw_overall = as.data.frame(matrix(aipw, ncol=2, byrow=TRUE))
names(aipw_overall)[names(aipw_overall) == 'V1'] = 'estimate'
names(aipw_overall)[names(aipw_overall) == 'V2'] = 'std.err'

ATE_table = rbind(aipw_overall, estimated_aipw_ate)

#save table
sink("CF_ATE_KP_VA_v6.txt")
# print tables with knitr package
ATE_table
sink()
closeAllConnections()

#import table after closing session
ATE_table = read.delim("CF_ATE_KP_VA_v6.txt", header=TRUE, sep="", dec=".")

#add columns for High, Low, and Error
ATE_table$Error = 1.96*ATE_table$std.err
ATE_table$High = ATE_table$estimate+ATE_table$Error
ATE_table$Low = ATE_table$estimate-ATE_table$Error

#rename estimate and std.err cols
ATE_table = ATE_table %>% 
  rename(ATE = estimate,
         SE = std.err)

#create row values for overall and each quartile
Quartile = c("Overall", "Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4") 
ATE_table$Quartile = Quartile

#move quartile col to first column
ATE_table = ATE_table[, c(6, 1, 2, 3, 4, 5)]

#export table into Word
library(rtf)
rtffile = RTF("filepath/ATE Table.doc")
addTable(rtffile, ATE_table)
done(rtffile)


#create AIPW Treatment Effect plot
library(ggplot2)

aipw = ggplot(ATE_table, aes(x=Quartile, y=ATE)) + 
  geom_pointrange(aes(ymin = Low, ymax = High)) + 
  ylab("Absolute Mortality Difference") +
  xlab("") +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks=c(0.0, 0.02, 0.04), 
                     labels=label_percent(accuracy = 0.1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 13, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(size = 13, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13))

#export figure into Word
rtffile = RTF("filepath/AIPW Plot2.doc")
addPlot(rtffile, plot.fun=print, width=7, height=5, res=500, aipw)
done(rtffile)


# Testing for equality using Wald test
waldtest_pvalue_aipw_ate <- wald.test(Sigma = diag(estimated_aipw_ate$std.err^2),
                                      b = estimated_aipw_ate$estimate,
                                      Terms = 1:num_tiles)$result$chi2[3]

# Regress each covariate on ntile assignment to means p
cov_means <- lapply(covariate_names, function(covariate) {
  lm(paste0(covariate, ' ~ 0 + ntile'), data = df_train)
})

# Extract the mean and standard deviation of each covariate per ntile
cov_table <- lapply(cov_means, function(cov_mean) {
  as.data.frame(t(coef(summary(cov_mean))[,c("Estimate", "Std. Error")]))
})
names(cov_table) <- covariate_names

vars <- names(cov_table)
est_table <- NULL
for (i in 1:length(vars)){
  est <- cov_table[[i]][1,]
  row.names(est) <- vars[i]
  est_table <- rbind(est_table,est)
}

# print tables with knitr package
sink("CF_means_KP_VA_v6.txt")
est_table
sink()
closeAllConnections()

#import table after closing session
est_table = read.delim("CF_means_KP_VA_v6.txt", header=TRUE, sep="", dec=".")

##Create heatmap table 

#create new var CV

#create function to get population sd
sdpop = function(x) {
  (mean((x-mean(x))^2))^0.5
}

#create CV column
est_table_cv = transform(est_table, CV=(apply(est_table, 1, sdpop))/(apply(est_table, 1, mean)))
est_table_cv

#sort descending order
est_table_cv = est_table_cv[order(est_table_cv$CV, decreasing = TRUE),]


#Create additional variance column
est_table_cv$Variance = ((apply(est_table_cv[,1:4], 1, FUN=max)) - (apply(est_table_cv[,1:4], 1, FUN=min)))/(rowMeans(est_table_cv[,1:4])) 
est_table_cv$Variance = as.numeric(est_table_cv$Variance)
est_table_cv$Variance = round(est_table_cv$Variance, digits=2)

##Format cols as percent
percent = function(x, digits = 1, format = "f", ...) {
  paste0(formatC(x*100, format = format, digits = digits, ...), "%")
}

install.packages("formattable")
library(formattable)

est_table_cv$ntile1 = formattable::percent(est_table_cv$ntile1, digits = 1)
est_table_cv$ntile2 = formattable::percent(est_table_cv$ntile2, digits = 1)
est_table_cv$ntile3 = formattable::percent(est_table_cv$ntile3, digits = 1)
est_table_cv$ntile4 = formattable::percent(est_table_cv$ntile4, digits = 1)
est_table_cv$CV = percent(est_table_cv$CV)

#make row names (aka Characteristic col) into an actual column and enter values
est_table_cv$Characteristic = c("Comorbid metastatic cancer",
                                "Shock",
                                "Acute hematologic dysfunction",
                                "Comorbid liver disease",
                                "Acute liver dysfunction",
                                "Comorbid lymphoma",
                                "Acute respiratory dysfunction",
                                "Comorbid coagulopathy",
                                "Comorbid weight loss",
                                "Age > 75",
                                "Comorbid anemia",
                                "Comorbid neurodegenerative disorders",
                                "Comorbid peptic ulcer disease",
                                "Comorbid valvular disease",
                                "Comorbid obesity",
                                "Comorbid fluid and electrolyte disorders",
                                "Comorbid pulmonary circulation disorders",
                                "Comorbid cardiac arrhythmias",
                                "Comorbid peripheral vascular disorders",
                                "Abnormal temperature",
                                "Comorbid congestive heart failure",
                                "Comorbid diabetes, uncomplicated",
                                "Abnormal respiratory rate",
                                "Comorbid solid tumor, without metastasis ",
                                "Acute lactate elevation",
                                "Comorbid diabetes, complicated",
                                "Comorbid renal failure",
                                "Acute renal dysfunction",
                                "Comorbid paralysis",
                                "Comorbid chronic pulmonary disease",
                                "Comorbid rheumatoid arthritis/collagen",
                                "Comorbid hypothyroidism",
                                "Comorbid hypertension",
                                "Abnormal WBC",
                                "Abnormal heart rate",
                                "Male sex")


#move to first column
est_table_cv = est_table_cv[, c(7, 1, 2, 3, 4, 5, 6)]


#create heatmap
library(gt)
library(scales)
library(readr)


##want color-coding to reflect differences in overall variance
#will calculate differences from the mean for each cell
library(matrixStats)
est_table_cv$mean = rowMeans(as.matrix(est_table_cv[,c(2,3,4,5)]))

est_table_diffs = est_table_cv %>% 
  select(Characteristic, ntile1, ntile2, ntile3, ntile4, CV, Variance, mean) %>% 
  mutate(
    diff1 = (ntile1-mean),
    diff2 = (ntile2-mean),
    diff3 = (ntile3-mean),
    diff4 = (ntile4-mean),
  )

min(est_table_diffs$diff1)
min(est_table_diffs$diff2)
min(est_table_diffs$diff3)
min(est_table_diffs$diff4)
#-0.1450888

max(est_table_diffs$diff1)
max(est_table_diffs$diff2)
max(est_table_diffs$diff3)
max(est_table_diffs$diff4)
#0.1783041


#rename columns
est_table_diffs = est_table_diffs %>% 
  rename("Quartile 1" = "ntile1",
         "Quartile 2" = "ntile2",
         "Quartile 3" = "ntile3",
         "Quartile 4" = "ntile4",
         "CV*" = "CV")


#create heatmap
hm = est_table_diffs %>% 
  gt %>% 
  data_color(columns = 9:12, target_columns = 2:5,
             colors = col_numeric(palette = c("#27AE60", "#F4D03F", "#F01908"),
                                  domain = c(-0.15,0.18))) %>% 
  data_color(columns = 7, 
             colors = col_numeric(palette = c("#27AE60", "#F4D03F", "#F01908"),
                                  domain = c(0.01, 2.9))) %>%  
  
  cols_hide(columns = c(mean, diff1, diff2, diff3, diff4)) %>% 
  
  tab_style(
    locations = cells_body(),
    style = cell_text(color="black"))


#export table into Word
hm |> gtsave("filepath/heatmap_new20230919.docx")

#save df_train data set as csv to then import in SAS to create table comparing quartiles
write.csv(df_train, 'filepath/df_train_ntiles.csv', row.names=FALSE, na= '')



