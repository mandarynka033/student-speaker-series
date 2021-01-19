library(dplyr)
library(ggplot2)
library(tidyr)

normie <- tibble(values = rnorm(10000))

ggplot(normie, aes(values)) + geom_histogram(bins = 40) + 
  theme_void()

#### Breast Cancer 1 ####
# https://archive.ics.uci.edu/ml/datasets/Breast+Cancer
df <- read.csv("breast-cancer.data", header = FALSE)
colnames(df) <- c("class", "age", "meno", "size", "nodes", "caps", "malig", "breast", "quad", "irradiant")
df$size <- factor(df$size,levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                                     "35-39", "40-44", "45-49", "50-54"))
ggplot(df) + geom_bar(aes(x = age))

ggplot(df) + geom_bar(aes(x = size))

ggplot(df) + geom_bar(aes(x = size)) + facet_wrap(vars(breast))

ggplot(df) + geom_bar(aes(x = size)) + facet_wrap(vars(malig))

barplot(table(df$quad))

#### Breast Cancer 2 ####
# https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Coimbra
cont <- read.csv("dataR2.csv")

hist(cont$Insulin, breaks = 30)

patient_lab <- c("Controls", "Patients")
names(patient_lab) <- c(1, 2)

ggplot(cont) + geom_histogram(aes(x = Insulin), bins = 40) +
  facet_wrap(vars(Classification), 
             labeller = labeller(Classification = patient_lab)) + 
  labs(x = "Insulin (µU/mL)", y = "Count",
       title = "Insulin measurements in healthy controls and breast cancer patients",
       caption = "Breast Cancer Coimbra Dataset, UCI Machine Learning Repository") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 14),
        plot.title = element_text(size = 18))


cont %>% select(Classification, Insulin) %>% 
  group_by(Classification) %>% 
  summarise(var(Insulin), mean(Insulin), n())



x <- cont %>% filter(Classification == 1) %>% select(Insulin)
y <- cont %>% filter(Classification == 2) %>% select(Insulin)

t.test(x$Insulin, y$Insulin)
wilcox.test(x$Insulin, y$Insulin)

set.seed(1729)
cont_r <- sample_n(cont, 80)


ggplot(cont_r) + geom_histogram(aes(x = Insulin), bins = 30) + 
  facet_wrap(vars(Classification), 
             labeller = labeller(Classification = patient_lab)) + 
  labs(x = "Insulin (µU/mL)", y = "Count",
       title = "Insulin measurements in healthy controls and breast cancer patients",
       caption = "Breast Cancer Coimbra Dataset, UCI Machine Learning Repository") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 14),
        plot.title = element_text(size = 18))

#### P Values
of_interest <- cont %>% select(Insulin, Classification) %>% 
  distinct(Insulin, .keep_all = TRUE)
p_values <- tibble(p_param = c(),
                   p_nonparam = c(),
                   nb_patients = c())

for (j in 1:100){
  set.seed(j)
  for (i in 20:nrow(of_interest)){
    sub_cont <- sample_n(of_interest, i)
    x <- filter(sub_cont, Classification == 1)$Insulin
    y <- filter(sub_cont, Classification == 2)$Insulin
    param_test <- t.test(x, y)
    nonparam_test <- wilcox.test(x, y)
    p_values <- p_values %>% 
      add_row(p_param = param_test$p.value,
              p_nonparam = nonparam_test$p.value,
              nb_patients = i)
  }}

p_plt <- p_values %>% group_by(nb_patients) %>% 
  summarise(param_aver = mean(p_param),
            nparam_aver = mean(p_nonparam), 
            .groups = "drop") %>% 
  pivot_longer(cols = c(param_aver, nparam_aver),
               names_to = "test_type", values_to = "value")

ggplot(p_plt, aes(nb_patients, value)) + 
  geom_line(aes(colour = test_type)) + 
  labs(x = "Number of Patients", y = "Average p-value",
       title = "Average p-value from 100 simulations by number of patients") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 14),
        plot.title = element_text(size = 18),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.width = unit(2,"cm") ) +
  scale_color_discrete(labels = c("Mann-Whitney", "t-test")) 

#### Survival ####

# http://www.stat.rice.edu/~sneeley/STAT553/Datasets/survivaldata.txt

trans <- read.csv("transplant.txt", sep = " ", header = FALSE) %>% select_if(~sum(!is.na(.)) > 0) %>% 
  rename(age = V1, censor = V2, time = V3)

hist(trans$time)
library(coin)
library(survival)
library(survminer)

plot(survfit(Surv(time, censor) ~ 1, data = trans))



trans <- trans %>% mutate(age_gap = ifelse(age < 45, 0, 1))
survdiff(Surv(time, censor) ~ age_gap, data = trans)

plt_intro <- ggsurvplot(fit = survfit(Surv(time, censor) ~ 1, 
                                data = trans), 
           conf.int = FALSE, legend = "none", 
           caption = "The Statistical Analysis of Failure Time Data,\n http://www.stat.rice.edu/~sneeley/STAT553/Datasets/survivaldata.txt", xlab = "Days", censor.size = 7,
           title = "Kaplan-Meier curve for the survival fuction after a heart transplant", ggtheme = theme_bw())

ggpar(plt_intro, font.title = 18, font.x = 12, font.y = 12)

ggsurvplot(fit = survfit(Surv(time, censor) ~ 1, data = filter(trans, age_gap == 0)))
ggsurvplot(fit = survfit(Surv(time, censor) ~ 1, data = filter(trans, age_gap == 1)))

survfit(Surv(time, censor) ~ 1, data = trans)


logrank_test(Surv(time, censor) ~ as.factor(age_gap), data = trans)

plot(survfit(Surv(time, censor) ~ age_gap, data = trans))

plt_groups <- ggsurvplot(survfit(Surv(time, censor) ~ age_gap, data = trans), 
           data = trans, conf.int = TRUE, pval = TRUE,
           legend.labs = c("<45", ">45"), 
           ggtheme = theme_bw(),
           xlab = "Days", censor.size = 7,
           title = "Comparison of survival probabilities in different age groups", legend = c(0.85, 0.85),
           legend.title = "Age group",
           conf.int.alpha=0.2)

ggpar(plt_groups, font.title = 18, font.x = 12, font.y = 12,
      font.legend = 14)

#### Bootstrapping ####
btsr <- cont %>% select(Classification, Insulin)
btsr_df <- tibble(final = c())

for (i in 1:1000){
  set.seed(i)
  data <- btsr %>% sample_n(50) %>% group_by(Classification) %>% 
    summarise(means = mean(Insulin), vars = var(Insulin)/n(), .groups = "drop") %>% 
    summarise(denom = sqrt(sum(vars)), numer = diff(means), final = numer / denom) %>% 
    select(final)
  btsr_df <- btsr_df %>% add_row(data)
}

library(scales)
cols <- hue_pal()(2)
colors <- c("Density estimate" = cols[2], "Normal distribution" = cols[1])

ggplot(btsr_df, aes(final)) + geom_histogram(aes(y=..density..), bins = 20) + 
  geom_density(size = 2, aes(color = cols[2])) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(btsr_df$final), sd = sd(btsr_df$final)), 
                aes(color = cols[1]), size = 2) + 
  theme_bw() + 
  labs(x = "t-statistic", y = "Density", 
       title = "Distribution of t-statistics for a Welch's t-test") + 
  theme(plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  scale_color_discrete(labels = c("Estimated density", 
                                  "Normal density"))
