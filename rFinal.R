#### Library Import ####
library(dplyr)
library(ggplot2)
library(tidyr)
library(coin)
library(survival)
library(survminer)
library(scales)
set.seed(1729)

#### ggplot aesthetics ####
th <- theme_bw() + 
  theme(strip.text.x = element_text(size = 14),
        plot.title = element_text(size = 24),
        legend.position = c(0.85, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.width = unit(2, "cm"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.caption = element_text(size = 14))

#### Visualisation of the shape of a normal distribution ####
normie <- tibble(values = rnorm(10000))

ggplot(normie, aes(values)) + geom_histogram(bins = 40, aes(y=..density..)) + 
  stat_function(fun = dnorm, size = 2, colour = cols[1],
                args = list(mean = mean(normie$values), 
                            sd = sd(normie$values))) +
  theme_void()

#### Breast Cancer data analysis ####
# https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Coimbra
cont <- read.csv("dataR2.csv")

# Labels for plotting
patient_lab <- c("Controls", "Patients")
names(patient_lab) <- c(1, 2)

ggplot(cont) + geom_histogram(aes(x = Resistin), bins = 40) +
  facet_wrap(vars(Classification), 
             labeller = labeller(Classification = patient_lab)) + 
  labs(x = "Resistin (ng/mL)", y = "Count",
       title = "Resistin in healthy controls and cancer patients",
       caption = "Breast Cancer Coimbra Dataset, UCI Machine Learning Repository") + th

# Check for means and variances between the groups
cont %>% select(Classification, Insulin) %>% 
  group_by(Classification) %>% 
  summarise(var(Insulin), mean(Insulin), n())

# Creating vectors for independence testing
x <- cont %>% filter(Classification == 1) %>% select(Insulin)
y <- cont %>% filter(Classification == 2) %>% select(Insulin)

# Parametric and non-parametric tests
t.test(x$Insulin, y$Insulin)
wilcox.test(x$Insulin, y$Insulin)

#### Average P Value by number of patients and test ####
of_interest <- cont %>% select(Resistin, Classification) %>% 
  distinct(Resistin, .keep_all = TRUE)
p_values <- tibble(p_param = c(),
                   p_nonparam = c(),
                   nb_patients = c())

# Conducting tests for different numbers of patients
for (j in 1:1000){
  set.seed(j)
  for (i in 20:nrow(of_interest)){
    sub_cont <- sample_n(of_interest, i)
    x <- filter(sub_cont, Classification == 1)$Resistin
    y <- filter(sub_cont, Classification == 2)$Resistin
    param_test <- t.test(x, y)
    nonparam_test <- wilcox.test(x, y)
    p_values <- p_values %>% 
      add_row(p_param = param_test$p.value,
              p_nonparam = nonparam_test$p.value,
              nb_patients = i)
  }}

# Averaging pvalues
p_plt <- p_values %>% group_by(nb_patients) %>% 
  summarise(param_aver = mean(p_param),
            nparam_aver = mean(p_nonparam), 
            .groups = "drop") %>% 
  pivot_longer(cols = c(param_aver, nparam_aver),
               names_to = "test_type", values_to = "value")

# Plotting
ggplot(p_plt, aes(nb_patients, value)) + 
  geom_line(aes(colour = test_type)) + 
  labs(x = "Number of Patients", y = "Average p-value",
       title = "Average p-value from 1000 simulations by number of patients") + 
  th + scale_color_discrete(labels = c("Mann-Whitney", "t-test")) 

#### Survival analysis ####
# http://www.stat.rice.edu/~sneeley/STAT553/Datasets/survivaldata.txt

trans <- read.csv("transplant.txt", sep = " ", header = FALSE) %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  rename(age = V1, censor = V2, time = V3)

# Plotting a survival curve for all patients
plt_intro <- ggsurvplot(fit = survfit(Surv(time, censor) ~ 1, 
                                data = trans), 
           conf.int = FALSE, legend = "none", 
           caption = "The Statistical Analysis of Failure Time Data,
           http://www.stat.rice.edu/~sneeley/STAT553/Datasets/survivaldata.txt",
           xlab = "Days", censor.size = 7,
           title = "Kaplan-Meier curve for survival data after a transplant",
           ggtheme = theme_bw())
ggpar(plt_intro, font.title = 24, font.x = 18, font.y = 18, font.caption = 16)

# Adding an indicator variable to split the data into subgroups based on age
trans <- trans %>% mutate(age_gap = ifelse(age < 45, 0, 1))

# A Kaplan Meier curve for subgroups 
plt_groups <- ggsurvplot(survfit(Surv(time, censor) ~ age_gap, data = trans), 
           data = trans, conf.int = TRUE, pval = TRUE,
           legend.labs = c("<45", ">45"), 
           ggtheme = theme_bw(),
           xlab = "Days", censor.size = 7,
           title = "Comparison of survival probabilities in different age groups",
           legend = c(0.85, 0.85),
           legend.title = "Age group",
           conf.int.alpha=0.2)
ggpar(plt_groups, font.title = 24, font.x = 18, font.y = 18, font.legend = 18)

#### Calculating pvalues for different subsets of the data ####
btsr <- cont %>% select(Classification, Insulin)
btsr_df <- tibble(final = c())

# Calculating the t statistic for different samples
for (i in 1:1000){
  data <- btsr %>% sample_n(50) %>% group_by(Classification) %>% 
    summarise(means = mean(Insulin), vars = var(Insulin)/n(), 
              .groups = "drop") %>% 
    summarise(denom = sqrt(sum(vars)), numer = diff(means), 
              final = numer / denom) %>% 
    select(final)
  btsr_df <- btsr_df %>% add_row(data)
}

# Colours for the plot
cols <- hue_pal()(2)
colors <- c("Density estimate" = cols[2], "Normal distribution" = cols[1])

# Plotting normal and estimated density over a histogram of the p values
ggplot(btsr_df, aes(final)) + geom_histogram(aes(y=..density..), bins = 20) + 
  geom_density(size = 2, aes(color = cols[2])) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(btsr_df$final), sd = sd(btsr_df$final)), 
                aes(color = cols[1]), size = 2) + 
  labs(x = "t-statistic", y = "Density", 
       title = "Distribution of t-statistics for a Welch's t-test") + th +
  scale_color_discrete(labels = c("Estimated density", 
                                  "Normal density")) 
