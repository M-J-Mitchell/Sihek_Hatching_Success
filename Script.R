##%##########################################################################%##
#                                                                              #
##%#                         Full Modelling Script                          #%##
#                                                                              #
#############################################################################%##
#### 1)   Necessary Packages & settings ####

library(eeptools)
library(lme4)
library(tidyverse)
library(corrr)
library(brms)
library(ggtext)
library(bayestestR)
library(ggdist)
library(beepr)
library(plotrix)
library(plyr)
library(loo)

options(future.globals.maxSize= 9891289600)

#### 2)   Functions & Search terms ####

ropeci = 1 # proportion of posterior which intersects with ROPE. 

# For more readable parameter names in plots
fn_labels <- function(string){
  string <- sub('b_','',string)
  string <- gsub('Dam.age.years','Maternal age',string)
  string <- gsub('dam_age_sq_sc',"Maternal age\u00B2",string)
  string <- gsub('Sire.age.years','Paternal age',string)
  string <- gsub('sire_age_sq_sc','Paternal age\u00B2',string)
  string <- gsub('Dam.Inbreeding','Maternal *f*',string)
  string <- gsub('Sire.Inbreeding','Paternal *f*',string)  
  string <- gsub('Pair.Inbreeding','Egg *f*', string)
  string <- gsub('Intercept','Intercept',string)
  string <- gsub('Artificial.incubationY','AI',string)
  string <- gsub("sd_Dam.Studbook.ID__Intercept", "SD (Dam ID intercept)", string)
  string <- gsub("sd_Sire.Studbook.ID__Intercept", "SD (Sire ID intercept)", string)
  string <- gsub("sd_Location__Intercept", "SD (Location intercept)", string)
  string <- gsub("sd_Management__Intercept", "SD (Management stage intercept)", string)
  string <- gsub(": $","",string) # remove colon at end of string (to remove "Mating:" if no subsequent variable)
  
  return(string)
  
}
#### 3)   Import data ####

data_raw <- read.csv("./Data/01_Data_ANON.csv", header = T)

## Add pair ID
data_raw$pair_ID<- paste(data_raw$Dam.Studbook.ID,data_raw$Sire.Studbook.ID, sep = "_")
data <- data_raw

#### 2.5) Egg summary numbers

# Number of egg records
nrow(data_raw)

# Years data is from
table(data_raw$Lay.Year)

# Number of institutions:
n_distinct(data_raw$Location)

# Number of eggs with unknown development
table(data_raw$Fertile)

#### 4)   Subset data to remove incomplete data and outliers ####

# Eggs laid by lone females and females in aviaries split from their mate
sum(is.na(data_raw$Fertile)) # 128
data <- subset(data, !is.na(Fertile))

# total developing/hatched eggs from whole dataset of eggs laid by pairs
table(data$Fertile)
table(data$Hatch.tot)

# Unknown parental hatching dates
length(which(data$Dam.DOH == "U" | data$Sire.DOH == "U")) # 76
data <- subset(data, !Dam.DOH == "U")
data <- subset(data, !Sire.DOH %in% c("U"))

# Unknown Lay date
length(which(data$Lay.Date=="U")) # 49
data <- subset(data, !Lay.Date %in% c("U"))

# Fertility is unknown
length(which(data$Fertile == "U")) #186
data <- subset(data, !Fertile %in% c("U"))
# drop unused factors (removes the U)
data<- data %>% mutate(Fertile = factor(Fertile))

# total eggs removed
982-543 # 439

# Remove pair with remarkably high kinship
length(which(data$Pair.Inbreeding >0.2))
data <- data %>% filter(Pair.Inbreeding < 0.25)

# Remove institution with only 2 records
table(data$Location)
data <- data %>% filter(Location != "Loc12")


#### 5)   Organise dates & ages ####

# Specify dates as dates 
data$Dam_DOH_use<- format(as.Date(data$Dam.DOH, format="%d/%m/%Y"))
data$Sire_DOH_use<- format(as.Date(data$Sire.DOH, format="%d/%m/%Y"))
data$Lay.date_use<- format(as.Date(data$Lay.Date, format="%d/%m/%Y"))

# Calculate parent ages in years and age difference in a pair
data$Dam.age.years<-round(age_calc(dob=as.Date(data$Dam_DOH_use), enddate=as.Date(data$Lay.date_use), units = "years", precise = TRUE), 0)
data$Sire.age.years<-round(age_calc(dob=as.Date(data$Sire_DOH_use), enddate=as.Date(data$Lay.date_use), units = "years", precise = TRUE), 0)

#### 6)   Set data types ####

data$Hatch.tot <- as.factor(data$Hatch.tot)
data$Dam.Inbreeding <- as.numeric((data$Dam.Inbreeding))
data$Sire.Inbreeding <- as.numeric((data$Sire.Inbreeding))
data$Pair.Inbreeding <- as.numeric((data$Pair.Inbreeding))
data$Artificial.incubation <- as.factor(data$Artificial.incubation)

data <- data %>% mutate(Sire.Studbook.ID = as.factor(Sire.Studbook.ID),
                        Dam.Studbook.ID = as.factor(Dam.Studbook.ID),
                        Location = as.factor(Location),
                        Lay.Year = as.factor(Lay.Year))

# Reduce dataset to only necessary columns

data <- data %>% select(Fertile,
                        Hatch.tot,
                        Sire.age.years, 
                        Dam.age.years,
                        Sire.Inbreeding,
                        Dam.Inbreeding,
                        Pair.Inbreeding,
                        Location,
                        Sire.Studbook.ID,
                        Dam.Studbook.ID,
                        Lay.Year,
                        Management,
                        Artificial.incubation,
                        pair_ID)

# round inbreeding coefficients
# unique(data$Sire.Inbreeding)
# data$Sire.Inbreeding <- round(data$Sire.Inbreeding, digits = 3)

# unique(data$Dam.Inbreeding)
# data$Dam.Inbreeding <- round(data$Dam.Inbreeding, digits = 3)


#### 7)   Summary stats ####

# Range of inbreeding coefficients in dataset
summary(data$Dam.Inbreeding)
summary(data$Sire.Inbreeding)
summary(data$Pair.Inbreeding)

# Number of pairs in usable dataset
length(unique(data$pair_ID))
pairs <- as.data.frame(unclass(table(data$pair_ID)))
summary(pairs$`unclass(table(data$pair_ID))`)
std.error(pairs$`unclass(table(data$pair_ID))`)

# Age range in dataset
summary(data$Dam.age.years)
summary(data$Sire.age.years)

#### 8)   Create datasets ####

## Egg viability dataset: 
## Duplicate data into a new set for analysis
data_dev <- data

# Scale and centre inputs to mean of zero & SD = 0.5
data_dev_sc <- data_dev %>% mutate(
  Dam.age.years = 0.5*scale(Dam.age.years, center=T, scale=T) ,
  Sire.age.years = 0.5*scale(Sire.age.years, center=T, scale=T),
  Dam.Inbreeding = 0.5*scale(Dam.Inbreeding, center=T, scale=T),
  Sire.Inbreeding = 0.5*scale(Sire.Inbreeding, center=T, scale=T),
  Pair.Inbreeding = 0.5*scale(Pair.Inbreeding, center=T, scale=T),
)

# Convert matrix columns to vectors (needs to be done as a result of using scale)
data_dev_sc <- data_dev_sc %>% mutate(Dam.age.years = as.numeric(Dam.age.years),
                                      Sire.age.years = as.numeric(Sire.age.years),
                                      Dam.Inbreeding = as.numeric(Dam.Inbreeding),
                                      Sire.Inbreeding = as.numeric(Sire.Inbreeding),
                                      Pair.Inbreeding = as.numeric(Pair.Inbreeding),
                                      Sire.Studbook.ID =as.factor(Sire.Studbook.ID),
                                      Dam.Studbook.ID =as.factor(Dam.Studbook.ID),
                                      Lay.Year=as.factor(Lay.Year),
                                      Location =as.factor(Location),
                                      Management = as.factor(Management)
)

## create variables to determine non-linear effects
data_dev_sc$sire_age_sq_sc <- data_dev_sc$Sire.age.years^2
data_dev_sc$dam_age_sq_sc <- data_dev_sc$Dam.age.years^2

# Prepare data for model 
data_dev_sc <-as_tibble(data_dev_sc)

data_dev_sc$Fertile <- as.character(data_dev_sc$Fertile)
data_dev_sc$Fertile[data_dev_sc$Fertile == "Y"] <- "1"
data_dev_sc$Fertile[data_dev_sc$Fertile == "N"] <- "0"
data_dev_sc$Fertile <- as.factor(data_dev_sc$Fertile)

data_dev_sc$Hatch.tot <- as.character(data_dev_sc$Hatch.tot)
data_dev_sc$Hatch.tot[data_dev_sc$Hatch.tot == "Y"] <- "1"
data_dev_sc$Hatch.tot[data_dev_sc$Hatch.tot == "N"] <- "0"
data_dev_sc$Hatch.tot <- as.factor(data_dev_sc$Hatch.tot)


#### 9)   Collinearity testing ####

data_corr_check <- data_dev %>% select(Dam.Inbreeding,Sire.Inbreeding,Pair.Inbreeding)

correlate(data_corr_check)

## model with all terms
model_dev_full <- brm(Fertile ~ 
                        Dam.age.years + dam_age_sq_sc +
                        Sire.age.years + sire_age_sq_sc +
                        Sire.Inbreeding + Pair.Inbreeding +
                        Dam.Inbreeding + 
                        (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                        (1|Management) + (1|Location),
                      data = data_dev_sc,
                      family = "bernoulli",
                      prior = NULL,
                      iter = 15000,
                      control=list(adapt_delta=0.99),
                      cores = 8,
                      sample_prior = TRUE,
                      seed = 13579)

## model only kinship
model_dev_pair <- brm(Fertile ~ 
                        Dam.age.years + dam_age_sq_sc +
                        Sire.age.years + sire_age_sq_sc +
                        Sire.Inbreeding + Pair.Inbreeding +
                        (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                        (1|Management) + (1|Location),
                      data = data_dev_sc,
                      family = "bernoulli",
                      prior = NULL,
                      iter = 15000,
                      control=list(adapt_delta=0.99),
                      cores = 8,
                      sample_prior = TRUE,
                      seed = 13579)

## model only dam
model_dev_dam <- brm(Fertile ~ 
                       Dam.age.years + dam_age_sq_sc +
                       Sire.age.years + sire_age_sq_sc +
                       Sire.Inbreeding +
                       Dam.Inbreeding +
                       (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                       (1|Management) + (1|Location),
                     data = data_dev_sc,
                     family = "bernoulli",
                     prior = NULL,
                     iter = 15000,
                     control=list(adapt_delta=0.99),
                     cores = 8,
                     sample_prior = TRUE,
                     seed = 13579)


## R2 values
bayes_R2(model_dev_full)
bayes_R2(model_dev_pair)
bayes_R2(model_dev_dam)

################################################################################
#### 10)  Egg viability models ####
#### 10a) Pair k ####
model_dev_pair <- brm(Fertile ~ 1 +
                       Dam.age.years + dam_age_sq_sc +
                       Sire.age.years + sire_age_sq_sc +
                       Sire.Inbreeding + Pair.Inbreeding +
                       (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                       (1|Management) + (1|Location),
                     data = data_dev_sc,
                     family = "bernoulli",
                     prior = NULL,
                     iter = 15000,
                     control=list(adapt_delta=0.99),
                     cores = 8,
                     sample_prior = TRUE,
                     seed = 13579,
                     save_pars = save_pars(all = TRUE))

bayes_R2(model_dev_pair)
plot(rope(model_dev_pair))

# describe posterior distributions # 1 devergent
describe_posterior(model_dev_pair, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
dev_pair_post <- describe_posterior(model_dev_pair, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_pair_post <- as.data.frame(dev_pair_post)
write.table(base:::format(dev_pair_post, digits=2), file = "./Supp_Info/Tables/01_dev_pair_post.txt", sep = ",", quote = FALSE, row.names = F)

#### 10b) Pair k: Egg viability model cross validation ####

# run reduced model based on undecided terms from full model
model_dev_pair_red <- brm(Fertile ~ 1 +
                            Sire.age.years + sire_age_sq_sc +
                            Sire.Inbreeding + Pair.Inbreeding +
                            (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID) +
                            (1|Management) + (1|Location),
                          data = data_dev_sc,
                          family = "bernoulli",
                          prior = NULL,
                          iter = 15000,
                          control=list(adapt_delta=0.99),
                          cores = 8,
                          sample_prior = TRUE,
                          seed = 13579,
                          save_pars = save_pars(all = TRUE))

describe_posterior(model_dev_pair_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_dev_pair_red)
plot(rope(model_dev_pair_red))

# save table
dev_pair_post_red <- describe_posterior(model_dev_pair_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_pair_post_red <- as.data.frame(dev_pair_post_red)
write.table(base:::format(dev_pair_post_red, digits=2), file = "./Supp_Info/Tables/09_dev_pair_post_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
dev_pair_loo <- loo(model_dev_pair, moment_match = TRUE, recompile = TRUE, 
                    save_psis = TRUE, reloo = TRUE)
dev_pair_red_loo <- loo(model_dev_pair_red, moment_match = TRUE, recompile = TRUE, 
                        save_psis = TRUE, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(dev_pair_loo, dev_pair_red_loo)

#### 10c) Dam f ####

model_dev_dam <- brm(Fertile ~ 1 +
                            Dam.age.years + dam_age_sq_sc +
                            Sire.age.years + sire_age_sq_sc +
                            Sire.Inbreeding + Dam.Inbreeding +
                            (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                            (1|Management) + (1|Location),
                          data = data_dev_sc,
                          family = "bernoulli",
                          prior = NULL,
                          iter = 15000,
                          control=list(adapt_delta=0.99),
                          cores = 8,
                          sample_prior = TRUE,
                          seed = 13579,
                          save_pars = save_pars(all = TRUE))

bayes_R2(model_dev_dam)
plot(rope(model_dev_dam))

# run reduced model based on undecided terms from full model
describe_posterior(model_dev_dam, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
dev_dam_post <- describe_posterior(model_dev_dam, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_dam_post <- as.data.frame(dev_dam_post)
write.table(format(dev_dam_post, digits=2), file = "./Supp_Info/Tables/02_dev_dam_post.txt", sep = ",", quote = FALSE, row.names = F)

#### 10d) Dam f: Egg viability model cross validation ####

model_dev_dam_red <- brm(Fertile ~ 1 +
                       Sire.age.years + sire_age_sq_sc +
                       Sire.Inbreeding +
                       (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                       (1|Management) + (1|Location),
                     data = data_dev_sc,
                     family = "bernoulli",
                     prior = NULL,
                     iter = 15000,
                     control=list(adapt_delta=0.99),
                     cores = 8,
                     sample_prior = TRUE,
                     seed = 13579,
                     save_pars = save_pars(all = TRUE))

describe_posterior(model_dev_dam_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_dev_dam_red)
plot(rope(model_dev_dam_red))

# save table
dev_dam_post_red <- describe_posterior(model_dev_dam_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_dam_post_red <- as.data.frame(dev_dam_post_red)
write.table(format(dev_dam_post_red, digits=2), file = "./Supp_Info/Tables/10_dev_dam_post_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
dev_dam_loo <- loo(model_dev_dam, moment_match = TRUE, recompile = TRUE, 
                   save_psis = TRUE, cores = 8, reloo = TRUE)
dev_dam_red_loo <- loo(model_dev_dam_red, moment_match = TRUE, recompile = TRUE, 
                       save_psis = TRUE, cores = 8, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(dev_dam_loo, dev_dam_red_loo)
                 
#### 11)  Egg viability model plots ####
#### 11a) Pair k (orange) ####

# Extract data from model
plot_stats_dev_pair <- as.data.frame(as.matrix(model_dev_pair)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l", "z_", "Intercept")))

# Add median and HDI to new df
plot_stats_dev_pair_long <- pivot_longer(plot_stats_dev_pair 
                                         %>% select(-b_Intercept), 
                                         cols=everything(), 
                                         names_to='variable')

# Retain order of variables for plot:
plot_stats_dev_pair_long$variable <- fct_relevel(plot_stats_dev_pair_long$variable, 
                                                 names(plot_stats_dev_pair))
plot_stats_dev_pair_long$y <- plot_stats_dev_pair_long$variable

fig_2a <- 
  plot(rope(model_dev_pair, range=rope_range(model_dev_pair))) +
  stat_pointinterval(data=plot_stats_dev_pair_long, 
                     aes(y=factor(variable,levels=names(plot_stats_dev_pair)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=1.5, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-13,13))+ # scale out to compare to other graph
  labs(title='(A)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Oranges", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));fig_2a

ggsave("./Figures/fig_2a_dev_pair.png", fig_2a, width = 9, height = 10, units = "cm")

#### 11b) Dam f (red) ####

# Extract data from model
plot_stats_dev_dam <- as.data.frame(as.matrix(model_dev_dam)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l", "z_", "Intercept")))

# Add median and HDI to new df
plot_stats_dev_dam_long <- pivot_longer(plot_stats_dev_dam 
                                    %>% select(-b_Intercept), 
                                    cols=everything(), 
                                    names_to='variable')

# Retain order of variables for plot:
plot_stats_dev_dam_long$variable <- fct_relevel(plot_stats_dev_dam_long$variable, 
                                            names(plot_stats_dev_dam))
plot_stats_dev_dam_long$y <- plot_stats_dev_dam_long$variable

fig_2b <- 
  plot(rope(model_dev_dam, range=rope_range(model_dev_dam))) +
  stat_pointinterval(data=plot_stats_dev_dam_long, 
                     aes(y=factor(variable,levels=names(plot_stats_dev_dam)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=1.5, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))+  
  scale_x_continuous(labels = function(x)x/2, limits = c(-13,13)) + # to convert posterior scale from 0.5 SD to SD
  labs(title='(B)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Reds", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));fig_2b

ggsave("./Figures/fig_2b_dev_dam.png", fig_2b, width = 9, height = 10, units = "cm")

#### 11c) Sire age on viability graph ####

viability_sire_age_graph_data<- 
  ddply(data_dev, c("Sire.age.years", "Sire.Studbook.ID"), 
        summarise, fertile=sum(Fertile=="Y"), 
        not.fertile=sum(Fertile=="N"))

viability_sire_age_graph_data$prop.fertile<- 
  viability_sire_age_graph_data$fertile/
  (viability_sire_age_graph_data$not.fertile+
     viability_sire_age_graph_data$fertile)


###calculate mean & standard deviation across pairs for each year...  
viability_sire_age_mean_viability<- 
  ddply(viability_sire_age_graph_data, "Sire.age.years", 
        summarise, mean.fertile=mean(prop.fertile), sd.fertile=sd(prop.fertile))

###plot 
ggplot(data=viability_sire_age_mean_viability, aes(x=Sire.age.years, y=mean.fertile))+
  geom_point()+
  geom_smooth(col = "Orange", fill = "Orange")+
  theme_classic()+
  labs(y= "Proportion of eggs viable", x = "Paternal age")+
  scale_x_discrete(limits=factor(0:22))+
  scale_y_continuous(limits = c(-0.2,1), breaks = c(0.0,0.25,0.5,0.75,1.0))

################################################################################
#### 12)  Egg viability AI models ####

data_dev_sc_AI <- data_dev_sc %>% filter(Artificial.incubation == "Y" | 
                                           Artificial.incubation == "N")
#### 12a) Pair k ####

model_dev_pair_AI <- brm(Fertile ~ 1 +
                      Dam.age.years + dam_age_sq_sc+
                      Sire.age.years + sire_age_sq_sc +
                      Sire.Inbreeding + Pair.Inbreeding +
                      Artificial.incubation +
                      (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                      (1|Management) + (1|Location),
                    data = data_dev_sc_AI,
                    family = "bernoulli",
                    prior = NULL,
                    iter = 15000,
                    control=list(adapt_delta=0.99),
                    cores = 8,
                    sample_prior = TRUE,
                    seed = 13579,
                    save_pars = save_pars(all = TRUE))

plot(rope(model_dev_pair_AI))
bayes_R2(model_dev_pair_AI)

# run reduced model based on undecided terms from full model
describe_posterior(model_dev_pair_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
dev_pair_post_ai <- describe_posterior(model_dev_pair_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_pair_post_ai <- as.data.frame(dev_pair_post_ai)
write.table(format(dev_pair_post_ai, digits=2), file = "./Supp_Info/Tables/03_dev_pair_post_ai.txt", sep = ",", quote = FALSE, row.names = F)

#### 12b) Pair k: Egg viability AI model cross validation; refit once each ####

model_dev_pair_AI_red <- brm(Fertile ~ 1 +
                           Sire.age.years + sire_age_sq_sc +
                           Sire.Inbreeding + Pair.Inbreeding +
                           (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                           (1|Management) + (1|Location),
                         data = data_dev_sc_AI,
                         family = "bernoulli",
                         prior = NULL,
                         iter = 15000,
                         control=list(adapt_delta=0.99),
                         cores = 8,
                         sample_prior = TRUE,
                         seed = 13579,
                         save_pars = save_pars(all = TRUE))

describe_posterior(model_dev_pair_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_dev_pair_AI_red)
plot(rope(model_dev_pair_AI_red))

# save table
dev_pair_post_ai_red <- describe_posterior(model_dev_pair_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_pair_post_ai_red <- as.data.frame(dev_pair_post_ai_red)
write.table(format(dev_pair_post_ai_red, digits=2), file = "./Supp_Info/Tables/11_dev_pair_post_ai.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
dev_pair_AI_loo <- loo(model_dev_pair_AI, moment_match = TRUE, recompile = TRUE,
                       save_psis = TRUE, cores = 8, reloo = TRUE)
dev_pair_AI_red_loo <- loo(model_dev_pair_AI_red, moment_match = TRUE, recompile = TRUE, 
                           save_psis = TRUE, cores = 8, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(dev_pair_AI_loo, dev_pair_AI_red_loo)

#### 12c) Dam f ####

model_dev_dam_AI <- brm(Fertile ~ 1 +
                      Dam.age.years + dam_age_sq_sc+
                      Sire.age.years + sire_age_sq_sc +
                      Sire.Inbreeding + 
                      Dam.Inbreeding +
                      Artificial.incubation +
                      (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                      (1|Management) + (1|Location),
                      data = data_dev_sc_AI,
                      family = "bernoulli",
                      prior = NULL,
                      iter = 15000,
                      control=list(adapt_delta=0.99),
                      cores = 8,
                      sample_prior = TRUE,
                      seed = 13579,
                      save_pars = save_pars(all = TRUE))

plot(rope(model_dev_dam_AI))
bayes_R2(model_dev_dam_AI)

# run reduced model based on undecided terms from full model
describe_posterior(model_dev_dam_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
dev_dam_post_ai <- describe_posterior(model_dev_dam_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_dam_post_ai <- as.data.frame(dev_dam_post_ai)
write.table(format(dev_dam_post_ai, digits=2), file = "./Supp_Info/Tables/04_dev_dam_post_ai.txt", sep = ",", quote = FALSE, row.names = F)

#### 12d) Dam f: Egg viability AI model cross validation; full refit once ####

model_dev_dam_AI_red <- brm(Fertile ~ 1 +
                               Sire.age.years + sire_age_sq_sc +
                               (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                               (1|Management) + (1|Location),
                             data = data_dev_sc_AI,
                             family = "bernoulli",
                             prior = NULL,
                             iter = 15000,
                             control=list(adapt_delta=0.99),
                             cores = 8,
                             sample_prior = TRUE,
                             seed = 13579,
                             save_pars = save_pars(all = TRUE))

describe_posterior(model_dev_dam_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_dev_dam_AI_red)
plot(rope(model_dev_dam_AI_red))

# save table
dev_dam_post_ai_red <- describe_posterior(model_dev_dam_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
dev_dam_post_ai_red <- as.data.frame(dev_dam_post_ai_red)
write.table(format(dev_dam_post_ai_red, digits=2), file = "./Supp_Info/Tables/12_dev_dam_post_ai_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
dev_dam_AI_loo <- loo(model_dev_dam_AI, moment_match = TRUE, 
                      recompile = TRUE, save_psis = TRUE, reloo = TRUE, cores = 8)
dev_dam_AI_red_loo <- loo(model_dev_dam_AI_red, moment_match = TRUE, 
                          recompile = TRUE, save_psis = TRUE, reloo = TRUE, cores = 8)

## compare models using expected log predictive density (ELPD)
loo_compare(dev_dam_AI_loo, dev_dam_AI_red_loo)

#### 13)  Egg viability AI plots ####
#### 13a) Pair k (green) ####
# Extract data from model
plot_stats_dev_pair_AI <- as.data.frame(as.matrix(model_dev_pair_AI)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l","z_", "Intercept")))

# Add median and HDI to new df
plot_stats_dev_pair_AI_long <- pivot_longer(plot_stats_dev_pair_AI 
                                       %>% select(-b_Intercept), 
                                       cols=everything(), 
                                       names_to='variable')

# Retain order of variables for plot:
plot_stats_dev_pair_AI_long$variable <- fct_relevel(plot_stats_dev_pair_AI_long$variable, 
                                               names(plot_stats_dev_pair_AI))
plot_stats_dev_pair_AI_long$y <- plot_stats_dev_pair_AI_long$variable

sup_fig_2a <- 
  plot(rope(model_dev_pair_AI, range=rope_range(model_dev_pair_AI))) +
  stat_pointinterval(data=plot_stats_dev_pair_AI_long, 
                     aes(y=factor(variable,levels=names(plot_stats_dev_pair_AI)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=2, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-16,16))+ # scale out to compare to other graph
  labs(title='(A)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Greens", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));sup_fig_2a

ggsave("./Supp_Info/Figures/sup_fig_2a_dev_pair_AI.png", sup_fig_2a, width = 9, height = 10, units = "cm")

#### 13b) Dam f (blue) ####

# Extract data from model
plot_stats_dev_dam_AI <- as.data.frame(as.matrix(model_dev_dam_AI)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l", "z_", "Intercept")))

# Add median and HDI to new df
plot_stats_dev_dam_AI_long <- pivot_longer(plot_stats_dev_dam_AI 
                                       %>% select(-b_Intercept), 
                                       cols=everything(), 
                                       names_to='variable')

# Retain order of variables for plot:
plot_stats_dev_dam_AI_long$variable <- fct_relevel(plot_stats_dev_dam_AI_long$variable, 
                                               names(plot_stats_dev_dam_AI))
plot_stats_dev_dam_AI_long$y <- plot_stats_dev_dam_AI_long$variable

sup_fig_2b <- 
  plot(rope(model_dev_dam_AI, range=rope_range(model_dev_dam_AI))) +
  stat_pointinterval(data=plot_stats_dev_dam_AI_long, 
                     aes(y=factor(variable,levels=names(plot_stats_dev_dam_AI)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=2, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-16,16))+ # scale out to compare to other graph
  labs(title='(B)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Blues", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));sup_fig_2b

ggsave("./Supp_Info/Figures/sup_fig_2b_dev_dam_AI.png", sup_fig_2b, width = 9, height = 10, units = "cm")

################################################################################
#### 20)  Hatching success models ####

## Remove records with unknown and negative viability outcomes
data_dev_hatch <- subset(data_dev_sc, !Fertile %in% c("0"))
data_dev_hatch_unsc <- subset(data_dev, !Fertile %in% c("N"))

#### 20a) Pair k ####

model_hatch_pair <- brm(Hatch.tot ~ 1 +
                     Dam.age.years + dam_age_sq_sc +
                     Sire.age.years + sire_age_sq_sc +
                     Sire.Inbreeding + Pair.Inbreeding +
                     (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                     (1|Management) + (1|Location),
                   data = data_dev_hatch,
                   family = "bernoulli",
                   prior = NULL,
                   iter = 15000,
                   control=list(adapt_delta=0.99),
                   cores = 8,
                   sample_prior = TRUE,
                   seed = 13579,
                   save_pars = save_pars(all = TRUE))

plot(rope(model_hatch_pair))
bayes_R2(model_hatch_pair)

# Posterior statistics
describe_posterior(model_hatch_pair, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
hat_pair_post <- describe_posterior(model_hatch_pair, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_pair_post <- as.data.frame(hat_pair_post)
write.table(format(hat_pair_post, digits=2), file = "./Supp_Info/Tables/05_hat_pair_post.txt", sep = ",", quote = FALSE, row.names = F)

#### 20b) Pair k: Hatching success model cross validation ####

model_hatch_pair_red <- brm(Hatch.tot ~ 1 +
                              Dam.age.years +
                              sire_age_sq_sc +
                              Sire.Inbreeding +
                              (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                              (1|Management) + (1|Location),
                        data = data_dev_hatch,
                        family = "bernoulli",
                        prior = NULL,
                        iter = 15000,
                        control=list(adapt_delta=0.99),
                        cores = 8,
                        sample_prior = TRUE,
                        seed = 13579,
                        save_pars = save_pars(all = TRUE))

describe_posterior(model_hatch_pair_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_hatch_pair_red)
plot(rope(model_hatch_pair_red))

# save table
hat_pair_post_red <- describe_posterior(model_hatch_pair_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_pair_post_red <- as.data.frame(hat_pair_post_red)
write.table(format(hat_pair_post_red, digits=2), file = "./Supp_Info/Tables/13_hat_pair_post_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
hatch_pair_loo <- loo(model_hatch_pair, moment_match = TRUE, recompile = TRUE, 
                      save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_pair_red_loo <- loo(model_hatch_pair_red, moment_match = TRUE, recompile = TRUE, 
                          save_psis = TRUE, cores = 8, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(hatch_pair_loo, hatch_pair_red_loo)

#### 20c) Dam f; 1 refit ####

model_hatch_dam <- brm(Hatch.tot ~ 1 +
                     Dam.age.years + dam_age_sq_sc +
                     Sire.age.years + sire_age_sq_sc +
                     Sire.Inbreeding + Dam.Inbreeding +
                     (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                     (1|Management) + (1|Location),
                   data = data_dev_hatch,
                   family = "bernoulli",
                   prior = NULL,
                   iter = 15000,
                   control=list(adapt_delta=0.99),
                   cores = 8,
                   sample_prior = TRUE,
                   seed = 13579,
                   save_pars = save_pars(all = TRUE))

plot(rope(model_hatch_dam))
bayes_R2(model_hatch_dam)

# run reduced model based on undecided terms from full model
describe_posterior(model_hatch_dam, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
hat_dam_post <- describe_posterior(model_hatch_dam, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_dam_post <- as.data.frame(hat_dam_post)
write.table(format(hat_dam_post, digits=2), file = "./Supp_Info/Tables/06_hat_dam_post.txt", sep = ",", quote = FALSE, row.names = F)

#### 20d) Dam f: Hatching success model cross validation ####

model_hatch_dam_red <- brm(Hatch.tot ~ 1 +
                             Dam.age.years + 
                             sire_age_sq_sc +
                             Sire.Inbreeding + 
                             (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                             (1|Management) + (1|Location),
                            data = data_dev_hatch,
                            family = "bernoulli",
                            prior = NULL,
                            iter = 15000,
                            control=list(adapt_delta=0.99),
                            cores = 8,
                            sample_prior = TRUE,
                            seed = 13579,
                            save_pars = save_pars(all = TRUE))

describe_posterior(model_hatch_dam_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_hatch_dam_red)
plot(rope(model_hatch_dam_red))

# save table
hat_dam_post_red <- describe_posterior(model_hatch_dam_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_dam_post_red <- as.data.frame(hat_dam_post_red)
write.table(format(hat_dam_post_red, digits=2), file = "./Supp_Info/Tables/14_hat_dam_post_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
hatch_dam_loo <- loo(model_hatch_dam, moment_match = TRUE, recompile = TRUE, 
                     save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_dam_red_loo <- loo(model_hatch_dam_red, moment_match = TRUE, recompile = TRUE, 
                         save_psis = TRUE, cores = 8, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(hatch_dam_loo, hatch_dam_red_loo)

#### 21)  Hatching success model plots ####

#### 21a) Pair k (orange) ####

# Extract data from model
plot_stats_hatch_pair <- as.data.frame(as.matrix(model_hatch_pair)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l","z_", "Intercept")))

# Add median and HDI to new df
plot_stats_hatch_pair_long <- pivot_longer(plot_stats_hatch_pair 
                                         %>% select(-b_Intercept), 
                                         cols=everything(), 
                                         names_to='variable')

# Retain order of variables for plot:
plot_stats_hatch_pair_long$variable <- fct_relevel(plot_stats_hatch_pair_long$variable, 
                                                 names(plot_stats_hatch_pair))
plot_stats_hatch_pair_long$y <- plot_stats_hatch_pair_long$variable

fig_3a <- 
  plot(rope(model_hatch_pair, range=rope_range(model_hatch_pair))) +
  stat_pointinterval(data=plot_stats_hatch_pair_long, 
                     aes(y=factor(variable,levels=names(plot_stats_hatch_pair)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=2, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-10,10)) + # scale out to compare to other graph
  labs(title='(A)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Oranges", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));fig_3a

ggsave("./Figures/fig_3a_hat_pair.png", fig_3a, width = 9, height = 10, units = "cm")

#### 21b) Dam f (red) ####

# Extract data from model
plot_stats_hatch_dam <- as.data.frame(as.matrix(model_hatch_dam)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l", "z_", "Intercept")))

# Add median and HDI to new df
plot_stats_hatch_dam_long <- pivot_longer(plot_stats_hatch_dam 
                                      %>% select(-b_Intercept), 
                                      cols=everything(), 
                                      names_to='variable')

# Retain order of variables for plot:
plot_stats_hatch_dam_long$variable <- fct_relevel(plot_stats_hatch_dam_long$variable, 
                                              names(plot_stats_hatch_dam))
plot_stats_hatch_dam_long$y <- plot_stats_hatch_dam_long$variable

fig_3b <- 
  plot(rope(model_hatch_dam, range=rope_range(model_hatch_dam))) +
  stat_pointinterval(data=plot_stats_hatch_dam_long, 
                     aes(y=factor(variable,levels=names(plot_stats_hatch_dam)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=2, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-10,10)) + # scale out to compare to other graph
  labs(title='(B)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Reds", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));fig_3b

ggsave("./Figures/fig_3b_hat_dam.png", fig_3b, width = 9, height = 10, units = "cm")

################################################################################
#### 22)  Hatching success AI models ####

data_dev_hatch_sc_AI <- data_dev_hatch %>% filter(Artificial.incubation == "Y" | 
                                                    Artificial.incubation == "N")

# remove eggs from first year females; -0.8 because data has been scaled
data_dev_hatch_sc_AI_nyf <- data_dev_hatch_sc_AI %>% filter(Dam.age.years > -0.8)

#### 22a) Pair k ####
model_hatch_pair_AI <- brm(Hatch.tot ~ 1 +
                     Dam.age.years + dam_age_sq_sc +
                     Sire.age.years + sire_age_sq_sc +
                     Sire.Inbreeding + Pair.Inbreeding +
                     Artificial.incubation +
                     (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                     (1|Management) + (1|Location),
                   data = data_dev_hatch_sc_AI_nyf,
                   family = "bernoulli",
                   prior = NULL,
                   iter = 15000,
                   control=list(adapt_delta=0.99),
                   cores = 8,
                   sample_prior = TRUE,
                   seed = 13579,
                   save_pars = save_pars(all = TRUE))

bayes_R2(model_hatch_pair_AI)
plot(rope(model_hatch_pair_AI))

# run reduced model based on undecided terms from full model
describe_posterior(model_hatch_pair_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
hat_pair_post_ai <- describe_posterior(model_hatch_pair_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_pair_post_ai <- as.data.frame(hat_pair_post_ai)
write.table(format(hat_pair_post_ai, digits=2), file = "./Supp_Info/Tables/07_hat_pair_post_ai.txt", sep = ",", quote = FALSE, row.names = F)

#### 22b) Pair k: hatching success cross validation ####

model_hatch_pair_AI_red <- brm(Hatch.tot ~ 1 +
                                 Dam.age.years + dam_age_sq_sc +
                                 Pair.Inbreeding +
                                 Artificial.incubation +
                                 (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                                 (1|Management) + (1|Location),
                               data = data_dev_hatch_sc_AI_nyf,
                               family = "bernoulli",
                               prior = NULL,
                               iter = 15000,
                               control=list(adapt_delta=0.99),
                               cores = 8,
                               sample_prior = TRUE,
                               seed = 13579,
                               save_pars = save_pars(all = TRUE))

describe_posterior(model_hatch_pair_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_hatch_pair_AI_red)
plot(rope(model_hatch_pair_AI_red))

# save table
hat_pair_post_ai_red <- describe_posterior(model_hatch_pair_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_pair_post_ai_red <- as.data.frame(hat_pair_post_ai_red)
write.table(format(hat_pair_post_ai_red, digits=2), file = "./Supp_Info/Tables/15_hat_pair_post_ai_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
hatch_pair_AI_loo <- loo(model_hatch_pair_AI, moment_match = TRUE, recompile = TRUE, 
                         save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_pair_AI_red_loo <- loo(model_hatch_pair_AI_red, moment_match = TRUE, recompile = TRUE, 
                             save_psis = TRUE, cores = 8, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(hatch_pair_AI_loo, hatch_pair_AI_red_loo)

#### 22c) Dam f ####

model_hatch_dam_AI <- brm(Hatch.tot ~ 1 +
                        Dam.age.years + dam_age_sq_sc +
                        Sire.age.years + sire_age_sq_sc +
                        Sire.Inbreeding + Dam.Inbreeding +
                        Artificial.incubation +
                        (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                        (1|Management) + (1|Location),
                      data = data_dev_hatch_sc_AI_nyf,
                      family = "bernoulli",
                      prior = NULL,
                      iter = 15000,
                      control=list(adapt_delta=0.99),
                      cores = 8,
                      sample_prior = TRUE,
                      seed = 13579,
                      save_pars = save_pars(all = TRUE))

bayes_R2(model_hatch_dam_AI)
plot(rope(model_hatch_dam_AI))

# run reduced model based on undecided terms from full model
describe_posterior(model_hatch_dam_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))

# save table
hat_dam_post_ai <- describe_posterior(model_hatch_dam_AI, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_dam_post_ai <- as.data.frame(hat_dam_post_ai)
write.table(format(hat_dam_post_ai, digits=2), file = "./Supp_Info/Tables/08_hat_dam_post_ai.txt", sep = ",", quote = FALSE, row.names = F)

#### 22d) Dam f: hatching success cross validation ####

model_hatch_dam_AI_red <- brm(Hatch.tot ~ 1 +
                                Dam.age.years + dam_age_sq_sc +
                                sire_age_sq_sc +
                                Artificial.incubation +
                                (1|Sire.Studbook.ID) + (1|Dam.Studbook.ID)+
                                (1|Management) + (1|Location),
                               data = data_dev_hatch_sc_AI_nyf,
                               family = "bernoulli",
                               prior = NULL,
                               iter = 15000,
                               control=list(adapt_delta=0.99),
                               cores = 8,
                               sample_prior = TRUE,
                               seed = 13579,
                               save_pars = save_pars(all = TRUE))

describe_posterior(model_hatch_dam_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
bayes_R2(model_hatch_dam_AI_red)
plot(rope(model_hatch_dam_AI_red))

# save table
hat_dam_post_ai_red <- describe_posterior(model_hatch_dam_AI_red, test = c("p_direction", "p_significance", "rope", 'equivalence_test'))
hat_dam_post_ai_red <- as.data.frame(hat_dam_post_ai_red)
write.table(format(hat_dam_post_ai_red, digits=2), file = "./Supp_Info/Tables/16_hat_dam_post_ai_red.txt", sep = ",", quote = FALSE, row.names = F)

## conduct loo cross validation on both full and reduced model
hatch_dam_AI_loo <- loo(model_hatch_dam_AI, moment_match = TRUE, recompile = TRUE, 
                        save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_dam_AI_red_loo <- loo(model_hatch_dam_AI_red, moment_match = TRUE, recompile = TRUE, 
                            save_psis = TRUE, reloo = TRUE)

## compare models using expected log predictive density (ELPD)
loo_compare(hatch_dam_AI_loo, hatch_dam_AI_red_loo)

#### 23)  Hatching success AI model plots ####
#### 23a) Pair k (green) ####
# Extract data from model
plot_stats_hatch_pair_AI <- as.data.frame(as.matrix(model_hatch_pair_AI)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l", "z_", "Intercept")))

# Add median and HDI to new df
plot_stats_hatch_pair_AI_long <- pivot_longer(plot_stats_hatch_pair_AI 
                                         %>% select(-b_Intercept), 
                                         cols=everything(), 
                                         names_to='variable')

# Retain order of variables for plot:
plot_stats_hatch_pair_AI_long$variable <- fct_relevel(plot_stats_hatch_pair_AI_long$variable, 
                                                 names(plot_stats_hatch_pair_AI))
plot_stats_hatch_pair_AI_long$y <- plot_stats_hatch_pair_AI_long$variable

fig_3c <- 
  plot(rope(model_hatch_pair_AI, range=rope_range(model_hatch_pair_AI))) +
  stat_pointinterval(data=plot_stats_hatch_pair_AI_long, 
                     aes(y=factor(variable,levels=names(plot_stats_hatch_pair_AI)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=2, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-15,15)) + # scale out to compare to other graph
  labs(title='(C)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Greens", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));fig_3c

ggsave("./Figures/fig_3c_hat_pair_AI.png", fig_3c, width = 9, height = 10, units = "cm")

#### 23b) Dam f (blue) ####

# Extract data from model
plot_stats_hatch_dam_AI <- as.data.frame(as.matrix(model_hatch_dam_AI)) %>% 
  select(!starts_with(c("r_", "sd_", "prior", "l", "z_", "Intercept")))

# Add median and HDI to new df
plot_stats_hatch_dam_AI_long <- pivot_longer(plot_stats_hatch_dam_AI 
                                         %>% select(-b_Intercept), 
                                         cols=everything(), 
                                         names_to='variable')

# Retain order of variables for plot:
plot_stats_hatch_dam_AI_long$variable <- fct_relevel(plot_stats_hatch_dam_AI_long$variable, 
                                                 names(plot_stats_hatch_dam_AI))
plot_stats_hatch_dam_AI_long$y <- plot_stats_hatch_dam_AI_long$variable

fig_3d <- 
  plot(rope(model_hatch_dam_AI, range=rope_range(model_hatch_dam_AI))) +
  stat_pointinterval(data=plot_stats_hatch_dam_AI_long, 
                     aes(y=factor(variable,levels=names(plot_stats_hatch_dam_AI)),
                         x=value, height=NULL,fill=NULL), point_interval=median_hdi,
                     .width=c(0.5, 0.95), normalize='xy',point_size=2, point_colour='Black', 
                     interval_colour='Dark Blue', shape=21, point_fill='Orange', alpha=1) +
  scale_y_discrete(labels=fn_labels, expand = expansion(add = c(-0,1.2)))  +
  scale_x_continuous(labels = function(x)x/2, # to convert posterior scale from 0.5 SD to SD
                     limits = c(-15,15)) + # scale out to compare to other graph
  labs(title='(D)', y = 'Variable', x = "Possible values (SD)")  +
  guides(fill = 'none')+
  theme_classic()+
  scale_fill_brewer(palette = "Blues", direction = -1)+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  theme(axis.text.y = element_markdown())+
  theme(plot.title = element_text(hjust = 0.05, vjust = -4, size = 15));fig_3d

ggsave("./Figures/fig_3d_hat_dam_AI.png", fig_3d, width = 9, height = 10, units = "cm")

#### 23c) Dam age and hatching success Ai graph ####

data_dev_hatch_AI_unsc <- data_dev %>% filter(Artificial.incubation == "Y" | 
                                                Artificial.incubation == "N",
                                                Fertile == "Y",
                                              Dam.age.years > 1)

dam_age_hatch_graph_data_ai<- 
  ddply(data_dev_hatch_AI_unsc, c("Dam.age.years", "Dam.Studbook.ID"), 
        summarise, hatch=sum(Hatch.tot=="Y"), 
        not.hatch=sum(Hatch.tot=="N"))

dam_age_hatch_graph_data_ai$prop.hatch<- 
  dam_age_hatch_graph_data_ai$hatch/(dam_age_hatch_graph_data_ai$not.hatch+dam_age_hatch_graph_data_ai$hatch)

###calculate mean & standard deviation across pairs for each year...  
dam.mean.hatch.age_ai<- ddply(dam_age_hatch_graph_data_ai, "Dam.age.years", 
                              summarise, mean.hatch=mean(prop.hatch), sd.hatch=sd(prop.hatch))

ggplot(data=dam.mean.hatch.age_ai, aes(x=Dam.age.years, y=mean.hatch))+
  geom_point()+
  geom_smooth(col = "Dark Green", fill = "Light Green")+
  theme_classic()+
  labs(y= "Proportion of eggs hatched", x = "Maternal age")+
  scale_x_continuous(limits=c(2,9), breaks = c(2,3,4,5,6,7,8,9))+
  scale_y_continuous(limits = c(0,1.3), breaks = c(0.0,0.25,0.5,0.75,1.0))

##%# End #####



#### Quick loo section ####

## conduct loo cross validation on both full and reduced model
dev_pair_loo <- loo(model_dev_pair, moment_match = TRUE, recompile = TRUE, 
                    save_psis = TRUE, reloo = TRUE)
dev_pair_red_loo <- loo(model_dev_pair_red, moment_match = TRUE, recompile = TRUE, 
                        save_psis = TRUE, reloo = TRUE)

## conduct loo cross validation on both full and reduced model
dev_dam_loo <- loo(model_dev_dam, moment_match = TRUE, recompile = TRUE, 
                   save_psis = TRUE, cores = 8, reloo = TRUE)
dev_dam_red_loo <- loo(model_dev_dam_red, moment_match = TRUE, recompile = TRUE, 
                       save_psis = TRUE, cores = 8, reloo = TRUE)

## conduct loo cross validation on both full and reduced model
dev_pair_AI_loo <- loo(model_dev_pair_AI, moment_match = TRUE, recompile = TRUE,
                       save_psis = TRUE, cores = 8, reloo = TRUE)
dev_pair_AI_red_loo <- loo(model_dev_pair_AI_red, moment_match = TRUE, recompile = TRUE, 
                           save_psis = TRUE, cores = 8, reloo = TRUE)

## conduct loo cross validation on both full and reduced model
dev_dam_AI_loo <- loo(model_dev_dam_AI, moment_match = TRUE, recompile = TRUE, 
                      save_psis = TRUE, reloo = TRUE, cores = 8)
dev_dam_AI_red_loo <- loo(model_dev_dam_AI_red, moment_match = TRUE, recompile = TRUE, 
                          save_psis = TRUE, reloo = TRUE, cores = 8)

## conduct loo cross validation on both full and reduced model
hatch_pair_loo <- loo(model_hatch_pair, moment_match = TRUE, recompile = TRUE, 
                      save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_pair_red_loo <- loo(model_hatch_pair_red, moment_match = TRUE, recompile = TRUE, 
                          save_psis = TRUE, cores = 8, reloo = TRUE)

## conduct loo cross validation on both full and reduced model
hatch_dam_loo <- loo(model_hatch_dam, moment_match = TRUE, recompile = TRUE, 
                     save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_dam_red_loo <- loo(model_hatch_dam_red, moment_match = TRUE, recompile = TRUE, 
                         save_psis = TRUE, cores = 8, reloo = TRUE)

## conduct loo cross validation on both full and reduced model
hatch_pair_AI_loo <- loo(model_hatch_pair_AI, moment_match = TRUE, recompile = TRUE, 
                         save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_pair_AI_red_loo <- loo(model_hatch_pair_AI_red, moment_match = TRUE, recompile = TRUE, 
                             save_psis = TRUE, cores = 8, reloo = TRUE)

## conduct loo cross validation on both full and reduced model
hatch_dam_AI_loo <- loo(model_hatch_dam_AI, moment_match = TRUE, recompile = TRUE, 
                        save_psis = TRUE, cores = 8, reloo = TRUE)
hatch_dam_AI_red_loo <- loo(model_hatch_dam_AI_red, moment_match = TRUE, recompile = TRUE, 
                            save_psis = TRUE, reloo = TRUE)

## all compares
loo_compare(dev_pair_loo, dev_pair_red_loo)
loo_compare(dev_dam_loo, dev_dam_red_loo)
loo_compare(dev_pair_AI_loo, dev_pair_AI_red_loo)
loo_compare(dev_dam_AI_loo, dev_dam_AI_red_loo)

loo_compare(hatch_pair_loo, hatch_pair_red_loo)
loo_compare(hatch_dam_loo, hatch_dam_red_loo)
loo_compare(hatch_pair_AI_loo, hatch_pair_AI_red_loo)
loo_compare(hatch_dam_AI_loo, hatch_dam_AI_red_loo)


#### Quick r2 section ####

bayes_R2(model_dev_pair)
bayes_R2(model_dev_dam)
bayes_R2(model_dev_pair_AI)
bayes_R2(model_dev_dam_AI)

bayes_R2(model_hatch_pair)
bayes_R2(model_hatch_dam)
bayes_R2(model_hatch_pair_AI)
bayes_R2(model_hatch_dam_AI)