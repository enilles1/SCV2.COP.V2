# - - - - - - - - - - - - - - - - - - - - - - - 
# Main code for COP estimation
# - - - - - - - - - - - - - - - - - - - - - - -

library(readr)
library(visreg)
library(mgcv)
library(splines)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(dplyr)

# Import data ------------------------------------------------

#setwd("~/Documents/GitHub/SCV2.COP/")
df_mu <- read_csv("data/df_mu_out.csv")
df_delta <- read_csv("data/df_delta_out.csv")
df_ba1 <- read_csv("data/df_ba1_out.csv")
df_ba2 <- read_csv("data/df_ba2_out.csv")
df_ba45 <- read_csv("data/df_ba45_out.csv")


# Model – run with each variant specific df -------------------------------

model_gam <- function(data_frame){
  # Format data
  data_frame$C19_pos3 <- as.factor(data_frame$C19_pos3)
  data_frame$vacc_dose_orig <- as.factor(data_frame$vacc_dose_orig)
  data_frame$setting <- as.factor(data_frame$setting)
  data_frame$gender_original <- as.factor(data_frame$gender_original)
  data_frame$nbr_hh_cat <- as.factor(data_frame$nbr_hh_cat)
  
  # Run model
  mgcv::gam(C19_pos3 ~  ns(s_titer_geom, df=2) + 
                              ns(PI_age_calc, df=4) + 
                              ns(date3, df=6) + 
                              vacc_dose_orig +  
                              dpso + 
                              ns(dpv, df=3) + 
                              setting + 
                              SO_location + 
                              nbr_hh_cat + 
                              gender_original, 
                              data = data_frame, 
                              family = "binomial")
}

# Generate models ------------------------------------------------------------

model_mu <- model_gam(df_mu)
model_delta <- model_gam(df_delta)
model_ba1 <- model_gam(df_ba1)
model_ba2 <- model_gam(df_ba2)
model_ba45 <- model_gam(df_ba45)

# XX EXPORT MODELS HERE, TO BE IMPORTED IN FINAL VERSION? XX

# Plot outputs ------------------------------------------------------------

xx1 <- expression(10^0); xx2 <- expression(10^1); xx3 <- expression(10^2)
xx4 <- expression(10^3); xx5 <- expression(10^4); xx6 <- expression(10^5)

# Plot function
plot_model <- function(model,title_name,data_frame){
  visreg(model, xvar = "s_titer_geom", data=data_frame, type = 'contrast', trans = exp, cond=list(s_titer_geom=-0.1), ylab = 'Relative risk',  xlab = NULL, line=list(col="darkblue"),overlay = T,  gg=TRUE, nn=1000, partial = F, rug = F)+
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4,5), labels = c(xx1, xx2, xx3, xx4, xx5, xx6) )+
    scale_y_continuous(trans="log10", breaks = c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1,2.5, 5, 10))+
    coord_cartesian(ylim = c(0.01, 10), xlim = c(0.01, 5))+
    theme_classic() +
    ggtitle(title_name)
}

# Plot variants
p_mu <- plot_model(model_mu,"Mu",df_mu)
p_delta <- plot_model(model_delta,"Delta",df_delta)
p_ba1 <- plot_model(model_ba1,"BA.1",df_ba1)
p_ba2 <- plot_model(model_ba2,"BA.2",df_ba2)
p_ba45 <- plot_model(model_ba45,"BA.4/5",df_ba45)

# Plot grid
grid.arrange(p_mu,p_delta,p_ba1,p_ba2,p_ba45, nrow = 1)

# Output figure
dev.copy(pdf,paste("plot.pdf",sep=""),width=14,height=3)
dev.off()
