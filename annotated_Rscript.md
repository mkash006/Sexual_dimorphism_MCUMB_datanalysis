# Loading required libraries
`ggplot2`, `emmeans`, `lme4`, and `dplyr` are the core libraries that we will use for calculations and plotting. The remaining libraries are used to format tables and figures for the publication associated with this dataset. 

```R
library("ggplot2")
library("emmeans")
library("lme4")
library("dplyr")
library("tibble")
library("writexl")
library("officer")
library("flextable")
```

# Description of the dataset
'cdsd_dat.csv' contains the following columns:

1. ID: Each entry in the ID column represents a unique identifier for an individual sample collected during the experiment. The ID codes follow a standardized naming convention that encodes metadata about the sample.
2. Block: Statistical replicates that represent shared ancestry between population pairs belonging to the two alternative selection regimes. 
3.  Selection: Categorical variable encoding data for the selection regime to which each morphological measurement belongs.
4. Treatment: Categorical variable describing the density treatment used during larval rearing.
5. Sex: Categorical variable describing the sex of each adult fly  
Column 6-11 contains morphological measurements for the trait listed in the column header for the same adult fly.

```R
dat<-read.csv("cdsd_dat.csv")
```

# Calculating and plotting principal components

The following code chunk uses `prcomp` (a base R function) to calculate principal components using a single correlation matrix (`scale = TRUE`), combining all the selection regimes, density treatments, sexes, and blocks for scaling. Calculated principal components are input into `ggplot` to create some useful visualizations of the dataset 

```R
# PCA using a single correlation matrix 
mypr_combined <- dat %>%            
  select(6:10) %>%                     
  prcomp(scale. = TRUE)                
# Flextable for eigenvalues and eigenvectors (TBA to supplementary materials)  
combinedanalysis_PCtable = pca_table(mypr_combined)

# figure 1 
combinedanalysis_plotdat = cbind(dat,mypr_combined$x)
plot_figure1 <- ggplot(data = combinedanalysis_plotdat ,
                 aes(Sex, PC1)) +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",size = 0.5, width = 0.3) +
  stat_summary(fun = "mean", size = 0.2) +
  ylab("PC1 (95% CIs)") + xlab("Treatment")  +
  theme_minimal(base_size = 18) + theme_bw()  +theme(axis.text.x = element_text(size = 12)) + facet_wrap(~Treatement)

# supplementary figure 
# plot illustrates the correlation matrix used in the combined analysis 
dat_with_pc1 <- dat %>%
    select(6:10) %>%
    mutate(PC1 = mypr_combined$x[, 1])
  dat_with_pc1 %>%
    cor(use = "complete.obs") %>%
    melt() %>%
    ggplot(aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    coord_fixed() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3)
```


# Body size correction and linear mixed effects modelling


The following code chunk uses custom functions (code annotations provided in the Custom functions section) to calculate residuals, which can be interpreted as body size corrected measures of all the morphological variables that are in our dataset. Another custom function call generates faceted plots similar to Figure 1 for each body size-corrected morphological variable

Code to generate a linear mixed effects model for composite body size (principal component 1) has also been provided along with the `emmeans` function to generate a multiple comparison table for the same. 

```R
# lmer and multiple comparisons for composite body size (principal component 1)
lmmcombined<-lmer(PC1~Treatement*Sex*Selection+(1|Block), combinedanalysis_plotdat)
anovatable<-car::Anova(lmmcombined, type=3)

multiple_compar <- emmeans(lmmcombined, pairwise ~ Sex*Selection|Treatement)
combinedlm_multiplectable<-summary(multiple_compar)
multiple_compar_df <- as.data.frame(summary(multiple_compar))

# dplyr pipe to remove the outlier MlfF1, which is an imagej scaling error, and perform multiple regressions between PC1 (calculated excluding the focal trait) and calculate residuals
traitwise_pc1_regressions <- dat %>%
    dplyr::filter(ID != "MLfF1") %>%
    body_size_correction(c("Femur","LW","RW","Thorax","Tibia","Wings")) # custom function to calculate residuals, which are essentially body size corrected measures since pc1 is composite body size

# using a custom function made to fit linear mixed effects models for body size-corrected morphological measurements
resid_data = traitwise_pc1_regressions$corrected_data
traitwise_lmer_pc1regressions <- run_traitwise_lmer(resid_data, c("Femur","LW","RW","Thorax","Tibia","Wings"))

# calling the plotting function to plot and means+-cis of calculated residuals
resid_data %>%
  mutate(Selection = if_else(Selection == "CU", "MCU", Selection)) %>%
  run_traitwise_plots(c("Femur","LW","RW","Thorax","Tibia","Wings")) -> traitwise_meanci_plots
combined_meanciplots <- wrap_plots(traitwise_meanci_plots, ncol = 2)
```


# Custom functions

The following contains annotated code for all the custom functions that are used in the main script. 

```R
#function to make eigenvalue and eigenvector tables. Takes a prcomp object as input and returns a data table of the proportion of variation explained by each PC. This table is used in a flextable pipeline to print a docx file for principal component loadings 
pca_table <- function(pr_obj){
  # % variance explained
  var_pct <- 100 * (pr_obj$sdev^2 / sum(pr_obj$sdev^2))
  
  # loadings matrix
  loadings <- pr_obj$rotation
  
  # assemble data frame
  df <- data.frame(Trait = rownames(loadings),
                   round(loadings, 3),
                   check.names = FALSE)
  colnames(df)[-1] <- paste0("PC", seq_along(var_pct))
  
  # add % variance as first row
  df <- bind_rows(
    data.frame(Trait = "% Variance", t(round(var_pct, 1))),
    df
  )
  df
}
# function to perform body size correction using pc1 (calculated excluding focal trait) - trait regression residuals. Produces regression plots and dataframe with residuals sotred as output
 body_size_correction <- function(data, trait_vars) {
    plot_list <- list()
    
    for (trait in trait_vars) {
      
      # Traits to include in PCA (exclude focal trait)
      traits_for_pca <- setdiff(trait_vars, trait)
      
      # Run PCA on correlation matrix
      pca <- prcomp(data[, traits_for_pca], scale. = TRUE)
      
      # Extract PC1 as composite body size
      data$PC1_composite <- pca$x[, 1]
      
      # Linear regression of focal trait ~ PC1
      mod <- lm(data[[trait]] ~ PC1_composite, data = data)
      
      # Store residuals as body-size-corrected trait
      data[[paste0(trait, "_resid")]] <- residuals(mod)
      
      # Build regression plot for this trait
      p <- ggplot(data, aes(x = PC1_composite, y = .data[[trait]])) +
        geom_point(alpha = 0.6, size = 0.1) +
        geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
        labs(
          title = paste("Body-size correction for", trait),
          x = paste("PC1 excluding", trait),
          y = trait
        ) +
        theme_minimal()
      
      plot_list[[trait]] <- p
    }   # ← FIX #1: closes the for-loop
    
    # Combine all regression panels into one plot
    final_plot <- wrap_plots(plot_list, ncol = 2)
    
    return(list(
      corrected_data = data,
      regression_panel = final_plot
    ))
  }     

# Function to perform mixed effect models in calculated residuals
run_traitwise_lmer <- function(corrected_output, trait_vars) {

  dat <- corrected_output
  resid_vars <- paste0(trait_vars, "_resid")

  model_list <- list()
  anova_list <- list()

  for (v in resid_vars) {

    fm <- as.formula(paste(v, "~ Treatement * Sex * Selection + (1|Block)"))

    mod <- lmer(fm, data = dat)
    model_list[[v]] <- mod

    anova_list[[v]] <- car::Anova(mod, type = 3)
  }

  return(list(
    models = model_list,
    anova_tables = anova_list
  ))
}

# function to plot means and confidence intervals of calculated residuals 
run_traitwise_plots <- function(corrected_output, trait_vars) {
  
  dat <- corrected_output
  resid_vars <- paste0(trait_vars, "_resid")
  
  plot_list <- list()
  
  for (v in resid_vars) {
    
    p <- ggplot(dat, aes(x = Selection, y = .data[[v]], shape = Sex)) +
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   size = 0.5, width = 0.3) +
      stat_summary(fun = "mean", size = 0.2) +
      ylab(paste0(v, " (95% CIs)")) +
      xlab("Sex") +
      theme_minimal(base_size = 18) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 12)) +
      facet_wrap(~Treatement)
    
    plot_list[[v]] <- p
  }
  
  return(plot_list)
}


```