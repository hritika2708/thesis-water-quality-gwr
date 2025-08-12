


"""**reading the file**"""

install.packages(c("sf","tidyverse","corrplot","miceadds","MuMIn","GWmodel"))

library(sf)
library(tidyverse)

river_data <- st_read("river-quality/Monthly_River_Quality_Data_SDCC.geojson")

"""# EDA"""

summary(river_data)

"""**creating bar graphs of cols**"""

library(ggplot2)
options(repr.plot.width = 12, repr.plot.height =12)

ggplot(river_data, aes(x = C_O_D__mg_l )) +
  geom_bar(fill = "black") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  labs(
    x = "Category",
    y = "Count",
    title = "Frequency of C_O_D__mg_l "
  )

options(repr.plot.width = 12, repr.plot.height =12)

ggplot(river_data, aes(x = Suspended_Solids_mg_l )) +
  geom_bar(fill = "black") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Category",
    y = "Count",
    title = "Frequency of Suspended_Solids_mg_l "
  )

options(repr.plot.width = 12, repr.plot.height =12)

ggplot(river_data, aes(x = Nitrate_mg_l_as_N )) +
  geom_bar(fill = "black") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  labs(
    x = "Category",
    y = "Count",
    title = "Frequency of Nitrate_mg_l_as_N "
      )

"""**finding if other columns of model have < or > symbol**"""

vars <- c(
  "C_O_D__mg_l",
  "Suspended_Solids_mg_l",
  "Nitrate_mg_l_as_N",
  "TON_mg_l_as_N",
  "Total_B_O_D__mg_l"
)

for (var in vars) {

  col_vals <- as.character(river_data[[var]])


  has_symbol <- grepl("[<>]", col_vals)

  if (any(has_symbol, na.rm = TRUE)) {
    matches <- unique(col_vals[has_symbol])

    cat("Column:", var, "\n")
    print(matches)
    cat("\n")
  }
}

"""based on the above graphs we can create indicator columns if needed

**finding missingness in variables**
"""

library(naniar)
vis_miss(river_data)
gg_miss_var(river_data)

"""we will drop variables with more than 1000 missing values as it is abpve the defined missingness threshold 
of the project
and also dropping sampling point because it isnt needed in analysis
"""

drop_vars <- c("Field7","Arsenic_µg_l","Total_Phosphorus_mg_l_as_P","OBJECTID","Y",
               "Sulphate_mg_l","Chloride_mg_l","Total_Hardness_mgCaCO3_l","X","Sample_Number",
               "Total_Alkalinity_mgCaCO3_l","Sampling_Point","Colour__Hazen", "B_O_D__mg_l")
df <- river_data %>% select(-all_of(drop_vars))

"""**converting character to numeric cols**"""

char_cols <- names(df)[sapply(df, is.character)]

for(col in char_cols) {
  new_col <- paste0("num_", col)
  df[[new_col]] <- as.numeric(df[[col]])
}

"""**finding correlation**"""

library(corrplot)
subdf <- sf::st_drop_geometry(df)
num_df <- subdf[ sapply(subdf, is.numeric) ]

corr_matrix <- cor(
  num_df,
  use    = "pairwise.complete.obs",
  method = "pearson"
)

print(corr_matrix)

corrplot(corr_matrix,
         method = "ellipse",
         type   = "upper",
         diag   = FALSE)

"""based on correlation with num COD we pick the model variables

correlation in model vars
"""

subdf <- sf::st_drop_geometry(df)[
        , c("num_C_O_D__mg_l",
                "num_Suspended_Solids_mg_l",
                "num_Nitrate_mg_l_as_N",
                "num_TON_mg_l_as_N",
                "num_Total_B_O_D__mg_l",
                "Conductivity__20_C__µS_cm",
                "pH"),
        drop = FALSE
      ]

corr_matrix <- cor(subdf, use = "pairwise.complete.obs", method = "pearson")

print(corr_matrix)

corrplot(corr_matrix,
         method = "ellipse",
         type   = "upper",
         tl.col = "black",
         addCoef.col = "black",
         diag   = FALSE)

#we drop TON because there is a complete correlation with Nitrate

"""# LR - 1 with all correlated model vars and mice"""

model_vars <- c("num_C_O_D__mg_l",
                "num_Suspended_Solids_mg_l",
                "num_Nitrate_mg_l_as_N",
                "num_Total_B_O_D__mg_l",
                "Conductivity__20_C__µS_cm",
                "pH")
df_mod <- df %>% select(all_of(model_vars))

df_mod <- st_set_geometry(df_mod, NULL)
sapply(df_mod, class)

na_counts <- sapply(df_mod, function(x) sum(is.na(x)))
print(na_counts)

library(mice)
imp <- mice(df_mod, m = 10, method = "cart",seed   = 123, maxit=50)

df1 <- complete(imp, 1)

keep        <- complete.cases(df1[, model_vars])
completed_cc <- df1[keep, ]

model1 <- lm(num_C_O_D__mg_l ~
               num_Suspended_Solids_mg_l +
               num_Nitrate_mg_l_as_N  +
               num_Total_B_O_D__mg_l +
               Conductivity__20_C__µS_cm +
               pH,
             data = completed_cc)

summary(model1)

library(MuMIn)
AICc(model1)

"""# GWR 1 with all correlated vars and mice"""

library(sf)
river_utm <- st_transform(river_data, 2157)

char_cols <- names(river_utm)[sapply(river_utm, is.character)]

for(col in char_cols) {
  new_col <- paste0("num_", col)
  river_utm[[new_col]] <- as.numeric(river_utm[[col]])
}

gutter <- river_utm %>%
  select(num_C_O_D__mg_l, num_Suspended_Solids_mg_l,num_Nitrate_mg_l_as_N,
         num_Total_B_O_D__mg_l, Conductivity__20_C__µS_cm, pH, geometry)
library(sp)
sp_gutter <- as(gutter, "Spatial")

library(GWmodel)
sapply(sp_gutter@data, class)

dat <- sp_gutter@data
md.pattern(dat)

imp <- mice(dat,
            m      = 10,
            method = 'cart',
            maxit  = 50,
            seed   = 123)

completed_dat <- complete(imp, 1)
sp_gutter_imp <- sp_gutter
sp_gutter_imp@data <- completed_dat

library(spgwr)
bw_cv <- gwr.sel(
  num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                        Conductivity__20_C__µS_cm +
                                      pH,
  data    = sp_gutter_imp,
  gweight = gwr.Gauss,
  longlat = FALSE
)

gwr_sp <- gwr(
  num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                      Conductivity__20_C__µS_cm +
                                      pH,
  data      = sp_gutter_imp,
  bandwidth = bw_cv,
  gweight   = gwr.Gauss,
  hatmatrix = TRUE,
  se.fit    = TRUE,
  longlat   = FALSE
)
print(gwr_sp)

"""# LR - 2 with all corr vars and no mice"""

model_vars <- c("num_C_O_D__mg_l",
                "num_Suspended_Solids_mg_l",
                "num_Nitrate_mg_l_as_N",
                "num_Total_B_O_D__mg_l",
                "Conductivity__20_C__µS_cm",
                "pH")
df_mod <- df %>% select(all_of(model_vars))
df_mod <- st_set_geometry(df_mod, NULL)

df_mod_cc <- df_mod[ complete.cases(df_mod[, model_vars]), ]

model_cc <- lm(num_C_O_D__mg_l ~
                 num_Suspended_Solids_mg_l +
                 num_Nitrate_mg_l_as_N +
                 num_Total_B_O_D__mg_l +
                 Conductivity__20_C__µS_cm +
                 pH,
               data = df_mod_cc)

summary(model_cc)

library(MuMIn)
MuMIn::AICc(model_cc)

"""# GWR - 2 with all corr vars and no mice"""

ok <- complete.cases(sp_gutter@data[, c("num_C_O_D__mg_l",
                                         "num_Suspended_Solids_mg_l",
                                         "num_Nitrate_mg_l_as_N",
                                         "num_Total_B_O_D__mg_l",
                                         "Conductivity__20_C__µS_cm",
                                         "pH")])
sp_gutter_clean <- sp_gutter[ok, ]

coords_clean <- sp::coordinates(sp_gutter_clean)
bw_cv <- gwr.sel(
  num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                        Conductivity__20_C__µS_cm +
                                      pH,
  data    = sp_gutter_clean,      
  gweight = gwr.Gauss,
  longlat = FALSE
)

gwr_sp <- gwr(
  num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                      Conductivity__20_C__µS_cm +
                                      pH,
  data      = sp_gutter_clean,
  bandwidth = bw_cv,
  gweight   = gwr.Gauss,
  hatmatrix = TRUE,
  se.fit    = TRUE,
  longlat   = FALSE
)
print(gwr_sp)

"""# plotting pairwise missingness"""

library(VIM)

vars <- c(
  "num_C_O_D__mg_l",
  "num_Suspended_Solids_mg_l",
  "num_Nitrate_mg_l_as_N",
  "num_Total_B_O_D__mg_l",
  "Conductivity__20_C__µS_cm",
  "pH"
)

sub_df <- sf::st_drop_geometry(df)[, vars]

pdf("pairwise_marginplots.pdf", width=7, height=7)

n <- length(vars)
for(i in seq_len(n-1)) {
  for(j in (i+1):n) {
    xvar <- vars[i]
    yvar <- vars[j]
    marginplot(
      sub_df[, c(xvar, yvar)],
      main = paste0(xvar, "  vs  ", yvar),
      xlab = xvar,
      ylab = yvar
    )
  }
}

dev.off()

"""# creating indicator variables

**COD indicator**
"""

df$i_lt_10 <- ifelse(

  grepl("^<\\s*10$", df$C_O_D__mg_l)
  |

  (grepl("^[0-9]+$", df$C_O_D__mg_l) & as.integer(df$C_O_D__mg_l) < 10),
  0L,
  1L
)

"""**Suspended solids indicator**"""

df$i_lt_5 <- ifelse(

  grepl("^<\\s*5$", df$Suspended_Solids_mg_l)
  |

  (grepl("^[0-9]+$", df$Suspended_Solids_mg_l) & as.integer(df$Suspended_Solids_mg_l) < 5),
  0L,
  1L
)

"""**Nitrate_mg_l_as_N indicator**"""

df$i_lt_01 <- with(df, ifelse(
  grepl("^<\\s*0\\.10$", Nitrate_mg_l_as_N)
    | (!is.na(suppressWarnings(as.numeric(Nitrate_mg_l_as_N)))
       & suppressWarnings(as.numeric(Nitrate_mg_l_as_N)) < 0.10),
  0L, 1L
))

"""# LR 3 - with indicator and all corr vars"""

model_vars <- c("num_C_O_D__mg_l",
                "num_Suspended_Solids_mg_l",
                "num_Nitrate_mg_l_as_N",

                "num_Total_B_O_D__mg_l",
                "Conductivity__20_C__µS_cm",
                "pH",
                "i_lt_5",
                "i_lt_10",
                "i_lt_01"
                )
df_mod <- df %>% select(all_of(model_vars))
df_mod <- st_set_geometry(df_mod, NULL)
sapply(df_mod, class)

library(mice)
imp <- mice(df_mod, m = 10, method = "pmm")

fit_list <- with(imp,
                 lm(num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                        Conductivity__20_C__µS_cm +
                                      pH + i_lt_5 + i_lt_10 +
                                       i_lt_01))

pooled <- pool(fit_list)
summary(pooled)

library(miceadds)

pooled_r2 <- pool.r.squared(fit_list)
print(pooled_r2)

library(car)
df <- complete(imp,1)
lm_cc <- lm(
  num_C_O_D__mg_l ~
    num_Suspended_Solids_mg_l +
    num_Nitrate_mg_l_as_N +
    num_Total_B_O_D__mg_l +
    Conductivity__20_C__µS_cm +
    pH +
    i_lt_5 +
    i_lt_10 +
    i_lt_01,
  data = df
)
summary(lm_cc)
vif_vals <- vif(lm_cc)
print(vif_vals)

library(MuMIn)
print(MuMIn::AICc(lm_cc))

"""# GWR 3 - with indicator and all corr vars

**creating indicator vars**
"""

river_utm$i_lt_10 <- ifelse(

  grepl("^<\\s*10$", river_utm$C_O_D__mg_l)
  |

  (grepl("^[0-9]+$", river_utm$C_O_D__mg_l) & as.integer(river_utm$C_O_D__mg_l) < 10),
  0L,
  1L
)

river_utm$i_lt_5 <- ifelse(

  grepl("^<\\s*5$", river_utm$Suspended_Solids_mg_l)
  |

  (grepl("^[0-9]+$", river_utm$Suspended_Solids_mg_l) & as.integer(river_utm$Suspended_Solids_mg_l) < 5),
  0L,
  1L
)

river_utm$i_lt_01 <- with(river_utm, ifelse(
  grepl("^<\\s*0\\.10$", Nitrate_mg_l_as_N)
    | (!is.na(suppressWarnings(as.numeric(Nitrate_mg_l_as_N)))
       & suppressWarnings(as.numeric(Nitrate_mg_l_as_N)) < 0.10),
  0L, 1L
))

gutter <- river_utm %>%
  select(num_C_O_D__mg_l, num_Suspended_Solids_mg_l,num_Nitrate_mg_l_as_N,
         num_Total_B_O_D__mg_l, Conductivity__20_C__µS_cm, pH,
         i_lt_5, i_lt_10,i_lt_01, geometry)
library(sp)
sp_gutter <- as(gutter, "Spatial")

library(GWmodel)
sapply(sp_gutter@data, class)
ok <- complete.cases(sp_gutter@data[, c("num_C_O_D__mg_l",
                                         "num_Suspended_Solids_mg_l",
                                         "num_Nitrate_mg_l_as_N",

                                         "num_Total_B_O_D__mg_l",
                                         "Conductivity__20_C__µS_cm",
                                         "pH",
                                         "i_lt_5",
                                         "i_lt_10",
                                         "i_lt_01")])
sp_gutter_clean <- sp_gutter[ok, ]

coords_clean <- sp::coordinates(sp_gutter_clean)
library(spgwr)

"""**checking duplicate coords**"""

sum( duplicated(coordinates(sp_gutter_clean)) )

"""**finding collinearity in the GWR Model**"""

fmla   = num_C_O_D__mg_l ~
                num_Suspended_Solids_mg_l +
                num_Nitrate_mg_l_as_N +
                num_Total_B_O_D__mg_l +
                Conductivity__20_C__µS_cm +
                pH +
                i_lt_5 +
                i_lt_10 +
                i_lt_01

X  <- model.matrix(fmla, data = sp_gutter_clean@data)

rk <- qr(X)$rank
nc <- ncol(X)
cat("Design matrix has rank", rk, "but", nc, "columns.\n")

lm_cc <- lm(fmla, data = sp_gutter_clean@data)

print(coef(lm_cc))

print(alias(lm_cc))

uv <- sapply(sp_gutter_clean@data[, all.vars(fmla)[-1]], function(x) length(unique(x)))
print(uv)

"""**we can see that the indicator values have only 1 unique value so they add collinearity**

# GWR 3 with MICE and indicator
"""

dat <- sp_gutter@data
md.pattern(dat)

imp <- mice(dat,
            m      = 10,
            method = 'pmm',
            maxit  = 50,
            seed   = 123)

completed_dat <- complete(imp, 1)
sp_gutter_imp <- sp_gutter
sp_gutter_imp@data <- completed_dat

library(spgwr)
bw_cv <- gwr.sel(
  num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                        Conductivity__20_C__µS_cm +
                                      pH + i_lt_5 +
                i_lt_10 +
                i_lt_01,
  data    = sp_gutter_imp,
  gweight = gwr.Gauss,
  longlat = FALSE
)

gwr_sp <- gwr(
  num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                      Conductivity__20_C__µS_cm +
                                      pH + i_lt_5 +
                i_lt_10 +
                i_lt_01,
  data      = sp_gutter_imp,
  bandwidth = bw_cv,
  gweight   = gwr.Gauss,
  hatmatrix = TRUE,
  se.fit    = TRUE,
  longlat   = FALSE
)
print(gwr_sp)



"""# LR 4 with COD,nitrate, and ph imputation"""

library(sf)
library(mice)
library(dplyr)

df_imp <- df %>%
  st_drop_geometry() %>%
  select(
    num_C_O_D__mg_l,
    pH,
    num_Nitrate_mg_l_as_N,
    Conductivity__20_C__µS_cm,
      num_Suspended_Solids_mg_l,
      num_Total_B_O_D__mg_l
  )

pred_mat <- make.predictorMatrix(df_imp)
pred_mat[,] <- 0
targets <- c("num_C_O_D__mg_l", "pH", "num_Nitrate_mg_l_as_N")
pred_mat[targets, "Conductivity__20_C__µS_cm"] <- 1

meth <- make.method(df_imp)
meth["Conductivity__20_C__µS_cm"] <- ""
meth["num_Suspended_Solids_mg_l"] <- ""
meth["num_Total_B_O_D__mg_l"] <- ""
meth[targets] <- "pmm"

set.seed(123)
imp <- mice(
  data            = df_imp,
  m               = 10,
  predictorMatrix = pred_mat,
  method          = meth,
  printFlag       = TRUE
)

"""**fit 1 complete case model**"""

model_vars <- c("num_C_O_D__mg_l",
                "num_Suspended_Solids_mg_l",
                "num_Nitrate_mg_l_as_N",

                "num_Total_B_O_D__mg_l",
                "Conductivity__20_C__µS_cm",
                "pH")

df1 <- complete(imp, 1)

keep        <- complete.cases(df1[, model_vars])
completed_cc <- df1[keep, ]


model1 <- lm(num_C_O_D__mg_l ~
               num_Suspended_Solids_mg_l +
               num_Nitrate_mg_l_as_N  +
               num_Total_B_O_D__mg_l +
               Conductivity__20_C__µS_cm +
               pH,
             data = completed_cc)


summary(model1)

library(MuMIn)
MuMIn::AICc(model1)

"""# GWR 4 with imputations"""

library(sf)
library(mice)
library(spgwr)

dat <- sp_gutter@data

pred <- make.predictorMatrix(dat)
pred[,] <- 0
targets <- c("num_C_O_D__mg_l", "pH", "num_Nitrate_mg_l_as_N")
pred[targets, "Conductivity__20_C__µS_cm"] <- 1

meth       <- make.method(dat)
meth[]     <- ""
meth[targets] <- "pmm"

set.seed(123)
imp <- mice(
  data            = dat,
  m               = 10,
  predictorMatrix = pred,
  method          = meth,
  printFlag       = TRUE
)

completed_df <- complete(imp, action = 1)

vars_all <- c(
  "num_C_O_D__mg_l",
  "num_Suspended_Solids_mg_l",
  "num_Nitrate_mg_l_as_N",
  "num_Total_B_O_D__mg_l",
  "pH",
  "Conductivity__20_C__µS_cm"
)
keep       <- complete.cases(completed_df[, vars_all])
completed_cc <- completed_df[keep, ]

coords   <- coordinates(sp_gutter)[keep, ]
proj4str <- proj4string(sp_gutter)

sp_gutter_imp_cc <- SpatialPointsDataFrame(
  coords      = coords,
  data        = completed_cc,
  proj4string = CRS(proj4str)
)

library(GWmodel)
fmla <- num_C_O_D__mg_l ~ num_Suspended_Solids_mg_l +
                                      num_Nitrate_mg_l_as_N +

                                      num_Total_B_O_D__mg_l +
                                        Conductivity__20_C__µS_cm +
                                      pH


bw_cv_cc <- bw.gwr(
  formula  = fmla,
  data     = sp_gutter_imp_cc,
  approach = "CV",
  kernel   = "gaussian",
  longlat  = FALSE
)
cat("Fixed bandwidth:", bw_cv_cc, "\n")

gwr_res_cc <- gwr(
  formula   = fmla,
  data      = sp_gutter_imp_cc,
  bandwidth = bw_cv_cc,
  gweight   = gwr.Gauss,
  hatmatrix = TRUE,
  se.fit    = TRUE,
  longlat   = FALSE
)

print(gwr_res_cc)

