###############################################################################
## Replication code for Cattaneo and Palomba (2025)
## "Leveraging Covariates in Regression Discontinuity Designs"
## Last modified: 2025-01-15
###############################################################################

rm(list=ls(all=TRUE))

##########################################
# Load stuff 
require(pacman)
pacman::p_load(ggplot2, wesanderson, dplyr, latex2exp,
               rdrobust, sandwich, lmtest)

theme_set(theme_bw())

##########################################
# Set paths
path <- "YOUR_PATH"
path.data <- paste0(path, "data/")
path.code <- paste0(path, "code/")
path.fig  <- paste0(path, "figures/")
path.out  <- paste0(path, "tables/")

# Auxiliary functions
source(paste0(path.code, "funs.R"))

################################################################################
## Prepare Head Start Data
################################################################################
data <- read.csv(paste0(path.data, "headstart.csv"))

data$census1960_pop_ind <- 1*(data$census1960_pop >= 10000L)

cutoff <- 59.1968  # head start cutoff
w.vals <- 2L       # number of heterogeneity categories

# create variables for analysis
y <- data$mort_age59_related_postHS  # outcome variable
x <- data$povrate60                  # running variable
w <- data$census1960_pop_ind         # covariate for heterogeneity
z <- cbind(data$census1960_pctblack, # covariates for efficiency
           data$census1960_pctsch1417, data$census1960_pctsch534,
           data$census1960_pctsch25plus, data$census1960_pop1417, 
           data$census1960_pop534, data$census1960_pop25plus,
           data$census1960_pcturban)


################################################################################
## 1) Create rdplots for all units and two sub-groups
################################################################################
# list to store all bandwidths
h.opt <- list()

# canonical RD
rd <- rdrobust(y, x, c=cutoff, rho=1)
h.opt[["all"]] <- rd$bws[1,1]

# heterogeneous RD
for (w.val in c(0, 1)) {
  rd <- rdrobust(y, x, c=cutoff, rho=1, subset = w == w.val)
  h.opt[[paste0("w", w.val)]] <- rd$bws[1,1]
}

#########################
# global plot (Figure 1a)
#########################

rdplot.res <- rdplot(y, x, c=cutoff, p=1, kernel="triangular", hide = TRUE)

toplot.bins <- data.frame(
  bin_mean_x = rdplot.res$vars_bins$rdplot_mean_bin,
  bin_mean_y = rdplot.res$vars_bins$rdplot_mean_y,
  treated = factor(ifelse(rdplot.res$vars_bins$rdplot_mean_bin >= cutoff, 1, 0)),
  group = "all units",
  type = "global",
  within.h = factor(ifelse(abs(rdplot.res$vars_bins$rdplot_mean_bin-cutoff)<=h.opt[["all"]],1,0))
)

toplot.poly <- data.frame(
  poly_y = rdplot.res$vars_poly$rdplot_y,
  poly_x = rdplot.res$vars_poly$rdplot_x,
  treated = factor(rep(c(0, 1), each = 500)),
  group = "all units",
  type = "global"
)

###############################################
# local plot with optimal bandwidth (Figure 1b)
###############################################

rdplot.res <- rdplot(y, x, c=cutoff, p=1, kernel="triangular",
                     subset= -h.opt[["all"]] + cutoff <= x & x <= h.opt[["all"]] + cutoff,
                     h=c(h.opt[["all"]], h.opt[["all"]]), hide=TRUE)

aux.bins <- data.frame(
  bin_mean_x = rdplot.res$vars_bins$rdplot_mean_bin,
  bin_mean_y = rdplot.res$vars_bins$rdplot_mean_y,
  treated = factor(ifelse(rdplot.res$vars_bins$rdplot_mean_bin >= cutoff, 1, 0)),
  group = "all units",
  type = "local",
  within.h = factor(1)
)

aux.poly <- data.frame(
  poly_y = rdplot.res$vars_poly$rdplot_y,
  poly_x = rdplot.res$vars_poly$rdplot_x,
  treated = factor(rep(c(0, 1), each = 500)),
  group = "all units",
  type = "local"
)

toplot.bins <- rbind(toplot.bins, aux.bins)
toplot.poly <- rbind(toplot.poly, aux.poly)

# sample sizes
n.all.l <- sum(x < cutoff, na.rm=TRUE)
n.0.l <- sum(x < cutoff & w == 0, na.rm=TRUE)
n.1.l <- sum(x < cutoff & w == 1, na.rm=TRUE)
n.all.r <- sum(x >= cutoff, na.rm=TRUE)
n.0.r <- sum(x >= cutoff & w == 0, na.rm=TRUE)
n.1.r <- sum(x >= cutoff & w == 1, na.rm=TRUE)

N.all.l <- sum(x >= cutoff - h.opt[["all"]] & x < cutoff, na.rm=TRUE)
N.0.l <- sum(x >= cutoff - h.opt[["w0"]] & x < cutoff & w == 0, na.rm=TRUE)
N.1.l <- sum(x >= cutoff - h.opt[["w1"]] & x < cutoff & w == 1, na.rm=TRUE)
N.all.r <- sum(x >= cutoff & x <= cutoff + h.opt[["all"]], na.rm=TRUE)
N.0.r <- sum(x >= cutoff & x <= cutoff + h.opt[["w0"]] & w == 0, na.rm=TRUE)
N.1.r <- sum(x >= cutoff & x <= cutoff + h.opt[["w1"]] & w == 1, na.rm=TRUE)

labs <- c("population <10k", "population $\\geq$ 10k")

for (w.val in c(0,1)) {

  #########################
  # global plot (Figure 2a)
  #########################
  
  rdplot.res <- rdplot(y, x, c=cutoff, p=1, kernel="triangular", hide=TRUE,
                       subset = w == w.val)
  
  aux.bins <- data.frame(
    bin_mean_x = rdplot.res$vars_bins$rdplot_mean_bin,
    bin_mean_y = rdplot.res$vars_bins$rdplot_mean_y,
    treated = factor(ifelse(rdplot.res$vars_bins$rdplot_mean_bin >= cutoff, 1, 0)),
    group = labs[w.val+1],
    type = "global",
    within.h = factor(ifelse(abs(rdplot.res$vars_bins$rdplot_mean_bin-cutoff)<=h.opt[[paste0("w",w.val)]],1,0))
  )
  
  aux.poly <- data.frame(
    poly_y = rdplot.res$vars_poly$rdplot_y,
    poly_x = rdplot.res$vars_poly$rdplot_x,
    treated = factor(rep(c(0, 1), each = 500)),
    group = labs[w.val+1],
    type = "global"
  )

  ###############################################
  # local plot with optimal bandwidth (Figure 2b)
  ###############################################
  
  rdplot.res <- rdplot(y, x, c=cutoff, p=1, kernel="triangular",
                       subset= -h.opt[[paste0("w", w.val)]] + cutoff<= x & 
                         x <= h.opt[[paste0("w", w.val)]] + cutoff & w == w.val,
                       h=c(h.opt[[paste0("w", w.val)]], h.opt[[paste0("w", w.val)]]), hide=TRUE)

  aux.bins.zoom <- data.frame(
    bin_mean_x = rdplot.res$vars_bins$rdplot_mean_bin,
    bin_mean_y = rdplot.res$vars_bins$rdplot_mean_y,
    treated = factor(ifelse(rdplot.res$vars_bins$rdplot_mean_bin >= cutoff, 1, 0)),
    group = labs[w.val+1],
    type = "local",
    within.h = factor(1)
  )
  
  aux.poly.zoom <- data.frame(
    poly_y = rdplot.res$vars_poly$rdplot_y,
    poly_x = rdplot.res$vars_poly$rdplot_x,
    treated = factor(rep(c(0, 1), each = 500)),
    group = labs[w.val+1],
    type = "local"
  )

  toplot.bins <- rbind(toplot.bins, aux.bins, aux.bins.zoom)
  toplot.poly <- rbind(toplot.poly, aux.poly, aux.poly.zoom)
}

toplot.bins$treated <- factor(toplot.bins$treated, levels = c(0,1),
                              labels=c("control", "treated"))
toplot.poly$treated <- factor(toplot.poly$treated, levels = c(0,1),
                              labels=c("control", "treated"))

# create and store Figure 1 and Figure 2
plotsGet(toplot.bins, toplot.poly, path.fig, h.opt)

################################################################################
## 2) Produce estimates for all units and sub-groups
################################################################################

#######################################
## Check "population" is pre-treatment
rdplot(w, x, c=cutoff)
summary(rdrobust(w, x, c=cutoff))


#######################################
## Prepare table
out <- matrix(NA,7,3)
rownames(out) <- c("$\\hat{\\tau}$",             # point estimate (bias-corrected)
                   "95\\% RCI",                  # robust bias-corrected confidence intervals
                   "CI length change (\\%)",     # CI length change wrt no covs case
                   "$p$-value",                  # p-value of t-test for null treatment effect
                   "$h$",                        # estimation MSE-optimal bandwidth
                   "$N_-\\,|\\, N_+$",           # number of obs to the left and right of cutoff
                   "\\% treatment effect")       # treatment effect compared with mean control outcome

colnames(out) <- c("canonical",
                   "with covariates for efficiency",
                   "with all covariates")

################################################################################
## 2.1) Analysis on all units

h.list <- list()

#################################
# 2.1a) RD on all units w/o covs
#################################

rd <- rdrobust(y, x, c=cutoff, rho=1)
h.list[["all"]] <- rd$bws[1,1]
mean.y0 <- meanControlsGet(y=y, x=x, c=cutoff, h=h.list[["all"]])

out[1,1] <- paste0("$", round(rd$coef[1],2), "$")
out[2,1] <- paste0("$[", round(rd$ci[3,1],2),",", round(rd$ci[3,2],2), "]$")
il <- rd$ci[3,2] - rd$ci[3,1]
out[3,1] <- "-"
out[4,1] <- round(rd$pv[3],3)
out[5,1] <- round(h.list[["all"]],3)
out[6,1] <- paste0("$", rd$N_h[1]," | ", rd$N_h[2], "$")
out[7,1] <- paste0("$", round(100*rd$coef[1]/mean.y0, 2), "$") 
  
#############################################
# 2.1b) RD on all units w covs for efficiency
#############################################

rd <- rdrobust(y, x, c=cutoff, covs=z, rho=1)
h.list[["covs_eff"]] <- rd$bws[1,1]
mean.y0 <- meanControlsGet(y=y, x=x, c=cutoff, z=z, h=h.list[["covs_eff"]])

out[1,2] <- paste0("$", round(rd$coef[1],2), "$")
out[2,2] <- paste0("$[", round(rd$ci[3,1],2),",", round(rd$ci[3,2],2), "]$")
il2 <- rd$ci[3,2] - rd$ci[3, 1]
out[3,2] <- paste0("$", round(((il2/il - 1) * 100),2))
out[4,2] <- round(rd$pv[3], 3)
out[5,2] <- round(h.list[["covs_eff"]], 3)
out[6,2] <- paste0("$", rd$N_h[1]," | ", rd$N_h[2], "$")
out[7,2] <- paste0("$", round(100*rd$coef[1]/mean.y0, 2), "$") 


#################################################
# 2.1c) RD on all units w all covs for efficiency
#################################################

rd <- rdrobust(y, x, c=cutoff, covs=cbind(z, w), rho=1)
mean.y0 <- meanControlsGet(y=y, x=x, c=cutoff, z=cbind(z, w), h=rd$bws[1,1])

out[1,3] <- paste0("$", round(rd$coef[1],2), "$")
out[2,3] <- paste0("$[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],2), "]$")
il3 <- rd$ci[3,2] - rd$ci[3,1]
out[3,3] <- paste0("$", round(((il3/il - 1) * 100),2))
out[4,3] <- round(rd$pv[3],3)
out[5,3] <- round(rd$bws[1],3)
out[6,3] <- paste0("$", rd$N_h[1]," | ", rd$N_h[2], "$")
out[7,3] <- paste0("$", round(100*rd$coef[1]/mean.y0, 2), "$") 

write.csv(out, paste0(path.out, "estimatesAll.csv"))

################################################################################
## 2.2) Heterogeneity analysis 

out.hte <- matrix(NA, 7, 5)

#########################################
# 2.2a-b) separate RD for each value of w
#########################################

for (w.val in c(0, 1)) {
  # without covariates for efficiency
  rd <- rdrobust(y, x, c=cutoff, rho=1, subset = w == w.val)
  mean.y0 <- meanControlsGet(y=y, x=x, c=cutoff, h=rd$bws[1,1], subset = w == w.val)
  
  out.hte[1,w.val+1] <- paste0("$", round(rd$coef[1],2), "$")
  out.hte[2,w.val+1] <- paste0("$[", round(rd$ci[3,1],2),",", round(rd$ci[3,2],2), "]$")
  il <- rd$ci[3,2] - rd$ci[3,1]
  out.hte[3,w.val+1] <- "-"
  out.hte[4,w.val+1] <- round(rd$pv[3],3)
  out.hte[5,w.val+1] <- round(rd$bws[1], 3)
  out.hte[6,w.val+1] <- paste0("$",rd$N_h[1],"\\, | \\,", rd$N_h[2],"$")
  out.hte[7,w.val+1] <- paste0("$", round(100*rd$coef[1]/mean.y0, 2), "$") 
  
  # with covariates for efficiency
  rd <- rdrobust(y, x, c=cutoff, covs=z, rho=1, subset = w == w.val)
  mean.y0 <- meanControlsGet(y=y, x=x, c=cutoff, z=z, h=rd$bws[1,1], subset = w == w.val)

  out.hte[1,w.val+4] <- paste0("$", round(rd$coef[1],2), "$")
  out.hte[2,w.val+4] <- paste0("$[", round(rd$ci[3,1],2),",", round(rd$ci[3,2],2), "]$")
  il2 <- rd$ci[3,2] - rd$ci[3,1]
  out.hte[3,w.val+4] <- paste0("$", round(((il2/il - 1) * 100), 2), "$")
  out.hte[4,w.val+4] <- round(rd$pv[3],3)
  out.hte[5,w.val+4] <- round(rd$bws[1], 3)
  out.hte[6,w.val+4] <- paste0("$", rd$N_h[1],"\\, | \\,", rd$N_h[2], "$")  
  out.hte[7,w.val+4] <- paste0("$", round(100*rd$coef[1]/mean.y0, 2), "$") 
}

aux.labs <- rownames(out)
aux.labs[1] <- "$\\widehat{\\tau}(z)$"
rownames(out.hte) <- c(aux.labs) 
aux.labs <- c("population $<$10k", "population $>$10k")

colnames(out.hte) <- c(aux.labs, " ", aux.labs)
out.hte[is.na(out.hte)] <- ""

write.csv(out.hte, paste0(path.out, "estimatesHte.csv"))

