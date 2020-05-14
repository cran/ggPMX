## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ggPMX)
library(nlmixr)

## -----------------------------------------------------------------------------
one.compartment <- function() {
    ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
    })
    model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
    })
}
nlmixr(one.compartment)

## ---- eval=FALSE--------------------------------------------------------------
#  fit <- nlmixr(one.compartment, theo_sd, est="saem", control=list(print=0))
#  saveRDS(fit,"fit.rds")

## -----------------------------------------------------------------------------
fit <- readRDS("fit.rds")
print(fit)

## -----------------------------------------------------------------------------
fit %>%
    pmx_nlmixr(vpc = FALSE) -> ## VPC is turned on by default, can turn off.
    ctr ## Assigned to controller 

## -----------------------------------------------------------------------------
ctr %>% pmx_plot_dv_pred
ctr %>% pmx_plot_npde_time
ctr %>% pmx_plot_vpc
ctr %>% pmx_plot_eta_box
ctr %>% pmx_plot_eta_matrix(
  shrink=list(size=3,hjust=1.5))

