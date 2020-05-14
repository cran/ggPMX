## ----load_package, echo=FALSE,warning=FALSE,message=FALSE---------------------
knitr::opts_chunk$set(out.width = "100%", warning = FALSE, message = FALSE)
library(ggPMX)
library(ggplot2)
library(xtable)
library(knitr)

theophylline <- file.path(system.file(package = "ggPMX"), "testdata", "theophylline")
work_dir <- file.path(theophylline, "Monolix")
input_data <- file.path(theophylline, "data_pk.csv")

ctr <- theophylline()

## ----illustrate_diagnostic, out.width='.25\\linewidth', fig.width=4, fig.height=4, fig.show='hold', fig.align='center', echo=FALSE----
ctr <- theophylline()
ctr %>% pmx_plot_dv_pred
ctr %>% pmx_plot_npde_time
ctr %>% pmx_plot_vpc
ctr %>% pmx_plot_eta_box
ctr %>% pmx_plot_eta_matrix(
  shrink=list(size=3,hjust=1.5))

## ----echo=FALSE, out.width='80%', fig.align='center'--------------------------
knitr::include_graphics("./ggPMX_arch.png")

## ----echo=FALSE---------------------------------------------------------------
theophylline_path <- file.path(system.file(package = "ggPMX"), "testdata", "theophylline")
work_dir          <- file.path(theophylline_path, "Monolix")
input_data_path   <- file.path(theophylline_path, "data_pk.csv")

input_data_theo   <- read.csv(input_data_path)
head(input_data_theo)

## ----echo=FALSE,results='asis',out.width='.9\\linewidth'----------------------

out <- rbind(
  c("sys", "Software used for model fittng (Monolix or nlmixr)", "mlx, mlx2018, nm"),
  c("config", "A pre-defined configuration is a set of default settings", "standing"),
  c("directory", "Path to the directory containing model output files", ""),
  c("input", "Path to input modeling dataset (dataset used for model fitting)", ""),
  c("dv", "Measurable variable name, as defined in the input modeling dataset", "DV, LIDV, LNDV, Y, etc."),
  c("dvid", "Endpoint (output) name, as defined in the input modeling dataset", "DVID, YTYPE, CMT, etc.")
)

colnames(out) <- c("Argument", "Description", "Values")

xt <- xtable(head(out), label = "tab:pmx_mandatory", caption = "Mandatory arguments of pmx() function")
print(xt, comment = F)

## ----init_ctr-----------------------------------------------------------------
theophylline_path <- file.path(system.file(package = "ggPMX"), "testdata", "theophylline")
work_dir          <- file.path(theophylline_path, "Monolix")
input_data_path   <- file.path(theophylline_path, "data_pk.csv")

ctr <- pmx(
  sys       = "mlx",
  config    = "standing",
  directory = work_dir,
  input     = input_data_path,
  dv        = "Y",
  dvid      = "DVID"
)

## -----------------------------------------------------------------------------
ctr <- theophylline()

## -----------------------------------------------------------------------------
ctr <- pmx_mlx(
  config    = "standing",
  directory = work_dir,
  input     = input_data_path,
  dv        = "Y",
  dvid      = "DVID"
)

## -----------------------------------------------------------------------------
mlxtran_path <- file.path(system.file(package = "ggPMX"), 
                          "testdata", "1_popPK_model", "project.mlxtran")

ctr <- pmx_mlxtran(file_name = mlxtran_path)

## ---- echo=F------------------------------------------------------------------
pkpd_path       <- file.path(system.file(package = "ggPMX"), "testdata", "pk_pd")
pkpd_work_dir   <- file.path(pkpd_path, "RESULTS")
pkpd_input_file <- file.path(pkpd_path, "pk_pd.csv")

input_data <- read.csv(pkpd_input_file)
head(input_data)

## -----------------------------------------------------------------------------
pkpd_path       <- file.path(system.file(package = "ggPMX"), "testdata", "pk_pd")
pkpd_work_dir   <- file.path(pkpd_path, "RESULTS")
pkpd_input_file <- file.path(pkpd_path, "pk_pd.csv")

ep <- pmx_endpoint(
  code  = "4",
  label = "some_label",
  unit  = "some_unit",
  file.code = "2", # will use predictions2.txt and finegrig2.txt
  trans = "log10"
)

ctr <- pmx_mlx(
  config    = "standing",
  directory = pkpd_work_dir,
  input     = pkpd_input_file,
  dv        = "dv",
  dvid      = "dvid",
  endpoint  = ep
)

## ---- echo=T, eval = FALSE----------------------------------------------------
#  pmx_mlx(
#    dvid = "YTYPE", ## use this column as obseration id
#    endpoint = 1,   ## select the first endpoint
#    ...)            ## other pmx parameters , config, input,etc..

## ----init_ctr_covar-----------------------------------------------------------
theophylline_path <- file.path(system.file(package = "ggPMX"), "testdata", "theophylline")
work_dir          <- file.path(theophylline_path, "Monolix")
input_data_path   <- file.path(theophylline_path, "data_pk.csv")

ctr <- pmx_mlx(
  config    = "standing",
  directory = work_dir,
  input     = input_data_path,
  dv        = "Y",
  dvid      = "DVID",
  cats      = c("SEX"),
  conts     = c("WT0", "AGE0"),
  strats    = c("STUD", "SEX")
)

## ----get_covar----------------------------------------------------------------
ctr %>% get_cats()
ctr %>% get_conts()
ctr %>% get_strats()
ctr %>% get_covariates()

## ----display_ctr--------------------------------------------------------------
ctr

## ----plot_lists---------------------------------------------------------------
ctr %>% plot_names()

## ----plot_types---------------------------------------------------------------
ctr %>% plots()

## ----datasets_list,echo=FALSE,results='asis'----------------------------------

out <- rbind(
  c("input", "Input modeling dataset"),
  c("estimates", "Estimated population parameters"),
  c("eta", "Random effects, their standard deviation and residual errors (to calculate shrinkage)"),
  c("predictions", "Observations and predictions at times of observations dataset"),
  c("finegrid", "Additional predictions (at times without observations)")
)

colnames(out) <- c("ggPMX dataset", "Description")
# knitr::kable(out)
# latex(head(out), file='', label='tab:ggPMX_datasets', caption='ggPMX datasets',where = "!htbp")
xt <- xtable(head(out), label = "tab:ggPMX_datasets", caption = "ggPMX datasets")
print(xt, comment = F)

## ----echo=FALSE---------------------------------------------------------------
theophylline <- file.path(system.file(package = "ggPMX"), "testdata", "theophylline")
work_dir <- file.path(theophylline, "Monolix")
input_data <- file.path(theophylline, "data_pk.csv")

ctr <- theophylline()

## -----------------------------------------------------------------------------
ctr %>% plot_names()

## ----fig.height=5, fig.width=8------------------------------------------------
ctr %>% pmx_plot_iwres_time

## ----basics_res, eval=F, out.width='.48\\linewidth', fig.height=4, fig.width=6, fig.show='hold', fig.align='center'----
#  ctr %>% pmx_plot_dv_pred
#  ctr %>% pmx_plot_dv_ipred
#  
#  ctr %>% pmx_plot_iwres_time
#  ctr %>% pmx_plot_npde_time
#  
#  ctr %>% pmx_plot_iwres_ipred
#  ctr %>% pmx_plot_abs_iwres_ipred
#  
#  ctr %>% pmx_plot_npde_pred

## ----basics_ebe_hist, eval=F , fig.height=3, fig.width=3, fig.show='hold', fig.align='center'----
#  ctr %>% pmx_plot_eta_hist
#  ctr %>% pmx_plot_eta_box

## ----basics_indiv, eval=F, fig.height=6, fig.width=6, fig.show='hold', fig.align='center'----
#  ctr %>% pmx_plot_individual(npage = 1)

## ----basics_qq, eval=F, fig.height=3, fig.width=3, fig.show='hold', fig.align='center'----
#  ctr %>% pmx_plot_npde_qq
#  ctr %>% pmx_plot_iwres_qq

## ----basics_matrix_plot, eval=F,  fig.height=6, fig.width=6, fig.show='hold', fig.align='center'----
#  ctr %>% pmx_plot_eta_matrix

## -----------------------------------------------------------------------------

theoph_path <- file.path(
  system.file(package = "ggPMX"), "testdata",
  "theophylline"
)
WORK_DIR <- file.path(theoph_path, "Monolix")
input_file <- file.path(theoph_path, "data_pk.csv")
vpc_file <- file.path(theoph_path, "sim.csv")

ctr <- pmx_mlx(
  config = "standing",
  directory = WORK_DIR,
  input = input_file,
  dv = "Y",
  cats = c("SEX"),
  conts = c("WT0", "AGE0"),
  strats = "STUD",
  settings = pmx_settings(
    use.labels=TRUE,
    cats.labels=list(
      SEX=c("0"="Male","1"="Female")
    )
  ),
  sim = pmx_sim(
    file = vpc_file,
    irun ="rep",
    idv="TIME"
  )
)



## ----fig.height=4, fig.width=6------------------------------------------------
ctr %>% pmx_plot_vpc

## ----fig.height=3.5, fig.width=5.5--------------------------------------------
ctr %>% pmx_plot_vpc(type ="scatter")

## ----fig.height=3.5, fig.width=5.5--------------------------------------------
ctr %>% pmx_plot_vpc(bin=pmx_vpc_bin(style = "kmeans",n=5))

## ----fig.height=7, fig.width=6------------------------------------------------
ctr %>% pmx_plot_vpc(strat.facet="SEX",facets=list(nrow=2))

## ----fig.height=7, fig.width=6------------------------------------------------
ctr %>% pmx_plot_vpc(
  strat.facet="SEX",
  facets=list(nrow=2),
  type="percentile",
  is.draft = FALSE,
  pi = pmx_vpc_pi(interval = c(0.1,0.9),
              median=list(color="green"),
              extreme= list(color="green")),
  obs = pmx_vpc_obs(color="blue",shape=18,size=2),
  ci = pmx_vpc_ci(interval = c(0.1,0.9),
              median=list(fill="red"))
)

## ----echo=FALSE---------------------------------------------------------------
theophylline <- file.path(system.file(package = "ggPMX"), "testdata", "theophylline")
work_dir <- file.path(theophylline, "Monolix")
input_data <- file.path(theophylline, "data_pk.csv")

ctr <- theophylline()

## ----eval=FALSE---------------------------------------------------------------
#  ctr %>% pmx_report(name='Diagnostic_plots2',
#                     save_dir = work_dir,
#                     format='both')

## ----eval=FALSE---------------------------------------------------------------
#  ctr %>% pmx_report(name='Diagnostic_plots1',
#                     save_dir = work_dir,
#                     format='report')

## ----eval=FALSE---------------------------------------------------------------
#  ctr %>% pmx_report(name='Diagnostic_plots3',
#                     save_dir = work_dir,
#                     format='report',
#                     template=file.path(work_dir,'Diagnostic_plots1.Rmd'))

## ----eval=F-------------------------------------------------------------------
#  ctr %>% pmx_plot_xx(list of options)

## ----eval=F-------------------------------------------------------------------
#  ctr %>% pmx_update(???xx???, list of options)

## ----eval=F-------------------------------------------------------------------
#  ctr %>% get_plot_config("xx")

## -----------------------------------------------------------------------------
ctr = theophylline()
bloq = list(show = FALSE, color = "blue", alpha = 0.2, size = 3, pch = 8)

## ----settings_example, fig.width=5, fig.height=4,eval=FALSE-------------------
#  
#  ## set one or more settings
#  my_settings <- pmx_settings(
#    is.draft   = FALSE,
#    use.abbrev = TRUE,
#    ...) ### set other settings parameters here
#  ctr <-
#    pmx_mlx(
#      ..., ## put here other pmx parametes
#      settings = my_settings
#    )

## ----settings_is_draft,fig.height=3, fig.width=3, fig.show='hold', fig.align='center'----

ctr <- theophylline(settings = pmx_settings(is.draft = FALSE))


## ----settings_get_abbrev------------------------------------------------------
ctr %>% get_abbrev

## ----settings_set_abbrev------------------------------------------------------
ctr %>% set_abbrev(TIME="TIME after the first dose")

## ----settings_use.abbrev------------------------------------------------------
ctr <- theophylline(settings=pmx_settings(use.abbrev = TRUE))
ctr %>% set_abbrev(TIME="Custom TIME axis")
ctr %>% pmx_plot_npde_time


## ----settings_use.finegrid,fig.height=9, fig.width=7--------------------------
ctr <- theophylline()
ctr %>% pmx_plot_individual(use.finegrid =FALSE)

## ----settings_color_scales_local----------------------------------------------
ctr <- theophylline()
ctr %>% pmx_plot_npde_time(strat.color="STUD")+ 
      ggplot2::scale_color_manual(
        "Study",
        labels=c("Study 1","Study 2"),
        values=c("1"="green","2"="blue"))
    


## ----settings_solor_scales,fig.height=5, fig.width=8--------------------------

ctr <- theophylline(
  settings=
    pmx_settings(
      color.scales=list(
        "Study",
        labels=c("Study 1","Study 2"),
        values=c("1"="orange","2"="magenta"))
    )
)

ctr %>% pmx_plot_npde_time(strat.color="STUD")

## ----settings_solor_scales_a,fig.height=5, fig.width=11-----------------------
ctr %>% pmx_plot_eta_box(strat.color="STUD")

## ----settings_cat_labels------------------------------------------------------


ctr <- theophylline(
  settings=
    pmx_settings(
      cats.labels=list(
        SEX=c("0"="M","1"="F"),
        STUD=c("1"="Study 1","2"="Study 2")
      ),
      use.labels = TRUE
    )
)

ctr %>% pmx_plot_npde_time(strat.facet=~SEX)

## ----settings_cat_labels2, fig.height=5, fig.width=8--------------------------


ctr <- theophylline(
  settings=
    pmx_settings(
      cats.labels=list(
        SEX=c("0"="M","1"="F"),
        STUD=c("1"="Study 1","2"="Study 2")
      ),
      use.labels = TRUE
    )
)

ctr %>% pmx_plot_npde_time(strat.facet=~SEX)

## ----settings_cat_labels3, fig.height=8, fig.width=8--------------------------
ctr %>% pmx_plot_eta_box(strat.facet =~SEX)


## ----plots_list,echo=FALSE,results='asis'-------------------------------------

out <- rbind(
  c("Scatter plot of NPDE vs population predictions", "SCATTER", "npde_pred"),
  c("Scatter plot of NPDE vs time", "SCATTER", "npde_time"),
  c("Scatter plot of IWRES vs time", "SCATTER", "iwres_time"),
  c("Scatter plot of observations vs population predictions", "SCATTER", "dv_pred"),
  c("Scatter plot of observations vs individual predictions", "SCATTER", "dv_ipred"),
  c("Scatter plot of absolute value of IWRES vs individual predictions", "SCATTER", "abs_iwres_ipred"),
  c("Scatter plot of IWRES vs individual predictions", "SCATTER", "iwres_ipred"),
  c("Plots of observations and model predictions per individual", "IND", "individual"),
  c("Histogram of EBE", "DIS", "ebe_hist"),
  c("Boxplot of EBE", "DIS", "ebe_box"),
  c("Distribution and quantile-quantile plot of IWRES", "QQ", "qq_iwres"),
  c("Distribution and correlation structure of RE (`ETA`)", "ETA_PAIRS", "eta_matrix"),
  c("Relationships between RE and categorical covariates", "ETA_COV", "eta_cats"),
  c("Relationships between RE and continuous covariates", "ETA_COV", "eta_conts"),
  c("Visual predictive check (VPC)", "VPC", "vpc")
)

colnames(out) <- c("Plot Name", "ggPMX type", "ggPMX name")
xt <- xtable(out, label = "tab:plots_list", caption = "List of all diagnostic plots")
print(xt, comment = F)

## ----functions_list,echo=FALSE,results='asis'---------------------------------

out <- rbind(
  c("1", "pmx, or pmx_mlx", "Creates a controller"),
  c("2", "plot_names or plots", "Lists controller plots"),
  c("3", "get_data", "Lists controller data"),
  c("4", "get_plot", "Prints a plot"),
  c("5", "set_plot", "Creates a new plot"),
  c("6", "pmx_update", "Updates an existing plot"),
  c("7", "pmx_filter", "Filters globally the data of the current session"),
  c("8", "pmx_copy", "Returns a deep copy of the controller")
)

colnames(out) <- c(" ", "Function name", "Description")

xt <- xtable(out, label = "tab:func_list", caption = "List of all `ggPMX` functions")
print(xt, comment = F)

## ----pmx_gpar_args------------------------------------------------------------
args(pmx_gpar)

## ----shrink_comp--------------------------------------------------------------
ctr %>% pmx_comp_shrink

## ----shrink_plot_box----------------------------------------------------------
ctr %>% pmx_plot_eta_box


## ----shrink_plot_hist,fig.height=8, fig.width=5-------------------------------
ctr %>% pmx_plot_eta_hist


## ----shrink_plot_no-----------------------------------------------------------
ctr %>%   pmx_plot_eta_box( is.shrink = FALSE) 


## ---- compute_var-------------------------------------------------------------
ctr %>% pmx_comp_shrink(  fun = "var")

## ---- shrink_plot_var---------------------------------------------------------
ctr %>% pmx_plot_eta_box( shrink=list(fun = "var"))

## ----shrink_comp_strat--------------------------------------------------------
ctr %>% pmx_comp_shrink(strat.facet = ~SEX)

## ----shrink_comp_strat_color--------------------------------------------------
ctr %>% pmx_comp_shrink(strat.color = "SEX")

## ----shrink_plot_strat, fig.width=9, fig.height=8-----------------------------
ctr %>% pmx_plot_eta_hist(is.shrink = TRUE, strat.facet = ~SEX,
                          facets=list(scales="free_y"))

## ----fig.width=12, fig.height=6-----------------------------------------------
ctr %>% pmx_plot_eta_box(is.shrink = TRUE, strat.facet = "SEX",
                          facets=list(scales="free_y",ncol=2))

