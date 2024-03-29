% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmx-plots-scatter.R
\name{residual_scatter}
\alias{residual_scatter}
\alias{pmx_plot_dv_pred}
\alias{pmx_plot_dv_ipred}
\alias{pmx_plot_iwres_ipred}
\alias{pmx_plot_abs_iwres_ipred}
\alias{pmx_plot_abs_iwres_time}
\alias{pmx_plot_iwres_time}
\alias{pmx_plot_npde_time}
\alias{pmx_plot_npde_pred}
\alias{pmx_plot_npd_time}
\alias{pmx_plot_npd_epred}
\alias{pmx_plot_npd_pred}
\alias{pmx_plot_cwres_time}
\alias{pmx_plot_cwres_cpred}
\alias{pmx_plot_cwres_pred}
\title{Scatter residual plots}
\usage{
residual_scatter(
  sim_blq,
  point,
  is.hline,
  hline,
  dname,
  bloq,
  filter,
  strat.facet,
  facets,
  strat.color,
  trans,
  pmxgpar,
  labels,
  axis.title,
  axis.text,
  ranges,
  is.smooth,
  smooth,
  is.band,
  band,
  is.draft,
  draft,
  is.identity_line,
  identity_line,
  scale_x_log10,
  scale_y_log10,
  color.scales,
  ...
)

pmx_plot_dv_pred(ctr, ...)

pmx_plot_dv_ipred(ctr, ...)

pmx_plot_iwres_ipred(ctr, ...)

pmx_plot_abs_iwres_ipred(ctr, ...)

pmx_plot_abs_iwres_time(ctr, ...)

pmx_plot_iwres_time(ctr, ...)

pmx_plot_npde_time(ctr, ...)

pmx_plot_npde_pred(ctr, ...)

pmx_plot_npd_time(ctr, ...)

pmx_plot_npd_epred(ctr, ...)

pmx_plot_npd_pred(ctr, ...)

pmx_plot_cwres_time(ctr, ...)

pmx_plot_cwres_cpred(ctr, ...)

pmx_plot_cwres_pred(ctr, ...)
}
\arguments{
\item{sim_blq}{\code{logical} if TRUE uses sim_blq values for plotting. Only for Monolix 2018 and later.}

\item{point}{\code{list} geom_point graphical parameters.}

\item{is.hline}{\code{logical} if TRUE add horizontal line y=0 ( TRUE by default).}

\item{hline}{\code{list} geom_hline graphical parameters.}

\item{dname}{\code{character} name of dataset to be used. User can create his own
dataset using \code{\link{set_data}} and pass it as dname to be plotted.}

\item{bloq}{\code{pmxBLOQ} object created by \code{\link{pmx_bloq}}.

\strong{pmx_update parameters}}

\item{filter}{\code{expression} filter which will be applied to plotting data.}

\item{strat.facet}{\code{formula} optional stratification parameter by facetting.
This split plot by strats(each strat in a facet)}

\item{facets}{\code{list} facet_wrap parameters.}

\item{strat.color}{\code{character} optional stratification parameter by grouping.
This will split the plot by group (color) of strat.}

\item{trans}{\code{character} define the transformation to apply on x or y or both variables}

\item{pmxgpar}{a object of class pmx_gpar possibly the output of the

\strong{pmx_gpar: Shared basic graphics parameters}}

\item{labels}{\code{list} list containing plot and/or axis labels: title, subtitle, x , y}

\item{axis.title}{\code{list} containing element_text attributes to customize
the axis title. (similar to ggplot2 axis.title theme)}

\item{axis.text}{\code{list} containing element_text attributes to customize
the axis text (similar to ggplot2 axis.text theme)}

\item{ranges}{\code{list} limits of x/y ranges}

\item{is.smooth}{\code{logical} if set to TRUE add smooth layer}

\item{smooth}{\code{list} geom_smooth graphical/smoothing fun parameters}

\item{is.band}{\code{logical} if TRUE add horizontal band}

\item{band}{\code{list} horizontal band parameters. geom_hline graphical parameters.}

\item{is.draft}{\code{logical} if TRUE add draft layer}

\item{draft}{\code{list} draft layer parameters. geom_text graphical parameters.}

\item{is.identity_line}{\code{logical} if TRUE add an identity line}

\item{identity_line}{\code{list}geom_abline graphical parameters.}

\item{scale_x_log10}{\code{logical} if TRUE use log10 scale for x axis.}

\item{scale_y_log10}{\code{logical} if TRUE use log10 scale for y axis.}

\item{color.scales}{\code{list} define scales parameter in case of strat.color \code{\link{pmx_settings}}}

\item{...}{others graphics parameters passed :
\itemize{
\item \code{\link{pmx_gpar}} internal function to customize shared graphical parameters
\item \code{\link{residual}} generic object for all residual (scatter) plots .
\item \code{\link{pmx_update}} function.
\item aess can be used to change time variable within the plot (e.g. aess = list(x="TADQBW"))
}

\strong{residual parameters}}

\item{ctr}{pmx controller}
}
\value{
ggplot2 object
}
\description{
Scatter residual plots

DV vs PRED plot

DV vs IPRED plot

IWRES vs IPRED plot

|IWRES| vs IPRED plot

|IWRES| vs TIME plot

IWRES vs TIME plot

NPDE vs TIME plot

NPDE vs PRED plot

NPD vs TIME plot

NPD vs EPRED plot

NPD vs PRED plot

CWRES vs TIME plot

CWRES vs CPRED plot

CWRES vs PRED plot
}
\examples{
\donttest{
# NOTES ######################
# examples are availables for all residual plots:
# - pmx_plot_abs_iwres_ipred
# - pmx_plot_dv_ipred
# - pmx_plot_dv_pred
# - pmx_plot_iwres_ipred
# - pmx_plot_iwres_time
# - pmx_plot_npde_time


# basic use  ---------------------------------------

ctr <- theophylline()
p <- ctr \%>\% pmx_plot_dv_pred()
## p is a ggplot2 object you can add any layer here
p + ggplot2::theme_minimal()

# update graphical parameter  ----------------------

## update labels
ctr \%>\% pmx_plot_dv_pred(
  labels = list(title = "DV versus PRED new title")
)

## remove draft
ctr \%>\% pmx_plot_dv_pred(is.draft = FALSE)

## remove horizontal line
ctr \%>\% pmx_plot_dv_pred(is.hline = FALSE)

## custom point aes and static parameters
## we can customize any geom_point parameter
ctr \%>\% pmx_plot_dv_pred(
  point = list(aes(alpha = DV), color = "green", shape = 4)
)



# stratification  ----------------------------------

## continuous stratification
ctr \%>\% pmx_plot_dv_pred(strat.color = ~WT0)
## categorical stratification
ctr \%>\% pmx_plot_dv_pred(strat.facet = ~SEX)
## using formula notation
ctr \%>\% pmx_plot_dv_pred(strat.facet = STUD~SEX)

# subsetting  --------------------------------------

## we can use any expression involving the data
ctr \%>\% pmx_plot_dv_pred(filter = DV > mean(DV) & PRED < median(PRED))
## filter and stratify
ctr \%>\% pmx_plot_dv_pred(filter = SEX == 1, strat.facet = ~SEX)


# transformation  --------------------------------------

## apply a log transformation in y
ctr \%>\% pmx_plot_dv_pred(trans = "log10_y")
}
}
\concept{residual}
