library(ggplot2)

select_covariate_effects <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  VERSION = "v2" # Last version with no VERSION variable definition (overwritten by load)
  load(paste0("results/", filename))
  filename = gsub("(.*)stanfit.Rdata$", "\\1", filename)
  out <- rstan::extract(fit)
  if (VERSION == "v3" || VERSION == "v2"){
    prepped_data = prepare_intervention_effects(out)
  } else {
    prepped_data = prepare_mobility_effects(out, formula_pooling, formula_partialpooling, region_to_country_map)
  }
  plot_covariate_effects(prepped_data[["data"]], prepped_data[["plot_labels"]], filename)
}

prepare_intervention_effects <- function (out){

  alpha <- data.frame(as.matrix(out$alpha))
  plot_labels <- c("School Closure",
                   "Self Isolation",
                   "Public Events",
                   "First Intervention",
                   "Lockdown", 'Social distancing \n encouraged')

  colnames(alpha) <- plot_labels
  first.intervention <- alpha[,c(1,2,3,5,6)] + alpha[,4]
  data1 <- bayesplot::mcmc_intervals_data(first.intervention,prob=.95,transformation=function(x) 1-exp(-x),point_est="mean")
  data1$type <- "First Intervention"

  data2 <- bayesplot::mcmc_intervals_data(alpha[,c(1,2,3,5,6)],  prob=.95,transformation=function(x) 1-exp(-x),point_est="mean")
  data2$type <- "Later Intervention"
  data <- rbind(data1,data2[1:5,])
  return(list(data=data, plot_labels=plot_labels))
}

prepare_mobility_effects <- function (out, formula_pooling, formula_partialpooling, region_to_country_map) {
  
  # stop("covariate effect plot not supported for mobility")

  plot_labels <- regmatches(formula_pooling, gregexpr("[a-zA-Z]+", formula_pooling))[[1]]
  name_labels_partial <- regmatches(formula_partialpooling, gregexpr("[a-zA-Z]+", formula_partialpooling))[[1]]
  
  alpha <- data.frame(as.matrix(out$alpha))
  plot_labels <- paste("global: ",plot_labels,sep="")
  colnames(alpha) <- plot_labels
  data_global <- bayesplot::mcmc_intervals_data(alpha,prob=.95,transformation=function(x) 1-exp(-x),point_est="mean")
  data_global$type <- "pooled"

  plot_labels_partial <- c()
  for(i in 1:length(region_to_country_map)){
    zone <- names(region_to_country_map)[i]
    zone_labels <- paste(paste0(zone, ": "), name_labels_partial, sep="")
    plot_labels_partial <- c(plot_labels_partial, zone_labels)
    alpha_partial <- data.frame(as.matrix(drop(out$alpha_state[,i,])))
    colnames(alpha_partial) <- zone_labels
    data_partial <- bayesplot::mcmc_intervals_data(alpha_partial,prob=.95,transformation=function(x) 1-exp(-x),point_est="mean")
    data_partial$type <- "partial"
    data_global <- rbind(data_global,data_partial)
  }

  plot_labels <- c(plot_labels, plot_labels_partial)
  # print(plot_labels)  
  return (list(data=data_global, plot_labels=plot_labels))

}

#data$type[1] = "First Intervention"
plot_covariate_effects <- function (data, plot_labels, filename) {

  levels(data$parameter) = gsub("t(", "", levels(data$parameter), fixed=TRUE)
  levels(data$parameter) = gsub(")", "", levels(data$parameter), fixed=TRUE)
  data$parameter = (as.character(data$parameter))

  no_point_est <- all(data$point_est == "none")
  x_lim <- range(c(data$ll, data$hh))
  x_range <- diff(x_lim)
  x_lim[1] <- x_lim[1] - 0.05 * x_range
  x_lim[2] <- x_lim[2] + 0.05 * x_range
  layer_vertical_line <- if (0 > x_lim[1] && 0 < x_lim[2]) {
    bayesplot::vline_0(color = "gray90", size = 0.5)
  } else {
    geom_blank(
      mapping = NULL, data = NULL,
      show.legend = FALSE, inherit.aes = FALSE)
  }
  args_outer <- list(mapping = aes_(x = ~ll, xend = ~hh, y = ~parameter, 
                                    yend = ~parameter)) #, color = bayesplot::get_color("mid"))
  args_inner <- list(mapping = aes_(x = ~l, xend = ~h, y = ~parameter, 
                                    yend = ~parameter), size = 2, show.legend = FALSE)
  args_point <- list(mapping = aes_(x = ~m, y = ~parameter), 
                     data = data, size = 4, shape = 21)

  args_point$color <- "blue" #get_color("dark_highlight")

  point_func <- geom_point
  layer_outer <- do.call(geom_segment, args_outer)
  layer_inner <- do.call(geom_segment, args_inner)
  layer_point <- do.call(point_func, args_point)

  data$parameter = factor(as.character(data$parameter),levels=plot_labels[order(plot_labels)[length(plot_labels):1]])
  # data = data[order(-data$m),]
  p = ggplot(data) + ggpubr::theme_pubr() +  geom_point(
      aes(x=m,y=parameter,colour=type),position = position_dodge(-.5)
    ) + geom_linerange(
      aes(xmin=ll,xmax=hh,y=parameter,colour=type), position = position_dodge(-.5)
    ) + 
    scale_x_continuous(
      breaks=seq(0,1,.25), labels = c(
        "0%\n(no effect on transmissibility)","25%","50%","75%","100%\n(ends transmissibility)"
      ), expand=c(0.005,0.005),expression(paste("Relative % reduction in  ",R[t]))
    ) + scale_colour_manual(name = "", values = c(("coral4"), ("seagreen"))) +
    geom_vline(xintercept=1,colour="darkgray") +
    scale_y_discrete("Governmental intervention\n") +
    #geom_vline(xintercept=0,colour="darkgray") + 
    theme(plot.margin = margin(0, 2, 0, .5, "cm"))
  #+ guides(fill=guide_legend(nrow=2))

  ggsave(filename = paste0("results/", filename, "covars-alpha-reduction.png"),
         p,height=4,width=8)
  write.csv(data, paste0("results/", filename, "covars-alpha-reduction.csv"))
  dir.create("web/figures/desktop/", showWarnings = FALSE, recursive = TRUE)
  dir.create("web/figures/mobile/", showWarnings = FALSE, recursive = TRUE)
  cowplot::save_plot(filename = paste0("web/figures/desktop/",  "covars-alpha-reduction.svg"), 
            p, base_height = 4, base_asp = 1.618 * 2 * 8/12)
  cowplot::save_plot(filename = paste0("web/figures/mobile/", "covars-alpha-reduction.svg"), 
            p, base_height = 4, base_asp = 1.1)
  
}

select_covariate_effects()