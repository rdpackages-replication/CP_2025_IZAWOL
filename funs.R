# function that computes the mean control outcome given a bandwidth
meanControlsGet <- function(y, x, c, h, z=NULL, subset=NULL) {
  not.missing <- !is.na(y) & !is.na(x)
  
  if (!is.null(z)) not.missing <- not.missing & rowSums(is.na(z)) == 0
  if (!is.null(subset)) not.missing <- not.missing & subset
  
  within.h.left <- x >= c - h & x < c
  y0 <- mean(y[not.missing & within.h.left])
  
  return(y0)
}


# ggplot code to generate figures
plotsGet <- function(toplot.bins, toplot.poly, path.fig, h.opt) {

  # graphic options
  text.size <- 14L
  line.wdt <- 1.8
  scatter.alpha <- 0.2
  dot.size <- 1.5
  my_palette <- wesanderson::wes_palette(n=2L, name="Darjeeling1")
  
  ######################################
  # Figure 1a) canonical RD: global plot
  ######################################

  p <- ggplot() +
    geom_line(data=subset(toplot.poly, group=="all units" & type=="global"),
              aes(x=poly_x, y=poly_y, group=treated, color=group), linewidth=line.wdt) +
    geom_point(data=subset(toplot.bins, group=="all units" & type=="global"),
               aes(x=bin_mean_x, y=bin_mean_y, group=treated, color=group),
               alpha=3*scatter.alpha, size=dot.size) +
    geom_vline(aes(xintercept = cutoff)) +
    scale_color_manual(values=c("black"), name=NULL) +
    xlab("poverty rate") + ylab("child mortality") +
    xlim(40, 70) +  ylim(0, 6) +
    theme(legend.position="none",
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size),
          axis.title.y = element_text(size = text.size),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  
  suppressWarnings(
    ggsave(paste0(path.fig, "rdplot_All.png"), plot=p, width=6, height=4, dpi="retina")
  )
  
  #####################################
  # Figure 1b) canonical RD: local plot
  #####################################
  p <- ggplot() +
    geom_point(data=subset(toplot.bins, group=="all units"),
               aes(x=bin_mean_x, y=bin_mean_y, group=treated, color=group, alpha=within.h),
               size=dot.size) +
    geom_line(data=subset(toplot.poly, group=="all units" & type=="local"),
              aes(x=poly_x, y=poly_y, group=treated, color=group), linewidth=line.wdt) +
    geom_vline(aes(xintercept = cutoff)) +
    geom_vline(aes(xintercept = cutoff + h.opt[["all"]]), linetype="dashed") + 
    geom_vline(aes(xintercept = cutoff - h.opt[["all"]]), linetype="dashed") + 
    geom_segment(
      aes(x = cutoff, xend = cutoff - h.opt[["all"]], y = 6, yend = 6),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")), 
      linewidth = 0.5, color = "black"
    ) +
    geom_segment(
      aes(x = cutoff - h.opt[["all"]], xend = cutoff, y = 6, yend = 6),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")), 
      linewidth = 0.5, color = "black"
    ) +
    annotate("text", x=cutoff - h.opt[["all"]]/2, y=5.6, label="h") +
    scale_color_manual(values=c("black"), name=NULL) +
    scale_alpha_manual(values=c(scatter.alpha,1)) +
    xlim(40, 70) + ylim(0, 6) +
    xlab("poverty rate") + ylab("child mortality") +
    theme(legend.position="none",
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size),
          axis.title.y = element_text(size = text.size),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  
  suppressWarnings(
    ggsave(paste0(path.fig, "rdplot_All_zoom.png"), plot=p, width=6, height=4, dpi="retina")
  )
  
  ##########################################
  # Figure 2a) heterogeneous RD: global plot
  ##########################################
  p <- ggplot() +
    geom_point(data=subset(toplot.bins, group!="all units" & type=="global"),
               aes(x=bin_mean_x, y=bin_mean_y, color=group, group=interaction(group, treated)),
               alpha=3*scatter.alpha, size=dot.size) +
    geom_line(data=subset(toplot.poly, group!="all units" & type=="global"),
              aes(x=poly_x, y=poly_y, color=group, group=interaction(treated, group)), linewidth=line.wdt) +
    geom_vline(aes(xintercept = cutoff)) +
    guides(colour = guide_legend(nrow=1L, byrow=TRUE)) +
    scale_color_manual(values=my_palette, name=NULL) +
    xlab("poverty rate") + ylab("child mortality") +
    xlim(40, 70) + ylim(0, 7) +
    theme(legend.position="bottom",
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size),
          axis.title.y = element_text(size = text.size),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.title=element_text(size=text.size),
          legend.text=element_text(size=text.size))  
  
  suppressWarnings(
    ggsave(paste0(path.fig, "rdplot_Hte.png"), plot=p, width=6, height=4, dpi="retina")
  )
  
  #########################################
  # Figure 2b) heterogeneous RD: local plot
  #########################################
  p <- ggplot() +
    geom_point(data=subset(toplot.bins, group!="all units"),
               aes(x=bin_mean_x, y=bin_mean_y, color=group,
                   group=interaction(group, treated), alpha=within.h),
               size=dot.size) +
    geom_line(data=subset(toplot.poly, group!="all units" & type=="local"),
              aes(x=poly_x, y=poly_y, color=group, group=interaction(treated, group)),
              linewidth=line.wdt) +
    geom_vline(aes(xintercept = cutoff)) +
    guides(colour = guide_legend(nrow=1L, byrow=TRUE)) +
    scale_color_manual(values=my_palette, name=NULL) +
    scale_alpha_manual(values=c(scatter.alpha,1)) +
    guides(alpha="none") +
    xlab("poverty rate") + ylab("child mortality") +
    xlim(40, 70) + ylim(0, 7) +
    geom_segment(
      aes(x = cutoff, xend = cutoff - h.opt[["w0"]], y = 6.75, yend = 6.75),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")), 
      linewidth = 0.5, color = my_palette[1]
    ) +
    geom_segment(
      aes(x = cutoff - h.opt[["w0"]], xend = cutoff, y = 6.75, yend = 6.75),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")), 
      linewidth = 0.5, color = my_palette[1]
    ) +
    annotate("text", x=cutoff - h.opt[["w0"]]/2, y=7, label=TeX("$h_0$"), color=my_palette[1]) +
    geom_vline(aes(xintercept = cutoff + h.opt[["w0"]]), linetype="dashed", color=my_palette[1]) + 
    geom_vline(aes(xintercept = cutoff - h.opt[["w0"]]), linetype="dashed", color=my_palette[1]) + 
    geom_segment(
      aes(x = cutoff, xend = cutoff - h.opt[["w1"]], y = 6.55, yend = 6.55),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")), 
      linewidth = 0.5, color = my_palette[2]
    ) +
    geom_segment(
      aes(x = cutoff - h.opt[["w1"]], xend = cutoff, y = 6.55, yend = 6.55),
      arrow = arrow(type = "closed", length = unit(0.05, "inches")), 
      linewidth = 0.5, color = my_palette[2]
    ) +
    annotate("text", x=cutoff - h.opt[["w0"]]/2, y=6.3, label=TeX("$h_1$"), color=my_palette[2]) +
    geom_vline(aes(xintercept = cutoff + h.opt[["w1"]]), linetype="dashed", color=my_palette[2]) + 
    geom_vline(aes(xintercept = cutoff - h.opt[["w1"]]), linetype="dashed", color=my_palette[2]) + 
    theme(legend.position="bottom",
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size),
          axis.title.y = element_text(size = text.size),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.title=element_text(size=text.size),
          legend.text=element_text(size=text.size))  
  suppressWarnings(
    ggsave(paste0(path.fig, "rdplot_Hte_zoom.png"), plot=p, width=6, height=4, dpi="retina")
  )
}
