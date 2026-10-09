# CODE FOR FIG 5.1

########EastDeltaDesignPlots.R for the PD design#############
# The sole purpose of this code is to plot the Sample Size function and the 
# Conditional Power function of any adaptive SSR trial created in East by 
# the CHW method. (Conditional power utilizes estimated delta, not delta_min)
#############################################################
# We have not developed R code to create the promising zone, obtain the 
# sample size rule inside the promising zone, or to simulate the operating 
# characteristics of designs that utilize conditional power based 
# on estimated delta. 
# The data for such designs is created in East and copied to a .xlsx file 
# and then saved as a .csv file called "Plot.Data-pd.csv". 
# Thus the values of z.L and z.U that define the promising zone are 
# hardcoded into the plotting routine below.
#############################################################

rm(list = ls())

# To prevent R from creating an Rplots.pdf file when this script is executed non-interactively.
options(device = function(...) {
  grDevices::pdf(NULL)
})

library(ggplot2)
library(latex2exp)
library(gridExtra)

## Paths and directories
# Sourcing the helper functions from the tools directory
tryCatch(
  source(file.path("tools", "helper.R"), local = .GlobalEnv),
  error = function(first.error) {
    source(
      file.path("..", "..", "tools", "helper.R"),
      local = .GlobalEnv
    )
  }
)

# Setting up the base directory for executing this code
base.dir <- get.base.dir(
  cur.dir = getwd(),
  sec.dir = "Chapter 5/Sec 5.2.1/"
)

# root <- getwd()
# base.dir <- file.path(root, "Adaptive with delta-est")
# prog.in  <- file.path(base.dir, "EstDeltaDesignPlots-pd.r")

#Read in the SSR and CP data pre-created in East and saved in Plot.Data-pd.csv
SSRdataZ <- read.csv(
  file = file.path(base.dir, "Plot.Data-pd.csv"),
  header = FALSE,
  sep = ","
)
names(SSRdataZ) <- c("z1", "SampSize", "CP")
#Read in z.L and z.U, the start and end of the promising zone, pre-computed in East
z.L <- 1.1298
z.U <- 2.0328

base_ss <- SSRdataZ$SampSize[which.min(abs(SSRdataZ$z1 - 1.13))]

SS.Plot <- ggplot(data = SSRdataZ, aes(x = z1, y = SampSize)) + 
  geom_line(linewidth=1.5) +
  coord_cartesian(xlim = c(-2, 3.5), ylim = c(80, 350)) +
  
  scale_x_continuous(
    breaks = seq(-2, 3.5, by = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(80, 360, by = 20)
  ) +
  
  xlab("Interim Test Statistic z1") + 
  ylab("Sample Size") +
  
  geom_vline(linetype="dashed", color="red", linewidth=1.5, xintercept = z.L) +
  geom_vline(linetype="dashed", color="red", linewidth=1.5, xintercept = z.U) +
  
  geom_segment(
    aes(x = z.L, xend = z.U, y = 170, yend = 170),
    #color = "red",
    linetype = "dotted",
    linewidth = 1.5
  ) +
  
  annotate("text", x=1.6, y=140, label="Promising Zone", size=7) +
  #annotate("text", x=1.6, y=130, label="Zone", size=7) +
  annotate(
    "text",
    label = TeX("$(z.L \\leq z_{1} \\leq z.U)$"),
    x = 1.6, y = 120, size = 7, parse = TRUE
  ) +
  
  theme(
    text = element_text(size = 22),
    
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    
    axis.text.x  = element_text(size = 22),
    axis.text.y  = element_text(size = 22),
    
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )




##############################################################################

CP.Plot <- ggplot(data = SSRdataZ, aes(x = z1, y = CP)) + 
  geom_line(linewidth=1.5) +
  coord_cartesian(xlim = c(-2, 3.5), ylim = c(0, 1)) +
  
  scale_x_continuous(
    breaks = seq(-2, 3.5, by = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1)
  ) +
  
  xlab("Interim Test Statistic z1") + 
 # ylab("Conditional Power") +
  ylab(TeX("Conditional Power Evaluated at $\\delta = 7$")) +
  
  geom_vline(linetype="dashed", color="red", linewidth=1.5, xintercept = z.L) +
  geom_vline(linetype="dashed", color="red", linewidth=1.5, xintercept = z.U) +
  
  annotate("text", x=1.6, y=0.2, label="Promising Zone", size=7) +
  #annotate("text", x=1.6, y=0.17, label="Zone", size=7) +
  annotate(
    "text",
    label = TeX("$(z.L \\leq z_{1} \\leq z.U)$"),
    x = 1.6, y = 0.14, size = 7, parse = TRUE
  ) +
  
  # theme(
  #   text = element_text(size = 16),
  #   
  #   axis.title.x = element_text(size = 25),
  #   axis.title.y = element_text(size = 25),
  #   
  #   axis.text.x  = element_text(size = 20),
  #   axis.text.y  = element_text(size = 20),
  #   
  #   panel.grid.major = element_line(color = "grey80"),
  #   panel.grid.minor = element_line(color = "grey90")
  # )

theme(
  text = element_text(size = 22),
  
  axis.title.x = element_text(size = 26),
  axis.title.y = element_text(size = 26),
  
  axis.text.x  = element_text(size = 22),
  axis.text.y  = element_text(size = 22),
  
  panel.grid.major = element_line(color = "grey80"),
  panel.grid.minor = element_line(color = "grey90")
)

# --- smooth connector inside the promising zone (hand-shaped with a control point) ---
conn_pts <- data.frame(
  x = c(z.L, .5*(z.L+z.U), z.U),     # middle x is the "shape" control
  y = c(0.484,  0.68, 0.803)       # tweak 0.66 up/down to change curvature
)

conn_x <- seq(min(conn_pts$x), max(conn_pts$x), length.out = 200)
conn_y <- spline(conn_pts$x, conn_pts$y, xout = conn_x, method = "natural")$y
conn_df <- data.frame(x = conn_x, y = conn_y)

CP.Plot <- CP.Plot +
  geom_line(
    data = conn_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    #color = "red",
    linetype = "dotted",
    linewidth = 1.5
  )


#arrange the plots in one column

combined_plot <- grid.arrange(SS.Plot, CP.Plot, ncol = 1)

#Save the plot directly in the screen folder for the textbook
# fig.dir <- "D:/Cyrus-Office/BookProject/All-Chapters/screens"
fig.dir <- base.dir

ggsave(
  filename = file.path(fig.dir, "PD_SampleSizeRuleEstDelta.png"),
  plot     = combined_plot,
  width    = 15,
  height   = 15,
  units    = "in",
  dpi      = 300
)







