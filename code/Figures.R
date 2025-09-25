# Figures
library(ragg)
library(ggplot2)
library(grid)

#------Main Text Figures-------

# Figure 1
Figure1 <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_JAS_disagreement, aes(color = TotalConcur), size = 1.5) +
  scale_color_distiller(limits = c(-1,13), palette = "YlOrRd", direction=+1,
                        breaks = c(0, 6, 12),
                        labels = c("Low", "Medium", "High")) +
  labs(color = "Disagreement", fill = NULL, x = "Longitude", y = "Latitude") +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "mm"))+
  geom_point(data = sites, 
             aes(x = LongNEW, y = LatNEW, fill = legend), size = 2)+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.7, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

ggsave(
  filename = "figure1.png",
  plot = Figure1,                        # your ggplot object
  width = 170, 
  height = 80,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# Figure2
Figure2 <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_ave_LST_Predict_long, aes(color = suitability), size = .5) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Habitat Suitability",        # Title for the legend
                        breaks = c(0.1, .45, 0.8),
                        labels = c("Low", "Medium", "High")) +  # Suitability color scale
  theme_minimal() +
  facet_wrap(~suitability_type)+
  labs(x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 10))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.5, y = 0.1, hjust = 0, gp = gpar(fontsize = 8)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

ggsave(
  filename = "figure2.png",
  plot = Figure2,                        # your ggplot object
  width = 170, 
  height = 90,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# Figure3
Figure3 <- ggplot(data2_summary, aes(x = Month, y = mean_value, fill = suitability_type))+
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9), color = "black")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")
  )+
  scale_fill_viridis_d(option = "viridis")+
  geom_text(aes(label = groups, group = suitability_type), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 4)+
  labs( x = "Month", y = "Average Habitat Suitability", fill = "Model") 

ggsave(
  filename = "figure3.png",
  plot = Figure3,                        # your ggplot object
  width = 170, 
  height = 120,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# Figure4
LogA <- ggplot(data_summary, aes(x = Month, y = -mean_value, color = Model, group = Model)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = groups, 
                y = -mean_value + ifelse(Model == "BRT", 0.032,
                                         ifelse(Model == "GAM", 0.022,
                                                ifelse(Model == "MaxEnt", -0.0182, -0.032)))),size = 3)+
  labs(x = "Month",
       y = "Model Performance (negative log score)") +
  scale_x_discrete(labels = c("7" = "July", "8" = "August","9" = "September")) +
  scale_color_viridis_d(option = "viridis", labels = c("BRT", "GAM", "MaxEnt", "RF")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13))

LogB <- ggplot(over_under_long, aes(x = Month, y = Bias, fill = Model)) +
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.9), color = "black") +
  labs(x = "Month", y = "Model Bias") +
  scale_fill_viridis_d(option = "viridis") +
  coord_cartesian(clip = "off") +   # allow legend to sit near edges
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 13),
    axis.text.y  = element_text(size = 8),
    axis.title.y = element_text(size = 13),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    
    # <<< place legend inside panel b >>>
    legend.position      = c(0.98, 0.25),   # (x, y) in [0,1] panel coords (right, near top)
    legend.justification = c(1, 1),         # anchor legend's top-right corner
    legend.background    = element_rect(fill = "white", color = "black"),
    legend.key.size      = unit(4, "mm"),
    legend.text          = element_text(size = 10),
    legend.title         = element_text(size = 13)
  )

Figure4 <- ggarrange(LogA, LogB, ncol = 2, labels = c("a)", "b)"))

ggsave(
  filename = "figure4.png",
  plot = Figure4,                        # your ggplot object
  width = 170, 
  height = 135,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# Figure5
Figure5 <- ggplot()+
  geom_errorbar(data = invasion_summary, aes(ymin = -(mean - se), ymax = -(mean + se), x = TrainingYrs, group = Model, color = Model), width = 0.2, stat = "identity")+
  geom_line(data = invasion_summary, aes(x = TrainingYrs, y = -mean, group = Model, color = Model))+
  theme_minimal()+
  scale_color_viridis_d(option = "viridis", labels = c("BRT","GAM", "MaxEnt", "RF"))+
  labs(x = "Training data", y = "Model Performance (negative log score)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position      = c(0.98, 0.25),   # (x, y) in [0,1] panel coords (right, near top)
        legend.justification = c(1, 1),         # anchor legend's top-right corner
        legend.background    = element_rect(fill = "white", color = "black"),
        legend.key.size      = unit(4, "mm"),
        legend.text          = element_text(size = 10),
        legend.title         = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))+
  scale_x_continuous(breaks = invasion_summary$TrainingYrs, labels = invasion_summary$TrainingYrs)

ggsave(
  filename = "figure5.png",
  plot = Figure5,                        # your ggplot object
  width = 170, 
  height = 135,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)



#------Supplemental Text Figures-------

# S1Figure
S1Figure <- ggplot()+
  geom_sf(data=boundary_AEA,aes(geometry=geometry), fill="grey90",color="black")+
  geom_point(data=sampling_ALB,
             aes(x = Longitude.y, y = Latitude.y, color = n_years), size = 1) +
  scale_color_distiller(limits = c(0,17), palette = "YlOrRd", direction=+1,
                        breaks = c(1, 4, 8, 12, 16),
                        labels = c("1", "4", "8", "12", "16")) +
  labs(x = "Longitude", y = "Latitude", color = "Sampling Years")+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "mm"))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.7, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

ggsave(
  filename = "S1Figure.png",
  plot = S1Figure,                        # your ggplot object
  width = 170, 
  height = 80,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# S2Figure
S2Figure <- ggplot()+
  geom_sf(data=boundary_AEA,aes(geometry=geometry), fill="grey90",color="black")+
  geom_point(data=first_ALB, 
             aes(x = Longitude.y, y = Latitude.y, color = Year), size = 1.5) + 
  scale_color_viridis_d(option = "viridis")+
  labs(x = "Longitude", y = "Latitude")+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.spacing.y  = unit(0.5, "mm"), legend.key.size   = unit(3, "mm"))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.7, y = 0.1, hjust = 0, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

ggsave(
  filename = "S2Figure.png",
  plot = S2Figure,                        # your ggplot object
  width = 170, 
  height = 90,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# S4Figure
IMP <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_imp), size = 1.25) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "% Impervious Surface") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
# DLST
DLST <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_DLST), size = 1.25) +
  scale_color_viridis_c(option = "inferno", direction = 1, name = "Day Surface Temperature") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
# NLST
NLST <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_NLST), size = 1.25) +
  scale_color_viridis_c(option = "inferno", direction = 1, name = "Night Surface Temperature") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))+
  annotation_custom(
    grob = textGrob("Source: gis.ny.gov", x = 0.6, y = 0.1, hjust = 0, gp = gpar(fontsize = 8)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
# EVI
EVI <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_EVI), size = 1.25) +
  scale_color_viridis_c(option = "mako", direction = 1, name = "EVI                                       ") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Suitability Score") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
# Land cover
LC <- ggplot() +
  geom_sf(data = boundary_AEA, fill = "grey90", color = "black", size = 0.2) +  # Background shape
  geom_sf(data = spatial_LST_Predict_predictors, aes(color = ave_LC), size = 1.25) +
  # scale_color_viridis(option = "plasma", name = "Landcover") +  # Suitability color scale
  theme_minimal() +
  labs( x = "Longitude", y = "Latitude", color = "Land Cover") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

S4Figure <- ggarrange(IMP, LC, EVI, DLST, NLST, ncol = 1, nrow = 5, labels = c("a)","b)","c)","d)","e)"))

ggsave(
  filename = "S4Figure.png",
  plot = S4Figure,                        # your ggplot object
  width = 170, 
  height = 340,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# S5Figure
S5Figure <- ggplot(ROCs_random,aes(x=FPR, y=TPR, color=Run, linetype = Data))+
  geom_line(alpha=.5,linewidth=.25)+
  geom_smooth(color="black")+
  geom_abline(slope=1,intercept = 0,linetype="dotted")+
  geom_label(data=meanAUC_random %>% filter(Data == "Test"),aes(x=.7,y=0.05,label=paste("Mean Testing AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  geom_label(data=meanAUC_random %>% filter(Data == "Train"),aes(x=.7,y=0.15,label=paste("Mean Training AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  scale_color_discrete(guide=NULL)+
  facet_wrap(Model~.)+
  xlab("False Positivity Rate")+
  ylab("True Positivity Rate")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 13))

ggsave(
  filename = "S5Figure.png",
  plot = S5Figure,                        # your ggplot object
  width = 170, 
  height = 150,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# S6Figure
S6Figure <- ggplot(ROCs_Lambda,aes(x=FPR, y=TPR, color=Run, linetype = Data))+
  geom_line(alpha=.5,linewidth=.25)+
  geom_smooth(color="black")+
  geom_abline(slope=1,intercept = 0,linetype="dotted")+
  geom_label(data=meanAUC_lambda %>% filter(Data == "Test"),aes(x=.7,y=0.05,label=paste("Mean Testing AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  geom_label(data=meanAUC_lambda %>% filter(Data == "Train"),aes(x=.7,y=0.15,label=paste("Mean Training AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  scale_color_discrete(guide=NULL)+
  xlab("False Positivity Rate")+
  ylab("True Positivity Rate")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 13))

ggsave(
  filename = "S6Figure.png",
  plot = S6Figure,                        # your ggplot object
  width = 170, 
  height = 150,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# S7Figure
S7Figure <- ggplot(ROCs_township,aes(x=FPR, y=TPR, color=Run, linetype = Data))+
  geom_line(alpha=.5,linewidth=.25)+
  geom_smooth(color="black")+
  geom_abline(slope=1,intercept = 0,linetype="dotted")+
  geom_label(data=meanAUC_township %>% filter(Data == "Test"),aes(x=.7,y=0.05,label=paste("Mean Testing AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  geom_label(data=meanAUC_township %>% filter(Data == "Train"),aes(x=.7,y=0.15,label=paste("Mean Training AUC = ",round(meanAUC,digits = 3))),
             color="black",size=3)+
  scale_color_discrete(guide=NULL)+
  facet_wrap(Model~.)+
  xlab("False Positivity Rate")+
  ylab("True Positivity Rate")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 13))

ggsave(
  filename = "S7Figure.png",
  plot = S7Figure,                        # your ggplot object
  width = 170, 
  height = 150,       # mm (double column wide, good aspect ratio)
  units = "mm",
  dpi = 600,                       # print-ready
  device = ragg::agg_png,          # crisp lines and text
  bg = "white"
)

# S8Figure
# S9Figure
