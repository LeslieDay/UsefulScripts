setwd("~/Desktop/Costa Lab/ExcelFiles/")
growthcurve <- read_excel("011621syntrophysummary.xlsx")

dH2ase_legend = bquote(paste(Delta, italic("H")[2], italic("ase")))
dH2ase_color = bquote("#88a0dc")
dH2ase_shape = bquote (2)

dfdh_legend = bquote(paste(Delta, italic("fdh")))
dfdh_color = bquote("#7c4b73")
dfdh_shape = bquote(16)

WTbarcode1_legend = bquote(paste("WT-barcode1"))
WTbarcode1_color = bquote("#edc775")
WTbarcode1_shape = bquote(23)

yaxis = bquote(paste(OD[600]))

growthcurve %>% 
  filter(media == "Formate" | media == "H2") %>% 
  ggplot(aes(x=time_hours, y=OD600, fill = condition, color = condition)) + 
  geom_point(size =3, aes(shape = condition)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0,50)) + 
  scale_y_log10() +
  geom_errorbar(aes(ymin = (OD600 - sd),
                    ymax = (OD600 + sd))) +
  scale_shape_manual(values = c(dfdh_shape, 17, 23,
                                1, dH2ase_shape,  5),
                     name = "Strain", 
                     labels = c(dfdh_legend, dH2ase_legend, "WT",
                                dfdh_legend, dH2ase_legend, "WT")) +
  scale_colour_manual(aesthetics = c("color", "fill"),
                      name = "Strain", 
                      labels = c(dfdh_legend, dH2ase_legend, "WT",
                                 dfdh_legend, dH2ase_legend, "WT"), 
                      values = c(dfdh_color, dH2ase_color, WTbarcode1_color,
                                 dfdh_color, dH2ase_color, WTbarcode1_color)) +
  theme_classic() +
  labs(y=yaxis, x= "Time (h)") +
  theme(legend.position = c(0.9, 0.7), text = element_text(size = 15))
