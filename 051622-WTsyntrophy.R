library(tidyverse)
growth.data <- read_excel("~/Desktop/Costa Lab/ExcelFiles/051922Wt-syntrophy.xlsx")
growth.summary <- growth.data %>% filter(condition != "ethanolamine_WT-syntrophy" | replicate != "B" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "E" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "F" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "G" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "H" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "I" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "J" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "K" &
                                           condition != "ethanolamine_WT-syntrophy" | replicate != "L" &
                                           condition != "ethanolamine_Scarbinolica" | replicate != "D" &
                                           condition != "acetoin_Scarbinolica" | replicate != "D" &
                                           condition != "ethanol_WT-syntrophy" | replicate != "D" &
                                           condition != "13PD_WT-syntrophy" | replicate != "D") %>% 
  group_by(strain, substrate,condition, hours) %>% 
  summarise(mean = mean(OD600), sd = sd(OD600))

growth.data2 <- read_excel("~/Desktop/Costa Lab/ExcelFiles/042122WtSynGrowth.xlsx")
growth.summary2 <- growth.data2 %>% 
  group_by(strain, substrate,condition, hours) %>% 
  summarise(mean = mean(OD600), sd = sd(OD600))

all.summary <- merge(growth.summary, growth.summary2, all = TRUE)

yaxis = bquote(paste(OD[600]))


#move legend 

growth.summary %>% filter(strain == "WT-syntrophy" & substrate != "acetoin") %>%
  ggplot(aes(x=hours, y=mean, group= condition, color = condition)) + 
  geom_point(aes(shape = condition, color = condition, fill = condition), size =3) +
  geom_line(size=1) +
  geom_errorbar(aes(hours, mean, ymin= mean-sd, ymax = mean+sd)) +
  scale_y_log10() +
  theme_classic() + 
  scale_x_continuous(limits=c(0,25)) +
  labs(y=yaxis, x= "Time (h)") + 
  scale_color_manual(breaks = c("ethanol_WT-syntrophy",
                                "13PD_WT-syntrophy",
                                "ethanolamine_WT-syntrophy"),
                     values = c("#2c6b67",
                                "#b0799a",
                                "#88a0dc"),
                     name = "Substrate",
                     labels = c("Ethanol",
                                "1,3-Propanediol",
                                "Ethanolamine")) +
  scale_fill_manual(breaks = c("ethanol_WT-syntrophy",
                               "13PD_WT-syntrophy",
                               "ethanolamine_WT-syntrophy"),
                    values = c("#2c6b67",
                               "#b0799a",
                               "#88a0dc"),
                    name = "Substrate",
                    labels = c("Ethanol",
                               "1,3-Propanediol",
                               "Ethanolamine"))+
  scale_shape_manual(values=c(0,
                              1,
                              2), 
                     breaks = c("ethanol_WT-syntrophy",
                                "13PD_WT-syntrophy",
                                "ethanolamine_WT-syntrophy"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Substrate", 
                     labels = c("Ethanol",
                                "1,3-Propanediol",
                                "Ethanolamine")) +
  theme(legend.position = c(0.1, 0.8), text = element_text(size = 15))

all.summary %>% 
  ggplot(aes(x=hours, y=mean, group= condition, color = condition)) + 
  geom_point() +
  geom_line() +
  geom_errorbar(aes(hours, mean, ymin= mean-sd, ymax = mean+sd)) + 
  scale_y_log10()  +
  theme(legend.position = c(0.1, 0.8), text = element_text(size = 15))

all.summary %>% filter(condition != "acetoin_WT-syntrophy") %>%
  ggplot(aes(x=hours, y=mean, group= condition, color = condition)) + 
  geom_point(aes(shape = condition, color = condition, fill = condition), size =3) +
  geom_line(size=1) +
  geom_errorbar(aes(hours, mean, ymin= mean-sd, ymax = mean+sd)) +
  scale_y_log10() +
  theme_classic() + 
  scale_x_continuous(limits=c(0,25)) +
  labs(y=yaxis, x= "Time (h)") + 
  scale_color_manual(breaks = c("ethanol_WT-syntrophy",
                                "13PD_WT-syntrophy",
                                "acetoin_Scarbinolica",
                                "ethanolamine_WT-syntrophy",
                                "ethanolamine_Scarbinolica",
                                "P.carb_EtOH",
                                "P.carb_13PD",
                                "MM901_EtOHNH2",
                                "MM901_EtOH",
                                "MM901_13PD"),
                     values = c("#2c6b67",
                                "#b0799a",
                                "#7c4b73",
                                "#88a0dc",
                                "#88a0dc",
                                "#2c6b67",
                                "#b0799a",
                                "#88a0dc",
                                "#2c6b67",
                                "#b0799a"),
                     name = "Condition",
                     labels = c("Syntrophy-Ethanol",
                                "Syntrophy-1,3-Propanediol",
                                "S.carbinolica-Acetoin",
                                "Syntrophy-Ethanolamine",
                                "S.carbinolica-Ethanolamine",
                                "S.carbinolica-Ethanol",
                                "S.carbinolica-1,3-Propanediol",
                                "M.maripaludis-Ethanolamine",
                                "M.maripaludis-Ethanol",
                                "M.maripaludis-1,3-Propanediol")) +
  scale_fill_manual(breaks = c("ethanol_WT-syntrophy",
                               "13PD_WT-syntrophy",
                               "acetoin_Scarbinolica",
                               "ethanolamine_WT-syntrophy",
                               "ethanolamine_Scarbinolica",
                               "P.carb_EtOH",
                               "P.carb_13PD",
                               "MM901_EtOHNH2",
                               "MM901_EtOH",
                               "MM901_13PD"),
                    values = c("#2c6b67",
                               "#b0799a",
                               "#7c4b73",
                               "#88a0dc",
                               "#88a0dc",
                               "#2c6b67",
                               "#b0799a",
                               "#88a0dc",
                               "#2c6b67",
                               "#b0799a"),
                    name = "Condition",
                    labels = c("Syntrophy-Ethanol",
                               "Syntrophy-1,3-Propanediol",
                               "S.carbinolica-Acetoin",
                               "Syntrophy-Ethanolamine",
                               "S.carbinolica-Ethanolamine",
                               "S.carbinolica-Ethanol",
                               "S.carbinolica-1,3-Propanediol",
                               "M.maripaludis-Ethanolamine",
                               "M.maripaludis-Ethanol",
                               "M.maripaludis-1,3-Propanediol"))+
  scale_shape_manual(values=c(0,
                              0,
                              1,
                              0,
                              1,
                              1,
                              1,
                              2,
                              2,
                              2), 
                     breaks = c("ethanol_WT-syntrophy",
                                "13PD_WT-syntrophy",
                                "acetoin_Scarbinolica",
                                "ethanolamine_WT-syntrophy",
                                "ethanolamine_Scarbinolica",
                                "P.carb_EtOH",
                                "P.carb_13PD",
                                "MM901_EtOHNH2",
                                "MM901_EtOH",
                                "MM901_13PD"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Condition", 
                     labels = c("Syntrophy-Ethanol",
                                "Syntrophy-1,3-Propanediol",
                                "S.carbinolica-Acetoin",
                                "Syntrophy-Ethanolamine",
                                "S.carbinolica-Ethanolamine",
                                "S.carbinolica-Ethanol",
                                "S.carbinolica-1,3-Propanediol",
                                "M.maripaludis-Ethanolamine",
                                "M.maripaludis-Ethanol",
                                "M.maripaludis-1,3-Propanediol")) +
  theme( text = element_text(size = 15))

Scarbinolica = bquote(paste(italic(S. ), italic(carbinolica)))
Mmaripaludis = bquote(paste(italic(M. ), italic(maripaludis)))
all.summary %>% filter(condition != "acetoin_WT-syntrophy") %>%
  ggplot(aes(x=hours, y=mean, group= condition, color = substrate)) + 
  geom_point(aes(shape = strain, color = substrate, fill = substrate), size =3) +
  geom_line(size=1) +
  geom_errorbar(aes(hours, mean, ymin= mean-sd, ymax = mean+sd)) +
  scale_y_log10() +
  theme_classic() + 
  scale_x_continuous(limits=c(0,25)) +
  labs(y=yaxis, x= "Time (h)") + 
  scale_color_manual(aesthetics = c("color", "fill"),
                     values = cols,
                     name = "Substrate",
                     label = c("Ethanol",
                     "1,3-Propanediol",
                     "Ethanolamine",
                     "Acetoin"))+
  scale_shape_manual(values = shapes,
                     name = "Strain",
                     label = c("Syntrophy",
                               Scarbinolica,
                               Mmaripaludis)) +
  theme(text = element_text(size = 15))

cols <- c("ethanol" = "#2c6b67",  "13PD" = "#b0799a", "ethanolamine" = "#88a0dc", "acetoin" = "black")
shapes <- c("WT-syntrophy" = 0, "Scarbinolica" = 1, "MM901" = 2)
