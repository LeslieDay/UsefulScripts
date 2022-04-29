# UsefulScripts
Place for little scripts I wrote to help with different data analysis / calculations. 


Use GenerationCalculator.py to calculate number of generations/doublings of a culture.
  - provide a tab delimited file containing inoculation OD, dilution factor, and final OD for 1 or more samples 
  - Example file.txt
 <img width="655" alt="Screen Shot 2022-02-02 at 4 49 37 PM" src="https://user-images.githubusercontent.com/92818902/152250830-27cbcbe2-bfa1-4be6-9831-da69346abe87.png">
 
 - to run 
 ```bash
 $python GenerationCalculator2.py file.txt 
 ```

R Script to generate ![image](https://user-images.githubusercontent.com/92818902/165995315-927f4e97-a130-4b21-a7a3-ed644c46a948.png)

```R
# script for analyzing fdh and H2 axenic growth on formate and H2 
# 011621syntrophygrowth file 
setwd("~/Desktop/Costa Lab/ExcelFiles/")
growthcurve <- read_excel("011621syntrophysummary.xlsx")

growthcurve %>% 
  filter(media == "Formate" | media == "H2") %>% 
  ggplot(aes(x=time_hours, y=OD600, fill = condition, color = condition)) + 
  geom_point(size =3, aes(shape = condition)) +
  geom_line(size=1 ) +
  scale_x_continuous(limits = c(0,50)) + 
  scale_y_log10() +
  geom_errorbar(aes(ymin = (OD600 - sd),
                    ymax = (OD600 + sd)),
                size = 1, color="black",
                width = 0.15) +
  scale_shape_manual(values = c(dfdh_shape, 17, 23,
                                1, dH2ase_shape,  5),
                     name = "Strain", 
                     labels = c(dfdh_legend, dH2ase_legend, "WT",
                                dfdh_legend, dH2ase_legend, "WT")) +
  scale_colour_manual(name = "Strain", 
                      labels = c(dfdh_legend, dH2ase_legend, "WT",
                                 dfdh_legend, dH2ase_legend, "WT"), 
                      values = c(dfdh_color, dH2ase_color, WTbarcode1_color,
                                 dfdh_color, dH2ase_color, WTbarcode1_color)) +
  scale_fill_manual(name = "Strain",
                    values = c(dfdh_color, dH2ase_color, WTbarcode1_color,
                               dfdh_color, dH2ase_color, WTbarcode1_color),
                    labels = c(dfdh_legend, dH2ase_legend, "WT",
                               dfdh_legend, dH2ase_legend, "WT")) + 
  theme_classic() +
  labs(x = "Time(h)")
```
