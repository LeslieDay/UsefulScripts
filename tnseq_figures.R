library(tidyverse)
library(ggthemes)
library(readxl)
#import all samples data from Desktop
all_samples <- read_excel("~/Desktop/Costa Lab/tnseq/tnseq_01/sample_data/all_samples.xlsx")

#list of essential genes from Sarmiento et. al. 2013
#make file with just genes that have EI from 1-3 in either lib1 or lib2 of rich medium McCm_T2 =556 genes 
EssentialGenes <- read_excel("Desktop/Costa Lab/tnseq/tnseq_01/essential_genes.xlsx")
McCm_T2_essential <- EssentialGenes %>% filter(lib1_McCm_T2 < 4 | lib2_McCm_T2 < 4)


#keep only CDS data 
CDS_tndata <- filter(all_samples, feature_type == "CDS")

S1_RM_0 = 19166143
S2_RM_4 = 11491039
S3_RM_9 = 11328098
S4_RM_19 = 10573408
S5_RM_28 = 10886591

#calculate fold change 
tndata_FC <- all_samples %>% mutate(RPKM_0 = RPK_gen_0/S1_RM_0*1000000,
                                   RPKM_4 = RPK_gen_4/S2_RM_4*1000000,
                                   RPKM_9 = RPK_gen_9/S3_RM_9*1000000,
                                   RPKM_19 = RPK_gen_19/S4_RM_19*1000000,
                                   RPKM_28 = RPK_gen_28/S5_RM_28*1000000,
                      RPKM_log2FC_4 = log2(RPKM_4/RPKM_0), 
                      RPKM_log2FC_9 = log2(RPKM_9/RPKM_0),
                      RPKM_log2FC_19 = log2(RPKM_19/RPKM_0),
                      RPKM_log2FC_28 = log2(RPKM_28/RPKM_0),
                      RPKM_log2FC_0 = 0)

tndata_long<- tndata_FC %>% pivot_longer(
  cols = CPM_gen_0:RPKM_log2FC_0,
  names_to = c("count_type", "sample"),
  names_pattern = "(...+_.*)_(.*)",
  values_to = "values",
  names_transform = list(sample = as.integer)
)
# remove all experimentally determined essential genes (known_essential_genes_S2) 
# and genes determined essential in lib 1&2 rich medium (McCm_T2_essential)
tn_nonessent <- tndata_FC %>% 
  filter(old_locus_tag %in% EssentialGenes$known_essential_genes_S2 == FALSE) %>% 
  filter(old_locus_tag %in% McCm_T2_essential$locus_tag == FALSE)



tndata_long %>% filter(count_type == "RPKM_log2FC") %>% 
  ggplot() +
  geom_line(aes(x=sample, y=values, fill = old_locus_tag),size = 0.05, alpha=1) +
  labs(x = "Generation", y = "Log2 Standardized Reads per Gene") +
  scale_x_continuous(breaks=c(0, 4, 9, 19, 28), limits=c(0,28)) +
  annotate("rect", xmin=0,xmax=28,ymin=c(-0.5, 0.5, -6),
           ymax=c(0.5, 2, -0.5), alpha= 0.1, fill=c("black", "blue", "red")) + 
  theme_classic()





###key genes and MMP numbers
# ehaO = MMP1462 = energy conserving hydrogenase A large subunit 
# fdhA = MMP0138 
# fdhB = MMP0139
# ehbN = MMP1153
# fruA = MMP1382
# frcA = MMP0820
# hmd  = MMP0127
# mtd = MMP0372
# vhuA = MMP1694
# vhcA = MMP0823

tndata_FC %>% filter(old_locus_tag == "MMP1297"|
                       old_locus_tag == "MMP1298") %>% 
  ggplot(aes(x=old_locus_tag,y=RPKM_log2FC_9,  fill = old_locus_tag)) + geom_col()

tndata_long %>% filter(count_type == "RPKM_log2FC") %>% filter(old_locus_tag == "MMP0372"|
                  old_locus_tag == "MMP1297"|
                  old_locus_tag == "MMP1298"|
                  old_locus_tag == "MMP1462"|
                  old_locus_tag == "MMP1153"|
                  old_locus_tag == "MMP1382"|
                  old_locus_tag == "MMP0820"|
                  old_locus_tag == "MMP1694"|
                  old_locus_tag == "MMP0823"|
                  old_locus_tag == "MMP0138") %>% 
  ggplot(aes(x=sample,y=values,  fill = old_locus_tag)) +
geom_point(aes(shape = old_locus_tag)) + 
  geom_line() +
  scale_color_discrete(breaks = c("MMP0138", "MMP0823", "MMP1382", "MMP0820", "MMP1694",
                                  "MMP1462", "MMP1153", "MMP1297","MMP1298", "MMP0372"),
                       name = "Gene",
                       labels = c("formate dehydrogenase A1", "vhcA", "fruA", "frcA", "vhuA",
                                  "ehaO", "ehbN", "formate dehydrogenase B2","formate dehydrogenase A2", "mtd")) +
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14,15, 16), 
                     breaks = c("MMP0138", "MMP0823", "MMP1382", "MMP0820", "MMP1694",
                                "MMP1462", "MMP1153", "MMP1297","MMP1298", "MMP0372"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Gene", 
                     labels = c("formate dehydrogenase A1", "vhcA", "fruA", "frcA", "vhuA",
                                "ehaO", "ehbN", "formate dehydrogenase B2","formate dehydrogenase A2", "mtd"))

# ehaO = MMP1462 = energy conserving hydrogenase A large subunit 
# fdhA1 = MMP0138 
# fdhB1 = MMP0139
# ehbN = MMP1153
# fruA = MMP1382
# frcA = MMP0820
# hmd  = MMP0127
# mtd = MMP0372
# vhuA = MMP1694
# vhcA = MMP0823

# fruA = MMP1382
# frcA = MMP0820
# vhuA = MMP1694
# vhcA = MMP0823


tndata_long %>% filter(count_type == "RPKM_log2FC") %>% filter(old_locus_tag == "MMP1297"|
                                                                 old_locus_tag == "MMP1298"|
                                                                 old_locus_tag == "MMP1382"|
                                                                 old_locus_tag == "MMP0820"|
                                                                 old_locus_tag == "MMP1694"|
                                                                 old_locus_tag == "MMP0823") %>% 
  ggplot(aes(x=sample,y=values,  color = old_locus_tag)) +
  geom_point(aes(shape = old_locus_tag, color = old_locus_tag, size =0.3)) + 
  geom_line(size=1) +
  labs(y="Log2 Fold Change Hits per Gene", x= "Generation")+
  scale_color_discrete(breaks = c("MMP0823", "MMP1382", "MMP0820", "MMP1694", "MMP1297","MMP1298"),
                       name = "Gene",
                       labels = c("hydrogenase vhcA", "hydrogenase fruA", "hydrogenase frcA",  
                                   "hydrogenase vhuA", "formate dehydrogenase B2", "formate dehydrogenase A2")) +
  scale_shape_manual(values=c(15,16,17,18,1,2), 
                     breaks = c("MMP0823", "MMP1382", "MMP0820", "MMP1694", "MMP1297","MMP1298"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Gene", 
                     labels = c("hydrogenase vhcA", "hydrogenase fruA", "hydrogenase frcA", 
                                "hydrogenase vhuA", "formate dehydrogenase B2", "formate dehydrogenase A2")) +
  scale_x_continuous(breaks=c(0, 4, 9, 19, 28), limits=c(0,28)) +
  annotate("rect", xmin=0,xmax=28,ymin=c(-0.5, 0.5, -6),
           ymax=c(0.5, 2, -0.5), alpha= 0.1, fill=c("black", "blue", "red")) + 
  theme_classic()

tndata_long %>% filter(count_type == "RPKM_log2FC") %>% filter(old_locus_tag == "MMP1297"|
                                                                 old_locus_tag == "MMP1298"|
                                                                 old_locus_tag == "MMP1382"|
                                                                 old_locus_tag == "MMP0820"|
                                                                 old_locus_tag == "MMP1694"|
                                                                 old_locus_tag == "MMP0823"|
                                                                 old_locus_tag == "MMP0127"|
                                                                 old_locus_tag == "MMP0372") %>% 
  ggplot(aes(x=sample,y=values,  color = old_locus_tag)) +
  geom_point(aes(shape = old_locus_tag, color = old_locus_tag, size =0.3)) + 
  geom_line(size=1) +
  labs(y="Log2 Fold Change Hits per Gene", x= "Generation")+
  scale_color_discrete(breaks = c("MMP0127","MMP0823", "MMP1382", "MMP0820", "MMP1694", "MMP1297","MMP1298","MMP0372"),
                       name = "Gene",
                       labels = c(hmd,vhcA, fruA, frcA, vhuA, fdhB, fdhA, mtd)) +
  scale_shape_manual(values=c(15,16,17,18,16,2,1,2), 
                     breaks = c("MMP0127","MMP0823", "MMP1382", "MMP0820", "MMP1694", "MMP1297","MMP1298","MMP0372"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Gene", 
                     labels = c(hmd,vhcA, fruA, frcA, vhuA, fdhB, fdhA, mtd)) +
  scale_x_continuous(breaks=c(0, 4, 9, 19, 28), limits=c(0,28)) +
  annotate("rect", xmin=0,xmax=28,ymin=c(-0.5, 0.5, -6),
           ymax=c(0.5, 2, -0.5), alpha= 0.1, fill=c("black", "blue", "red")) + 
  theme_classic()

# fdhB2 = MMP0139
## fdhA1 = MMP1298
# fdhB1 = MMP1297
# ehbN = MMP1153
## fruA = MMP1382
# frcA = MMP0820
## hmd  = MMP0127
# mtd = MMP0372
## vhcA = MMP0823



# vhuA = MMP1694


## hmd  = MMP0127
## fdhA2 = MMP0138 
## vhcA = MMP0823
## fruA = MMP1382
## frcA = MMP0820
## vhuA = MMP1694
## ehaO = MMP1462 = energy conserving hydrogenase A large subunit 


ehbN = bquote(paste(italic(ehb),N))
hmd = bquote(paste(italic(hmd)))
vhcA = bquote(paste(italic(vhc),A))
fruA = bquote(paste(italic(fru),A))
frcA = bquote(paste(italic(frc),A))
vhuA = bquote(paste(italic(vhu),A))
fdhA2 = bquote(paste(italic(fdh),A2))
fdhB2 = bquote(paste(italic(fdh),B2))
fdhA1 = bquote(paste(italic(fdh),A1))
fdhB1 = bquote(paste(italic(fdh),B1))
mtd = bquote(paste(italic(mtd)))
ehaO = bquote(paste(italic(eha),O))
yaxis = bquote(paste(log[2], ~Fold, ~Change, ~RPKM))

hmd_color = bquote("#2c6b67")
fdh_color = bquote("#7c4b73")
dH2ase_color = bquote("#88a0dc")
dmtd_color = bquote("#b0799a")
eha_color = bquote("#e9b109")

tndata_long %>% filter(count_type == "RPKM_log2FC") %>% filter(old_locus_tag == "MMP1298"|
                                                                 old_locus_tag == "MMP1382"|
                                                                 old_locus_tag == "MMP0820"|
                                                                 old_locus_tag == "MMP1694"|
                                                                 old_locus_tag == "MMP0823"|
                                                                 old_locus_tag == "MMP0127"|
                                                                 old_locus_tag == "MMP0372"|
                                                                 old_locus_tag == "MMP1462"|
                                                                 old_locus_tag == "MMP0138"|
                                                                 old_locus_tag == "MMP1153") %>% 
  ggplot(aes(x=sample,y=values,  color = old_locus_tag)) +
  geom_point(aes(shape = old_locus_tag, color = old_locus_tag, fill = old_locus_tag), size =3) + 
  geom_line(size=1) +
  labs(y=yaxis, x= "Generation")+
  scale_color_manual(breaks = c("MMP0127",
                                  "MMP0138",
                                  "MMP0823", 
                                  "MMP1382", 
                                  "MMP0820", 
                                  "MMP1694", 
                                  "MMP1462",
                                  "MMP1153",
                                  "MMP1298",
                                  "MMP0372"),
                       values = c(hmd_color,
                                  fdh_color,
                                  dH2ase_color,
                                  dH2ase_color,
                                  dH2ase_color,
                                  dH2ase_color,
                                  eha_color,
                                  eha_color,
                                  fdh_color,
                                  dmtd_color),
                       name = "Gene",
                       labels = c(hmd,
                                  fdhA2,
                                  vhcA, 
                                  fruA, 
                                  frcA, 
                                  vhuA, 
                                  ehaO,
                                  ehbN, 
                                  fdhA1, 
                                  mtd)) +
  scale_fill_manual(breaks = c("MMP0127",
                               "MMP0138",
                               "MMP0823", 
                               "MMP1382", 
                               "MMP0820", 
                               "MMP1694", 
                               "MMP1462",
                               "MMP1153",
                               "MMP1298",
                               "MMP0372"),
                    values = c(hmd_color,
                               fdh_color,
                               dH2ase_color,
                               dH2ase_color,
                               dH2ase_color,
                               dH2ase_color,
                               eha_color,
                               eha_color,
                               fdh_color,
                               dmtd_color),
                    name = "Gene",
                    labels = c(hmd,
                               fdhA2,
                               vhcA, 
                               fruA, 
                               frcA, 
                               vhuA, 
                               ehaO,
                               ehbN, 
                               fdhA1, 
                               mtd))+
  scale_shape_manual(values=c(0,
                              1,
                              2,
                              6,
                              25,
                              17,
                              23,
                              5,
                              16,
                              15), 
                     breaks = c("MMP0127",
                                "MMP0138",
                                "MMP0823", 
                                "MMP1382", 
                                "MMP0820", 
                                "MMP1694",
                                "MMP1462",
                                "MMP1153",
                                "MMP1298",
                                "MMP0372"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Gene", 
                     labels = c(hmd,
                                fdhA2,
                                vhcA, 
                                fruA, 
                                frcA, 
                                vhuA, 
                                ehaO,
                                ehbN, 
                                fdhA1, 
                                mtd)) +
  scale_x_continuous(limits=c(0,30)) + 
  theme_classic()

tndata_long %>% filter(count_type == "RPKM_log2FC") %>% filter(old_locus_tag == "MMP1298"|
                                                                 old_locus_tag == "MMP1382"|
                                                                 old_locus_tag == "MMP0820"|
                                                                 old_locus_tag == "MMP1694"|
                                                                 old_locus_tag == "MMP0823"|
                                                                 old_locus_tag == "MMP0127"|
                                                                 old_locus_tag == "MMP0372"|
                                                                 old_locus_tag == "MMP1462"|
                                                                 old_locus_tag == "MMP0138"|
                                                                 old_locus_tag == "MMP1153"|
                                                                 old_locus_tag == "MMP0713") %>% 
  ggplot(aes(x=sample,y=values,  color = old_locus_tag)) +
  geom_point(aes(shape = old_locus_tag, color = old_locus_tag, fill = old_locus_tag), size =3) + 
  geom_line(size=1) +
  scale_color_manual(breaks = c("MMP0127",
                                "MMP0138",
                                "MMP0823", 
                                "MMP1382", 
                                "MMP0820", 
                                "MMP1694", 
                                "MMP1462",
                                "MMP1153",
                                "MMP1298",
                                "MMP0372",
                                "MMP0713"),
                     values = c(hmd_color,
                                fdh_color,
                                dH2ase_color,
                                dH2ase_color,
                                dH2ase_color,
                                dH2ase_color,
                                eha_color,
                                eha_color,
                                fdh_color,
                                dmtd_color,
                                "red"),
                     name = "Gene",
                     labels = c(hmd,
                                fdhA2,
                                vhcA, 
                                fruA, 
                                frcA, 
                                vhuA, 
                                ehaO,
                                ehbN, 
                                fdhA1, 
                                mtd,
                                "ior")) +
  scale_fill_manual(breaks = c("MMP0127",
                               "MMP0138",
                               "MMP0823", 
                               "MMP1382", 
                               "MMP0820", 
                               "MMP1694", 
                               "MMP1462",
                               "MMP1153",
                               "MMP1298",
                               "MMP0372",
                               "MMP0713"),
                    values = c(hmd_color,
                               fdh_color,
                               dH2ase_color,
                               dH2ase_color,
                               dH2ase_color,
                               dH2ase_color,
                               eha_color,
                               eha_color,
                               fdh_color,
                               dmtd_color,
                               "red"),
                    name = "Gene",
                    labels = c(hmd,
                               fdhA2,
                               vhcA, 
                               fruA, 
                               frcA, 
                               vhuA, 
                               ehaO,
                               ehbN, 
                               fdhA1, 
                               mtd,
                               "ior"))+
  scale_shape_manual(values=c(0,
                              1,
                              2,
                              6,
                              25,
                              17,
                              23,
                              5,
                              16,
                              15,
                              2), 
                     breaks = c("MMP0127",
                                "MMP0138",
                                "MMP0823", 
                                "MMP1382", 
                                "MMP0820", 
                                "MMP1694",
                                "MMP1462",
                                "MMP1153",
                                "MMP1298",
                                "MMP0372",
                                "MMP0713"),
                     
                     ## NOTE THAT THESE NAMES AND LABELS MUST MATCH SCALE_COLOR TO PRINT ON SAME LEGEND
                     name = "Gene", 
                     labels = c(hmd,
                                fdhA2,
                                vhcA, 
                                fruA, 
                                frcA, 
                                vhuA, 
                                ehaO,
                                ehbN, 
                                fdhA1, 
                                mtd,
                                "ior")) +
  scale_x_continuous(limits=c(0,30)) + 
  theme_classic()
