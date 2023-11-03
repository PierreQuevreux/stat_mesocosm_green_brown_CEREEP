### PACKAGES ### ----
sessionInfo()
# R version 4.3.1
# attached packages
# rstudioapi_0.15.0  lubridate_1.9.3    vegan_2.6-4        permute_0.9-7      factoextra_1.0.7   FactoMineR_2.9    
# multcompView_0.1-9 lsmeans_2.30-0     pscl_1.5.5.1       blmeco_1.4         multcomp_1.4-25    TH.data_1.1-2     
# survival_3.5-5     mvtnorm_1.2-3      emmeans_1.8.8      car_3.1-2          carData_3.0-5      lmerTest_3.1-3    
# MASS_7.3-60        nlme_3.1-162       lme4_1.1-34        Matrix_1.6-1.1     tibble_3.2.1       dplyr_1.1.3       
# Rmisc_1.5.1        plyr_1.8.9         lattice_0.21-8     reshape2_1.4.4     Unicode_15.0.0-1   RColorBrewer_1.1-3
# viridis_0.6.4      viridisLite_0.4.2  scales_1.2.1       gridExtra_2.3      cowplot_1.1.1      ggplot2_3.4.3     
# loaded via a namespace (and not attached):
# tidyselect_1.2.0     fastmap_1.1.1        digest_0.6.33        estimability_1.4.1   timechange_0.2.0    
# lifecycle_1.0.3      cluster_2.1.4        magrittr_2.0.3       compiler_4.3.1       rlang_1.1.1         
# tools_4.3.1          utf8_1.2.3           htmlwidgets_1.6.2    scatterplot3d_0.3-44 abind_1.4-5         
# withr_2.5.1          numDeriv_2016.8-1.1  grid_4.3.1           fansi_1.0.4          xtable_1.8-4        
# colorspace_2.1-0     flashClust_1.01-2    cli_3.6.1            generics_0.1.3       minqa_1.2.6         
# stringr_1.5.0        splines_4.3.1        parallel_4.3.1       vctrs_0.6.3          boot_1.3-28.1       
# sandwich_3.0-2       ggrepel_0.9.3        glue_1.6.2           nloptr_2.0.3         codetools_0.2-19    
# DT_0.30              stringi_1.7.12       gtable_0.3.4         munsell_0.5.0        pillar_1.9.0        
# htmltools_0.5.6.1    R6_2.5.1             leaps_3.1            arm_1.13-1           Rcpp_1.0.11         
# coda_0.19-4          mgcv_1.8-42          zoo_1.8-12           pkgconfig_2.0.3     
# Graph packages
library(ggplot2)
library(cowplot)
library(gridExtra)
library(scales)
library(viridis)
library(RColorBrewer)
library(Unicode)
# Data management
library(reshape2)
library(Rmisc)
library(dplyr)
library(tibble)
# Stat packages
library(lme4) # mixed models
library(nlme) # random effects
library(MASS) # for quasi-Poisson mixed model
library(lmerTest)
library(car) # Anova function (type II anova)
library(emmeans) # for post hoc Tukey tests with lmer
library(multcomp) # to dispay the results of the Tukey test
library(blmeco) # to test the dispersion of Poisson distribution
library(Rmisc)
library(pscl)
library(lsmeans) # for post hoc Tukey tests with lmer
library(lmerTest)
library(multcompView)
library(FactoMineR) # PCA analysis
library(factoextra) # to represent eigen values of PCA
library(vegan)
# Miscellaneous packages
library(lubridate)
library(rstudioapi) # to set the working directory

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
path_figures="Figures/"
path_data="Data/"

### PLOT OPTIONS ### ----
### colour scales ###
colour_OM<-scale_fill_manual(values=c("lightgray","dimgray"))
colour_light<-scale_colour_manual(values=c("black","dimgray"))
colour_chloro<-scale_fill_manual(values=c("green1","green4"))
colour_zoo<-scale_fill_manual(values=c("darkorange1","darkorange4"))
colour_fish<-scale_fill_manual(values=c("red","darkred"))
colour_sediment<-scale_fill_manual(values=c("lightgray","dimgray"))
colour_heterotrophs<-scale_fill_manual(values=c("darkgoldenrod1","darkgoldenrod4"))
colour_bacteria<-scale_fill_manual(values=c("hotpink1","hotpink4"))

### themes ###
theme<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour='grey'),
        panel.grid.major.y = element_line(colour='grey'),
        text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        axis.title.x = element_blank(),
        legend.key=element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))

### axes ###
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

### labels ###
label_algae<-expression("Concentration (cell mL"^-1*")")
label_chloro<-expression("Chlorophyll "*italic(a)*" (\u03BCg L"^{-1}*")")
label_zoo<-"Individuals per 24 L"
label_OD<-"Optical density (OD)"

coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

### OUTPUT FORMATTING ### ----
get_stat_table<-function(stat){
  stat_out<-stat
  for(i in 1:ncol(stat_out)){
    stat_out[,i]<-as.character(stat_out[,i])
  }
  for(i in c(1,3)){
    stat_out[stat[,i]>=10,i]<-round(stat[stat[,i]>=10,i],2)
    stat_out[stat[,i]<10,i]<-round(stat[stat[,i]<10,i],3)
  }
  stat_out[stat_out[,3]=="0",3]<-"<.001"
  stat_out<-t(stat_out)
  stat_out<-stat_out[c(3,1,2),]
  return(stat_out)
}
  
### LOAD DATA ### ----
# treatments (table summarising the characteristics of neach mesocosm) #
treatment<-read.table(paste0(path_data,"treatment.csv"),sep=";",header=TRUE)
treatment[1:ncol(treatment)]<-lapply(treatment[1:ncol(treatment)], as.factor) 
levels(treatment$fish)<-c("Fish -","Fish +")
levels(treatment$OM)<-c("OM -","OM +")
levels(treatment$light)<-c("Light -","Light +")
treatment$treatment<-paste(treatment$fish,treatment$light,treatment$OM,sep=" ")

# BBE data (spectroscopy of phytoplankton pigments measured by the BBE probe) #
BBE<-read.table(paste0(path_data,"BBE.txt"),sep=';',header=T)
BBE$mesocosm<-as.factor(BBE$mesocosm)
BBE<-merge(BBE,treatment,by="mesocosm")

# Dry zooplankton biomass #
zoodry<-read.table(paste0(path_data,"zoo_dry.csv"),sep=';',header=T)
zoodry$mesocosm<-as.factor(zoodry$mesocosm)
zoodry$concentration<-(zoodry$after-zoodry$before)/zoodry$volume*1000 # conversion into liters
zoodry$date<-"01/10/15"
zoodry<-merge(zoodry,treatment,by="mesocosm")

# fish growth #
fish<-read.table(paste0(path_data,"fish.csv"),sep=";",header=TRUE)
fish$mesocosm<-as.factor(fish$mesocosm)
fish$growth<-(fish$final_mass - fish$initial_mass)/fish$duration
fish$relative_growth<-(fish$final_mass - fish$initial_mass)/fish$duration/fish$initial_mass
fish<-merge(fish,treatment,by="mesocosm")

# Sediments (sediments collected at the end of the experiment) #
sediment<-read.table(paste0(path_data,"sediment.csv"),sep=';',header=T)
sediment$mesocosm<-as.factor(sediment$mesocosm)
sediment$mass<-sediment$poids_total_sec-sediment$poids_pot # difference between the total dry weight and the weight of the jar
sediment$date<-"09/11/15"
sediment<-merge(sediment,treatment,by="mesocosm")

# Cytometry concentration (data as cell concentration) #
cyto<-read.table(paste0(path_data,"cyto.txt"),sep=';',header=T)
cyto$mesocosm<-as.factor(cyto$mesocosm)
cyto$date<-as.factor(cyto$date)
cyto$date=factor(cyto$date,levels(cyto$date)[c(2,3,1,4)])
cyto<-merge(cyto,treatment,by="mesocosm")

# Cytometry count (data as cell count) #
cyto_count<-read.table(paste0(path_data,"cyto_count.txt"),sep=';',header=T)
cyto_count$mesocosm<-as.factor(cyto_count$mesocosm)
cyto_count$sample_date<-as.factor(cyto_count$sample_date)
cyto_count$sample_date=factor(cyto_count$sample_date,levels(cyto_count$sample_date)[c(2,3,1,4)])
cyto_count<-merge(cyto_count,treatment,by="mesocosm")

# DOC (Dissolved organic carbon measured by the Shimadzu) #
DOC<-read.table(paste0(path_data,"DOC.txt"),sep=';',header=T)
DOC$mesocosm<-as.factor(DOC$mesocosm)
DOC<-DOC[DOC$DOC<25,] # removal of aberrant points
DOC<-merge(DOC,treatment,by="mesocosm")
DOC$day<-ymd(DOC$date) # turn the dates into the appropriate format
Sys.setenv(LANGUAGE="en")
Sys.setlocale("LC_TIME", "English")

# Phytoplankton taxa abundance (for PCA analysis, microscope identification data) #
phyto<-read.table(paste0(path_data,"algues_species.csv"),sep=';',header=T)
### per species #
phyto_1<-t(as.matrix(phyto[,c(2,6:ncol(phyto))]))
temp<-phyto_1[1,1:(ncol(phyto_1)-1)]
phyto_1<-as.data.frame(phyto_1[2:nrow(phyto_1),1:(ncol(phyto_1)-1)])
names(phyto_1)<-temp
phyto_1<-phyto_1 %>% mutate_at(1:ncol(phyto_1), as.numeric)
phyto_1$mesocosm<-c(1:36)
phyto_1$mesocosm<-as.factor(phyto_1$mesocosm)
phyto_1<-merge(treatment,phyto_1,by="mesocosm")
### aggregated per family #
phyto_2<-melt(phyto, id.vars = c("family_old","species","family_1","family_2","mixotrophic"),
              variable.name = "mesocosm", 
              value.name = "concentration")
phyto_2<-aggregate(concentration~family_1+mesocosm, data=phyto_2, sum)
phyto_2<-dcast(data=phyto_2,mesocosm~family_1,value.var="concentration")
phyto_2<-phyto_2[,c(2:ncol(phyto_2),1)]
phyto_2<-phyto_2 %>% mutate_at(1:ncol(phyto_2), as.numeric)
phyto_2$mesocosm<-c(1:36)
phyto_2$mesocosm<-as.factor(phyto_2$mesocosm)
phyto_2<-merge(treatment,phyto_2,by="mesocosm")

# Shannon-Weaver diversity
diversity<-read.table(paste0(path_data,"algues_diversity.csv"),sep=';',header=T)
diversity<-merge(diversity,treatment,by="mesocosm")

# Zooplankton taxa abundance (for PCA analysis, microscope identification data) #
zoo1<-read.table(paste0(path_data,"zoo_2015_08_19.txt"),sep=';',header=T)
zoo2<-read.table(paste0(path_data,"zoo_2015_09_25.csv"),sep=';',header=T)
zoo1$mesocosm<-as.factor(zoo1$mesocosm)
zoo2$mesocosm<-as.factor(zoo2$mesocosm)
zoo1<-subset(zoo1,select=-c(conochilus,ceratium,hour,counting_date))
zoo1$brachionus<-0
zoo1$anuraeopsis<-0
zoo1$hexarthra<-0
zoo1$bosmine<-0
zoo1$large_lecane<-0
zoo2<-subset(zoo2,select=-c(filinia,oligochete))
zoo2$ascomorpha<-0
zoo1$cladocerae<-zoo1$chydoridae
zoo2$cladocerae<-zoo2$chydoridae+zoo2$bosmine
zoo1$copepode<-zoo1$nauplius+zoo1$cyclopide
zoo2$copepode<-zoo2$nauplius+zoo2$cyclopide
zoo1$rotifer<-zoo1$lecane+zoo1$lepadella+zoo1$bdelloide+zoo1$polyarthra+zoo1$ascomorpha+zoo1$keratella
zoo2$rotifer<-zoo2$lecane+zoo2$lepadella+zoo2$bdelloide+zoo2$polyarthra+
  zoo2$brachionus+zoo2$anuraeopsis+zoo2$hexarthra+zoo2$large_lecane
#zoo1<-zoo1[,c("mesocosm","date","chaoborus","cladocerae","copepode","rotifer")]
#zoo2<-zoo2[,c("mesocosm","date","chaoborus","cladocerae","copepode","rotifer")]
zoo2<-zoo2[,names(zoo1)]
zoo<-rbind(zoo1,zoo2)
rm(zoo1)
rm(zoo2)
zoo<-merge(zoo,treatment,by="mesocosm")
zoo$date<-as.factor(zoo$date)
levels(zoo$date)<-c("19/08/2015","25/09/2015")

# Ecoplates (functional diversity of bacteria according to substrat degradation) #
eco<-read.table(paste0(path_data,"ecoplate.txt"),sep=';',header=T)
eco$mesocosm<-as.factor(eco$mesocosm)
eco$day<-as.factor(eco$day)
eco<-dcast(eco,substrat+mesocosm+date+compartment+series+name+family~day,value.var="absorbance")
eco$absorbance<-eco$'7'-eco$'1' # absorbance 7 days later
eco<-eco[,c(1:7,14)]
eco<-merge(eco,treatment,by="mesocosm")
eco<-eco[eco$family!="blank",]
eco$compartment<-factor(eco$compartment,levels=c("pelagic","benthic"))

# CNP (stoichiometry of the seston) #
CNP<-read.table(paste0(path_data,"CNP.csv"),sep=';',header=T)
CNP$mesocosm<-as.factor(CNP$mesocosm)
CNP<-merge(CNP,treatment,by="mesocosm")

# Multipar (multiparameter probe: chlorophyll, pH, oxygen...) #
multi<-read.table(paste0(path_data,"multipar.txt"),sep=";",header=TRUE)
multi$depth<-as.factor(multi$depth)
multi<-merge(multi,treatment,by="mesocosm")

# Dissolved nutrients
nutri<-read.table(paste0(path_data,"nutrients.csv"),sep=';',header=T)
nutri$mesocosm<-as.factor(nutri$mesocosm)
nutri$date<-as.factor(nutri$date)
nutri$date=factor(nutri$date,levels(nutri$date)[c(3,1,2)])
levels(nutri$date)<-c("07/23/2015","10/01/2015","10/21/2015")
nutri$N=nutri$NO2*14/(14+2*16)+ # calculate N mass (removes O and H)
        nutri$NO3*14/(14+3*16)+
        nutri$NH4*14/(14+4)
nutri$P=nutri$PO4*31/(31+4*16) # calculate P mass
nutri$NP<-nutri$N/nutri$P
nutri<-nutri[nutri$remark=="ok",] # remove aberrant data points (measures are not very relyable)
nutri<-nutri[nutri$P<100,]
data_N<-summarySE(nutri, measurevar="N", groupvars=c("date","mesocosm"),na.rm=TRUE)
data_N<-data_N[,c(1,2,4)]
data_P<-summarySE(nutri, measurevar="P", groupvars=c("date","mesocosm"),na.rm=TRUE)
data_P<-data_P[,c(1,2,4)]
data_NP<-summarySE(nutri, measurevar="NP", groupvars=c("date","mesocosm"),na.rm=TRUE)
data_NP<-data_NP[,c(1,2,4)]
data<-expand.grid(mesocosm=1:36,
                  date=levels(nutri$date))
data<-merge(data,treatment,by="mesocosm",all=TRUE)
nutri<-merge(data_N,data,by=c("date","mesocosm"),all=TRUE)
nutri<-merge(nutri,data_P,by=c("date","mesocosm"),all=TRUE)
nutri<-merge(nutri,data_NP,by=c("date","mesocosm"),all=TRUE)
rm(data_N,data_P,data_NP,data)
nutri<-nutri[nutri$date=="07/23/2015" | nutri$date=="10/21/2015",]
nutri$date<-droplevels(nutri$date)
levels(nutri$date)<-c("Jul","Oct")
#nutri$date<-mdy(nutri$date) # turn the dates into the appropriate format

# seston
seston<-read.table(paste0(path_data,"seston_dry.csv"),sep=';',header=T)
sediment$mesocosm<-as.factor(sediment$mesocosm)
seston<-merge(seston,treatment,by="mesocosm")
seston$concentration<-(seston$after - seston$before)/seston$volume*1000 # conversion into liters
seston$date<-"01/10/15"

### RESULT STORAGE STRUCTURES ### ----
stat_list<-list()
graph_list<-list()

############# ----
# MAIN TEXT # ----
############# ----
### GENERAL FOOD WEB DESCRIPTION # ----
# Chlorophyll BBE # ----
data<-BBE[BBE$taxon_BBE=="green_algae",]

model<-lmer(data=data, log(BBE)~fish*light*OM+(1|mesocosm)+(1|date), na.action = "na.fail", REML=FALSE)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$chloro=get_stat_table(stat)

# Tukey test
model<-lmer(data=BBE[BBE$taxon_BBE=="green_algae",], log(BBE)~fish*light+(1|mesocosm)+(1|date), na.action = "na.fail", REML=FALSE)
summary(model) # to get the quantitative effect of significant treatments
tukey<-emmeans(model, specs = pairwise~fish:light)
tukey<-cld(tukey$emmeans,
           alpha=0.05,
           Letters=letters)
tukey$y<-c(50,75,75,75)

# figure
data<-summarySE(data, measurevar="BBE", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)
p1<-ggplot(data=data)+
  geom_boxplot(aes(light,BBE,fill=OM))+
  geom_text(data=tukey,aes(light,y,label=.group),size=7,fontface = "bold")+
  facet_wrap(~fish)+
  colour_chloro+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(label_chloro)

# displays significant effects on the graph
label<-data.frame(label="fish ***\nlight ***\nfish:light ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.5,0.1,0.45,0.2,hjust=0)

graph_list$chloro=graph

# Zooplankton dry mass # ----
model<-lm(data=zoodry,log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$zoo=get_stat_table(stat)

p1<-ggplot(data=zoodry)+
  geom_boxplot(aes(light,concentration,fill=OM))+
  facet_wrap(~fish)+
  colour_zoo+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Zooplankton (mg "*L^{-1}*")"))

graph_list$zoo=p1

# Fish growth # ----
model<-lm(data=fish,growth~light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:3,c("F value","Df","Pr(>F)")]
stat_list$fish=get_stat_table(stat)

# Tukey test
model<-lm(data=fish, growth~light)
summary(model) # to get the quantitative effect of significant treatments
tukey<-emmeans(model, specs = ~light)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=c(1,2)
tukey$y=c(0.012,0.02)

p1<-ggplot(data=fish)+
  geom_boxplot(aes(light,growth,fill=OM),position="dodge")+
  geom_text(data=tukey,aes(x,y+0.002,label=.group),size=7,fontface = "bold")+
  colour_fish+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Fish growth (g "*day^{-1}*")"))

# displays significant effects on the graph
label<-data.frame(label="light **",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.7,0.05,0.3,0.2,hjust=0)

graph_list$fish=graph

# Sediments # ---------------------------------
model<-lm(data=sediment,log(mass)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$sediment=get_stat_table(stat)

# Tukey test
model<-lm(data=sediment,log(mass)~fish)
summary(model) # to get the quantitative effect of significant treatments
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(0.3,0.5)

p1<-ggplot(data=sediment)+
  geom_boxplot(aes(light,mass,fill=OM))+
  geom_text(data=tukey,aes(x,y+0.02,label=.group),size=7,fontface = "bold")+
  geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=0.01)+
  facet_wrap(~fish)+
  colour_sediment+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Sediments (mg)"))

# displays significant effects on the graph
label<-data.frame(label="fish ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.7,0.03,0.3,0.2,hjust=0)

graph_list$sediment=graph

# Bacteria # ----
model<-glmer(bacteria~fish*light*OM+(1|mesocosm)+(1|sample_date), offset=log(volume_bacteria),family=poisson,data=cyto_count)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$bacteria=get_stat_table(stat)

# Tukey test
model<-glmer(bacteria~fish+(1|mesocosm)+(1|sample_date), offset=log(volume_bacteria),family=poisson,data=cyto_count)
summary(model) # to get the quantitative effect of significant treatments
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(2900,3400)

# figure
data<-cyto[cyto$taxon_cyto=="bacteria",]
data<-summarySE(data, measurevar="cyto", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)

p1<-ggplot(data=data)+
  geom_boxplot(aes(light,cyto,fill=OM))+
  geom_text(data=tukey,aes(x,y+200,label=.group),size=7,fontface = "bold")+
  geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=200)+
  facet_wrap(~fish)+
  colour_bacteria+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Bacteria (cell µL"^{-1}*")"))

# displays significant effects on the graph
label<-data.frame(label="fish **",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.75,0.03,0.25,0.2,hjust=0)

graph_list$bacteria=graph

# Heterotroph # ----
model<-glmer(hetero~fish*light*OM+(1|mesocosm)+(1|sample_date), offset=log(volume_hetero),family=poisson,data=cyto_count)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$hetero=get_stat_table(stat)

# Tukey test
model<-glmer(hetero~fish+(1|mesocosm)+(1|sample_date), offset=log(volume_hetero),family=poisson,data=cyto_count)
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(80,60)

# figure
databis<-cyto[cyto$taxon_cyto=="heterotroph",]
databis<-summarySE(databis, measurevar="cyto", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)

p1<-ggplot(data=databis)+
  geom_boxplot(aes(light,cyto,fill=OM))+
  geom_text(data=tukey,aes(x,y+5,label=.group),size=7,fontface = "bold")+
  geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=5)+
  facet_wrap(~fish)+
  colour_heterotrophs+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Heterotroph (cell µL"^{-1}*")"))+
  ylim(15,100)

# displays significant effects on the graph
label<-data.frame(label="fish *",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.75,0.78,0.25,0.2,hjust=0)

graph_list$hetero=graph

# DOC # ----
# continiuous time
model<-lmer(log(DOC)~time*fish*light*OM+(1|mesocosm), na.action = "na.fail",REML=FALSE,data=DOC)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model) # to get the quantitative effect of significant treatments
stat<-Anova(model,type=2,test.statistic="Chisq")
stat<-get_stat_table(stat)
stat<-stat[,colnames(stat)%in%c("time:fish","time:light","time:OM","time:fish:light","time:fish:OM","time:light:OM","time:fish:light:OM")]
colnames(stat)<-c("fish","light","OM","fish:light","fish:OM","light:OM","fish:light:OM")
stat_list$DOC=stat

# discrete time
model<-lmer(data=DOC, log(DOC)~fish*light*OM+(1|mesocosm)+(1|date), na.action = "na.fail", REML=FALSE)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
#stat_list<-c(stat_list,list(chloro=get_stat_table(stat)))

# figure
p1<-ggplot(DOC[DOC$DOC<25,])+ # remove the 5 aberrant extra points
  geom_point(aes(day,DOC),alpha=0.5)+
  geom_smooth(aes(day,DOC,color=light,linetype=OM),linewidth=2,method='lm')+
  facet_wrap(~fish)+
  scale_colour_manual(values=c("darkblue","deepskyblue"))+
  scale_linetype_manual(values=c("22","solid"))+
  theme+theme(legend.key.width= unit(1, 'cm'),
              legend.position=c(0.07,0.8),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  guides(colour=guide_legend(override.aes=list(fill=NA)),
         linetype=guide_legend(override.aes=list(fill=NA))) + 
  xlab("")+
  ylab(expression("DOC (mg "*L^{-1}*")"))

# displays significant effects on the graph
label<-data.frame(label="fish ***\nlight ***\nOM ***\nfish:light ***\nfish:OM *\nlight:OM *\nfish:light:OM **",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.4,0.4,0.35,0.6,hjust=0)

graph_list$DOC=graph

# Figure food web # ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 3)) +
  draw_plot(graph_list$chloro, 0, 1, 1, 1)+
  draw_plot(graph_list$zoo, 0, 2, 1, 1)+
  draw_plot(graph_list$fish, 1, 2, 1, 1)+
  draw_plot(graph_list$hetero, 2, 2, 1, 1)+
  draw_plot(graph_list$bacteria, 2, 1, 1, 1)+
  draw_plot(graph_list$sediment, 2, 0, 1, 1)+
  draw_plot(graph_list$DOC, 0, 0, 2, 1)+
  draw_plot_label(c("A","B","C","D","E","F","G"),
                  c(0,0,1,2,2,2,0),
                  c(2,3,3,3,2,1,1), size = 30)
ggsave(paste(path_figures,"figure_food_web.pdf",sep=""), graph, width = 18, height = 16, device=cairo_pdf)

# Export stat table # ----
stat_table<-NULL
names_data<-c("Chlorophyll","Zooplankton","Fish growth","Sediments","Bacteria","Heterotrophs","DOC")
for(i in 1:length(stat_list)){
  data<-as.data.frame(stat_list[[i]])
  data<-rownames_to_column(data,"row_names")
  data<-cbind(Variable=c("",names_data[i],""),data)
  stat_table<-bind_rows(stat_table,data)
}
write.table(stat_table,paste0(path_figures,"stat_foodweb.csv"),row.names=FALSE)
### FUNCTIONAL RESPONSE OF COMPARTMENTS # ----
# PERMANOVA phytoplankton # ----
### per species #
data<-phyto_1[,c((ncol(treatment)+1):ncol(phyto_1))]
adonis2(data ~ fish*light*OM, data=phyto_1, permutations = 10000, method="bray")

### aggregated per family #
data<-phyto_2[,c((ncol(treatment)+1):ncol(phyto_2))]
data$none<-NULL
adonis2(data ~ fish*light*OM, data=phyto_2, permutations = 10000, method="bray")

# PCA phytoplankton # ----
### per species #
data<-phyto_1[,c((ncol(treatment)+1):ncol(phyto_1))]
data_pca<-PCA(data, scale.unit = TRUE, ncp = 5, graph = TRUE)
data_pca$eig

fviz_eig(data_pca, addlabels = TRUE) +
  theme+theme(plot.title=element_blank())

fviz_pca_biplot(data_pca,
                label = "var") +
  geom_point(aes(shape=factor(phyto_1$light),
                 colour=factor(phyto_1$fish)),size=3) +
  scale_color_manual(values=c("black","red"))+
  scale_size_manual(values=c(2,4))+
  theme+theme(axis.title.x = element_text())+
  ggtitle("PCA phytoplakton per species")

### aggregated per family #
data<-phyto_2[,c((ncol(treatment)+1):ncol(phyto_2))]
data$none<-NULL
data_pca<-PCA(data, scale.unit = TRUE, ncp = 5, graph = TRUE)
data_pca$eig

fviz_eig(data_pca, addlabels = TRUE) +
  theme+theme(plot.title=element_blank())

p1<-fviz_pca_biplot(data_pca,
                label = "var",
                col.var ="black") +
  geom_point(aes(shape=factor(phyto_2$light),
                 colour=factor(phyto_2$fish)),size=3) +
  scale_color_manual(values=c("cadetblue","orangered"))+
  scale_size_manual(values=c(2,4))+
  theme+theme(axis.title.x = element_text(),
              legend.box="horizontal",
              legend.position=c(0.75,0.92),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ggtitle("PCA phytoplakton")

label<-data.frame(label="fish *\nfish:light **",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.5,0.15,0.5,0.2,hjust=0)

graph_list$PCA_phyto=graph

# PERMANOVA zooplankton # ----
data<-zoo[,c(3:18)]
adonis2(data ~ fish*light*OM, data=zoo, permutations = 1000, method="bray", strata=zoo$date)

data_rel<-decostand(data, method = "total") # convert into relative abundances
data_dist<-as.matrix(vegdist(data_rel, method = "bray")) # distance matrix
NMDS<-metaMDS(data_dist,
              distance = "bray",
              k = 2,
              maxit = 999, 
              trymax = 500,
              wascores = TRUE)
stressplot(NMDS)
data_plot<-scores(NMDS) %>% as_tibble()
data_plot$mesocosm<-zoo$mesocosm
data_plot$date<-zoo$date
data_plot<-merge(data_plot,treatment,by="mesocosm")

ggplot(data=data_plot,aes(NMDS1,NMDS2,colour=fish))+
  geom_point(size=3)+
  stat_ellipse()+
  facet_wrap(~date, scale="free")+
  scale_colour_manual(values=c("blue","red"))+
  theme

# PCA zooplankton # ----
zoo$FL<-as.factor(paste(zoo$fish,zoo$light))
zoo$treatment<-as.factor(zoo$treatment)
#data<-zoo[,c("cladocerae","copepode","rotifer")]
data<-zoo[,c(3:11,13:14)]
data<-zoo[,c("cladocerae","copepode","chaoborus","rotifer")]
data_pca<-PCA(data, scale.unit = TRUE, ncp = 5, graph = TRUE)
data_pca$eig

fviz_eig(data_pca, addlabels = TRUE) +
  theme+theme(plot.title=element_blank())

p1<-fviz_pca_biplot(data_pca,
                label = "var",
                col.var ="black") +
  geom_point(aes(colour=factor(zoo$fish)),size=3) +
  scale_color_manual(values=c("cadetblue","orangered"))+
  scale_size_manual(values=c(2,4))+
  theme+theme(axis.title.x = element_text(),
              legend.position=c(0.15,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ggtitle("PCA zooplakton")

label<-data.frame(label="fish ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.02,0.6,0.3,0.2,hjust=0)

graph_list$PCA_zoo=graph

# Functional diversity of bacteria (radar chart) # ----
data<-eco
data<-summarySE(data, measurevar="absorbance", groupvars=c("family","treatment","compartment"),na.rm=TRUE)
data<-data[,c("absorbance","treatment","compartment","family")]
levels(data$compartment)<-c("Pelagic","Benthic")
data$family<-as.factor(data$family)
levels(data$family)<-c("amines","amino\nacids","carbohydrate","carboxylic\nacids","phenolic\nacids","polymer")

p1<-ggplot(data=data,aes(family, absorbance))+
  geom_polygon(aes(group=treatment, colour=treatment), fill = NA, size = 2, show.legend = FALSE)+
  geom_line(aes(group=treatment, colour=treatment), linewidth = 2)+
  facet_wrap(~compartment)+
  coord_radar(start=pi/6)+
  scale_colour_manual(values = c("cyan2","cadetblue","dodgerblue4","blue4",
                                 "goldenrod2","darkorange","orangered","firebrick3"))+
  theme+theme(legend.title=element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(vjust=0.5,size=15),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_line(colour="black",linetype = "dashed"),
              panel.background = element_blank(),
              strip.background = element_rect(colour=NA,fill=NA))+
  ylab("Optical density (OD)")

graph_list$ecoplate=p1

# Figure functional diversity # ----
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(graph_list$PCA_phyto, 0, 1, 1, 1)+
  draw_plot(graph_list$PCA_zoo, 1, 1, 1, 1)+
  draw_plot(graph_list$ecoplate, 0, 0, 2, 1)+
  draw_plot_label(c("A","B","C"),
                  c(0,1,0),
                  c(2,2,1), size = 30)
ggsave(paste0(path_figures,"figure_fundiv.pdf"), graph, width = 12, height = 10, device=cairo_pdf)

### PHYSIOLOGICAL RESPONSE OF PHYTOPLANKTON # ----
# Chlorophyll concentration in seston (~cellular concentration of chlorophyll) # ----
data<-merge(BBE[BBE$taxon_BBE=="green_algae",],seston[,c("mesocosm","concentration","date")],by=c("mesocosm","date"))
data$ratio<-data$BBE/data$concentration

model<-lm(data=data,ratio~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$chloro_seston=get_stat_table(stat)

p1<-ggplot(data=data)+
  geom_boxplot(aes(light,ratio,fill=OM))+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position=c(0.9,0.1),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Chlorophyll "*italic(a)*" in seston (\u03BCg " *mg^{-1}*")"))

graph_list$chloro_seston=p1

# Mixotrophs # ----
data<-melt(phyto[phyto$mixotrophic=="yes",], id.vars = c("family_old","species","family_1","family_2","mixotrophic"),
           variable.name = "mesocosm", 
           value.name = "concentration")
data$mesocosm<-as.factor(data$mesocosm)
levels(data$mesocosm)<-c(1:36)
data<-merge(data,treatment,by="mesocosm")
data<-aggregate(data=data, concentration~mesocosm+OM+light+fish,FUN=sum)

model<-lm(data=data,concentration~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$mixo=get_stat_table(stat)

p1<-ggplot(data=data)+
  geom_boxplot(aes(light,log(concentration),fill=OM))+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position=c(0.1,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Mixotrophs (Log cell "*mL^{-1}*")"))

graph_list$mixo=p1

# C:N of seston # ----
model<-lm(data=CNP,C.N~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$CN=get_stat_table(stat)

p1<-ggplot(data=CNP)+
  geom_boxplot(aes(light,C.N,fill=OM))+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position=c(0.1,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Seston C:N ratio"))

graph_list$CN=p1

# C:P of seston # ----
model<-lm(data=CNP,C.P~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:3,c("F value","Df","Pr(>F)")]
stat_list$CP=get_stat_table(stat)

# Tukey test
model<-lm(data=CNP, C.P~fish)
summary(model) # to get the quantitative effect of significant treatments
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(900,1120)

# figure
p1<-ggplot(data=CNP)+
      geom_boxplot(aes(light,C.P,fill=OM))+
      geom_text(data=tukey,aes(x,y+50,label=.group),size=7,fontface = "bold")+
      geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=30)+
      facet_wrap(~fish)+
      colour_OM+
      theme+theme(legend.position=c(0.1,0.9),
                  legend.background=element_blank(),
                  legend.box.background=element_rect(colour="black",fill="white"))+
      ylab(expression("Seston C:P ratio"))

# displays significant effects on the graph
label<-data.frame(label="fish ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.75,0.03,0.25,0.2,hjust=0)

graph_list$CP=graph

# Figure physiological response of phytoplankton # ----
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(graph_list$chloro_seston, 0, 1, 1, 1)+
  draw_plot(graph_list$mixo, 1, 1, 1, 1)+
  draw_plot(graph_list$CN, 0, 0, 1, 1)+
  draw_plot(graph_list$CP, 1, 0, 1, 1)+
  draw_plot_label(c("A","B","C","D"),
                  c(0,1,0,1),
                  c(2,2,1,1), size = 30)
ggsave(paste(path_figures,"figure_physio.pdf",sep=""), graph, width = 14, height = 10, device=cairo_pdf)

########################## ----
# SUPPORTING INFORMATION # ----
########################## ----

### Multipar # ----
# Statistic table and figures # ----
stat_list<-list()

# Chlorophyll # ----
model<-lmer(log(chlorophyll)~fish*light*OM+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$chloro_multipar=get_stat_table(stat)

# Tukey test
model<-lmer(log(chlorophyll)~fish*light+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
tukey<-emmeans(model, specs = pairwise~fish:light)
tukey<-cld(tukey$emmeans,
           alpha=0.05,
           Letters=letters)
tukey$y<-c(35,55,55,55)

# time series
data<-summarySE(BBE[BBE$taxon_BBE=="green_algae",], measurevar="BBE", groupvars=c("date","fish","light"),na.rm=TRUE)
data$date<-dmy(data$date) # turn the dates into the appropriate format
Sys.setenv(LANGUAGE="en")

p1<-ggplot(data=data)+
      geom_line(aes(date,BBE,colour=light,linetype=fish),size=2)+
      geom_errorbar(aes(date,ymin=BBE-sd, ymax=BBE+sd,colour=light), width=5)+
      scale_colour_manual(values=c("green4","green1"))+
      scale_linetype_manual(values=c("22","solid"))+
      theme+theme(legend.key.width= unit(1, 'cm'),
                  legend.position=c(0.9,0.8),
                  legend.background=element_blank(),
                  legend.box.background=element_rect(colour="black",fill="white"))+
      ylab(label_chloro)

graph_list$chloro_TS=p1

# chlorophyll multipar
data<-summarySE(multi, measurevar="chlorophyll", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)
p1<-ggplot(data=data)+
  geom_boxplot(aes(light,chlorophyll,fill=OM))+
  geom_text(data=tukey,aes(light,y,label=.group),size=7,fontface = "bold")+
  facet_wrap(~fish)+
  colour_chloro+
  theme+theme(legend.position=c(0.1,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(label_chloro)

# displays significant effects on the graph
label<-data.frame(label="fish ***\nlight ***\nfish:light *\ndepth ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.65,0.15,0.35,0.25,hjust=0)

graph_list$chloro_multi=graph

# chloro stratification
data<-summarySE(multi, measurevar="chlorophyll", groupvars=c("light","fish","mesocosm","depth"),na.rm=TRUE)
p1<-ggplot(data=data)+
  geom_boxplot(aes(depth,chlorophyll,fill=light))+
  facet_wrap(~fish)+
  scale_fill_manual(values=c("green4","green1"))+
  theme+theme(legend.position=c(0.12,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"),
              axis.title.x=element_text())+
  xlab("Depth (m)")+
  ylab(label_chloro)

# displays significant effects on the graph
label<-data.frame(label="fish ***\nlight ***\nfish:light *\ndepth ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.65,0.15,0.35,0.25,hjust=0)

graph_list$chloro_depth=graph

# final graph
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(graph_list$chloro_TS, 0, 1, 2, 1)+
  draw_plot(graph_list$chloro_multi, 0, 0, 1, 1)+
  draw_plot(graph_list$chloro_depth, 1, 0, 1, 1)+
  draw_plot_label(c("A","B","C"),
                  c(0,0,1),
                  c(2,1,1), size = 30)
ggsave(paste(path_figures,"supp_chloro.pdf",sep=""), graph, width = 14, height = 12, device=cairo_pdf)

# Temperature # ----
model<-lmer(temperature~fish*light*OM+depth+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$temperature=get_stat_table(stat)

# figure
data<-summarySE(data=multi[multi$temperature>10,],measurevar="temperature",groupvars=c("depth","date"))
data<-na.omit(data)
data$date<-mdy(data$date) # turn the dates into the appropriate format
Sys.setenv(LANGUAGE="en")

p1<-ggplot(data=data)+
      geom_line(aes(date,temperature,colour=depth),size=2)+
      scale_colour_manual(values=c("cadetblue1","cadetblue3","cadetblue4"),
                          name="Depth (m)")+
      theme+theme(legend.position=c(0.15,0.15),
                  legend.title=element_text(),
                  legend.background=element_blank(),
                  legend.box.background=element_rect(colour="black",fill="white"))+
      ylab(expression("Temperature (°C)"))

# displays significant effects on the graph
label<-data.frame(label="depth ***\ndate ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.65,0.8,0.30,0.2,hjust=0)

graph_list$temperature=graph

# Oxygen # ----
model<-lmer(oxygen~fish*light*OM+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$oxygen=get_stat_table(stat)

# Tukey test
model<-lmer(oxygen~fish+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(14.75,16.2)

# figure
data<-summarySE(multi, measurevar="oxygen", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)
p1<-ggplot(data=data)+
      geom_boxplot(aes(light,oxygen,fill=OM))+
      geom_text(data=tukey,aes(x,y+0.2,label=.group),size=7,fontface = "bold")+
      geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=0.2)+
      facet_wrap(~fish)+
      colour_OM+
      theme+theme(legend.position=c(0.1,0.9),
                  legend.background=element_blank(),
                  legend.box.background=element_rect(colour="black",fill="white"))+
      ylab(expression("Oxygen concentration (mg L"^{-1}*")"))

# displays significant effects on the graph
label<-data.frame(label="fish ***\ndepth ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.65,0.15,0.35,0.25,hjust=0)

graph_list$oxygen=graph

# pH -------------
model<-lmer(pH~fish*light*OM+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$pH=get_stat_table(stat)

# Tukey test
model<-lmer(pH~fish+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(10.3,10.5)

# figure
data<-summarySE(multi, measurevar="pH", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)
p1<-ggplot(data=data)+
  geom_boxplot(aes(light,pH,fill=OM))+
  geom_text(data=tukey,aes(x,y+0.1,label=.group),size=7,fontface = "bold")+
  geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=0.1)+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position=c(0.1,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab("pH")

# displays significant effects on the graph
label<-data.frame(label="fish ***\ndepth *",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.7,0.1,0.3,0.2,hjust=0)

graph_list$pH=graph

# Turbidity -------------
data<-multi[multi$turbidity<20,]
model<-lmer(log(turbidity)~fish*light*OM+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=data)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$turbidity=get_stat_table(stat)

# Tukey test
model<-lmer(log(turbidity)~fish*light+depth+(1|mesocosm)+(1|date), na.action = "na.omit", REML=FALSE, data=multi)
tukey<-emmeans(model, specs = ~fish*light)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$y<-c(0.8,1.2,1.2,1.2)

# figure
data<-multi
data$turbidity<-log(data$turbidity)
data<-summarySE(data, measurevar="turbidity", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)
p1<-ggplot(data=data)+
      geom_boxplot(aes(light,turbidity,fill=OM))+
      geom_text(data=tukey,aes(light,y,label=.group),size=7,fontface = "bold")+
      facet_wrap(~fish)+
      colour_OM+
      theme+theme(legend.position=c(0.1,0.9),
                  legend.background=element_blank(),
                  legend.box.background=element_rect(colour="black",fill="white"))+
      ylab("Turbidity (log)")

# displays significant effects on the graph
label<-data.frame(label="fish ***\nlight ***\nfish:light **\ndepth ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.6,0.15,0.35,0.3,hjust=0)

graph_list$turbidity=graph

# Figure water physico-chemistry # ----
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(graph_list$temperature, 0, 1, 1, 1)+
  draw_plot(graph_list$oxygen, 1, 1, 1, 1)+
  draw_plot(graph_list$pH, 0, 0, 1, 1)+
  draw_plot(graph_list$turbidity, 1, 0, 1, 1)+
  draw_plot_label(c("A","B","C","D"),
                  c(0,1,0,1),
                  c(2,2,1,1), size = 30)
ggsave(paste(path_figures,"supp_multipar.pdf",sep=""), graph, width = 14, height = 10, device=cairo_pdf)

### Dissolved nutrients # ----
p1<-ggplot(data=nutri)+
  geom_boxplot(aes(date,N),fill="lightgrey")+
  theme+
  ylab(expression("Nitrogen (\u03BCg L"^{-1}*")"))

p2<-ggplot(data=nutri)+
  geom_boxplot(aes(date,P),fill="lightgrey")+
  theme+
  ylab(expression("Phosphorus (\u03BCg L"^{-1}*")"))

p3<-ggplot(data=nutri)+
  geom_boxplot(aes(date,NP),fill="lightgrey")+
  theme+
  ylab("N:P ratio")

graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  draw_plot_label(c("A","B","C"),
                  c(0,1,2),
                  c(1,1,1), size = 30)
ggsave(paste(path_figures,"supp_nutrient.pdf",sep=""), graph, width = 14, height = 4.5, device=cairo_pdf)

### Seston # ----
model<-lm(data=seston,log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statisticf="F")
stat_list$seston=get_stat_table(stat)

# Tukey test
model<-lm(data=seston,log(concentration)~fish*light)
tukey<-emmeans(model, specs = pairwise~fish*light)
tukey<-cld(tukey[[1]],
           alpha=0.05,
           Letters=letters)
tukey$y=c(7,12,12,12)

# figure
p1<-ggplot(data=seston)+
  geom_boxplot(aes(light,concentration,fill=OM))+
  geom_text(data=tukey,aes(light,y,label=.group),size=7,fontface = "bold")+
  facet_wrap(~fish)+
  colour_chloro+
  theme+theme(legend.position=c(0.15,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Seston dry concentration (mg "*L^{-1}*")"))

# displays significant effects on the graph
label<-data.frame(label="fish *\nlight *\nfish:light *",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,-0.02,0.5,0.4,0.2,hjust=0)

graph_list$seston=graph

ggsave(paste(path_figures,"supp_seston.pdf",sep=""), graph, width = 6, height = 5, device=cairo_pdf)

### Phytoplankton cytometry # ----
# Pico phytoplankton # ----
model<-glmer(pico~fish*light*OM+(1|mesocosm)+(1|sample_date), offset=log(volume_phyto),family=poisson,data=cyto_count)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$pico=get_stat_table(stat)

# figure
data<-cyto[cyto$taxon_cyto=="picophytoplankton",]
data<-summarySE(data, measurevar="cyto", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)

p1<-ggplot(data=data)+
  geom_boxplot(aes(light,cyto,fill=OM))+
  facet_wrap(~fish)+
  colour_chloro+
  theme+theme(legend.position=c(0.1,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Pico phytoplankton (cell µL"^{-1}*")"))

graph_list$pico=p1

# Nano phytoplankton # ----
model<-glmer(nano~fish*light*OM+(1|mesocosm)+(1|sample_date), offset=log(volume_phyto),family=poisson,data=cyto_count)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statistic="Chisq")
stat_list$nano=get_stat_table(stat)

# Tukey test
model<-glmer(nano~fish*light+(1|mesocosm)+(1|sample_date), offset=log(volume_phyto),family=poisson,data=cyto_count)
tukey<-emmeans(model, specs = pairwise~fish*light)
tukey<-cld(tukey$emmeans,
           alpha=0.05,
           Letters=letters)
tukey$y<-c(3.1,5.5,5.5,5.5)

# figure
data<-cyto[cyto$taxon_cyto=="nanophytoplankton",]
data<-summarySE(data, measurevar="cyto", groupvars=c("OM","light","fish","mesocosm"),na.rm=TRUE)

p1<-ggplot(data=data)+
  geom_boxplot(aes(light,cyto,fill=OM))+
  geom_text(data=tukey,aes(light,y,label=.group),size=7,fontface = "bold")+
  facet_wrap(~fish)+
  colour_chloro+
  theme+theme(legend.position=c(0.1,0.9),
              legend.background=element_blank(),
              legend.box.background=element_rect(colour="black",fill="white"))+
  ylab(expression("Nano phytoplankton (cell µL"^{-1}*")"))

# displays significant effects on the graph
label<-data.frame(label="fish ***\nlight ***\nfish:light ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.6,0.1,0.4,0.2,hjust=0)

graph_list$nano=graph

# Figure phytoplankton cytometry # ----
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(graph_list$pico, 0, 0, 1, 1)+
  draw_plot(graph_list$nano, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"),
                  c(0,1),
                  c(1,1), size = 30)
ggsave(paste(path_figures,"supp_cyto.pdf",sep=""), graph, width = 14, height = 6, device=cairo_pdf)

### Phytoplankton identification # ----
data<-melt(phyto[phyto$family_1!="none",c(3,6:ncol(phyto))], id.vars = c("family_1"),
           variable.name = "mesocosm", 
           value.name = "concentration")
levels(data$mesocosm)<-c(1:36)
data<-merge(treatment,data,by="mesocosm")
data$family_1<-as.factor(data$family_1)
data<-aggregate(data=data,concentration~family_1+mesocosm+fish+light+OM+treatment,FUN=sum)

databis<-aggregate(data=data,concentration~family_1+treatment,FUN=mean)

# Global abundance of the main taxa
p1<-ggplot(data=databis)+
  geom_bar(aes(treatment,concentration,fill=family_1),stat="identity")+
  scale_fill_brewer(palette="Paired")+
  theme+theme(axis.title.y=element_blank(),
              axis.title.x=element_text())+
  coord_flip()+
  xlab("")+
  ylab(label_algae)

graph_list$phyto_taxa=p1

# extract taxa colours
names<-levels(data$family_1)
colours<-brewer.pal(length(names),"Paired")

# Chlorophyceae
model<-lm(data=data[data$family_1=="Chlorophyceae",],log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$phyto_chlorophyceae=get_stat_table(stat)

p1<-ggplot(data=data[data$family_1=="Chlorophyceae",])+
  geom_boxplot(aes(light,concentration,fill=OM))+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position="none",
              plot.title=element_text(colour=colours[names=="Chlorophyceae"]))+
  y_axis_log10+
  ylab(label_algae)+
  ggtitle("Chlorophyceae")

graph_list$phyto_chlorophyceae=p1

# Cyanobacteria
model<-lm(data=data[data$family_1=="Cyanobacteria",],log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$phyto_cyanobacteria=get_stat_table(stat)

p1<-ggplot(data=data[data$family_1=="Cyanobacteria",])+
  geom_boxplot(aes(light,concentration,fill=OM))+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position="none",
              plot.title=element_text(colour=colours[names=="Cyanobacteria"]))+
  y_axis_log10+
  ylab(label_algae)+
  ggtitle("Cyanobacteria")

graph_list$phyto_cyanobacteria=p1

# Dinophyceae
model<-lm(data=data[data$family_1=="Dinophyceae",],log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$phyto_dinophyceae=get_stat_table(stat)

# Tukey test
model<-lm(data=data[data$family_1=="Dinophyceae",],log(concentration)~fish)
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(10^2.5,1e3)
tukey$.group<-as.factor(tukey$.group)
levels(tukey$.group)<-c("a","b")

p1<-ggplot(data=data[data$family_1=="Dinophyceae",])+
  geom_boxplot(aes(light,concentration,fill=OM))+
  geom_text(data=tukey,aes(x,y+c(150,500),label=.group),size=5,fontface = "bold")+
  geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=0.1)+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position="none",
              plot.title=element_text(colour=colours[names=="Dinophyceae"]))+
  y_axis_log10_short+
  ylab(label_algae)+
  ggtitle("Dinophyceae")

# displays significant effects on the graph
label<-data.frame(label="fish ***",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.7,0.65,0.25,0.2,hjust=0)

graph_list$phyto_dinophyceae=graph

# Trebouxiophyceae
model<-lm(data=data[data$family_1=="Trebouxiophyceae",],log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$phyto_trebouxiophyceae=get_stat_table(stat)

p1<-ggplot(data=data[data$family_1=="Trebouxiophyceae",])+
  geom_boxplot(aes(light,concentration,fill=OM))+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position="none",
              plot.title=element_text(colour=colours[names=="Trebouxiophyceae"]))+
  y_axis_log10+
  ylab(label_algae)+
  ggtitle("Trebouxiophyceae")

graph_list$phyto_trebouxiophyceae=p1

# Zygnematophyceae
model<-lm(data=data[data$family_1=="Zygnematophyceae",],log(concentration)~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$phyto_zygnematophyceae=get_stat_table(stat)

# Tukey test
model<-lm(data=data[data$family_1=="Zygnematophyceae",],log(concentration)~fish*light)
tukey<-emmeans(model, specs = pairwise~fish:light)
tukey<-cld(tukey$emmeans,
           alpha=0.05,
           Letters=letters)
tukey$y<-c(10^3.1,10^3.6,10^3.6,10^3.6)

p1<-ggplot(data=data[data$family_1=="Zygnematophyceae",])+
  geom_boxplot(aes(light,concentration,fill=OM))+
  geom_text(data=tukey,aes(light,y,label=.group),size=5,fontface = "bold")+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position="none",
              plot.title=element_text(colour=colours[names=="Zygnematophyceae"]))+
  y_axis_log10+
  ylab(label_algae)+
  ggtitle("Zygnematophyceae")

# displays significant effects on the graph
label<-data.frame(label="fish *\nlight *\nfish:light *",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.6,0.15,0.35,0.2,hjust=0)

graph_list$phyto_zygnematophyceae=graph

# Shannon-Weaver diversity
model<-lm(data=diversity,diversity~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[1:7,c("F value","Df","Pr(>F)")]
stat_list$phyto_shannon=get_stat_table(stat)

# Tukey test
model<-lm(data=diversity,diversity~fish)
tukey<-emmeans(model, specs = ~fish)
tukey<-cld(tukey,
           alpha=0.05,
           Letters=letters)
tukey$x=1.5
tukey$xmin=0.5
tukey$xmax=2.5
tukey$y=c(2,2.6)
tukey$.group<-as.factor(tukey$.group)
levels(tukey$.group)<-c("a","b")

p1<-ggplot(data=diversity)+
  geom_boxplot(aes(light,diversity,fill=OM))+
  geom_text(data=tukey,aes(x,y+0.1,label=.group),size=5,fontface = "bold")+
  geom_errorbarh(data=tukey,aes(xmin=xmin,xmax=xmax,y=y),height=0.1)+
  facet_wrap(~fish)+
  colour_OM+
  theme+theme(legend.position="none")+
  ylab("Shannon index")+
  ggtitle("")

# displays significant effects on the graph
label<-data.frame(label="fish *",x=1,y=1)
label<-ggplot(data=label)+
  geom_label(aes(x=x,y=y,label=label),hjust=0,size=6)+
  theme_void()

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(label,0.75,0.65,0.25,0.2,hjust=0)

graph_list$phyto_shannon=graph

# final graph
p1<-p1+theme(legend.position="right")
legend<-get_legend(p1)

graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 3)) +
  draw_plot(graph_list$phyto_taxa, 0.1, 2, 2.9, 1)+
  draw_plot(graph_list$phyto_chlorophyceae, 0, 1, 1, 1)+
  draw_plot(graph_list$phyto_cyanobacteria, 1, 1, 1, 1)+
  draw_plot(graph_list$phyto_dinophyceae, 2, 1, 1, 1)+
  draw_plot(graph_list$phyto_trebouxiophyceae, 0, 0, 1, 1)+
  draw_plot(graph_list$phyto_zygnematophyceae, 1, 0, 1, 1)+
  draw_plot(graph_list$phyto_shannon, 2, 0, 1, 1)+
  draw_plot(legend, 2.8, 1.95, 0.2, 0.2)+
  draw_plot_label(c("A","B","C","D","E","F","G"), c(0,0,1,2,0,1,2), c(3,2,2,2,1,1,1), size = 30)
ggsave(paste(path_figures,"supp_phytoplankton_family.pdf",sep=""),graph, width = 20, height = 14, device = cairo_pdf)

### Zooplankton identification # ----
# Copepoda # ----
# global
model<-glmmPQL(data=zoo, copepode~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$copepoda=get_stat_table(stat)
model<-glmmPQL(data=zoo, copepode~light+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
summary(model)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], copepode~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$copepoda_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], copepode~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$copepoda_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,copepode,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Copepoda")

graph_list$copepoda=p1

# Cladocera # ----
# global
model<-glmmPQL(data=zoo, cladocerae~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$cladocera=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], cladocerae~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$cladocera_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], cladocerae~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$cladocera_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,cladocerae,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Cladocera")

graph_list$cladocera=p1

# Chaoborus # ----
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,chaoborus,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  ylab(label_zoo)+
  ggtitle("Chaoborus")

graph_list$chaoborus=p1

# Rotifers # ----
# global
model<-glmmPQL(data=zoo, rotifer~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$rotifer=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], rotifer~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$rotifer_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], rotifer~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$rotifer_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,rotifer,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Rotifers")

graph_list$rotifer=p1

# Figure aggregated taxa # ----
p1<-p1+theme(legend.position="right")
legend<-get_legend(p1)

graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 2)) +
  draw_plot(graph_list$copepoda, 0, 1, 1, 1)+
  draw_plot(graph_list$cladocera, 1, 1, 1, 1)+
  draw_plot(graph_list$chaoborus, 0, 0, 1, 1)+
  draw_plot(graph_list$rotifer, 1, 0, 1, 1)+
  draw_plot(legend, 2, 1.4, 0.2, 0.2)+
  draw_plot(legend, 2, 0.4, 0.2, 0.2)+
  draw_plot_label(c("A","B","C","D"),
                  c(0,1,0,1),
                  c(2,2,1,1), size = 30)
ggsave(paste(path_figures,"supp_zoo_aggregated.pdf",sep=""), graph, width = 14, height = 12, device=cairo_pdf)

# Cyclopids # ----
# global
model<-glmer(data=zoo, cyclopide~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, cyclopide~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$cyclopidae=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], cyclopide~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$cyclopidae_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], cyclopide~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$cyclopidae_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,cyclopide,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Cyclopidae")

graph_list$cyclopidae=p1

# Nauplii # ----
# global
model<-glmer(data=zoo, nauplius~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, nauplius~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$nauplii=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], nauplius~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$nauplii_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], nauplius~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$nauplii_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,nauplius,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Nauplii")

graph_list$nauplii=p1

# Chydoridae # ----
# global
model<-glmer(data=zoo, chydoridae~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, chydoridae~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$chydoridae=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], chydoridae~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$chydoridae_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], chydoridae~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$chydoridae_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,chydoridae,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Chydoridae")

graph_list$chydoridae=p1

# Lecane # ----
# global
model<-glmer(data=zoo, lecane~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, lecane~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$lecane=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], lecane~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$lecane_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], lecane~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$lecane_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,lecane,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Lecane")

graph_list$lecane=p1

# Lepadella # ----
# global
model<-glmer(data=zoo, lepadella~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, lepadella~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$lepadella=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], lepadella~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$lepadella_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], lecane~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$lepadella_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,lepadella,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Lepadella")

graph_list$lepadella=p1

# Bdelloide # ----
# global
model<-glmer(data=zoo, bdelloide~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, bdelloide~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$bdelloide=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], bdelloide~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$bdelloide_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], bdelloide~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$bdelloide_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,bdelloide,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Bdelloide")

graph_list$bdelloide=p1

# Polyarthra # ----
# global
model<-glmer(data=zoo, polyarthra~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, polyarthra~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$polyarthra=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], polyarthra~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$polyarthra_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], polyarthra~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$polyarthra_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,polyarthra,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Polyarthra")

graph_list$polyarthra=p1

# Keratella # ----
# global
model<-glmer(data=zoo, keratella~fish*light*OM+date+(1|mesocosm), family = poisson)
dispersion_glmer(model)
model<-glmmPQL(data=zoo, keratella~fish*light*OM+date, random = ~1|mesocosm, family = quasipoisson(link='log'))
plot(model)
qqnorm(resid(model))
qqline(resid(model))
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$keratella=get_stat_table(stat)

# per sample date
model<-glm(data=zoo[zoo$date=="19/08/2015",], keratella~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$keratella_aug=get_stat_table(stat)

model<-glm(data=zoo[zoo$date=="25/09/2015",], keratella~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$keratella_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo)+
  geom_boxplot(aes(light,keratella,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Keratella")

graph_list$keratella=p1

# Brachionus # ----
model<-glm(data=zoo[zoo$date=="25/09/2015",], brachionus~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$brachionus_sep=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo[zoo$date=="25/09/2015",])+
  geom_boxplot(aes(light,brachionus,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Brachionus")

graph_list$brachionus=p1

# Anuraeopsis # ----
model<-glm(data=zoo[zoo$date=="25/09/2015",], anuraeopsis~fish*light*OM, family = quasipoisson(link='log'))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statisticf="Chisq")
stat_list$anuraeopsis=get_stat_table(stat)

# figure
p1<-ggplot(data=zoo[zoo$date=="25/09/2015",])+
  geom_boxplot(aes(light,anuraeopsis,fill=OM))+
  facet_grid(date~fish,scale="free")+
  colour_zoo+
  theme+theme(legend.position="none")+
  y_axis_log10+
  ylab(label_zoo)+
  ggtitle("Anuraeopsis")

graph_list$anuraeopsis=p1

# Figures # ----
p1<-p1+theme(legend.position="right")
legend<-get_legend(p1)

# large zooplankton
graph<-ggdraw(xlim = c(0, 3.2), ylim = c(0, 1)) +
  draw_plot(graph_list$cyclopidae, 0, 0, 1, 1)+
  draw_plot(graph_list$nauplii, 1, 0, 1, 1)+
  draw_plot(graph_list$chydoridae, 2, 0, 1, 1)+
  draw_plot(legend, 3, 0.4, 0.2, 0.2)+
  draw_plot_label(c("A","B","C"),
                  c(0,1,2),
                  c(1,1,1), size = 30)
ggsave(paste(path_figures,"supp_zoo_large.pdf",sep=""), graph, width = 20, height = 7, device=cairo_pdf)

# small rotifers
graph<-ggdraw(xlim = c(0, 3.2), ylim = c(0, 1)) +
  draw_plot(graph_list$lecane, 0, 0, 1, 1)+
  draw_plot(graph_list$lepadella, 1, 0, 1, 1)+
  draw_plot(graph_list$bdelloide, 2, 0, 1, 1)+
  draw_plot(legend, 3, 0.4, 0.2, 0.2)+
  draw_plot_label(c("A","B","C"),
                  c(0,1,2),
                  c(1,1,1), size = 30)
ggsave(paste(path_figures,"supp_rotifer_small.pdf",sep=""), graph, width = 20, height = 7, device=cairo_pdf)

# large rotifers
graph<-ggdraw(xlim = c(0, 3.2), ylim = c(0, 1)) +
  draw_plot(graph_list$polyarthra, 0, 0, 1, 1)+
  draw_plot(graph_list$keratella, 1, 0, 1, 1)+
  draw_plot(graph_list$anuraeopsis, 2, 0, 1, 0.5)+
  draw_plot(graph_list$brachionus, 2, 0.5, 1, 0.5)+
  draw_plot(legend, 3, 0.4, 0.2, 0.2)+
  draw_plot_label(c("A","B","C","D"),
                  c(0,1,2,2),
                  c(1,1,1,0.5), size = 30)
ggsave(paste(path_figures,"supp_rotifer_large.pdf",sep=""), graph, width = 20, height = 8, device=cairo_pdf)

### Functional diversity of bacteria - Ecoplates # ----
# amines # ----
# pelagic
model<-lm(data=eco[eco$family=="amines" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+substrat+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,6:9),c("F value","Df","Pr(>F)")]
stat_list$eco_amine_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$family=="amines" & eco$compartment=="benthic",],
          absorbance~fish*light*OM+substrat)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,5:8),c("F value","Df","Pr(>F)")]
stat_list$eco_amine_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$family=="amines",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("Amines")

graph_list$eco_amine=p1

# amino acids # ----
# pelagic
model<-lm(data=eco[eco$family=="amino acids" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+substrat+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,6:9),c("F value","Df","Pr(>F)")]
stat_list$eco_aminoacid_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$family=="amino acids" & eco$compartment=="benthic",],
          absorbance~fish*light*OM+substrat)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,5:8),c("F value","Df","Pr(>F)")]
stat_list$eco_aminoacid_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$family=="amino acids",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("Amino acids")

graph_list$eco_aminoacid=p1

# carbohydrate # ----
model<-lm(data=eco[eco$family=="carbohydrate" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+substrat+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,6:9),c("F value","Df","Pr(>F)")]
stat_list$eco_carbohydrate_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$family=="carbohydrate" & eco$compartment=="benthic",],
          absorbance~fish*light*OM+substrat)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,5:8),c("F value","Df","Pr(>F)")]
stat_list$eco_carbohydrate_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$family=="carbohydrate",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("Carbohydrate")

graph_list$eco_carbohydrate=p1

# carboxylic acids # ----
# pelagic
model<-lm(data=eco[eco$family=="carboxylic acids" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+substrat+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,6:9),c("F value","Df","Pr(>F)")]
stat_list$eco_carboacid_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$family=="carboxylic acids" & eco$compartment=="benthic",],
          absorbance~fish*light*OM+substrat)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,5:8),c("F value","Df","Pr(>F)")]
stat_list$eco_carboacid_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$family=="carboxylic acids",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("Carboxylic acids")

graph_list$eco_carboacid=p1

# phenolic acids # ----
# pelagic
model<-lm(data=eco[eco$family=="phenolic acids" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+substrat+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,6:9),c("F value","Df","Pr(>F)")]
stat_list$eco_phenol_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$family=="phenolic acids" & eco$compartment=="benthic",],
          absorbance~fish*light*OM+substrat)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,5:8),c("F value","Df","Pr(>F)")]
stat_list$eco_phenol_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$family=="phenolic acids",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("Phenolic acids")

graph_list$eco_phenol=p1

# polymer # ----
# pelagic
model<-lm(data=eco[eco$family=="polymer" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+substrat+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,6:9),c("F value","Df","Pr(>F)")]
stat_list$eco_polymer_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$family=="polymer" & eco$compartment=="benthic",],
          absorbance~fish*light*OM+substrat)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:3,5:8),c("F value","Df","Pr(>F)")]
stat_list$eco_polymer_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$family=="polymer",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("Polymer")

graph_list$eco_polymer=p1

# Figure substrate families # ----
p1<-p1+theme(legend.position="right")
legend<-get_legend(p1)

graph<-ggdraw(xlim = c(0, 3.25), ylim = c(0, 2)) +
  draw_plot(graph_list$eco_amine, 0, 1, 1, 1)+
  draw_plot(graph_list$eco_aminoacid, 1, 1, 1, 1)+
  draw_plot(graph_list$eco_carbohydrate, 2, 1, 1, 1)+
  draw_plot(graph_list$eco_carboacid, 0, 0, 1, 1)+
  draw_plot(graph_list$eco_phenol, 1, 0, 1, 1)+
  draw_plot(graph_list$eco_polymer, 2, 0, 1, 1)+
  draw_plot(legend, 3, 1.4, 0.25, 0.2)+
  draw_plot(legend, 3, 0.4, 0.25, 0.2)+
  draw_plot_label(c("A","B","C","D","E","F"),
                  c(0,1,2,0,1,2),
                  c(2,2,2,1,1,1), size = 30)
ggsave(paste(path_figures,"supp_ecoplate_substrats.pdf",sep=""), graph, width = 18, height = 12, device=cairo_pdf)

# alpha-cyclodextrine # ----
# pelagic
model<-lm(data=eco[eco$substrat=="E1" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:7),c("F value","Df","Pr(>F)")]
stat_list$eco_alphacyclo_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$substrat=="E1" & eco$compartment=="benthic",],
          absorbance~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:7),c("F value","Df","Pr(>F)")]
stat_list$eco_alphacyclo_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$substrat=="E1",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle(expression(alpha*"-cyclodextrine"))

graph_list$eco_alphacyclo=p1

# D-cellobiose # ----
# pelagic
model<-lm(data=eco[eco$substrat=="G1" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:7),c("F value","Df","Pr(>F)")]
stat_list$eco_cellobiose_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$substrat=="G1" & eco$compartment=="benthic",],
          absorbance~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:7),c("F value","Df","Pr(>F)")]
stat_list$eco_cellobiose_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$substrat=="G1",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("D-cellobiose")

graph_list$eco_cellobiose=p1

# L-asparagine # ----
# pelagic
model<-lm(data=eco[eco$substrat=="B4" & eco$compartment=="pelagic",],
          absorbance~fish*light*OM+date)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:7),c("F value","Df","Pr(>F)")]
stat_list$eco_asparagine_pel=get_stat_table(stat)

# benthic
model<-lm(data=eco[eco$substrat=="B4" & eco$compartment=="benthic",],
          absorbance~fish*light*OM)
par(mfrow=c(2,2))
plot(model)
summary(model)
stat<-Anova(model,type=2,test.statistic="F")
stat<-stat[c(1:7),c("F value","Df","Pr(>F)")]
stat_list$eco_asparagine_ben=get_stat_table(stat)

# figure
p1<-ggplot(data=eco[eco$substrat=="B4",])+
  geom_boxplot(aes(light,absorbance,fill=OM))+
  facet_grid(compartment~fish,scale="free")+
  colour_bacteria+
  theme+theme(legend.position="none")+
  ylab(label_OD)+
  ggtitle("L-asparagine")

graph_list$eco_asparagine=p1

# Figure substrate DOC treatment # ----
p1<-p1+theme(legend.position="right")
legend<-get_legend(p1)

graph<-ggdraw(xlim = c(0, 3.25), ylim = c(0, 1)) +
  draw_plot(graph_list$eco_carboacid, 0, 0, 1, 1)+
  draw_plot(graph_list$eco_phenol, 1, 0, 1, 1)+
  draw_plot(graph_list$eco_polymer, 2, 0, 1, 1)+
  draw_plot(legend, 3, 0.4, 0.25, 0.2)+
  draw_plot_label(c("A","B","C"),
                  c(0,1,2),
                  c(1,1,1), size = 30)
ggsave(paste(path_figures,"supp_ecoplate_substrats_DOC.pdf",sep=""), graph, width = 18, height = 6, device=cairo_pdf)
