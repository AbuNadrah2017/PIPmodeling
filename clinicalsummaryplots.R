options(warn=-1)
library(gtools)
library("plyr")

library(ggside)
library(ggcorrplot)
library(ggstatsplot)
library(ggplot2)
library(PMCMRplus)




clinicalModel = NULL

# train dataset processed
source('Tables 1.R')

train_data = data


data_omit_NA = data



for(i in 2:ncol(data_omit_NA)){
  
  num_classes <- c("integer", "numeric", "double")
  str_classes <- c("factor", "character")
  
  
  if(class(data_omit_NA[,i]) %in% num_classes){
    
    ggbetweenstats(
      data  = data_omit_NA,
      x     = response,
      y     = names(data_omit_NA[i]),
      type = "parametric",
      title = names(data_omit_NA[i])) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
    
  }
    
    
  if(coarse) {
    num_classes <- c("integer", "numeric", "double")
    str_classes <- c("factor", "character")
    out <- ifelse(out %in% num_classes, "numeric", out)
    out <- ifelse(out %in% str_classes, "string", out)
    out <- ifelse(out %in% c("numeric", "string"), out, "other")
  }
  
  
}


## for reproducibility
set.seed(124)

c = ggbarstats(data_omit_NA,
               
               x = Treatment,
               
               y = response,
               
               bf.message = FALSE,
               
               proportion.test = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  scale_fill_discrete(name = "Treatment") +
  
  
  labs(title="Treatment",  x="response", y = "Proportions",
       caption = NULL) +
  
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 



## Gender
d = ggbarstats(data_omit_NA,
               
               x = Sex,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  scale_fill_discrete(name = "Gender") +
  
  labs(title="Gender",  x="response", y = "Proportions",
       caption = NULL) +
  
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 




## Line of ICB
e = ggbarstats(data_omit_NA,
               
               x = Line,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  scale_fill_discrete(name = "Line of ICB") +
  
  labs(title="Line of ICB",  x="response", y = "Proportions",
       caption = NULL) +
 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 





## Primary.melanoma.site.JW
f = ggbarstats(data_omit_NA,
               
               x = `Primary melanoma`,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  labs(title="Primary melanoma site JW",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Primary melanoma site JW") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 








## "Cutaneous primary YES/NO"
g = ggbarstats(data_omit_NA,
               
               x = `Cutaneous primary`,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  labs(title="Cutaneous primary YES/NO",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Cutaneous primary YES/NO") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 



## Previous MAPKi Yes/No
h = ggbarstats(data_omit_NA,
               
               x = `Previous MAPKi`,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  labs(title="Previous MAPKi Yes/No",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Previous MAPKi Yes/No") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 




## M Stage at Entry
i= ggbarstats(data_omit_NA,
              
              x = `M Stage`,
              
              y = response,
              
              bf.message = FALSE,
              
              package      = "wesanderson",
              
              palette      = "Darjeeling2",
              
              ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  labs(title="M Stage at Entry",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "M Stage at Entry") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 





## "Baseline LDH"
j =  ggbarstats(data_omit_NA,
                
                x = `Baseline LDH`,
                
                y = response,
                
                bf.message = FALSE,
                
                package      = "wesanderson",
                
                palette      = "Darjeeling2",
                
                ggtheme      = ggthemes::theme_tufte(base_size = 12) ) + 
  
  labs(title="Baseline LDH",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Baseline LDH") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 



## "ECOG"  LDH"
k = ggbarstats(data_omit_NA,
               
               x = ECOG,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) +
  
  labs(title="ECOG",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "ECOG") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 








## "Brain metastasis"
l = ggbarstats(data_omit_NA,
               
               x = `Brain metastasis`,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) +
  
  labs(title="Brain metastasis",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Brain metastasis") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 



## "Brain metastasis"
m = ggbarstats(data_omit_NA,
               
               x = `Lung metastasis`,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) +
  
  labs(title="Lung metastasis",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Lung metastasis") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 



## "Brain metastasis"
n = ggbarstats(data_omit_NA,
               
               x = `Liver metastasis`,
               
               y = response,
               
               bf.message = FALSE,
               
               package      = "wesanderson",
               
               palette      = "Darjeeling2",
               
               ggtheme      = ggthemes::theme_tufte(base_size = 12) ) +
  
  labs(title="Liver metastasis",  x="response", y = "Proportions",
       caption = NULL) +
  scale_fill_discrete(name = "Liver metastasis") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
        legend.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(colour = "black", face = "bold", size = 12),
        axis.text = element_text(colour = "black", face='bold', size = 12),
        axis.text.y = element_text(colour = "black", face='bold', size = 12),
        axis.text.x = element_text(colour = "black", face='bold', size = 12),
        axis.title.y = element_text(colour = "black", face='bold', size = 12),
        axis.title.x = element_text(colour = "black", face='bold', size = 12),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", colour = "black",  size = 12)    # Subtitle customization
  ) 





### Age
my_comparisons = alply(combinations(2 ,2, 1:2),1)

o=ggbetweenstats(
  data  = data_omit_NA,
  x     = response,
  y     = Age,
  type = "parametric",
  title = "Age at start of treatment") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))



p=ggbetweenstats(
  data  = data_omit_NA,
  x     = response,
  y     = `Neutro Lympho ratio`,
  type = "parametric",
  title = "Neutro Lympho ratio") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))



p=ggbetweenstats(
  data  = data_omit_NA,
  x     = response,
  y     = `Hemoglobin`,
  type = "parametric",
  title = "Hemoglobin") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))




svg(paste0("RNAClincal",  ".svg"), width=22, height=20)
gridExtra::grid.arrange(c,
             d, e, f,
             g, h, i, j, k,l, m, n, o,p,
             nrow=4,
             ncol=4)
dev.off()





#https://cran.r-project.org/web/packages/ggstatsplot/readme/README.html
#https://reposhub.com/python/data-validation/IndrajeetPatil-ggstatsplot.html
#https://indrajeetpatil.github.io/ggstatsplot/index.html
#https://www.fatalerrors.org/a/0tp10T4.html
