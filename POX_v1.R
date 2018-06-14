# Alexander Kutschera 23.05.2018

library(ggplot2)
library(XLConnect)
library(plyr)
library(reshape2)
library(grid)

inputfile <- "ALK180523_1.xlsx"
layout <- "ALK180523_1_layout.xlsx"



### reading the layout
raw_plate_layout <- readWorksheetFromFile(layout, sheet=1)
plate_layout <- mutate(raw_plate_layout,
                       Row=as.numeric(match(toupper(substr(well, 1, 1)), LETTERS)),
                       Column=as.numeric(substr(well, 2, 5)))

raw_plate_layout <- raw_plate_layout[complete.cases(raw_plate_layout),] #delete rows with NAs

elicitors <- unique(raw_plate_layout$elicitor)
genotypes <- unique(raw_plate_layout$genotype)

### reading the data
rawdata <- readWorksheetFromFile(inputfile, sheet=1, header = TRUE)

### calculations

sdom.calculate <- function(x){
  xsdom <- sd(x)/sqrt(length(x));
  return(xsdom)
}

divide.by.control <- function(x, dataset_mean){ #x is one column of dataset_mean
  values <- as.numeric(x[3:nrow(dataset_mean)])
  genotype <- x[1]
  blank <- as.numeric(dataset_mean[,dataset_mean[1,]==genotype & dataset_mean[2,]=="control"][3:nrow(dataset_mean)])
  normed_mean = values/blank
  return(normed_mean)
}

#ggplots formatting
annotate_type <- function(x, data_info, type){
  type <- data_info[,colnames(data_info) == x["name"]][rownames(data_info) == type ]
  return(type)
}

#ggplots formatting
add_sd <- function(x, data_info, type){
  type <- data_info[,colnames(data_info) == x["name"]][rownames(data_info) == type ]
  return(type)
}

dataset_mean <- matrix(nrow = 1, ncol = length(elicitors)*length(genotypes), dimnames = NULL)
dataset_raw <- matrix(nrow = 96/(length(elicitors)*length(genotypes)), ncol = length(elicitors)*length(genotypes), dimnames = NULL)
dataset_sd <- matrix(nrow = 1, ncol = length(elicitors)*length(genotypes), dimnames = NULL)
dataset_sdom <- matrix(nrow = 1, ncol = length(elicitors)*length(genotypes), dimnames = NULL)
dataset_info <- matrix(nrow = 2, ncol = length(elicitors)*length(genotypes), dimnames = NULL)

raw_graph <- matrix(nrow = length(elicitors)*length(genotypes)*(96/(length(elicitors)*length(genotypes))), ncol = 3, dimnames = NULL)

columnames <- c()
i = 1
i2 = 1
for(eli in elicitors){
  for(gen in genotypes){
    wells <- raw_plate_layout$well[raw_plate_layout$elicitor == eli & raw_plate_layout$genotype == gen]
    current_set <- rawdata$Mean[rawdata$Well %in% wells]
    
    dataset_mean[,i] <- mean(rawdata$Mean[rawdata$Well %in% wells]) #record mean
    dataset_sdom[,i] <- sdom.calculate(rawdata$Mean[rawdata$Well %in% wells]) #record SDOM
    dataset_sd[,i] <- sd(rawdata$Mean[rawdata$Well %in% wells]) #record SD
    dataset_info[1,i] <- gen
    dataset_info[2,i] <- eli
    
    bg_wells <- raw_plate_layout$well[raw_plate_layout$elicitor == "control" & raw_plate_layout$genotype == gen]
    current_background <- rawdata$Mean[rawdata$Well %in% bg_wells]
    
    dataset_raw[,i] <- current_set
    
    raw_graph[i2:(i2+length(current_set)-1),1] <- rep(eli, length(current_set))
    raw_graph[i2:(i2+length(current_set)-1),2] <- rep(gen, length(current_set))
    raw_graph[i2:(i2+length(current_set)-1),3] <- current_set/mean(current_background)

    
    i = i+1
    i2 = i2 + length(current_set)
    columnames <- append(columnames, paste(gen, eli, sep = " "))
    
    
  }
}
dataset_mean <- rbind(dataset_info, dataset_mean)
dataset_raw <- rbind(dataset_info, dataset_raw)
#  dataset_sd <- rbind(dataset_info, dataset_sd)

colnames(dataset_mean) <- columnames
colnames(dataset_raw) <- columnames
colnames(dataset_sd) <- columnames
colnames(dataset_sdom) <- columnames
colnames(dataset_info) <- columnames
rownames(dataset_info) <- c("genotype", "elicitor")

mean_fold_induction_values <- apply(dataset_mean, 2, divide.by.control, dataset_mean = dataset_mean)

mean_fold_induction <- dataset_mean #dublicate dataset_mean to use frost two rows

mean_fold_induction[3:nrow(dataset_mean),] <- mean_fold_induction_values
mean_fold_induction <- rbind(mean_fold_induction, dataset_sd)

rownames(mean_fold_induction) <- c("genotype", "elicitor", "value", "sd")
t_mean_fold_induction <- t(mean_fold_induction)
num_values <- t_mean_fold_induction[,3:4]
class(num_values) <- "numeric" # sonst werden max und se nicht alszahlen erkannt...

df_induction <- data.frame(
  eli = t_mean_fold_induction[,2],
  gen = factor(t_mean_fold_induction[,1], levels = genotypes),
  fold_induction = num_values[,1],
  se = num_values[,2]
)

bpmax <- ggplot(df_induction, aes(fill=gen, y=fold_induction, x=eli)) + facet_wrap(~gen, ncol = 1) + ggtitle("Fold induction")
bpmax <- bpmax + geom_bar(position="dodge", stat="identity")
bpmax <- bpmax + geom_errorbar(aes(ymin=fold_induction-se, ymax=fold_induction+se),
                               width=.2,                    # Width of the error bars
                               position=position_dodge(.9))
bpmax <- bpmax + labs(x="", y="L/Lmax") + theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
                                                legend.title=element_blank(),
                                                panel.background = element_rect(fill = "white", colour = "grey90"),
                                                panel.grid.major = element_line(color = "grey90"),
                                                strip.background = element_rect(fill = "grey90", colour = NA),
                                                panel.grid.minor = element_line(colour = "grey90", size = 0.25))


show(bpmax)

num_raw_graph <- raw_graph[,3]
class(num_raw_graph) <- "numeric" # sonst werden max und se nicht alszahlen erkannt...
raw_graph_bp <- data.frame(
  eli = raw_graph[,1],
  gen = factor(raw_graph[,2], levels = genotypes),
  induction = num_raw_graph
)

bp <- ggplot(raw_graph_bp, aes(fill=gen, y=induction, x=eli)) + ggtitle("Fold induction")
bp <- bp + geom_boxplot()
bp <- bp + labs(x="", y="L/Lmax") + theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
                                                legend.title=element_blank(),
                                                panel.background = element_rect(fill = "white", colour = "grey90"),
                                                panel.grid.major = element_line(color = "grey90"),
                                                strip.background = element_rect(fill = "grey90", colour = NA),
                                                panel.grid.minor = element_line(colour = "grey90", size = 0.25))

show(bp)
### plotting the layout
layout_plot <- ggplot(data=plate_layout, aes(x=Column, y=Row)) +
  geom_point(size=10) +
  geom_point(size=8, aes(colour = elicitor)) +
  geom_point(size=3, aes(shape = genotype)) +
  scale_alpha_discrete(range = c(0.4, 0.05)) +
  labs(title="Plate Layout") +
  scale_color_discrete(na.value="white")

layout_plot <- layout_plot +
  coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.8, 12.2), ylim=c(0.6, 8.4)) +
  scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
  scale_x_continuous(breaks=seq(1, 12))

layout_plot <- layout_plot + theme_bw() + guides(colour = guide_legend(title = "Elicitors", ncol = 2, byrow = TRUE), alpha = guide_legend(title = "Genotypes", ncol = 2, byrow = TRUE)) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.key=element_blank()
  )

show(layout_plot)

