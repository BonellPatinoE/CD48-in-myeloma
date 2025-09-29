library(tidyverse)
library(dslabs)
library(dplyr)
library(ggrepel)
library(ggthemes)
library(ggplot2)
library(reshape2)
library(gprofiler2)
library(ggpubr)
library(naniar)
library(ggpubr)

## Loading the coMMpass data set ##

coMMpass<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
coMMpass<-data.frame(coMMpass)

colnames(coMMpass) 
coMMpass1<-coMMpass

rownames(coMMpass1) <- coMMpass1[,1]
coMMpass1[,1] <- NULL
head(coMMpass1)

## gconvert to get gene names ##

coMMpass2<-gconvert(row.names(coMMpass1), organism = "hsapiens", target="HGNC", mthreshold = 1, filter_na = FALSE)
coMMpass2

coMMpass3<-coMMpass1
rownames(coMMpass3)<-NULL
coMMpass3$Gene_name<-coMMpass2$name
nrow(coMMpass3)

coMMpass4<-coMMpass3 %>% distinct(Gene_name, .keep_all = TRUE)

coMMpass4<-coMMpass4%>%filter(!is.na(Gene_name))

row.names(coMMpass4)<-coMMpass4$Gene_name
coMMpass5<-coMMpass4  
coMMpass5
coMMpass5$Gene_name<-NULL

coMMpass6<-t(coMMpass5)
coMMpass6<-data.frame(coMMpass6)

class(coMMpass6)
colnames(coMMpass6)
coMMpass6$CD70

coMMpass7<-log2(coMMpass6+1)
coMMpass7[1:10, 1:10]

coMMpass6%>%select(CD48, KDM6A)

coMMpass6[1:10, 1:10]
ncol(coMMpass6)
nrow(coMMpass6)
log2(coMMpass6)+1

# Define your dataset 'coMMpass7' and the gene of interest 'CD48'
gene_of_interest <- "CD48"  # Gene of interest

# Initialize an empty list to store correlation results
correlation_results <- list()

# Iterate over all genes except the gene of interest
for (gene in names(coMMpass7)[-which(names(coMMpass7) == gene_of_interest)]) {
  # Perform correlation test
  correlation_test <- cor.test(coMMpass7[, gene_of_interest], coMMpass7[, gene])
  
  # Extract correlation coefficient and p-value
  correlation_coefficient <- correlation_test$estimate
  p_value <- correlation_test$p.value
  
  # Store results in a list
  correlation_results[[gene]] <- c(Correlation = correlation_coefficient, P_Value = p_value)
}

# Convert the list to a data frame
correlation_table <- do.call(rbind, correlation_results)

# Print the correlation table
print(correlation_table)

## Correlation plots ##

coMMpass7%>%ggplot(aes(CD48))+
  geom_histogram()

coMMpass7%>%ggplot(aes(UCK2))+
  geom_histogram()

coMMpass7%>%ggplot(aes(IL12A))+
  geom_histogram()

ggscatter(coMMpass6, x = "CD48", y = "MEF2B", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD48 log2", ylab = "MEF2B log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs MEF2B (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "IRF4", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "IRF4 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs IRF4 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "PAX5", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "PAX5 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs PAX5 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "SPI1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "SPI1 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs SPI1 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "TAF1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "TAF1 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs TAF1 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "TBP", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "TBP log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs TBP (Pearson R value)")


ggscatter(coMMpass6, x = "CD48", y = "YY1", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "YY1 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs YY1 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "TCF12", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "TCF12 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs TCF12 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "IRF4", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "CD48 log2", ylab = "IRF4 log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs IRF4 (Pearson R value)")

ggscatter(coMMpass6, x = "CD48", y = "MEF2C", color="darkblue",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightblue"),
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "CD48 log2", ylab = "MEF2C log2")+
  yscale("log2", .format = TRUE)+
  xscale("log2", .format = TRUE)+
  ggtitle("CD48 vs MEF2C (Pearson R value)")

comparison<-coMMpass6%>%select(c("CD48", "CD38", "TNFRSF17"))
comparison<-log2(comparison)+1
comparison$Sample_ID<-rownames(comparison)
rownames(comparison)<-NULL

# Convert to long format for ggplot
long_data <- melt(comparison, id.vars = "Sample_ID", variable.name = "Gene", value.name = "Expression")


# Create plot
ggplot(long_data, aes(x = Gene, y = Expression, fill= Gene, group=Gene)) +
  geom_boxplot()+
  yscale("log2")
