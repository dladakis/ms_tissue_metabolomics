library(tidyverse)
library(readxl)
library(GGally)
library(lubridate)
library(lme4)
library(ggpubr)
library(anytime)
library(directlabels)
library(nlme)
library(impute)
library(janitor)
library(pathifier)
library(Numero)
library(WGCNA)
library(geepack)
library(GSA)
library(rlang)
library(ggbeeswarm)
library(DescTools)
library(patchwork)
library(ggrepel)
library(gee)
library(fgsea)
library(RColorBrewer)
library(ggsci)
library(scales)
library(rstatix)
library(AnnotationDbi)
library(GO.db)
library(ReactomePA)
library(gplots)
library(org.Hs.eg.db)
library(clusterProfiler)
set.seed(10271996)


#Prepare data----
demographics <- read_excel("brain_tissue_demographics.xlsx")
df_full <- read_excel("JHOP-03-20MD CDT (1).xlsx", sheet = 2, col_names = F)
col_16mg <- which(df_full[8,] == 16) #Identify Sample with 16mg instead of 24mg
sample_16mg <- df_full[[1,col_16mg]]

df_subset <- df_full[ -c(2:10), c(2,14:ncol(df_full))]
df_subset_wide <- as.data.frame(t(df_subset))
colnames(df_subset_wide) <- df_subset_wide[1,]
colnames(df_subset_wide)[1:2]<- c("id", "group")
df_subset_wide <- df_subset_wide[-1,]

df_subset_wide[3:ncol(df_subset_wide)] <- lapply(df_subset_wide[3:ncol(df_subset_wide)], as.numeric)

df <- as.data.frame(t(df_subset_wide[,-2])) %>% row_to_names(row_number = 1)
df[1:ncol(df)] <- lapply(df[1:ncol(df)], as.numeric)
df[sample_16mg] <- df[sample_16mg]*24/16 

path <- df_full[11:nrow(df_full), 1:13] %>% row_to_names(row_number = 1) %>% clean_names() %>% mutate(metabolite = janitor::make_clean_names(biochemical))

## Remove Metabolites with NAs >= 30%
df <- df[rowMeans(is.na(df)) < 0.3, ]

#KNN imputation of missing values
df_imputed <- impute.knn(as.matrix(df), k = 3)
df_imputed <- as.data.frame(df_imputed$data)
df_original <- df_imputed
df_imputed[1:ncol(df_imputed)] <- lapply(df_imputed[1:ncol(df_imputed)], log) #log-transformation



df <- as.data.frame(t(df_imputed)) %>% rownames_to_column(var = "id")


df <- right_join(df_subset_wide[,1:2], df) %>% mutate(normals = group == "Control NAWM") %>% relocate(normals, .after = group)


#Fix names of samples
df <- df %>% mutate(id = str_replace_all(id, "\\.", "-"),
                    ms = as.factor(normals == F)) %>% relocate(ms, .after = normals)
df <- df %>% mutate(id = str_replace(id, "12-11-5", "12-78")) %>% 
        mutate(id = str_replace(id, "13-67", "13-47")) %>% 
        mutate(id = str_replace(id, "09-67", "9-67")) %>% 
        mutate(id = str_replace(id, "C-", "C"))

df <- df %>% mutate(subjects = ifelse(str_detect(id, " "), str_extract(id, "[:graph:]*(?= )"), id)) %>% relocate(subjects, .after = id) %>%
        mutate(subjects = str_replace(subjects, "-", "_"))
df <- df %>% mutate(group = ifelse(id %in% c("12-78", "13-47 PLA4 CA", "9-34 PLA1", "9-67 MRI13 CA", "13-47 PCA8 CA"), "ca_lesion",
                                   ifelse(id %in% c("13-15 PA 15-5", "13-47 PLA8 C1", "12-78 14-4", "13-15 MRI13", "9-34 MRI2"), "ci_lesion",
                                          ifelse(id %in% c("9-67 PLA61", "13-47 PLA4 cone", "13-47 PLA13-8 cone"), "core_lesion", 
                                                 ifelse(id %in% c("9-34 PLA1 NAWM", "13-15 MRI13 NAWM", "12-78 PLA11-5 NAWM", "9-67 MRI13 WM"), "nawm", "control"))))) %>%
        relocate(group, .after = ms) %>% left_join(demographics) %>% relocate(gender:pmi, .after = subjects)

#remove controls with no demographics
df <- df %>% filter(!(str_detect(id, "P WM")))

df$subjects <- factor(df$subjects)

#Scale metabolites
scale <- function(x){(x-mean(x))/sd(x)}
df[9:ncol(df)] <- lapply(df[9:ncol(df)], scale)




# WGCNA -------------------------------------------------------------------

datExpr0 <- df
datExpr0 <- datExpr0[, -c(2:8)] %>% column_to_rownames(var = "id")
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 50, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

allTraits = df[,c(1,2,6:8)]

traitRows = match(rownames(datExpr0), allTraits$id)
datTraits = allTraits[traitRows, -2];
rownames(datTraits) = allTraits[traitRows, 1];
datTraits <- datTraits %>% mutate(group = ifelse(normals == F, "MS", "HC")) %>% dplyr::select(-2)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = labels2colors(datTraits[,-c(1,2)]);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = "group",
                    main = "Sample dendrogram and trait heatmap")



# sample network based on squared Euclidean distance note that we
# transpose the data
A = adjacency(t(datExpr), type = "distance")
# this calculates the whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized connectivity
Z.k = scale(k)

# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -3
remove.samples = Z.k < thresholdZ.k | is.na(Z.k)
sum(remove.samples)



# Choose a set of soft-thresholding powers
powers = c(1:20)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#power of 16 for signed network
net = blockwiseModules(datExpr, power = 16, TOMType = "signed",
                       networkType = "signed", minModuleSize = 6,
                       corType = "bicor",
                       mergeCutHeight = 0.25,
                       numericLabels = T, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)
table(net$colors)
length(unique(net$colors))
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#moduleTraitCor = cor(MEs, datTraits, use = "p");
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# Define variable weight containing the weight column of datTrait
#weight = as.data.frame(datTraits$weight_g);
#names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
#geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
#names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = names(datExpr),
                       moduleColor = moduleColors)
geneInfo0 <- geneInfo0 %>% rename(biochemical = substanceBXH)
geneInfo0 <- right_join(path[,c(4,2)], geneInfo0)

# Order modules by their significance for weight
modOrder = order(1:ncol(MMPvalue));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                               MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}


#Arranged based on color and then module membership
all_colors <- c()
for(i in unique(geneInfo0$moduleColor)[order(unique(geneInfo0$moduleColor))]){
        var <- paste0("MM.", i)
        color <- geneInfo0 %>% filter(moduleColor == i) %>% 
                arrange(desc(abs(!!sym(var)))) %>% 
                mutate(MM_sign = ifelse(!!sym(var) < 0, "-", "+")) %>%
                relocate(MM_sign, .after = moduleColor)
        all_colors <- bind_rows(all_colors, color)
}

all_colors <- all_colors[, c(colnames(all_colors)[1:4], colnames(all_colors)[5:ncol(all_colors)][order(colnames(all_colors)[5:ncol(all_colors)])])]



#Only >=0.8 module membership metabolites
ordered_colors <- unique(geneInfo0$moduleColor)[order(unique(geneInfo0$moduleColor))]
all_colors_short <- c()
for(i in ordered_colors[which(ordered_colors != "grey")]){
        var <- paste0("MM.", i)
        color <- geneInfo0 %>% filter(moduleColor == i) %>% 
                filter(abs(!!sym(var)) >= 0.8) %>%
                arrange(desc(abs(!!sym(var)))) %>% 
                mutate(MM_sign = ifelse(!!sym(var) < 0, "-", "+")) %>%
                relocate(MM_sign, .after = moduleColor)
        all_colors_short <- bind_rows(all_colors_short, color)
}

all_colors_short <- all_colors_short[, c(colnames(all_colors_short)[1:4], colnames(all_colors_short)[5:ncol(all_colors_short)][order(colnames(all_colors_short)[5:ncol(all_colors_short)])])]

metabolomics_MEs <- MEs
metabolomics_MEs <- metabolomics_MEs %>% rownames_to_column(var = "id")

final <- MEs %>% rownames_to_column(var = "id")

final <- left_join(final, df[,c(1,2,3,4,7,8)])
final <- final %>% mutate(ms = ifelse(ms == T, 1, 0))
final$group <- factor(final$group, levels = c("control", "nawm", "ca_lesion" , "ci_lesion", "core_lesion"))
final <- final %>% mutate(group_2 = ifelse(group %in% c("ca_lesion", "ci_lesion", "core_lesion"), "ms_lesion", as.character(group)))
final$group_2 <- factor(final$group_2, levels = c("control", "nawm", "ms_lesion"))

final <- final %>% arrange(subjects)

module_names <- data.frame(ME = paste0("ME", sort(unique(all_colors$moduleColor))),
                           module = str_to_title(sort(unique(all_colors$moduleColor))),
                           module_composition = c("Nucleotide and Histidine metabolism", "Dipeptides", "Unsaturated Fatty Acids, Endocannabinoids & Lysophospholipids",
                                                  "Benzoate Metabolism & Food Component", "Phosphatidylethanolamines & Phosphatidylcholines", "Hexocylceramides & Plasmalogens",
                                                  NA, "Nucleotide Metabolism", "Energy Metabolites", "Nonspecific", "Lysophospholipids", "Hexosylceramides",
                                                  "Acyl Carnitines", "Sphingosines", "Monoacylglycerols", "Amino Acids Metabolism", "Sphingomyelins and Ceramides"))



all_colors_first5 <- all_colors_short %>% group_by(moduleColor) %>% slice_head(n = 5)
all_colors_first3 <- all_colors_short %>% group_by(moduleColor) %>% slice_head(n = 3) %>% ungroup() %>% mutate(ME = paste0("ME",moduleColor)) %>% left_join(module_names) %>% select(module, module_composition, biochemical, sub_pathway)

#write_csv(all_colors_first3, "top_3_metabolites_per_module.csv")



# GEE - 3groups - HC-Periplaque-MS -------------------------

gee_results_3_groups_adj <- c()
for(i in colnames(final)[str_detect(colnames(final), "ME")]){
        formula <- as.formula(paste0(i, "~ group_2 + gender + age"))
        gee_model <- geeglm(formula, data = final, id = subjects, corstr = "unstructured")
        coef <- summary(gee_model)$coefficients[2:3,]
        coef <- cbind(i, coef) %>% mutate(lwr = Estimate - 1.96*Std.err,
                                          upr = Estimate + 1.96*Std.err)
        gee_results_3_groups_adj <- rbind(gee_results_3_groups_adj, coef)
}

gee_results_3_groups_adj <- gee_results_3_groups_adj %>% rename(p_value = `Pr(>|W|)`,
                                                                module = i) %>%
        rownames_to_column(var = "comparison")  %>%
        mutate(comparison = str_replace_all(comparison, "group_2", ""))

gee_results_3_groups_adj <- gee_results_3_groups_adj %>% mutate(sign = ifelse(p_value < 0.05, "*", "ns"),
                                                                p_adj = p.adjust(gee_results_3_groups_adj$p_value, method = "BH"),
                                                                p_adj_sign = ifelse(p_adj < 0.05, "*", "ns"),
                                                                direction = ifelse(Estimate < 0 & p_value < 0.05, intToUtf8(8595), 
                                                                                   ifelse(Estimate > 0 & p_value < 0.05, intToUtf8(8593), "ns")),
                                                                moduleColor = str_extract(module, "(?<=ME)[:graph:]*"),
                                                                estimate_ci = paste0(round(Estimate, 2), " (", round(lwr,2), " to ", round(upr,2), ")")) %>%
        arrange(moduleColor, comparison)

gee_results_ms_lesions_adj <- gee_results_3_groups_adj %>% filter(str_detect(comparison, "ms_lesion"))
gee_results_3_groups_adj <- gee_results_3_groups_adj %>% mutate(comparison = str_replace_all(comparison, "[:digit:]", ""))

module_size <- as.data.frame(table(all_colors$moduleColor)) %>% rename(moduleColor = Var1)
gee_results_3_groups_adj <- gee_results_3_groups_adj %>% left_join(module_size) %>% arrange(comparison)


gee_results_table<- gee_results_3_groups_adj %>% mutate(p_value = ifelse(p_value < 0.001, p_value, 
                                                                         ifelse(p_value >= 0.001 & p_value < 0.01, round(p_value, 3), round(p_value,2))),
                                                        p_adj = ifelse(p_adj < 0.001, p_adj, 
                                                                       ifelse(p_adj >= 0.001 & p_adj < 0.01, round(p_adj, 3), round(p_adj,2)))) %>%
        select(comparison, module, Freq, estimate_ci, p_value, p_adj) %>%  pivot_wider(names_from = comparison, values_from = c(estimate_ci, p_value, p_adj), names_vary = "slowest")

sort(unique((gee_results_3_groups_adj %>% filter(sign == "*"))$moduleColor))
sort(unique((gee_results_3_groups_adj %>% filter(p_adj_sign == "*"))$moduleColor))


gee_results_table <- gee_results_table %>% mutate(module = str_replace_all(module, "ME", ""),
                                                  module = str_to_title(module)) %>%
        left_join(module_names) %>% relocate(module_composition, .after = module)

all_colors_except_grey <- sort(unique(module_names$module))[unique(module_names$module) != "Grey"]
colors_order <- c(all_colors_except_grey,"Grey")
all_colors <- all_colors %>% left_join(module_names %>% mutate(moduleColor = str_extract(ME, "(?<=ME)[:alnum:]+")) %>% select(-ME)) %>% relocate(c(module,module_composition), .before = 1) %>% arrange(factor(module, levels = colors_order)) %>%
        select(-c(MM_sign, moduleColor)) %>% rename(metabolite = biochemical)
# write_csv(all_colors, "all_modules_with_metabolites.csv")



# Supplementary table 2 ----
# write_csv(gee_results_table, "gee_results_cntrl_periplaque_lesions_unstructured.csv")


#ANOVA 5 groups----
aov_results_5_groups <- c()
for(i in colnames(final)[str_detect(colnames(final), "ME")]){
        formula <- as.formula(paste0(i, "~ group"))
        aov_model <- aov(formula, data = final)
        tukey <- cbind(as.data.frame(TukeyHSD(aov_model)$group), i) %>% clean_names()
        aov_results_5_groups <- rbind(aov_results_5_groups, tukey)
}

aov_results_5_groups <- aov_results_5_groups %>% 
        rownames_to_column(var = "comparison") %>%
        mutate(sign = ifelse(p_adj < 0.05 & p_adj >= 0.01, "*",
                             ifelse(p_adj < 0.01 & p_adj >=0.001, "**", 
                                    ifelse(p_adj < 0.001, "***", "ns"))))

aov_results_5_groups <- aov_results_5_groups %>% mutate(dif_ci = paste0(round(diff, 2), " (", round(lwr, 2), " to ", round(upr, 2), ")"))
aov_results_5_groups <- aov_results_5_groups %>% rename(ME = i) %>% left_join(module_names)
aov_results_5_groups <- aov_results_5_groups %>% mutate(comparison = str_replace_all(comparison, "[:digit:]", ""))
aov_results_5_groups <- aov_results_5_groups %>% mutate(group_1 = str_extract(comparison, "[:graph:]+(?=-)"),
                                                        group_2 = str_extract(comparison, "(?<=-)[:graph:]+")) %>% 
        relocate(sign, .after = group_2)

aov_results_5_groups <- aov_results_5_groups %>% mutate(module = ifelse(ME == "MEgrey", "grey", module)) %>% select(ME, comparison, group_1, group_2, dif_ci, p_adj, sign) %>% left_join(module_names) %>% relocate(module:module_composition, .before = 1) %>% select(-ME)
aov_results_5_groups_sign <- aov_results_5_groups %>% filter(p_adj <0.05)
sort(unique(aov_results_5_groups_sign$ME))

#Supplementary table 3----
# write_csv(aov_results_5_groups, "tukey_all_groups.csv") 



#snRRNA seq data - Cell analysis----
df_names <- df %>% left_join(final[c("id","group_2")]) %>% clean_names() %>% arrange(subjects)


names_old <- read_excel("snRNAseq_for_metabolomics.xlsx", sheet = 1, range = "A1:C21") %>% clean_names()
metabolon_names <- read_excel("tissue_metabolon.xlsx", sheet = 2, range = "N1:AR6", col_names = F)
metabolon_names <- as.data.frame(t(metabolon_names[c(1,6), ]))
colnames(metabolon_names) <- c("nbb_number", "id")
metabolon_names <- metabolon_names %>% mutate(id = str_replace_all(id, "\\.", "-"),
                                              id = str_replace(id, "12-11-5", "12-78"),
                                              id = str_replace(id, "13-67", "13-47"),
                                              id = str_replace(id, "09-67", "9-67"),
                                              id = str_replace(id, "C-", "C")) %>%
        filter(!(str_detect(id, "P WM")))

names_old <- names_old %>% left_join(metabolon_names)
sum(names_old$id %in% final$id)

##Cell count ----
cells <- read_excel("snRNAseq_for_metabolomics.xlsx", sheet = 6, col_names = F)
cell_types <- cells[ ,1:2] %>% row_to_names(row_number = 1) %>% clean_names() %>% rename(cell = sample_number) %>% mutate(cell = janitor::make_clean_names(cell))
cells <- as.data.frame(t(cells[-2])) %>% row_to_names(row_number = 1) %>% clean_names()
get_numbers <- function(x) {as.numeric(str_extract(x, "[:digit:]+"))}
cells[] <- lapply(cells, get_numbers)
cells <- cells %>% rowwise() %>% mutate(dc = dc_1 + dc_2,
                                        mims = mims_foamy + mims_iron,
                                        imm = sum(c_across(homeostatic_microglia:plasmablasts)),
                                        opc = sum(c_across(homeostatic_cl0_2:odc_cl6)),
                                        astro = sum(c_across(non_reactive_cl0:microtubules_cl4)),
                                        oligo = sum(c_across(c(ol_0:ol_4, ol_6, ol_7))),
                                        stressed_ol = ol_6 + ol_7,
                                        oligo_not_stressed = sum(c_across(ol_0:ol_4)),
                                        astro_non_reactive = sum(c_across(non_reactive_cl0:non_reactive_cl3)))

cells <- cells %>% left_join(names_old) %>% relocate(nbb_number:id, .before = 1)

cells_modules <- final %>% right_join(cells)

cells_to_test <- colnames(cells)[-c(1:4,6,7,11,13,22:24, 29:36,39:42)]

cells <- cells %>% rowwise %>% mutate(total = sum(c_across(all_of(cells_to_test))))

cells_ratio <- cells %>% mutate(across(5:46, ~ .x/total))

cells_metabolite <- df_names %>% right_join(cells) 
cell_ratio_metabolite <- df_names %>% right_join(cells_ratio) 
sum(cells_to_test %in% colnames(cell_ratio_metabolite))

## Spearman individual metabolites----
cell_ratio_metabolite <- cell_ratio_metabolite %>% filter(sample_number != 6)
metabolite_results_spearman <- c()
for(i in colnames(cell_ratio_metabolite)[9:684]){
        for(j in cells_to_test) {
                spearman <- cor.test(cell_ratio_metabolite[[i]], cell_ratio_metabolite[[j]], method = "spearman")
                est <- as.data.frame(cbind(i, j, est = as.numeric(spearman$estimate), p_value = as.numeric(spearman$p.value)))
                metabolite_results_spearman <- rbind(metabolite_results_spearman, est)
        }
}
metabolite_results_spearman[3:4] <- lapply(metabolite_results_spearman[3:4], as.numeric)



metabolite_results_spearman <- metabolite_results_spearman %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>%
        rename(metabolite = i, 
               cell = j) %>%
        left_join(cell_types) %>%
        left_join(path[c("sub_pathway", "biochemical", "metabolite")]) %>% 
        relocate(sub_pathway:biochemical, .before = 1) %>%
        relocate(type, .before = cell)
metabolite_results_spearman <- metabolite_results_spearman %>% mutate(type = ifelse(cell %in% c("dc", "mims"), "imm",
                                                                                    ifelse(cell %in% c("stressed_ol", "oligo_not_stressed"), "ol",
                                                                                           ifelse(cell == "astro_non_reactive", "astro", type))))


##MSEA----
pathways <- split(path$biochemical, path$sub_pathway) #Metabolon pathways with metabolites

set.seed(10271996)
msea_results <- c()
for(i in unique(metabolite_results_spearman$cell)) {
        tryCatch({
                df_temp <- metabolite_results_spearman %>% filter(cell == i)
                ranks <- df_temp$est
                names(ranks) <- df_temp$biochemical
                ranks <- ranks[order(ranks)]
                msea_spearman  <- fgsea(pathways = pathways, 
                                        stats    = ranks,
                                        minSize  = 4,
                                        maxSize  = 500,
                                        nPermSimple  = 10000)
                msea_spearman <- cbind(cell = i,msea_spearman)
                msea_results <- rbind(msea_results, msea_spearman)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

msea_results <- msea_results %>% mutate(qval = p.adjust(pval, method = "BH"))
msea_results <- msea_results %>% mutate(pval_color = ifelse(qval < 0.05 & NES < 0, 1,
                                                            ifelse(pval < 0.05 & NES < 0, 2,
                                                                   ifelse(pval > 0.05 & NES < 0, 3, 
                                                                          ifelse(qval < 0.05 & NES > 0, 6,
                                                                                 ifelse(pval < 0.05 & NES > 0, 5, 4))))))


msea_results <- msea_results %>% left_join((path %>% rename(pathway = sub_pathway) %>% distinct(pathway, .keep_all = T))[c("pathway", "super_pathway")])

msea_results <- msea_results %>% mutate(cell_clean = str_replace_all(cell, c("stressed_cl1" = "Stressed OPC",
                                                                             "premyelinating_cl7" = "Premyelinating OPC",
                                                                             "odc_cl6" = "Pre-oligo",
                                                                             "immune_like_cl3_5" = "Immune-like OPC",
                                                                             "homeostatic_cl0_2" = "Homeostatic OPC",
                                                                             "hmgb1_senescent_cl4" = "Senescent OPC",
                                                                             "stressed_ol" = "Stressed Oligo",
                                                                             "oligo_not_stressed" = "Oligo",
                                                                             "t_cells" = "T cells",
                                                                             "stressed_microglia" = "Stressed microglia",
                                                                             "plasmablasts" = "Plasmablasts",
                                                                             "mims" = "MIMS",
                                                                             "microglia_high_myein_genes" = "Imm (high-myelin RNA)",
                                                                             "macrophages" = "Macrophages",
                                                                             "homeostatic_microglia" = "Homeostatic microglia",
                                                                             "dc" = "Mono/moDC",
                                                                             "reactive_cl1_5" = "Reactive/Stressed astrocytes",
                                                                             "perinodal_cl7" = "Perinodal astrocytes",
                                                                             "microtubules_cl4" = "Senescent astrocytes",
                                                                             "inflammatory_cl6" = "Inflammatory astrocytes",
                                                                             "astro_non_reactive" = "Non-reactive astrocytes",
                                                                             "other" = "Other")))



#Figure 2 - Bubble Heat Map Plot----
cell_order <- unique((msea_results  %>% left_join(cell_types) %>%
                              mutate(type = ifelse(cell %in% c("dc", "mims"), "imm",
                                                   ifelse(cell %in% c("stressed_ol", "oligo_not_stressed"), "ol",
                                                          ifelse(cell == "astro_non_reactive", "astro", type)))) %>%
                              arrange(type, cell))$cell_clean)

cell_order <- cell_order[c(1:10,12:14,11,15:22)]




text_colors <- c(rep("#33286f", 7), rep("#7d161a", 7))
msea_results_backup <- msea_results
msea_results <- msea_results %>% mutate(pathway = ifelse(pathway == "Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)", "Polyunsaturated Acyl Carnitine",
                                                         ifelse(pathway == "Long Chain Polyunsaturated Fatty Acid (n3 and n6)", "Long Chain Polyunsaturated Fatty Acid", 
                                                                ifelse(pathway == "Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)", "Long Chain Saturated Acyl Carnitine",
                                                                       ifelse(pathway == "Hexosylceramides (HCER)", "Hexosylceramides",
                                                                              ifelse(pathway == "Phosphatidylcholine (PC)", "Phosphatidylcholines",
                                                                                     ifelse(pathway == "Phosphatidylethanolamine (PE)", "Phosphatidylethanolamines", pathway)))))))

msea_results <- msea_results %>% mutate(pathway = ifelse(str_detect(pathway, "Long Chain|glycerol|Endocan|phospholipid|alogen"), paste0(pathway,"s"), pathway))
pathway_order <- c("Pantothenate and CoA Metabolism", "Monoacylglycerols", "Lysophospholipids", "Long Chain Polyunsaturated Fatty Acids", "Long Chain Monounsaturated Fatty Acids", 
                   "Long Chain Saturated Acyl Carnitines", "Endocannabinoids", "Plasmalogens", "Phosphatidylethanolamines", "Phosphatidylcholines", "Diacylglycerols", "Sphingosines",
                   "Sphingomyelins", "Ceramides")
length(pathway_order); sum(unique(msea_results$pathway) %in% pathway_order)
pathway_order[!(pathway_order %in% unique(msea_results$pathway))]

figure_2 <- ggplot(msea_results %>% filter(pathway %in% pathway_order), aes(x = factor(cell_clean, level = cell_order), y = factor(pathway, level = pathway_order))) + 
        geom_point(aes(color = factor(pval_color), size = abs(NES)*5)) + 
        scale_color_manual(name = "Significance",
                           breaks = factor(c(3,4,2,5,1,6)),
                           labels = c("Decreasing p>0.05", "Increasing p>0.05", "Decreasing p<0.05", "Increasing p<0.05", "Decreasing q<0.05", "Increasing q<0.05"),
                           values = c("#bbb7e2", "#f39aa7", "#5055af", "#e22a34", "#33286f", "#7d161a")) +
        theme_bw() +
        scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 17),
                         position = "top") +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 38)) +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 0, size = 14),
              axis.text.y = element_text(color = text_colors, size = 14),
              axis.text = element_text(face = "bold"),
              legend.position = "bottom",
              axis.text.x.top = element_text(vjust = 0.5),
              legend.title = element_text(face = "bold", size = 13),
              legend.text = element_text(size = 12)) +
        guides(size = "none", color = guide_legend(override.aes = list(size = 4))) +
        scale_size_continuous(range = c(2,7)) +
        xlab(element_blank()) +
        ylab(element_blank()) + 
        geom_vline(xintercept = c(5.5, 14.5, 16.5), linetype = 2) +
        geom_hline(yintercept = 7.5) 

pdf("figure_2b.pdf", width = 17, height = 6.8)
figure_2
dev.off()




#WGCNA - Cell analysis----

#module data set with cell ratios
cell_ratio_modules <- final %>% right_join(cells_ratio) 
sum(cells_to_test %in% colnames(cell_ratio_metabolite))
cell_ratio_modules <- cell_ratio_modules %>% filter(sample_number != 6)
## Spearman using Module Scores----
module_results_spearman <- c()
for(i in colnames(cell_ratio_modules)[str_which(names(cell_ratio_modules), "^ME")]){
        for(j in cells_to_test) {
                spearman <- cor.test(cell_ratio_modules[[i]], cell_ratio_modules[[j]], method = "spearman")
                est <- as.data.frame(cbind(i, j, est = as.numeric(spearman$estimate), p_value = as.numeric(spearman$p.value)))
                module_results_spearman <- rbind(module_results_spearman, est)
        }
}
module_results_spearman[3:4] <- lapply(module_results_spearman[3:4], as.numeric)

module_results_spearman <- module_results_spearman %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>%
        rename(module = i, 
               cell = j) %>%
        left_join(cell_types) %>%
        mutate(type = ifelse(cell %in% c("dc", "mims"), "imm",
                             ifelse(cell %in% c("stressed_ol", "oligo_not_stressed"), "ol",
                                    ifelse(cell == "astro_non_reactive", "astro", type))))

module_results_spearman <- module_results_spearman %>% 
        mutate(module = str_to_title(str_extract(module, "(?<=ME)[:alpha:]+"))) %>%
        left_join(all_colors %>% distinct(module, .keep_all = T) %>% select(module, module_composition)) %>% mutate(title = paste0(module, " - ", module_composition))


module_results_spearman <- module_results_spearman %>% mutate(pval_color = ifelse(fdr < 0.05 & est < 0, 1,
                                                                                  ifelse(p_value < 0.05 & est < 0, 2,
                                                                                         ifelse(p_value > 0.05 & est < 0, 3, 
                                                                                                ifelse(fdr < 0.05 & est > 0, 6,
                                                                                                       ifelse(p_value < 0.05 & est > 0, 5, 4))))))

module_results_spearman<- module_results_spearman %>% mutate(cell_clean = str_replace_all(cell, c("stressed_cl1" = "Stressed OPC",
                                                                                                  "premyelinating_cl7" = "Premyelinating OPC",
                                                                                                  "odc_cl6" = "Pre-oligo",
                                                                                                  "immune_like_cl3_5" = "Immune-like OPC",
                                                                                                  "homeostatic_cl0_2" = "Homeostatic OPC",
                                                                                                  "hmgb1_senescent_cl4" = "Senescent OPC",
                                                                                                  "stressed_ol" = "Stressed Oligo",
                                                                                                  "oligo_not_stressed" = "Oligo",
                                                                                                  "t_cells" = "T cells",
                                                                                                  "stressed_microglia" = "Stressed microglia",
                                                                                                  "plasmablasts" = "Plasmablasts",
                                                                                                  "mims" = "MIMS",
                                                                                                  "microglia_high_myein_genes" = "Imm (high-myelin RNA)",
                                                                                                  "macrophages" = "Macrophages",
                                                                                                  "homeostatic_microglia" = "Homeostatic microglia",
                                                                                                  "dc" = "Mono/moDC",
                                                                                                  "reactive_cl1_5" = "Reactive/Stressed astrocytes",
                                                                                                  "perinodal_cl7" = "Perinodal astrocytes",
                                                                                                  "microtubules_cl4" = "Senescent astrocytes",
                                                                                                  "inflammatory_cl6" = "Inflammatory astrocytes",
                                                                                                  "astro_non_reactive" = "Non-reactive astrocytes",
                                                                                                  "other" = "Other")))
# pdf("figure_2_module_correlations.pdf", width = 17, height = 10)
# ggplot(module_results_spearman, aes(x = factor(cell_clean, level = cell_order), y = factor(title, level = pathway_module_order))) + 
#         geom_point(aes(color = factor(pval_color), size = abs(est)*5)) + 
#         scale_color_manual(name = "Significance",
#                            breaks = factor(c(3,4,2,5,1,6)),
#                            labels = c("Decreasing p>0.05", "Increasing p>0.05", "Decreasing p<0.05", "Increasing p<0.05", "Decreasing q<0.05", "Increasing q<0.05"),
#                            values = c("#b5bbf5", "#f2aaaa", "#4952f5", "#e03838", "#340e75", "#800101")) +
#         theme_bw() +
#         scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 17),
#                          position = "top") +
#         scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 38)) +
#         theme(panel.grid = element_blank(),
#               axis.text.x = element_text(angle = 90, hjust = 0, size = 14),
#               axis.text.y = element_text(size = 10),
#               axis.text = element_text(face = "bold"),
#               legend.position = "bottom",
#               axis.text.x.top = element_text(vjust = 0.5),
#               legend.title = element_text(face = "bold", size = 13),
#               legend.text = element_text(size = 12),
#               title = element_text(size = 15, face = "bold")) +
#         guides(size = "none", color = guide_legend(override.aes = list(size = 4))) +
#         scale_size_continuous(range = c(2,7)) +
#         xlab(element_blank()) +
#         ylab(element_blank()) + 
#         geom_vline(xintercept = c(5.5, 14.5, 16.5), linetype = 2) +
#         ggtitle("Spearman's Correlations - Module Scores and Cell Ratios")
# dev.off()



##MSEA using WGCNA pathways----
pathways_module <- split(all_colors$metabolite, all_colors$module) #Metabolon pathways with metabolites

set.seed(10271996)
msea_results_module <- c()
for(i in unique(metabolite_results_spearman$cell)) {
        tryCatch({
                df_temp <- metabolite_results_spearman %>% filter(cell == i)
                ranks <- df_temp$est
                names(ranks) <- df_temp$biochemical
                ranks <- ranks[order(ranks)]
                msea_spearman  <- fgsea(pathways = pathways_module, 
                                        stats    = ranks,
                                        minSize  = 4,
                                        maxSize  = 500,
                                        nPermSimple  = 10000)
                msea_spearman <- cbind(cell = i,msea_spearman)
                msea_results_module <- rbind(msea_results_module, msea_spearman)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

msea_results_module <- msea_results_module %>% mutate(qval = p.adjust(pval, method = "BH"))
msea_results_module <- msea_results_module %>% mutate(pval_color = ifelse(qval < 0.05 & NES < 0, 1,
                                                                          ifelse(pval < 0.05 & NES < 0, 2,
                                                                                 ifelse(pval > 0.05 & NES < 0, 3, 
                                                                                        ifelse(qval < 0.05 & NES > 0, 6,
                                                                                               ifelse(pval < 0.05 & NES > 0, 5, 4))))))

msea_results_module <- msea_results_module %>% mutate(cell_clean = str_replace_all(cell, c("stressed_cl1" = "Stressed OPC",
                                                                                           "premyelinating_cl7" = "Premyelinating OPC",
                                                                                           "odc_cl6" = "Pre-oligo",
                                                                                           "immune_like_cl3_5" = "Immune-like OPC",
                                                                                           "homeostatic_cl0_2" = "Homeostatic OPC",
                                                                                           "hmgb1_senescent_cl4" = "Senescent OPC",
                                                                                           "stressed_ol" = "Stressed Oligo",
                                                                                           "oligo_not_stressed" = "Oligo",
                                                                                           "t_cells" = "T cells",
                                                                                           "stressed_microglia" = "Stressed microglia",
                                                                                           "plasmablasts" = "Plasmablasts",
                                                                                           "mims" = "MIMS",
                                                                                           "microglia_high_myein_genes" = "Imm (high-myelin RNA)",
                                                                                           "macrophages" = "Macrophages",
                                                                                           "homeostatic_microglia" = "Homeostatic microglia",
                                                                                           "dc" = "Mono/moDC",
                                                                                           "reactive_cl1_5" = "Reactive/Stressed astrocytes",
                                                                                           "perinodal_cl7" = "Perinodal astrocytes",
                                                                                           "microtubules_cl4" = "Senescent astrocytes",
                                                                                           "inflammatory_cl6" = "Inflammatory astrocytes",
                                                                                           "astro_non_reactive" = "Non-reactive astrocytes",
                                                                                           "other" = "Other")))

msea_results_module <- msea_results_module %>% left_join(all_colors %>% distinct(module, .keep_all = T) %>% select(module, module_composition) %>% rename(pathway = module)) %>% mutate(title = paste0(pathway, " - ", module_composition))
pathway_module_order <- sort(unique(msea_results_module$title), decreasing = T)



pdf("figure_2_module_clusters_19_samples.pdf", width = 17, height = 10)
ggplot(msea_results_module , aes(x = factor(cell_clean, level = cell_order), y = factor(title, level = pathway_module_order))) +
        geom_point(aes(color = factor(pval_color), size = abs(NES)*5)) +
        scale_color_manual(name = "Significance",
                           breaks = factor(c(3,4,2,5,1,6)),
                           labels = c("Decreasing p>0.05", "Increasing p>0.05", "Decreasing p<0.05", "Increasing p<0.05", "Decreasing q<0.05", "Increasing q<0.05"),
                           values = c("#b5bbf5", "#f2aaaa", "#4952f5", "#e03838", "#340e75", "#800101")) +
        theme_bw() +
        scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 17),
                         position = "top") +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 38)) +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 0, size = 14),
              axis.text.y = element_text(size = 10),
              axis.text = element_text(face = "bold"),
              legend.position = "bottom",
              axis.text.x.top = element_text(vjust = 0.5),
              legend.title = element_text(face = "bold", size = 13),
              legend.text = element_text(size = 12),
              title = element_text(size = 15, face = "bold")) +
        guides(size = "none", color = guide_legend(override.aes = list(size = 4))) +
        scale_size_continuous(range = c(2,7)) +
        xlab(element_blank()) +
        ylab(element_blank()) +
        geom_vline(xintercept = c(5.5, 14.5, 16.5), linetype = 2) +
        ggtitle("MSEA with WGCNA pathways")
dev.off()




#Figure 1 ----------------------------------------------------

final$group <- factor(final$group, levels = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"))

brown <- ggplot(final, aes(x = factor(group), y = MEbrown, shape = factor(group))) + 
        stat_summary(fun = "median", fun.min = "median", fun.max = "median", geom = "crossbar",show.legend = F, alpha = 0.1, width = 0.5, size = 0.4, color = "grey20") +
        geom_beeswarm(aes(fill = factor(group), color = factor(group)), cex = 3.7, alpha = 0.7, size = 5) +
        geom_beeswarm(aes(color = factor(group)),cex = 3.7, size = 5, alpha = 1, stroke = 1) +
        scale_x_discrete(limits = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = "Pathological Stage",
                              breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = "Pathological Stage",
                             breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = "Pathological Stage",
                           breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab("Module Score") +
        theme_classic() +
        theme(legend.position = "bottom",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0))) +
        geom_bracket(y.position = 0.45, xmin = 1, xmax = 5, size = 0.5, label.size = 8,
                     label = "**", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.35, xmin = 2, xmax = 5, size = 0.5, label.size = 8,
                     label = "*", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        labs(title = "Unsaturated FAs, Endocannabinoids\n& Lysophospholipids") 



green <- ggplot(final, aes(x = factor(group), y = MEgreen, shape = factor(group))) + 
        stat_summary(fun = "median", fun.min = "median", fun.max = "median", geom = "crossbar",show.legend = F, alpha = 0.1, width = 0.5, size = 0.4, color = "grey20") +
        geom_beeswarm(aes(fill = factor(group), color = factor(group)), cex = 3.7, alpha = 0.7, size = 5) +
        geom_beeswarm(aes(color = factor(group)),cex = 3.7, size = 5, alpha = 1, stroke = 1) +
        scale_x_discrete(limits = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = "Pathological Stage",
                              breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = "Pathological Stage",
                             breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = "Pathological Stage",
                           breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab(element_blank()) +
        theme_classic() +
        theme(legend.position = "bottom",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0))) +
        geom_bracket(y.position = 0.45, xmin = 1, xmax = 5, size = 0.5, label.size = 8,
                     label = "*", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        labs(title = "Dipeptides") 


lightcyan <- ggplot(final, aes(x = factor(group), y = MElightcyan, shape = factor(group))) + 
        stat_summary(fun = "median", fun.min = "median", fun.max = "median", geom = "crossbar",show.legend = F, alpha = 0.1, width = 0.5, size = 0.4, color = "grey20") +
        geom_beeswarm(aes(fill = factor(group), color = factor(group)), cex = 3.7, alpha = 0.7, size = 5) +
        geom_beeswarm(aes(color = factor(group)),cex = 3.7, size = 5, alpha = 1, stroke = 1) +
        scale_x_discrete(limits = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = "Pathological Stage",
                              breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = "Pathological Stage",
                             breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = "Pathological Stage",
                           breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab(element_blank()) +
        theme_classic() +
        theme(legend.position = "bottom",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0))) +
        geom_bracket(y.position = 0.55, xmin = 1, xmax = 5, size = 0.5, label.size = 8,
                     label = "*", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        labs(title = "Nucleotide Metabolism") 


pink <- ggplot(final, aes(x = factor(group), y = MEpink, shape = factor(group))) + 
        stat_summary(fun = "median", fun.min = "median", fun.max = "median", geom = "crossbar",show.legend = F, alpha = 0.1, width = 0.5, size = 0.4, color = "grey20") +
        geom_beeswarm(aes(fill = factor(group), color = factor(group)), cex = 3.7, alpha = 0.7, size = 5) +
        geom_beeswarm(aes(color = factor(group)),cex = 3.7, size = 5, alpha = 1, stroke = 1) +
        scale_x_discrete(limits = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = "Pathological Stage",
                              breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = "Pathological Stage",
                             breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = "Pathological Stage",
                           breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab(element_blank()) +
        theme_classic() +
        theme(legend.position = "bottom",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0))) +
        geom_bracket(y.position = 0.7, xmin = 1, xmax = 5, size = 0.5, label.size = 8,
                     label = "***", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.6, xmin = 1, xmax = 4, size = 0.5, label.size = 8,
                     label = "**", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.5, xmin = 1, xmax = 3, size = 0.5, label.size = 8,
                     label = "**", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.4, xmin = 1, xmax = 2, size = 0.5, label.size = 8,
                     label = "*", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        labs(title = "Lysophospholipids") 






salmon <- ggplot(final, aes(x = factor(group), y = MEsalmon, shape = factor(group))) + 
        stat_summary(fun = "median", fun.min = "median", fun.max = "median", geom = "crossbar",show.legend = F, alpha = 0.1, width = 0.5, size = 0.4, color = "grey20") +
        geom_beeswarm(aes(fill = factor(group), color = factor(group)), cex = 3.7, alpha = 0.7, size = 5) +
        geom_beeswarm(aes(color = factor(group)),cex = 3.7, size = 5, alpha = 1, stroke = 1) +
        scale_x_discrete(limits = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = "Pathological Stage",
                              breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = "Pathological Stage",
                             breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = "Pathological Stage",
                           breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab(element_blank()) +
        theme_classic() +
        theme(legend.position = "bottom",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0))) +
        geom_bracket(y.position = 0.55, xmin = 1, xmax = 5, size = 0.5, label.size = 8,
                     label = "**", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.45, xmin = 1, xmax = 3, size = 0.5, label.size = 8,
                     label = "*",inherit.aes = F, tip.length = 0, vjust = 0.5) +
        labs(title = "Sphingosines")


yellow <- ggplot(final, aes(x = factor(group), y = MEyellow, shape = factor(group))) + 
        stat_summary(fun = "median", fun.min = "median", fun.max = "median", geom = "crossbar",show.legend = F, alpha = 0.1, width = 0.5, size = 0.4, color = "grey20") +
        geom_beeswarm(aes(fill = factor(group), color = factor(group)), cex = 3.7, alpha = 0.7, size = 5) +
        geom_beeswarm(aes(color = factor(group)),cex = 3.7, size = 5, alpha = 1, stroke = 1) +
        scale_x_discrete(limits = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = "Pathological Stage",
                              breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = "Pathological Stage",
                             breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = "Pathological Stage",
                           breaks = c("control", "nawm", "ci_lesion", "ca_lesion", "core_lesion"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab(element_text("Module Score")) +
        theme_classic() +
        theme(legend.position = "bottom",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0))) +
        geom_bracket(y.position = 0.85, xmin = 1, xmax = 5, size = 0.5, label.size = 8,
                     label = "***", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.75, xmin = 1, xmax = 3.95, size = 0.5, label.size = 8,
                     label = "***", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.75, xmin = 4.05, xmax = 5, size = 0.5, label.size = 8,
                     label = "***", inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.65, xmin = 1, xmax = 2.95, size = 0.5, label.size = 8,
                     label = "***",inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.65, xmin = 3.05, xmax = 5, size = 0.5, label.size = 8,
                     label = "***",inherit.aes = F, tip.length = 0, vjust = 0.5) +
        geom_bracket(y.position = 0.55, xmin = 2, xmax = 5, size = 0.5, label.size = 8,
                     label = "***",inherit.aes = F, tip.length = 0, vjust = 0.5) +
        labs(title = "Sphingomyelins & Ceramides")


figure_1_palette_F <- (yellow | salmon | green) / plot_spacer() / (brown + pink + lightcyan) + 
        plot_layout(heights = c(10,1,10), guides = 'collect') &
        theme(legend.position= "right",
              legend.direction = "vertical",
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20, face = "bold"),
              axis.title.y = element_text(size = 16, face = "bold"),
              axis.text = element_text(size = 15),
              axis.text.x = element_text(face = "bold"),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 15, face = "bold"))

pdf("figure_1_F_legend.pdf", width = 23, height = 14)
figure_1_palette_F
dev.off()




#Supplementary color map - ALL ----
path_alphabetical_order <- order(unique(msea_results$pathway), decreasing = T)
head(unique(msea_results$pathway)[path_alphabetical_order])

pdf("supplementary_figure_color_map.pdf", width = 15.75, height = 15.75)
ggplot(msea_results, aes(x = factor(cell_clean, level = cell_order), y = factor(pathway, level = unique(msea_results$pathway)[path_alphabetical_order]))) +
        geom_point(aes(color = factor(pval_color), size = abs(NES))) +
        scale_color_manual(name = "Significance",
                           breaks = factor(c(3,4,2,5,1,6)),
                           labels = c("Decreasing p>0.05", "Increasing p>0.05", "Decreasing p<0.05", "Increasing p<0.05", "Decreasing q<0.05", "Increasing q<0.05"),
                           values = c("#b5bbf5", "#f2aaaa", "#4952f5", "#e03838", "#340e75", "#800101")) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
              legend.position = "bottom") +
        scale_x_discrete(position = "top") +
        guides(size = "none", color = guide_legend(override.aes = list(size = 4))) +
        xlab(element_blank()) +
        ylab(element_blank()) +
        geom_vline(xintercept = c(5.5, 14.5, 16.5), linetype = 2)
dev.off()

#dummy plot_for_biorender ----
bio_render <- ggplot(msea_results %>% filter(pathway %in% sample(unique(msea_results$pathway), 12)) %>% filter(cell_clean %in% sample(unique(msea_results$cell_clean), 16)), aes(x = cell_clean, y = factor(pathway))) + 
        geom_point(aes(color = factor(pval_color), size = abs(NES)*1000)) + 
        scale_color_manual(name = "Significance",
                           breaks = factor(c(3,4,2,5,1,6)),
                           labels = c("Decreasing p>0.05", "Increasing p>0.05", "Decreasing p<0.05", "Increasing p<0.05", "Decreasing q<0.05", "Increasing q<0.05"),
                           values = c("#b5bbf5", "#f2aaaa", "#4952f5", "#e03838", "#340e75", "#800101")) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              legend.position = "none",
              title = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=2)) +
        xlab(element_blank()) +
        ylab(element_blank()) 

#Supplementary table 4----
msea_results_2 <- msea_results
msea_results_2$cell_clean <- factor(msea_results_2$cell_clean, levels = cell_order)
msea_results_2 <- msea_results_2 %>% arrange(cell_clean, pathway)
msea_results_2$NES <- round(msea_results_2$NES, 2)
msea_table <- msea_results_2 %>% select(pathway, cell_clean, NES, pval, qval) %>% pivot_wider(names_from = cell_clean, values_from = c(NES, pval, qval), names_vary = "slowest")
# write_csv(msea_table, "supplementary_table_4_msea_results.csv")


#Supplememtary table 1 - GEE Individual Metabolites  ------------------------------------------------------------
gee_single_met_3_groups <- c()
for(i in colnames(df_names)[9:684]){
        gee_formula <- as.formula(paste0(i, "~ group_2 + gender + age"))
        gee_model <- geeglm(gee_formula, data = df_names, id = subjects, corstr = "unstructured")
        coef <- summary(gee_model)$coefficients[2:3,]
        coef <- cbind(i, coef) %>% mutate(lwr = Estimate - 1.96*Std.err,
                                          upr = Estimate + 1.96*Std.err)
        gee_single_met_3_groups <- rbind(gee_single_met_3_groups, coef)
}
#check if the same
gee_single_met_3_groups <- gee_single_met_3_groups %>% rename(p_value = `Pr(>|W|)`) %>%
        rownames_to_column(var = "comparison") #%>%
#        mutate(lwr = Estimate-1.96*Std.err,
#               upr=Estimate+1.96*Std.err)
gee_single_met_3_groups <- gee_single_met_3_groups %>% 
        mutate(significance = ifelse(p_value < 0.05, "*", "ns"),
               FDR_p_value = p.adjust(gee_single_met_3_groups$p_value, method = "BH"),
               FDR_significance = ifelse(FDR_p_value < 0.05, "*", "ns"))

gee_single_met_3_groups <- gee_single_met_3_groups %>% rename(metabolite = i) %>%
        left_join(path[c("sub_pathway", "biochemical", "metabolite")]) %>%
        relocate(sub_pathway:biochemical, .before = 1)

gee_single_met_3_groups <- gee_single_met_3_groups %>% 
        mutate(estimate_ci = paste0(round(Estimate, 2), " (", round(lwr,2), " to ", round(upr,2), ")"))


gee_single_met_3_groups <- gee_single_met_3_groups %>% mutate(sign = ifelse(p_value < 0.05, "*", "ns"),
                                                              p_adj = p.adjust(gee_single_met_3_groups$p_value, method = "BH"),
                                                              p_adj_sign = ifelse(p_adj < 0.05, "*", "ns"),
                                                              direction = ifelse(Estimate < 0 & p_value < 0.05, intToUtf8(8595), 
                                                                                 ifelse(Estimate > 0 & p_value < 0.05, intToUtf8(8593), "ns")),
                                                              estimate_ci = paste0(round(Estimate, 2), " (", round(lwr,2), " to ", round(upr,2), ")")) 


gee_single_met_3_groups <- gee_single_met_3_groups %>% mutate(comparison = str_replace_all(comparison, "[:digit:]", ""),
                                                              comparison = str_replace(comparison, "group_", "")) %>% 
        arrange(comparison, p_value)



gee_single_met_results_table<- gee_single_met_3_groups %>% select(comparison, sub_pathway, biochemical, estimate_ci, p_value, p_adj) %>%  pivot_wider(names_from = comparison, values_from = c(estimate_ci, p_value, p_adj), names_vary = "slowest") %>%
        rename(metabolite = biochemical)


# write_csv(gee_single_met_results_table, "supplementary_table_1_gee_single_met_results_table.csv")




