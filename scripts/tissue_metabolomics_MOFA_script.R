# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "C:/Users/dladaki1/Miniconda3/envs/mofa_env/Scripts/python.exe")
library(reticulate)
reticulate::use_condaenv(condaenv = "mofa_env")
reticulate::use_python("C:/Users/dladaki1/Miniconda3/envs/mofa_env/Scripts/python.exe")
py_config()

# 1Introduction -----------------------------------------------------------
# This vignette show how to use MOFA+ on the bulk multi-omics data set that was used in the first publication of MOFA and the original vignette of the MOFA package.

# Briefly, the data consists of four omics including DNA methylation, RNA-seq, somatic mutations and drug response data from blood for N=200 patients with Chronic Lymphocytic Leukemia (CLL). The data set is explained in detail in this article and is publicly available here

# 2Load libraries and data ------------------------------------------------
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(msigdbr)

# read in the data --------------------------------------------------------
# read in the LUT of the samples
# save the table of the metadata
LUT_sample <- read_tsv("out/table/LUT_sample_metabolite_RNA.tsv") %>% 
  mutate(id_fix = str_replace_all(id_fix,pattern = "-",replacement = ".")) %>% 
  # add the clumn sample to avoid issue
  mutate(sample = sample_number_fix) %>% 
  mutate(pathological_stage = factor(pathological_stage)) %>% 
  # recode the order of the samples by degree
  # notice that this will affect the correlation estimate with the metadata
  mutate(pathological_stage = factor(pathological_stage,levels = c("Control WM","Periplaque","CI","CA","core"))) 

#Remove sample no 6 as it was not the same tissue for transcriptomics and metabolomics
id_6 <- (LUT_sample %>% filter(sample_number == 6))$id_fix
LUT_sample <- LUT_sample %>% filter(sample_number != 6)

# read in the scaled metabolomic dataset
# boxplot(read_tsv("out/table/df_metabolites_scale_adj_top.tsv")[,-1])
# df_metabolic2 <- read_tsv("out/table/df_metabolites_scale_adj_top_featureScale.tsv")
# boxplot(df_metabolic2[,-1])
# load the recommended processing for the metabolomics
df_metabolic <- read_tsv("out/table/df_metabolites_fix_adj_top_featureScale.tsv")

#Remove sample no 6
df_metabolic <- df_metabolic %>% select(-id_6)

boxplot(df_metabolic[,-1])
# change the order of the samples and change the name of the column
df_metabolic_fix <- df_metabolic %>%
  column_to_rownames() %>% 
  # change the order of the samples
  dplyr::select(LUT_sample$id_fix) %>% 
  # rename the samples
  setNames(LUT_sample$sample_number_fix) %>% 
  as.matrix()

# read in the scaled transcriptomic dataset I only kept the top 5000
# boxplot(read_tsv("out/table/mRNA_matrix_counts_Sample_norm_scaled.tsv")[,-1])
# boxplot(read_tsv("out/table/df_rna_scale_adj_top.tsv")[,-1])
# this is the feature scaling that Giulia si recommending, done on the CPM matrix
# df_rna <- read_tsv("out/table/df_rna_scale_adj_top_featureScale_Giulia_CPM.tsv")
# boxplot(df_rna[,-1])
# this is the features scaling from the VST transformed data
df_rna <- read_tsv("out/table/df_rna_scale_adj_top_featureScale.tsv")

#Remove sample no 6
df_rna <- df_rna %>% select(-s6)

boxplot(df_rna[,-1])

# change the order of the samples and change the name of the column
df_rna_fix <- df_rna %>%
  column_to_rownames() %>% 
  # change the order of the samples
  dplyr::select(LUT_sample$sample_number_fix) %>% 
  as.matrix()

# build the list of layers ------------------------------------------------
# put the layers in a list
list_layers <- list(mRNA = df_rna_fix,metabolic = df_metabolic_fix)
str(list_layers)
list_layers$mRNA[1:10,1:10]

# metadata of the dataset
LUT_sample

# 3Create the MOFA obejct and train the model -----------------------------
# Create the MOFA object
MOFAobject <- create_mofa(list_layers)
MOFAobject

# 3.1Plot data overview ---------------------------------------------------
# Visualise the number of views (rows) and the number of groups (columns) exist, what are their corresponding dimensionalities and how many missing information they have (grey bars).
plot_data_overview(MOFAobject)

# 3.2.1Data options -------------------------------------------------------
# Important arguments:
# scale_groups: scale groups to the same total variance? Default is FALSE
# scale_views: scale views to the same total variance? Default is FALSE
# views: views names
# groups: groups names
data_opts <- get_default_data_options(MOFAobject)
data_opts

# 3.2.2Model options ------------------------------------------------------
# Important arguments:
# num_factors: number of factors
# likelihoods: likelihood per view (options are “gaussian”, “poisson”, “bernoulli”). By default they are inferred automatically.
# spikeslab_factors: use spike-slab sparsity prior in the factors? default is FALSE.
# spikeslab_weights: use spike-slab sparsity prior in the weights? default is TRUE.
# ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
# ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.
model_opts <- get_default_model_options(MOFAobject)
# ideal number of factor is 5.
model_opts$num_factors <- 4

model_opts

# 3.2.3Training options ---------------------------------------------------
# Important arguments:
# maxiter: number of iterations. Default is 1000.
# convergence_mode: “fast”, “medium” (default), “slow”. For exploration, the fast mode is good enough.
# seed: random seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 44

train_opts

# 3.3Train the MOFA model -------------------------------------------------
# Prepare the MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)
# Train the model: this should take ~2min

# NOTE: The software has evolved since the original publication and the results will not be 100% identical to the original publication, please use the pretrained model if you are running through the vignette for the fist time
MOFAobject <- run_mofa(MOFAobject, outfile="out/object/MOFA2_Sample_adj_top_featureScale.hdf5")
saveRDS(MOFAobject,"out/object/MOFA2_Sample_adj_top_featureScale.rds")
MOFAobject <- readRDS("out/object/MOFA2_Sample_adj_top_featureScale.rds")
# MOFAobject <- readRDS("out/object/MOFA2_CLL.rds")

# 4Overview of the trained MOFA model -------------------------------------

# 4.1Slots ----------------------------------------------------------------
# The MOFA object consists of multiple slots where relevant data and information is stored. For descriptions, you can read the documentation using ?MOFA. The most important slots are:
# data: input data used to train the model (features are centered at zero mean)
# samples_metadata: sample metadata information
# expectations: expectations of the posterior distributions for the Weights and the Factors
slotNames(MOFAobject)

# Data:
names(MOFAobject@data)
dim(MOFAobject@data$mRNA$group1)
dim(MOFAobject@data$metabolic$group1)

# Factor and Weight values (expectations of the posterior distributions):
names(MOFAobject@expectations)

# Dimensionality of the factor matrix: 200 samples, 15 factors
dim(MOFAobject@expectations$Z$group1)

# Dimensionality of the mRNA Weight matrix: 5000 features, 15 factors
dim(MOFAobject@expectations$W$mRNA)

# 4.2Add sample metadata to the model -------------------------------------
# The sample metadata must be provided as a data.frame and it must contain a column sample with the sample IDs. Make sure that the samples in the metadata match the samples in the model

# Sanity check
stopifnot(all(sort(LUT_sample$sample)==sort(unlist(samples_names(MOFAobject)))))

# Add sample metadata to the model
samples_metadata(MOFAobject) <- LUT_sample
samples_metadata(MOFAobject)

# 4.3Correlation between factors ------------------------------------------
# A good sanity check is to verify that the Factors are largely uncorrelated. In MOFA there are no orthogonality constraints such as in Principal Component Analysis, but if there is a lot of correlation between Factors this suggests a poor model fit. Reasons? Perhaps you used too many factors or perhaps the normalisation is not adequate.
pdf("out/image/01_MOFA2_Sample_adj_top_featureScale_CorFactor.pdf",width = 5,height = 4)
plot_factor_cor(MOFAobject)
dev.off()

# 4.4Plot variance decomposition ------------------------------------------

# 4.4.1Variance decomposition by Factor -----------------------------------
# The most important insight that MOFA generates is the variance decomposition analysis. This plot shows the percentage of variance explained by each factor across each data modality (and group, if provided). It summarises the sources of variation from a complex heterogeneous data set in a single figure.
plot_variance_explained(MOFAobject, max_r2=NULL)
ggsave("out/image/02_MOFA2_Sample_adj_top_featureScale_VarDecom.pdf",width = 4,height = 6)

# What insights from the data can we learn just from inspecting this plot?

# Factor 1,2,3 captures a source of variability that is present across all data modalities. Thus, its etiology is likely to be something very important for the disease
# Factor 4 captures a source of variation that is exclusive to the mRNA data.
# Factor 6 captures a source of variation that is exclusive to the metabolic data.

# Factor 1 captures a source of variability that is present across all data modalities. Thus, its etiology is likely to be something very important for the disease
# Factor 2 captures a very strong source of variation that is exclusive to the drug response data.
# Factor 3 captures variation that is present across multiple data modalities, except for DNA methylation. This is likely to be important too.
# Factor 5 is capturing some co-variation between the mRNA and the drug response assay.

# 4.4.2Total variance explained per view ----------------------------------
# A reasonable question is whether the model is providing a good fit to the data. For this we can plot the total variance explained (using all factors). The resulting values will depend on the nature of the data set, the number of samples, the number of factors, etc. Some general guidelines:
# Noisy data sets with strong non-linearities will result in small amounts of variance explained (<10%).
# The higher the number of samples the smaller the total variance explained
# The higher the number of factors, the higher the total variance explained.
# MOFA is a linear and sparse model. This is helpful to prevent overfitting, but it will never explain 100% of the variance, even if using a lot of Factors.
# In this data set, using only K=15 factors the model explains up to ~54% of the variation in the Drug response and ~42% of the variation in the mRNA data. This is quite remarkable for a linear model.
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
ggsave("out/image/03_MOFA2_Sample_adj_top_featureScale_VarLayer.pdf",width = 4,height = 4)

calculate_variance_explained(MOFAobject)
# 5Characterisation of Factor 1 -------------------------------------------
# There are a few systematic strategies to characterise the molecular etiology underlying the MOFA Factors and to relate them to the sample covariates:
# Association analysis between the sample metadata and the Factor values.
# Inspection of factor values.
# Inspection of the feature weights.
# Gene set enrichment analysis on the mRNA weights.

# 5.1Association analysis -------------------------------------------------
# Let’s test the association between MOFA Factors and Gender, survival outcome (dead vs alive) and age:
pdf("out/image/04_MOFA2_Sample_adj_top_featureScale_FactorCovariate.pdf",width = 4,height = 4)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("gender","age","pmi","ms","pathological_stage"),
                                  plot="log_pval"
)
dev.off()

correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("gender","age","pmi","ms","pathological_stage"),
                                  plot="log_pval", 
                                  return_data = T)
p_value <- 10^(-4.32628)
# pull the factor values per sample
df_factors <- get_factors(MOFAobject,
                          factors = "all",
                          as.data.frame = TRUE
)

# add the info from the LUT
left_join(df_factors,LUT_sample,by = "sample") %>% 
  ggplot(aes(y=value,x=pmi)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5) +
  facet_wrap(~factor) +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("out/image/04_MOFA2_Sample_adj_top_featureScale_FactorCovariate_scatter_pmi.pdf",width = 8,height = 6)

# add the info from the LUT
left_join(df_factors,LUT_sample,by = "sample") %>% 
  mutate(pathological_stage = factor(pathological_stage,levels = c("Control WM","Periplaque","CI","CA","core"))) %>% 
  ggplot(aes(y=value,x=pathological_stage)) + 
  # geom_smooth(method = "lm") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),alpha=0.5) +
  facet_wrap(~factor) +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("out/image/04_MOFA2_Sample_adj_top_featureScale_FactorCovariate_boxplot_PathologicalStage.pdf",width = 8,height = 6)

# Most Factors don’t have a clear association with any of the covariates. Only Factor 11 has a small association with survival outcome. We will go back to associations with clinical information at the end of the vignette.

# 5.2Plot factor values ---------------------------------------------------
# How do we interpret the factor values?
# Each factor captures a different source of variability in the data. Mathematically, each Factor is defined by a linear combination of the input features. As the data is centered prior to running MOFA, each Factor ordinates cells along a one-dimensional axis that is centered at zero. Samples with different signs manifest opposite phenotypes along the inferred axis of variation, with higher absolute value indicating a stronger effect. Note that the interpretation of MOFA factors is analogous to the interpretation of the principal components in PCA.
plot_factor(MOFAobject, 
            factors = c(1), 
            color_by = "Factor1"
)
ggsave("out/image/05_MOFA2_Sample_adj_top_featureScale_Factor1.pdf",width = 3,height = 3)

plot_factor(MOFAobject, 
            factors = c(1), 
            color_by = "pathological_stage"
)
ggsave("out/image/05_MOFA2_Sample_adj_top_featureScale_Factor1_factor.pdf",width = 3,height = 3)

# plot the top 7 factors
plot_factor(MOFAobject, 
            factors = c(1:4)
)

# 5.3Plot feature weights -------------------------------------------------
# How do we interpret the weights?
# The weights provide a score for each feature on each factor. Features with no association with the corresponding factor are expected to have values close to zero, whereas features with strong association with the factor are expected to have large absolute values. The sign of the weights indicates the direction of the effect: a positive weights indicates that the feature has higher levels in the cells with positive factor values, and vice-versa.

# 5.3.1Plot feature weights for somatic mutations -------------------------
# By looking at the variance explained plot, we saw that Factor 1 captures variation in all data modalities. Out of all omics, the somatic mutation data is a good place to start, as somatic mutations are very sparse, easy to interpret and any change in the DNA is likely to have downstream consequences to all other molecular layers. Let’s plot the weights:
plot_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
ggsave("out/image/06_MOFA2_Sample_adj_top_featureScale_Factor1_mRNA.pdf",width = 6,height = 3)

plot_weights(MOFAobject,
             view = "metabolic",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
ggsave("out/image/06_MOFA2_Sample_adj_top_featureScale_Factor1_metabolites.pdf",width = 6,height = 3)

# Notice that most features lie at zero, indicating that most features have no association with Factor 1. There is however one gene that clearly stands out: IGHV (immunoglobulin heavy chain variable region). This is the main clinical marker for CLL.

# An alternative visualistion to the full distribution of weights is to do a line plot that displays only the top features with the corresponding weight sign on the right:
plot_top_weights(MOFAobject,
                 view = "mRNA",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave("out/image/07_MOFA2_Sample_adj_top_featureScale_Factor1_mRNA.pdf",width = 4,height = 3)

plot_top_weights(MOFAobject,
                 view = "metabolic",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave("out/image/07_MOFA2_Sample_adj_top_featureScale_Factor1_metabolic.pdf",width = 6,height = 3)

# check for specific features
set.seed(1)
plot_factor(MOFAobject, 
            factors = 1, 
            color_by 
            = "SLC6A16",
            add_violin = TRUE,
            dodge = TRUE
)

set.seed(1)
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "ceramide_d18_1_14_0_d16_1_16_0",
            add_violin = TRUE,
            dodge = TRUE
)


# We can also plot Factor values coloured by other covariates, for example Gender. As shown above, this variable has no association with Factor 1:
plot_factor(MOFAobject,
            factors = 1,
            color_by = "ms",
            dodge = TRUE,
            add_violin = TRUE
)

plot_factor(MOFAobject,
            factors = 1,
            color_by = "pathological_stage",
            dodge = TRUE,
            add_violin = TRUE
) + scale_fill_viridis_d(name = "Pathological Stage",
                         option = "F",
                         end = 0.9, 
                         direction = -1) +
        xlab("Stage")



# 5.3.2Plot gene weights for mRNA expression ------------------------------
# From the variance explained plot we know that Factor 1 drives variation across all data modalities. Let’s visualise the mRNA expression changes that are associated with Factor 1:
plot_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10
)

# 5.3.3Plot molecular signatures in the input data ------------------------
# In this case we have a large amount of genes that have large positive and negative weights. Genes with large positive values will be more expressed in the samples with IGHV mutation, whereas genes with large negative values will be more expressed in the samples without the IGHV mutation. Let’s verify this. The function plot_data_scatter generates a scatterplot of Factor 1 values (x-axis) versus expression values (y-axis) for the top 4 genes with largest positive weight. Samples are coloured by IGHV status:
plot_data_scatter(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 9,
                  sign = "positive",
                  color_by = "ceramide_d18_1_14_0_d16_1_16_0"
) + 
  labs(y="RNA expression")
# if I want to explore the coexpression or correlation of individula factors I can extract the following
# extract the weights per factor
weights <- get_weights(MOFAobject, 
                       views = "all", 
                       factors = "all", 
                       as.data.frame = TRUE 
)

top_mrna <- weights %>%
  filter(factor=="Factor1") %>% 
  arrange(desc(abs(value))) %>% 
  filter(view == "mRNA") %>% 
  dplyr::slice(1:5) %>% 
  pull(feature) %>% 
  as.character()

top_metabolite <- weights %>%
  filter(factor=="Factor1") %>% 
  arrange(desc(abs(value))) %>% 
  filter(view == "metabolic") %>% 
  dplyr::slice(1:5) %>% 
  pull(feature) %>% 
  as.character()

# extrac the expresison values
# pick the top weight of factor 1
left_join(MOFAobject@data$mRNA %>% 
            data.frame() %>% 
            rownames_to_column("gene") %>% 
            filter(gene%in%top_mrna) %>% 
            pivot_longer(names_to = "sample",values_to = "exp",-gene),
          MOFAobject@data$metabolic %>% 
            data.frame() %>% 
            rownames_to_column("metabolite") %>% 
            filter(metabolite%in%top_metabolite) %>% 
            pivot_longer(names_to = "sample",values_to = "exp",-metabolite),
          by = "sample",suffix = c(".rna",".metabolite")
) %>% 
  ggplot(aes(x=exp.rna,y=exp.metabolite)) +
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5) +
  theme_bw() +
  theme(strip.background = element_blank())+
  facet_grid(gene~metabolite,scales = "free")
ggsave("out/image/07_MOFA2_Sample_adj_top_featureScale_Factor1_corr_top_mrna_vs_metabol.pdf",width = 15,height = 15)

# plot_data_heatmap has an interesting argument to “beautify” the heatmap: denoise = TRUE. Instead of plotting the (noisy) input data, we can plot the data reconstructed by the model, where noise has been removed:
pdf("out/image/07_MOFA2_Sample_adj_top_featureScale_Factor1_heatmap_top_mrna.pdf",width = 8,height = 6)
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 20,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = T,
                  scale = "row",
                  annotation_samples = "pathological_stage"
)
dev.off()

# plot_data_heatmap(MOFAobject,
#                   view = "mRNA",
#                   factor = 1,
#                   features = as.character(df_gene_fix$feature)[1:5],
#                   denoise = TRUE,
#                   cluster_rows = FALSE, cluster_cols = FALSE,
#                   show_rownames = TRUE, show_colnames = T,
#                   scale = "row",
#                   annotation_samples = "pathological_stage"
# )

pdf("out/image/07_MOFA2_Sample_adj_top_featureScale_Factor1_heatmap_top_metabol.pdf",width = 10,height = 6)
plot_data_heatmap(MOFAobject,
                  view = 2,
                  factor = 1,
                  features = 20,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = T,
                  scale = "row",
                  annotation_samples = "pathological_stage",
                  legend = F
) + theme(legend.position= "right",
        legend.direction = "horizontal")
dev.off()

# 9Gene set enrichment analysis (GSEA) ------------------------------------
# In addition to exploring the individual weights for each factor, we can use enrichment analysis to look for significant associations of factors to genesets. Here, we use the Reactome genesets for illustrations, which is contained in the MOFAdata package. For more details on how the GSEA works we encourage the users to read the GSEA vignette

# 9.1Load Reactome gene set annotations. ----------------------------------
# read the gene x term matrix for the enrichment analysis
# reactome_t2g <- readRDS("out/object/reactome_t2g.rds")
# annotation <- readRDS("out/object/reactome_t2g.rds")
# # kegg_t2g <- readRDS("out/object/kegg_t2g.rds")
# annotation <- readRDS("out/object/kegg_t2g.rds")
# # bp_t2g <- readRDS("out/object/bp_t2g.rds")
# annotation <- readRDS("out/object/bp_t2g.rds")

list_annotation <- list(reactome = readRDS("out/object/reactome_t2g.rds"),
                        kegg = readRDS("out/object/kegg_t2g.rds"),
                        bp = readRDS("out/object/bp_t2g.rds")
                        )

# 9.2Run enrichment analysis ----------------------------------------------
# These are the steps for doing Gene Set Enrichment Analysis (GSEA) with MOFA:
# (1) Define your gene set matrix: this can be specified as a binary matrix where rows are gene sets and columns are genes. A value of 1 indicates that gene j belongs to pathway i. A value of 0 indicates elsewise.
# (2) Select a gene set statistic: the statistic used to quantify the scores at the pathway level. Must be one of the following: mean.diff (difference in the average weight between foreground and background genes) or rank.sum (difference in the sum of ranks between foreground and background genes).
# (3) Select a statistical test: the statistical test used to compute the significance of the gene set statistics under a competitive null hypothesis. Must be one of the following: parametric (a simple and very liberal parametric t-test), cor.adj.parametric (parametric t-test adjusted by the correlation between features), permutation (unparametric, the null distribution is created by permuting the weights. This option is computationally expensive, but it preserves the correlation structure between features in the data.).

# An important consideration when running GSEA is that MOFA contains positive and negative weights. There will be cases where the genes with negative weights all belong to a specific pathway but genes with positive weights belong to other pathways. If this is true, doing GSEA with all of them together could dilute the signal. Hence, we recommend the user to do GSEA separately for (+) and (-) weights, and possibly also jointly with all weights.

# annotation <- list_annotation$reactome

# loop with all the annotation
list_GSEA <- pmap(list(list_annotation,names(list_annotation)),function(annotation,name){
  print(name)
  annotation[1:10,1:10]
  
  # GSEA on positive weights, with default options
  res.positive <- run_enrichment(MOFAobject, 
                                 feature.sets = annotation, 
                                 view = "mRNA",
                                 sign = "positive"
  )
  
  # GSEA on negative weights, with default options
  res.negative <- run_enrichment(MOFAobject, 
                                 feature.sets = annotation, 
                                 view = "mRNA",
                                 sign = "negative"
  )
  
  # no difference
  res.all <- run_enrichment(MOFAobject,
                            feature.sets = annotation,
                            view = "mRNA",
                            sign = "all"
  )
  
  # save the GSEA analysis
  saveRDS(object = list(res.positive=res.positive,
                        res.negative=res.negative,
                        res.all=res.all),
          file = paste0("out/object/MOFA2_Sample_adj_top_featureScale_GSEA_RNA_",name,".rds"))
  
  # The enrichment analysis returns a list of 5 elements:
  # feature.sets: the feature set matrix filtered by the genes that overlap with the MOFA model.
  # pval: the nominal p-values.
  # pval.adj: the FDR-adjusted p-values.
  # feature.statistics: the feature statistics (i.e. the weights).
  # set.statistics: matrices with the gene set statistics.
  # sigPathways: list with significant pathways per factor at a specified FDR threshold
  names(res.positive)
  names(res.negative)
  names(res.all)
  
  return(list(res.positive = res.positive,
              res.negative = res.negative,
              res.all = res.all))
})

# 9.2.1Plot enrichment analysis results -----------------------------------
# Plot an overview of the number of significant pathways per factor.
# It seems that most of the Factors do not have clear gene set signatures. A clear exception is Factor 5, which has a very strong enrichment for genes with positive weights.

pmap(list(list_GSEA,names(list_GSEA)),function(obj,name){
  # heatamp from the NES for each factor 
  plot_enrichment_heatmap(obj$res.positive)
  plot_enrichment_heatmap(obj$res.negative)
  plot_enrichment_heatmap(obj$res.all)
  dev.off()
  # save the plots for the top signatures
  # Let’s plot the GSEA results for Factor 5. It seems that this Factor is capturing differences in the stress response of the blood cells.
  pdf(paste0("out/image/08_MOFA2_Sample_adj_top_featureScale_Factor1_GSEA_",name,"_mRNA.pdf"),width = 15,height = 8)
  print(plot_enrichment(obj$res.positive, factor = 1, max.pathways = 15) + ggtitle(paste0("res.positive mRNA ",name)))
  print(plot_enrichment(obj$res.negative, factor = 1, max.pathways = 15) + ggtitle(paste0("res.negative mRNA ",name)))
  print(plot_enrichment(obj$res.all, factor = 1, max.pathways = 15) + ggtitle(paste0("res.all mRNA ",name)))
  dev.off()
  
  # It is always advised to not rely only on the p-values and to visualise which genes are driving the enrichment within each pathways. There are problematic cases where a single gene is driving the enrichment results in very small pathways.
  pdf(paste0("out/image/08_MOFA2_Sample_adj_top_featureScale_Factor1_GSEA_",name,"_mRNA_gene.pdf"),width = 15,height = 8)
  
  print(plot_enrichment_detailed(
    enrichment.results = obj$res.positive,
    factor = 1, 
    max.pathways = 10
  ) + ggtitle(paste0("res.positive mRNA ",name)))
  
  print(plot_enrichment_detailed(
    enrichment.results = obj$res.negative,
    factor = 1, 
    max.pathways = 10
  ) + ggtitle(paste0("res.negative mRNA ",name)))
  
  print(plot_enrichment_detailed(
    enrichment.results = obj$res.all,
    factor = 1, 
    max.pathways = 10
  ) + ggtitle(paste0("res.all mRNA ",name)))
  dev.off()
  
})

# 9.1Load Reactome gene set annotations. ----------------------------------
# read the gene x term matrix for the enrichment analysis
metabolites_t2g <- readRDS("out/object/metabolites_t2g.rds")

# sample GSEA on metabolites ----------------------------------------------
# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = metabolites_t2g, 
                               view = "metabolic",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = metabolites_t2g, 
                               view = "metabolic",
                               sign = "negative"
)

# all
res.all <- run_enrichment(MOFAobject, 
                          feature.sets = metabolites_t2g, 
                          view = "metabolic",
                          sign = "all"
)

saveRDS(object = list(res.positive=res.positive,
                      res.negative=res.negative,
                      res.all=res.all),
        file = "out/object/MOFA2_Sample_adj_top_featureScale_GSEA_METABOL_costume.rds")

# The enrichment analysis returns a list of 5 elements:
# feature.sets: the feature set matrix filtered by the genes that overlap with the MOFA model.
# pval: the nominal p-values.
# pval.adj: the FDR-adjusted p-values.
# feature.statistics: the feature statistics (i.e. the weights).
# set.statistics: matrices with the gene set statistics.
# sigPathways: list with significant pathways per factor at a specified FDR threshold
names(res.positive)
names(res.negative)
names(res.all)

# 9.2.1Plot enrichment analysis results -----------------------------------
# Plot an overview of the number of significant pathways per factor.
# It seems that most of the Factors do not have clear gene set signatures. A clear exception is Factor 5, which has a very strong enrichment for genes with positive weights.
plot_enrichment_heatmap(res.positive)
plot_enrichment_heatmap(res.negative)
plot_enrichment_heatmap(res.all)

# Let’s plot the GSEA results for Factor 5. It seems that this Factor is capturing differences in the stress response of the blood cells.
pdf("out/image/09_MOFA2_Sample_adj_top_featureScale_Factor1_GSEA_metabol_costume.pdf",width = 15,height = 8)
print(plot_enrichment(res.positive, factor = 1, max.pathways = 15) +
        ggtitle("res.positive matabol costume"))
print(plot_enrichment(res.negative, factor = 1, max.pathways = 15) +
        ggtitle("res.negative matabol costume"))
print(plot_enrichment(res.all, factor = 1, max.pathways = 15) +
        ggtitle("res.all matabol costume"))
dev.off()

# It is always advised to not rely only on the p-values and to visualise which genes are driving the enrichment within each pathways. There are problematic cases where a single gene is driving the enrichment results in very small pathways.
pdf("out/image/09_MOFA2_Sample_adj_top_featureScale_Factor1_GSEA_metabol_costume_feature.pdf",width = 15,height = 8)
print(plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 1, 
  max.pathways = 10) +
    ggtitle("res.positive matabol costume"))
print(plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 1, 
  max.pathways = 10) +
    ggtitle("res.negative matabol costume"))
print(plot_enrichment_detailed(
  enrichment.results = res.all,
  factor = 1, 
  max.pathways = 10) +
    ggtitle("res.all matabol costume"))
dev.off()

# save the objetc ---------------------------------------------------------
saveRDS(MOFAobject,file = "out/object/MOFA2_Sample_adj_top_featureScale_final.rds")


#Figure 3 ----
#Violin plot with MOFA factor 1----
mofa_factors <- data.frame(get_factors(MOFAobject)$group1) %>% rownames_to_column(var = "sample")
mofa_metadata <- data.frame(MOFAobject@samples_metadata)
mofa_factors <- mofa_factors %>% left_join(mofa_metadata)
cor.test(mofa_factors$Factor1, as.numeric(mofa_factors$pathological_stage))
#anova
anova_factor_1 <- aov(Factor1 ~ pathological_stage, mofa_factors)
summary(anova_factor_1)
TukeyHSD(anova_factor_1)
DunnettTest(mofa_factors$Factor1, mofa_factors$pathological_stage)

factor_1_violin <- ggplot(mofa_factors, aes(x = pathological_stage, y = Factor1)) +
        geom_violin(aes(fill = pathological_stage), alpha = 0.4, scale = "width") +
        geom_beeswarm(aes(fill = pathological_stage, shape = pathological_stage), size = 3) +
        scale_x_discrete(limits = c("Control WM", "Periplaque", "CI", "CA", "core"),
                         labels = c("Control\nWM", "Periplaque\nWM", "CI edge", "CA edge", "Core")) + 
        scale_color_viridis_d(name = NULL,
                              breaks = c("Control WM", "Periplaque", "CI", "CA", "core"),
                              labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                              option = "F",
                              end = 0.9, 
                              direction = -1) +
        scale_fill_viridis_d(name = NULL,
                             breaks = c("Control WM", "Periplaque", "CI", "CA", "core"),
                             labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                             option = "F",
                             end = 0.9, 
                             direction = -1) +
        scale_shape_manual(name = NULL,
                           breaks = c("Control WM", "Periplaque", "CI", "CA", "core"),
                           labels = c("Control WM", "Periplaque WM", "Chronic Inactive Edge", "Chronic Active Edge", "Core"),
                           values = c(21, 24, 23, 22, 25)) +
        xlab(element_blank()) + 
        ylab("Factor 1") +
        theme_bw() + 
        theme(panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20, face = "bold"),
              axis.title.y = element_text(size = 16, face = "bold"),
              axis.text = element_text(size = 15),
              axis.text.x = element_text(face = "bold"),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size = 14),
              legend.position = "none") 

# pdf("factor_1_violin.pdf", height = 7, width = 8)
# factor_1_violin
# dev.off()



#Pathway analysis plots form MOFA factor 1 ----
bp_pos <- data.frame(list_GSEA$bp$res.positive$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "BP", direction = "+")
bp_neg <- data.frame(list_GSEA$bp$res.negative$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "BP", direction = "-")
bp <- bind_rows(bp_pos, bp_neg) %>% rownames_to_column(var = "pathway") %>% arrange(desc(direction), desc(Factor1))
bp <- bp %>% mutate(pathway = str_replace(pathway, "GOBP_", "")) %>%
        mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
        mutate(pathway = str_to_title(pathway))
bp$pathway <- factor(bp$pathway, levels = bp$pathway)

kegg_pos <- data.frame(list_GSEA$kegg$res.positive$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "KEGG", direction = "+")
kegg_neg <- data.frame(list_GSEA$kegg$res.negative$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "KEGG", direction = "-")
kegg <- bind_rows(kegg_pos, kegg_neg) %>% rownames_to_column(var = "pathway")

reactome_pos <- data.frame(list_GSEA$reactome$res.positive$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "Reactome", direction = "+")
reactome_neg <- data.frame(list_GSEA$reactome$res.negative$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "Reactome", direction = "-")
reactome <- bind_rows(reactome_pos, reactome_neg) %>% rownames_to_column(var = "pathway")


bp_plot <- ggplot(bp, aes(x = -log10(Factor1), y = pathway, group = pathway, color = direction, fill = direction)) +
        geom_bar(stat="identity",width=0.15) +
        geom_point(size = 5) +
        scale_color_manual(name = NULL,
                           breaks = c("+", "-"),
                           labels = c("upregulated pathways", "downregulated pathways"),
                           values = c("#D6604D", "#2166AC")) +
        scale_fill_manual(name = NULL,
                          breaks = c("+", "-"),
                          labels = c("upregulated pathways", "downregulated pathways"),
                          values = c("#D6604D", "#2166AC")) +
        geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.9, color = "gray30") +
        xlab(expression(paste("\u2013", log[10],"(q-value)"))) +
        ylab(element_blank()) +
        theme_classic() +
        theme(panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20, face = "bold"),
              axis.title.x = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 13),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size = 14),
              legend.position = "none")  

# pdf("bp_pathways.pdf", height = 8, width = 10)
# bp_plot
# dev.off()


metabolic_pos <- data.frame(res.positive$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "metabolic", direction = "+")
metabolic_neg <- data.frame(res.negative$pval.adj) %>% filter(Factor1 < 0.05) %>% slice_min(Factor1, n=10) %>% mutate(database = "metabolic", direction = "-")
metabolic <- bind_rows(metabolic_pos, metabolic_neg) %>% rownames_to_column(var = "pathway") %>% arrange(desc(direction), desc(Factor1))

metabolic$pathway <- factor(metabolic$pathway, levels = metabolic$pathway)


metabolic_plot <- ggplot(metabolic, aes(x = -log10(Factor1), y = pathway, group = pathway, color = direction, fill = direction)) +
        geom_bar(stat="identity",width=0.08) +
        geom_point(size = 5) +
        scale_color_manual(breaks = c("+", "-"),
                           values = c("#D6604D", "#2166AC")) +
        scale_fill_manual(breaks = c("+", "-"),
                          values = c("#D6604D", "#2166AC")) +
        geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.9, color = "gray30") +
        xlab(expression(paste("\u2013", log[10],"(q-value)"))) +
        ylab(element_blank()) +
        theme_classic() +
        theme(panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title = element_text(hjust = 0.5, margin=margin(0,0,30,0), size = 20, face = "bold"),
              axis.title.x = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 13),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size = 14),
              legend.position = "none") 
bp_plot + metabolic_plot + plot_layout(ncol = 1)

# pdf("mofa_pathways.pdf", height = 11, width = 8)
# bp_plot / metabolic_plot
# dev.off()
