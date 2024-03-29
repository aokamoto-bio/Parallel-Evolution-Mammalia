## Script to generate manuscript supplemental figures
## Alexander Okamoto
## June 15, 2023

#load packages
library(data.table)
library(tidyverse)
library(nlme)
library(phytools)
library("Hmisc", lib.loc="~/R-4.2.1")
library(RERconverge) #for permulations
library(parallel) #for running code in parallel
library(ggplot2) #for plotting
library(ggpubr) #for final plot
library(ggpmisc) #for adding equation to plot
library(grid)
library(eulerr) #add a Venn diagram
library(rphylopic) #add animal sillouettes 
library(cowplot) #a better ggarrange

#create new directory for results 
try(system("mkdir Big_Loop_Output2"))

print("Number of cores")
detectCores()

#load trait data
traits_df <- read.table(file = "~/Polygenic_Inference/Big_Loop_Output2/traits_df.tsv", header = TRUE)

#load the SNP data
SNP_df_mammal_cons <- fread(file = "SNP_df_mammal_cons.tsv", header = TRUE,  showProgress= TRUE)
SNP_df_primates_cons <- fread(file = "SNP_df_primates_cons.tsv", header = TRUE,  showProgress= TRUE)
#load the monophytic SNP data
#this means SNPs with an allele fixed in one lineage and something else in all other species
SNP_df_mammal_monophyletic <- fread(file = "SNP_df_mammal_monophyletic.tsv", header = TRUE) %>% as.data.frame()
#load variant data, the filter is just to reduced the dataset size quickly
quality_variants_matched <- fread(file = "quality_variants_matched.tsv", header = TRUE,  showProgress= TRUE) %>% filter(minor_AF > 0.01)
SNP_df_primates_cons_merged <- merge(x = quality_variants_matched, y = SNP_df_primates_cons, by.x = 'Row.names', by.y = 'rsID')
SNP_df_mammal_cons_merged <- merge(x = quality_variants_matched, y = SNP_df_mammal_cons, by.x = 'Row.names', by.y = 'rsID')
#get primate unique sets
SNP_primates_unique <- SNP_df_primates_cons_merged[!SNP_df_primates_cons_merged$Row.names %in% SNP_df_mammal_cons_merged$Row.names, ]

#load the phenotype data
ZIMS_data_full <- fread(file = "ZIMS Medical Reference Data 2022 - Sheet1.tsv", header = TRUE,  showProgress= TRUE)

#load pantheria data
#this is used to provide body mass data
pantheria <-fread("Pantheria_data.csv", header=TRUE, showProgress = TRUE)
#make a couple edits so that R will handle the data properly
pantheria$species[1496] <- "Giraffa_tippelskirchi"
pantheria[pantheria == -999] <- NA

#load CC phenotype data
#take phenotype means combining the sexes for each strain
JAX_CC_Pheno_Data <- read.csv("JAX_CC_Pheno_Data.csv") %>% group_by(strain) %>% summarise_at(traits_df$trait_abrev, .funs = mean, na.rm = TRUE)
#modify CC names so phenotype and genotype match
JAX_CC_Pheno_Data$strain <- sub("/", ".", JAX_CC_Pheno_Data$strain)
JAX_CC_Pheno_Data$strain <- sub("J", "", JAX_CC_Pheno_Data$strain)

#tidy the data for easier subsetting later on
JAX_CC_Pheno_Data_tidy <- JAX_CC_Pheno_Data %>% pivot_longer(cols = traits_df$trait_abrev, 
                                                             names_to = "trait", 
                                                             values_to = "trait_mean",
                                                             values_drop_na = TRUE)

#load CC genotype data
CC_genotypes <- read.csv("CC_genotypes.csv")

#reduce SNP dataset to informative columns for this analysis
SNP_df_mus <- SNP_df_mammal_cons_merged %>% select(Row.names, variant, chr, pos, ref, alt)
#merge CC data with mammal SNPs
#filter out CC rows with indels or no information
SNP_df_mouse_merged_filtered <- merge(x = SNP_df_mus, y = CC_genotypes, by.x = 'Row.names', by.y = 'rsID') %>% filter(observed != 'A') %>% filter(observed != 'T') %>% filter(observed != 'G') %>% filter(observed != 'C') %>% filter(observed !="")

#load in mammals data
mammals_info <- fread(file = "241_mammals_info.csv", header = TRUE) %>% select(1:3)

#load info on the 241 mammals
species_list <- mammals_info  %>% select(Species) %>% filter(Species != 'Homo_sapiens') 
#make primate filtered subsets
primates <- mammals_info[,1:3]  %>% filter(Order == 'PRIMATES')  %>% filter(Species != 'Homo_sapiens')

#clean and expand pantheria data
pantheria_body_mass <- pantheria %>% filter(species %in% species_list$Species) %>% filter(`5-1_AdultBodyMass_g` != 'NA')
#add additional body size estimates
body_size_pheno_data <- data.frame(species = species_list$Species)
body_size_pheno_data <- merge(body_size_pheno_data, pantheria_body_mass[,c(5, 7)], all = TRUE)
colnames(body_size_pheno_data) <- c("species", "body_mass_g")

#supplement additional data
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Balaenoptera_bonaerensis')] <- 8350000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Bos_indicus')] <- 280000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Bos_mutus')] <- 546250
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Camelus_ferus')] <- 475000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Capra_aegagrus')] <- 53750
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Ceratotherium_simum_cottoni')] <- 1900000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Chlorocebus_sabaeus')] <- 3975
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Cricetulus_griseus')] <- 40
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Ellobius_lutescens')] <- 71
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Equus_przewalskii')] <- 250000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Eubalaena_japonica')] <- 57595750
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Eulemur_flavifrons')] <- 1850
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Fukomys_damarensis')] <- 153.25 
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Galeopterus_variegatus')] <- 1450
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Ictidomys_tridecemlineatus')] <- 190
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Murina_feae')] <- 4.8
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Mus_pahari')] <- 28
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Myotis_davidii')] <- 5.95
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Nannospalax_galili')] <- 162.5
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Neomonachus_schauinslandi')] <- 204000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Neophocaena_asiaeorientalis')] <- 56000
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Piliocolobus_tephrosceles')] <- 9700
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Rhinolophus_sinicus')] <- 9.9
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Spermophilus_dauricus')] <- 223.8 
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Spilogale_gracilis')] <- 492.5
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Tonatia_saurophila')] <- 28.56
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Tupaia_chinensis')] <- 125
body_size_pheno_data$body_mass_g[which(body_size_pheno_data$species=='Vicugna_pacos')] <- 44400
body_size_pheno_data <- body_size_pheno_data %>% filter(body_mass_g > 0)
body_size_pheno_data$rank <- rank(body_size_pheno_data$body_mass_g)

#load species tree
AllSpeciesTree <- read.tree("species_tree.nwk")
#modify tips to match Zoonomia data
AllSpeciesTree$tip.label[AllSpeciesTree$tip.label == "Plecturocebus_donacophilus"] <- "Callicebus_donacophilus"
AllSpeciesTree$tip.label[AllSpeciesTree$tip.label == "Monachus_schauinslandi"] <- "Neomonachus_schauinslandi"

#load function to analyze traits
analyze_trait <- function(trait_species, trait_values, trait_GWAS, dataset){
  
  #create a table to store output results
  #the rank is here is often redundant but useful after filtering species (e.g. reducing set to primates)
  results_table <- data.frame(species = trait_species, trait_rank = rank(trait_values), pos_ct = 0, neg_ct = 0, total = 0, pos_prop = 0)
  
  #process GWAS info
  #read in GWAS summary stats and filter by pval and MAF
  GWAS_sum_stats <- trait_GWAS %>% dplyr::filter(pval < 5e-8 & minor_AF >0.01)
  bps <- c("A", "C", "G", "T")
  #merge GWAS summary statistics with SNP data, remove indel matches, and create alt effect column
  full_data_set <- merge(x = GWAS_sum_stats, y = dataset, by = 'variant') %>% filter(ref %in% bps, alt %in% bps) %>% mutate(alt_effect = sign(beta)) %>% as.data.frame()
  
  results_table$n <- nrow(full_data_set)
  
  #this checks for alleles not matching a human chromosome 
  for (animal in trait_species){
    #calculate trait increasing proportion of alleles
    neg_count <- length(which(full_data_set[animal] == full_data_set$ref & full_data_set$alt_effect == 1)) + length(which(full_data_set[animal] == full_data_set$alt & full_data_set$alt_effect == -1))
    results_table[results_table$species == animal, 'neg_ct'] <- neg_count
    pos_count <- length(which(full_data_set[animal] == full_data_set$ref & full_data_set$alt_effect == -1)) + length(which(full_data_set[animal] == full_data_set$alt & full_data_set$alt_effect == 1))
    results_table[results_table$species == animal, 'pos_ct'] <- pos_count
    results_table[results_table$species == animal, 'total'] <- pos_count + neg_count
    results_table[results_table$species == animal, 'pos_prop'] <- pos_count/(pos_count + neg_count)
    #calculate score
    results_table$pos_prop_rank <- rank(results_table$pos_prop)
  }
  return(results_table)
}

#get causative SNPs for a given dataset using wilcoxon rank sum test
causative_SNP_test <- function(s, pheno_ranks, species_list, SNP_set, test_dir){
  ref_ranks <- pheno_ranks[as.character(SNP_set[s, 'ref']) == SNP_set[s, species_list]]
  alt_ranks <- pheno_ranks[as.character(SNP_set[s,'alt']) == SNP_set[s, species_list]]
  utest <- wilcox.test(ref_ranks, alt_ranks, exact = FALSE, alternative = test_dir)
  
  return(list(utest$p.value, mean(ref_ranks), mean(alt_ranks), as.numeric(utest$statistic)))
}
#_________________________________


#fix value for RBC figure
for (i in c(1:7, 10:11)){
  #for (i in c(1:2)){
  
  #make sure file exists:
  if(file.exists(traits_df$GWAS_file[i]) == F){
    break
  }
  print(paste("Running ", traits_df$trait[i]))
  #load phenotypic data
  #check if phenotype data comes from ZIMS:
  if(traits_df$pheno_source[i] == "ZIMS"){
    #check if data needs serum filter:
    if(traits_df$Serum[i]){
      ZIMS_temp <- ZIMS_data_full %>% filter(Test == traits_df$ZIMS_ID[i]) %>% filter(`Sample Type` == "Serum")
      ZIMS_temp$rank <- rank(ZIMS_temp$Mean)
      pheno_data_temp <- data.frame(species =  ZIMS_temp$Species, pheno_rank = ZIMS_temp$rank)
      
    }else{
      ZIMS_temp <- ZIMS_data_full %>% filter(Test == traits_df$ZIMS_ID[i])
      ZIMS_temp$rank <- rank(ZIMS_temp$Mean)
      pheno_data_temp <- data.frame(species =  ZIMS_temp$Species, pheno_rank = ZIMS_temp$rank)
    }
  }
  if(traits_df$pheno_source[i] == "PANTHERIA_PLUS"){
    pheno_data_temp <- data.frame(species =  body_size_pheno_data$species, pheno_rank = body_size_pheno_data$rank)
  }
  #keep track of number of mammals with phenotypic data for the trait
  traits_df$phenodatspecn[i] <- length(unique(pheno_data_temp$species))
  pheno_data_temp_primates <- pheno_data_temp[which(pheno_data_temp$species %in% primates$Species),]
  #keep track of number of primates with phenotypic data for the trait
  traits_df$phenodatspecn_prim[i] <- length(unique(pheno_data_temp_primates$species))
  
  #read in GWAS summary statistics 
  GWAS_sum_stats <- fread(file = traits_df$GWAS_file[i], header = TRUE, showProgress= F) %>% filter(minor_AF >0.01)
  GWAS_sum_stats_sig <- GWAS_sum_stats %>% filter(pval < 5e-8)
  #get results for mammals, primates, and primate unique
  trait_results_mammals <- analyze_trait(trait_species = pheno_data_temp$species, trait_values = pheno_data_temp$pheno_rank, GWAS_sum_stats, SNP_df_mammal_cons_merged)
  trait_results_primates <- analyze_trait(trait_species = pheno_data_temp_primates$species, trait_values = pheno_data_temp_primates$pheno_rank, GWAS_sum_stats, SNP_df_primates_cons_merged)
  trait_results_primates_unique <- analyze_trait(trait_species = pheno_data_temp_primates$species, trait_values = pheno_data_temp_primates$pheno_rank, GWAS_sum_stats, SNP_primates_unique)
  
  #run PGLS
  
  #clean tree if needed
  species_to_drop <- setdiff(mammals_info$Species, trait_results_mammals$species)
  #remove species from tree with no phenotypic data available
  FilteredTree <- drop.tip(AllSpeciesTree, species_to_drop)
  #extract the species in the appropriate order 
  PGLS_FilteredTree <- FilteredTree$tip.label
  
  #make brownian evolution model 
  PGLS_model <- corBrownian(1, phy = FilteredTree, form = ~PGLS_FilteredTree)
  trait_results_mammals <- merge(x= trait_results_mammals, y = mammals_info, by.x = "species", by.y = "Species") %>% as.data.frame()
  trait_results_mammals$Group <- "Other"
  trait_results_mammals$Group[which(trait_results_mammals$Order=="PRIMATES")] <- "Primates"
  
  #generate model
  PGLS_posprop_model <- gls(pos_prop_rank ~ trait_rank, trait_results_mammals, correlation = PGLS_model, na.action = na.omit)
  PGLS_posprop_model_sum <- summary(PGLS_posprop_model)
  #summary(PGLS_posprop_model) #to visualize manually
  
  #generate plot
  PGLS_plot_mam <- ggplot(trait_results_mammals, aes(x=trait_rank, y=pos_prop_rank, shape = Group)) + geom_point(color="blue", size = 0.5, stroke = 0.2)  + 
    labs(x = paste(traits_df$trait_abrev[i], " Rank", sep = ""), y = "Score Rank \n (Mammals)") + 
    scale_shape_manual(values = c("Primates"= 15, "Other"= 1)) + 
    theme(legend.position = "none", panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.title=element_text(size = 8), 
          axis.text=element_text(size = 8),
          text = element_text(family = "Arial"), 
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm")
    ) 
  
  if(PGLS_posprop_model_sum$tTable[2,4] < 0.05){
    PGLS_plot_mam <- PGLS_plot_mam + geom_abline(slope = coef(PGLS_posprop_model)[[2]], intercept = coef(PGLS_posprop_model)[[1]], color ='black')
    
  }
  
  #PGLS for primates
  species_to_drop_primates <- setdiff(mammals_info$Species, trait_results_primates$species)
  #remove species from tree with no phenotypic data available
  FilteredTree_primates <- drop.tip(AllSpeciesTree, species_to_drop_primates)
  #extract the species in the appropriate order 
  PGLS_FilteredTree_primates <- FilteredTree_primates$tip.label
  
  #make brownian evolution model 
  PGLS_model_primates <- corBrownian(1, phy = FilteredTree_primates, form = ~PGLS_FilteredTree_primates)
  
  #add family and order data 
  trait_results_primates <- merge(x= trait_results_primates, y = mammals_info, by.x = "species", by.y = "Species") %>% as.data.frame()
  #trait_results_primates_unique <- merge(x= trait_results_primates_unique, y = mammals_info, by.x = "species", by.y = "Species") %>% as.data.frame()
  
  #generate model
  PGLS_posprop_model_primates <- gls(pos_prop_rank ~ trait_rank, trait_results_primates, correlation = PGLS_model_primates, na.action = na.omit)
  PGLS_posprop_model_sum_primates <- summary(PGLS_posprop_model_primates)
  #summary(PGLS_posprop_model) #to visualize manually
  
  #generate plot
  PGLS_plot_primates <- ggplot(trait_results_primates, aes(x=trait_rank, y=pos_prop_rank)) + geom_point(color="red", size = 0.5, stroke = 0.2)  + 
    labs(x = paste(traits_df$trait_abrev[i], " Rank", sep = ""), y = "Score Rank \n (Primate)") + 
    theme(legend.position = "none", panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.title=element_text(size = 8), 
          axis.text=element_text(size = 8),
          text = element_text(family = "Arial"), 
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm")
    ) 
  
  if(PGLS_posprop_model_sum_primates$tTable[2,4] < 0.05){
    PGLS_plot_primates <- PGLS_plot_primates + geom_abline(slope = coef(PGLS_posprop_model_primates)[[2]], intercept = coef(PGLS_posprop_model_primates)[[1]], color ='black')
    
  }
  
  # #primates unique
  # 
  # #generate model
  # PGLS_posprop_model_primates_unique <- gls(pos_prop_rank ~ trait_rank, trait_results_primates_unique, correlation = PGLS_model_primates, na.action = na.omit)
  # PGLS_posprop_model_sum_primates_unique <- summary(PGLS_posprop_model_primates_unique)
  # #summary(PGLS_posprop_model) #to visualize manually
  # 
  # #generate plot
  # PGLS_plot_primates_unique <- ggplot(trait_results_primates_unique, aes(x=trait_rank, y=pos_prop_rank)) + geom_point(color = "pink") + 
  #   labs(x = "RBC Rank", y = "Score Rank (Primate Unique)") + 
  #   geom_abline(slope = coef(PGLS_posprop_model_primates_unique)[[2]], intercept = coef(PGLS_posprop_model_primates_unique)[[1]], color ='black') + 
  #   theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=6), axis.text=element_text(size=6))
  # 
  # print("done PGLS")
  
  
  #COLLABORATIVE CROSS MOUSE ANALYSIS
  
  #create dataframe for specific phenotype
  JAX_CC_Pheno_Data_tidy_reduced <- JAX_CC_Pheno_Data_tidy %>% dplyr::filter(trait == traits_df$trait_abrev[i])
  JAX_CC_Pheno_Data_tidy_reduced$pheno_rank <- rank(JAX_CC_Pheno_Data_tidy_reduced$trait_mean)
  
  #analyze trait in CC mice strains
  CC_results <- analyze_trait(trait_species = JAX_CC_Pheno_Data_tidy_reduced$strain, trait_values = JAX_CC_Pheno_Data_tidy_reduced$pheno_rank, trait_GWAS = GWAS_sum_stats, dataset = SNP_df_mouse_merged_filtered)
  
  #test the correlation
  CC_corr_test <- cor.test(CC_results$trait_rank, CC_results$pos_prop_rank, method = "spearman", alternative = "greater")
  
  #plot results 
  #plot results 
  CC_plot <- ggplot(data = CC_results, aes(x= trait_rank, y = pos_prop_rank)) + geom_point(size = 0.5) + 
    labs( x = paste(traits_df$trait_abrev[i], " Rank", sep = ""), y= "Score Rank \n (CC mice)") +
    theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), axis.title=element_text(size = 8), 
          axis.text=element_text(size = 8), 
          text = element_text(family = "Arial"), 
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm")
    ) 
  #+ add_phylopic(uuid="6b2b98f6-f879-445f-9ac2-2c2563157025", x= 15, y= 3, ysize=4, alpha = 1, color = "black")
  if(CC_corr_test$p.value < 0.05){
    CC_plot <- CC_plot + geom_smooth(method = "lm", se = F, color ='black', linewidth = 0.5) 
    
  }
  
  
  
  
  print("done CC analysis")
  
  
  
  #remove species from tree with no phenotypic data available
  MonoFilteredTree <- drop.tip(AllSpeciesTree, "Homo_sapiens")
  
  #match nodes from all species tree (minus human) and filtered free 
  nodes_df <- as.data.frame(matchNodes(MonoFilteredTree, FilteredTree))
  colnames(nodes_df) <- c("all_tree_node", "reduced_tree_node")
  
  #create a dataframe of results for each node 
  # count total monophyletic SNPs, number overlapping GWAS significant SNPs, and the ratio
  phylo_results_df <- nodes_df %>% dplyr::filter(!is.na(reduced_tree_node)) %>% add_column(monoSNP_n = NA, GWAS_overlap = NA, trait_change = NA, causative_n = NA, branch_length = NA)
  if(phylo_results_df$reduced_tree_node[1] == length(FilteredTree$tip.label) + 1){
    phylo_results_df <- phylo_results_df[2:nrow(phylo_results_df),]  
  }    
  #BRANCH FIXATION TEST
  
  #get phenotype ranks organized
  pheno_ranks <- pheno_data_temp$pheno_rank
  names(pheno_ranks) <-  pheno_data_temp$species
  
  #add rsID data to GWAS summary stats
  GWAS_sum_stats_sig_merged <- merge(x = GWAS_sum_stats_sig, y = quality_variants_matched, by = "variant", all.y = F)
  GWAS_sum_stats_sig_merged <- GWAS_sum_stats_sig_merged[!duplicated(GWAS_sum_stats_sig_merged)]
  
  #make function to calculate the amount of phenotypic change at a node:
  get_edge_change2 <- function(tree, values, node){
    #so the structure of the tree edges data frame is the numeric ID of the clade outside first, then the ID of the second clade which is contained in the first
    
    #get the lists of species in the clade of interest
    clade1 <- tree$tip.label[getDescendants(tree, node = node)[(getDescendants(tree, node = node ) <= length(tree$tip.label))]]
    #clade_outer <- tree$tip.label[getDescendants(tree, node = edge[1])[(getDescendants(tree, node = edge[1] ) <=length(tree$tip.label))]]
    #clade2 <- setdiff(clade_outer, clade1)
    
    if(length(clade1) > 1){
      #now get the phenotypic data
      clade1_phenos <- values[clade1]
      fit<-fastAnc(tree, x = values, vars=F, CI=F)
      
      anc_val <- fit[as.character(getParent(tree = tree, node = node))]
      node_val <- fit[as.character(node)]
      #calculate change in trait value from LCA and decescendants 
      clade_diff <- node_val - anc_val
      
    }
    else{
      return(NA)
    }
    return(clade_diff)
  }
  
  #get tree branch lengths
  distances <- dist.nodes(MonoFilteredTree)
  #get important values for each node in the filtered tree
  
  for (x in 1:length(phylo_results_df$all_tree_node)){
    #for (x in 1:25){ #testing
    #get the relevant node
    node <- phylo_results_df$all_tree_node[x]
    #identify the SNPs monophyletic on this clade's branch
    clade_monoSNPs <- SNP_df_mammal_monophyletic[SNP_df_mammal_monophyletic$monophyletic_node==node,] %>% pull(rsID)
    #store the number of monophyletic SNPs in the results df
    phylo_results_df$monoSNP_n[x] <- length(clade_monoSNPs)
    #find how many of those monophyletic SNPs are GWAS significant for that trait
    clade_overlap <- intersect(GWAS_sum_stats_sig_merged$Row.names, clade_monoSNPs) 
    #store result
    phylo_results_df$GWAS_overlap[x] <- length(clade_overlap)
    #calculate change on basal branch for this clade 
    #phylo_results_df$trait_change[x] <- get_edge_change(tree = FilteredTree, values = pheno_ranks, node = phylo_results_df$reduced_tree_node[x])
    phylo_results_df$trait_change[x] <- get_edge_change2(tree = FilteredTree, values = pheno_ranks, node = phylo_results_df$reduced_tree_node[x])
    #now test number of SNPs which seem causative
    #example species
    clade_species <- MonoFilteredTree$tip.label[getDescendants(MonoFilteredTree, node = node)[(getDescendants(MonoFilteredTree, node = node) <= length(MonoFilteredTree$tip.label))]]
    #get SNPs for the focal node and merge with GWAS summary stats
    causative_monoSNPs <- SNP_df_mammal_monophyletic[SNP_df_mammal_monophyletic$monophyletic_node==node,] %>% 
      inner_join(y = GWAS_sum_stats_sig_merged, by = c('rsID' = "Row.names")) 
    #identify SNPs for which the ref or alt allele matches the focal clade
    ref_monoSNPs <- causative_monoSNPs[which(causative_monoSNPs[clade_species[1]] == causative_monoSNPs$ref),]
    alt_monoSNPs <- causative_monoSNPs[which(causative_monoSNPs[clade_species[1]] == causative_monoSNPs$alt),]
    
    #first for positive clade trait change
    if (phylo_results_df$trait_change[x] > 0){
      ref_matches <- ref_monoSNPs %>% dplyr::filter(beta < 0)
      alt_matches <- alt_monoSNPs %>% dplyr::filter(beta > 0)
      all_matches <- rbind(ref_matches, alt_matches)
      ref_fails <- ref_monoSNPs %>% dplyr::filter(beta > 0)
      alt_fails <- alt_monoSNPs %>% dplyr::filter(beta < 0)
      all_fails <- rbind(ref_fails, alt_fails)
      phylo_results_df$causative_n[x] <- nrow(all_matches)
      phylo_results_df$fails_n[x] <- nrow(all_fails)
    }
    #now for negative trait change
    if (phylo_results_df$trait_change[x] < 0){
      ref_matches <- ref_monoSNPs %>% dplyr::filter(beta > 0)
      alt_matches <- alt_monoSNPs %>% dplyr::filter(beta < 0)
      all_matches <- rbind(ref_matches, alt_matches)
      ref_fails <- ref_monoSNPs %>% dplyr::filter(beta < 0)
      alt_fails <- alt_monoSNPs %>% dplyr::filter(beta > 0)
      all_fails <- rbind(ref_fails, alt_fails)
      phylo_results_df$causative_n[x] <- nrow(all_matches)
      phylo_results_df$fails_n[x] <- nrow(all_fails)
    }
    
    #add branch length data
    phylo_results_df$branch_length[x] <- distances[getParent(tree = MonoFilteredTree, node = node), node]
    
  }
  
  
  #sum causative and mismatch
  phylo_results_df$ref_alt_match <- phylo_results_df$causative_n + phylo_results_df$fails_n
  phylo_results_df$causative_ratio <- phylo_results_df$causative_n/phylo_results_df$ref_alt_match
  #filter out problematic rows
  phylo_results_df <-phylo_results_df %>% dplyr::filter(monoSNP_n > 0, !is.na(trait_change), trait_change != 0)
  #convert to ranks and calculate correlation
  phylo_results_df$abs_trait_change <- abs(phylo_results_df$trait_change)
  
  #add residuals of 
  fit <- lm(causative_n ~ fails_n, data=phylo_results_df) 
  phylo_results_df$residual <- fit$residuals
  
  #fit GWAS against 
  gwasfit <- lm(GWAS_overlap ~ branch_length, data=phylo_results_df) 
  phylo_results_df$GWASresidual <- gwasfit$residuals
  
  #it's important to know which lineages work and which don't 
  #so let's figure that out! 
  
  phylo_results_df$species <- NA
  phylo_results_df$Order <- "Deep"
  
  for (x in 1:length(phylo_results_df$all_tree_node)){
    #get the relevant node
    node <- phylo_results_df$all_tree_node[x]
    #get species
    clade_species <- MonoFilteredTree$tip.label[getDescendants(MonoFilteredTree, node = node)[(getDescendants(MonoFilteredTree, node = node) <= length(MonoFilteredTree$tip.label))]]
    phylo_results_df$species[x] <- paste(clade_species, collapse = ", ")
    Orders <- unique(mammals_info$Order[which(mammals_info$Species %in% clade_species)])
    if(length(Orders) == 1){
      phylo_results_df$Order[x] <- Orders
    }
  } #end add species loop
  
  
  phylo_results_df_reduced <- phylo_results_df  %>% dplyr::filter(abs_trait_change > quantile(phylo_results_df$abs_trait_change, probs = 0.90))
  
  
  phylo_results_df$Order <- factor(phylo_results_df$Order)
  
  
  clades_wilcox_corr <- wilcox.test(phylo_results_df$abs_trait_change[which(phylo_results_df$residual > 0)], phylo_results_df$abs_trait_change[which(phylo_results_df$residual <0)], alternative = "greater")
  
  #plot by residual
  phylo_results_df$residual_sign <- "Negative"
  phylo_results_df$residual_sign[which(phylo_results_df$residual > 0)] <- "Positive"
  phylo_results_df$residual_sign <- factor(phylo_results_df$residual_sign, levels = c("Positive", "Negative"))
  
  
  #plot and save results 
  #optimal X axis is probably the absolute value of the trait change (calculated using the FastAnc function for both focal node and parent instead of the mean of descendant species)
  #this should capture the extremity of the change in phenotype even if the actually magnitude is meaningless
  monoPlot <- ggplot(phylo_results_df, aes(y = abs_trait_change, x = residual_sign)) + geom_boxplot(fill = "gray", linewidth = 0.1, outlier.size = 0.0001) + 
    labs(y = "| Trait Change |", x = "Residual of sign matches \n to mismatches", title = NULL )+ 
    theme_bw() +
    theme(axis.line = element_line(colour = "black"), 
          axis.title=element_text(size = 8), 
          axis.text=element_text(size = 8), 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          text = element_text(family = "Arial"),
          axis.text.x = element_text(size = 6),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm")
    ) +
    annotate(geom="text", x=1.5, y=(max(phylo_results_df$abs_trait_change) + 3), label=paste("Wilcoxon-rank sum, p < ", round(clades_wilcox_corr$p.value, digits = 3), sep = ""),
             color="black", size = 6/.pt)
  # 
  monoPlot
  composite_figure <- plot_grid(PGLS_plot_mam, PGLS_plot_primates, CC_plot, monoPlot, nrow = 2, ncol =2, labels = "AUTO", label_size = 8, align = "v")
  
  
  #save plot
  ggsave(filename = paste("Figure_S", i, "_", traits_df$trait_abrev[i], ".png", sep = ""), device = "png", plot = composite_figure, path = "~/Polygenic_Inference/Big_Loop_Output2", width = 7, height = 7)
  ggsave(filename = paste("Figure_S", i, "_", traits_df$trait_abrev[i], ".tiff", sep = ""), device = "tiff", plot = composite_figure, path = "~/Polygenic_Inference/Big_Loop_Output2", width = 11, height = 11, units = "cm", dpi = 1200)
  
}
