
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")


source("/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/envs/r_env_jgteach/renv/activate.R")


.libPaths("/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/envs/r_env_jgteach/renv/library/linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu")
library(tidyverse)
library(hierfstat)
library(boot)
library(devtools)
library("JGTeach")
library(vcfR)
library(VGAM)
library(data.table)
library(gaston)

set.seed(123)

#### Get correct file from command line arguments
args <- commandArgs(trailingOnly = TRUE)
file_index <- as.numeric(args[1])

Simulation_directory <-"/scratch/jntoko/vcf_output/stepping/"
Neutral_file_name <- list.files(path = Simulation_directory, pattern = "*.vcf", full.names = TRUE)
file <- Neutral_file_name[file_index]

replicate_number <- file_index

### Declare and create output directories
rep_dir <- file.path("/scratch/jntoko/results_GW_final_stepping", paste0("Rep_", replicate_number))
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)



outdir <- rep_dir

outdir_alleles <- file.path(outdir, "alleles")
if (!dir.exists(outdir_alleles)) {
  dir.create(outdir_alleles, recursive = TRUE)
}
out_dir_allele_maf <- file.path(outdir, "allele_maf")
if (!dir.exists(out_dir_allele_maf)) {
  dir.create(out_dir_allele_maf, recursive = TRUE)
}

outdir_effects <- file.path(outdir, "effects")
if (!dir.exists(outdir_effects)) {
  dir.create(outdir_effects, recursive = TRUE)
}

outdir_components <- file.path(outdir, "components")
if (!dir.exists(outdir_components)) {
  dir.create(outdir_components, recursive = TRUE)
}

outdir_dif_maf_qst <- file.path(outdir, "maf_qst")
if (!dir.exists(outdir_dif_maf_qst)) {
  dir.create(outdir_dif_maf_qst, recursive = TRUE)
}

### define constants
ns <- 10  # Number of sires
nd <- 10  # Number of dams
no <- 2  # Number of offspring per pair
n_loci_qtl_list <- c(1,2,5,10,100,1000)  # Number of QTLs
maplength_quanti <- 50  # Map length for drop along ped
np <- 8

# variable replicate
replicate_number <- file_index

compute_maf_dos <- function(dosage_matrix, threshold = 0.05, return_maf = FALSE) {
  # Compute allele frequency from dosage (0, 1, 2)
  allele_freq <- colMeans(dosage_matrix, na.rm = TRUE) / 2
  
  # Compute minor allele frequency (MAF)
  maf <- pmin(allele_freq, 1 - allele_freq)
  
  # Apply filtering
  keep_loci <- maf >= threshold
  filtered_matrix <- dosage_matrix[, keep_loci, drop = FALSE]
  
  # Optional: return MAF vector too
  if (return_maf) {
    return(list(filtered_matrix = filtered_matrix, maf = maf))
  } else {
    return(filtered_matrix)
  }
}

generate_effect_matrix <- function(geno_matrix, effect_vector, h = rep(0.5, length(effect_vector))) {
  effect_matrix <- matrix(0, nrow = nrow(geno_matrix), ncol = ncol(geno_matrix))
  
  
  # Homozygous variant allele (aa): 2 * effect
  aa_idx <- which(geno_matrix == 2, arr.ind = TRUE)
  effect_matrix[aa_idx] <- 2 * effect_vector[aa_idx[, 2]]
  
  
  # Heterozygous (Aa): h * 2 * effect
  het_idx <- which(geno_matrix == 1, arr.ind = TRUE)
  effect_matrix[het_idx] <- (h[het_idx[, 2]]) * 2 * effect_vector[het_idx[, 2]]
  
  # Homozygous ancestral allele (AA): remains 0
  return(effect_matrix)
}

# Function to compute epistatic phenotype
compute_epistatic_phenotype <- function(dos_matrix, additive_effects, epi_effect_matrix) {
  n_loci <- ncol(dos_matrix)
  n_individuals <- nrow(dos_matrix)
  
  phenotype <- rep(0, n_individuals)
  
  # --- Epistatic part ---
  for (i in 1:(n_loci-1)) {
    for (j in (i+1):n_loci) {
      interaction_term <- dos_matrix[, i] * dos_matrix[, j]
      phenotype <- phenotype + interaction_term * epi_effect_matrix[i, j]
    }
  }
  
  # --- Additive part ---
  additive_contribution <- as.vector(dos_matrix %*% additive_effects)
  phenotype <- phenotype + additive_contribution
  
  return(phenotype)
}

compute_allele_freq <- function(geno_matrix) {
  allele_freqs <- apply(geno_matrix, 2, function(col) {
    sum(col, na.rm = TRUE) / (2 * sum(!is.na(col)))
  })
  return(allele_freqs)
}

compute_allele_freq_by_pop <- function(geno_matrix, pop_vector) {
  pops <- unique(pop_vector)
  
  # Create a matrix to store results
  freq_matrix <- matrix(NA, nrow = length(pops), ncol = ncol(geno_matrix))
  rownames(freq_matrix) <- pops
  colnames(freq_matrix) <- colnames(geno_matrix)
  
  for (i in seq_along(pops)) {
    pop <- pops[i]
    pop_rows <- which(pop_vector == pop)
    submatrix <- geno_matrix[pop_rows, , drop = FALSE]
    
    freq_matrix[i, ] <- apply(submatrix, 2, function(col) {
      sum(col, na.rm = TRUE) / (2 * sum(!is.na(col)))
    })
  }
  return(freq_matrix)
}

save_allele_freqs <- function(outdir_alleles, geno_matrix, pop_vector, replicate_number, prefix = "maf", n_loci_qtl = NULL, maf = FALSE) {
  # 1. Compute global and per-population allele frequencies
  global_allele_freq <- compute_allele_freq(geno_matrix)
  allele_freq_per_pop <- compute_allele_freq_by_pop(geno_matrix, pop_vector)
  
  # 2. Determine number of loci if not provided
  if (is.null(n_loci_qtl)) {
    n_loci_qtl <- ncol(geno_matrix)
  }
  
  # 3. Turn global into data.frame
  df_global <- as.data.frame(t(global_allele_freq))
  df_global$Pop <- "Global"
  df_global$Rep <- replicate_number
  
  # 4. Turn per population into data.frame
  df_pop <- as.data.frame(allele_freq_per_pop)
  df_pop$Pop <- rownames(df_pop)
  df_pop$Rep <- replicate_number
  
  # 5. Combine and reorder columns
  df_all <- rbind(df_global, df_pop)
  df_all <- df_all[, c("Rep", "Pop", setdiff(names(df_all), c("Rep", "Pop")))]
  
  # 6. File name logic
  if (maf) {
    file_name <- file.path(outdir_alleles, paste0(prefix, "_", n_loci_qtl, "_rep", replicate_number, "_allele_freq.csv"))
    dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
    write.table(df_all, file = file_name, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
  } else {
    file_name <- file.path(outdir_alleles, paste0(prefix, "_", n_loci_qtl, "_allele_freq.csv"))
    dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
    write.table(df_all, file = file_name, sep = ",", row.names = FALSE, col.names = !file.exists(file_name), append = TRUE)
  }
}

generate_epistatic_phenotypes <- function(dos_matrix, normal_effect, l_effect) {
  n_loci_qtl <- ncol(dos_matrix)
  
  if (n_loci_qtl <= 1) {
    return(list(
      epi_norm_phenotype = rep(0, nrow(dos_matrix)),
      epi_l_phenotype = rep(0, nrow(dos_matrix))
    ))
  }
  
  # 1. Create epistatic interaction matrix (0.5 everywhere except diagonal)
  epi_effect_matrix <- matrix(0.5, nrow = n_loci_qtl, ncol = n_loci_qtl)
  diag(epi_effect_matrix) <- 0  # no self-interactions
  
  # 2. Convert matrix into a long-format dataframe
  locus_pairs <- combn(1:n_loci_qtl, 2)
  epi_effect_values <- apply(locus_pairs, 2, function(pair) {
    i <- pair[1]
    j <- pair[2]
    epi_effect_matrix[i, j]
  })
  
  epi_effect_df <- data.frame(
    Locus1 = locus_pairs[1,],
    Locus2 = locus_pairs[2,],
    EpistaticEffect = epi_effect_values
  )
  
  # 3. Compute epistatic phenotypes for each type of additive effect
  epi_norm_phenotype    <- compute_epistatic_phenotype(dos_matrix, normal_effect, epi_effect_matrix)
  epi_l_phenotype       <- compute_epistatic_phenotype(dos_matrix, l_effect, epi_effect_matrix)
  
  # 4. Return everything as a named list
  return(list(
    epi_norm_phenotype       = epi_norm_phenotype,
    epi_l_phenotype          = epi_l_phenotype
  ))
}

compute_epistatic_phenotype <- function(dos_matrix, additive_effects, epi_effect_matrix) {
  n_loci <- ncol(dos_matrix)
  n_individuals <- nrow(dos_matrix)
  
  phenotype <- rep(0, n_individuals)
  
  # --- Epistatic part ---
  for (i in 1:(n_loci-1)) {
    for (j in (i+1):n_loci) {
      interaction_term <- dos_matrix[, i] * dos_matrix[, j]
      phenotype <- phenotype + interaction_term * epi_effect_matrix[i, j]
    }
  }
  
  # --- Additive part ---
  additive_contribution <- as.vector(dos_matrix %*% additive_effects)
  phenotype <- phenotype + additive_contribution
  
  return(phenotype)
}


####---- Simulation ----####
if (!file.exists(file)) stop("File not found: ", file)
## -- data import -- ##
vcf <- read.vcfR(file)
gt <- extract.gt(vcf)

# any missing value 
is_missing <- is.na(gt) | gt == "." | gt == "./." | gt == ".|."
gt_clean <- gt
gt_clean[is_missing] <- NA  # Replace missing values with NA

# Filter for only biallelic SNPs
vcf_data <- vcf@fix
biallelic <- !grepl(",", vcf_data[,"ALT"])
gt_clean <- gt_clean[biallelic, ]

# polymorphic SNPs
is_polymorphic <- function(row) {
  alleles <- unique(na.omit(row))
  length(alleles) > 1
}
poly_sites <- apply(gt_clean, 1, is_polymorphic)
gt_clean <- gt_clean[poly_sites, ]

# recode genotype
geno_map <- c("0|0" = "11", "0|1" = "12", "1|0" = "21", "1|1" = "22")
gt_recoded <- apply(gt_clean, c(1,2), function(x) {
  if (is.na(x)) return(NA)
  if (x %in% names(geno_map)) return(geno_map[[x]])
  return(NA)
})

data <- gt_recoded[complete.cases(gt_recoded), ]
data <- as.matrix(data)
data <- t(data)


## sampling loci for FST (MAF and without MAF)
sampled_loci_idx <- sample(ncol(data), size = 3000, replace = FALSE)
sampled_data <- data[,sampled_loci_idx ]
sampled_data <- apply(sampled_data, c(1, 2), as.numeric)
sampled_data <- biall2dos(sampled_data, diploid = TRUE)
founder_data <- data[,-sampled_loci_idx ]


compute_maf <- function(genotypes) {
  gts <- na.omit(genotypes)
  alleles <- unlist(strsplit(gts, split = ""))
  allele_freqs <- table(alleles) / length(alleles)
  if (length(allele_freqs) < 2) return(0)
  return(min(allele_freqs))
}

maf_values <- apply(data, 2, compute_maf)
maf_thresh <- 0.05


data_maf <- data[,maf_values >= maf_thresh]

sampled_maf_idx <- sample(ncol(data_maf), size = 3000, replace = FALSE)
data_founders_maf <- data_maf[, sampled_maf_idx]
data_founders_maf <- apply(data_founders_maf, c(1,2), as.numeric)

data_founders_maf <- biall2dos(data_founders_maf, diploid = TRUE)
founder_data <- data[,-sampled_maf_idx]

# Compare MAF Fst vs non-MAF Qst, MAF Fst vs MAF post sampling QST
# Non MAF FST vs non MAF Qst, non MAF Fst vs MAF post sampling QST

##
pop_vector <- rep(paste0("pop", 1:8), each = 20)
fst_results <- list(
  non_maf_fst = fs.dosage(sampled_data, pop_vector)$Fs[2,np+1],
  maf_fst = fs.dosage(data_founders_maf, pop_vector)$Fs[2,np+1])


#-------------------------------------------------------

#Pedigree-------
nft <- np * (ns + nd)
sire <- rep(1:ns, each = nd * no)
dam <- rep(ns + 1:nd, each = no, ns)

sire <- rep(0:(np - 1) * (nd + ns), each = (nd * ns * no)) + sire
dam <- rep(0:(np - 1) * (nd + ns), each = (nd * ns * no)) + dam

sire <- c(rep(NA, nft), sire)
dam <- c(rep(NA, nft), dam)
nt <- length(sire)
nf <- -c(1:nft)
ni <- ns + nd
pop_P <- rep(1:np, each = (ns + nd))
pop_F1 <- rep(1:np, each = (ns * nd * no))
ped <- data.frame(ind = 1:nt, sire = sire, dam = dam)
n_ind <- no * ns * nd * np #total number of individuals F1
indperpop <- 1000 #Original population size (before sampling)


#-----------------

#--------- Generate F1 generation
founder_data <- founder_data[, sample(ncol(founder_data), replace = FALSE)] # randomly shuffle 

# Split data into two half
n_loci <- ncol(founder_data)
half_loci <- floor(n_loci / 2)
neutral_data <- founder_data[, 1:half_loci]
qtl_data <- founder_data[, (half_loci + 1):n_loci]

combined_data <- cbind(neutral_data, qtl_data)
combined_data <- apply(combined_data, 2, as.numeric)

# create founders data
# Correct locus counts after final split
nl <- ncol(neutral_data)
nl_quanti <- ncol(qtl_data)
combined_nl <- nl + nl_quanti
combined_founders <- array(0, dim = c(nft, combined_nl, 2))

# No need to extract allele pairs; directly assign as founders
tmp1<-as.matrix(combined_data%/%10-1)



tmp2<-as.matrix(combined_data%%10-1)
combined_founders[,,1]<-tmp1
combined_founders[,,2]<-tmp2
# Keep genotype data in dosage format
nfounders_combined <- rbind(combined_founders[, , 1], combined_founders[, , 2])

genos.F1_combined <- drop.along.ped(ped, founders.genotypes = nfounders_combined, nloc = combined_nl, maplength = 1000)


#Extract the neutral and quantitative genotypes from the combined result
genos.F1_neutral <- genos.F1_combined[, 1:nl, ]

genos.F1_quanti <- genos.F1_combined[, (nl + 1):combined_nl, ]

dos.F1_neutral <- genos.F1_neutral[, , 1] + genos.F1_neutral[, , 2]
bed.F1_neutral <- as.bed.matrix(dos.F1_neutral)


dos.F1_quanti <- genos.F1_quanti[, , 1] + genos.F1_quanti[, , 2]
bed.F1_quanti <- as.bed.matrix(dos.F1_quanti)
#------------------------------------------------------------------------
qtl_values <- c(1, 2, 5, 10 , 100, 1000)
#Building phenotype
#Extract F1 individuals from quantitative genotype matrix
dos.Just.F1_quanti <- dos.F1_quanti[(nft+1):dim(dos.F1_quanti)[1],]


original_dos_Just_F1_quanti <- dos.Just.F1_quanti

#JEROME'S PORTION OF THE SCRIPT: 
nboot <- 1000
nrep<-nboot

#to modify if needs be
Fst<-fst_results$non_maf_fst
Fst_maf <- fst_results$maf_fst

get.rand.fst<-function(dos,pop){
  nl<-ncol(dos)
  np<-length(table(pop))
  x<-sample(nl,size=nl,replace=TRUE)
  hierfstat::fst.dosage(dos[,x],pop=pop)[np+1]
}

fst.hat <- sapply(1:nrep, function(i) get.rand.fst(sampled_data, rep(1:np, each=ni)))
fst.hat_maf <- sapply(1:nrep, function(i) get.rand.fst(data_founders_maf, rep(1:np, each=ni)))
fst_hat_list <- list( "unfiltered" = fst.hat,
               "maf_filtered" = fst.hat_maf)

dos_subset <- original_dos_Just_F1_quanti
# QST computation
for(n_loci_qtl in qtl_values){
  
  repeat {
    # Sample the QTL loci: if more loci than desired, randomly choose n_loci_qtl loci.
    if(ncol(dos_subset) >= n_loci_qtl) {
      
      sampled_cols <- sample(ncol(dos_subset), n_loci_qtl, replace = FALSE)
      
      dos_subset_loci <- dos_subset[, sampled_cols, drop = FALSE]
      
    } else if(ncol(dos_subset) < n_loci_qtl) {
      warning("Not enough QTL loci available; using all available loci.")
      # Optionally, could sample with replacement:
      # sampled_cols <- sample(ncol(dos_subset), n_loci_qtl, replace = TRUE)
      # dos_subset <- dos_subset[, sampled_cols, drop = FALSE]
    }
    
    dataf_list <- list()
    
    h_ref_dom <- rep(1, n_loci_qtl)
    h_ref_rec <- rep(0.01, n_loci_qtl)
    
    geno_matrix <- dos_subset_loci
    
    normal_effect <- rnorm(n_loci_qtl, mean = 2, sd = 1)
    l_effect <- rgamma(n_loci_qtl, shape = 0.7, scale = 1.5)
    
    # dominance phenotype
    effect_matrix_dom <- generate_effect_matrix(geno_matrix, normal_effect, h = h_ref_dom)
    effect_matrix_l <- generate_effect_matrix(geno_matrix, l_effect, h = h_ref_rec)
    
    # Fully recessive phenotype
    effect_matrix_rec <- generate_effect_matrix(geno_matrix, normal_effect, h = h_ref_rec)
    effect_matrix_l_rec <- generate_effect_matrix(geno_matrix, l_effect, h = h_ref_rec)
    
    # compute phenotypes
    dom_normal <- rowSums(effect_matrix_dom)
    dom_l <- rowSums(effect_matrix_l)
    rec_normal <- rowSums(effect_matrix_rec)
    rec_l <- rowSums(effect_matrix_l_rec)
    
    # Additive phenotype
    add_norm <- rowSums(geno_matrix * normal_effect)
    add_l <- rowSums(geno_matrix * l_effect)
    
    csv_normal <- paste0(outdir_effects, "/normal_effects_", n_loci_qtl, "_loci.csv")
    
    csv_l <- paste0(outdir_effects,"/l_effects_", n_loci_qtl, "_loci.csv")
    normal_effect_df <- t(as.data.frame(normal_effect))
    l_effect_df <- t(as.data.frame(l_effect))
    
    colnames(normal_effect_df) <- paste0("Locus_", seq_along(normal_effect))
    colnames(l_effect_df) <- paste0("Locus_", seq_along(l_effect))
    
    write.table(l_effect_df, file = csv_l, row.names = FALSE, append = file.exists(csv_l), sep = ",",
                col.names = !file.exists(csv_l))
    write.table(normal_effect_df, file = csv_normal, row.names = FALSE, append = file.exists(csv_normal), sep = ",",
                col.names = !file.exists(csv_normal))
    
    
    # Compute allele frequency for all DF and store them in csv files
    #number of populations
    populations <- unique(pop_F1)
    pop_F1_vector <- rep(paste0("Pop_", 1:8), each = 200)
    
    global_allele_freq <- compute_allele_freq(geno_matrix)
    allele_freq_per_pop <- compute_allele_freq_by_pop(geno_matrix, pop_F1_vector)
    save_allele_freqs(outdir_alleles,geno_matrix, pop_F1_vector, replicate_number = replicate_number, prefix = "unfiltered", 
                      n_loci_qtl = n_loci_qtl)
    
    if (length(unique(add_norm)) > 1) break
  }
    # Epistatic phenotype
    epi <- generate_epistatic_phenotypes(dos_matrix = geno_matrix, 
                                         normal_effect = normal_effect, 
                                         l_effect = l_effect)
    
    epi_norm_phenotype <- epi$epi_norm_phenotype
    epi_l_phenotype <- epi$epi_l_phenotype
    
    
    #DATA
    individuals <- 1:n_ind
    populations <- rep(1:np, each = ns * nd * no)
    
    sires <- rep(rep(1:ns, each = nd * no), np)
    dams <- rep(rep(1:nd, each = no), ns * np)
    
    dataf_list$add_normal <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = add_norm, Distribution = "Additive Normal" ) 
    dataf_list$add_l <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = add_l, Distribution = "Additive L" ) 
    dataf_list$dom_normal <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = dom_normal, Distribution = "Dominance Normal" ) 
    dataf_list$dom_l <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = dom_l, Distribution = "Dominance L" ) 
    dataf_list$rec_normal <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = rec_normal, Distribution = "Recessive Normal" ) 
    dataf_list$rec_l <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = rec_l, Distribution = "Recessive L" ) 
    
    if (sum(epi_norm_phenotype) != 0) {
      dataf_list$epi_norm <- data.frame(
        Individual   = individuals,
        Population   = factor(populations),
        Sire         = factor(sires),
        Dam          = factor(dams),
        Y            = epi_norm_phenotype,
        Distribution = "Epistatic Normal"
      )
    }
    
    if (sum(epi_l_phenotype) != 0) {
      dataf_list$epi_l <- data.frame(
        Individual   = individuals,
        Population   = factor(populations),
        Sire         = factor(sires),
        Dam          = factor(dams),
        Y            = epi_l_phenotype,
        Distribution = "Epistatic L"
      )
    }
    
    for(name in names(dataf_list)){
      
      
      for (dos_type in names(fst_hat_list)){
        fst.hat <- fst_hat_list[[dos_type]]
        
        if(dos_type == "unfiltered"){
          Fst = fst_results$non_maf_fst
        } else {
          Fst = fst_results$maf_fst
        }
        
        dataf <- dataf_list[[name]]
        distribution_type <- unique(dataf$Distribution)
        
        
        my.anova <- anova(aov(Y ~Population/(Sire + Dam), data = dataf))
        MSs <- my.anova[,3]
        DFs <- my.anova[,1]
        MSsire <- MSs[2]
        MSdam <- MSs[3]
        MSpop <- MSs[1]
        MSwithin <- MSs[4]
        DFsire <- DFs[2]
        DFdam <- DFs[3]
        DFpop <- DFs[1]
        DFwithin <- DFs[4]
        
        # To get MSpopNeutral
        x.hat<-fst.hat/(1-fst.hat)
        
        MSsire.hat <- MSsire/DFsire*rchisq(nrep, DFsire)
        MSdam.hat <- MSdam/DFdam*rchisq(nrep, DFdam)
        MSwithin.hat <- MSwithin/DFwithin*rchisq(nrep, DFwithin)
        
        # estimated var components dam and sire
        sigA.dam<- (MSdam-MSwithin)/(ns*no)
        sigA.sire <- (MSsire-MSwithin)/(nd*no)
        
        sigA.dam.hat <- (MSdam.hat-MSwithin.hat)/(ns*no)
        sigA.sire.hat <- (MSsire.hat-MSwithin.hat)/(nd*no)
        
        sigA.hat<-(sigA.dam.hat+sigA.sire.hat)/2
        VA.hat<-4*sigA.hat
        
        #only place to get variablity for fst I think
        MSpopNeutral<-MSwithin+sigA.dam*ns*no+sigA.sire*nd*no+4*sigA.dam*(ns*nd*no*Fst/(1-Fst))+4*sigA.sire*(nd*ns*no*Fst/(1-Fst))
        MSpopneutral.hat<-MSpopNeutral/DFpop*rchisq(nrep,df=DFpop)  
        
        sigA.pop<-(MSpop-sigA.dam*no*ns-sigA.sire*no*nd-MSwithin)/(no*ns*nd)
        
        QST<-sigA.pop/(sigA.pop+4*sigA.sire+4*sigA.dam)
        
        sigA.pop.neutral.hat<-(MSpopneutral.hat-sigA.dam.hat*no*ns-sigA.sire.hat*no*nd-MSwithin.hat)/(no*ns*nd)
        
        QST.neutral.hat<-sigA.pop.neutral.hat/(sigA.pop.neutral.hat+4*sigA.dam.hat+4*sigA.sire.hat)
        
        
        tmp.qstwg<-data.frame(fst.star=fst.hat,qstneut.star=QST.neutral.hat,sApopneut.star=sigA.pop.neutral.hat,
                              sASire.star=sigA.sire.hat,sAdam.star=sigA.dam.hat,sAwithin.star=MSwithin.hat)
        
        
     
        five.num.qstgw<-quantile(QST.neutral.hat-QST- fst.hat + Fst,c(0.025,0.25,0.5,0.75,0.975), na.rm = TRUE)
        pval.qstgw.neg<-(sum(QST.neutral.hat-fst.hat <= QST- Fst )+1)/(nrep+1)
        pval.qstgw.pos<-(sum(QST.neutral.hat - fst.hat >= QST - Fst)+1)/(nrep+1)
        
        p_value <- sum(abs(QST.neutral.hat-fst.hat) >= abs(QST- Fst)) / length(QST.neutral.hat-fst.hat)
        
        results_df <- data.frame(replicate_number = replicate_number,
                                 estimated_FST = Fst,
                                 estimated_QST = QST,
                                 p_value_pos = pval.qstgw.pos,
                                 p_value_neg = pval.qstgw.neg,
                                 p_value = p_value)
        file_csv_name <-  paste0(outdir,"/",dos_type, "_",distribution_type,"_MethodsWG_", n_loci_qtl, "_loci.csv")
        write.table(results_df, file = file_csv_name, append = TRUE, sep = ",", col.names = !file.exists(file_csv_name), row.names = FALSE)
        
      }
      
      
    }
    
    
    
}

### MAF post sampling phenotype
n_loci_qtl <- 1000
repeat {
  

  repeat {
    sampled_cols <- sample(ncol(dos_subset), n_loci_qtl, replace = FALSE)
    dos_subset_loci <- dos_subset[, sampled_cols, drop = FALSE]
    
    # MAF filtering
    maf_data <- compute_maf_dos(dos_subset_loci, return_maf = TRUE)
    maf_geno_matrix<-maf_data$filtered_matrix
    
    
    # Check how many loci remain after filtering
    if (ncol(maf_geno_matrix) > 1) break
  }
  
  nbr_loci_retained <- ncol(maf_geno_matrix)
  
  dataf_list <- list()
  
  # dominance coef
  h_ref_dom <- rep(1, nbr_loci_retained)
  h_ref_rec <- rep(0.01, nbr_loci_retained)
  
  # effect size
  normal_effect <- rnorm(nbr_loci_retained, mean = 2, sd = 1)
  l_effect <- rgamma(nbr_loci_retained, shape=0.7, scale = 1.5)
  
  # Dominance phenotype
  effect_matrix_dom <- generate_effect_matrix(maf_geno_matrix, normal_effect, h = h_ref_dom)
  effect_matrix_l <- generate_effect_matrix(maf_geno_matrix, l_effect, h = h_ref_rec)
  
  # Fully recessive phenotype
  effect_matrix_rec <- generate_effect_matrix(maf_geno_matrix, normal_effect, h = h_ref_rec)
  effect_matrix_l_rec <- generate_effect_matrix(maf_geno_matrix, l_effect, h = h_ref_rec)
  
  # compute phenotypes
  dom_normal <- rowSums(effect_matrix_dom)
  dom_l <- rowSums(effect_matrix_l)
  rec_normal <- rowSums(effect_matrix_rec)
  rec_l <- rowSums(effect_matrix_l_rec)
  
  
  # Additive phenotype
  add_norm <- rowSums(maf_geno_matrix * normal_effect)
  add_l <- rowSums(maf_geno_matrix * l_effect)
  
  # Compute allele frequency for all DF and store them in csv files
  #number of populations
  populations <- unique(pop_F1)
  pop_F1_vector <- rep(paste0("Pop_", 1:8), each = 200)
  
  
  global_allele_freq <- compute_allele_freq(maf_geno_matrix)
  allele_freq_per_pop <- compute_allele_freq_by_pop(maf_geno_matrix, pop_F1_vector)
  save_allele_freqs(out_dir_allele_maf,maf_geno_matrix, pop_F1_vector,
                    replicate_number = replicate_number,
                    prefix = paste0("maf_filtered_", ncol(maf_geno_matrix), "_loci"), 
                    n_loci_qtl = ncol(maf_geno_matrix), maf = TRUE)
  
  if (length(unique(add_norm)) > 1) break
}


# Epistatic phenotype
epi <- generate_epistatic_phenotypes(dos_matrix = maf_geno_matrix, 
                                     normal_effect = normal_effect, 
                                     l_effect = l_effect)

epi_norm_phenotype <- epi$epi_norm_phenotype
epi_l_phenotype <- epi$epi_l_phenotype

#DATA
individuals <- 1:n_ind
populations <- rep(1:np, each = ns * nd * no)

sires <- rep(rep(1:ns, each = nd * no), np)
dams <- rep(rep(1:nd, each = no), ns * np)

dataf_list$add_normal <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = add_norm, Distribution = "Additive Normal" ) 
dataf_list$add_l <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = add_l, Distribution = "Additive L" ) 
dataf_list$dom_normal <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = dom_normal, Distribution = "Dominance Normal" ) 
dataf_list$dom_l <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = dom_l, Distribution = "Dominance L" ) 
dataf_list$rec_normal <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = rec_normal, Distribution = "Recessive Normal" ) 
dataf_list$rec_l <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = rec_l, Distribution = "Recessive L" ) 

if (sum(epi_norm_phenotype) != 0) {
  dataf_list$epi_norm <- data.frame(
    Individual   = individuals,
    Population   = factor(populations),
    Sire         = factor(sires),
    Dam          = factor(dams),
    Y            = epi_norm_phenotype,
    Distribution = "Epistatic Normal"
  )
}

if (sum(epi_l_phenotype) != 0) {
  dataf_list$epi_l <- data.frame(
    Individual   = individuals,
    Population   = factor(populations),
    Sire         = factor(sires),
    Dam          = factor(dams),
    Y            = epi_l_phenotype,
    Distribution = "Epistatic L"
  )
}

for(name in names(dataf_list)){
  
  
  
  for (dos_type in names(fst_hat_list)){
    fst.hat <- fst_hat_list[[dos_type]]
    
    if(dos_type == "unfiltered"){
      Fst = fst_results$non_maf_fst
    } else {
      Fst = fst_results$maf_fst
    }
    
    dataf <- dataf_list[[name]]
    distribution_type <- unique(dataf$Distribution)
    
    
    my.anova <- anova(aov(Y ~Population/(Sire + Dam), data = dataf))
    MSs <- my.anova[,3]
    DFs <- my.anova[,1]
    MSsire <- MSs[2]
    MSdam <- MSs[3]
    MSpop <- MSs[1]
    MSwithin <- MSs[4]
    DFsire <- DFs[2]
    DFdam <- DFs[3]
    DFpop <- DFs[1]
    DFwithin <- DFs[4]
    
    # To get MSpopNeutral
    x.hat<-fst.hat/(1-fst.hat)
    
    MSsire.hat <- MSsire/DFsire*rchisq(nrep, DFsire)
    MSdam.hat <- MSdam/DFdam*rchisq(nrep, DFdam)
    MSwithin.hat <- MSwithin/DFwithin*rchisq(nrep, DFwithin)
    
    # estimated var components dam and sire
    sigA.dam<- (MSdam-MSwithin)/(ns*no)
    sigA.sire <- (MSsire-MSwithin)/(nd*no)
    
    sigA.dam.hat <- (MSdam.hat-MSwithin.hat)/(ns*no)
    sigA.sire.hat <- (MSsire.hat-MSwithin.hat)/(nd*no)
    
    sigA.hat<-(sigA.dam.hat+sigA.sire.hat)/2
    VA.hat<-4*sigA.hat
    
    #only place to get variablity for fst I think
    MSpopNeutral<-MSwithin+sigA.dam*ns*no+sigA.sire*nd*no+4*sigA.dam*(ns*nd*no*Fst/(1-Fst))+4*sigA.sire*(nd*ns*no*Fst/(1-Fst))
    MSpopneutral.hat<-MSpopNeutral/DFpop*rchisq(nrep,df=DFpop)  
    
    sigA.pop<-(MSpop-sigA.dam*no*ns-sigA.sire*no*nd-MSwithin)/(no*ns*nd)
    
    QST<-sigA.pop/(sigA.pop+4*sigA.sire+4*sigA.dam)
    
    sigA.pop.neutral.hat<-(MSpopneutral.hat-sigA.dam.hat*no*ns-sigA.sire.hat*no*nd-MSwithin.hat)/(no*ns*nd)
    
    QST.neutral.hat<-sigA.pop.neutral.hat/(sigA.pop.neutral.hat+4*sigA.dam.hat+4*sigA.sire.hat)
    
    tmp.qstwg<-data.frame(fst.star=fst.hat,qstneut.star=QST.neutral.hat,sApopneut.star=sigA.pop.neutral.hat,
                          sASire.star=sigA.sire.hat,sAdam.star=sigA.dam.hat,sAwithin.star=MSwithin.hat)
    
    
    
    five.num.qstgw<-quantile(QST.neutral.hat-QST- fst.hat + Fst,c(0.025,0.25,0.5,0.75,0.975), na.rm = TRUE)
    pval.qstgw.neg<-(sum(QST.neutral.hat-fst.hat <= QST- Fst )+1)/(nrep+1)
    pval.qstgw.pos<-(sum(QST.neutral.hat - fst.hat >= QST - Fst)+1)/(nrep+1)
    
    p_value <- sum(abs(QST.neutral.hat-fst.hat) >= abs(QST- Fst)) / length(QST.neutral.hat-fst.hat)
    
    results_df <- data.frame(replicate_number = replicate_number,
                             estimated_FST = Fst,
                             estimated_QST = QST,
                             p_value_pos = pval.qstgw.pos,
                             p_value_neg = pval.qstgw.neg,
                             p_value = p_value)
    file_csv_name <-  paste0(outdir_dif_maf_qst,"/",dos_type, "_",distribution_type,"_MethodsWG_", n_loci_qtl, "_loci.csv")
    write.table(results_df, file = file_csv_name, append = TRUE, sep = ",", col.names = !file.exists(file_csv_name), row.names = FALSE)
  }
}
