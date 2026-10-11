################################################################################
# 16S rRNA amplicon analysis (DADA2 + phyloseq) - Practice 11
# ------------------------------------------------------------------------------
# Fully commented version of "analisis1.txt".
#  - Every line of code carries a comment that explains what it does.
#  - The original Spanish object names were kept so the script stays recognisable.
#  - Lines changed with respect to the original are tagged [CHANGED];
#    lines that did not exist in the original are tagged [ADDED].
#  - Blocks that were unnumbered in the original received numbers 14-17.
#  - Lines of the form "# 1.1 Name ----" are RStudio code sections (outline).
################################################################################


##################################
# 0. WORKSPACE PREPARATION       #
##################################

# 0.1 Paths and output folders ----

# Original: set the working directory (the author was already inside an RStudio project).
# Not needed here because every path below is absolute.
# setwd("")

# [CHANGED] Root folder of the bioinformatics work
ruta_base <- "/home/samuel/Documents/science/computational-biology/bioinformatics"

# [ADDED] Raw-data folder of this practice (contains the sub-folders fasta, fastq, metadata and mothur)
ruta_raw <- file.path(ruta_base, "raw-data", "practice-11")

# [CHANGED] Folder with the raw paired-end FASTQ files (the "ruta" object of the original script)
ruta <- file.path(ruta_raw, "fastq")

# [ADDED] Folder where the generated data (filtered reads, R objects) are saved
ruta_datos <- file.path(ruta_base, "raw-data")

# [ADDED] Results folder of this practice
ruta_resultados <- file.path(ruta_base, "results", "practice-11")

# [ADDED] Sub-folder for figures (PNG)
dir_figuras <- file.path(ruta_resultados, "figuras")

# [ADDED] Sub-folder for tables (CSV)
dir_tablas <- file.path(ruta_resultados, "tablas")

# [ADDED] Sub-folder for the phylogenetic fit and tree files
dir_filogenia <- file.path(ruta_resultados, "filogenia")

# [ADDED] Folder for the generated R objects (.rds); kept with the generated data
dir_objetos <- file.path(ruta_datos, "objetos_rds")

# [CHANGED] Metadata table, located in the "metadata" sub-folder (tab-delimited mothur file)
ruta_metadatos <- file.path(ruta_raw, "metadata", "mouse.dpw.metadata")

# [CHANGED] SILVA v132 training set for genus-level classification, expected in the "fasta" sub-folder
ruta_silva_train <- file.path(ruta_raw, "database", "silva", "silva_nr_v132_train_set.fa.gz")

# [CHANGED] SILVA v132 file for exact species-level matching, expected in the "fasta" sub-folder
ruta_silva_especies <- file.path(ruta_raw, "database", "silva","silva_species_assignment_v132.fa.gz")

# [ADDED] Vector with every folder that must exist before the analysis starts
carpetas <- c(dir_figuras, dir_tablas, dir_filogenia, dir_objetos)

# [ADDED] Loop over each folder of the vector
for (d in carpetas) {
  # [ADDED] Create the folder (and any missing parent folder) only if it does not exist yet
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  # [ADDED] End of the loop
}

# [ADDED] Helper that saves a ggplot object as a PNG file inside the figures folder
guardar_fig <- function(grafico, nombre, ancho = 8, alto = 6) {
  # Build the file name, save the plot, and set size (inches) and resolution (dpi)
  ggsave(filename = file.path(dir_figuras, paste0(nombre, ".png")),
         plot = grafico, width = ancho, height = alto, dpi = 300)
  # End of the helper
}

# [ADDED] Empty table that will log the number of ASVs and samples after each filtering step
registro_filtros <- data.frame()

# [ADDED] Helper that appends one row (step name, ASVs left, samples left) to the log
registrar_paso <- function(paso, objeto) {
  # "<<-" writes to the global log; rbind adds the new row at the bottom
  registro_filtros <<- rbind(registro_filtros,
                             data.frame(Paso = paso,                  # name of the step
                                        ASVs = ntaxa(objeto),         # ASVs left in the phyloseq object
                                        Muestras = nsamples(objeto))) # samples left in the phyloseq object
  # End of the helper
}


# 0.2 Packages and random seed ----

# Ejecutar desde la raíz del proyecto, con renv activado
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
# renv::init()   # solo si el proyecto todavía no usa renv

# CRAN
cran_packages <- c(
  "tidyverse", "devtools", "knitr", "qtl", "bookdown", "gridExtra",
  "xtable", "kableExtra", "remotes", "viridis", "dendextend",
  "ggpubr", "phangorn"
)

# Bioconductor (se instalan con el prefijo "bioc::")
bioc_packages <- c(
  "phyloseq", "dada2", "DECIPHER", "DESeq2", "genefilter", "philr"
)

# GitHub (usuario/repo)
git_packages <- c("twbattaglia/btools", "gmteunisse/fantaxtic")

# Nombres de todos los paquetes tal como se cargan con library()
pkgs <- c(cran_packages, bioc_packages, basename(git_packages))

# Paquetes que faltan en la biblioteca del proyecto
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))

if (length(missing_pkgs) > 0) {
  
  # 1) INSTALACIÓN: solo lo que falta (renv configura BiocManager y repos)
  to_install <- c(
    intersect(cran_packages, missing_pkgs),
    paste0("bioc::", intersect(bioc_packages, missing_pkgs)),
    git_packages[basename(git_packages) %in% missing_pkgs]
  )
  renv::install(to_install)
  
  # La sesión queda apuntando a la carpeta temporal de renv (.renv-staging),
  # así que hay que reiniciar R antes de cargar nada.
  message(
    "\nInstalación terminada. Reinicia R (Ctrl+Shift+F10 en RStudio) ",
    "y vuelve a ejecutar este script para cargar los paquetes y guardar el renv.lock."
  )
  
} else {
  
  # 2) CARGA: todos los paquetes; FALSE indica uno que no se pudo cargar
  loaded <- sapply(pkgs, require, character.only = TRUE)
  
  if (all(loaded)) {
    # 3) Guardar versiones en renv.lock solo si todo cargó bien
    renv::snapshot(type = "all", prompt = FALSE)
  } else {
    warning("No se pudieron cargar: ", paste(names(loaded)[!loaded], collapse = ", "),
            "\nNo se actualizó renv.lock.")
  }
}

# Fix the random seed so that stochastic steps (e.g., tree optimisation) are reproducible
set.seed(2026)

# 0.3 Loading the sequence files ----

# Forward reads: file names containing "_R1_001.fastq", sorted alphabetically (so pairs stay aligned)
SecSentido       <- sort(list.files(ruta, pattern = "_R1_001.fastq"))

# Reverse reads: file names containing "_R2_001.fastq", sorted alphabetically
SecContraSentido <- sort(list.files(ruta, pattern = "_R2_001.fastq"))

# Sample names: split each forward file name at "_" and keep the first piece
NombresMuestras  <- sapply(strsplit(SecSentido, "_"), `[`, 1)

# Full path of each forward file (avoids ambiguity about where the files are)
SecSentido       <- file.path(ruta, SecSentido)

# Full path of each reverse file
SecContraSentido <- file.path(ruta, SecContraSentido)

# Safeguard: stop here with a clear message if the folder is wrong, nothing matched,
# the forward/reverse files are not paired, or a path does not point to a real file
stopifnot(dir.exists(ruta),
          length(SecSentido) > 0,
          length(SecSentido) == length(SecContraSentido),
          all(file.exists(SecSentido), file.exists(SecContraSentido)))

# Quick look at what was found (expected: "F3D0" "F3D1" "F3D141" "F3D142" "F3D143" "F3D144")
head(NombresMuestras)
head(SecSentido, 1)


###############################################
# 1. INSPECTION OF SAMPLE QUALITY             #
###############################################

# 1.1 Quality profiles ----

# Quality profile (per-cycle Phred scores) of the first three forward files
PerfilCalidadSentido<- plotQualityProfile(SecSentido[1:3])

# Quality profile of the first three reverse files
PerfilCalidadContraSentido<- plotQualityProfile(SecContraSentido[1:3])

# Display the forward profiles
PerfilCalidadSentido

# Display the reverse profiles
PerfilCalidadContraSentido

# [ADDED] Save the forward profiles as PNG
guardar_fig(PerfilCalidadSentido, "01_perfil_calidad_sentido")

# [ADDED] Save the reverse profiles as PNG
guardar_fig(PerfilCalidadContraSentido, "01_perfil_calidad_contrasentido")


############################
# 2. FILTERING AND TRIMMING #
############################

# 2.1 Quality filtering ----

# [CHANGED] Folder for the "clean" reads, now inside raw-data (generated data) instead of the FASTQ folder
ruta_filtradas <- file.path(ruta_datos, "filtradas")

# Create that folder if it does not exist ("-d" tests for a directory); [CHANGED] recursive = TRUE added
if(!file_test("-d", ruta_filtradas)) dir.create(ruta_filtradas, recursive = TRUE)

# Output file names for the filtered forward reads (sample name + "_F_filt.fastq.gz")
filtFs <- file.path(ruta_filtradas, paste0(NombresMuestras, "_F_filt.fastq.gz"))

# Output file names for the filtered reverse reads (sample name + "_R_filt.fastq.gz")
filtRs <- file.path(ruta_filtradas, paste0(NombresMuestras, "_R_filt.fastq.gz"))

# Quality control proper. On Windows multithread must be FALSE; on Linux it could be TRUE
out <- filterAndTrim(SecSentido, filtFs,       # input and output files, forward reads
                     SecContraSentido, filtRs, # input and output files, reverse reads
                     maxN=0,                   # discard reads with any ambiguous base (N)
                     maxEE=c(2,4),             # maximum expected errors: 2 (forward) and 4 (reverse)
                     truncQ=2,                 # truncate each read at the first base with quality <= 2
                     multithread=FALSE)        # single-threaded execution

# Show the structure of the summary table (reads in / reads out per file)
str(out)

# Show the first rows of that table
head(out)

# [ADDED] Save the filtering summary
write.csv(out, file.path(dir_tablas, "filtrado_lecturas.csv"))


#####################################
# 3. REDUCTION OF REDUNDANCY        #
#####################################

# 3.1 Dereplication ----

# Collapse identical forward reads into unique sequences with their abundances
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Collapse identical reverse reads into unique sequences with their abundances
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the dereplicated forward objects with the sample names
names(derepFs) <- NombresMuestras

# Name the dereplicated reverse objects with the sample names
names(derepRs) <- NombresMuestras


#########################
# 4. ERROR MODEL        #
#########################

# 4.1 Learning error rates ----

# Learn the substitution error model of the forward reads from the filtered files
errF <- learnErrors(filtFs, multithread=TRUE)

# Learn the substitution error model of the reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

# [CHANGED] Plot the forward error model (observed rates as points, fitted model as a line) and keep the plot
p_errF <- plotErrors(errF)

# [CHANGED] Plot the reverse error model and keep the plot
p_errR <- plotErrors(errR)

# Display the forward error plot
p_errF

# Display the reverse error plot
p_errR

# [ADDED] Save the forward error plot
guardar_fig(p_errF, "04_errores_sentido", ancho = 10, alto = 8)

# [ADDED] Save the reverse error plot
guardar_fig(p_errR, "04_errores_contrasentido", ancho = 10, alto = 8)


############################################################
# 5. INFERENCE OF AMPLICON SEQUENCE VARIANTS (ASVs)        #
############################################################

# 5.1 Denoising ----

# Infer forward ASVs; pool = TRUE analyses all samples jointly (more sensitive to rare variants, slower)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)

# Infer reverse ASVs with the same settings
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)


##########################
# 6. MERGING THE READS   #
##########################

# 6.1 Pair merging ----

# Merge each forward ASV with its reverse mate where the overlap is compatible
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)


##########################################
# 7. BUILDING A SEQUENCE TABLE           #
##########################################

# 7.1 Sequence table ----

# Samples x ASVs count matrix, leaving out any sample whose name contains "Mock" (synthetic control)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

# Distribution of ASV lengths (number of ASVs per sequence length in base pairs)
table(nchar(getSequences(seqtabAll)))


############################
# 8. CHIMERA REMOVAL       #
############################

# 8.1 De novo chimera removal ----

# Remove sequences that look like chimeras of two more abundant sequences (consensus across samples)
seqtabNoC <- removeBimeraDenovo(seqtabAll)

# Distribution of ASV lengths after chimera removal
table(nchar(getSequences(seqtabNoC)))

# Dimensions: samples x ASVs after chimera removal
dim(seqtabNoC)

# Fraction of total reads that survived chimera removal (original note: "I kept 43.68%")
sum(seqtabNoC)/sum(seqtabAll)

# [ADDED] Save the chimera-free sequence table (generated data)
saveRDS(seqtabNoC, file.path(dir_objetos, "seqtabNoC.rds"))


#############################
# 9. TAXONOMIC ASSIGNMENT   #
#############################

# 9.1 Classification against SILVA ----

# Assign taxonomy (Kingdom to Genus) to every ASV with the naive Bayesian classifier and the SILVA v132 training set
taxaTab <- assignTaxonomy(seqtabNoC, ruta_silva_train, multithread=TRUE)

# Add species names by exact matching (100 % identity) against the SILVA species file
taxaTabExtra <- addSpecies(taxaTab, ruta_silva_especies)

# Copy of the genus-level table to be printed
taxa.print <- taxaTab

# Copy of the table that includes species, to be printed
taxa.print2 <- taxaTabExtra

# Remove sequence row names for display only (the full DNA sequences make the table unreadable)
rownames(taxa.print) <- NULL

# Show the first rows of the genus-level table
head(taxa.print)

# Same cosmetic step for the table with species
rownames(taxa.print2) <- NULL

# Show the first rows of the table with species
head(taxa.print2)

# Export the determinations to CSV (row names were removed above, so sequences are not kept in the file)
write.csv(taxa.print, file.path(dir_tablas, "Determinaciones.csv"), row.names = FALSE)

# Export the table with species (row numbers are written because row.names is left at its default)
write.csv(taxa.print2, file.path(dir_tablas, "DeterminacionesB.csv"))

# [ADDED] Save the classification tables with the sequences as row names (generated data)
saveRDS(taxaTab, file.path(dir_objetos, "taxaTab.rds"))
saveRDS(taxaTabExtra, file.path(dir_objetos, "taxaTabExtra.rds"))

# 9.2 Completeness of the classification ----

# [ADDED] Number and percentage of ASVs that received a label at each taxonomic rank
completitud_taxonomica <- data.frame(Rango = colnames(taxaTabExtra),                     # rank names (Kingdom ... Species)
                                     ASVs_asignados = colSums(!is.na(taxaTabExtra)),     # ASVs with a non-missing label
                                     Porcentaje = round(100 * colSums(!is.na(taxaTabExtra)) / nrow(taxaTabExtra), 2), # share of all ASVs
                                     row.names = NULL)                                   # use plain row numbers

# [ADDED] Save the completeness table
write.csv(completitud_taxonomica, file.path(dir_tablas, "completitud_taxonomica.csv"), row.names = FALSE)


##############################################
# 10. BUILDING A PHYLOGENETIC TREE           #
##############################################

# 10.1 Sequence alignment ----

# Extract the ASV sequences (the column names of the sequence table)
seqs <- getSequences(seqtabNoC)

# Use the sequences as their own names so the names propagate to the tree
names(seqs) <- seqs

# Multiple sequence alignment with DECIPHER (no anchoring, silent)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

# 10.2 Starting tree and maximum-likelihood fit ----

# Convert the alignment to a phangorn "phyDat" object
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")

# Pairwise maximum-likelihood distances between sequences
dm <- dist.ml(phangAlign)

# Neighbour-joining tree used as the starting tree
treeNJ <- NJ(dm)

# Likelihood object for that tree and alignment (default substitution model)
fit = pml(treeNJ, data=phangAlign)

# Update the model: 4 gamma rate categories and a starting proportion of invariant sites of 0.2
fitGTR <- update(fit, k=4, inv=0.2)

# Optimise the GTR model (rates, invariant sites, gamma shape) and the tree topology (stochastic search)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# [CHANGED] Save the printed summary of the fit (plain text, hence .txt instead of the original data.csv)
capture.output(fitGTR, file = file.path(dir_filogenia, "ajuste_GTR.txt"))

# [ADDED] Save the complete fitted object (generated data)
saveRDS(fitGTR, file.path(dir_objetos, "fitGTR.rds"))


#######################################
# 11. ASSEMBLY OF THE PHYLOSEQ OBJECT #
#######################################

# 11.1 Phyloseq object ----

# [CHANGED] Read the sample metadata; the first column provides the row names (sample names)
samdf <- read.csv(ruta_metadatos, header=TRUE, row.names = 1)

# Show the structure of the metadata
str(samdf)

# Check which sequence-table samples have metadata (TRUE/FALSE per sample)
rownames(seqtabNoC) %in% rownames(samdf)

# Check that every sample of the original table has metadata (single TRUE/FALSE)
all(rownames(seqtabAll) %in% rownames(samdf))

# Build the phyloseq object from the four components
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), # counts (samples in rows, ASVs in columns)
               sample_data(samdf),                        # sample metadata
               tax_table(taxaTab),                        # taxonomy (genus level; species not included)
               phy_tree(fitGTR$tree))                     # maximum-likelihood tree

# Safeguard: drop any sample called "Mock" (it was already excluded when the sequence table was built)
ps <- prune_samples(sample_names(ps) != "Mock", ps)

# Print a summary of the object (original note: 755 taxa)
ps

# 11.2 Renaming the ASVs ----
# Short names (e.g., ASV21) are easier to use in tables and plots than the full DNA sequence.
# The sequences are stored in the "refseq" slot, so refseq(ps) always recovers them.

# Create a DNA string set whose elements are the current taxa names (the sequences)
dna <- Biostrings::DNAStringSet(taxa_names(ps))

# Name each sequence with itself
names(dna) <- taxa_names(ps)

# Attach the sequences to the phyloseq object (refseq slot)
ps <- merge_phyloseq(ps, dna)

# Replace the sequence names with short generic names ASV1, ASV2, ...
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Print the object again to confirm the change
ps

# [ADDED] Log the starting point of the filtering cascade
registrar_paso("Phyloseq object (all ASVs)", ps)

# [ADDED] Save the complete phyloseq object (generated data)
saveRDS(ps, file.path(dir_objetos, "ps_completo.rds"))

# First rows of the taxonomy table, piped into a styled, scrollable HTML table
head(tax_table(ps)) %>%
  kable(format = "html", col.names = colnames(tax_table(ps))) %>%  # build the HTML table with the rank names as headers
  kable_styling() %>%                                              # apply the default styling
  kableExtra::scroll_box(width = "100%", height = "320px")         # make it scrollable


################################################
# 12. QUALITY CONTROL OF THE 16S ANALYSIS      #
################################################

# 12.1 Prevalence ----

# For each ASV, count in how many samples it is present (count > 0)
prevdf = apply(X = otu_table(ps),                                  # counts table
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),# 1 = rows if ASVs are rows, otherwise 2 = columns
               FUN = function(x){sum(x > 0)})                      # number of samples with at least one read

# Add total abundance and taxonomy to the prevalence values
prevdf = data.frame(Prevalence = prevdf,                           # samples in which each ASV occurs
                    TotalAbundance = taxa_sums(ps),                # total reads of each ASV
                    tax_table(ps))                                 # taxonomy of each ASV

# Open the table in the RStudio viewer ([CHANGED] only when running interactively)
if (interactive()) view(prevdf)

# Mean and summed prevalence per phylum (ddply splits by Phylum and applies the function to each piece)
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) -> dfprev

# Open the phylum-level table ([CHANGED] only when running interactively)
if (interactive()) view(dfprev)

# [ADDED] Save the prevalence tables
write.csv(prevdf, file.path(dir_tablas, "prevalencia_ASV.csv"))
write.csv(dfprev, file.path(dir_tablas, "prevalencia_filo.csv"), row.names = FALSE)

# 12.2 Filtering cascade ----

# Phyla to remove (rare phyla) plus NA (ASVs without phylum assignment)
filterPhyla = c("Armatimonadetes", "Chlamydiae", "Entotheonellaeota", "FBP", "Nitrospirae", "Verrucomicrobia", NA)

# Keep only ASVs whose phylum is not in that list (the outer parentheses also print the result; original note: 744 taxa)
(psd1 = subset_taxa(ps, !Phylum %in% filterPhyla))

# [ADDED] Log this step
registrar_paso("Removal of rare and unassigned phyla", psd1)

# Labels that identify organelle or eukaryotic sequences (not bacterial/archaeal 16S)
filterPhyla2 <- c("Chloroplast", "Mitochondria", "Eukaryota")

# Remove ASVs labelled as such at the Kingdom rank
psd1 <- subset_taxa(psd1, !Kingdom %in% filterPhyla2)

# ... at the Phylum rank
psd1 <- subset_taxa(psd1, !Phylum %in% filterPhyla2)

# ... at the Class rank
psd1 <- subset_taxa(psd1, !Class %in% filterPhyla2)

# ... at the Order rank
psd1 <- subset_taxa(psd1, !Order %in% filterPhyla2)

# ... at the Family rank
psd1 <- subset_taxa(psd1, !Family %in% filterPhyla2)

# ... at the Genus rank
psd1 <- subset_taxa(psd1, !Genus %in% filterPhyla2)

# Print the object (original note: 705 taxa)
psd1

# [ADDED] Log this step
registrar_paso("Removal of chloroplast, mitochondria and eukaryotes", psd1)

# Keep ASVs whose mean count across samples is above 1e-5 (a very permissive abundance threshold)
psd2 <- filter_taxa(psd1, function(x) mean(x) > 1e-5, TRUE)

# Print the object (original note: 705)
psd2

# [ADDED] Log this step
registrar_paso("Mean-abundance filter (> 1e-5)", psd2)

# Keep ASVs seen more than 2 times in more than 5 % of the samples
psd3 <- filter_taxa(psd2, function(x) sum(x > 2) > (0.05*length(x)), TRUE)

# Print the object (original note: 681)
psd3

# [ADDED] Log this step
registrar_paso("Prevalence filter (> 2 reads in > 5 % of samples)", psd3)

# Keep only samples with more than 1000 reads
psd4 = prune_samples(sample_sums(psd3) > 1000, psd3)

# Print the object (original note: 681 taxa, without the seedling sample)
psd4

# [ADDED] Log this step
registrar_paso("Removal of samples with <= 1000 reads", psd4)

# [ADDED] Save the filtering log
write.csv(registro_filtros, file.path(dir_tablas, "registro_filtrado.csv"), row.names = FALSE)

# [ADDED] Save the filtered objects (generated data): psd1 keeps all samples, psd4 only those above 1000 reads
saveRDS(psd1, file.path(dir_objetos, "psd1.rds"))
saveRDS(psd4, file.path(dir_objetos, "psd4.rds"))

# Restrict the prevalence table to the phyla still present in psd4
prevdf1a = subset(prevdf, Phylum %in% get_taxa_unique(psd4, "Phylum"))

# Extract the taxonomy table of psd4
prueba<-tax_table(psd4)

# Convert it to a data frame
prueba2<-as.data.frame(prueba)

# Open it in the viewer ([CHANGED] only when running interactively)
if (interactive()) view(prueba2)

# Original note: the following plots are not informative because there is only one sample per material.
# Kept commented, exactly as in the original:
# ggplot(prevdf1a, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
#   # Agregamos una línea para nuestro umbral
#   geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
#   scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#   facet_wrap(~Phylum) + theme(legend.position="none")

# Data frame with the sequencing depth (total reads) of each retained sample
sample_sum_df <- data.frame(sum = sample_sums(psd4))

# Histogram of sequencing depth, also commented out in the original:
# ggplot(sample_sum_df, aes(x = sum)) +
#   geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
#   ggtitle("Distribution of sample sequencing depth") +
#   xlab("Read counts") +
#   theme(axis.title.y = element_blank())


##############################
# 13. RAREFACTION CURVES     #
##############################

# 13.1 Rarefaction curves ----

# Names of the helper scripts (from the phyloseq-extended repository) that are loaded remotely
scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")

# Full GitHub (raw) URL of each script
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

# Loop over the URLs (requires an internet connection)
for (url in urls) {
  # Download and run the script, which defines its functions in the global environment (e.g., ggrare)
  source(url)
  # End of the loop
}

# Rarefaction curves of psd4: richness as a function of the number of reads drawn (step of 100 reads), coloured by material
p <- ggrare(psd4, step = 100, color = "material", se = TRUE)

# Split the plot into one panel per material and print it
(p <- p + facet_wrap(~material))

# Same curves for psd1 (the object that still includes the sample below 1000 reads)
p2<-ggrare(psd1, step = 100, color = "material", se = TRUE)

# [ADDED] Save both rarefaction plots
guardar_fig(p, "13_rarefaccion_psd4")
guardar_fig(p2, "13_rarefaccion_psd1")


#####################
# 14. USEFUL TABLES #
#####################

# 14.1 Count, taxonomy and metadata tables ----

# COUNT TABLE: relates ASVs to samples and to the number of reads assigned to each.
# Here the number of reads is proportional to how often an ASV was observed.

# First rows of the psd4 count table as a styled, scrollable HTML table
head(otu_table(psd4)) %>%
  kable(format = "html", col.names = colnames(otu_table(psd4))) %>%  # HTML table with ASV names as headers
  kable_styling() %>%                                                # default styling
  kableExtra::scroll_box(width = "100%", height = "350px")           # scrollable box

# Count table that INCLUDES the seedling sample (psd1)
tb1<-otu_table(psd1)

# Open it in the viewer ([CHANGED] only when running interactively)
if (interactive()) view(tb1)

# Convert it to a data frame
df_1<-as.data.frame(tb1)

# [CHANGED] Save it as CSV in the tables folder
write.csv(df_1, file.path(dir_tablas, "tabla_taxones_psd1.csv"))

# Count table that EXCLUDES the seedling sample (psd4)
tb2<-otu_table(psd4)

# Convert it to a data frame
df_2<-as.data.frame(tb2)

# [CHANGED] Save it as CSV in the tables folder
write.csv(df_2, file.path(dir_tablas, "tabla_taxones_psd4.csv"))

# TAXONOMY TABLE: first rows for psd4 as a styled, scrollable HTML table
head(tax_table(psd4)) %>%
  kable(format = "html", col.names = colnames(tax_table(psd4))) %>%  # HTML table with rank names as headers
  kable_styling() %>%                                                # default styling
  kableExtra::scroll_box(width = "100%", height = "320px")           # scrollable box

# Taxonomy table of psd1
tb3<-tax_table(psd1)

# Convert it to a data frame
df_3<-as.data.frame(tb3)

# Open it in the viewer ([CHANGED] only when running interactively)
if (interactive()) view(df_3)

# [ADDED] Save the taxonomy of psd1 as CSV
write.csv(df_3, file.path(dir_tablas, "tabla_taxonomia_psd1.csv"))

# METADATA TABLE: extract the sample data of psd4 as a plain data frame
as(sample_data(psd4), "data.frame") -> metad

# Show the metadata as a styled, scrollable HTML table
metad %>%
  kable(format = "html", col.names = colnames(metad)) %>%   # HTML table with variable names as headers
  kable_styling() %>%                                       # default styling
  kableExtra::scroll_box(width = "100%", height = "400px")  # scrollable box

# PHYLOGENETIC TREE: tree only (no sample annotations), branches ordered "ladderized" to the left
p_arbol <- plot_tree(psd4, method = "treeonly", ladderize = "left")

# Display the tree
p_arbol

# [ADDED] Save the tree figure (tall page because of the number of tips)
guardar_fig(p_arbol, "14_arbol_filogenetico_psd4", ancho = 8, alto = 12)

# SUMMARY OF THE FILTERED PHYLOSEQ OBJECT: load the microbiome package
library(microbiome)

# Summary statistics (reads per sample, ASVs, sparsity, etc.) of psd4
summarize_phyloseq(psd4)

# Melt the phyloseq object into one long table (one row per ASV x sample) with counts, taxonomy and metadata
df <- psmelt(psd4)

# First rows of the long table as a styled, scrollable HTML table
head(df) %>%
  kable(format = "html", col.names = colnames(df)) %>%      # HTML table with column names as headers
  kable_styling() %>%                                       # default styling
  kableExtra::scroll_box(width = "100%", height = "400px")  # scrollable box

# Open the long table in the viewer ([CHANGED] only when running interactively)
if (interactive()) view(df)

# [ADDED] Save the long table
write.csv(df, file.path(dir_tablas, "psd4_formato_largo.csv"), row.names = FALSE)


###########################
# 15. RELATIVE ABUNDANCE  #
###########################

# 15.1 Percentage transformation ----

# Convert the counts of each sample into percentages (each sample sums to 100)
psd4r  = transform_sample_counts(psd4, function(x) 100*x / sum(x) )

# Keep ASVs whose summed relative abundance (in percentage points, across samples) is above 0.01.
# Original comment says "below 1 %", but with data in percentages the threshold 0.01 means 0.01 %, not 1 %.
(psd4r.filtrado = filter_taxa(psd4r, function(x) sum(x) > 0.01, TRUE))

# Shorter alias for the filtered relative-abundance object
psd5<- psd4r.filtrado

# [ADDED] Save the relative-abundance object (generated data)
saveRDS(psd4r.filtrado, file.path(dir_objetos, "psd4r_filtrado.rds"))


####################
# 16. ALPHA DIVERSITY #
####################

# 16.1 Richness plots (phyloseq) ----

# Commented out in the original (bar version with custom colours):
# plot_richness(psd4, color = "material", x = "material", measures  = c("Observed", "Chao1", "Shannon")) +
#   geom_bar(aes(fill = material), alpha=.7) +
#   scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) +
#   scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))

# [CHANGED] Default richness plot (all phyloseq measures) by material; the plot is kept in an object
p_rich1 <- plot_richness(psd4, x="material")

# Display it
p_rich1

# [CHANGED] Observed richness, Chao1 and Shannon, coloured by material; the plot is kept in an object
p_rich2 <- plot_richness(psd4, color = "material", measures  = c("Observed", "Chao1", "Shannon"))

# Display it
p_rich2

# Same three measures with the material on the x axis and as colour
p = plot_richness(psd4, x="material", color="material", measures=c("Observed","Chao1", "Shannon"))

# [CHANGED] Enlarge and soften the points; the result is kept in an object
p_rich3 <- p + geom_point(size=5, alpha=0.7)

# Display it
p_rich3

# [ADDED] Save the three richness plots
guardar_fig(p_rich1, "16_riqueza_todas_medidas")
guardar_fig(p_rich2, "16_riqueza_observada_chao1_shannon")
guardar_fig(p_rich3, "16_riqueza_puntos_por_material")

# 16.2 Taxonomic composition plots (microViz) ----
# In the original these plots sit inside the alpha-diversity block; they describe composition.

# Install line kept as a reminder (check that 0.9.3 is the latest release):
# devtools::install_github("david-barnett/microViz@0.9.3")
# Load microViz (comp_barplot, tax_transform, comp_heatmap)
library(microViz)

# [CHANGED] Bar plot of the composition at Phylum level, labelled by material, with horizontal bars (more readable)
p_comp <- psd4 %>%
  comp_barplot(tax_level = "Phylum",   # taxonomic rank to display
               label = "material") +   # sample label shown on the axis
  coord_flip()                         # horizontal bars

# Display it
p_comp

# [ADDED] Save the composition plot
guardar_fig(p_comp, "16_composicion_filo")

# [CHANGED] Aggregate to Class, apply the centred log-ratio (CLR) transformation and draw a heatmap; keep the object
hm_clr <- psd4r.filtrado %>%
  tax_transform("clr", rank = "Class") %>%                                # CLR transformation at Class level
  comp_heatmap(colors = heat_palette(sym = TRUE), name = "CLR")           # symmetric colour palette, legend titled "CLR"

# Display the heatmap
hm_clr

# [ADDED] Save the heatmap (not a ggplot object, so a PNG device is opened explicitly)
png(file.path(dir_figuras, "16_heatmap_clr_clase.png"), width = 2400, height = 1800, res = 300)

# [ADDED] Draw the heatmap on that device
print(hm_clr)

# [ADDED] Close the device so the file is written (invisible() suppresses the "null device" message)
invisible(dev.off())

# 16.3 Alpha-diversity tables ----

# Phylogenetic diversity (Faith's PD) and species richness (SR) per sample, from the btools package (psd4)
alpha_pd <- estimate_pd(psd4)

# Print it
alpha_pd

# Same measures for psd1 (includes the seedling sample)
alpha_pd2 <- estimate_pd(psd1)

# Print it
alpha_pd2

# Open the psd4 count table in the viewer ([CHANGED] only when running interactively)
if (interactive()) view(otu_table(psd4))

# Combine the sample metadata with the alpha-diversity values
data <- cbind(sample_data(psd4), alpha_pd)

# Open the combined table ([CHANGED] only when running interactively)
if (interactive()) view(data)

# Commented out in the original (ANOVA of PD by material):
# # Y calculamos un ANOVA
# psd5.anova <- aov(PD ~ material, data)
# library(xtable)
# psd5.anova.table <- xtable(psd5.anova)
# kable(psd5.anova.table, caption = "Tabla ANOVA", digits = 5, format = "markdown")

# All alpha-diversity indices available in the microbiome package
tab <-microbiome::alpha(psd4, index = "all")

# Open the table ([CHANGED] only when running interactively)
if (interactive()) view(tab)

# Quick text version of the first rows
kable(head(tab))

# First rows as a styled, scrollable HTML table (2 decimals)
head(tab) %>%
  kable(format = "html", col.names = colnames(tab), digits = 2) %>%  # HTML table, 2 decimals
  kable_styling() %>%                                                # default styling
  kableExtra::scroll_box(width = "100%", height = "300px")           # scrollable box

# RICHNESS indices (microbiome package)
tab2 <- richness(psd4)

# First rows as a styled, scrollable HTML table
head(tab2) %>%
  kable(format = "html", col.names = colnames(tab2), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# DOMINANCE indices (microbiome package)
tab3 <- dominance(psd4, index = "all")

# First rows as a styled, scrollable HTML table
head(tab3) %>%
  kable(format = "html", col.names = colnames(tab3), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# RARITY indices (microbiome package)
tab4 <- rarity(psd4, index = "all")

# First rows as a styled, scrollable HTML table
head(tab4) %>%
  kable(format = "html", col.names = colnames(tab4), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# [ADDED] Save the alpha-diversity tables
write.csv(alpha_pd, file.path(dir_tablas, "alfa_pd_psd4.csv"))
write.csv(alpha_pd2, file.path(dir_tablas, "alfa_pd_psd1.csv"))
write.csv(tab, file.path(dir_tablas, "alfa_microbiome_todos.csv"))
write.csv(tab2, file.path(dir_tablas, "alfa_riqueza.csv"))
write.csv(tab3, file.path(dir_tablas, "alfa_dominancia.csv"))
write.csv(tab4, file.path(dir_tablas, "alfa_rareza.csv"))

# 16.4 Shannon box plot ----

# [ADDED] boxplot_alpha() is provided by the microbiomeutilities package (load it if it is installed)
require(microbiomeutilities)

# Box plot of the Shannon index by material, with one fill colour per material
p.shannon <- boxplot_alpha(psd4,
                           index = "shannon",                                           # diversity index to plot
                           x_var = "material",                                          # grouping variable
                           fill.colors = c(Suelo="cyan4", tejidoVegetal="deeppink4"))   # soil and plant tissue colours

# Apply a minimal theme and larger text to the plot
p.shannon <- p.shannon + theme_minimal() +
  labs(x="Material", y="Shannon diversity") +    # axis titles
  theme(axis.text = element_text(size=12),       # axis tick labels
        axis.title = element_text(size=16),      # axis titles
        legend.text = element_text(size=12),     # legend labels
        legend.title = element_text(size=16))    # legend title

# Display the final plot
p.shannon

# [ADDED] Save the Shannon box plot
guardar_fig(p.shannon, "16_shannon_boxplot", ancho = 6, alto = 5)


#####################
#  17. BETA DIVERSITY #
#####################

# 17.1 Bray-Curtis dissimilarity ----

# Rarefy the raw counts to an even depth (default: the smallest library) without replacement; rngseed = TRUE fixes the seed (TRUE acts as the value 1)
OTU_filt_rar = rarefy_even_depth(otu_table(psd4), rngseed = TRUE, replace = FALSE)

# Store the rarefied counts as an ordinary data frame (samples in rows)
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar))

# [CHANGED] Build a phyloseq object with the rarefied counts; the original used taxaTab and samdf, whose ASV names
# (full sequences) no longer match the renamed ASVs, so the taxonomy and metadata of psd4 are used instead
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, tax_table(psd4), sample_data(psd4))

# [ADDED] Load vegan, which provides vegdist()
library(vegan)

# Bray-Curtis dissimilarity between samples, converted from a "dist" object to a full matrix
dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = "bray"))

# Print the matrix
dist_bc

# [ADDED] Save the matrix and the rarefied object (generated data)
write.csv(dist_bc, file.path(dir_tablas, "distancia_braycurtis.csv"))
saveRDS(data_phylo_filt_rar, file.path(dir_objetos, "psd4_rarefied.rds"))

# 17.2 Community bar plot ----

# Stacked bar plot of the counts of psd4, filled by phylum
p_bar <- plot_bar(psd4, fill="Phylum")

# Display it
p_bar

# [ADDED] Save the bar plot
guardar_fig(p_bar, "17_barras_filo")


#####################
# SESSION INFO      #
#####################

# [ADDED] Save the R and package versions used, for reproducibility
capture.output(sessionInfo(), file = file.path(ruta_resultados, "sessionInfo.txt"))