##################################
# 0. PREPARACIÓN PARA EL TRABAJO #
##################################

# Especificación del directorio de trabajo (ya estoy en un proycto)
# setwd("")

# Especificación de la ruta donde se encuentran los datos
ruta<-"ruta"

# Cargar paquetes

cran_packages <- c( "tidyverse", "devtools", "knitr", "qtl", "bookdown", "grid", "gridExtra", "xtable", "kableExtra", "remotes", "viridis", "dendextend")
bioc_packages <- c("phyloseq", "dada2", "DECIPHER", "phangorn", "ggpubr", "DESeq2", "genefilter", "philr")
git_packages <- c("btools", "fantaxtic")


sapply(c(cran_packages, bioc_packages, git_packages), require, character.only = TRUE)

# set seed
set.seed(2026)



# Cargar las secuencias
SecSentido <- sort(list.files(ruta, pattern="_1.fastq"))
SecContraSentido <- sort(list.files(ruta, pattern="_2.fastq"))

# Extracción del nombre de las muestras
NombresMuestras <- sapply(strsplit(SecSentido, "_"), `[`, 1)

# Especificamos la ruta completo para evitar errores de ambigüedad
SecSentido <- file.path(ruta, SecSentido)
SecContraSentido <- file.path(ruta, SecContraSentido)


###############################################
# 1. INSPECCIÓN DE LA CALIDAD DE LAS MUESTRAS #
###############################################

PerfilCalidadSentido<- plotQualityProfile(SecSentido[1:3])
PerfilCalidadContraSentido<- plotQualityProfile(SecContraSentido[1:3])

PerfilCalidadSentido
PerfilCalidadContraSentido

############################
# 2. FILTRADO Y RECORTADO  #
############################

# Creamos un directorio para poner las reads "limpias"
ruta_filtradas <- file.path(ruta, "filtradas") 
if(!file_test("-d", ruta_filtradas)) dir.create(ruta_filtradas)
filtFs <- file.path(ruta_filtradas, paste0(NombresMuestras, "_F_filt.fastq.gz"))
filtRs <- file.path(ruta_filtradas, paste0(NombresMuestras, "_R_filt.fastq.gz"))

# Y finalmente procedemos con el control de calidad
# En Windows multithread=FALSE
out <- filterAndTrim(SecSentido, filtFs, SecContraSentido, filtRs,maxN=0, 
                     maxEE=c(2,4), 
                     truncQ=2,
                     multithread=FALSE) 

str(out)
head(out)


#####################################
# 3. DISMINUCIÓN DE LA REDUNDANCIA  #
#####################################
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Acá simplemente agregamos los nombres de las muestras al objeto de-replicado
names(derepFs) <- NombresMuestras
names(derepRs) <- NombresMuestras

#########################
# 4. MODELO DE ERRORES  #
#########################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Graficamos los errores para cada par
plotErrors(errF)
plotErrors(errR)


############################################################
# 5. INFERENCIA DE LAS "AMPLICON SEQUENCE VARIANTS (ASVs)  #
############################################################

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

##########################
# 6. UNION DE LAS READS  #
##########################

# Unimos las reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

##########################################
# 7. ENSAMBLE DE UNA TABLA DE SECUENCIAS #
##########################################
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

############################
# 8. REMOSCIÓN DE QUIMERAS #
############################
seqtabNoC <- removeBimeraDenovo(seqtabAll)
table(nchar(getSequences(seqtabNoC)))
dim(seqtabNoC)
sum(seqtabNoC)/sum(seqtabAll) # me quedé con el 43.68%

#############################
# 9. ASIGNACIÓN DE TAXONES  #
#############################
taxaTab <- assignTaxonomy(seqtabNoC, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#asignación de especie a los ASV
taxaTabExtra <- addSpecies(taxaTab, "silva_species_assignment_v132.fa.gz")
# ver las determinaciones
taxa.print <- taxaTab
taxa.print2 <- taxaTabExtra
# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

rownames(taxa.print2) <- NULL
head(taxa.print2)

#exporta el cuadro de determinaciones a un CSV
write.csv(taxa.print,"Determinaciones.csv", row.names = FALSE)
write.csv(taxa.print2,"DeterminacionesB.csv")

##############################################
# 10. CONSTRUCCION DE UN ARBOL FILOGENÉTICO  #
##############################################
# Alineación de secuencias               #
##########################################

seqs <- getSequences(seqtabNoC)

# Este paso propaga los nombres de las secuencias al árbol
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

########################################
# Construcción de un árbol de partida  #
########################################

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

capture.output(fitGTR, file = "data.csv")

#######################################
# 11. ENSAMBLADO DEL OBJETO PHYLOSEQ  #
#######################################

samdf <- read.csv("metadatos.csv", header=TRUE, row.names = 1) 
str(samdf)
rownames(seqtabNoC) %in% rownames(samdf)
all(rownames(seqtabAll) %in% rownames(samdf))

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxaTab),
               phy_tree(fitGTR$tree))


# Remover potenciales muestras sintéticas
ps <- prune_samples(sample_names(ps) != "Mock", ps) 

ps #755 Taxa

# Renombrado de los taxa
## Es más conveniente usar nombres cortos para nuestros ASV (por ejemplo, ASV21) en lugar 
## de la secuencia de ADN completa cuando se trabaja con algunas de las tablas y 
## visualizaciones de phyloseq, pero queremos mantener las secuencias de ADN completas para otros fines.
## Se almacenan las secuencias de ADN de nuestros ASV en la ranura "refseq" del objeto phyloseq y 
## luego se cambia el nombre de nuestros taxones a una cadena corta. De esa forma, los nombres breves 
## de los nuevos taxones aparecerán en tablas y diagramas, y aún podremos recuperar las secuencias 
## de ADN correspondientes a cada ASV según sea necesario con refseq(ps).

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)

### Se reemplazan las secuencias por un nombre genérico (ASV)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))  

ps

head(tax_table(ps)) %>%
  kable(format = "html", col.names = colnames(tax_table(ps))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "320px")

################################################
# 12. CONTROL DE CALIDAD DEL ANÁLISIS DEL 16S  #
################################################

prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Le agregamos la taxonomía
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
view(prevdf)
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) -> dfprev
view(dfprev)

# filtrado

# Definimos taxa a filtrar
filterPhyla = c("Armatimonadetes", "Chlamydiae", "Entotheonellaeota", "FBP", "Nitrospirae", "Verrucomicrobia", NA)

# Procedemos a filtrar
(psd1 = subset_taxa(ps, !Phylum %in% filterPhyla))

# 744 taxa

filterPhyla2 <- c("Chloroplast", "Mitochondria", "Eukaryota")
psd1 <- subset_taxa(psd1, !Kingdom %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Phylum %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Class %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Order %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Family %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Genus %in% filterPhyla2)

psd1 #705 taxa

# Filtramos taxa de acuerdo a un umbral de número medio de _read counts_, 
# en este caso 1e-5
psd2 <- filter_taxa(psd1, function(x) mean(x) > 1e-5, TRUE)
psd2 # 705

# También podemos remover taxa que no se observe más de X veces en al menos 5% de las muestras
psd3 <- filter_taxa(psd2, function(x) sum(x > 2) > (0.05*length(x)), TRUE)
psd3 # 681

# Y finalmente filtrar muestras con menos de 1000 reads
psd4 = prune_samples(sample_sums(psd3) > 1000, psd3)

psd4 # 681 sin plántula



prevdf1a = subset(prevdf, Phylum %in% get_taxa_unique(psd4, "Phylum"))
prueba<-tax_table(psd4)
prueba2<-as.data.frame(prueba)
view(prueba2)

# DEBIDO A QUE SOLO TENEMOS UNA MUESTRA DE CADA TIPO, ESTOS GRÁFICOS NO SON ÚTILE

# ggplot(prevdf1a, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
#   # Agregamos una línea para nuestro umbral
#   geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
#   scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#   facet_wrap(~Phylum) + theme(legend.position="none")



sample_sum_df <- data.frame(sum = sample_sums(psd4))

# ggplot(sample_sum_df, aes(x = sum)) + 
#   geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
#   ggtitle("Distribution of sample sequencing depth") + 
#   xlab("Read counts") +
#   theme(axis.title.y = element_blank()) 


##############################
# 13. CURVAS DE RAREFACCIÓN  #
##############################

# Primero cargamos algunos scripts de manera remota
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
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}

# Y graficamos
p <- ggrare(psd4, step = 100, color = "material", se = TRUE)
(p <- p + facet_wrap(~material))


p2<-ggrare(psd1, step = 100, color = "material", se = TRUE)



#####################
# TABLAS UTILES     #
#####################
# DE CUENTAS

# La tabla de cuentas relaciona el nombre de las taxa con las muestras y 
# con el número de reads mapeadas en contra de ellas. 
# Acá el número de reads es directamente proporcional al número de veces que 
# se observa un taxon.

head(otu_table(psd4)) %>%
  kable(format = "html", col.names = colnames(otu_table(psd4))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "350px")


# Tablas de cuentas que incluye a la plantula en formato csv
tb1<-otu_table(psd1)
view(tb1)
df_1<-as.data.frame(tb1)
write.csv(df_1, "tabla_taxones_psd1.csv")

# Tablas de cuentas que no incluye a la plantula en formato csv
tb2<-otu_table(psd4)
df_2<-as.data.frame(tb2)
write.csv(df_2, "tabla_taxones_psd4.csv")


# TABLA DE TAXONOMIA
head(tax_table(psd4)) %>%
  kable(format = "html", col.names = colnames(tax_table(psd4))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "320px")

tb3<-tax_table(psd1)
df_3<-as.data.frame(tb3)
view(df_3)


## TABLA DE METADATOS

# Esta es la metadata
as(sample_data(psd4), "data.frame") -> metad

metad %>%
  kable(format = "html", col.names = colnames(metad)) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "400px")

# ARBOL FILOGENÉTICO

plot_tree(psd4, method = "treeonly", ladderize = "left")


## RESUMEN DEL OBJETO PHYLOSEQ FILTRADO
library(microbiome)

summarize_phyloseq(psd4)

df <- psmelt(psd4)


head(df) %>%
  kable(format = "html", col.names = colnames(df)) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "400px")

view(df)
## PORCENTAJE DE TAXA
# Transformamos las cuentas en porcentaje
psd4r  = transform_sample_counts(psd4, function(x) 100*x / sum(x) )



# Filtramos las taxa con una abundancia inferior al 1%
(psd4r.filtrado = filter_taxa(psd4r, function(x) sum(x) > 0.01, TRUE))

psd5<- psd4r.filtrado

####################
# DIVERSIDAD ALFA  #
####################

# plot_richness(psd4, color = "material", x = "material", measures  = c("Observed", "Chao1", "Shannon")) + 
#   geom_bar(aes(fill = material), alpha=.7) + 
#   scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) + 
#   scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))

plot_richness(psd4, x="material")
plot_richness(psd4, color = "material", measures  = c("Observed", "Chao1", "Shannon"))

p = plot_richness(psd4, x="material", color="material", measures=c("Observed","Chao1", "Shannon"))
p + geom_point(size=5, alpha=0.7)


# devtools::install_github("david-barnett/microViz@0.9.3") # check 0.9.3 is the latest release
library(microViz)

psd4 %>%
  comp_barplot(tax_level = "Phylum",
               label = "material") + 
  coord_flip() # horizontal bars are often more readable
#> Registered S3 method overwritten by 'seriation':
#>   method         from 
#>   reorder.hclust vegan
#>   


psd4r.filtrado %>% 
  tax_transform("clr", rank = "Class") %>% 
  comp_heatmap(colors = heat_palette(sym = TRUE), name = "CLR")


# Guardamos un dataframe con las medidas de diversidad alfa
alpha_pd <- estimate_pd(psd4)
alpha_pd

alpha_pd2 <- estimate_pd(psd1)
alpha_pd2


view(otu_table(psd4))

# Combinamos la metadata con alpha.diversity
data <- cbind(sample_data(psd4), alpha_pd) 
view(data)
# # Y calculamos un ANOVA
# psd5.anova <- aov(PD ~ material, data) 
# library(xtable)
# psd5.anova.table <- xtable(psd5.anova)
# kable(psd5.anova.table, caption = "Tabla ANOVA", digits = 5, format = "markdown")


# Indices de diversidad alfa
tab <-microbiome::alpha(psd4, index = "all")
view(tab)
kable(head(tab))

head(tab) %>%
  kable(format = "html", col.names = colnames(tab), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# RIQUEZA
tab2 <- richness(psd4)
head(tab2) %>%
  kable(format = "html", col.names = colnames(tab2), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# DOMINANCIA
tab3 <- dominance(psd4, index = "all")
head(tab3) %>%
  kable(format = "html", col.names = colnames(tab3), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# RAREZA
tab4 <- rarity(psd4, index = "all")
head(tab4) %>%
  kable(format = "html", col.names = colnames(tab4), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")


p.shannon <- boxplot_alpha(psd4, 
                           index = "shannon",
                           x_var = "material",
                           fill.colors = c(Suelo="cyan4", tejidoVegetal="deeppink4"))


p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Material", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon

#####################
#  DIVERSIDAD BETA  #
#####################

OTU_filt_rar = rarefy_even_depth(otu_table(psd4), rngseed = TRUE, replace = FALSE) # rarefy the raw data using Phyloseq package
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar)) # create a separated file
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, taxaTab, samdf) # create a phyloseq object

dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = "bray")) 
dist_bc


plot_bar(psd4, fill="Phylum")