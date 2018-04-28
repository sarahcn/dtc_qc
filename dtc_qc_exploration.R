# Purpose: evaluate allele freqs and rates of missingness among potentially clinically relevant SNPs in DTC "raw" data
# SN, April 2018


##### Contents
# Identify DTC array manfiests
# Identify ACMG genes
# Add allele freq information

##### 

#####
# Identify DTC array manfiests
#####

library(dplyr)
library(readr)

# Using the table of arrays here: https://dna-explained.com/2017/04/11/autosomal-dna-transfers-which-companies-accept-which-tests/
# From Dec 2010 to Aug 2017, appears that 23andMe was using Illumina OmniExpress as a basis (w/ or w/o custom content)
# Use current Illumina Omni Express manifest as buest guess of what SNPs most 23andMe customers were tested on (recognizing it's imperfect)

# I've previously downloaded Omni Express current manifest from https://support.illumina.com/array/array_kits/humanomniexpress-24-beadchip-kit/downloads.html 

# read in manifest, w/ Illumina-provided annotation (gene name and pathogenic vs. benign, etc.)
# (run while connected to GAC network drives)
fn <- "T:/data/SNP_annotation/Illumina/HumanOmniExpress-24v1-1/HumanOmniExpress-24v1-1_A/HumanOmniExpress-24v1-1_A.annotated.txt"
ann <- read_tsv(file=fn, col_names=TRUE)
dim(ann) # 713014      8
print(ann, width=Inf)
# # A tibble: 713,014 x 8
# Name         Chr MapInfo Alleles `Transcript(s)`        `Gene(s)`     `In-exon` `Mutation(s)`
# <chr>      <int>   <int> <chr>   <chr>                  <chr>         <chr>     <chr>        
#   1 rs11063263    12  191619 [A/G]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 2 rs4980929     12  193818 [T/C]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 3 rs7975313     12  197841 [T/C]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 4 rs7294904     12  199532 [T/C]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 5 rs868249      12  207886 [T/C]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 6 rs7974165     12  216234 [A/G]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 7 rs1106983     12  218588 [A/C]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 8 rs1106984     12  218719 [T/C]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 9 rs2011738     12  219988 [A/G]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# 10 rs7135126     12  221203 [A/G]   NM_001170738,NM_015232 IQSEC3,IQSEC3 NA        Silent,Silent
# ... with 713,004 more rows

ann %>% select(`Mutation(s)`) %>% unique() %>% head(20)
# entries such as silent, synonymous, missense, etc. 
# use 'missense' as a proxy (albeit limited) for potentially clinically relevant

# add column that flags any missense vars - simplify names with parens
names(ann)[c(5:6, 8)] <- c("transcripts", "genes", "mutations") # couldn't get 'rename' to work
names(ann)
ann <- mutate(ann, missense=ifelse(grepl("Missense", mutations), TRUE, FALSE))

table(ann$missense)
# FALSE   TRUE 
# 700850  12164

#####
# Identify ACMG genes
#####

# Found a flat text file of ACMG 59 genes (from github repo of paper on ACMG penantrance)
acmg_fn <- "https://raw.githubusercontent.com/jamesdiao/2017-ACMG-penetrance/master/Supplementary_Files/download_output.txt"
acmg <- read_tsv(file=acmg_fn, col_names=TRUE)
dim(acmg); head(acmg)acmg[1] 59  6
# # A tibble: 6 x 6
#   gene  name         chrom    start      end downloaded
#   <chr> <chr>        <chr>    <int>    <int> <lgl>     
# 1 BRCA1 NM_007294    17    41196311 41277500 TRUE      
# 2 BRCA2 NM_000059    13    32889616 32973809 TRUE      
# 3 TP53  NM_000546    17     7571719  7590868 TRUE      
# 4 STK11 NM_000455    19     1205797  1228434 TRUE      
# 5 MLH1  NM_000249    3     37034840 37092337 TRUE      
# 6 MSH2  NM_001258281 2     47630205 47710367 TRUE  

length(unique(acmg$gene)) # 59
sort(acmg$gene)

# can we match up to Illumina manifest on gene name, or should we do chrom and position?
select(ann, genes) %>% unique() %>% head()
# # A tibble: 6 x 1
#   genes                                   
#   <chr>                                       
# 1 IQSEC3,IQSEC3                               
# 2 IQSEC3,IQSEC3,LOC574538                     
# 3 IQSEC3                                      
# 4 NA                                          
# 5 SLC6A12,SLC6A12,SLC6A12,SLC6A12             
# 6 SLC6A12,SLC6A12,SLC6A12,SLC6A12,LOC101929384

# multiple entries per column; we need row-wise grep
tmp <- select(ann, genes) %>% unique() %>% head()

# stingr::str_detect is vectorised over both 'string' and 'pattern'. brilliant
library(stringr)
head(stringr::str_detect(tmp, acmg$gene))

# test positive control
test <- ann %>% 
  filter(grepl("BRCA1", genes)) %>%  
  select(genes) 

# 23 rows with "BRCA1" in genes col
dim(test)
#[1] 23  1

# str_detect(string, pattern)
stringr::str_detect(test, acmg$gene)  %>% table()
# FALSE  TRUE 
#    58     1 # this maps back onto acmg$gene. we want logicals mapping back to test

stringr::str_detect(acmg$gene, test) 
# Error in UseMethod("type") : 
#   no applicable method for 'type' applied to an object of class "c('tbl_df', 'tbl', 'data.frame')"

# what I don't understand is if this can be vectorised both ways
stringr::str_detect(unlist(acmg$gene), unlist(test)) %>% table()
# FALSE 
#    59 
# Warning message:
# In stri_detect_regex(string, pattern, opts_regex = opts(pattern)) :
#   longer object length is not a multiple of shorter object length

# try to get the outcome i want with just one row
# build with for loops then see if there is some combo of apply, row-wise, vectorized fcns, etc.
tmp <- rbind(tmp, test)
tmp$acmg <- FALSE
tmp$acmg.genes <- NA
acmg$flag <- FALSE

for (i in 1:nrow(tmp)) {
  print(i)
  acmg.tmp <- acmg
  
  # loop over 59 acmg genes
  for (g in acmg.tmp$gene) {
    acmg.tmp$flag[acmg.tmp$gene %in% g] <- grepl(g, tmp[i,"genes"])
  }
  
  # if there was a match, annotate it
  if(sum(acmg.tmp$flag) > 0){
    genes.lab <- paste(acmg.tmp$gene[acmg.tmp$flag], collapse=", ")
    tmp$acmg[i] <- TRUE
    tmp$acmg.genes[i] <- genes.lab
  }
}

table(tmp$acmg)
# FALSE  TRUE 
#     6    46 

tmp %>%
  filter(acmg) %>%
  group_by(acmg.genes) %>%
  summarize(count=n())
# # A tibble: 1 x 2
#   acmg.genes count
#   <chr>      <int>
# 1 BRCA1         46
    
# beginning to think it may be better to use Granges and match on chrom and position rather than these names!
# i'm sure there's a more efficient way vs the nested for loops above
# but if i first limit Illumina annotation to the 12K missense, and maybe further by chroms with ACMG genes, 
# then it might not be to bad (from wall clock time, setting aside code efficiency) :/

# make into a function, with more generic names
flag_gene_list <- function(input, gene_list){
  
  # prepare input lists
  input$flag <- FALSE
  input$flag.genes <- NA
  gene_list$flag <- FALSE

  for (i in 1:nrow(input)) {
    # print(i)
    gene_list_tmp <- gene_list
    
    # loop over genes in gene list
    for (g in gene_list_tmp$gene) {
      gene_list_tmp$flag[gene_list_tmp$gene %in% g] <- grepl(g, input[i,"genes"])
    }
    
    # if there was a match, annotate it
    if(sum(gene_list_tmp$flag) > 0){
      genes.lab <- paste(gene_list_tmp$gene[gene_list_tmp$flag], collapse=", ")
      input$flag[i] <- TRUE
      input$flag.genes[i] <- genes.lab
    }
  }
}

# check time for 1K rows
sort(unique(acmg$chrom)) # 19 chroms. wouldn't help mutch

ann.msns <- ann %>% filter(missense)
dim(ann.msns)
# [1] 12164     9

system.time(ann.msns %>%
              head(n=1000) %>%
              flag_gene_list(gene_list=acmg)) # 82 sec. not great. 

(82*12)/60 # 16 minutes to annotate the 12K list.

# revisit str detect. what if I collapse acmg gene list down to one long string
acmg_str <- paste(acmg$gene, collapse="|")
  
str_detect(acmg_str, tmp$genes) # No, all false

# rewrite flag_gene_list as something that might work row-wise (i.e. only for loop is over ACMG list)
# gene_list needs "gene" and "flag" cols
# input needs to be gene cols
flag_gene_list_byrow <- function(input, gene_list){
  
  gene_list_tmp <- gene_list
  
  # loop over genes in gene list
  for (g in gene_list_tmp$gene) {
    gene_list_tmp$flag[gene_list_tmp$gene %in% g] <- grepl(g, input)
  }
  
  # if there was a match, annotate 
  genes.lab <- NA
  if(sum(gene_list_tmp$flag) > 0){
    genes.lab <- paste(gene_list_tmp$gene[gene_list_tmp$flag], collapse=", ")
    # input_row$flag <- TRUE
    # input_row$flag.genes <- genes.lab
    # return(input_row) # not sure how to collect this information
  }
  
  # simply return genes.lab
  return(genes.lab)
}

flag_gene_list_byrow(input=tmp[1,"genes"], gene_list=acmg) #NA
# test with a BRCA1 row
flag_gene_list_byrow(input=tmp[20,"genes"], gene_list=acmg) #BRCA1. yay!

# try rowwise
ann.msns %>%
  head(n=10) %>%
  rowwise() %>%
  mutate(acmg.gene=flag_gene_list_byrow(genes, gene_list=acmg)) %>%
  group_by(acmg.gene) %>%
  summarize(count=n())

# system.time the 1000
system.time(test.out <- ann.msns %>%
              head(n=1000) %>%
              rowwise() %>%
              mutate(acmg.gene=flag_gene_list_byrow(genes, gene_list=acmg)))
# 4 sec. MUCH IMPROVED

system.time(ann.msns <- ann %>% 
              filter(missense) %>%
              rowwise() %>%
              mutate(acmg.gene=flag_gene_list_byrow(genes, gene_list=acmg)))
# 57 sec

ann.msns %>% print(width=Inf) # 12,154 rows
ann.msns %>%
  filter(!is.na(acmg.gene)) %>% 
  group_by(missense, acmg.gene, genes) %>%
  summarize(count=n()) 

# # A tibble: 116 x 4
# # Groups:   missense, acmg.gene [?]
#    missense acmg.gene genes                           count
#    <lgl>    <chr>     <chr>                           <int>
#  1 TRUE     APC       ANAPC5,ANAPC5                       1
#  2 TRUE     APC       APC,APC,APC                         1
#  3 TRUE     APC       APCDD1                              1
#  4 TRUE     APC       APCDD1L                             1
#  5 TRUE     APC       MSH5-SAPCD1,SAPCD1                  4
#  6 TRUE     APC       MSH5,MSH5,MSH5,MSH5-SAPCD1,MSH5     6
#  7 TRUE     APC       SNAPC3                              1
#  8 TRUE     APC       SNAPC4                              1
#  9 TRUE     APOB      APOB                               12
# 10 TRUE     APOB      APOBEC2                             1

# oops. grepl is too liberal. can there be a whole word search equivalent in grep?
# or do i neet to strplit the 'genes' column?

# looks likw whole word searches might have to assume multiple entries in genes column (i.e. with trailing whitespace or commas)
# instead rewrite function to strsplit (or dplyr equivalent) genes

flag_gene_list_byrow <- function(input, gene_list){
  
  gene_list_tmp <- gene_list
  input_expand <-stringr::str_split(input, ",", simplify=TRUE) 
  
  # now I can use %in% without a for loop!
  gene_list_tmp$flag <- gene_list_tmp$gene %in% input_expand
  
  # if there was a match, annotate 
  genes.lab <- NA
  if(sum(gene_list_tmp$flag) > 0){
    genes.lab <- paste(gene_list_tmp$gene[gene_list_tmp$flag], collapse=", ")
  }
  
  # simply return genes.lab
  return(genes.lab)
}

# rerun
system.time(ann.msns <- ann %>% 
              filter(missense) %>%
              rowwise() %>%
              mutate(acmg.gene=flag_gene_list_byrow(genes, gene_list=acmg)))
# 3 sec

ann.msns %>% print(width=Inf) # 12,154 rows
ann.msns %>%
  filter(!is.na(acmg.gene)) %>% 
  group_by(missense, acmg.gene, genes) %>%
  summarize(count=n())
# # A tibble: 36 x 4
# # Groups:   missense, acmg.gene [?]
#    missense acmg.gene genes                               count
#    <lgl>    <chr>     <chr>                               <int>
#  1 TRUE     APC       APC,APC,APC                             1
#  2 TRUE     APOB      APOB                                   12
#  3 TRUE     ATP7B     ATP7B,ATP7B,ATP7B                       6
#  4 TRUE     BRCA1     BRCA1,BRCA1,BRCA1,BRCA1,BRCA1,BRCA1    12
#  5 TRUE     BRCA2     BRCA2                                  14
#  6 TRUE     CACNA1S   CACNA1S                                 1
#  7 TRUE     COL3A1    COL3A1                                  3
#  8 TRUE     DSC2      DSC2,DSC2                               1
#  9 TRUE     DSG2      DSG2                                    1
# 10 TRUE     DSG2      DSG2,LOC100652770                       1
# # ... with 26 more rows
# Warning message:
# Grouping rowwise data frame strips rowwise nature 

# from 16 minutes down to 3 sec, and now I'm actually getting the right answers

sum(!is.na(ann.msns$acmg.gene))
# 104. not sure if that will give enough power to detect differences in missingness
# remember caveat there that I'm taking Illumina gene annotation at face value, vs. searching by position

# save this interim product
write.table(ann.msns[!is.na(ann.msns$acmg.gene),], 
            file="omniex_missense_acmg.txt", row.names=FALSE, quote=FALSE, sep="\t",
            col.names=TRUE)

#####
# Add allele freq information
#####

# found MAF info from Illumina website:
# ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/ProductFiles/HumanOmniExpress-24/v1-1/humanomniexpress-24-v-1-1-population-reports-maf-copy-numbers.zip

maf_fn <- "T:/data/SNP_annotation/Illumina/HumanOmniExpress-24v1-1/HumanOmniExpress-24v1-1_A/Population Reports (MAF, Copy Numbers)/HumanOmniExpress-24v1-1_A_PopulationReport_MAF_012015.txt"
maf <- read_tsv(file=maf_fn, col_names = TRUE)

dim(maf)
# [1] 713014      8

# would really help if they specified the minor allele!
# merge in the freqs, at the very least

names(maf)
# [1] "Name"                   "Chr"                    "Position"              
# [4] "Sample Group All"       "Sample Group CEU_Unrel" "Sample Group CHB_Unrel"
# [7] "Sample Group JPT_Unrel" "Sample Group YRI_Unrel"

(names(maf) <- gsub("Sample Group ", "maf_", names(maf)))
# [1] "Name"          "Chr"           "Position"      "maf_All"       "maf_CEU_Unrel"
# [6] "maf_CHB_Unrel" "maf_JPT_Unrel" "maf_YRI_Unrel"

ann.mrg <- merge(ann.msns[!is.na(ann.msns$acmg.gene),],
                 maf[,c(1,4:8)], by="Name", all.x=TRUE, all.y=FALSE)

dim(ann.mrg); head(ann.mrg)

summary(ann.mrg$maf_CEU_Unrel)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.008929 0.086064 0.127232 0.500000 

library(ggplot2)
ggplot(ann.mrg, aes(x=maf_CEU_Unrel)) + geom_histogram() + scale_x_log10()  

with(ann.mrg, table(maf_CEU_Unrel < 0.01)) 
# FALSE  TRUE 
#    48    56 

# so over half are MAF < 1% in CEU.

# write over interim product
write.table(ann.mrg, file="omniex_missense_acmg.txt", 
            row.names=FALSE, quote=FALSE, sep="\t", col.names=TRUE)       

# put on GitHub for Bastian (and general transparency)
