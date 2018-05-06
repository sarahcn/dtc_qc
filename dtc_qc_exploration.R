# Purpose: evaluate allele freqs and rates of missingness among potentially clinically relevant SNPs in DTC "raw" data
# SN, April 2018


##### Contents
# Identify DTC array manfiests
# Identify ACMG genes
# Add allele freq information
# Incorporate openSNP summary data - missingness

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

#####
# Incorporate openSNP summary data - missingness
#####
# 5/5/18

library(readr)

# Bastian parsed openSNP raw data for the 104 selected snps (json format)

# my first time dealing with json format - see if this will work
# install.packages("rjson")
library(rjson)
dat <- rjson::fromJSON(file="opensnp_genotypes_processed.json")
class(dat) # list
names(dat)
# [1] "rs3219497"  "rs3219484"  "rs11591147" "rs11583680" "rs562556"   "rs28362277"
# [7] "rs505151"   "rs3850625"  "rs4523540"  "rs3730238"  "rs3766871"  "rs34967813"
# [13] "rs12742148" "rs1042034"  "rs1042031"  "rs1801703"  "rs1801701"  "rs676210"  
# [19] "rs533617"   "rs568413"   "rs13306187" "rs12713843" "rs679899"   "rs1367117" 
# [25] "rs4987188"  "rs2020908"  "rs2020912"  "rs1800255"  "rs2271683"  "rs1516446" 
# [31] "rs2340917"  "rs1799977"  "rs2020873"  "rs7626962"  "rs1805124"  "rs459552"  
# [37] "rs17604693" "rs2076299"  "rs28763966" "rs6929069"  "rs34738426" "rs1805123" 
# [43] "rs1073123"  "rs1799939"  "rs3026760"  "rs17158558" "rs11570112" "rs11570060"
# [49] "rs3729989"  "rs2959656"  "rs2071312"  "rs1046116"  "rs2428140"  "rs766173"  
# [55] "rs144848"   "rs1799944"  "rs1799951"  "rs4987048"  "rs4987117"  "rs1799954" 
# [61] "rs11571660" "rs28897743" "rs169547"   "rs169548"   "rs11571746" "rs4987047" 
# [67] "rs1801426"  "rs17071686" "rs3092905"  "rs7334118"  "rs732774"   "rs1061472" 
# [73] "rs2277447"  "rs1801243"  "rs363821"   "rs140647"   "rs140597"   "rs140586"  
# [79] "rs12324002" "rs25403"    "rs16967494" "rs17882252" "rs28897695" "rs16942"   
# [85] "rs2227945"  "rs4986852"  "rs16941"    "rs4986848"  "rs799917"   "rs1800709" 
# [91] "rs4986850"  "rs1799950"  "rs28897674" "rs28897673" "rs1893963"  "rs2230234" 
# [97] "rs2278792"  "rs13306510" "rs1800321"  "rs12392549" "rs1864423"  "rs3729986" 
# [103] "VG13S52444" "rs1800328" 

# nice! the 104 SNPs
# look at one snp
names(dat[[1]])
# [1] "absolute_values"        "relative_values"        "number_of_observations"

dat[[1]]
# $`absolute_values`
# $`absolute_values`$`CC`
# [1] 1015
# 
# $`absolute_values`$C
# [1] 7
# 
# $`absolute_values`$CT
# [1] 2
# 
# $`absolute_values`$GG
# [1] 2
# 
# 
# $relative_values
# $relative_values$`CC`
# [1] 0.9892788
# 
# $relative_values$C
# [1] 0.006822612
# 
# $relative_values$CT
# [1] 0.001949318
# 
# $relative_values$GG
# [1] 0.001949318
# 
# 
# $number_of_observations
# [1] 1026

# how can we calcualte missingness?
dat[[2]] # oh, now I see there are some snps that have a relative count of "--"

# also looks like we need to deal with strand issues. throw out where genotype appears to be haploid?

# read in my original table
acmg <- read_tsv("omniex_missense_acmg.txt", col_names = TRUE)
dim(acmg); names(acmg)
# [1] 104  15
# [1] "Name"          "Chr"           "MapInfo"       "Alleles"       "transcripts"  
# [6] "genes"         "In-exon"       "mutations"     "missense"      "acmg.gene"    
# [11] "maf_All"       "maf_CEU_Unrel" "maf_CHB_Unrel" "maf_JPT_Unrel" "maf_YRI_Unrel"

# how to extract the named element "$relative_values$`--`" from each list item from json
head(unlist(dat))
# rs3219497.absolute_values.CC  rs3219497.absolute_values.C rs3219497.absolute_values.CT 
# 1.015000e+03                 7.000000e+00                 2.000000e+00 
# rs3219497.absolute_values.GG rs3219497.relative_values.CC  rs3219497.relative_values.C 
# 2.000000e+00                 9.892788e-01                 6.822612e-03

# looks like i can just unlist and grep for name
dat.unlst <- unlist(dat)
head(dat.unlst[grep("relative_values.--", names(dat.unlst))]) # yup, that's what I want

miss <- dat.unlst[grep("relative_values.--", names(dat.unlst))]

rs.nms <- gsub(".relative_values.--", "", names(miss)); head(rs.nms)
miss.freq <- data.frame(rsID=rs.nms, miss.frac=miss)
dim(miss.freq); head(miss.freq)
# 73 2

summary(miss.freq$miss.frac)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003388 0.0009881 0.0027192 0.0121796 0.0060362 0.3538406

# overall, the MCRs look reasonable. 
sum(miss.freq$miss.frac > 0.05) # only 3 with MCR > 5%

# tack onto acmg genes
acmg$opensnp.miss.frac <- miss.freq$miss.frac[match(acmg$Name, miss.freq$rsID)]
with(acmg, table(acmg.gene, is.na(opensnp.miss.frac)))
# acmg.gene FALSE TRUE
# APC         1    0
# APOB        4    8
# ATP7B       3    3
# BRCA1      10    2
# BRCA2      11    3
# CACNA1S     1    0
# COL3A1      1    2
# DSC2        1    0
# DSG2        2    0
# DSP         4    1
# FBN1        4    2
# GLA         1    0
# KCNH2       1    0
# LDLR        0    1
# MEN1        2    0
# MLH1        1    1
# MSH2        1    0
# MSH6        2    0
# MUTYH       1    1
# MYBPC3      4    0
# MYH11       0    1
# MYL2        1    0
# OTC         1    1
# PCSK9       4    1
# PKP2        1    0
# RB1         2    0
# RET         2    1
# RYR2        1    2
# SCN5A       2    0
# TMEM43      1    0
# TNNT2       1    1
# TP53        1    0
# TSC1        1    0

# unclear if we should say the other 104 - 73 are missing=0 or NA. 
# safer to say NA until I confirm with Bastian the other interpretation
# I'm also just noticing that APOL1 isn't on ACMG list

# make a histogram
library(ggplot2)
ggplot(acmg, aes(x=opensnp.miss.frac)) + 
  geom_histogram(binwidth=0.01, fill="darkblue", colour="white") +
  theme_bw() + 
  ggtitle("Fraction of missing calls, over openSNP users")
ggsave("frac_misscall.png")

# before saving annotated list of SNPs, adding in the number of openSNP observations
obs <- dat.unlst[grepl("number_of_observations", names(dat.unlst))]

rs.nms <- gsub(".number_of_observations", "", names(obs)); head(rs.nms)
num.obs <- data.frame(rsID=rs.nms, nobs=obs)
dim(num.obs) # 104 2
summary(num.obs$nobs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7    1070    1277    1867    2936    2953 

acmg$opensnp.numobs <- num.obs$nobs[match(acmg$Name, num.obs$rsID)]
acmg %>% print(width=Inf)

# look at the SNPs with MCR > 5%. how many observations do they have?
acmg %>% filter(opensnp.miss.frac > 0.05) %>% 
  select(Name, acmg.gene, maf_All, opensnp.miss.frac, opensnp.numobs)
# # A tibble: 3 x 5
# Name       acmg.gene maf_All opensnp.miss.frac opensnp.numobs
# <chr>      <chr>       <dbl>             <dbl>          <dbl>
# 1 rs11583680 PCSK9      0.0804            0.354            2747
# 2 rs2959656  MEN1       0.216             0.0610           1033
# 3 rs3026760  RET        0                 0.0663            996

# no, all have high # of observations.
