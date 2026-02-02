library(dplyr)
library(data.table)
library(stringr)
library(glue)
library(Matrix)

pu_dir = './out_pileup'

vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
    rename(CHROM = `#CHROM`) %>%
    mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_')) %>%
    mutate(CHROM = str_remove(CHROM, 'chr')) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
    filter(CHROM != 'X')

vcf_phased = vcf_pu %>% mutate(GT = '1|0', cM = 0)

# pileup count matrices
AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))
barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F)$V1

# convert to dataframe
DP = as.data.frame(Matrix::summary(DP)) %>%
    mutate(
        cell = barcodes[j],
        snp_id = vcf_pu$snp_id[i]
    ) %>%
    select(-i, -j) %>%
    rename(DP = x) %>%
    select(cell, snp_id, DP)

AD = as.data.frame(Matrix::summary(AD)) %>%
    mutate(
        cell = barcodes[j],
        snp_id = vcf_pu$snp_id[i]
    ) %>%
    select(-i, -j) %>%
    rename(AD = x) %>%
    select(cell, snp_id, AD)

df_allele = DP %>% left_join(AD, by = c("cell", "snp_id")) %>%
    mutate(AD = ifelse(is.na(AD), 0, AD))

# attach genotype info
df_allele = df_allele %>% inner_join(
    vcf_phased %>% select(snp_id, CHROM, POS, REF, ALT, GT, cM),
    by = 'snp_id')



# out
saveRDS(df_allele, 'df_allele.rds')


