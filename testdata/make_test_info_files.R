library(tidyverse)
library(glue)

set.seed(123)

chromosomes = 1:3
n_info_files_per_chromosome = 2
n_variants_per_chromosome = 100

header = c(
    '##fileformat=VCFv4.2',
    '##filedate=20250222',
    '##source=Minimac v4.1.6',
    '##phasing=full',
    '##contig=<ID=chr{chromosome}>',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Alternate Allele Frequency">',
    '##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">',
    '##INFO=<ID=AVG_CS,Number=1,Type=Float,Description="Average Call Score">',
    '##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">',
    '##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">',
    '##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed">',
    '##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped">'
)

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO'

for (chromosome in chromosomes) {

    # Generate a list of variants.
    variants = tibble(
        id = glue("chr{chromosome}_{i}", i=1:n_variants_per_chromosome),
        chr=chromosome,
        position = 1:n_variants_per_chromosome,
        ref = sample(c("A", "C", "G", "T"), n_variants_per_chromosome, replace = TRUE),
        alt = sample(c("A", "C", "G", "T"), n_variants_per_chromosome, replace = TRUE),
        type = sample(c("IMPUTED", "TYPED"), prob=c(0.9, 0.1), n_variants_per_chromosome, replace = TRUE),
    )

    for (j in 1:n_info_files_per_chromosome) {
        info = variants %>%
            # Add random data.
            mutate(
                QUAL=".",
                FILTER=".",
                AF=runif(n_variants_per_chromosome, 0, 1),
                MAF=abs(AF - sample(c(0, 1), n_variants_per_chromosome, prob=c(0.9, 0.1), replace=TRUE)),
                AVG_CS=runif(n_variants_per_chromosome, 0, 1),
                R2=runif(n_variants_per_chromosome, 0, 1),
                ER2=runif(n_variants_per_chromosome, 0, 1),
            )

        # Final output format
        info =
            info %>%
            # Formatting expected by VCF
            mutate(across(c(AF, MAF, AVG_CS, R2, ER2), ~sprintf('%0.3f', .))) %>%
            mutate(
                chr = glue("chr{chr}"),
                INFO=glue("AF={AF};MAF={MAF};AVG_CS={AVG_CS};R2={R2};ER2={ER2};{type}")
            ) %>%
            select(
                `#CHROM`=chr,
                POS = position,
                ID = id,
                REF = ref,
                ALT = alt,
                QUAL,
                FILTER,
                INFO
            )

        # Save it with the header.
        outfile = "info_chr{chromosome}_{j}.txt.gz" %>% glue()
        header %>% sapply(glue) %>% write_lines(outfile)
        info %>% write_tsv(outfile, append = TRUE, col_names = TRUE)
    }
}
