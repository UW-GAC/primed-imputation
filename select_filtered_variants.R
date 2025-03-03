library(argparser)
library(tidyverse)

# Parse command line arguments
argp <- arg_parser("Filter imputation server results") %>%
    add_argument("--input", help="Input variant files; one line per variant with format chr:pos:ref:alt", nargs=Inf) %>%
    add_argument("--output", help="Output variant file")

argv <- parse_args(argp)
print(argv)

(input_files = argv$input)

# Read input files
input = list()
for (file in input_files) {
    input[[file]] = read_tsv(file, col_names="SNP", col_types=cols(col_character())) %>%
        separate_wider_delim(SNP, delim=":", names=c("chr"), too_many="drop", cols_remove=FALSE)
}

# Make sure there is only one chromosome per file.
stopifnot(all(sapply(input, function(x) length(unique(x$chr))) == 1))

files_by_chromosome = input %>%
    lapply(slice_head, n=1) %>%
    bind_rows(.id="file") %>%
    select(-SNP) %>%
    group_by(chr) %>%
    summarise(files = list(file))
files_by_chromosome %>% unnest(files)

# Now go through the files for each chromosome and find the variants that are in all of them.
variants = list()
for (chromosome in files_by_chromosome$chr) {
    files = files_by_chromosome$files[files_by_chromosome$chr == chromosome] %>% unlist()
    variants[[chromosome]] = input[files] %>%
        bind_rows(.id="file") %>%
        group_by(SNP) %>%
        summarise(n_files = n()) %>%
        filter(n_files == length(files)) %>%
        select(-n_files)
    # # This was the co-pilot suggestion. Maybe it is faster?
    # variants[[chromosome]] = Reduce(intersect, lapply(files, function(file) input[[file]]$chr))
}

variants = variants %>% bind_rows()

# Write out.
write_tsv(variants, argv$output, col_names=FALSE)
