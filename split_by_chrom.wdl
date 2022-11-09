version 1.0

import "imputation_server_submit.wdl" as tasks

workflow chrom_split {
     input {
          File vcf_file
     }

     call tasks.split_by_chrom {
          input: vcf_file = vcf_file
     }

     output {
          Array[File] chrom_files = split_by_chrom.chrom_files
     }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
     }
}
