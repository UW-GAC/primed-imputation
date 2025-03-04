version 1.0

workflow imputation_server_filter {
     input {
        Array[File] info_files
        Float r2_minimum
     }

    scatter (info_file in info_files) {
        call filter_info_file {
            input: info_file = info_file,
                    r2_minimum = r2_minimum
        }
    }

     output {
          File variant_id_file = intersect_by_chromosome.filtered_id_file
     }

    call intersect_by_chromosome {
        input: filtered_files = filter_info_file.filtered_info_file
    }

     meta {
          author: "Adrienne Stilp"
          email: "amstilp@uw.edu"
     }
}


task filter_info_file {
    input {
        File info_file
        Float r2_minimum
    }

    Int disk_gb = ceil(size(info_file, "GB")) + 2
    String output_file = "subset_" + basename(info_file)

    command {
        bcftools query \
            -f '%CHROM:%POS:%REF:%ALT\n' \
            -e "INFO/R2 < ~{r2_minimum}" \
            ~{info_file} \
            -o ~{output_file}
    }

    output {
        File filtered_info_file = output_file
    }

    runtime {
        docker: "staphb/bcftools:1.16"
        disks: "local-disk ${disk_gb} SSD"
    }
}


task intersect_by_chromosome {
    input {
        Array[File] filtered_files
    }

    command {
        Rscript /usr/local/primed-imputation/intersect_by_chromosome.R \
            --input ${sep=' ' filtered_files} \
            --output filtered_ids.txt
    }

    output {
        File filtered_id_file = "filtered_ids.txt"
    }

    runtime {
        docker: "uwgac/primed-imputation-filter:0.1.0"
    }
}
