version 1.0

workflow imputation_server_filter {
     input {
        File info_file
        Float r2_minimum
     }

     call filter_info_file {
          input: info_file = info_file,
                 r2_minimum = r2_minimum
     }

     output {
          File filtered_info_file = filter_info_file.filtered_info_file
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

    command {
        bcftools view \
            -e "INFO/R2 < ~{r2_minimum}" \
            ~{info_file} \
            -o filtered.txt.gz
    }

    output {
        File filtered_info_file = "filtered.txt.gz"
    }

    runtime {
        docker: "staphb/bcftools:1.16"
        disks: "local-disk ${disk_gb} SSD"
    }
}
