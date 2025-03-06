version 1.0

workflow imputation_server_filter {
     input {
        Array[Array[File]] info_files
        Float r2_minimum
     }

    scatter(files in info_files) {

        scatter (info_file in files) {
            call filter_info_file {
                input: info_file = info_file,
                        r2_minimum = r2_minimum
            }

        }

        call intersect {
            input: vcf_files = filter_info_file.filtered_info_file,
                    index_files = filter_info_file.filtered_index_file
        }

    }

    call concat_files {
        input: files = intersect.filtered_id_file
    }

    output {
        File variant_file = concat_files.output_file
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
            -G \
            -e "INFO/R2 < ~{r2_minimum}" \
            ~{info_file} \
            -o subset.vcf.gz

        bcftools index subset.vcf.gz
    }

    output {
        File filtered_info_file = "subset.vcf.gz"
        File filtered_index_file = "subset.vcf.gz.csi"
    }

    runtime {
        docker: "staphb/bcftools:1.20"
        disks: "local-disk ${disk_gb} SSD"
    }
}


task intersect {
    input {
        Array[File] vcf_files
        Array[File] index_files
    }

    Int n_files = length(vcf_files)

    command <<<
        set -e -o pipefail

        # We need to create a file that lists the VCF files and their indexes.
        # Index file paths cannot be assumed due to localization, so we will explicitly specify them.
        # Format expected by bcftools:
        # <VCF_FILE>##idx##<INDEX_FILE>

        echo "writing input file"
        VCF_ARRAY=(~{sep=" " vcf_files}) # Load array into bash variable
        INDEX_ARRAY=(~{sep=" " index_files}) # Load array into bash variable
        for idx in ${!VCF_ARRAY[*]}
        do
            echo "${VCF_ARRAY[$idx]}##idx##${INDEX_ARRAY[$idx]}"
        done > files.txt

        # Perform the intersection
        bcftools isec \
            -n =2 \
            -l files.txt \
            -w 1 |
            bcftools query \
            -f '%CHROM:%POS:%REF:%ALT\n' \
            > filtered_ids.txt
    >>>

    output {
        File filtered_id_file = "filtered_ids.txt"
    }

    runtime {
        docker: "staphb/bcftools:1.20"
    }

}

task concat_files {
    input {
        Array[File] files
    }

    Int disk_gb = ceil(size(files, "GB")*2.5) + 5

    command <<<
        cat ~{sep=' ' files} > concat.txt
    >>>

    output {
        File output_file = "concat.txt"
    }

    runtime {
        docker: "staphb/bcftools:1.20"
        disks: "local-disk ${disk_gb} SSD"
    }
}
