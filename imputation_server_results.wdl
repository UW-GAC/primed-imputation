version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-checks/main/validate_genotype_model.wdl" as validate

workflow imputation_server_results {
     input {
          String hostname
          String token
          String job_id
          String password
          Int? disk_gb
          String sample_set_id
          String source_dataset_id
          String source_genotypes
          String refpanel
          Float r2_filter
          String model_url
          String workspace_name
          String workspace_namespace
          Boolean import_tables = true
          Boolean overwrite = true
     }

     String reference_assembly = if (refpanel == "topmed-r2" || refpanel == "1000g-phase3-low" || refpanel == "1000g-phase3-deep") then "GRCh38" else "GRCh37"
     String quality_filter = if (r2_filter == 0) then "None" else "r2=" + r2_filter

     call results { 
          input: hostname = hostname,
                 token = token,
                 job_id = job_id,
                 password = password,
                 disk_gb = disk_gb
     }

    call imputation_data_model {
          input: imputed_files = results.imputed,
               sample_set_id = sample_set_id,
               source_dataset_id = source_dataset_id,
               source_genotypes = source_genotypes,
               reference_panel = refpanel,
               reference_assembly = reference_assembly,
               quality_filter = quality_filter
    }

    call validate.validate {
          input: table_files = imputation_data_model.table_files,
               model_url = model_url,
               workspace_name = workspace_name,
               workspace_namespace = workspace_namespace,
               overwrite = overwrite,
               import_tables = import_tables
    }
    
     output {
          Array[File] imputed = results.imputed
          File md5 = results.md5
          File qc_report = results.qc_report
          Array[File] qc_stats = results.qc_stats
          Array[File] log = results.log
          File validation_report = validate.validation_report
          Array[File]? tables = validate.tables
     }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
     }
}

task results {
     input {
          String hostname
          String token
          String job_id
          String password
          Int disk_gb = 10
     }

     command {
          mkdir ~/.imputationbot
          printf -- "-  hostname: %s\n   token: %s\n" ${hostname} ${token} > ~/.imputationbot/imputationbot.instances
          imputationbot download ${job_id} --password ${password}
     }

     output {      
          #Array[File] imputed = glob("${job_id}/*/*")
          Array[File] imputed = glob("${job_id}/local/*.gz")
          File md5 = "${job_id}/local/results.md5"
          File qc_report = "${job_id}/qcreport/qcreport.html"
          Array[File] qc_stats = glob("${job_id}/statisticDir/*.txt")
          Array[File] log = glob("${job_id}/logfile/*.log")
     }

     runtime {
          docker: "uwgac/primed-imputation:0.2.0"
          disks: "local-disk ${disk_gb} SSD"
     }
}


task imputation_data_model {
     input {
          Array[String] imputed_files
          String sample_set_id
          String source_dataset_id
          String source_genotypes
          String reference_panel
          String reference_assembly
          String imputation_software = "Minimac4"
          String phasing_software = "Eagle2.4"
          String quality_filter
     }

     command <<<
        Rscript -e "\
        library(dplyr); \
        library(stringr); \
        dat <- tibble(field='sample_set_id', value='~{sample_set_id}'); \
        dat <- bind_rows(dat, tibble(field='source_dataset_id', value='~{source_dataset_id}')); \
        dat <- bind_rows(dat, tibble(field='source_genotypes', value='~{source_genotypes}')); \
        dat <- bind_rows(dat, tibble(field='reference_panel', value='~{reference_panel}')); \
        dat <- bind_rows(dat, tibble(field='reference_assembly', value='~{reference_assembly}')); \
        dat <- bind_rows(dat, tibble(field='imputation_software', value='~{imputation_software}')); \
        dat <- bind_rows(dat, tibble(field='phasing_software', value='~{phasing_software}')); \
        dat <- bind_rows(dat, tibble(field='quality_filter', value='~{quality_filter}')); \
        readr::write_tsv(dat, 'imputation_dataset_table.tsv'); \
        parse_array <- function(x) unlist(strsplit(x, split=' ', fixed=TRUE)); \
        files <- parse_array('~{sep=' ' imputed_files}'); \
        chr <- str_extract(files, 'chr[:alnum:]+[:punct:]'); \
        chr <- sub('chr', '', chr, fixed=TRUE); \
        chr <- sub('.', '', chr, fixed=TRUE); \
        file_type <- ifelse(grepl('vcf', files), 'VCF', ifelse(grepl('info', files), 'quality metrics', 'supporting file')); \
        dat <- tibble(file_path = files, chromosome = chr, file_type = file_type); \
        writeLines(dat[['file_path']], 'files.txt'); \
        readr::write_tsv(dat, 'imputation_file_table.tsv'); \
        "
        while read f; do
            echo $f
            gsutil ls -L $f | grep "md5" | awk '{print $3}' > md5_b64.txt
            echo "b64 checksum: "; cat md5_b64.txt
            python3 -c "import base64; import binascii; print(binascii.hexlify(base64.urlsafe_b64decode(open('md5_b64.txt').read())))" | cut -d "'" -f 2 >> md5_hex.txt
            echo "hex checksum: "; cat md5_hex.txt
        done < files.txt
        Rscript -e "\
        dat <- readr::read_tsv('imputation_file_table.tsv'); \
        md5_hex <- readLines('md5_hex.txt'); \
        dat <- dplyr::mutate(dat, md5sum=md5_hex); \
        readr::write_tsv(dat, 'imputation_file_table.tsv'); \
        "
     >>>

     output {
        Map[String, File] table_files = {
            "imputation_dataset": "imputation_dataset_table.tsv",
            "imputation_file": "imputation_file_table.tsv"
        }
     }

     runtime {
          docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.18.0"
     }
}
