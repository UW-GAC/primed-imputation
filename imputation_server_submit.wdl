version 1.0

workflow imputation_server_submit {
     input {
          String hostname
          String token
          Array[File] vcf_files
          Boolean multi_chrom_file
          String build
          String refpanel
          String population
          String password
     }

     if (multi_chrom_file) {
          call split_by_chrom {
               input: vcf_file = vcf_files[0]
          }
     }
     Array[File] inp_files = select_first([split_by_chrom.chrom_files, vcf_files])

     call submit { 
          input: hostname = hostname,
                 token = token,
                 vcf_files = inp_files,
                 build = build,
                 refpanel = refpanel,
                 population = population,
                 password = password
     }

     output {
          String job_id = submit.job_id
     }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
     }
}

task split_by_chrom {
     input {
          File vcf_file
     }

     String vcf_basename = basename(vcf_file, ".vcf.gz")

     command {
          bcftools index ${vcf_file}
          bcftools query -f '%CHROM\n' ${vcf_file} | sort -u > chroms.txt
          while read -r c; do
               bcftools view --regions "$c" -Oz -o ${vcf_basename}".$c.vcf.gz" ${vcf_file}
          done < chroms.txt
     }

     output {
          Array[File] chrom_files = glob("*.vcf.gz")
     }

     runtime {
        docker: "staphb/bcftools:1.16"
    }
}

task submit {
     input {
          String hostname
          String token
          Array[File] vcf_files
          String build
          String refpanel
          String population
          String password
          Boolean meta_imputation = true
     }

     String server = hostname + "/api/v2/jobs/submit/" + if (hostname == "https://imputation.biodatacatalyst.nhlbi.nih.gov") then "imputationserver" else "minimac4"

     command {
          #mkdir ~/.imputationbot
          #printf -- "-  hostname: %s\n   token: %s\n" ${hostname} ${token} > ~/.imputationbot/imputationbot.instances
          #imputationbot impute --file ${sep=' ' vcf_files} --build ${build} --refpanel ${refpanel} --population ${population} --password ${password} > tmp
          #grep -o \'job.*\' tmp | sed "s/'//g" > job_id.txt

          curl ${server} \
               -X "POST" \
               -H "X-Auth-Token: ${token}" \
               -F "files=@${sep='" -F "files=@' vcf_files}" \
               -F "build=${build}" \
               -F "refpanel=apps@${refpanel}" \
               -F "population=${population}" \
               -F "meta=${true='yes' false ='no' meta_imputation}" \
               -F "password=${password}" \
               > job_id.json
          cat job_id.json | json id > job_id.txt
     }

     output {
          String job_id = read_string("job_id.txt")
     }

     runtime {
          docker: "uwgac/primed-imputation:0.1.0"
     }
}
