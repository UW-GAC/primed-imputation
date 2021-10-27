version 1.0

workflow imputation_server_submit {
     input {
          String hostname
          String token
          Array[File] vcf_files
          String refpanel
          String population
          String password
     }

     call submit { 
          input: hostname = hostname,
                 token = token,
                 vcf_files = vcf_files,
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

task submit {
     input {
          String hostname
          String token
          Array[File] vcf_files
          String refpanel
          String population
          String password
     }

     command {
          mkdir ~/.imputationbot
          printf -- "-  hostname: %s\n   token: %s\n" ${hostname} ${token} > ~/.imputationbot/imputationbot.instances
          imputationbot impute --file ${vcf_files} --refpanel ${refpanel} --population ${population} --password ${password} > tmp
          grep -o \'job.*\' tmp | sed "s/'//g" > job_id.txt
     }

     output {
          String job_id = read_string("job_id.txt")
     }

     runtime {
          docker: "uwgac/primed-imputation:0.1.0"
     }
}
