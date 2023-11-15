version 1.0

workflow imputation_server_results {
     input {
          String hostname
          String token
          String job_id
          String password
          Int? disk_gb
     }

     call results { 
          input: hostname = hostname,
                 token = token,
                 job_id = job_id,
                 password = password,
                 disk_gb = disk_gb
     }

     output {
          Array[File] imputed = results.imputed
          File md5 = results.md5
          File qc_report = results.qc_report
          Array[File] qc_stats = results.qc_stats
          Array[File] log = results.log
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
          docker: "uwgac/primed-imputation:0.1.0"
          disks: "local-disk ${disk_gb} SSD"
     }
}
