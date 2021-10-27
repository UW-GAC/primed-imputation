version 1.0

workflow bdc_imputation_submit {
     input {
          String token
     }

     call submit { 
          input: token = token 
     }

     output {
          File output_file = submit.output_file
     }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
     }
}
       
task submit {
     input {
          String token
     }

     command {
          mkdir ~/.imputationbot
          printf -- "-  hostname: https://imputation.biodatacatalyst.nhlbi.nih.gov\n   token: %s\n" ${token} > ~/.imputationbot/imputationbot.instances
          imputationbot instances
          imputationbot refpanels
     }

     output {      
          File output_file = stdout()
     }

     runtime {
          docker: "uwgac/primed-imputation:0.1.0"
     }
}
