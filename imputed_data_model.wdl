version 1.0

import "imputation_server_results.wdl" as tasks

workflow run_data_model {
     input {
          Array[String] imputed_files
          String sample_set_id
          String source_dataset_id
          String source_genotypes
          String reference_panel
          String reference_assembly
          String quality_filter
     }

    call tasks.imputation_data_model {
        input: imputed_files = imputed_files,
        sample_set_id = sample_set_id,
        source_dataset_id = source_dataset_id,
        source_genotypes = source_genotypes,
        reference_panel = reference_panel,
        reference_assembly = reference_assembly,
        quality_filter = quality_filter
    }
    
    output {
        Map[String, File] table_files = imputation_data_model.table_files
    }

    meta {
        author: "Stephanie Gogarten"
        email: "sdmorris@uw.edu"

    }
}
