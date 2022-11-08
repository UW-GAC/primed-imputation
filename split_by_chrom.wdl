version 1.0

workflow chrom_split {
     input {
          Array[File] vcf_files
     }

     call split_by_chrom {
          input: vcf_file = vcf_files[0]
     }

     output {
          Array[File] chrom_files = split_by_chrom.chrom_files
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
          File chrom_list = "chroms.txt"
     }

     runtime {
        docker: "staphb/bcftools:1.16"
    }
}
