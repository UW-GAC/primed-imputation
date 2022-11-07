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

     command {
          bcftools query -f '%CHROM\n' ${vcf_file} | sort -u > chroms.txt
          while read -r c; do
               bcftools view --regions "$c" -Oz -o "chr$c.vcf.gz"
          done < chroms.txt
     }

     output {
          Array[File] chrom_files = glob("*.vcf.gz")
     }

     runtime {
        docker: "xbrianh/xsamtools:v0.5.2"
    }
}
