version 1.0

## Version 07-2-2021
##
## This WDL workflow runs Annovar and a Whitelist Filter on the ouput VCFs from the Mutect2 Workflow.
##
##
## ** ANNOVAR **
## Annovar functionally annotates genetic variants detected from diverse genomes.
## Given a list of variants with chromosome, start position, end position, reference nucleotide
## and observed nucleotides, Annovar can perform gene-based annotation, region-based annotation,
## filter-based annotation, and more.
##
## See ANNOVAR documentation to fully understand functionality:
## https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
##
## annovar_zip: the zipped folder with all of the needed files to run Annovar
##              NOTE: This file path is set on Terra - The file must be in the Workspace's bucket
## annovar_vcf_input: the Tables/sample column containing the vcf output files from a run of Mutect2
##                    NOTE: This is set on Terra (ex. this.filtered_vcf)
## annovar_protocols: the specificed protocols needed to run annovar (default = refGene,cosmic70)
##                    NOTE: You must add the needed file paths to annovar_data_sources
## annovar_operation: the specified operations needed to run annovar (default = g,f)
##                    NOTE: They must match up with annovar_protocols
## ref_name: the reference name needed for annovar to run (default = hg38)
## annovar_docker: the docker image to be used in the Annovar task
##
## ** WHITELIST_FILTER **
## WhitelistFilter filters annovar's output based on only relevant data to our lab's whitelist.
##
## You can find the R script code with comments on github: https://github.com/charliecondon/Annovar_Whitelist_Filter_WDL
##
## sample_id: set to the corresponding sample id for a given run
##			  NOTE: This is set on Terra (ex. this.sample_id)
## whitelist_filter_zip: the zipped folder with all of the needed files, including the .R Script file, needed to run WhitelistFilter
##                       NOTE: This file path is set on Terra - The file must be in the Workspace's bucket
## txt_input: the txt input file that was an output of annovar
## whitelist_filter_docker: the docker image to be used in the Annovar task
##
## ** WDL OUTPUTS **
##  Four CSV files: 
## 		- one with the whitelist filter applied 
##		- one ready for manual review
##		- one with check information for debugging 
##		- one with all the pre-whitelist variants listed
##
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Charlie Condon
## Contact <ccondon@vols.utk.edu>

workflow AnnovarAndWhitelistFilter {
    input {
      File annovar_vcf_input
      String sample_id

      String annovar_docker = "perl:5.34.0"
      String whitelist_filter_docker = "ccondon/whitelist_filter:latest"
    }

    call Annovar {
      input:
        annovar_docker = annovar_docker,
        vcf_input = annovar_vcf_input
    }

    File whitelist_filter_annovar_txt_input = annovar_annotated_file_table
    call WhitelistFilter {
      input:
        whitelist_filter_docker = whitelist_filter_docker,
        sample_id = sample_id,
        txt_input = whitelist_filter_annovar_txt_input
    }

    output {
      File annovar_annotated_file_vcf = Annovar.annovar_output_file_vcf
      File annovar_annotated_file_table = Annovar.annovar_output_file_table # this is the .txt file that we need for the R script

      File whitelist_filter_output_wl = WhitelistFilter.whitelist_filter_output_wl_csv
      File whitelist_filter_output_manual_review = WhitelistFilter.whitelist_filter_output_manual_review_csv
      File whitelist_filter_output_check = WhitelistFilter.whitelist_filter_output_check_csv
      File whitelist_filter_output_allvariants = WhitelistFilter.whitelist_filter_output_allvariants_csv
    }
}

task Annovar {
    input {
      Float memory = 4.0
      Int annovar_disk_space = 300
      Int cpu = 1
      Int preemptible = 1
      Int maxRetries = 0
      String annovar_docker

      File vcf_input
      File annovar_zip

      String ref_name = "hg38"
      String annovar_protocols = "refGene,cosmic70"
      String annovar_operation = "g,f"
    }

    command {
      set -euo pipefail

      cp ~{annovar_zip} .
      unzip annovar_files.zip

      chmod +x annovar_files/convert2annovar.pl
      chmod +x annovar_files/table_annovar.pl
      chmod +x annovar_files/annotate_variation.pl
      chmod +x annovar_files/coding_change.pl
      chmod +x annovar_files/retrieve_seq_from_fasta.pl
      chmod +x annovar_files/variants_reduction.pl

      perl annovar_files/table_annovar.pl ${vcf_input} annovar_files \
        -buildver ${default="hg38" ref_name} \
        -out "annovar_out" \
        -remove \
        -protocol ${default="refGene,cosmic70" annovar_protocols} \
        -operation ${default="g,f" annovar_operation} \
        -nastring . -vcfinput
    }

    runtime {
      docker: annovar_docker
      memory: memory + " GiB"
      disk: "local-disk " + annovar_disk_space + " HDD"
      cpu: cpu
      preemptible: preemptible
      maxRetries: maxRetries
    }

    output {
      File annovar_output_file_vcf = "annovar_out.hg38_multianno.vcf"
      File annovar_output_file_table = "annovar_out.hg38_multianno.txt"
    }
}

task WhitelistFilter {
    input {
      Float memory = 10.0
      Int whitelist_filter_disk_space = 300
      Int cpu = 1
      Int preemptible = 1
      Int maxRetries = 0
      String whitelist_filter_docker

      File txt_input
      String sample_id
      File whitelist_filter_zip
    }

    command {
      set -euo pipefail
      
      cp ~{whitelist_filter_zip} .
      unzip whitelist_filter_files.zip
      mv whitelist_filter_files/* .
      
      cp ~{txt_input} .

      Rscript whitelist_filter_rscript.R --args ~{sample_id}
    }

    runtime {
      docker: whitelist_filter_docker
      memory: memory + " GiB"
      disk: "local-disk " + whitelist_filter_disk_space + " HDD"
      cpu: cpu
      preemptible: preemptible
      maxRetries: maxRetries
    }

    output {
      File whitelist_filter_output_check_csv = sample_id + ".annovar.varsOI.check.csv"
      File whitelist_filter_output_allvariants_csv = sample_id + ".annovar.varsOI.allvariants.csv"
      File whitelist_filter_output_wl_csv = sample_id + ".annovar.varsOI.wl.csv"
      File whitelist_filter_output_manual_review_csv = sample_id + ".annovar.varsOI.manualreview.csv"
    }
}
