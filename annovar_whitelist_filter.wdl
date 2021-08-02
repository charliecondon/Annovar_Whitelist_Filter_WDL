version 1.0

## Version 8-2-2021
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
## run_whitelist: if true, the WhitelistFilter task is run
## sample_id: set to the corresponding sample id for a given run
##			  NOTE: This is set on Terra (ex. this.sample_id)
## whitelist_filter_zip: the zipped folder with all of the files needed to run WhitelistFilter
##                       NOTE: This file path is set on Terra - The file must be in the Workspace's bucket
## txt_input: the txt input file that was an output of annovar
## whitelist_filter_docker: the docker image to be used in the Annovar task
##
## ** WDL OUTPUTS **
##  Four CSV files:
## 		- one with the whitelist filter applied
##		- one ready for manual review
##		- one with variant count information for debugging
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
      Boolean? run_whitelist

      String annovar_docker = "perl:5.34.0"
      String whitelist_filter_docker = "ccondon/whitelist_filter:latest"
    }

    Boolean run_whitelist_or_default = select_first([run_whitelist, true])

    call Annovar {
      input:
        annovar_docker = annovar_docker,
        sample_id = sample_id,
        vcf_input = annovar_vcf_input
    }

    if (run_whitelist_or_default) {
      File whitelist_filter_annovar_txt_input = annovar_annotated_file_table
      call WhitelistFilter {
        input:
          whitelist_filter_docker = whitelist_filter_docker,
          sample_id = sample_id,
          txt_input = whitelist_filter_annovar_txt_input
      }
    }

    output {
      File annovar_annotated_file_vcf = Annovar.annovar_output_file_vcf
      File annovar_annotated_file_table = Annovar.annovar_output_file_table # this is the .txt file that we need for the R script

      File? whitelist_filter_output_wl = WhitelistFilter.whitelist_filter_output_wl_csv
      File? whitelist_filter_output_manual_review = WhitelistFilter.whitelist_filter_output_manual_review_csv
      File? whitelist_filter_output_varcount = WhitelistFilter.whitelist_filter_output_varcount_csv
      File? whitelist_filter_output_allvariants = WhitelistFilter.whitelist_filter_output_allvariants_csv
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

      String sample_id
      String file_prefix = sample_id + ".annovar_out"
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
        -out ${file_prefix} \
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
      File annovar_output_file_vcf = file_prefix + ".hg38_multianno.vcf"
      File annovar_output_file_table = file_prefix + ".hg38_multianno.txt"
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

      unzip ~{whitelist_filter_zip}
      mv whitelist_filter_files/* .

      wget https://raw.githubusercontent.com/charliecondon/Annovar_Whitelist_Filter_WDL/main/whitelist_filter_rscript.R

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
      File? whitelist_filter_output_varcount_csv = sample_id + ".annovar.varsOI.varcount.csv"
      File? whitelist_filter_output_allvariants_csv = sample_id + ".annovar.varsOI.allvariants.csv"
      File? whitelist_filter_output_wl_csv = sample_id + ".annovar.varsOI.wl.csv"
      File? whitelist_filter_output_manual_review_csv = sample_id + ".annovar.varsOI.manualreview.csv"
    }
}
