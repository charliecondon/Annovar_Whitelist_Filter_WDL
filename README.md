# Annovar_Whitelist_Filter_WDL

WDL that runs Annovar and a filtering R Script on the output of Mutect2

This WDL workflow can also be found on dockstore: https://dockstore.org/my-tools/registry.hub.docker.com/ccondon/whitelist_filter

The Annovar task uses an official perl image: https://hub.docker.com/_/perl

The WhitelistFilter task uses a custom docker image with R 4.1.0 and needed packages installed: https://hub.docker.com/repository/docker/ccondon/whitelist_filter

# annovar_whitelist_filter.wdl

** ANNOVAR **

Annovar functionally annotates genetic variants detected from diverse genomes.
Given a list of variants with chromosome, start position, end position, reference nucleotide
and observed nucleotides, Annovar can perform gene-based annotation, region-based annotation,
filter-based annotation, and more.

See ANNOVAR documentation to fully understand functionality: https://annovar.openbioinformatics.org/en/latest/user-guide/startup/

annovar_zip: the zipped folder with all of the needed files to run Annovar
  - NOTE: This file path is set on Terra - The file must be in the Workspace's bucket

annovar_vcf_input: the Tables/sample column containing the vcf output files from a run of Mutect2
  - NOTE: This is set on Terra (ex. this.filtered_vcf)

annovar_protocols: the specificed protocols needed to run annovar (default = refGene,cosmic70)
  - NOTE: You must add the needed file paths to annovar_data_sources

annovar_operation: the specified operations needed to run annovar (default = g,f)
  - NOTE: They must match up with annovar_protocols

ref_name: the reference name needed for annovar to run (default = hg38)

annovar_docker: the docker image to be used in the Annovar task



** WHITELIST_FILTER **

WhitelistFilter filters annovar's output based on only relevant data to our lab's whitelist.

You can find the R script code with comments on github: https://github.com/charliecondon/Annovar_Whitelist_Filter_WDL

sample_id: set to the corresponding sample id for a given run
  - NOTE: This is set on Terra (ex. this.sample_id)

whitelist_filter_zip: the zipped folder with all of the needed files, including the .R Script file, needed to run WhitelistFilter
  - NOTE: This file path is set on Terra - The file must be in the Workspace's bucket

txt_input: the txt input file that was an output of annovar

whitelist_filter_docker: the docker image to be used in the Annovar task


** WDL OUTPUTS **

Four CSV files: 
  - one with the whitelist filter applied 
  - one ready for manual review
  - one with check information for debugging 
  - one with all the pre-whitelist variants listed
