FROM ubuntu
MAINTAINER Felix Francis <felixfrancier@gmail.com>

# Install all the software needed to run the pipeline
RUN apt-get -qq update
RUN apt-get install -qqy python3-setuptools python3-docutils python3-flask
RUN easy_install3 snakemake
 
 
# clone the most recent version of the pipeline
WORKDIR /home/user/
RUN git clone https://github.com/dalloliogm/pipeline_play.git
 
# download the input files
WORKDIR /home/user/pipeline_play
RUN tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 22:23862589-25055718 > myregion.vcf
 
 
# Execute the pipeline
RUN snakemake