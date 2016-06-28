FROM ubuntu:latest

MAINTAINER Felix Francis <felixfrancier@gmail.com> 

# Install all the softwares & dependencies required to run the pipeline
RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y git
RUN apt-get install -y python-pip                                                       # v 2.7.11+
#RUN apt-get install python-dev build-essential 
RUN pip install --upgrade pip shutilwhich numpy pandas biopython
RUN pip install Cython                                                                  # v 0.24
RUN pip install primer3-py                                                              # v 0.5.1
RUN pip install networkx                                                                # v 1.11

# Download and extract multiplx software
RUN wget --content-disposition http://bioinfo.ut.ee/download/dl.php?file=24
RUN tar xzf multiplx_linux_64_20101011.tar.gz
RUN rm -rf multiplx_linux_64_20101011.tar.gz

# Install BLASTn 
RUN apt-get install -y ncbi-blast+
RUN PATH=$PATH:~/opt/bin

# Download TA codes and sample genome, polymorphism vcf files
RUN git clone https://github.com/drmaize/ThermoAlign.git
RUN mv -t ./ ThermoAlign/sample_genome/ ThermoAlign/sample_vcf/ ThermoAlign/TA_codes/
RUN rm -rf ThermoAlign/
