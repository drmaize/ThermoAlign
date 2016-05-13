
FROM ubuntu:latest
MAINTAINER Felix Francis <felixfrancier@gmail.com>

# Install all the softwares & dependencies required to run the pipeline
RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y python-pip
#RUN apt-get install python-dev build-essential 
RUN pip install --upgrade pip shutilwhich biopython numpy pandas
RUN pip install Cython
RUN pip install primer3-py
RUN pip install networkx

# Download and extract multiplx software
RUN wget --content-disposition http://bioinfo.ut.ee/download/dl.php?file=24
RUN tar xzf multiplx_linux_64_20101011.tar.gz
RUN rm -rf multiplx_linux_64_20101011.tar.gz

# Install BLASTn 
RUN apt-get install -y ncbi-blast+
RUN PATH=$PATH:~/opt/bin