FROM ubuntu:latest
MAINTAINER Felix Francis <felixfrancier@gmail.com>

# Install all the softwares & dependencies required to run the pipeline
RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y python-pip
#RUN apt-get install python-dev build-essential 
RUN pip install --upgrade pip
RUN pip install shutilwhich
RUN pip install biopython
RUN pip install numpy
RUN pip install pandas
RUN pip install Cython
RUN pip install primer3-py
RUN pip install networkx