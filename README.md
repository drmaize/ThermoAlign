ThermoAlign:
================================================
ThermoAlign is a a genome-aware oligonucleotide design algorithm embedded within a distributable tool (i.e. consolidating analysis of thermodynamics and sequence alignment). ThermoAlign determines the priming specificity of oligonucleotides based on genome-wide analysis of the thermodynamics of hybridization for full-length and relevant primer-template interactions in the presence of mismatches, to design and select template specific oligonucleotides for any species with a reference genome sequence. It is also capable of using prior information on genomic variants during the selection of candidate oligonucleotides. In addition to the design of primer pairs for standard, single locus PCR amplification, ThermoAlign uses a directed graph analysis approach to design multiplex sets of tiled primer pairs for targeted re-sequencing applications.


Running ThermoAlign:
================================================
ThermoAlign is released under a GNU GPLv3 open source license. The [source code](https://github.com/drmaize/ThermoAlign/tree/master/TA_codes) can be found [here.](https://github.com/drmaize/ThermoAlign/tree/master/TA_codes)

The required parameters need to be used in the [parameters.py](https://github.com/drmaize/ThermoAlign/blob/master/TA_codes/parameters.py) file.

The shell script [pipeline.sh](https://github.com/drmaize/ThermoAlign/blob/master/TA_codes/pipeline.sh) would run all the required ThermoAlign components.

Alternatively, individual components of ThermoAlign may also be separately run in the following order:
1) TRS.py; 2) UOD.py 3) PSE.py 4) PPS.py



Requirements for directly running the source code
================================================
* Linux/Unix
* [Python  2.7](http://python.org/)
* [NumPy   1.9.2](http://www.numpy.org/)
* [pandas  0.18.1](http://pandas.pydata.org/)
* [Cython  0.24](http://cython.org/)
* [primer3-py  0.5.1](https://pypi.python.org/pypi/primer3-py)
* [BLAST   2.2.31+](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [MultiPLX    2.0](http://bioinfo.ut.ee/download/dl.php?file=24)
* [networkx    1.11](https://networkx.github.io/)


ThermoAlign Docker images:
================================================
Beyond the source code, the following Docker images are made available: 
* (i) TA_1.0.0_s is a sample run version containing a small set of sample files that can be used to test ThermoAlign
* (ii) TA_1.0.0_d is a general distributable version which requires user supplied files
* (iii) TA_1.0.0_Zm3 is a maize ready version containing all components required for running ThermoAlign as described in this study

    
These docker images help automate deployment of ThermoAlign inside docker containers. It is an efficient way to port ThermoAlign across systems and operating systems.


Installing Docker:
================================================
Docker installation can be done by visiting the [official Docker installation page](https://docs.docker.com/engine/installation/) and following the instructions tailored for your operating system.

Basic [docker commands](https://goo.gl/TfU9AY) to run software containers such as ThermoAlign can be found [here](https://goo.gl/TfU9AY) or on [docs.docker.com](https://docs.docker.com/).


To run a docker image as a container that can be accessed via bash:

    docker run -t -i drmaize/thermoalign:TA_1.0.0_s /bin/bash
                            OR
    docker run -t -i drmaize/thermoalign:TA_1.0.0_d /bin/bash
                            OR
    docker run -t -i drmaize/thermoalign:TA_1.0.0_Zm3 /bin/bash

This will generate a new ThermoAlign container everytime its used.


Build your own ThermoAlign docker image:
================================================

You may directly use the [ThermoAlign docker reference file](https://github.com/drmaize/ThermoAlign/blob/master/Dockerfile) to create your own ThermoAlign docker images.
    
    docker build https://github.com/drmaize/ThermoAlign.git

    
Support:
================================================
For support using ThermoAlign, please visit the ThermoAlign Google Group.







   
                                    ##### END OF README #####    

