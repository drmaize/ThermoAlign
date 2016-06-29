ThermoAlign:
================================================
ThermoAlign is a a genome-aware oligonucleotide design algorithm embedded within a distributable tool (i.e. consolidating analysis of thermodynamics and sequence alignment). ThermoAlign determines the priming specificity of oligonucleotides based on genome-wide analysis of the thermodynamics of hybridization for full-length and relevant primer-template interactions in the presence of mismatches, to design and select template specific oligonucleotides for any species with a reference genome sequence. It is also capable of using prior information on genomic variants during the selection of candidate oligonucleotides. In addition to the design of primer pairs for standard, single locus PCR amplification, ThermoAlign uses a directed graph analysis approach to design multiplex sets of tiled primer pairs for targeted sequencing of large segments of a genome.


Running ThermoAlign:
================================================
ThermoAlign is released under a GNU GPLv3 open source license. The [source code](https://github.com/drmaize/ThermoAlign/tree/master/TA_codes) can be found [here](https://github.com/drmaize/ThermoAlign/tree/master/TA_codes)


Requirements for directly running the source code
================================================
Linux/Unix
Python  2.7
NumPy   1.9.2
pandas  0.18.1
Cython  0.24
primer3-py  0.5.1
BLAST   2.2.31+
MultiPLX    2.0
networkx    1.11


ThermoAlign Docker images:
================================================
Beyond the source code, the following Docker images are made available: (i) TA_1.0.0_d is a general distributable version which requires user supplied files; (ii) TA_1.0.0_s is a sample run version containing a small set of sample files that can be used to test ThermoAlign; (iii) TA_1.0.0_Zm3 is a maize ready version containing all components required for running ThermoAlign as described in this study.

These docker images help automate deployment of ThermoAlign inside docker containers. It is an efficient way to port ThermoAlign across systems and operating systems.

You may also directly use the [ThermoAlign docker reference file](https://github.com/drmaize/ThermoAlign/blob/master/Dockerfile) to create your own ThermoAlign docker images.


Installing Docker:
================================================
Docker installation can be done by visiting the [official Docker installation page](https://docs.docker.com/engine/installation/) and following the instructions tailored for your operating system.


