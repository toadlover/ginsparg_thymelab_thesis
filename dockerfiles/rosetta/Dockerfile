FROM ubuntu:bionic-20220315



RUN apt-get update

RUN apt-get install -y python
RUN apt-get install -y python3
#RUN apt-get install -y python-is-python3

RUN apt-get install -y git

RUN apt-get install -y g++

RUN apt-get install -y gcc-6

RUN apt-get install -y build-essential

RUN apt-get install -y zlib1g-dev

RUN git clone https://toadlover:[My personal access token]@github.com/RosettaCommons/rosetta.git

WORKDIR rosetta/source/src

WORKDIR ..

#checkout the branch that has the faster discovery code (will eventually push to master)

RUN git pull

#RUN git checkout ligand_discovery_profiling_november_2023

#RUN git checkout ligand_discovery_method_optimization_jan_2024_ag

RUN git checkout updated_ligand_discovery_search_code_ag_may_2024

RUN echo hi

RUN python --version
RUN python3 --version

#RUN which python3
#RUN alias python=python3
#RUN source ~/.bashrc
#RUN source ~/.bash_aliases
#RUN which python
#RUN ls /usr/bin/env/python3

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.6 3

RUN ./scons.py -j 1 mode=release identify_ligand_motifs extract_motifs ligand_discovery_search_protocol remove_duplicate_motifs

RUN rm -drf tests src external

WORKDIR ..

RUN rm -drf .git

WORKDIR /

RUN apt-get install nano

#RUN /main/source/bin/ligand_discovery_search_protocol.linuxgccrelease
