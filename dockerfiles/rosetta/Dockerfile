FROM ubuntu:bionic-20220315

#package installs
RUN apt-get update

RUN apt-get install -y python
RUN apt-get install -y python3
RUN apt-get install -y git
RUN apt-get install -y g++
RUN apt-get install -y gcc-6
RUN apt-get install -y build-essential
RUN apt-get install -y zlib1g-dev

#install Rosetta from github
#Anyone else using this needs to use their own login and token
RUN git clone https://github.com/RosettaCommons/rosetta.git

#move to src to source to install
WORKDIR rosetta/source

#checkout the branch that has the faster discovery code (will eventually push to main)

RUN git pull

#checkout most recent branch of code, in final version, will not have to do this
#RUN git checkout updated_ligand_discovery_search_code_ag_may_2024

#get python3.6 for scons.py installer
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.6 3

#install Rosetta wit hscons.py, create executables for motifs-related operations
RUN ./scons.py -j 1 mode=release identify_ligand_motifs extract_motifs ligand_discovery_search_protocol remove_duplicate_motifs

#remove extra data that is not needed for motifs discovery operations
#this helps keep memory overhead down on container size
RUN rm -drf tests src external

#move up
WORKDIR ..

#remove the large .git, since we won't need any of the github functionality and it is large
RUN rm -drf .git

#move to top of container
WORKDIR /

#get nano too
RUN apt-get install nano
