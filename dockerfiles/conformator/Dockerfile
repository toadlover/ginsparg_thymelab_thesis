FROM ubuntu:bionic-20220315

#package installs
RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y git
RUN apt-get install -y g++
RUN apt-get install -y gcc
RUN apt-get install -y build-essential
RUN apt-get install -y zlib1g-dev

#get copy of conformator from conformator_for_container repo
#this is a personal private repository, and anyone trying to use this should use other means to get their copy of conformator into the container
#this repository also contains Rosetta python packages to use the molfile_to_params.py script to generate params files
RUN git clone https://github.com/toadlover/conformator_for_container.git

#move into conformator_for_container repo
WORKDIR /conformator_for_container

#unzip conformator
RUN tar -xvzf conformator_1.2.1_ubuntu-16.04-64bit.tar.gz

#move into conformator directory
WORKDIR /conformator_for_container/conformator_1.2.1

#activate conformator with license key 
#register for conformator here: https://software.zbh.uni-hamburg.de/
RUN ./conformator --license {conformator license key}

#get wget to get anaconda
RUN apt-get install wget

# Install base utilities
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

#prepare conda and pip
RUN conda update conda
RUN conda install pip

#build openbabel
#get cmake
#https://askubuntu.com/questions/356014/package-cmake-has-no-installation-candidate-ubuntu-12-04
RUN apt-get update -y
RUN apt-get update
RUN apt-get install cmake -y

#build openbabel
#https://open-babel.readthedocs.io/en/latest/Installation/install.html#compiling-open-babel
#openbabel package is in conformator_for_container repo
WORKDIR /conformator_for_container
RUN tar -xzvf openbabel-openbabel-2-4-0.tar.gz
RUN mkdir build
WORKDIR /conformator_for_container/build
RUN cmake ../openbabel-openbabel-2-4-0
RUN make
RUN make install
