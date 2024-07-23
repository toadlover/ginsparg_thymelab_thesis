FROM ubuntu:bionic-20220315
  
RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y git

#RUN apt-get install -y glibc

RUN apt-get install -y g++
RUN apt-get install -y gcc
RUN apt-get install -y build-essential
RUN apt-get install -y zlib1g-dev


RUN git clone https://toadlover:{my personal access token}@github.com/toadlover/conformator_for_container.git
WORKDIR /conformator_for_container
RUN tar -xvzf conformator_1.2.1_ubuntu-16.04-64bit.tar.gz
WORKDIR /conformator_for_container/conformator_1.2.1

#register for conformator here: https://software.zbh.uni-hamburg.de/
RUN ./conformator --license {conformator license key}

RUN apt-get install wget
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
#RUN sha256sum Anaconda3-2022.05-Linux-x86_64.sh
#RUN bash Anaconda3-2022.05-Linux-x86_64.sh --quiet

# Install base utilities
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

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

WORKDIR /conformator_for_container
RUN tar -xzvf openbabel-openbabel-2-4-0.tar.gz
RUN mkdir build
WORKDIR /conformator_for_container/build
RUN cmake ../openbabel-openbabel-2-4-0
RUN make
RUN make install


#RUN python setup.py build_ext -I/usr/local/include/openbabel3 -L/usr/local/lib
#RUN python setup.py install


#RUN pip install openbabel


#THIS IS THE WRONG PYBEL, NOT THE EXTENSION OF OPENBABEL
#RUN pip install pybel

#RUN conda install -c bioconda pybel



#RUN conda init bash
#RUN conda update anaconda
#RUN conda create --name pybel_env python=3
#RUN conda init /bin/bash
#RUN conda activate
#RUN conda activate pybel_env
#RUN conda install pybel
