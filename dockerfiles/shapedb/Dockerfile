FROM ubuntu:bionic-20220315

#Package installs
RUN apt-get update
RUN apt-get install -y python3.8
RUN apt-get install -y git
RUN apt-get install -y g++
RUN apt-get install -y gcc
RUN apt-get install -y build-essential
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y nano
RUN apt install -y  autoconf automake libtool ghostscript liblemon-dev libeigen3-dev libann-dev bmagic libcgicc-dev libgoogle-perftools-dev libglpk-dev coinor-* libjsoncpp-dev cmake libboost-dev swig python-dev libxml2-dev  libcairo2-dev libboost-all-dev libcurl4-openssl-dev

#get fastcgi
RUN git clone https://github.com/FastCGI-Archives/fcgi2.git
WORKDIR /fcgi2
RUN autoreconf -i
RUN ./configure
RUN make
RUN make install
WORKDIR /

#get openbabel
RUN git clone https://github.com/openbabel/openbabel.git
WORKDIR openbabel
RUN mkdir build
WORKDIR build
RUN cmake .. -DPYTHON_BINDINGS=1 -DRUN_SWIG=1  -DWITH_MAEPARSER=0 -DWITH_COORDGEN=0
RUN make -j12
RUN make install
WORKDIR /

#get smina
RUN git clone http://git.code.sf.net/p/smina/code smina
WORKDIR smina
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make -j12
RUN cp libsmina.a /usr/local/lib/
WORKDIR /

#get pharmit/shapedb
RUN git clone http://git.code.sf.net/p/pharmit/code pharmit
WORKDIR pharmit/src
#RUN ls
RUN mkdir build
#WORKDIR build
RUN cp /smina/src/lib/SminaConverter.h /pharmit/src

RUN sed 's/#include \"parsing.h\"/#include \"\/smina\/src\/lib\/parsing.h\"/g' SminaConverter.h > output
RUN sed 's/#include \"PDBQTUtilities.h\"/#include \"\/smina\/src\/lib\/PDBQTUtilities.h\"/g' output > output2

RUN cp output2 SminaConverter.h
WORKDIR build
RUN cmake .. -DSMINALIB=/usr/local/lib/libsmina.a -DSMINA_DIR=$HOME/git/smina -DLEMON_DIR=/usr/lib/cmake/ 
RUN make shapedb -j12

WORKDIR /opt
ENV PATH="/opt/miniconda3/bin:${PATH}"
ARG PATH="/opt/miniconda3/bin:${PATH}"
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /opt/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN /opt/miniconda/bin/conda --version
RUN /opt/miniconda/bin/conda install -c conda-forge openbabel -y
RUN /opt/miniconda/bin/conda install numpy
