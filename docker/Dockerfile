FROM java:7


#########################################
# samtools 0.1.19
#########################################

# MAINTAINER Kathrin Klee
RUN apt-get update && apt-get -y install samtools=0.1.19-1  && apt-get clean


#########################################
# bwa 0.7.12
#########################################

# Install compiler and perl stuff
RUN apt-get install --yes build-essential gcc-multilib apt-utils zlib1g-dev wget

# Get source code
WORKDIR /tmp

RUN wget -q http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2

RUN tar xjf bwa-0.7.12.tar.bz2

WORKDIR /tmp/bwa-0.7.12

# Patch Makefile
RUN sed -i 's/CFLAGS=\\t\\t-g -Wall -Wno-unused-function -O2/CFLAGS=-g -Wall -Wno-unused-function -O2 -static/' Makefile

# Compile
RUN make

RUN cp -p bwa /usr/local/bin

# Cleanup
RUN rm -rf /tmp/bwa-0.7.12

RUN apt-get clean

WORKDIR /root


#########################################
# shapeit.v2
#########################################

RUN mkdir /tmp/shapeit
WORKDIR /tmp/shapeit
RUN wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
RUN tar -zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
RUN cp -p bin/shapeit /usr/local/bin







# RUN apt-get clean
# RUN apt-get remove --yes --purge build-essential gcc-multilib apt-utils zlib1g-dev wget
