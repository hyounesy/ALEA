###########################################################################
# Dockerfile
#
# Software:         ALEA
# Software Version: 1.9.6
# Description:      A toolbox for Allele-specific Epigenomics Analysis
# Base Image:       java:7
# Build Cmd:        docker build -t /path/to/your/dockerfile-directory/
# Pull Cmd:         docker pull taskkoike/alea:1.9.6
# Run Cmd:          docker run -it taskkoike/alea:1.9.6
###########################################################################

FROM java:7

MAINTAINER Tasuku KOIKE <44215022@nodai.ac.jp>


###########################################################################
# Preparation
###########################################################################

RUN apt-get update && \
  apt-get install -y \
  wget \
  gzip \
  bzip2 \
  g++ \
  make \
  pkg-config \
  ncurses-dev \
  zlib1g-dev \
  build-essential \
  gcc-multilib \
  apt-utils \
  vim


###########################################################################
# samtools 0.1.19
###########################################################################

WORKDIR /tmp
RUN wget -q https://github.com/samtools/samtools/archive/0.1.19.tar.gz && \
  tar -xf 0.1.19.tar.gz -C /opt

WORKDIR /opt/samtools-0.1.19
RUN make all

RUN rm -f /tmp/0.1.19.tar.gz


###########################################################################
# samtools 1.3.1 & tabix & gbzip
###########################################################################

ENV SAMTOOLS_INSTALL_DIR=/opt/samtools-1.3.1
WORKDIR /tmp
RUN wget -q https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
  tar -xjf samtools-1.3.1.tar.bz2

WORKDIR /tmp/samtools-1.3.1
RUN ./configure --enable-plugins --prefix=$SAMTOOLS_INSTALL_DIR && \
  make all all-htslib && \
  make install install-htslib

RUN rm -f /tmp/samtools-1.3.1.tar.bz2 && \
  rm -rf /tmp/samtools-1.3.1


###########################################################################
# bwa 0.7.12
###########################################################################

WORKDIR /tmp
RUN wget -q http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2 && \
  tar -xjf bwa-0.7.12.tar.bz2 -C /opt

WORKDIR /opt/bwa-0.7.12
# Patch Makefile and Compile
RUN sed -i 's/CFLAGS=\\t\\t-g -Wall -Wno-unused-function -O2/CFLAGS=-g -Wall -Wno-unused-function -O2 -static/' Makefile && \
  make

RUN rm -f /tmp/bwa-0.7.12.tar.bz2


###########################################################################
# bowtie 1.1.2
###########################################################################

WORKDIR /tmp
RUN wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip && \
  unzip bowtie-1.1.2-linux-x86_64.zip -d /opt && \
  rm -f /tmp/bowtie-1.1.2-linux-x86_64.zip


###########################################################################
# bowtie2-2.2.9
###########################################################################

WORKDIR /tmp
RUN wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip && \
  unzip bowtie2-2.2.9-linux-x86_64.zip -d /opt && \
  rm -f /tmp/bowtie2-2.2.9-linux-x86_64.zip


###########################################################################
# bedGraphToBigWig
###########################################################################

ENV BEDGRAPHTOBIGWIG_VERSION=287
ENV BEDGRAPHTOBIGWIG_URL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v${BEDGRAPHTOBIGWIG_VERSION}/bedGraphToBigWig

WORKDIR /tmp
RUN wget -q ${BEDGRAPHTOBIGWIG_URL} && \
  chmod +x bedGraphToBigWig

RUN mkdir /opt/bedGraphToBigWig && \
  cp -p /tmp/bedGraphToBigWig /opt/bedGraphToBigWig/ && \
  rm -f /tmp/bedGraphToBigWig


###########################################################################
# shapeit.v2
###########################################################################

WORKDIR /tmp
RUN wget -q https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz && \
  mkdir /opt/shapeit && \
  tar -xf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz -C /opt/shapeit && \
  rm -f /tmp/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz


###########################################################################
# VCFtools 0.1.14
###########################################################################

ENV VCFTOOLS_INSTALL_DIR=/opt/vcftools-0.1.14
WORKDIR /tmp
RUN wget -q https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz && \
  tar -xf vcftools-0.1.14.tar.gz -C /opt

ENV PERL5LIB=$VCFTOOLS_INSTALL_DIR/src/perl
WORKDIR /opt/vcftools-0.1.14
RUN ./configure --prefix=$VCFTOOLS_INSTALL_DIR && \
  make && \
  make install

RUN rm -f /tmp/vcftools-0.1.14.tar.gz


###########################################################################
# Bismark 0.16.3
###########################################################################

WORKDIR /tmp
RUN wget -q http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.16.3.tar.gz && \
  tar -xf bismark_v0.16.3.tar.gz -C /opt && \
  rm -f bismark_v0.16.3.tar.gz


###########################################################################
# tophat2-2.1.1
###########################################################################

WORKDIR /tmp
RUN wget -q wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && \
	tar -xf tophat-2.1.1.Linux_x86_64.tar.gz -C /opt && \
	rm -f tophat-2.1.1.Linux_x86_64.tar.gz


###########################################################################
# STAR-2.4.0i
###########################################################################

WORKDIR /tmp
RUN wget -q https://github.com/alexdobin/STAR/archive/STAR_2.4.0i.tar.gz && \
  tar -xf /opt/star/STAR_2.4.0i.tar.gz -C /opt && \
  make STAR -C /opt/star/STAR-STAR_2.4.0i/source && \
  rm -f STAR_2.4.0i.tar.gz


###########################################################################
# ALEA 1.9.6
###########################################################################

ENV ALEA_INSTALL_DIR=/opt/alea
ENV ALEA_VERSION=1.9.6
ENV ALEA_FILENAME=alea.${ALEA_VERSION}.tar.gz
ENV ALEA_URL=https://github.com/hyounesy/ALEA/raw/Task/dist/${ALEA_FILENAME}

WORKDIR /tmp
RUN wget -q $ALEA_URL && \
  tar -xf ${ALEA_FILENAME} -C /opt && \
  rm -f /tmp/${ALEA_FILENAME}


###########################################################################
# Create links
###########################################################################

RUN ln -s /opt/samtools-0.1.19/samtools $ALEA_INSTALL_DIR/bin/ && \
  ln -s $SAMTOOLS_INSTALL_DIR/bin/tabix $ALEA_INSTALL_DIR/bin/tabix && \
  ln -s $SAMTOOLS_INSTALL_DIR/bin/bgzip $ALEA_INSTALL_DIR/bin/bgzip && \
  ln -s /opt/bwa-0.7.12/bwa $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-align-l $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-align-s $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-build $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-build-l $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-build-s $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-inspect $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-inspect-l $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie-1.1.2/bowtie-inspect-s $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2 $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-align-l $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-align-s $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-build $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-build-l $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-build-s $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-inspect $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-inspect-l $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bowtie2-2.2.9/bowtie2-inspect-s $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bedGraphToBigWig/bedGraphToBigWig $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/shapeit/bin/shapeit $ALEA_INSTALL_DIR/bin/ && \
  ln -s $VCFTOOLS_INSTALL_DIR/bin/* $ALEA_INSTALL_DIR/bin/


###########################################################################
# Clean
###########################################################################

WORKDIR /
RUN apt-get clean && \
  apt-get remove -y --force-yes --purge \
  wget \
  bzip2 \
  g++ \
  make \
  pkg-config \
  ncurses-dev \
  zlib1g-dev \
  build-essential \
  gcc-multilib \
  apt-utils


###########################################################################
# ENTRYPOINT
###########################################################################

#ENTRYPOINT ["/opt/alea/bin/alea"]
