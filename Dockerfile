############################################################
# Dockerfile to build MHC Imputation tools
# Based on Ubuntu 18.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:18.04

# File Author 
LABEL maintenair="Ruth Nanjala rnanjala@icipe.org"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

################## BEGIN INSTALLATION ######################

# Install wget
RUN apt-get update && apt-get install -y \
  autoconf \
  build-essential \
  git \
  libncurses5-dev \
  pkg-config \
  unzip \
  wget curl \
  python python-dev \
  libbz2-dev \
  liblzma-dev \
  zlib1g-dev &&\
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#Install IMPUTE4
RUN wget https://www.dropbox.com/sh/k6b34fzw9w4s8bg/AADRaxyaRFsDn6P5InT5qMiga/impute4.1.2_r300.3?dl=0

# Install SNP2HLA
RUN wget http://software.broadinstitute.org/mpg/snp2hla/data/SNP2HLA_package_v1.0.3.tar.gz && \
  tar -xzvf SNP2HLA_package_v1.0.3.tar.gz && \
  cd SNP2HLA_package_v1.0.3 && \
  make && \
  mv ./bin/SNP2HLA_package_v1.0.3 /usr/local/bin/SNP2HLA

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
  tar -xvf htslib-1.9.tar.bz2 && \
  cd htslib-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf htslib-1.9*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
  tar -xvf samtools-1.9.tar.bz2 && \
  cd samtools-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf samtools-1.9*

# Install VCFTools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
  tar -xvf vcftools-0.1.16.tar.gz && \
  cd vcftools-0.1.16 && \
  ./configure && \
  make && \
  make install

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
  tar -xvf bcftools-1.9.tar.bz2 && \
  cd bcftools-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf bcftools-1.9*

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
  tar -zxvf bedtools-2.28.0.tar.gz && \
  cd bedtools2 && \
  make && \
  mv bin/bedtools /usr/local/bin/ && \
  cd .. && rm -r bedtools2

# Install minimac3
RUN wget ftp://share.sph.umich.edu/minimac3/Minimac3.v2.0.1.tar.gz  && \
  tar -xzvf Minimac3.v2.0.1.tar.gz  && \
  cd Minimac3/  && \
  make  && \
  mv ./bin/Minimac3 /usr/local/bin/minimac3

# Install minimac4
RUN wget http://debian.mirror.ac.za/debian/pool/main/libs/libstatgen/libstatgen0_1.0.14-5_amd64.deb && \
  dpkg -i libstatgen0_1.0.14-5_amd64.deb && \
  wget https://github.com/statgen/Minimac4/releases/download/v1.0.2/minimac4-1.0.2-Linux.deb && \
  dpkg -i minimac4-1.0.2-Linux.deb

# Install PLINK2
# there is an undocumented stable url (without the date)
RUN wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O plink.zip && \
  unzip plink.zip -d /usr/local/bin/ && \
  rm -f plink.zip


# Install R and R packages
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
  /bin/bash ~/miniconda.sh -b -p /opt/conda && \
  rm ~/miniconda.sh && \
  echo "conda activate base" >> ~/.bashrc
RUN conda clean --all --yes && \
  conda install -y -c bioconda r-ggplot2 r-dplyr r-plyr r-tidyr r-data.table r-reshape2 r-optparse r-sm
RUN conda clean --all --yes && \
  conda install -y -c conda-forge r-ggsci
RUN conda clean --all --yes && \
  conda install -c bioconda bioconductor-hibag

# Change ownership
RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

# Run the specified command within the container
CMD ["/bin/bash","-i"]
