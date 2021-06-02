############################################################
# Dockerfile to build Genotype imputation tools
# Based on Ubuntu 16.04
############################################################
# Set the base image to Ubuntu
FROM ubuntu:16.04
# File Author / Maintainer
LABEL maintainer="Ruth Nanjala rnanjala@icipe.org"
LABEL description="This is custom Docker Image for \
MHC Imputation tools and their dependencies."
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
# Install IMPUTE2
RUN wget http://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
  tar -zxvf impute_v2.3.2_x86_64_static.tgz && \
  mv impute_v2.3.2_x86_64_static/impute2 /usr/local/bin/impute2 && \
  mkdir /opt/impute2/example -p && \
  mv impute_v2.3.2_x86_64_static/Example/* /opt/impute2/example && \
  rm -rf impute_v2.3.2_x86_64_static impute_v2.3.2_x86_64_static.tgz
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
  make && \
  mv ./bin/Minimac3 /usr/local/bin/minimac3 && \
  cd .. && rm -r Minimac3
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
# Install Eagle
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
  gunzip Eagle_v2.4.1.tar.gz && \
  tar xvf Eagle_v2.4.1.tar && \
  mv Eagle_v2.4.1/eagle /usr/local/bin/ && \
  rm -rf Eagle_v2.4.1
# Install R and R packages
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
  /bin/bash ~/miniconda.sh -b -p /opt/conda && \
  rm ~/miniconda.sh && \
  echo "conda activate base" >> ~/.bashrc
RUN conda clean --all --yes && \
  conda install -y -c conda-forge r-base r-ggsci && \
  conda install -y -c bioconda r-ggplot2 r-dplyr r-plyr r-tidyr r-data.table r-reshape2 r-optparse r-sm
# HIBAG tool
RUN conda clean --all --yes && \
  conda install -c bioconda bioconductor-hibag
#Install DEEPHLA
RUN git clone https://github.com/tatsuhikonaito/DEEP-HLA.git && \
  mv DEEP-HLA /usr/local/bin/
#install csh
RUN conda clean --all --yes && \
  conda install -c conda-forge tcsh
#install svn  
RUN conda clean --all --yes && \
  conda install -c anaconda svn
# Install SNP2HLA
# chmod a+x /usr/local/bin/SNP2HLA_package_v1.0.3 && \
RUN git clone https://github.com/nanjalaruth/MHC-Imputation-Accuracy.git && \
  cd MHC-Imputation-Accuracy/ && \
  tar -xvzf SNP2HLA_package_v1.0.3.tar.gz && \ 
  mv -f SNP2HLA_package_v1.0.3/MakeReference/* /usr/local/bin && \
  mv -f SNP2HLA_package_v1.0.3/SNP2HLA/* /usr/local/bin && \
  chmod -R a+rwx /usr/local/bin/* && \  
  rm -fr SNP2HLA_package_v1.0.3.tar.gz SNP2HLA_package_v1.0.3 
WORKDIR /usr/local/bin/
RUN wget http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip && \
  unzip plink-1.07-x86_64.zip && \
  mv -f plink-1.07-x86_64/plink . && \
  rm -fr plink-1.07-x86_64.zip plink-1.07-x86_64 && \
  wget http://faculty.washington.edu/browning/beagle/recent.versions/beagle_3.0.4_05May09.zip && \
  unzip beagle_3.0.4_05May09.zip && \
  mv -f beagle.3.0.4/beagle.jar . && \
  rm -fr beagle_3.0.4_05May09.zip beagle.3.0.4 && \
  wget http://faculty.washington.edu/browning/beagle_utilities/linkage2beagle.jar && \
  wget http://faculty.washington.edu/browning/beagle_utilities/beagle2linkage.jar && \
  wget https://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar && \
  mv beagle.18May20.d20.jar beagle5.jar && \
  wget https://faculty.washington.edu/browning/beagle_utilities/vcf2beagle.jar && \
  wget https://faculty.washington.edu/browning/beagle_utilities/beagle2vcf.jar && \
  wget https://faculty.washington.edu/browning/beagle_utilities/transpose.jar
#install nano
RUN conda clean --all --yes && \
  conda install -c conda-forge nano
#install java
RUN conda clean --all --yes && \
  conda install -c bioconda java-jdk 
#MACH 
RUN wget http://csg.sph.umich.edu/abecasis/mach/download/mach.1.0.18.Linux.tgz && \
  tar -xvzf mach.1.0.18.Linux.tgz && \
  mv executables/mach1 /usr/local/bin/ && \
  rm -fr mach.1.0.18.Linux.tgz examples executables README
#plink 1.9
RUN conda clean --all --yes && \
  conda install -c bioconda plink
#python 3.6
RUN conda clean --all --yes && \
  conda install -c anaconda python=3.6
#pandas 0.25.3
RUN conda clean --all --yes && \
  conda install pandas=0.25.3
#perl 5.26.2
RUN conda clean --all --yes && \
  conda install -c anaconda perl=5.26.2
#rbase
RUN conda clean --all --yes && \
  conda install -y -c conda-forge r-base=3.6
#pyliftover 0.4
RUN conda clean --all --yes && \
  conda install -c bioconda pyliftover=0.4
#CookHLA
RUN git clone https://github.com/WansonChoi/CookHLA.git && \
  cd CookHLA && \
  rm -fr example 1000G_REF img README.md CookHLA_LINUX.yml CookHLA_OSX.yml dependency && \
  mv -f * /usr/local/bin && \
  chmod -R a+rwx /usr/local/bin/CookHLA.py && \ 
  chmod -R a+rwx /usr/local/bin/CookHLA_lab.py && \
  chmod -R a+rwx /usr/local/bin/CookHLA_lab_bglv5.py && \
  chmod -R a+rwx /usr/local/bin/MakeGeneticMap/* && \
  chmod -R a+rwx /usr/local/bin/measureAcc/__init__.py && \
  chmod -R a+rwx /usr/local/bin/measureAcc/__main__.py && \
  chmod -R a+rwx /usr/local/bin/measureAcc/measureAccuracy.py && \
  chmod -R a+rwx /usr/local/bin/measureAcc/NomenCleaner/* && \
  chmod -R a+rwx /usr/local/bin/measureAcc/src/* && \
  chmod -R a+rwx /usr/local/bin/measureAcc/data/* && \
  chmod -R a+rwx /usr/local/bin/src/* && \
  cd .. && \
  rm -fr CookHLA
RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu
USER ubuntu
CMD ["/bin/bash","-i"]
