Bootstrap: docker
From: ubuntu:18.04 
Stage: base-build

%post
# DAb-seq Dockerfile
# Ben Demaree 1.26.2021

# start with ubuntu base

### install linux dependencies
apt-get -y update
apt-get -y install python3
apt-get -y install python3-pip
ln -s /usr/bin/python3 /usr/bin/python & ln -s /usr/bin/pip3 /usr/bin/pip
apt-get install -y wget
apt-get install -y pigz
apt-get install -y libncurses5-dev
apt-get install -y zlib1g-dev
apt-get install -y libbz2-dev
apt-get install -y liblzma-dev
apt-get install -y libcurl4-openssl-dev
apt-get install -y unzip
apt-get install -y git
apt-get install -y default-jdk

### install bioinformatics programs
mkdir -p /dabseq/programs
cd /dabseq/programs

# htslib
wget -q --show-progress --progress=bar:force:noscroll https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -xjf htslib-1.9.tar.bz2
mkdir htslib
mkdir -p /dabseq/programs/htslib-1.9
cd /dabseq/programs/htslib-1.9
./configure --prefix=/dabseq/programs/htslib
make
make install
PATH="$PATH:/dabseq/programs/htslib/bin"

# samtools
mkdir -p /dabseq/programs
cd /dabseq/programs
wget -q --show-progress --progress=bar:force:noscroll https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
tar -xjf samtools-1.8.tar.bz2
mkdir samtools
mkdir -p /dabseq/programs/samtools-1.8
cd /dabseq/programs/samtools-1.8
./configure --prefix=/dabseq/programs/samtools
make
make install
PATH="$PATH:/dabseq/programs/samtools/bin"

# bcftools
mkdir -p /dabseq/programs
cd /dabseq/programs
wget -q --show-progress --progress=bar:force:noscroll https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -xjf bcftools-1.9.tar.bz2
mkdir bcftools
mkdir -p /dabseq/programs/bcftools-1.9
cd /dabseq/programs/bcftools-1.9
./configure --prefix=/dabseq/programs/bcftools
make
make install
PATH="$PATH:/dabseq/programs/bcftools/bin"

# gatk
mkdir -p /dabseq/programs
cd /dabseq/programs
wget -q --show-progress --progress=bar:force:noscroll https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip
unzip gatk-4.1.3.0.zip
PATH="$PATH:/dabseq/programs/gatk-4.1.3.0"

# bowtie2
wget -q --show-progress --progress=bar:force:noscroll https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip
unzip bowtie2-2.3.4.1-linux-x86_64.zip
PATH="$PATH:/dabseq/programs/bowtie2-2.3.4.1-linux-x86_64"

# bedtools
mkdir bedtools
mkdir -p /dabseq/programs/bedtools
cd /dabseq/programs/bedtools
wget -q --show-progress --progress=bar:force:noscroll https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
PATH="$PATH:/dabseq/programs/bedtools"
chmod 777 /dabseq/programs/bedtools/bedtools

# bbmap
mkdir -p /dabseq/programs
cd /dabseq/programs
wget -q --show-progress --progress=bar:force:noscroll https://sourceforge.net/projects/bbmap/files/BBMap_38.57.tar.gz
gunzip BBMap_38.57.tar.gz
tar -xf BBMap_38.57.tar
PATH="$PATH:/dabseq/programs/bbmap"

# cleanup installation files
rm *.zip
rm *.bz2
rm *.tar
rm -rf htslib-1.9/
rm -rf samtools-1.8/

### install missing python modules
pip3 install numpy==1.19.2
pip3 install scipy==1.5.2
pip3 install pandas==1.1.3
pip3 install tables==3.6.1
pip3 install regex==2020.10.15
pip3 install h5py==2.10.0
pip3 install matplotlib==3.3.2
pip3 install seaborn==0.11.0
pip3 install statsmodels==0.12.0
pip3 install scikit-allel==1.3.2
pip3 install cutadapt==2.10
pip3 install slackclient==2.9.2
pip3 install ipython>=7.4

# set input directories to fully open - so any user can write
mkdir /input
chmod -R 777 /input

# get pipeline files from dab-seq repo
mkdir -p /dabseq/pipeline
cd /dabseq/pipeline
git clone https://github.com/AbateLab/DAb-seq.git

%environment
export PATH="$PATH:/dabseq/programs/htslib/bin"
export PATH="$PATH:/dabseq/programs/samtools/bin"
export PATH="$PATH:/dabseq/programs/bcftools/bin"
export PATH="$PATH:/dabseq/programs/gatk-4.1.3.0"
export PATH="$PATH:/dabseq/programs/bowtie2-2.3.4.1-linux-x86_64"
export PATH="$PATH:/dabseq/programs/bedtools"
export PATH="$PATH:/dabseq/programs/bbmap"

Bootstrap: docker
From: base-build 
Stage: human-build

%post

# the following section is for builds with a HUMAN reference
########################################################################################################################
########################################################################################################################

mkdir -p /dabseq/programs
cd /dabseq/programs

# idtseek
git clone https://github.com/tommyau/itdseek
PATH="$PATH:/dabseq/programs/itdseek"

# snpeff
wget -q --show-progress --progress=bar:force:noscroll http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
PATH="$PATH:/dabseq/programs/snpEff"
PATH="$PATH:/dabseq/programs/snpEff/scripts"
# install snpEff database for hg19
snpEff download hg19

### download genome reference files
mkdir -p /dabseq/references
cd /dabseq/references

# get clinvar db
wget -q --show-progress --progress=bar:force:noscroll https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2020/clinvar_20200329.vcf.gz
gunzip clinvar_20200329.vcf.gz
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' clinvar_20200329.vcf > clinvar_20200329.chr.vcf
rm clinvar_20200329.vcf
bgzip -f -@ 16 clinvar_20200329.chr.vcf
tabix clinvar_20200329.chr.vcf.gz

# get hg19 fasta and pre-built indices
wget -q --show-progress --progress=bar:force:noscroll ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O hg19.fasta.gz
gunzip hg19.fasta.gz
# genrate using 'gatk CreateSequenceDictionary'
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/6yzq1n5fwtp06hf/hg19.dict?dl=0 -O hg19.dict
# generate using 'samtools faidx'
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/jbfyt5uz7jnfrm1/hg19.fasta.fai?dl=0 -O hg19.fasta.fai
# compress with 'tar -czvf [archive_name] [dir_to_compress]'
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/asqiuaiyzvqnkj0/hg19_bt2.tar.gz?dl=0 -O hg19_bt2.tar.gz
gunzip hg19_bt2.tar.gz
tar -xf hg19_bt2.tar
rm hg19_bt2.tar
mv hg19_bt2/*.bt2 .

########################################################################################################################
########################################################################################################################

%environment
export PATH="$PATH:/dabseq/programs/itdseek"
export PATH="$PATH:/dabseq/programs/snpEff"
export PATH="$PATH:/dabseq/programs/snpEff/scripts"

Bootstrap: docker
From: base-build 
Stage: hiv-build

%post

# the following section is for builds with a HUMAN and HIV (HXB2) reference
########################################################################################################################
########################################################################################################################

mkdir -p /dabseq/programs
cd /dabseq/programs

# idtseek
git clone https://github.com/tommyau/itdseek
PATH="$PATH:/dabseq/programs/itdseek"

# snpeff
wget -q --show-progress --progress=bar:force:noscroll http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
PATH="$PATH:/dabseq/programs/snpEff"
PATH="$PATH:/dabseq/programs/snpEff/scripts"
# install snpEff database for hg19
snpEff download hg19

### download genome reference files
mkdir -p /dabseq/references
cd /dabseq/references

# get clinvar db
wget -q --show-progress --progress=bar:force:noscroll https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2020/clinvar_20200329.vcf.gz
gunzip clinvar_20200329.vcf.gz
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' clinvar_20200329.vcf > clinvar_20200329.chr.vcf
rm clinvar_20200329.vcf
bgzip -f -@ 16 clinvar_20200329.chr.vcf
tabix clinvar_20200329.chr.vcf.gz

# get hg19/hxb2 fasta and pre-built indices
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/hby5s2elzu4do60/hg19_hxb2.fasta.gz?dl=0 -O hg19_hxb2.fasta.gz
gunzip hg19_hxb2.fasta.gz
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/3o2vrm8hqr67pbu/hg19_hxb2.dict?dl=0 -O hg19_hxb2.dict
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/pp87jvi6rnpm9t2/hg19_hxb2.fasta.fai?dl=0 -O hg19_hxb2.fasta.fai
wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/qnjazwnt08u0z14/hg19_hxb2_bt2.tar.gz?dl=0 -O hg19_hxb2_bt2.tar.gz
gunzip hg19_hxb2_bt2.tar.gz
tar -xf hg19_hxb2_bt2.tar
rm hg19_hxb2_bt2.tar
mv hg19_hxb2_bt2/*.bt2 .

########################################################################################################################
########################################################################################################################
%environment
export PATH="$PATH:/dabseq/programs/itdseek"
export PATH="$PATH:/dabseq/programs/snpEff"
export PATH="$PATH:/dabseq/programs/snpEff/scripts"
%runscript
cd /dabseq/references
exec /bin/bash "$@"
%startscript
cd /dabseq/references
exec /bin/bash "$@"
