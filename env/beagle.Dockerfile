FROM openjdk:8-jre
MAINTAINER Mike Sauria <mike.sauria@jhu.edu>

RUN wget https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar -O /beagle.28Jun21.220.jar && \
echo "#!/bin/bash" > /usr/local/sbin/beagle && \
echo "java -Xmx98g -jar /beagle.28Jun21.220.jar \$*" >> /usr/local/sbin/beagle && \
chmod a+rx /usr/local/sbin/beagle

RUN apt-get --allow-releaseinfo-change update && \
apt-get install -y libbz2-dev libvcflib-tools libvcflib-dev procps autoconf automake make gcc \
    perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev && \
rm -rf /var/lib/apt/lists/* && \
    wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2 && \
tar -xjvf bcftools.tar.bz2 && \
cd bcftools-1.3.1 && \
make && \
make prefix=/usr/local/bin install && \
mv /usr/local/bin/bin/bcftools /usr/bin/bcftools && \
rm -rf /usr/local/bin/bin && \
cd /usr/local/bin && \
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
tar -vxjf htslib-1.9.tar.bz2 && \
cd htslib-1.9 && \
make && \
mv bgzip ../ && \
cd ../ && \
rm -rf htslib-1.9
