FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y build-essential wget unzip python2.7 \
					   python-dev git python-pip bats awscli curl \
					   libcurl4-openssl-dev make gcc zlib1g-dev

# Set the default langage to C
ENV LC_ALL C

# Use /share as the working directory
RUN mkdir /share
WORKDIR /share

# Add files
RUN mkdir /usr/map_viruses
ADD requirements.txt /usr/map_viruses

# Install python requirements
RUN pip install -r /usr/map_viruses/requirements.txt && rm /usr/map_viruses/requirements.txt


# Install DIAMOND v0.9.10
RUN cd /usr/map_viruses && \
	wget -q https://github.com/bbuchfink/diamond/releases/download/v0.9.10/diamond-linux64.tar.gz && \
	tar xzf diamond-linux64.tar.gz && \
	mv diamond /usr/bin/ && \
	rm diamond-linux64.tar.gz


# Install the SRA toolkit
RUN cd /usr/local/bin && \
	wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz && \
	tar xzf sratoolkit.2.8.2-ubuntu64.tar.gz && \
	ln -s /usr/local/bin/sratoolkit.2.8.2-ubuntu64/bin/* /usr/local/bin/ && \
	rm sratoolkit.2.8.2-ubuntu64.tar.gz


# Add the run script to the PATH
ADD map_viruses.py /usr/map_viruses
ADD make_viral_db.py /usr/map_viruses
ADD lib /usr/map_viruses/lib
RUN cd /usr/map_viruses && \
	chmod +x map_viruses.py && \
	ln -s /usr/map_viruses/map_viruses.py /usr/bin/  && \
	ln -s /usr/map_viruses/make_viral_db.py /usr/bin/


# Add a wrapper to help execution via SciLuigi
RUN apt-get install -y python3-pip
RUN pip3 install bucket_command_wrapper==0.1.0 


# Run tests and then remove the folder
ADD tests /usr/map_viruses/tests
RUN bats /usr/map_viruses/tests/ && rm -r /usr/map_viruses/tests/
