FROM mgibio/samtools-cwl:1.0.0
MAINTAINER John Garza <johnegarza@wustl.edu>

LABEL \
    description="Versioned steps to produce a clinvar vcf for use in workflows"

RUN apt-get update -y && apt-get install -y \
    curl \
    libnss-sss \
    python-pip

RUN pip install --upgrade pip
RUN pip install unidecode

COPY ncbi_to_vcf.py /usr/bin/ncbi_to_vcf.py

