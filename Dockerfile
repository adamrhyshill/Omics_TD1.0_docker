FROM continuumio/miniconda3:4.9.2
COPY ./requirements.txt /requirements.txt
RUN pip install -r requirements.txt
RUN git clone https://github.com/adamrhyshill/Omics_TD1.0_docker.git
WORKDIR /Omics_TD1.0_docker
