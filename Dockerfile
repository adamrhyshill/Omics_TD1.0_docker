FROM continuumio/miniconda3:main
COPY ./requirements.txt /requirements.txt
RUN pip install -r requirements.txt
RUN git clone https://github.com/adamrhyshill/Omics_TD1.0_docker.git
WORKDIR /Omics_TD1.0_docker
