FROM continuumio/miniconda3:4.9.2
RUN apt-get update -y --allow-releaseinfo-change
RUN apt-get install gcc -y
RUN apt-get install g++ -y
COPY ./requirements.txt /requirements.txt
RUN pip install -r requirements.txt --ignore-installed ruamel.yaml
RUN git clone https://github.com/adamrhyshill/Omics_TD1.0_docker.git /Omics
WORKDIR /Omics
CMD ["gunicorn", "-b :8043","index:server"]
