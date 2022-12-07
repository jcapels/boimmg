FROM informaticsmatters/rdkit-python3-debian

USER root

#update linux and accessories
RUN apt-get -y update
RUN apt-get -y upgrade

#get pip
RUN apt-get install -y python-dev build-essential
RUN apt-get install python3-pip

COPY ./requirements.txt .
RUN pip3 install -r requirements.txt

#clean output; -y acepts everything
RUN apt-get -y autoclean
RUN apt-get -y clean