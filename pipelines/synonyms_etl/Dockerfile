FROM python:3.8-slim-buster


#update linux and accessories
RUN apt-get -y update
RUN apt-get -y upgrade

#get pip
RUN apt-get install -y python-dev build-essential
RUN apt-get install -y python-pip
RUN python -m pip install --upgrade pip
RUN apt-get install -y libxrender1
RUN apt install libxext6
RUN apt-get install -y libxrender-dev

# create directory in container
WORKDIR /code

# copy requirements to directory created
COPY pipelines/synonyms_etl/requirements.txt /code/requirements.txt


# run and update requirements
RUN python -m pip install --no-cache-dir -r /code/requirements.txt
RUN pip install -U flask


# copy all folders in local root to the image
COPY ./ .

# install setup.py
RUN python "setup.py" install
RUN python "pipelines/synonyms_etl/setup_database_credentials.py"

# run script test_etl.py with python
CMD ["python","pipelines/synonyms_etl/etl_pipeline.py"]
