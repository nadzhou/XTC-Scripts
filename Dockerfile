FROM continuumio/anaconda3

WORKDIR /usr/src/app

COPY xtc_read.py .
COPY requirements.txt .

RUN conda install --file requirements.txt

CMD ["python", "./xtc_read.py"]
