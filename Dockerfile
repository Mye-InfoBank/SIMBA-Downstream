FROM rocker/r2u

WORKDIR /app
COPY requirements.* .

RUN ./requirements.sh