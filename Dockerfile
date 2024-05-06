FROM continuumio/miniconda3:24.3.0-0

WORKDIR /app
COPY environment.yml ./

RUN conda env create -f environment.yml

COPY *.py ./
COPY dgea dgea
COPY data data
COPY tree tree

ENV PORT=8080

CMD conda run -n shiny uvicorn app:app --host 0.0.0.0 --port $PORT