FROM rocker/r2u

WORKDIR /app
COPY requirements.* .

RUN ./requirements.sh

COPY *.py ./
COPY data /app/data

ENV PORT=8080

CMD uvicorn app:app --host 0.0.0.0 --port $PORT