# SIMBA🦁 Downstream toolkit
This repository contains the source code of the python shiny app which is recommended to be used for the downstream analysis of the output of the SIMBA🦁 pipeline. Available instances can be found at at [https://exbio.wzw.tum.de/cost/overview/](https://exbio.wzw.tum.de/cost/overview/).

## Offline Usage
### Docker
The app can be run using the docker image `bigdatainbiomedicine/simba-downstream`. The image can be pulled from the docker hub using the following command:
```bash
docker run -p <port>:8080 bigdatainbiomedicine/simba-downstream
```
Assuming that the port 1234 is used, the app can be accessed at [http://localhost:1234](http://localhost:1234).

### Local
The app can also be run locally. The following steps are required to run the app locally:
1. Clone the repository:
```bash
git clone https://github.com/Mye-InfoBank/SIMBA-Downstream.git
```
2. Install the required packages:
```bash
pip install -r requirements.txt
```
3. Run the app:
```bash
cd src
uvicorn app:app --host 0.0.0.0 --port <port>
```
Assuming that the port 1234 is used, the app can be accessed at [http://localhost:1234](http://localhost:1234).
