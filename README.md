# CarrierSeq
Carrier/Pore-Maintainer Sequencing Analysis for Nanopore Sequencing

## About
...

## Getting Started

### Building Your Own Docker Image

1. Download & install Docker - https://www.docker.com/
2. Start Docker and increase the threads/memory setting to your preference. (The included script is set to 6 CPUs and 18 GB ram.)
3. Save Dockerfile to your directory 
4. ```$ cd to/your/directory```
5. ```$ docker build -t <name-your-image> .```

### Or Pull from Docker Hub

1. Start docker
2. run ```$ docker pull mojarro/carrierseq:latest```
