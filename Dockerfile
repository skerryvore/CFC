FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y gfortran libgdal-dev ssh libopenmpi-dev

WORKDIR /home
CMD gfortran --version




