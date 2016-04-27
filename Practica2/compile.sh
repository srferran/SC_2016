#!/bin/sh

g++ -O3 main.cpp pattmatch.cpp -o main -Wall -W -ansi -pedantic -Dcimg_use_vt100 -I/usr/X11R6/include  -lm -L/usr/X11R6/lib -lpthread -lX11 -I/usr/local/cuda/include -L/usr/local/cuda/lib -lOpenCL

