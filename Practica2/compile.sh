#!/bin/sh

g++ -O3 main.cpp pattmatch.cpp -o main -Wall -W -ansi -pedantic -Dcimg_use_vt100 -I/home/ferran/AMDAPPSDK-3.0/include/CL  -lm -L/home/ferran/AMDAPPSDK-3.0/lib/x86_64 -lpthread -lX11 -lOpenCL