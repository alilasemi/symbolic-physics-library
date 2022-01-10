#!/bin/bash

cd generated
gcc -fPIC -shared -O3 e_and_cv.c -o e_and_cv.so
gcc -fPIC -shared -O3 e_s.c -o e_s.so
gcc -fPIC -shared -O3 cv_s.c -o cv_s.so
gcc -fPIC -shared -O3 wdot.c -o wdot.so
