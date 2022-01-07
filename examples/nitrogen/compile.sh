#!/bin/bash

cd generated
gcc -fPIC -shared -O3 e.c -o e.so
gcc -fPIC -shared -O3 cv.c -o cv.so
gcc -fPIC -shared -O3 e_s.c -o e_s.so
gcc -fPIC -shared -O3 cv_s.c -o cv_s.so
