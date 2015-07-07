#!/bin/bash
#

compiler='g95'
rm -f example LDC3.o ldc3.mod
$compiler -O3 -c LDC3.f90
$compiler -O3 -o example example.f90 LDC3.o
