#!/bin/bash
gfortran -c d1mach.f
gfortran -c drc3jj.f
gfortran -c fdump.f
gfortran -c i1mach.f
gfortran -c j4save.f
gfortran -c spline.for
gfortran -c splint.for
gfortran -c xercnt.f
gfortran -c xerhlt.f
gfortran -c xermsg.f
gfortran -c xerprn.f
gfortran -c xgetua.f
gfortran -c xersve.f
g++ -c cholesky_new_unbin_L_1_2048_planck_new.c 
g++ -o cholesky_new_unbin_L_1_2048_planck_new.exe cholesky_new_unbin_L_1_2048_planck_new.o d1mach.o drc3jj.o fdump.o i1mach.o j4save.o spline.o splint.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o -lgfortran
rm *.o
