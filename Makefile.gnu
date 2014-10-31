
FC90 = gfortran
FC77 = gfortran

OPTS = -O2 -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -ffree-line-length-none

LIB = -L/etc/lapack-3.5.0 -llapack -lblas

