#!/bin/bash

# Select Fortran compiler.
f90=gfortran

# Choose compiler options.
case "$f90" in
mpif90*|mpgfortran*|gfortran*)
  if [ "$1" = -d ] ; then
    opts="-Wall -Wextra -fimplicit-none -O0 -fbounds-check -g -pg -pedantic\
       -fbacktrace -fcray-pointer"
  else
    opts="-O2 -fprotect-parens -march=native -fcray-pointer"
  fi ;;
mpifort*|mpiifort*|ifort*)
  if [ "$1" = -d ] ; then
    opts="-check all -check noarg_temp_created -fp-stack-check -g -O0\
       -implicitnone -std95 -traceback -warn all,nounused -debug all -ftrapuv\
       -assume buffered_io"
  else
    opts="-O3 -no-prec-div -no-prec-sqrt -funroll-loops -no-fp-port -ip\
       -complex-limited-range -assume protect_parens -assume buffered_io"
  fi ;;
mpnagfor*|nagfor*)
  if [ "$1" = -d ] ; then
    opts="-g -C=all -colour -f95 -strict95 -u -v -gline -nan\
       -no_underflow_warning -w=uda -O0"
  else
    opts="-O4 -Ounsafe"
  fi ;;
*)
  echo "Compiler '$f90' not configured."
  exit 1 ;;
esac

# Compile.
$f90 $opts -o vib_levels vib_levels.f90 -llapack 2>&1
rm -f *.mod *.o
