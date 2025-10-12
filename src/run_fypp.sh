#!/bin/sh -e
fypp -F -DWITH_REAL32  elliptic_nome_fortran_real.fypp elliptic_nome_fortran_real32.f90
fypp -F -DWITH_REAL64  elliptic_nome_fortran_real.fypp elliptic_nome_fortran_real64.f90
fypp -F -DWITH_REAL128 elliptic_nome_fortran_real.fypp elliptic_nome_fortran_real128.f90
fypp -F                elliptic_nome_fortran.fypp      elliptic_nome_fortran.f90
