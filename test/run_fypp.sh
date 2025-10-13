#!/bin/sh -e
fypp -F -DWITH_REAL32  check.fypp check_real32.f90
fypp -F -DWITH_REAL64  check.fypp check_real64.f90
fypp -F -DWITH_REAL128 check.fypp check_real128.f90
