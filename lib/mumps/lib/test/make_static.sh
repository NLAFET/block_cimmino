#!/bin/bash

rm -f *.o
ar -x ../libdmumps.a
ar -rcs ../libdmumps_s.a *.o
