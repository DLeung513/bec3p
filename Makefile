all:	bec3p bec3pc

bec3p:	bec3p.f90 parameters3.f90
	f95 parameters3.f90 bec3p.f90 -o bec3p

bec3pc:	bec3p.c++ parameters3.h
	c++ -O3 bec3p.c++ -o bec3pc
