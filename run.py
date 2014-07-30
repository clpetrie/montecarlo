import os

os.system('gfortran -c random.f90')
os.system('gfortran -c mcinteg.f90')
os.system('gfortran random.o mcinteg.o -o mcinteg')
os.system('./mcinteg')
