# CC = gfortran 
CC = ifort 

all: fourier_series_2d

fourier_series_2d :  fourier_series_2d.F90
	$(CC) -O3 -o fourier_series_2d fourier_series_2d.F90

run: fourier_series_2d
	./fourier_series_2d
