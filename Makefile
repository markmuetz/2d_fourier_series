CC = gfortran 

all: fourier_series_2d

fourier_series_2d :  fourier_series_2d.F90
	$(CC) -o fourier_series_2d fourier_series_2d.F90

run: fourier_series_2d
	./fourier_series_2d
