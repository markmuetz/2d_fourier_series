CC = gfortran 

all: data/a.bin fourier_series_2d test

data/a.bin :  gen_array.py
	python gen_array.py

fourier_series_2d :  fourier_series_2d.F90
	$(CC) -o fourier_series_2d fourier_series_2d.F90

test: fourier_series_2d
	./fourier_series_2d
