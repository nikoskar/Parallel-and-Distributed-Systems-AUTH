v3serial: v3_serial.o mmio.o
	gcc v3_serial.o mmio.o -o v3serial

v3cilk: v3_cilk.o mmio.o
	gcc v3_cilk.o mmio.o -o v3cilk -fcilkplus

v3omp: v3_openmp.o mmio.o
	gcc v3_openmp.o mmio.o -o v3omp -fopenmp

v4serial: v4_serial.o mmio.o
	gcc v4_serial.o mmio.o -o v4serial

v4cilk: v4_cilk.o mmio.o
	gcc v4_cilk.o mmio.o -o v4cilk -fcilkplus

v3_serial.o: v3_serial.c
	gcc -c v3_serial.c

v3_cilk.o: v3_cilk.c
	gcc -c v3_cilk.c -fcilkplus

v3_openmp.o: v3_openmp.c
	gcc -c v3_openmp.c -fopenmp

v4_serial.o: v4_serial.c
	gcc -c v4_serial.c

v4_cilk.o: v4_cilk.c
	gcc -c v4_cilk.c -fcilkplus

v4_openmp.o: v4_openmp.c
	gcc -c v4_openmp.c -fopenmp

v4_pthreads.o: v4_pthreads.c
	gcc -c v4_pthreads.c -pthread

mmio.o: mmio.c mmio.h
	gcc -c mmio.c

clean:
	rm *.o v3serial v3cilk v3omp v4serial v4cilk v4omp v4pthreads
