
output: graph.o mmio.o
	gcc graph.o mmio.o -o output

graph.o: graph.c
	gcc -c graph.c

coo2csr_2: coo2csr_2.c
	gcc coo2csr_2 -o coo2csr_2

example: example.c
	gcc example.c -o example

mmio.o: mmio.c
	gcc -c mmio.c

to_coo: to_coo.c
	gcc to_coo.c -o to_coo

read_mtx: readmtx.c
	gcc readmtx.c -o readmtx

clean: 
	rm *.o

cleanall:
	rm 
