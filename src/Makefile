CC	= nvcc

%.o: %.c
	$(CC) -Xcompiler -fpic -g -c -o $@ $<

all: lib

test: matrix.o
	$(CC) -Xcompiler -fpic test.c matrix.o -lm -g -o test

clean:
	rm *.o *.so

lib: lib.c matrix.o
	$(CC) -Xcompiler -fpic -g -shared lib.c matrix.o -lm -g -o lib.so
