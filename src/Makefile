CC	= nvcc

%.o: %.c
	$(CC) -Xcompiler -fpic -g -c -o $@ $<

all: lib

test: matrix.o
	$(CC) -Xcompiler -fPIC test.c matrix.o -lm -g -o test

rebuild: clean all

clean:
	rm *.o *.so

lib: lib.c matrix.o device.o
	$(CC) -Xcompiler -fpic -g -shared lib.c device.o matrix.o -lm -g -o lib.so
