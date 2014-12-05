CC	= nvcc

build:
	$(CC) -Xcompiler -fpic -shared src/matrix.c -lm -g -o matrix.so
