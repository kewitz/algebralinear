CC	= nvcc
gen = -gencode arch=compute_30,code=sm_30
flags = -Xcompiler -fPIC -g
#flags = -Xcompiler -fPIC -g -m 64

all: lib kernels.so

%.o: %.c
	$(CC) $(flags) -c -o $@ $< $(gen)

kernels.so:
	$(CC) $(flags) -shared kernels.cu -lm -g -o kernels.so -lcuda $(gen)

rebuild: clean all

clean:
	rm *.o *.so

lib: matrix.o device.o
	$(CC) $(flags) -shared lib.c *.o -lm -g -o lib.so -lcuda $(gen)
