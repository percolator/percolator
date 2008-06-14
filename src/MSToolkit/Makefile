ZLIB = adler32.o compress.o crc32.o deflate.o gzio.o inffast.o inflate.o infback.o inftrees.o trees.o uncompr.o zutil.o
RAMP = ramp_base64.o ramp.o
MSTOOLKIT = Spectrum.o MSReader.o MSObject.o 
CC = g++
GCC = gcc
CFLAGS = -O3 -static -I./zLib -I./RAMP -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

all:  MSTDemo.cpp MSSingleScan.cpp MSConvertFile.cpp $(ZLIB) $(RAMP) $(MSTOOLKIT)
	ar rcs libmstoolkit.a $(ZLIB) $(RAMP) $(MSTOOLKIT)
	$(CC) $(CFLAGS) MSTDemo.cpp -L. -lmstoolkit -o MSTDemo
	$(CC) $(CFLAGS) MSSingleScan.cpp -L. -lmstoolkit -o MSSingleScan
	$(CC) $(CFLAGS) MSConvertFile.cpp -L. -lmstoolkit -o MSConvertFile

demo:  MSTDemo.cpp MSSingleScan.cpp MSConvertFile.cpp $(ZLIB) $(RAMP) $(MSTOOLKIT) 
	$(CC) $(CFLAGS) $(ZLIB) $(RAMP) $(MSTOOLKIT) MSTDemo.cpp -o MSTDemo
	$(CC) $(CFLAGS) $(ZLIB) $(RAMP) $(MSTOOLKIT) MSSingleScan.cpp -o MSSingleScan
	$(CC) $(CFLAGS) $(ZLIB) $(RAMP) $(MSTOOLKIT) MSConvertFile.cpp -o MSConvertFile

lib: $(ZLIB) $(RAMP) $(MSTOOLKIT)
	ar rcs libmstoolkit.a $(ZLIB) $(RAMP) $(MSTOOLKIT)

clean:
	rm *.o MSTDemo MSSingleScan MSConvertFile libmstoolkit.a

# zLib objects

adler32.o : ./zLib/adler32.c
	$(GCC) $(CFLAGS) zLib/adler32.c -c

compress.o : ./zLib/compress.c
	$(GCC) $(CFLAGS) zLib/compress.c -c

crc32.o : ./zLib/crc32.c
	$(GCC) $(CFLAGS) zLib/crc32.c -c

deflate.o : ./zLib/deflate.c
	$(GCC) $(CFLAGS) zLib/deflate.c -c

gzio.o : ./zLib/gzio.c
	$(GCC) $(CFLAGS) zLib/gzio.c -c

inffast.o : ./zLib/inffast.c
	$(GCC) $(CFLAGS) zLib/inffast.c -c

inflate.o : ./zLib/inflate.c
	$(GCC) $(CFLAGS) zLib/inflate.c -c

infback.o : ./zLib/infback.c
	$(GCC) $(CFLAGS) zLib/infback.c -c

inftrees.o : ./zLib/inftrees.c
	$(GCC) $(CFLAGS) zLib/inftrees.c -c

trees.o : ./zLib/trees.c
	$(GCC) $(CFLAGS) zLib/trees.c -c

uncompr.o : ./zLib/uncompr.c
	$(GCC) $(CFLAGS) zLib/uncompr.c -c

zutil.o : ./zLib/zutil.c
	$(GCC) $(CFLAGS) zLib/zutil.c -c

#RAMP objects

ramp_base64.o : ./RAMP/ramp_base64.cpp
	$(CC) $(CFLAGS) RAMP/ramp_base64.cpp -c

ramp.o : ./RAMP/ramp.cpp
	$(CC) $(CFLAGS) RAMP/ramp.cpp -c

#MSToolkit objects

Spectrum.o : Spectrum.cpp
	$(CC) $(CFLAGS) Spectrum.cpp -c

MSReader.o : MSReader.cpp
	$(CC) $(CFLAGS) MSReader.cpp -c

MSObject.o : MSObject.cpp
	$(CC) $(CFLAGS) MSObject.cpp -c

