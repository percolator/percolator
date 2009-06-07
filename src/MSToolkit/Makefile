ZLIB = adler32.o compress.o crc32.o deflate.o gzio.o inffast.o inflate.o infback.o inftrees.o trees.o uncompr.o zutil.o
RAMP = ramp_base64.o ramp.o
MSTOOLKIT = Spectrum.o MSReader.o MSObject.o
MSMAT = msmat.o msmat_base.o msmat_header_parser.o LazyMatrix.o crawutils.o
SQLITE = sqlite3.o 
CC = g++ -g
GCC = gcc -g
BOOST_LOCATION=/usr/local/include/boost
MSMAT_FLAGS=-D_NO_CEPHES  -DOLD_BOOST
#MSMAT_FLAGS=-D_NO_CEPHES  -DOLD_BOOST
CFLAGS = -O3 -static -I $(BOOST_LOCATION) -I. -I./zLib -I./RAMP -I./SQLite -I./msmat -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 $(MSMAT_FLAGS)


LIBS = -lm -lpthread -ldl

nosqlite: $(ZLIB) $(RAMP) $(MSTOOLKIT) $(MSMAT)
	ar rcs libmstoolkitlite.a $(ZLIB) $(RAMP) $(MSTOOLKIT) $(MSMAT)
	$(CC) $(CFLAGS) MSSingleScan.cpp -L. -lmstoolkitlite -o MSSingleScan
	$(CC) $(CFLAGS) MSConvertFile.cpp -L. -lmstoolkitlite -o MSConvertFile


all:  $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) $(MSMAT)
	ar rcs libmstoolkit.a $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) $(MSMAT)
	$(CC) $(CFLAGS) MSTDemo.cpp -L. -lmstoolkit $(LIBS) -o MSTDemo
	$(CC) $(CFLAGS) MSSingleScan.cpp -L. -lmstoolkit $(LIBS) -o MSSingleScan
	$(CC) $(CFLAGS) MSConvertFile.cpp -L. -lmstoolkit $(LIBS) -o MSConvertFile
	$(CC) $(CFLAGS) MSWarper.cpp -L. -lmstoolkit $(LIBS) -o MSWarper
	$(CC) $(CFLAGS) msmat_test.cpp -L. -lmstoolkit $(LIBS) -o msmat_test

demo:  $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) $(MSMAT)
	$(CC) $(CFLAGS) $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) MSTDemo.cpp $(LIBS) -o MSTDemo
	$(CC) $(CFLAGS) $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) MSSingleScan.cpp $(LIBS) -o MSSingleScan
	$(CC) $(CFLAGS) $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) MSConvertFile.cpp $(LIBS) -o MSConvertFile


lib: $(ZLIB) $(RAMP) $(MSTOOLKIT) $(SQLITE) $(MSMAT)
	ar rcs libmstoolkit.a $(ZLIB) $(RAMP) $(SQLITE) $(MSTOOLKIT) $(MSMAT)

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

#SQLite object
sqlite3.o : ./SQLite/sqlite3.c
	$(GCC) $(CFLAGS) SQLite/sqlite3.c -c

#MSToolkit objects

Spectrum.o : Spectrum.cpp
	$(CC) $(CFLAGS) Spectrum.cpp -c

MSReader.o : MSReader.cpp
	$(CC) $(CFLAGS) MSReader.cpp -c

MSObject.o : MSObject.cpp
	$(CC) $(CFLAGS) MSObject.cpp -c

#MSMAT objects

msmat.o :  ./msmat/msmat.cpp ./msmat/msmat.H
	$(CC) -o $@ $(CFLAGS) -c $<

msmat_base.o : ./msmat/msmat_base.cpp ./msmat/msmat_base.H 
	$(CC) -o $@ $(CFLAGS) -c $<

msmat_header_parser.o : ./msmat/msmat_header_parser.cpp ./msmat/msmat_header_parser.H 
	$(CC) -o $@ $(CFLAGS) -c $<

LazyMatrix.o : ./msmat/LazyMatrix.cpp ./msmat/LazyMatrix.H 
	$(CC) -o $@ $(CFLAGS) -c $<

crawutils.o : ./msmat/crawutils.cpp ./msmat/crawutils.H 
	$(CC) -o $@ $(CFLAGS) -c $<
