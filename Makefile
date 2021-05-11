.PHONY=clean

CXX?=g++

LIB=-lz
INCLUDE+=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -Ibonsai -I. -Isrc
OPT+=-std=c++17 -O3 -march=native -fopenmp
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts
EXTRA+=-DNOCURL

OFS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
OBJ=$(OFS)
OBJLD=$(patsubst %.o,%.ldo,$(OFS))
OBJF=$(patsubst %.o,%.fo,$(OFS))

all: dashing2 dashing2-ld dashing2-f32
obh: echo $(OBJ)

dashing2: $(OBJ) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-ld: $(OBJLD) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJLD) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-f32: $(OBJF) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJF) -o $@ $(LIB) $(EXTRA) libBigWig.a
readbw: src/bwsketch.o test/readbw.o
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< test/readbw.o -o readbw $(LIB) $(EXTRA) libBigWig.a
%.o: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA)
%.ldo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="long double"
%.fo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="float"

libBigWig.a: $(wildcard libBigWig/*.c) $(wildcard libBigWig/*.h)
	cd libBigWig && $(MAKE) && cp libBigWig.a ..

clean:
	rm -f dashing2 dashing2-ld dashing2-f32 libBigWig.a $(OBJ) $(OBJLD) $(OBJF)
