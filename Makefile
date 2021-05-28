.PHONY=clean

CXX?=g++

LIB=-lz
INCLUDE+=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -Ibonsai -I. -Isrc
OPT+=-std=c++17 -O3 -march=native -fopenmp
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts -pedantic
EXTRA+=-DNOCURL
CXXFLAGS+=-std=c++17

OFS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
OBJ=$(OFS)
OBJLD=$(patsubst %.o,%.ldo,$(OFS))
OBJF=$(patsubst %.o,%.fo,$(OFS))

all: dashing2
unit: readfx readbw readbed
obh: echo $(OBJ)

all3d: dashing2 dashing2-f dashing2-ld


OBJFS=src/enums.cpp src/counter.cpp src/fastxsketch.cpp src/merge.cpp src/bwsketch.cpp src/bedsketch.cpp
LIBOBJ=$(patsubst %.cpp,%.o,$(OBJFS))
FLIBOBJ=$(patsubst %.cpp,%.fo,$(OBJFS))
LDLIBOBJ=$(patsubst %.cpp,%.ldo,$(OBJFS))

dashing2: $(OBJ) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-ld: $(OBJLD) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJLD) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-f: $(OBJF) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJF) -o $@ $(LIB) $(EXTRA) libBigWig.a
read%: test/read%.o $(LIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(LIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
read%-ld: test/read%.ldo $(LDLIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(LDLIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a -DDSKETCH_FLOAT_TYPE="long double"
read%-f: test/read%.fo $(FLIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(FLIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a -DSKETCH_FLOAT_TYPE="float"
%.o: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA)
%.ldo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="long double"
%.fo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="float"

libBigWig.a: $(wildcard libBigWig/*.c) $(wildcard libBigWig/*.h)
	cd libBigWig && sed -i '' 's/HAVE_CURL:/#/' Makefile && $(MAKE) && cp libBigWig.a ..

test: readfx readbw

clean:
	rm -f dashing2 dashing2-ld dashing2-f libBigWig.a $(OBJ) $(OBJLD) $(OBJF) readfx readfx-f readfx-ld readbw readbw readbw-f readbw-ld
