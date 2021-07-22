.PHONY=clean

CXX?=g++

LIB=-lz
INCLUDE+=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -Ibonsai -I. -Isrc
OPT+=-std=c++17 -O3 -march=native -fopenmp
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts -pedantic -Wno-shift-count-overflow
EXTRA+=-DNOCURL
CXXFLAGS+=-std=c++17

OFS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
OBJ=$(OFS)
OBJLD=$(patsubst %.o,%.ldo,$(OFS))
OBJF=$(patsubst %.o,%.fo,$(OFS))
OBJDBG=$(patsubst %.o,%.do,$(OFS))
OBJADD=$(patsubst %.o,%.sano,$(OFS))
OBJLTO=$(patsubst %.o,%.lto,$(OFS))
OBJ0=$(patsubst %.o,%.0,$(OFS))

all: dashing2 dashing2-f dashing2-ld
unit: readfx readbw readbed
obh: echo $(OBJ)

all3d: dashing2 dashing2-f dashing2-ld
SEDSTR=
ifeq ($(shell uname -s ),Darwin)
    SEDSTR = " '' "
endif


OBJFS=src/enums.cpp src/counter.cpp src/fastxsketch.cpp src/merge.cpp src/bwsketch.cpp src/bedsketch.cpp src/fastxsketchbyseq.cpp src/bwreduce.cpp
LIBOBJ=$(patsubst %.cpp,%.o,$(OBJFS))
LIB0=$(patsubst %.cpp,%.0,$(OBJFS))
DLIBOBJ=$(patsubst %.cpp,%.do,$(OBJFS))
GLIBOBJ=$(patsubst %.cpp,%.go,$(OBJFS))
FLIBOBJ=$(patsubst %.cpp,%.fo,$(OBJFS))
LDLIBOBJ=$(patsubst %.cpp,%.ldo,$(OBJFS))

dashing2: $(OBJ) libBigWig.a $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a -DNDEBUG -flto
dashing2-0: $(OBJ0) libBigWig.a $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ0) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-d: $(OBJDBG) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJDBG) -o $@ $(LIB) $(EXTRA) libBigWig.a -O0
dashing2-d0: $(OBJDBG) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJDBG) -o $@ $(LIB) $(EXTRA) libBigWig.a -O0
dashing2-add: $(OBJADD) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJADD) -o $@ $(LIB) $(EXTRA) libBigWig.a -fsanitize=address
dashing2-g: $(GLIBOBJ) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-ld: $(OBJLD) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJLD) -o $@ $(LIB) $(EXTRA) libBigWig.a -DNDEBUG -flto
dashing2-f: $(OBJF) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJF) -o $@ $(LIB) $(EXTRA) libBigWig.a -DNDEBUG -flto
read%: test/read%.o $(LIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(LIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
read%-ld: test/read%.ldo $(LDLIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(LDLIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a -DDSKETCH_FLOAT_TYPE="long double"
read%-f: test/read%.fo $(FLIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(FLIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a -DSKETCH_FLOAT_TYPE="float"
%: test/%.cpp $(LIBOBJ)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< $(LIBOBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a -DSKETCH_FLOAT_TYPE="float"
	# $(wildcard src/*.h)
%.o: %.cpp
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DNDEBUG -flto
%.0: %.cpp
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DNDEBUG -O1
%.lto: %.cpp
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DNDEBUG -flto
%.do: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -O0
%.sano: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -flto -fsanitize=address
%.go: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -pg
%.ldo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="long double" -DNDEBUG -flto
%.fo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="float" -DNDEBUG



libBigWig.a: $(wildcard libBigWig/*.c) $(wildcard libBigWig/*.h)
	cd libBigWig && sed -i $(SEDSTR) 's/HAVE_CURL:/#/' Makefile && $(MAKE) && cp libBigWig.a ..

test: readfx readbw

clean:
	rm -f dashing2 dashing2-ld dashing2-f libBigWig.a $(OBJ) $(OBJLD) $(OBJF) readfx readfx-f readfx-ld readbw readbw readbw-f readbw-ld src/*.0 src/*.do src/*.fo src/*.go src/*.ldo
