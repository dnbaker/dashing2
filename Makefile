CXX?=g++

LIB=-lz
INCLUDE+=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -Ibonsai
OPT+=-std=c++17 -O3 -march=native
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts
EXTRA+=-DNOCURL

OBJ=src/d2.o src/bedsketch.o
OBJLD=src/d2.ldo src/bedsketch.ldo
OBJF=src/d2.fo src/bedsketch.fo
dashing2: $(OBJ) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-ld: $(OBJLD) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJLD) -o $@ $(LIB) $(EXTRA) libBigWig.a
dashing2-f32: $(OBJF) libBigWig.a
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $(OBJ) -o $@ $(LIB) $(EXTRA) libBigWig.a
%.o: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA)
%.ldo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="long double"
%.fo: %.cpp $(wildcard src/*.h)
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(LIB) $(EXTRA) -DSKETCH_FLOAT_TYPE="float"

libBigWig.a: $(wildcard libBigWig/*.c) $(wildcard libBigWig/*.h)
	cd libBigWig && $(MAKE) && cp libBigWig.a ..

all: dashing2 dashing2-ld dashing2-f32
