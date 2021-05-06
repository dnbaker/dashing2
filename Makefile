CXX?=g++

LIB=-lz
INCLUDE+=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -Ibonsai
OPT+=-std=c++17 -O3 -march=native
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts
EXTRA+=-DNOCURL

%: src/%.cpp
	$(CXX) $(INCLUDE) $(OPT) $(WARNING) $(MACH) $< -o $@ $(LIB) $(EXTRA) libBigWig.a

libBigWig.a: libBigWig/libBigWig.a
	cd libBigWig && $(MAKE) && cp libBigWig.a ..
