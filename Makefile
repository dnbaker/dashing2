# Special target that is not associated with a file but with a command to be executed.
.PHONY=clean

#Default C++ & C compiler
CXX?=g++
CC?=gcc

#Cache size default value & flag
CACHE_SIZE?=4194304
CACHE_SIZE_FLAG:=-DD2_CACHE_SIZE=${CACHE_SIZE}
GIT_VERSION?=v2.1.19


# If on M1, use -target arm64-apple-macos11 -mmacosx-version-min=11.0
# Otherwise, use march=native
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),arm)
	TARGET_FLAG=-target arm64-apple-macos11 -mmacosx-version-min=11.0
else
	TARGET_FLAG=-march=native
endif

LIB=-lz # -lfmt, library to link against
INC=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -Ibonsai -I. -Isrc -Ifmt/include #includes
OPT+= -O3 \ # additions to the compiler flags
	$(TARGET_FLAG) \ 
	-fopenmp -pipe $(CACHE_SIZE_FLAG)

OPTMV:=$(OPT)
CXXSTD?=-std=c++20
OPT+= $(CXXSTD)
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts -pedantic -Wno-array-bounds # -Wno-shift-count-overflow
EXTRA+=-DNOCURL -DDASHING2_VERSION=\"$(GIT_VERSION)\" -DFMT_HEADER_ONLY
CXXFLAGS+= $(CXXSTD)
CFLAGS+= -std=c11

D2SRC=$(wildcard src/*.cpp)
OFS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) $(patsubst %.c,%.o, $(wildcard src/*.c))
OBJ=$(OFS) src/osfmt.o
OBJLD=$(patsubst %.o,%.ldo,$(OFS)) src/osfmt.o
OBJF=$(patsubst %.o,%.fo,$(OFS)) src/osfmt.o
OBJF64=$(patsubst %.o,%.f64o,$(OFS)) src/osfmt.o
OBJLD64=$(patsubst %.o,%.ld64o,$(OFS)) src/osfmt.o
OBJ64=$(patsubst %.o,%.64o,$(OFS)) src/osfmt.o
OBJDBG=$(patsubst %.o,%.do,$(OFS)) src/osfmt.o
OBJADD=$(patsubst %.o,%.sano,$(OFS)) src/osfmt.o
OBJLTO=$(patsubst %.o,%.lto,$(OFS)) src/osfmt.o
OBJ0=$(patsubst %.o,%.0,$(OFS)) src/osfmt.o
OBJV=$(patsubst %.o,%.vo,$(OFS)) src/osfmt.o
OBJG=$(patsubst %.o,%.gobj,$(OFS)) src/osfmt.o
OBJW=$(patsubst %.o,%.wo,$(OFS)) src/osfmt.o
OBJNLTO=$(patsubst %.o,%.nlto,$(OFS)) src/osfmt.o

# Preprocessing commands for AVX2 and AVX512BW specific object files
D2SRCSTATICAVX2=$(patsubst %.cpp,%.static-avx2.o,$(D2SRC))
D2SRCSTATICAVX512BW=$(patsubst %.cpp,%.static-avx512bw.o,$(D2SRC))

libdashing2: $(OBJ)
	ar rcs $@ $^

%.o: %.cpp # compiles  .cpp files to .o files
	$(CXX) $(INC) $(OPT) $(WARNING) $(MACH) $< -c -o $@ $(EXTRA) -DNDEBUG -O3
%.o: %.c # compiles .c files to .o files
	$(CC) $(INC) $(OPTMV) $(WARNING) $(MACH) $< -c -o $@ $(EXTRA) -DNDEBUG -O3 -std=c11
src/osfmt.o: fmt/src/os.cc
	$(CXX) -I fmt/include $(OPT) $(WARNING) $< -c -o $@ $(EXTRA)

clean:
	rm -f libdashing2 $(OBJ)
