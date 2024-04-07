.PHONY: clean all

# Default C++ & C compiler
CXX?=g++
CC?=gcc

# Cache size default value & flag
CACHE_SIZE?=4194304
CACHE_SIZE_FLAG:=-DD2_CACHE_SIZE=${CACHE_SIZE}
GIT_VERSION?=v2.1.19

# Target architecture flag for amd64 Linux
TARGET_FLAG=-m64

LIB=-lz
INC=-IlibBigWig -Ibonsai/include -Ibonsai -Ibonsai/hll -Ibonsai/hll/include -I. -Isrc -Ifmt/include
OPT+= -O3 \
    $(TARGET_FLAG) \
    -fopenmp -pipe $(CACHE_SIZE_FLAG)

OPTMV:=$(OPT)
CXXSTD?=-std=c++20
OPT+= $(CXXSTD)
WARNING+=-Wall -Wextra -Wno-unused-function -Wno-char-subscripts -pedantic -Wno-array-bounds
EXTRA+=-DNOCURL -DDASHING2_VERSION=\"$(GIT_VERSION)\" -DFMT_HEADER_ONLY

# Object files
OFS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) $(patsubst %.c,%.o, $(wildcard src/*.c))
OBJ=$(OFS) src/osfmt.o

# Rule to compile .cpp files to .o files
%.o: %.cpp
	$(CXX) $(INC) $(OPT) $(WARNING) $< -c -o $@ $(EXTRA) -DNDEBUG

# Rule to compile .c files to .o files
%.o: %.c
	$(CC) $(INC) $(OPTMV) $(WARNING) $< -c -o $@ $(EXTRA) -DNDEBUG -std=c11

src/osfmt.o: fmt/src/os.cc
	$(CXX) -I fmt/include $(OPT) $(WARNING) $< -c -o $@ $(EXTRA)

# Main target to build the static library
all: libdashing2.a

libdashing2.a: $(OBJ)
	ar rcs $@ $^

clean:
	rm -f libdashing2.a $(OBJ)
