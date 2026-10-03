# lbm-squirmer -- build with `make` (OpenMP enabled if the compiler supports it)
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
OMPFLAGS ?= -fopenmp
SRC      := src/main.cpp $(wildcard src/utils/*.cpp)
HDR      := $(wildcard src/utils/*.hpp)
TARGET   := bin/lbm_squirmer

all: $(TARGET)

$(TARGET): $(SRC) $(HDR)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(SRC) -o $@

test: $(TARGET)
	python3 -m pytest tests -v

test-all: $(TARGET)
	python3 -m pytest tests -v --runslow

clean:
	rm -rf bin tests/output

.PHONY: all test test-all clean
