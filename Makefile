CXX=g++
CXXFLAGS=-fopenmp -O3 -I./include -I./../Libraries/eigen
JSONFLAGS=-I./../Libraries/json/include
BUILD_DIR=build
LIB_DIR=$(BUILD_DIR)/lib
BIN_DIR=$(BUILD_DIR)/bin

.PHONY: all clean

all: $(BIN_DIR)/_cyclicPerturbation.out $(BIN_DIR)/_spectralFunction.out $(BIN_DIR)/_classicalSpectralFunction.out

$(LIB_DIR):
mkdir -p $(LIB_DIR)

$(BIN_DIR):
mkdir -p $(BIN_DIR)

$(LIB_DIR)/2spin.o: src/2spin.cpp | $(LIB_DIR)
$(CXX) -c $< $(CXXFLAGS) -o $@

$(LIB_DIR)/cyclicPerturbation.o: src/cyclicPerturbation.cpp | $(LIB_DIR)
$(CXX) -c $< $(CXXFLAGS) $(JSONFLAGS) -o $@

$(BIN_DIR)/_cyclicPerturbation.out: $(LIB_DIR)/cyclicPerturbation.o $(LIB_DIR)/2spin.o | $(BIN_DIR)
$(CXX) -o $@ -fopenmp -O3 $^

$(LIB_DIR)/spectralFunction.o: src/spectralFunction.cpp | $(LIB_DIR)
$(CXX) -c $< $(CXXFLAGS) $(JSONFLAGS) -o $@

$(BIN_DIR)/_spectralFunction.out: $(LIB_DIR)/spectralFunction.o $(LIB_DIR)/2spin.o | $(BIN_DIR)
$(CXX) -o $@ -fopenmp -O3 $^

$(LIB_DIR)/classical2spin.o: src/classical2spin.cpp | $(LIB_DIR)
$(CXX) -c $< $(CXXFLAGS) -o $@

$(LIB_DIR)/classicalSpectralFunction.o: src/classicalSpectralFunction.cpp | $(LIB_DIR)
$(CXX) -c $< $(CXXFLAGS) $(JSONFLAGS) -o $@

$(BIN_DIR)/_classicalSpectralFunction.out: $(LIB_DIR)/classicalSpectralFunction.o $(LIB_DIR)/classical2spin.o | $(BIN_DIR)
$(CXX) -o $@ -fopenmp -O3 $^

clean:
rm -rf $(BUILD_DIR)
