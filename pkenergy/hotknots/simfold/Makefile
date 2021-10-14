CXX = g++
RANLIB = ranlib
CXXFLAGS = -g -I./include -I./src/common -I./src/simfold -Wno-deprecated -fPIC -Wno-write-strings
LFLAGS = -lm

SOURCES = $(wildcard src/common/*.cpp) $(wildcard src/simfold/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = libsimfold.a

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@

clean:
	rm -rf $(TARGET) $(OBJECTS)
