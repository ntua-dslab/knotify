CXX = g++
CXXFLAGS = -g -I./include -I../simfold/include -I../simfold/src/common -I../simfold/src/simfold -Wno-deprecated -Wno-write-strings -fPIC
LFLAGS = -lm -fPIC

SOURCES = $(wildcard src/*.cpp ../simfold/src/common/*.cpp ../simfold/src/simfold/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = libpkenergy.so

$(TARGET): $(OBJECTS)
	$(CXX) -shared -o $@ $(OBJECTS) $(LFLAGS)

clean:
	rm -rf $(TARGET) $(OBJECTS)
