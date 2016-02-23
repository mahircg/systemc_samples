CXX = 	 g++
HEADER = -I/cygdrive/c/Users/Mahircan/Desktop/Foo/systemc-2.3.1/include
LIB = 	 -L/cygdrive/c/Users/Mahircan/Desktop/Foo/systemc-2.3.1/lib-cygwin
LIBRARY = -lsystemc
SRC =  fifo.h consumer.h producer.h main.cpp
TARGET = fifo
CXX_FLAGS = -fpermissive -std=c++0x

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXX_FLAGS) $(HEADER) $(LIB) -o $(TARGET) $(SRC) $(LIBRARY)

clean: 
	$(RM) $(TARGET)
