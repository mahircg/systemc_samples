CXX = 	 g++
HEADER = -I/cygdrive/c/Users/Mahircan/Desktop/Foo/systemc-2.3.1/include
LIB = 	 -L/cygdrive/c/Users/Mahircan/Desktop/Foo/systemc-2.3.1/lib-cygwin
LIBRARY = -lsystemc
SRC =  main.cpp $(TARGET).cpp $(TARGET).h $(TARGET_TB).cpp $(TARGET_TB).h
TARGET = alu_16b
TARGET_TB = alu_16b_tb
CXX_FLAGS = -fpermissive -std=c++0x

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXX_FLAGS) $(HEADER) $(LIB) -o $(TARGET) $(SRC) $(LIBRARY)

clean: 
	$(RM) $(TARGET) $(TARGET)_trace.vcd
