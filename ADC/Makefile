CXX = 	g++ 
HEADER = -I$(SYSTEMC)include -I$(SYSTEMC_AMS)include
LIB = 	 -L$(SYSTEMC)lib-linux64 -L$(SYSTEMC_AMS)lib-linux64	 
LIBRARY = -lsystemc -lsystemc-ams
PV_SRC = pv_main.cpp
SRC =  sc_main.cpp
TUV_SRC = tuv_ams_library.cpp
TARGET = exercise2
CXX_FLAGS = -fpermissive -std=c++0x

all: $(TARGET)

$(TARGET): pv_main.o tuv_ams_library.o sc_main.o
	$(CXX) $(CXX_FLAGS) -Wl,-rpath $(LIB) -o $(TARGET)  pv_main.o sc_main.o pv_main.o tuv_ams_library.o $(LIBRARY) 

sc_main.o: $(SRC)
	$(CXX) $(CXX_FLAGS) $(HEADER) -c $(SRC) 

tuv_ams_library.o: $(TUV_SRC)
	$(CXX) $(CXX_FLAGS) $(HEADER) -c $(TUV_SRC)

pv_main.o:$(PV_SRC)
	$(CXX) $(CXX_FLAGS) $(HEADER) -c $(PV_SRC)

clean: 
	$(RM) $(TARGET) *.o $(TARGET)_in.vcd $(TARGET)_out.vcd 
