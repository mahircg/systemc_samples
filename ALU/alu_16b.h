#ifndef H_ALU_16B
#define H_ALU_16B
#include<systemc.h>

#define MIN_INT -32768
#define MAX_INT 32767

SC_MODULE(alu_16b)
{

private:
	sc_int<16> data1_reg, data2_reg, acc_reg;
	sc_bv<8> stat_reg;

public:
	sc_in<sc_logic> clk,rst;
	sc_in<sc_uint<3>> op_sel;
	sc_in<sc_int<16>> data1,data2;
	sc_out<sc_int<16>> acc;
	sc_out<sc_bv<8>> stat;



	SC_CTOR(alu_16b)
	{
		data1_reg = data1_reg = acc_reg = 0;
		stat_reg = false;
		SC_METHOD(alu_process);
		sensitive << clk.pos();
	}

	void alu_process();
	
};

#endif