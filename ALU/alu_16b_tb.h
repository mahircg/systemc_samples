#ifndef H_ALU_16B_TB
#define H_ALU_16B_TB

#include<systemc.h>



SC_MODULE(alu_16b_tb)
{

	sc_out<sc_logic> clk, rst;
	sc_out<sc_uint<3>> op_sel;
	sc_out<sc_int<16>> data1, data2;
	sc_in<sc_int<16>> acc;
	sc_in<sc_bv<8>> stat;
	
	bool clk_val;
	void clk_process();
	void tb_process();

	SC_CTOR(alu_16b_tb)
	{
		clk_val = false;
		SC_THREAD(clk_process);
		SC_THREAD(tb_process);
	}

};


#endif