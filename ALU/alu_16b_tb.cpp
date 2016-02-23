#include"alu_16b_tb.h"
#include<iostream>

using namespace std;

void alu_16b_tb::clk_process()
{
	while (true)
	{
		clk_val = !clk_val;
		clk.write(static_cast<sc_logic>(clk_val));
		wait(1, SC_NS);
	}
}

void alu_16b_tb::tb_process()
{
	rst.write(SC_LOGIC_0);		//keep reset for 10 ns
	wait(10, SC_NS);
	rst.write(SC_LOGIC_1);

	data1 = 0x00ff;
	data2= 0xff00;
	op_sel = 0;		// 0x00ff & 0xff00
	


	wait(2, SC_NS);
	op_sel = 1;		// 0x00ff | 0xff00


	wait(2, SC_NS);
	op_sel = 2;		// 0x00ff ^ 0xff00


	wait(2, SC_NS);	 // 0x00ff ror 1
	op_sel = 3;


	wait(2, SC_NS);  // 0x00ff rol 1
	op_sel = 4;


	wait(2, SC_NS); // 0x00ff + 0xff00 
	op_sel = 5;


	wait(2, SC_NS); // 8*4
	data1 = 0x0008;
	data2 = 0x0004;
	op_sel = 6;

	wait(2, SC_NS); // 8*-4 
	data1 = 0x0008;
	data2 = 0xFFFC;

	wait(2, SC_NS); // -8*-4
	data1 = 0xFFF8;
	data2 = 0xFFFC;


	wait(2, SC_NS);	// -8 / -4
	op_sel = 7;

	wait(2, SC_NS); // -8/ 0
	data1 = 0xFFF8;
	data2 = 0x0000;

	wait(2, SC_NS);

	rst.write(SC_LOGIC_0);

}