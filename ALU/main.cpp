#include "alu_16b.h"
#include "alu_16b_tb.h"


int sc_main(int argc, char* argv[])
{
	sc_signal<sc_logic> clk, rst;
	sc_signal<sc_int<16>> data1, data2, acc;
	sc_signal<sc_uint<3>> op_sel;
	sc_signal<sc_bv<8>> stat;

	alu_16b alu("ALU_16B");
	alu.clk(clk);
	alu.rst(rst);
	alu.data1(data1);
	alu.data2(data2);
	alu.stat(stat);
	alu.acc(acc);
	alu.op_sel(op_sel);

	alu_16b_tb alu_tb("ALU_16B_TB");
	alu_tb.clk(clk);
	alu_tb.rst(rst);
	alu_tb.data1(data1);
	alu_tb.data2(data2);
	alu_tb.stat(stat);
	alu_tb.acc(acc);
	alu_tb.op_sel(op_sel);

	sc_trace_file *vcd_trace = sc_create_vcd_trace_file("alu_16b_trace");
	vcd_trace->set_time_unit(1, SC_NS);

	sc_trace(vcd_trace, clk, "CLK");
	sc_trace(vcd_trace, rst, "RST");
	sc_trace(vcd_trace, data1, "DATA1_IN");
	sc_trace(vcd_trace, data2, "DATA2_IN");
	sc_trace(vcd_trace, acc, "ACC_OUT");
	sc_trace(vcd_trace, op_sel, "OP_SEL");
	sc_trace(vcd_trace, stat, "STAT_REG");

	sc_start(40, SC_NS);
	sc_stop();

	sc_close_vcd_trace_file(vcd_trace);

	return 0;

}