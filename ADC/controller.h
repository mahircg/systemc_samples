#ifndef H_REGFILE
#define H_REGFILE

#include<systemc.h>
#include "adc_pv.h"


SC_MODULE(controller)
{
private:
	sc_bv<32> regfile[24];
	sc_event reg_updated;
	sc_uint<4> active_line;
	bool start = false;
	bool reset = false;
	int offset, gain;
	bool resolution = false;	//12-bits by default
public:
	// User interface signals,from TLM target/control logic
	sc_in<sc_uint<5>> reg_addr;
	sc_in<bool> reg_ctrl;
	sc_in<sc_bv<32>> reg_data_in;
	sc_out<sc_bv<32>> reg_data_out;
	sc_out<sc_int<32>> adc_gain;
	sc_out<sc_int<32>> adc_offset;
	sc_in<sc_bv<ADC_OUTPUT_PREC>> adc_data;
	sc_out<sc_uint<4>> in_select;

	SC_CTOR(controller)
	{
		
		for (auto& i : regfile)
			i = 0;
		active_line = 0;
		SC_METHOD(reg_process);
		sensitive << reg_ctrl << reg_addr << reg_data_in;
		SC_METHOD(read_process);
		sensitive << adc_data;
		SC_THREAD(control_process);
	}

	void read_process()
	{
		regfile[active_line + 7] = adc_data.read();
	}

	void control_process()
	{
		while (true)
		{
			wait(reg_updated);
			reset = regfile[0][0];
			if (reset == true)
			{
				for (auto& i : regfile)
					i = 0;
				active_line = 0;
			}
			else
			{
				
				start = regfile[0][1];
				resolution = regfile[1][4];
				active_line = regfile[2].range(15, 0).to_uint();
				offset = regfile[6].to_int();
				gain = regfile[5].to_int();
				adc_gain.write(gain);
				adc_offset.write(offset);
				in_select.write(active_line);
			}
		}
	}

	void reg_process()
	{
		bool rw = reg_ctrl.read();
		unsigned int data;
		unsigned short addr = reg_addr.read();

		if (rw == true)
		{
			data = reg_data_in.read().to_uint();
			regfile[addr] = data;
			reg_updated.notify();
		}
		else
		{
			if (start == true)
				data = regfile[addr].to_uint();
			else
				data = 0;
			reg_data_out.write(data);
		}

	}
};

#endif