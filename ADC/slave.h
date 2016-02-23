#ifndef H_SLAVE
#define H_SLAVE

#include"adc.h"
#include"adc_target.h"
#include"controller.h"
#include"conv_tdf_sc.h"

SCA_TDF_MODULE(analog_encoder)
{
private:
	sca_time timestep;
public:
	sc_vector<sca_tdf::sca_in<double>> in;
	sca_tdf::sca_de::sca_in<sc_uint<4>> in_select;
	sca_tdf::sca_out<double> out;

	analog_encoder(sc_core::sc_module_name n, sca_time tstep = sca_time(1, SC_NS)) : timestep(tstep)
	{
		in.init(3);
	}

	void set_attributes()
	{
		set_timestep(timestep);
	}
	void initialize()
	{
		
	}

	void processing()
	{
		unsigned index = in_select.read();
		out.write(in[index].read());
	}


};

SC_MODULE(slave)
{
private:
	sc_core::sc_signal<sc_bv<32>> reg_data_in;
	sc_core::sc_signal<sc_bv<32>> reg_data_out;
	sc_core::sc_signal<sc_uint<5>> reg_addr;
	sc_core::sc_signal<bool> reg_ctrl;
	sc_core::sc_signal<sc_uint<4>> in_select;
	sc_core::sc_signal<double> controller_gain;
	sc_core::sc_signal<double> controller_offset;
	sc_core::sc_signal<sc_bv<ADC_OUTPUT_PREC>> adc_out_de;
	sc_core::sc_signal<sc_int<32>> adc_gain;
	sc_core::sc_signal<sc_int<32>> adc_offset;
	sca_tdf::sca_signal<double> encoder_out;

	adc i_adc;
	controller i_controller;
	sca_tdf::sca_signal<double> adc_sig_buffer;
	analog_encoder i_encoder;
	

public:
	sc_vector<sca_tdf::sca_in<double>> in;
	adc_target i_target;


	SC_CTOR(slave) :i_controller("Controller"), i_target("ADC_Target"), i_adc("ADC",sca_time(1,SC_NS)), i_encoder("Encoder")
	{
		
		
		in.init(3);
		int j = 0;
		for (auto& i : in)
		{
			i_encoder.in[j](i);
			j += 1;
		}

		i_encoder.out(encoder_out);
		i_encoder.in_select(in_select);
		

		//ADC signal bindings

		i_adc.gain_in(adc_gain);
		i_adc.offset_in(adc_offset);
		i_adc.out(adc_out_de);
		i_adc.in(encoder_out);

		//TLM target bindings
		i_target.reg_addr(reg_addr);
		i_target.reg_ctrl(reg_ctrl);
		i_target.reg_data_in(reg_data_in);
		i_target.reg_data_out(reg_data_out);

		//Controller bindings
		i_controller.reg_addr(reg_addr);
		i_controller.reg_ctrl(reg_ctrl);
		i_controller.reg_data_in(reg_data_out);
		i_controller.reg_data_out(reg_data_in);
		i_controller.in_select(in_select);
		i_controller.adc_data(adc_out_de);
		i_controller.adc_gain(adc_gain);
		i_controller.adc_offset(adc_offset);
		
		

	}


};


#endif