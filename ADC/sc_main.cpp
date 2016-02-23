#include"slave.h"
#include"adc_initiator.h"
#include"tuv_ams_library.h"
#include<thread>


adc_initiator i_initiator("ADC_Initiator");
bool init = false;


SCA_TDF_MODULE(freq_gen)
{
	sca_tdf::sca_out<double> freq;

	SCA_CTOR(freq_gen)
	{}

	void set_attributes()
	{
		set_timestep(sca_time(1, SC_NS));
	}

	void initialize()
	{
		freq.initialize(0);
	}

	void processing()
	{
		freq.write(1000000);
	}
};



SCA_TDF_MODULE(sine_to_cosine)
{
	sca_tdf::sca_in<double> in;
	sca_tdf::sca_out<double> out;

	SCA_CTOR(sine_to_cosine)
	{}

	void set_attributes()
	{
		out.set_rate(1);
	}

	void initialize()
	{
		
	}

	void processing()
	{
		
		out.write(1 - in.read());
	}

};

extern void pv_main();

void ADC_Configuration(unsigned int addr, int value)
{
	
	configure(i_initiator, addr, value);
}

unsigned int ADC_GetStatus()
{

	return get_status(i_initiator);

}

unsigned int ADC_GetData(unsigned short index)
{
	return get_data(i_initiator, index);
}


int sc_main(int argc, char** argv)
{
	
	sc_set_time_resolution(1, SC_NS);
	slave i_slave("ADC_Slave");
	freq_gen i_gen("Freq_gen");
	


	//Bind initiator and target sockets
	i_initiator.socket.bind(i_slave.i_target.socket);
	
	sca_tdf::sca_signal<double> const_freq;
	sca_tdf::sca_signal<double> sig_sin_out;
	
	//Initialize sine generator
	i_gen.freq(const_freq);
	TUV_ams_lib::bb::sine Sine("Sine",0, 1, 0.0, 0.0, false, true, 1);
	Sine.set_timestep(1, sc_core::SC_NS); 
	Sine.freq_con(&const_freq);
	Sine.out(sig_sin_out);
	sca_util::sca_trace_file* atf = sca_util::sca_create_vcd_trace_file("adc_in");
	

	sca_tdf::sca_signal<double> saw_signal;
	//Initialize saw-tooth generator
	TUV_ams_lib::bb::saw_gen Saw("Saw",1000000,1);
	Saw.out(saw_signal);
	Saw.set_timestep(1, SC_NS);

	

	sca_tdf::sca_signal<double> cos_signal;
	//Initialize cosine generator
	sine_to_cosine Cos("Cos");
	Cos.in(sig_sin_out);
	Cos.out(cos_signal);

	
	i_slave.in[0](sig_sin_out);
	i_slave.in[1](cos_signal);
	i_slave.in[2](saw_signal);

	sca_util::sca_trace(atf, sig_sin_out, "Sin");
	sca_util::sca_trace(atf, cos_signal, "Cos");
	sca_util::sca_trace(atf, saw_signal, "Saw");

	//call user main function
	
	thread user_main_thread(pv_main);
	
	
	sc_start(100, SC_US);

	//user_main_thread.join();

	sca_util::sca_close_vcd_trace_file(atf);

	return 0;


}
