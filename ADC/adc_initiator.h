#ifndef H_ADC_INITIATOR
#define H_ADC_INITIATOR

#include "adc_pv.h"
#include <tlm.h>
#include <systemc.h>
#include <tlm_utils/simple_initiator_socket.h>

using namespace tlm;
using namespace tlm_utils;

SC_MODULE(adc_initiator)
{

	friend void configure(adc_initiator &i_adc,unsigned int addr, int value);
	friend unsigned int get_status(adc_initiator &i_adc);
	friend unsigned int get_data(adc_initiator &i_adc, unsigned short index);

private:
	int reg_data;
	sc_bv<32> reg_addr;
	bool reg_ctrl;
	int data;
	uint64 addr;
	bool rw;
	tlm_command cmd;
	tlm_generic_payload *payload;

public:


	tlm_utils::simple_initiator_socket<adc_initiator> socket;
	SC_CTOR(adc_initiator) : socket("TLM_Initiator")
	{
		
		reg_data = 0;
		reg_addr = 0;
		reg_ctrl = false;
	}


	void comm_process()
	{
		sc_time delay=sc_time(0,SC_NS);	
		payload = new tlm_generic_payload();
		addr = reg_addr.to_uint64();
	    rw = reg_ctrl;
		SC_REPORT_INFO("TLM2", "Command received on TLM initiator");
		if (rw == 1)
		{
			data = reg_data; // store data to be written
			cmd = TLM_WRITE_COMMAND;
		}
		else
			cmd = TLM_READ_COMMAND;

		payload->set_command(cmd);
		payload->set_address(addr);
		payload->set_data_ptr(reinterpret_cast<unsigned char*>(&data));
		payload->set_data_length(4);
		payload->set_streaming_width(4);	//no streaming
		payload->set_byte_enable_ptr(0); 
		payload->set_dmi_allowed(false); 
		payload->set_response_status(tlm::TLM_INCOMPLETE_RESPONSE); 

		socket->b_transport(*payload, delay);
		

		if (payload->is_response_error())
		{
			SC_REPORT_INFO("TLM2", "Error on response");
			cout << payload->get_response_string() << endl;
			reg_data = 0;
		}
		else if (rw == false)
		{
			memcpy(&data, payload->get_data_ptr(), payload->get_data_length());
			reg_data=data;
		}
		delete payload;
		payload = nullptr;
		

	}
};

#endif

void configure(adc_initiator &i_adc,unsigned int addr, int value)
{
	i_adc.reg_ctrl = true;
	i_adc.reg_addr = addr;
	i_adc.reg_data = value;
	i_adc.comm_process();
	cout << "Control register " << hex << addr << " set to " << hex << value << endl;
	return;
}

unsigned int get_status(adc_initiator &i_adc)
{
	i_adc.reg_ctrl = false;
	i_adc.reg_addr = ADC_CHSR;
	i_adc.comm_process();
	return i_adc.reg_data;
}

unsigned int get_data(adc_initiator &i_adc, unsigned short index)
{
	i_adc.reg_ctrl = false;
	i_adc.reg_addr = ADC_CDR0 + index * 4;
	i_adc.comm_process();
	return i_adc.reg_data;
}
