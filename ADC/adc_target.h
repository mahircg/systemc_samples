#ifndef H_TARGET
#define H_TARGET

#include <tlm.h>
#include <systemc.h>
#include <tlm_utils/simple_target_socket.h>
#include "adc_pv.h"
#include "adc.h"

using namespace tlm;
using namespace tlm_utils;

SC_MODULE(adc_target)
{
private:
	int buffer;
public:
	simple_target_socket<adc_target> socket;
	sc_in<sc_bv<32>> reg_data_in;
	sc_out<sc_bv<32>> reg_data_out;
	sc_out<sc_uint<5>> reg_addr;
	sc_out<bool> reg_ctrl;
	
	SC_CTOR(adc_target) : socket("target_socket")
	{

		socket.register_b_transport(this, &adc_target::b_transport);
	}

	virtual void b_transport(tlm_generic_payload& trans, sc_time& delay)
	{
		tlm_command command = trans.get_command();
		uint64 adr = trans.get_address();
		unsigned int len = trans.get_data_length();
		unsigned char* ptr = trans.get_data_ptr();
		if (adr<ADC_BASE || adr > ADC_WPMR)
		{
			SC_REPORT_ERROR("TLM-2", "Address out of range");
			trans.set_response_status(TLM_ADDRESS_ERROR_RESPONSE);
		}
		else if (len > ADC_PAYLOAD_LENGTH || len == 0)
		{
			
			SC_REPORT_ERROR("TLM2", "Invalid payload length");
			trans.set_response_status(TLM_BURST_ERROR_RESPONSE);
		}
		else
		{
#pragma region HANDLE WRITE
			if (command == TLM_WRITE_COMMAND)
			{
				SC_REPORT_INFO("TLM2","Write command received");
				memcpy(&buffer, ptr, len);
				switch (adr)
				{
				case ADC_CR:
					reg_ctrl.write(true);
					reg_addr.write(0);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_CR write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_MR:
					reg_ctrl.write(true);
					reg_addr.write(1);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_MR write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_CHER:
					reg_ctrl.write(true);
					reg_addr.write(2);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_CHER write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_CHDR:
					reg_ctrl.write(true);
					reg_addr.write(3);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_CHDR write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_WPMR:
					reg_ctrl.write(true);
					reg_addr.write(23);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_WPMR write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_COR:
					reg_ctrl.write(true);
					reg_addr.write(6);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_COR write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_CGR:
					reg_ctrl.write(true);
					reg_addr.write(5);
					reg_data_out.write(buffer);
					SC_REPORT_INFO("TLM2", "ADC_CGR write completed");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				default:
					SC_REPORT_ERROR("TLM2", "Cannot write read-only register");
					trans.set_response_status(TLM_GENERIC_ERROR_RESPONSE);
					break;
					
				}
			}
#pragma endregion
#pragma region HANDLE READ
			else if (command == TLM_READ_COMMAND)
			{
				SC_REPORT_INFO("TLM2", "Read command received");
				switch (adr)
				{
				case ADC_MR:
					reg_addr.write(1);
					reg_ctrl.write(false);
					buffer=reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_MR");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_CHSR:
					reg_addr.write(4);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CHSR");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_CGR:
					reg_addr.write(5);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CGR");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_COR:
					reg_addr.write(6);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_COR");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_CDR0:
					reg_addr.write(7);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR0");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR1:
					reg_addr.write(8);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR1");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR2:
					reg_addr.write(9);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR2");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR3:
					reg_addr.write(10);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR3");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR4:
					reg_addr.write(11);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR4");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR5:
					reg_addr.write(12);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR5");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR6:
					reg_addr.write(13);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR6");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR7:
					reg_addr.write(14);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR7");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR8:
					reg_addr.write(15);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR8");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR9:
					reg_addr.write(16);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR9");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR10:
					reg_addr.write(17);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR10");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR11:
					reg_addr.write(18);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR11");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR12:
					reg_addr.write(19);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR12");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR13:
					reg_addr.write(20);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR13");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR14:
					reg_addr.write(21);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR14");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;

				case ADC_CDR15:
					reg_addr.write(22);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_CDR15");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				case ADC_WPMR:
					reg_addr.write(23);
					reg_ctrl.write(false);
					buffer = reg_data_in.read().to_int();
					memcpy(ptr, &buffer, len);
					SC_REPORT_INFO("TLM2", "Target sends ADC_WPMR");
					trans.set_response_status(TLM_OK_RESPONSE);
					break;
				default:
					SC_REPORT_ERROR("TLM2", "Cannot read write-only register");
					trans.set_response_status(TLM_GENERIC_ERROR_RESPONSE);
					break;
				}
			}
#pragma endregion
		}
	}

	
};


#endif
