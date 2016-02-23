#ifndef H_COUNTER_13B_TB
#define H_COUNTER_13B_TB

#include<systemc.h>

SC_MODULE(counter_13b_tb)
{
  sc_out<sc_logic> clk,rst,cnt_en,ud_ctrl;
  sc_in<sc_logic> ovf_intr,unf_intr;
  sc_in<sc_uint<13>> cnt_out;
  sc_signal<bool> clk_signal;

  SC_CTOR(counter_13b_tb)
  {
    SC_THREAD(tb_process);
    SC_THREAD(clk_process);
	clk_signal = false;
  }

  void tb_process(void);
  void clk_process(void);
  
};

#endif
