#ifndef H_COUNTER_13B_H
#define H_COUNTER_13B_H

#include <systemc.h>

#define COUNTER_MAX 8192

SC_MODULE(counter_13b)
{
  sc_in<sc_logic> clk,rst,cnt_en,ud_ctrl;
  sc_out<sc_logic> ovf_intr,unf_intr;
  sc_out<sc_uint<13>> cnt_out;
  sc_uint<13> counter;

 

 SC_CTOR(counter_13b):counter(0)
  {
    SC_METHOD(count_process);
    sensitive<<clk.pos();
  }

  void count_process();
  
};

#endif
