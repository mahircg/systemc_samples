#include<iostream>
#include"counter_13b_tb.h"
#include<systemc.h>

using namespace std;


void counter_13b_tb::clk_process(void)
{
  while(true)
    {
	  clk_signal.write(!clk_signal.read());
	  clk.write(static_cast<sc_logic>(clk_signal.read()));
      wait(1,SC_NS);
    }
}

 void counter_13b_tb::tb_process(void)
{
  //reset for first 100 ns
  rst=SC_LOGIC_0;


  wait(10,SC_NS);

  //deassert reset,start counting up
  rst = SC_LOGIC_0;
  
  cnt_en = SC_LOGIC_1;
  ud_ctrl = SC_LOGIC_1;
  
  wait(8196*2,SC_NS);

  ud_ctrl = SC_LOGIC_0;

  wait(8196*2,SC_NS);

  
}
