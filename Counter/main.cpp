#include "counter_13b.h"
#include "counter_13b_tb.h"
#include <systemc.h>



int sc_main(int argc,char* argv[])
{
  sc_set_time_resolution(1,SC_PS);

  sc_signal<sc_logic> clk,rst,cnt_en,ud_ctrl,ovf_intr,unf_intr;
  sc_signal<sc_uint<13>> cnt;
  
  counter_13b counter("Counter-13Bit");
  counter.clk(clk);
  counter.rst(rst);
  counter.cnt_en(cnt_en);
  counter.ud_ctrl(ud_ctrl);
  counter.ovf_intr(ovf_intr);
  counter.unf_intr(unf_intr);
  counter.cnt_out(cnt);
  
  counter_13b_tb counter_tb("Counter-Testbench");
  counter_tb.clk(clk);
  counter_tb.rst(rst);
  counter_tb.cnt_en(cnt_en);
  counter_tb.ud_ctrl(ud_ctrl);
  counter_tb.ovf_intr(ovf_intr);
  counter_tb.unf_intr(unf_intr);
  counter_tb.cnt_out(cnt);
  

  sc_trace_file *vcd_trace = sc_create_vcd_trace_file("counter_13b_trace");
  vcd_trace->set_time_unit(1,SC_NS);

  sc_trace(vcd_trace,clk,"CLK");
  sc_trace(vcd_trace,rst,"RST");
  sc_trace(vcd_trace,cnt_en,"CNT_EN");
  sc_trace(vcd_trace,ud_ctrl,"UD_CTRL");
  sc_trace(vcd_trace,ovf_intr,"OVERFLOW");
  sc_trace(vcd_trace,unf_intr,"UNDERLOW");
  sc_trace(vcd_trace,cnt,"CNT_OUT");

  sc_start(8196*4+100,SC_NS);
  sc_stop();

  sc_close_vcd_trace_file(vcd_trace);
  
  return(0);
}
