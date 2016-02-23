#include "counter_13b.h"


void counter_13b::count_process(void)
{
  if(rst==SC_LOGIC_0)
    {
      if(cnt_en ==SC_LOGIC_1)
	{
	  ovf_intr = (sc_logic)(counter == COUNTER_MAX - 1  && ud_ctrl == SC_LOGIC_1);
	  unf_intr = (sc_logic)(counter == 0 && ud_ctrl == SC_LOGIC_0);
	  if(ud_ctrl ==SC_LOGIC_1)
	    counter+=1;
	  else
	    counter-=1;
	}
      cnt_out.write(counter);
    }
  else
    {
      counter =0;
      unf_intr =(SC_LOGIC_0);
      ovf_intr =(SC_LOGIC_0);
      cnt_out = 0;
    }
}
