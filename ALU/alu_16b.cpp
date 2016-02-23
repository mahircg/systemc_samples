#include"alu_16b.h"


void alu_16b::alu_process()
{

	if (rst == SC_LOGIC_1)
	{
		data1_reg = data1.read();
		data2_reg = data2.read();
		stat_reg = 0x00;
		switch (op_sel.read())
		{
		case 0:
			acc_reg = data1_reg & data2_reg;
			stat_reg[2] = (acc_reg < 0);
			stat_reg[0] = (acc_reg == 0);
			stat_reg[3] = (data1_reg > 0 && data2_reg > 0 && acc_reg < 0);
			stat_reg[4] = (data1_reg < 0 && data2_reg < 0 && acc_reg > 0);
			break;

		case 1:
			acc_reg = data1_reg | data2_reg;
			stat_reg[2] = (acc_reg < 0);
			stat_reg[0] = (acc_reg == 0);
			stat_reg[3] = (data1_reg > 0 && data2_reg > 0 && acc_reg < 0);
			stat_reg[4] = (data1_reg < 0 && data2_reg < 0 && acc_reg > 0);
			break;

		case 2:
			acc_reg = data1_reg ^ data2_reg;
			stat_reg[2] = (acc_reg < 0);
			stat_reg[0] = (acc_reg == 0);
			stat_reg[3] = (data1_reg > 0 && data2_reg > 0 && acc_reg < 0);
			stat_reg[4] = (data1_reg < 0 && data2_reg < 0 && acc_reg > 0);
			break;

		case 3:
			acc_reg = static_cast<sc_bv<16>>(data1_reg).rrotate(1);
			stat_reg[2] = (acc_reg < 0);
			stat_reg[0] = (acc_reg == 0);
			stat_reg[3] = (data1_reg > 0  && acc_reg < 0);
			stat_reg[4] = (data1_reg < 0  && acc_reg > 0);
			break;
		case 4:
			acc_reg = static_cast<sc_bv<16>>(data1_reg).lrotate(1);
			stat_reg[2] = (acc_reg < 0);
			stat_reg[0] = (acc_reg == 0);
			stat_reg[3] = (data1_reg > 0  && acc_reg < 0);
			stat_reg[4] = (data1_reg < 0  && acc_reg > 0);
			break;
		case 5:
			acc_reg = data1_reg + data2_reg;
			stat_reg[0] = (acc_reg == 0);		
			stat_reg[3] = (data1_reg > 0 && data2_reg > 0 && acc_reg < 0); //overflow		
			stat_reg[4] = (data1_reg < 0 && data2_reg < 0 && acc_reg > 0); //underlow
			stat_reg[2] = (acc_reg < 0);					//negative flag
			stat_reg[1] = (acc_reg > data1_reg);			//carry flag
			break;
		case 6:
			acc_reg = data1_reg * data2_reg;
			stat_reg[0] = (acc_reg == 0);
			
			if (data1_reg > 0 && data2_reg > 0)
				stat_reg[3] = (acc_reg < data1_reg | acc_reg< data2_reg);
			else if (data1_reg > 0 && data2_reg < 0)
				stat_reg[4] = (acc_reg > data2_reg);
			else if (data1_reg < 0 && data2_reg > 0)
				stat_reg[4] = (acc_reg > data1_reg);
			else if (data1_reg < 0 && data2_reg < 0)
				stat_reg[4] = (acc_reg < 0);
			else
				stat_reg[4] = stat_reg[3] = SC_LOGIC_0;

			stat_reg[2] = (acc_reg < 0);					//negative flag
			stat_reg[1] = false;	
			break;
		case 7:
			if (data2_reg != 0)
			{
				acc_reg = data1_reg / data2_reg;
				stat_reg[0] = (acc_reg == 0);
				stat_reg[2] = (acc_reg < 0);					//negative flag
			}
			else
			{
				acc_reg = 0;
				stat_reg[3] = true;							//set both underflow and overflow in case of divide-by-zero
				stat_reg[4] = true;
			}
			break;
		}
		acc.write(acc_reg);
		stat.write(stat_reg);
	}
	else
	{
		stat_reg = 0x0000;
		data1_reg = data1_reg = acc_reg = 0;
		acc.write(acc_reg);
		stat.write(stat_reg);
	}
}