#ifndef H_PRODUCER
#define H_PRODUCER
#include"fifo.h"

const unsigned short DATA_OUT_SIZE=31;

SC_MODULE(producer)
{
	
	sc_port<writer_interface<sc_uint<5>>> out;
	sc_uint<5> out_data_1[DATA_OUT_SIZE + 1], out_data_2[DATA_OUT_SIZE + 1];
	
	

	SC_CTOR(producer) 
	{
		for (int i = 0; i <= DATA_OUT_SIZE; i++)
		{
			out_data_1[i] = i;
			out_data_2[i] = DATA_OUT_SIZE-i;
		}
			
		SC_THREAD(write_data_proc_1);
		SC_THREAD(write_data_proc_2);
	}

	void write_data_proc_1()
	{
		unsigned short index = 0;
		do
		{
			
			out->write_b(out_data_1[index]);
			cout << this->name() <<" (process 1) "<< "---" << "writing: " << out_data_1[index] << endl;
			index += 1;
			wait(1, SC_SEC);
		} while (index <= DATA_OUT_SIZE);

	}

	void write_data_proc_2()
	{
		unsigned short index = 0;
		do
		{
			/*From SystemC Documentation,page 100 : */
			/*---Process instances sensitive to the event will not be resumed or triggered until the process that called notify has suspended or returned.---*/
			/*So,if producer and consumer start write_nb and read_nb,respectively, at the same time, consumer has to wait until FIFO gets full(or writer exits)*/
			/*This explains the behaviour when producer writes data without any delay*/
			out->write_b(out_data_2[index]);
			cout << this->name() << " (process 2) " << "---" << "writing: " << out_data_2[index] << endl;
			index += 1;
		} while (index <= DATA_OUT_SIZE);

	}



	
};

#endif

