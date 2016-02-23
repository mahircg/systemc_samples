#ifndef H_CONSUMER
#define H_CONSUMER
#include"fifo.h"

SC_MODULE(consumer)
{
	sc_port < reader_interface<sc_uint<5>>> in;
	sc_uint<5> in_data;

	SC_CTOR(consumer)
	{
		SC_THREAD(read_data);
	}

	void read_data()
	{
		while (true)
		{
			in->read_b(in_data);
			cout <<this->name()<<"---"<<"Data:"<< in_data << endl;
			cout << this->name() << "---" << "Number of elements:" << (unsigned)(in->get_ndata()) << endl;
		}
	}
};


#endif