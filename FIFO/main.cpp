#include"fifo.h"
#include"producer.h"
#include"consumer.h"

int sc_main(int argc, char* argv[])
{
	generic_fifo<sc_uint<5>, 16> *fifo = new generic_fifo<sc_uint<5>, 16>("FIFO");
	producer *prod = new producer("Producer");
	consumer *cons = new consumer("Consumer");
	prod->out(*fifo);
	cons->in(*fifo);
	sc_start();

	return 0;
}