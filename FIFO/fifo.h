#ifndef H_FIFO
#define H_FIFO

#include<systemc.h>


template<class Type>
class writer_interface: virtual public sc_interface
{
public:
	virtual void write_nb(Type) = 0;
	virtual void write_b(Type)=0;
	virtual void reset() = 0;
};

template<class Type>
class reader_interface : virtual public sc_interface
{
public:
	virtual void read_nb(Type&) = 0;
	virtual void read_b(Type&) = 0;
	virtual unsigned short get_ndata() = 0;
};


template<class Type, unsigned short size>
class generic_fifo: public sc_channel,public writer_interface<Type>,public reader_interface<Type>
{
private:
	sc_event read_event, write_event;
	unsigned short front,back,num_elements;
	Type fifo[size];
	sc_semaphore *fifo_lock;
	
public:
	generic_fifo(sc_module_name name) : sc_channel(name), front(0), back(0), num_elements(0)
	{
		fifo_lock = new sc_semaphore(1);
	}

	~generic_fifo()
	{
		delete fifo_lock;
	}

	void write_nb(Type);
	void write_b(Type);
	void read_nb(Type&);
	void read_b(Type&);
	void reset();
	unsigned short get_ndata();

};

template<class Type, unsigned short size>
void generic_fifo<Type, size>::write_nb(Type data)
{

	if (num_elements == size)
	{
		if (fifo_lock->trywait() == -1)
			fifo_lock->post();
		wait(read_event);
		
	}
	
	if (back == size)
		back = 0;
	fifo[back] = data;
	back += 1;
	num_elements += 1;
	write_event.notify();
	

}

template<class Type, unsigned short size>
void generic_fifo<Type, size>::write_b(Type data)
{
	fifo_lock->wait();
	write_nb(data);
	if (fifo_lock->get_value()==0)
		fifo_lock->post();
}

template<class Type, unsigned short size>
void generic_fifo<Type, size>::read_nb(Type& data)
{
	
	if (front == back)
	{
		
		if (fifo_lock->trywait() == -1)
			fifo_lock->post();
		wait(write_event);
		
		
	}
	if (front == size)
		front = 0;
	data = fifo[front];
	front += 1;
	num_elements -= 1;
	read_event.notify();
	
	
}

template<class Type, unsigned short size>
void generic_fifo<Type, size>::read_b(Type& data)
{
	fifo_lock->wait();
	read_nb(data);
	if (fifo_lock->get_value()==0)
		fifo_lock->post();

}

template<class Type, unsigned short size>
void generic_fifo<Type, size>::reset()
{
	cout << "FIFO is reset"<<endl<<endl;
	front = back = 0;
	num_elements = 0;
}

template<class Type, unsigned short size>
unsigned short generic_fifo<Type, size>::get_ndata()
{
	return num_elements;
}

#endif H_FIFO