#include"adc_pv.h"
#include<iostream>
#include<iomanip>
#include <chrono>
#include <thread>
#include<exception>

using namespace std;

void pv_main()
{

	this_thread::sleep_for(chrono::seconds(1));
	

	//start conversion right now
	ADC_Configuration(ADC_CR, 2);

	//set resolution to 12 bits
	ADC_Configuration(ADC_MR, 0);

	//set gain to x1
	ADC_Configuration(ADC_CGR, 1);

	this_thread::sleep_for(chrono::seconds(2));

	//switch to second channel
	ADC_Configuration(ADC_CHER, 1);

	this_thread::sleep_for(chrono::seconds(2));

	//switch to third channel
	ADC_Configuration(ADC_CHER, 2);


	this_thread::yield();
}