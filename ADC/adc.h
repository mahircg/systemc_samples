#ifndef H_ADC
#define H_ADC

#include<systemc-ams.h>


//Mostly copied from TU-Vien library,only slight modifications made
SCA_TDF_MODULE(adc) {


private:
	double gain;
	double offset;
	sca_core::sca_time timestep;
	sca_util::sca_trace_file* atf;
public:

	sca_tdf::sca_in<double> in;	         // Input port
	sca_tdf::sca_de::sca_out<sc_bv<12> > out;            // Output port
	sca_tdf::sca_de::sca_in<sc_int<32>> offset_in;
	sca_tdf::sca_de::sca_in<sc_int<32>> gain_in;

	//variablen


	double u_ref;		        // converter's reference voltage
	double lsb;			// least significant bit
	long max;			// maximal output value
	int nl_mode;			// mode selection for liniarity error

	sc_bv<12> bv_erg;		// 16- bit Vektor for conv. value
	int first_inl;		// Merker fuer ersten Aufruf von INL



	

	void processing()
	{
		double analog;
		double temp;
		long erg;

		analog = in.read();	                
		offset = offset_in.read();
		gain = gain_in.read();

		temp = analog + offset;
		temp = temp*gain;
		

		erg = roundValue(temp / lsb);	// calculate digital value and rounding
		erg = maximumValue(erg);    	// output limitation
		bv_erg = erg;	        		// save as bBitvektor

		out.write(bv_erg);          	// write result to output port


	}

	void initialize(){

		// objekt of class "noiseGenerator"
		
		// analog value according to the maximal digital output of adc
		max = (long)(pow(2.0, 12 - 1.0));
		// 1bit lost because of sign
		lsb = u_ref / (max);	// calculate lsb of the adc

		first_inl = 0;          // to note the first call of INL()


		
	}

	void set_attributes()
	{
		set_timestep(timestep);
	}

	long roundValue(double val){
		return (long)lround(val);
	}



	long maximumValue(long value){

		if (abs(value) > max - 1){

			if (value > 0)
				value = max - 1;	//by positiv results => pos. maximum - 1 (according to the maximum analog value Uref-Uls)


			else
				value = -max;	//by negativ results => neg. Maximum
			// also -Uref
		}
		return value;
	}

	

	// Konstruktor (siehe Dokumentation)
	adc(sc_core::sc_module_name n, sca_core::sca_time tstep = sca_core::sca_time(100, SC_NS)) :timestep(tstep), offset_in("offset"), gain_in("gain"), out("out")
	{

		offset = gain = 0;
		u_ref=2.0,
		atf = sca_util::sca_create_vcd_trace_file("adc_out");
		sca_util::sca_trace(atf, out, "ADC_Out");
	}

	~adc()
	{
		sca_util::sca_close_vcd_trace_file(atf);
	}
};

#endif