/*******************************************************************************************************
Copyright (c) 2008, Institute of Computer Technology (University of Technology Vienna)
Authors:
        Markus Damm <damm@ict.tuwien.ac.at>
        Jan Haase <haase@ict.tuwien.ac.at>
        Jiong Ou <ou@ict.tuwien.ac.at>
        Joseph Wenninger <wenninger@ict.tuwien.ac.at>
For the license please read LICENSE_BB.txt
**********************************************************************************************************/

#include <queue>
#include <systemc-ams.h>

#ifndef ANDRES_BB_LIB_H
#define ANDRES_BB_LIB_H 1


typedef double BB_DBL; //for visual c++, but perhaps it should be long double on windows?
#if defined(__APPLE__) && defined(__MACH__)
typedef unsigned int uint;
#endif

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <math.h>
#define log2(x) (log((BB_DBL)x)/log(2.0))
//acosh and asinh are from http://www.devx.com/vb2themax/Tip/19026
#define acosh(x) (log(x + sqrt(x * x - 1.0)))
#define asinh(x) (log(x + sqrt(x * x + 1.0)))
typedef unsigned int uint;
#endif

using namespace std;
using std::string;

namespace TUV_ams_lib {
  namespace bb {



    /*********************************Multiplier***********************************/

    SCA_TDF_MODULE(mul) {

      sca_tdf::sca_in<double>  in1;
      sca_tdf::sca_in<double>  in2;
      sca_tdf::sca_out<double>  out;

    private:
      void set_attributes() ;

      void initialize() ;

      void processing();


    public:


      mul(sc_core::sc_module_name nm);
    };


    /***************************************** Divider ****************************************************/
    // This module divides the value of input in1 from the value of input in2 and writes the result to the output port.
    SCA_TDF_MODULE(divi) {

      sca_tdf::sca_in<double>  in1;
      sca_tdf::sca_in<double>  in2;
      sca_tdf::sca_out<double>  out;

    private:
      void processing();

    public:
      divi(sc_core::sc_module_name nm);
    };


    /************************************Adder**********************************************************/

    SCA_TDF_MODULE(add) {

      sca_tdf::sca_in<double>  in1;
      sca_tdf::sca_in<double>  in2;
      sca_tdf::sca_out<double>  out;

    private:
      double gain ;
      int data_rate;

      void set_attributes() ;

      void initialize() ;

      void processing();


    public:
      add(sc_core::sc_module_name nm,int rate=1);
    };


    /***********************************Subtractor*********************************/
    SCA_TDF_MODULE(sub) {

      sca_tdf::sca_in<double>  in1;
      sca_tdf::sca_in<double>  in2;
      sca_tdf::sca_out<double>  out;

    private:
      void set_attributes() ;

      void initialize() ;

      void processing();


    public:
      sub(sc_core::sc_module_name nm);
                                        };

/************************************************* logarithm *******************************************/

SCA_TDF_MODULE(logarithm) {
 public:
  sca_tdf::sca_in<double> in;                        //input port
  sca_tdf::sca_out<double> out;                      //output port

 private:
 string base;

  void initialize();
  void processing();

 public:
  logarithm(sc_core::sc_module_name n, string _base);
                                                     };


/***********************************************tangnes******************************************/
SCA_TDF_MODULE(tang) {

 public:

  sca_tdf::sca_in<double> in;                        //input port
  sca_tdf::sca_out<double> out;                      //output port

 private:
  void processing();

 public:
  tang(sc_core::sc_module_name n);
                                  };

/******************************************** XNOR*********************************************/
SCA_TDF_MODULE(xnor) {

  sca_tdf::sca_in<bool> in1;
  sca_tdf::sca_in<bool> in2;
  sca_tdf::sca_out<bool> out;
 private:
    void processing();
 public:
    xnor(sc_core::sc_module_name n);
};

/****************************************Delay MODULE *********************************************/
// This class delays the input bit stream for "_delay" takts. The number of delay and the init value of output port can be set with parameter "_delay"
// and "_init", respectively.
SCA_TDF_MODULE(d_flipflop) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<bool> out;

 private:
  void set_attributes();
  void initialize();
  void processing();

  int delay;
  bool initial;
 public:
  // The constructor takes the inital value of the output port and the number of delay. They are set to true and 1 by default.
  d_flipflop(sc_core::sc_module_name n,bool _init=true,int _delay=1);
};

/********************************** Compare****************************** ***************************/
// This class compares the input of the module with a threshold value then output "true" if input is larger than the threshold, and vice verse.

SCA_TDF_MODULE(compare) {

  sca_tdf::sca_in<double>  in;
  sca_tdf::sca_out<bool>  out;

 private:
  double threshold;
  void processing();

 public:
  // The constructor takes the value of the threshold.It is 0.0 by default.
  compare(sc_module_name,double _threshold=0.0);
};




/******************************** NRZ **********************************/

SCA_TDF_MODULE(nrz) {

  sca_tdf::sca_in<bool>  in;
  sca_tdf::sca_out<double>  out;

 private:
  double vout;
  void processing();

 public:
  nrz(sc_module_name, double _vout=1.0);
};





    /***********************************Saturation*********************************/

    SCA_TDF_MODULE(saturation) {

      sca_tdf::sca_in<double>  in;
      sca_tdf::sca_out<double>  out;

    private:

      double OUT_MAX ;
      double OUT_MIN  ;


      void set_attributes() ;
      void initialize() ;
      void processing();

    public:
      saturation(sc_core::sc_module_name nm, double _OUT_MIN = 0.0 , double _OUT_MAX = 5.0);
                                                                                               };


    /**********************************Dead Zone**********************************/

    SCA_TDF_MODULE(deadzone) {

      sca_tdf::sca_in<double>  in;
      sca_tdf::sca_out<double>  out;

    private:

      double high;
      double low;


      double internal;

      void set_attributes() ;

      void initialize() ;

      void processing();


    public:
   deadzone(sc_module_name, double _high = 0.1, double _low = -0.1);

                                                                       };

    /********************************coulomb************************************/

    SCA_TDF_MODULE(coulomb) {

      sca_tdf::sca_in<double>  in;
      sca_tdf::sca_out<double>  out;

    private:

      double gain;
      double y;
      double internal;

      void set_attributes() ;

      void initialize() ;

      void processing();

    public:
      coulomb(sc_module_name, double _gain = 1, double _y = -0.1);

                                                                     };



/************************** INOUT**************************/

SCA_TDF_MODULE(inout) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<bool> out;
 private:
  void processing();
 public:
  inout(sc_core::sc_module_name n);
};




/*************************Offset**************************************/

SCA_TDF_MODULE(offset) {

  sca_tdf::sca_in<double>  in    ;
  sca_tdf::sca_out<double>  out ;

 private:
  double ofs ;

  void set_attributes() ;

  void initialize() ;

  void processing();

 public:
  offset(sc_module_name, double _offset = 1.0);

};

/**************************Integrator***************************/

SCA_TDF_MODULE(integ) {

  sca_tdf::sca_in<double> in;	// input port
  sca_tdf::sca_out<double> out;	// output port
 private:
  sca_tdf::sca_ltf_nd ltf_1;	        // The Laplace-Transform module
  double initial;	        // the initial value of output


  sca_vector<double> A,B;	// Vectors for the Laplace-Transform module

  void processing();

 public:
  integ(sc_core::sc_module_name n, double _initial=0.0);
};

/******************************Up Sampler**************************************/

SCA_TDF_MODULE(upsample) {

  sca_tdf::sca_in<double>  in     ;
  sca_tdf::sca_out<double>  out   ;

 private:
  int rate ;


  void set_attributes() ;

  void initialize() ;

  void processing();


 public:
  upsample(sc_module_name, int _rate = 1 );

};

/**************************Down Sampler***************************/
// This module decreases the data rate of the input signal. It reads _rate values and writes one value to the output. With the parameter _sel it is possible to selcet the _selth data which you want to output. It is an integer number between 1 and _rate.
SCA_TDF_MODULE(downsample) {

  sca_tdf::sca_in<double>  in     ;
  sca_tdf::sca_out<double>  out   ;

 private:
  int rate ;

  int sel;

  void set_attributes() ;

  void initialize() ;

  void processing();

 public:
  // The constructor takes the value of parameter _sel and the data rate of the input port.They are set to 1 and 1 by default, respectively.
  downsample(sc_core::sc_module_name nm,int _sel=1,int _rate=1);
};

/*********************************************************Fork with one output delayed one sample*********************************************************/

SCA_TDF_MODULE(o2t) {        //sub module for modulator which split one token to two tokens for the delay of half symbold uration //we have one symbol duration delay

  sca_tdf::sca_in<bool>  in;
  sca_tdf::sca_out<bool>  out;

 private:
  bool inp;
  void processing();
  void set_attributes();
 public:
  o2t(sc_core::sc_module_name n);
};
/************************************************Reverse Fork***************************************************************/
SCA_TDF_MODULE(t2o) {                                            //sub module to combine two tokens to one token

  sca_tdf::sca_in<bool>  in;
  sca_tdf::sca_out<bool>  out;

 private:
  bool inp;
  void processing();
  void set_attributes();
 public:
  t2o(sc_core::sc_module_name n);
};

/******************************************Nonlinearity***************************************************/
SCA_TDF_MODULE(nonlinearity) {

 public:
  sca_tdf::sca_in<double> in;
  sca_tdf::sca_out<double> out;

 private:
  double a, c, in_max, out_max;

  void processing();
 public:
  nonlinearity(sc_core::sc_module_name n, double gain, double ip3);
};



/***************************************** parallel to serial converter**********************************/

template <class T, int N>

  SCA_TDF_MODULE(p2s) {

  sca_tdf::sca_in<T> in[N];	                // input ports
  sca_tdf::sca_out<T> out;	                // output port

 private:
  int in_rate;
  int out_rate;
  T *data;

  void set_attributes()
  {
    for (int i=0;i<N;i++)
      {
	in[i].set_rate(in_rate);
      }
    out_rate=in_rate*N;
    out.set_rate(out_rate);
	data= new T[out_rate];
  }

  void processing()
  {

    int k=0;
    for (int i = 0; i < N; i++)
      {
	for(int j=0;j<in_rate;j++)
	  {
	    data[k]=in[i].read(j);
	    k++;
	  }
      }

    for(int j=0;j<out_rate;j++)
      {
	out.write(data[j],j);
      }
  }

 public:
  p2s(sc_core::sc_module_name n, int _in_rate=1)
    {
      in_rate=_in_rate;
    }
};


/***************************************serial to paralle converter*****************************************/

template <class T,int N>

  SCA_TDF_MODULE(s2p) {

  sca_tdf::sca_in<T> in;	        // input
  sca_tdf::sca_out<T> out[N];	// output

 private:
  int out_rate;
  int in_rate;
  T *symbol;

  void initialize(){};

  void set_attributes()
  {
    for (int i=0;i<N;i++)
      {
	out[i].set_rate(out_rate);
      }
    in_rate=N*out_rate;
    in.set_rate(in_rate);

    symbol= new T[N*out_rate];
  }

  void processing()
  {

    for (int i = 0; i < in_rate ; i++)
      {
	symbol[i]=in.read(i);
      }
    int k=0;
    for (int i=0;i<N;i++)
      {
	for(int j=0;j<out_rate;j++)
	  {
	    out[i].write(symbol[k],j);
	    k++;
	  }
      }
  }
 public:
  s2p(sc_core::sc_module_name n, int _out_rate=1)  {
    out_rate=_out_rate;
  }
};

/******************Uniformly distributed random Numbers****************/
    // This class generates uniformly distributed random numbers. The parameters _min and _max are used to set the upper and lower bound of the output values.
    // The data rate can be increased by the parameter _rate.
    // class definition : noise_uniform(sc_core::sc_module_name nm,double _min,double _max)
    // interfaces :       sca_tdf::sca_out<double> out

    SCA_TDF_MODULE(noise_uniform) {

      sca_tdf::sca_out<double>  out;                   // output port

    private:

      double min ;                                // variable to store value of lower bound of the output
      double max ;                                // variable to store value of upper bound of the output

      int rate;                                   // variable for data rate
      double rnd;                                 // variable to store calculated random numbers

      void set_attributes() ;
      void initialize() ;
      void processing() ;

    public:
      // The constructor takes the value of the upper and lower bound of the output value,
      // as well as the output data rate of the data token. They are by default 0.0, 1.0 and 1,respectively.
      noise_uniform(sc_core::sc_module_name, double _min = 0.0, double _max = 1.0, int _rate = 1) ;
    };
    /******************Uniformly distributed random Bits****************/
    // This class generates a uniformly distributed random sequence of bits on its output.
    // Class definition: rand_bool(sc_core::sc_module_name nm, int _rate)
    // Interfaces:       sca_tdf::sca_out<bool> out

    SCA_TDF_MODULE(rand_bool) {

      sca_tdf::sca_out<bool>  out;


    private:

      double rnd;                                // variable for random number
      int rate;                                  // variable for data rate

      void set_attributes();
      void initialize() ;
      void processing() ;

    public:
      // The constructor takes the value of the output data rate of the data token. It is 1 by default.
      rand_bool(sc_core::sc_module_name nm, int _rate=1) ;
    };


 /******************Gausssian distributed random numbers****************/
    // This class generates Gaussian distributed random numbers. The parameters _mean and _variance are used to
    // set the mean and variance of the random numbers. The data rate can be increased with the parameter _rate.
    // Class definition: noise_gauss(sc_core::sc_module_name nm, double _variance, double _mean);
    // Interfaces: 	 sca_tdf::sca_out<double> out ;

    SCA_TDF_MODULE(noise_gauss) {

      sca_tdf::sca_out<double>  out;

    private:

      double mean;
      double variance;
      int rate;

      void set_attributes();
      void initialize();
      void processing();

    public:
      // The constructor takes the value of the variance, mean of random number and data rate of the output data rate of the data token.
      // They are by default 1.0, 0.0 and 1,respectively.
      noise_gauss(sc_core::sc_module_name nm, double _variance = 1, double _mean = 0.0, int _rate = 1);
    };


    /************************** Sine Generator ***************************/
    // The output of this class is a sine wave. The frequency, amplitude, offset and the initial phase of the output
    // can be set with accordant parameters. You can also change the frequency and/or amplitude of the sine wave
    // during simulation using "freq_in" and "amp_in" ports. In this case parameter "freq_c" (frequency configuration)
    // and/or "amp_c" (amplitude configuration) should be set to "true". Otherwise "freq_in" and/or "amp_in" are/or not
    // available. Meanwhile you must add the "&" symbol before name of signals which are connected to these two ports.

    // Class definition: sine(sc_core::sc_module_name n, double freq_def, double amp_def, double _phi, double _offset, bool amp_c, bool freq_c,int datarate)
    // Interfaces: 	 sca_tdf::sca_out<double> out
    //                   sca_tdf::sca_in<double>* freq_in
    //                   sca_tdf::sca_in<double>* amp_in


    SCA_TDF_MODULE(sine) {
    public:
      sca_tdf::sca_out<double> out;                                            // output port
      sca_tdf::sca_in<double>* freq_in;                                        // input of frequency
      sca_tdf::sca_in<double>* amp_in;                                         // input of amplitude

      void freq_con (sca_tdf::sca_signal<double>* f_in);                       // function to refer value of signal "f_in" to the optional port "freq_in"
      void amp_con (sca_tdf::sca_signal<double>* a_in);                        // function to refer value of signal "a_in" to the optional port "ampl_in"

    private:
      bool freq_recon, ampl_recon;
      int rate;		                                                  // The datarate of the output
      double stepsize;	                                                  // the steps between the output token (i.e. if one token hase the
   			                                                  // value sin(x), then the next token has the value sin(x+stepsize)
      double actval;	                                                  // state variable to store the actual "x" value, if whe don't want
  			                                                  // to use sca_get_time()
      double sample_time;                                                 // variable to store sample time, which is the time between two tokens

      double freq_reg;
      double amp_reg;
      double offset;
      double phi;

      void set_attributes() ;
      void initialize() ;
      void processing();

    public:
      // The constructor takes the value of frequency, amplitute, initial phase and offset of the expected output signal.They are set to 1000,1.0,0.0,0.0 by default,
      // respectively. The reconfigurability of amplitude and frequency can be set per the parameter "amp_c" and "freq_c". Besides, one can set the data rate of output
      // port with "datarate".
      sine(sc_core::sc_module_name n, double freq_def=1000., double amp_def=1.0,double _phi=0.0,double _offset=0.0, bool amp_c=false, bool freq_c=false,int datarate=1);
    };


    /************************** Square Wave Generator ***************************/
    // The output of this class is a square wave. With the parameter "freq" and "amp" user can define the frequency and amplitude of the output signal, respectively.
    // The parameter "ofs" sets the offset of the output signal.With "d_rate" user can change the data rate of output port.
    SCA_TDF_MODULE(sqr_gen) {
      sca_tdf::sca_out<double> out;

    private:
      int rate;
      double stepsize;
      double f;
      double actval;
      double amp;
      double ofs;

      void set_attributes();
      void initialize();
      void processing();

    public:
      // The constructor takes the value of frequency, amplitute and offset of the expected output signal. Besides, user can set the data rate of output port with "d_rate".
      sqr_gen(sc_core::sc_module_name n, double freq, double _amp, double _ofs=0.0, int d_rate=1);
    };


 /************************** Triangle Wave Generator ***************************/
    // The output of this class is a triangle wave. With the parameter "freq" and "amp" user can define the frequency and amplitude of the output signal, respectively.
    // The parameter "ofs" sets the offset of the output signal.With "d_rate" user can change the data rate of output port.
    SCA_TDF_MODULE(tri_gen) {

      sca_tdf::sca_out<double> out;	// output is the filtered wave

    private:
      int rate;
      double stepsize;
      double actval;
      double amp;
      double ofs;
      double f;


      void set_attributes();
      void initialize();
      void processing();

    public:
      // The constructor takes the value of frequency, amplitute and offset of the expected output signal. Besides, user can set the data rate of output port with "d_rate".
      tri_gen(sc_core::sc_module_name n, double freq, double _amp, double _ofs=0.0, int d_rate=1);
    };

    /************************** Sawtooth Wave Generator ***************************/
    // The output of this class is a Sawtooth wave. With the parameter "freq" and "amp" user can define the frequency and amplitude of the output signal, respectively.
    // The parameter "ofs" sets the offset of the output signal.With "d_rate" user can change the data rate of output port.
    SCA_TDF_MODULE(saw_gen) {
      sca_tdf::sca_out<double> out;

    private:

      int rate;
      double f;
      double stepsize;
      double actval;
      double amp;
      double ofs;

      void set_attributes();
      void initialize();
      void processing();

    public:
      // The constructor takes the value of frequency, amplitute and offset of the expected output signal. Besides, user can set the data rate of output port with "d_rate".
      saw_gen(sc_core::sc_module_name n, double freq,double _amp, double _ofs=0.0, int d_rate=1);
    };






/********************************wave generator not supported************************************/



      SCA_TDF_MODULE(wave) {

	sca_tdf::sca_out<double> out;	// input double (wave)

      private:
	sqr_gen* sqr_sub;
	tri_gen* tri_sub;
	saw_gen* saw_sub;

	string w_type;

      public:
	void set_time_step(double sp,char unit);

	wave(sc_core::sc_module_name n, string type, double freq ,double amp, double ofs=0.0, int rate=1);
      };


/***************************FIRST Order Lowpass filter************************************/

SCA_TDF_MODULE(lp) {

  sca_tdf::sca_in<double> in;	// input double (wave)
  sca_tdf::sca_out<double> out;	// output is the filtered wave
 private:
  sca_tdf::sca_ltf_nd ltf_1;	// The Laplace-Transform module
  double freq_cutoff;	// the cutoff-frequency of the lowpass
  double out_tmp;

  void initialize();

  sca_vector<double> A,B;	// Vectors for the Laplace-Transform module

  void processing();

 public:
  lp(sc_core::sc_module_name n, double freq_cut);
};



#if  ( (!defined( _MSC_VER)) && (!defined(__APPLE__)) && (!defined(__MACH__)) )

/**********************************butterworth analog filter************************************************/

SCA_TDF_MODULE(btworth) {

	sca_tdf::sca_in<double> in;			// Input double (wave)
   	sca_tdf::sca_out<double> out;		// Output is the filtered wave

	sca_tdf::sca_ltf_nd ltf;

	string type;				//Filter type
	double gp;				//Passband gain
	double gs;				//Stopband gain
	double wp;				//Passband frequency
	double ws;				//Stopband frequency
	double wc;				//half-power frequency

	sca_vector<double> A;			//Coefficients of Butterworth polynomial
	sca_vector<double> B;			//
	sca_vector<double> poly_coef;		//Buffer for the Coefficients of Butterworth polynomial
	sca_vector<double> poly_dummy;		//Buffer for the Coefficients of Butterworth polynomial

	sca_vector<double> poles_unsorted;	//Unsorted Butterworth poles
	sca_vector<double> poles_sorted;

	void initialize();

	void processing();
	int order();
	void calc_poles(int n);
	void sort_poles(int n);
	void make_real(int n);
	void calc_char_poly(int n);
	int poly_multi(int m, int n);
	btworth(sc_core::sc_module_name nm, string _type, double _gp, double _gs, double _fp, double _fs);
 };

#endif


/***********************************chebyshev analog filter*************************************/

SCA_TDF_MODULE(chebyshev) {

 public:
  sca_tdf::sca_in<double> in;
  sca_tdf::sca_out<double> out;

  sca_tdf::sca_ltf_nd ltf_1;

 private:
  sca_vector<double> A,B;	// Vectors for the Laplace-Transform module
  string type;                  // Filter type:"lowpass", "highpass"
  int n;                        // Filter order
  double r,gs,wp,ws;            // Variables to define amplitude response of filter.wp:passband frequency in ; gp:minimum passband gain in dB by wp;ws: stopband frequency; gp:maximum stopband gain in dB by ws;
  double wp_tmp,ws_tmp;         // Tmp variables used by desciding wp and ws

  bool back_door;               // Variable for debuging

  void processing();

 public:
  chebyshev(sc_core::sc_module_name n, string _type="lowpass", double _ratio=2., double _gs=-20., double _wp=10, double _ws=20);
};

/*******************************************Gaussian puls shaping**********************************/
SCA_TDF_MODULE(gauss_shaping) {

  sca_tdf::sca_in<double> in;	// input double (wave)
  sca_tdf::sca_out<double> out;	// output is the filtered wave

 private:
  sca_tdf::sca_ltf_nd ltf_1;    	// The Laplace-Transform module
  double freq_cut;  	// the cutoff-frequency of the lowpass
  double bts;
  double ts;

  void initialize();
  sca_vector<double> A,B;
  void processing();

 public:
  gauss_shaping(sc_core::sc_module_name n, double _bts, double _ts);
};
/********************************Gain controller************************************/

SCA_TDF_MODULE(gain) {

  sca_tdf::sca_in<double>  in    ;
  sca_tdf::sca_out<double>  out ;
 private:
  double g ;
  double offset ;
  int rate;
  void set_attributes() ;

  void initialize() ;

  void processing();

 public:

  gain(sc_module_name, double _gain = 1.0, double _offset = 0.0, int data_rate=1);

};

/**************************************Mixer****************************************/

SCA_TDF_MODULE(mixer) {
 public:
  sca_tdf::sca_in<double> sig_in;
  sca_tdf::sca_in<double> lo_in;
  sca_tdf::sca_out<double> out;
 private:
  sca_tdf::sca_signal<double> sig_nonlinear;
  sca_tdf::sca_signal<double> sig_gout;
  bool ideal;
void processing();
  gain* i_gain;
  mul* ideal_mixer;
  nonlinearity* mixer_nonlinearity;

 public:
  mixer(sc_core::sc_module_name n, double _gain, double ip3, bool _ideal);
};

/**********************************phase detector/comparator***********************************/

SCA_TDF_MODULE(phc) {

 public:

  sca_tdf::sca_in<double> in_ref;        // input port
  sca_tdf::sca_in<double> in_vco;
  sca_tdf::sca_out<double> out;	// output port

 private:

  int rate;		// The datarate of the output
  double gain;          // Gain

  void set_attributes() ;


  void processing();

 public:

  phc(sc_core::sc_module_name n, int datarate, double _gain);
};



/*************************************peak detector***************************8************/
SCA_TDF_MODULE(peak_dt) {
 public:
  sca_tdf::sca_in<double> in;	                // input ports
  sca_tdf::sca_out<double> out;

 private:

  double peak;                                  //variable to store peak value
  inline double p_dt(double sig);

  void initialize();

  void processing();

 public:
  peak_dt(sc_core::sc_module_name n);
};

/****************************************Signal Splitter*******************************************/

template <class T, int N>

SCA_TDF_MODULE(splitter) {

  sca_tdf::sca_in<T> in;	                                 // input ports
  sca_tdf::sca_out<T> out[N];	                         // output port

  void processing()
  {
    for(int i=0; i<N; i++)  		                 // write <rate> data token for one symbol
      {
	out[i].write(in.read());
      }
  }

  splitter(sc_core::sc_module_name n){}
};


   /**************************** -45 grad shift modul***************************** *****/

SC_MODULE(nshift_45) {

 public:
   sca_tdf::sca_in<double> in;	// input double (wave)
   sca_tdf::sca_out<double> out;	// output is the filtered wave

 private:
   sca_eln::sca_tdf::sca_vsource* vin;	// SDF to voltage converter
   sca_eln::sca_tdf::sca_vsink* vout;	// and vice versa

   sca_eln::sca_c* c1;		// a capacitor
   sca_eln::sca_r* r1;		// a resistor

   sca_eln::sca_node n1,n2;	// two electrical nodes
   sca_eln::sca_node_ref gnd;	// and a reference


   double r;
   double c;

 public:
   nshift_45(sc_core::sc_module_name n, double w_if,double delta_R, double delta_C);
} ;


/********************************* +45 grad shift modul *************************************/

SC_MODULE(pshift_45) {

 public:
   sca_tdf::sca_in<double> in;	// input double (wave)
   sca_tdf::sca_out<double> out;	// output is the filtered wave

 private:
   sca_eln::sca_tdf::sca_vsource* vin;	// SDF to voltage converter
   sca_eln::sca_tdf::sca_vsink* vout;	// and vice versa

   sca_eln::sca_c* c1;		// a capacitor
   sca_eln::sca_r* r1;		// a resistor

   sca_eln::sca_node n1,n2;	// two electrical nodes
   sca_eln::sca_node_ref gnd;	// and a reference


   double r;
   double c;

 public:
   pshift_45(sc_core::sc_module_name n, double w_if,double delta_R, double delta_C);
} ;

/****** **********************************shifter ***********************************************/

SC_MODULE(shift_45) {
 public:
   sca_tdf::sca_in<double>   in;
   sca_tdf::sca_out<double>  out;

   // private:
   sca_tdf::sca_signal<double> gain_out;
   pshift_45* i_pshift;
   gain* i_gain;
   nshift_45* i_nshift;

 public:
   shift_45(sc_core::sc_module_name n,string mode, double w_if, double delta_r=0.0, double delta_c=0.0);
};

/**************************************Digital VCO**************************************************/


SCA_TDF_MODULE(d_vco) {

 public:
  sca_tdf::sca_in<double> in;        // input port
  sca_tdf::sca_out<double> out;	// output port

 private:

  int rate;		// The datarate of the output
  double freq;          // Central frequency [Hz]
  double kvco;          // Sensitivity [Hz/V]
  double gain;          // Gain

  int actsamples;
  int samples;
  double sample_time;
  double actval;

  void set_attributes() ;

  void initialize() ;

  void processing();

 public:

  d_vco(sc_core::sc_module_name n, double freq_data, int datarate, double _kvco, double _gain);
};


  /**************************analog VCO**************************************************/

SCA_TDF_MODULE(a_vco) {

 public:
  sca_tdf::sca_in<double> in;        // input port
  sca_tdf::sca_out<double> out;	// output port

 private:

  int rate;		// The datarate of the output
  double freq;          // Central frequency [Hz]
  double kvco;          // Sensitivity [Hz/V]
  double gain;          // Gain

  double sample_time;
  double stepsize;
  double actval;


  void set_attributes() ;

  void initialize() ;

  void processing();

 public:

  a_vco(sc_core::sc_module_name n, double freq_data, int datarate, double _kvco, double _gain);
};

/**********************************PLL****************************************************/

SC_MODULE(pll) {
 public:
  sca_tdf::sca_in<double>   ref;
  sca_tdf::sca_out<double>  vcoo, lpo;

 private:
  sca_tdf::sca_signal<double>  pco;

  phc* phc_sub;
  lp*  lp_sub;
  a_vco* vco_sub;

 public:
  pll(sc_core::sc_module_name n,double phc_gain,double lp_fc,double vco_freq,double kvco,double vco_gain,int vco_out_rate=1,int phc_out_rate=1);
};

/**********************************LNA****************************************/
SCA_TDF_MODULE(lna){

   public:
   // Ports
   sca_tdf::sca_in<double> in;    // Input
   sca_tdf::sca_out<double> out;  // Output

   private:
   double gain;              // Gain in dB
   double ip3;               // Third Input Intercept Point in dBm

   // Coefficients of output polynomial v = a*i - b*i*i - c*i*i*i
   double a;
   double b;
   double c;

   // Variablen for calculation of maximal input and output
   double input_max;
   double output_max;

   bool ideal;               //  ideal lna or not, true --> ideal

   #if _LNA_SUB_PROT
   // Wird zum Protokollieren verwendet
   ofstream prot;
   #endif

   public:
   // Constructor
   // name, gain in dB, ip3 in dBm, ideal (bool)
   lna(sc_core::sc_module_name n, double _gain, double _ip3, bool _ideal);

   // Destruktor
   ~lna();

   private:
   // Funktions which will be called in every simulation cycles
   void processing();

};
/********************************air***************************************************************/

SCA_TDF_MODULE(air) {

  sca_tdf::sca_in<double> in;	// input double (wave)
  sca_tdf::sca_out<double> out;	// output is the filtered wave

 private:
  int rate;
  double mean;
  double variance;
  double max;
  double min;
  double gain;
  string mode;

  void set_attributes();
  void initialize();
  void processing();

 public:
  air(sc_core::sc_module_name n, double atten, string n_art,double a,double b, int d_rate);
};

/************************** Differential encoder ***************************/

SC_MODULE(d_coding) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<bool> out;

 private:
  sca_tdf::sca_signal<bool> sig_xnor;
  sca_tdf::sca_signal<bool> sig_out;

  inout* inout_sub;
  d_flipflop* d_fp_sub;
  xnor* xnor_sub;

 public:
  d_coding(sc_core::sc_module_name n);
};


/************************** Differential decoder ***************************/
SC_MODULE(d_decoding) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<bool> out;
 private:
  sca_tdf::sca_signal<bool> sig_dff;
  d_flipflop* d_fp_sub;
  xnor* xnor_sub;
 public:
  d_decoding(sc_core::sc_module_name n);
};
/************************** AM(Amplitude Modulation) Modulator ***************************/   //TODO:AM demodulator
      // This class modulates the input signal to a carrier. The carrier frequency, amplitude and initial phase can be set with parameters. The parameter _offset adds
      // a constant value to the input signal. It is possible to increase the output data rate with the _rate parameter.
      SCA_TDF_MODULE(mod_am) {

	sca_tdf::sca_in<double>   in;
	sca_tdf::sca_out<double>   out;

      private:
	double ampl;
	double offset;
	double delta_T;

	double freq;
	int rate;

	void set_attributes();
	void initialize()  ;
	void processing();

      public:
	// The constructor takes the value of frequency, amplitute and offset of the carrier. Besides, user can set the data rate of output port with "d_rate".
	mod_am(sc_module_name, double _freq, double _ampl = 1.0, double _offset= 0.0, int _rate = 1);
      };


/************************** BPSK Modulator ***************************/
      // This class modulates an input bitstream to a carrier using binary phase shift keying. Frequency of the carrier and the data rate of output port can be set with
      // parameters.

      SCA_TDF_MODULE(bpsk) {

	sca_tdf::sca_in<bool> in;
	sca_tdf::sca_out<double> out;

      private:
	int out_rate;
	double sample_time;
	double actval;

	double freq;
	double step_size;
	double inp;

	void set_attributes();
	void initialize();
	void processing();

      public:
	// The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
	bpsk(sc_core::sc_module_name n, double _freq=1000,int _out_rate=1);
      };

/************************** BPSK Demodulator ***************************/
// This class demodulates a BPSK modulated signal. Frequency of the carrier and the data rate of input port can be set with
// parameters.

 SCA_TDF_MODULE(bpsk_de) {

  sca_tdf::sca_in<double> in;                                          // in and output port
  sca_tdf::sca_out<bool> out;

 private:
  int in_rate;
  double freq;
  sca_tdf::sca_signal<double> sig_sine;                                // signals to connect sub modules
  sca_tdf::sca_signal<double> sig_mix;
  sca_tdf::sca_signal<double> sig_lp;
  sca_tdf::sca_signal<double> sig_dd;
void processing();
  void set_attributes();
  sine* sine_sub;                                                 // declare sub modules
  mixer* mixer_sub;                   /*mixer* mixer_sub;*/
  lp* lp_sub;
  downsample* decider;
  compare* comp_sub;

 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the input port. They are 1000.0 and 1 by default, respectively.
  bpsk_de(sc_core::sc_module_name n, double _freq=1000,int _in_rate=1);
};

 /************************** BASK Modulator ***************************/
      // This class modulates an input bit stream to a carrier. Frequency, amplitude and initial phase of the carrier can be set with parameters. Parameter "_ampl1"
      // and "_ampl0" define the carrier amplitude for binary value 1 and 0, respectively.
      SCA_TDF_MODULE(mod_bask) {

	sca_tdf::sca_in<bool>  in;
	sca_tdf::sca_out<double>  out;

      private:
	double freq;
	double ampl1;
	double ampl0;
	double phi;
	int rate;

	double carrier;
	double delta_T;

	void set_attributes();
	void initialize();
	void processing();
      public:
	// The constructor takes the value of frequency, amplitutes and initial phase of the carrier.Data rate of output port can be set per "d_rate".They are set to
	// 1.0,0.0,0.0 and 1 by default, respectively.
	mod_bask(sc_core::sc_module_name nm, double _freq, double _ampl1 = 1.0, double _ampl0 = 0.0, double _phi = 0.0, int _rate = 1);
      };
/************************** BASK Demodulator ***************************/
      // This class demodulates a BASK modulated high frequency signal to a bit stream. Frequency of carrier signal can be set with parameter "_freq". A threshold signal level has to be set in
      // order to let the module work correctly. If detected signal is larger than the threshold value the module will output "true" and vice verse.
      SCA_TDF_MODULE(demod_bask) {

	sca_tdf::sca_in<double>  in;
	sca_tdf::sca_out<bool>  out;

      private:
	double freq;

	double level;
	int rate;
	double internal ;
	double internal2 ;
	sca_vector<double> A,B,S;
	sca_tdf::sca_ltf_nd ltf;

	void set_attributes() ;
	void initialize() ;
	void processing();
      public:
	// The constructor takes the value of frequency of the carrier, threshold signal level and the data rate of the output port.
	demod_bask(sc_core::sc_module_name nm, double _levle, double _freq = 10.0e3, int _rate = 1);
      };




/********************************quadrature mixer for receiver ************************************************/

SCA_TDF_MODULE(q_mixer_re) {

 public:

  sca_tdf::sca_in<double> in;             // input of frequency
  sca_tdf::sca_out<double> i_out;         // output of I
  sca_tdf::sca_out<double> q_out;         // output of Q

  sca_tdf::sca_in<double>* freq_in;       // pointer for optional input port of frequency

  void freq_con (sca_tdf::sca_signal<double>* f_in);

 private:

  double stepsize;	// the steps between the output token (i.e. if one token hase the
   			// value sin(x), then the next token has the value sin(x+stepsize)
  double actval;	// state variable to store the actual "x" value, if whe don't want
  			// to use sca_get_time()
  double sample_time;
  double freq;
  bool freq_recon;
  int rate;
  double theta;
  double e;
  double amp;

  void set_attributes();
  void initialize();
  void processing() ;

 public:
  q_mixer_re(sc_core::sc_module_name n,double _freq,double _amp,int data_rate=1,bool f_config=false,double _theta=0.0,double _e=0.0);
};


/**************************************quadrature mixer for transmitter****************************/

SCA_TDF_MODULE(q_mixer_tr) {

 public:

  sca_tdf::sca_out<double> out;             // output of frequency
  sca_tdf::sca_in<double> i_in;         // input of I
  sca_tdf::sca_in<double> q_in;         // input of Q

  sca_tdf::sca_in<double>* freq_in;       // pointer for optional input port of frequency

  void freq_con (sca_tdf::sca_signal<double>* f_in);

 private:

  double stepsize;	// the steps between the output token (i.e. if one token hase the
   			// value sin(x), then the next token has the value sin(x+stepsize)
  double actval;	// state variable to store the actual "x" value, if whe don't want
  			// to use sca_get_time()
  double sample_time;
  double freq;
  bool freq_recon;
  int rate;
  double theta;
  double e;
  double amp;

  void set_attributes();

  void initialize()  ;

  void processing() ;
 public:
  q_mixer_tr(sc_core::sc_module_name n, double _freq, double _amp,int data_rate=1,bool f_config=false,double _theta=0.0, double _e=0.0);
};


/****************************************** qam mapper for philip*************************************/

SCA_TDF_MODULE(qam_map_phi) {

  sca_tdf::sca_in<bool>   in;
  sca_tdf::sca_out<double> out_i;
  sca_tdf::sca_out<double> out_q;
  sca_tdf::sca_in<bool> reset;
  sca_tdf::sca_out<bool> valid;

 private:
  sc_bv<8> symbol;

  int rate_in;
  int bit_num;
  int wait;
  void set_attributes();

  void initialize()  ;

  void processing() ;

  double qam4(sc_bv<1>);
  double qam16(sc_bv<2>);
  double qam64(sc_bv<3>);
  double qam256(sc_bv<4>);

  double i_o;
  double q_o;

 public:
  qam_map_phi(sc_core::sc_module_name nm, int dim_cons);
};


/****** *********************************qam mapper ************************************/

SCA_TDF_MODULE(qam_map) {

  sca_tdf::sca_in<bool>   in;
  sca_tdf::sca_out<double> out_i;
  sca_tdf::sca_out<double> out_q;

 private:
  sc_bv<8> symbol;

  int rate_in;
  int bit_num;
  int wait;
  void set_attributes();

  void initialize()  ;

  void processing() ;

  double qam4(sc_bv<1>);
  double qam16(sc_bv<2>);
  double qam64(sc_bv<3>);
  double qam256(sc_bv<4>);

  double i_o;
  double q_o;

 public:
  qam_map(sc_core::sc_module_name nm, int _rate);
};






/****** **************************qam demapper**************************** ********/

SCA_TDF_MODULE(qam_demap) {

  sca_tdf::sca_in<double> in_i;
  sca_tdf::sca_in<double> in_q;
  sca_tdf::sca_out<bool>   out;

 private:

  int rate_in;
  int rate_out;

  double delta;                           //error tolerance
  void set_attributes();

  void initialize() ;

  void processing() ;

  sc_uint<1> qam4(int);
  sc_uint<2> qam16(int);
  sc_uint<3> qam64(int);
  sc_uint<4> qam256(int);

  int error_correction(double);

 public:
  qam_demap(sc_core::sc_module_name nm, int _rate=4, double _delta=1 );
};


/****************************** QPSK modulator*******************************/
// This class modulates an input bit stream to a carrier using quadrature phase shift keying. Frequency of the carrier and data rate of output port can be set with parameters.
SC_MODULE(qpsk) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<double> out;

 private:
  sca_tdf::sca_signal<bool> sig_dcode;                                                     // signal for connecting sub module
  sca_tdf::sca_signal<bool> sig_i;
  sca_tdf::sca_signal<bool> sig_q;
  sca_tdf::sca_signal<double> sig_n_i;
  sca_tdf::sca_signal<double> sig_n_q;

  s2p<bool,2>* s2p_sub;                                                               // declare sub module
  nrz* nrz_i_sub;
  nrz* nrz_q_sub;
  q_mixer_tr* mixer_sub;

 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  qpsk(sc_core::sc_module_name n, double _freq, int rate);
};

/****************************** QPSK demodulator*******************************/
// This class demodulates a QPSK modulated signal to a bit stream. Frequency of the carrier and data rate of input port can be set with parameters.
SC_MODULE(qpsk_de) {

  sca_tdf::sca_in<double> in;
  sca_tdf::sca_out<bool> out;

 private:

  sca_tdf::sca_signal<double> sig_i;                                                        // signal for connecting sub module
  sca_tdf::sca_signal<double> sig_q;
  sca_tdf::sca_signal<double> sig_di;
  sca_tdf::sca_signal<double> sig_dq;
  sca_tdf::sca_signal<double> sig_lp_i;
  sca_tdf::sca_signal<double> sig_lp_q;
  sca_tdf::sca_signal<bool> sig_ci;
  sca_tdf::sca_signal<bool> sig_cq;

  p2s<bool,2>* p2s_sub;                                                                // declare sub module
  lp* lp_i_sub;
  lp* lp_q_sub;
  downsample* decider_i;
  downsample* decider_q;
  q_mixer_re* mixer_sub;
  compare* comp_i_sub;
  compare* comp_q_sub;

 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  qpsk_de(sc_core::sc_module_name n, double _freq, int rate);
};




/****************************** OQPSK modulator*******************************/
// This class modulates an input bit stream to a carrier using offset quadrature phase shift keying. Frequency of the carrier and data rate of output port can be set with
// parameters

SC_MODULE(oqpsk) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<double> out;

 private:
  sca_tdf::sca_signal<bool> sig_dcode;                                                     // sub signals
  sca_tdf::sca_signal<bool> sig_i;
  sca_tdf::sca_signal<bool> sig_i2;
  sca_tdf::sca_signal<bool> sig_q;
  sca_tdf::sca_signal<bool> sig_q2;
  sca_tdf::sca_signal<bool> sig_dq;
  sca_tdf::sca_signal<double> sig_n_i;
  sca_tdf::sca_signal<double> sig_n_q;

  d_flipflop* delay_sub;                                                              // sub module
  s2p<bool,2>* s2p_sub;
  nrz* nrz_i_sub;
  nrz* nrz_q_sub;
  q_mixer_tr* mixer_sub;
  o2t* o2t_sub;
  o2t* o2t_i_sub;

 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  oqpsk(sc_core::sc_module_name n, double _freq, int rate);
};


/****************************** OQPSK demodulator*******************************/
// This class demodulates a OQPSK modulated signal to a bit stream. Frequency of the carrier and data rate of input port can be set with parameters.


SC_MODULE(oqpsk_de) {

  sca_tdf::sca_in<double> in;
  sca_tdf::sca_out<bool> out;

 private:

  sca_tdf::sca_signal<double> sig_i;                                                           //sub signals
  sca_tdf::sca_signal<double> sig_q;
  sca_tdf::sca_signal<double> sig_di;
  sca_tdf::sca_signal<double> sig_dq;
  sca_tdf::sca_signal<double> sig_lp_i;
  sca_tdf::sca_signal<double> sig_lp_q;
  sca_tdf::sca_signal<bool> sig_ci;
  sca_tdf::sca_signal<bool> sig_cq;
  sca_tdf::sca_signal<bool> sig_deq;
  sca_tdf::sca_signal<bool> sig_dei;
  sca_tdf::sca_signal<bool> sig_qto;

  sca_tdf::sca_signal<bool> sig_ito;
  p2s<bool,2>* p2s_sub;                                                                 // sub module
  lp* lp_i_sub;
  lp* lp_q_sub;
  downsample* decider_i;
  downsample* decider_q;
  q_mixer_re* mixer_sub;
  compare* comp_i_sub;
  compare* comp_q_sub;
  t2o* t2o_q_sub;
  t2o* t2o_i_sub;
  d_flipflop* delay_isub;

  d_flipflop* delay_qsub;

 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  oqpsk_de(sc_core::sc_module_name n, double _freq, int rate);
};


/****************************** DBPSK modulator *******************************/
// This class modulates an input bit stream to a carrier using differential binary phase shift keying. Frequency of the carrier and data rate of output port can be set with
// parameters.
SC_MODULE(dbpsk) {

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<double> out;
 private:
  sca_tdf::sca_signal<bool> sig_dcode;                                                 //signal to connect sub module

  d_coding* d_code_sub;                                                           // declare sub module
  bpsk* bpsk_sub;
 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  dbpsk(sc_core::sc_module_name n, double _freq, int rate);
};

/****************************** DBPSK demodulator *******************************/
// This class modulates an input bit stream to a carrier using differential binary phase shift keying. Frequency of the carrier and data rate of output port can be set with
// parameters.
 SC_MODULE(dbpsk_de) {

   sca_tdf::sca_in<double> in;
   sca_tdf::sca_out<bool> out;

 private:
   sca_tdf::sca_signal<bool> sig_dem;                                                 //signal to connect sub module

   d_decoding* d_decode_sub;                                                     // declare sub module
   bpsk_de* bpsk_de_sub;

 public:
   // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
   dbpsk_de(sc_core::sc_module_name n, double _freq, int rate);
 };




/****************************** DQPSK modulator*******************************/
// This class modulates an input bit stream to a carrier using differential quadrature phase shift keying. Frequency of the carrier and data rate of output port can be set with
// parameters.
SC_MODULE(dqpsk) {                          //one symbol duration delay

  sca_tdf::sca_in<bool> in;
  sca_tdf::sca_out<double> out;

 private:
  sca_tdf::sca_signal<bool> sig_i;                                                                 //sub signals
  sca_tdf::sca_signal<bool> sig_i_encode;
  sca_tdf::sca_signal<bool> sig_q;
  sca_tdf::sca_signal<bool> sig_q_encode;
  sca_tdf::sca_signal<double> sig_nrz_i;
  sca_tdf::sca_signal<double> sig_nrz_q;

  s2p<bool,2>* s2p_sub;                                                                       //sub module
  d_coding* encode_i_sub;
  d_coding* encode_q_sub;
  nrz* nrz_i_sub;
  nrz* nrz_q_sub;
  q_mixer_tr* mixer_sub;

 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  dqpsk(sc_core::sc_module_name n, double _freq, int rate);
};
/****************************** DQPSK demodulator*******************************/
// This class demodulates a DQPSK modulated signal to a bit stream. Frequency of the carrier and data rate of input port can be set with parameters.
SC_MODULE(dqpsk_de) {

  sca_tdf::sca_in<double> in;
  sca_tdf::sca_out<bool> out;

 private:

  sca_tdf::sca_signal<double> sig_i;                                                              //sub signals
  sca_tdf::sca_signal<double> sig_q;
  sca_tdf::sca_signal<double> sig_di;
  sca_tdf::sca_signal<double> sig_dq;
  sca_tdf::sca_signal<double> sig_lp_i;
  sca_tdf::sca_signal<double> sig_lp_q;
  sca_tdf::sca_signal<bool> sig_ci;
  sca_tdf::sca_signal<bool> sig_cq;
  sca_tdf::sca_signal<bool> sig_decq;
  sca_tdf::sca_signal<bool> sig_deci;
  sca_tdf::sca_signal<bool> sig_qto;

  sca_tdf::sca_signal<bool> sig_ito;                                                             //sub module
  p2s<bool,2>* p2s_sub;
  lp* lp_i_sub;
  lp* lp_q_sub;
  downsample* decider_i;
  downsample* decider_q;
  q_mixer_re* mixer_sub;
  compare* comp_i_sub;
  compare* comp_q_sub;
  d_decoding* decode_isub;
  d_decoding* decode_qsub;


 public:
  // The constructor takes the value of frequency of the carrier and the data rate of the output port. They are 1000.0 and 1 by default, respectively.
  dqpsk_de(sc_core::sc_module_name n, double _freq, int rate);
};

/************************** M-FSK Modulator ***************************/
      // The M-FSK modulator transmits signal through discrete frequency changes of a carrier wave. The lowest carrier frequency and the frequency interval between two
      // adjacent carriers can be set per the parameter _freq_basis and _freq_interval, respectively. With the template N you can set the bit width of input port. The
      // frequencies will be than assigned to each symbol automatically. For instance, by setting N=2, _freq_basis=1000 and _freq_interval=1000, the frequency 1000 will
      // be assigned to symbol "00" , the frequency 2000 will be assigned to symbol "01" and so on. It is also possible to increase the output data rate with the _out_rate
      // parameter.
      template <int N>                                                           // TODO:change the name of module to mfsk
	SCA_TDF_MODULE(fsk) {
      public:
	sca_tdf::sca_in<sc_bv<N> > in;
	sca_tdf::sca_out<double> out;

      private:
	double* step_size;

	int out_rate;
	double sample_time;
	double actval;

	double freq_actval;
	double freq_basis;

	void initialize()
	{
	  step_size = new double[(int)pow(2.0,(double)N)];	                // create an array to store values of stepsize for each symbol
	  sample_time = out.get_time().to_seconds();                               // get sample time of output port

	  for(int i=0;i<(int)pow(2.0,(double)N);i++)                            // calculate the step_size of carrier wave for each symbol and save them into array stepsize;(remove -1 to solve the problem)
	    {
	      step_size[i]=2.*M_PI*sample_time*(freq_basis+i*freq_actval);
	    }
	}

	void processing()
	{
	  sc_uint<N> in_bv;
	  int in_int=0;

	  in_bv=in.read();                                                     // read input port and convert the value to a integer value

	  in_int=(int)in_bv;
	  for(int i=0; i<out_rate; i++)  		                       // write <rate> data token for one symbol
	    {
	      out.write(sin( actval ),i);
	      actval+=step_size[in_int];
	    }
	}
      public:
	// The constructor takes the value of the lowest carrier frequency,the frequency interval and the data rate of the output port.
	fsk(sc_core::sc_module_name n, double _freq_basis=1000, double _freq_interval=1000, int _out_rate=1)
	  {
	    freq_basis=_freq_basis;
	    freq_actval=_freq_interval;
	    out_rate=_out_rate;
	    actval=0;
	  }
      };

      /************************** M-PSK Modulator ***************************/
      // The PSK modulator conveys data by changing, or modulating, the phase of a reference signal (the carrier wave).  The carrier frequency can be set per the
      // parameter _freq. With the template N you can set the bit width of input port. The phases will be than assigned to each symbol automatically.
      template <int N>                                                           //TODO:change the name of module to mpsk
	SCA_TDF_MODULE(psk) {
      public:
	sca_tdf::sca_in<sc_bv<N> > in;
	sca_tdf::sca_out<double> out;

      private:
	double* phi;

	int out_rate;
	double sample_time;
	double actval;
	double freq;
	double step_size;
	void initialize()
	{
	  phi = new double[(int)pow(2.0,(double)N)];                            // create an array to store phase for each symbol

	  sample_time = out.get_time().to_seconds();                               // read sample time
	  step_size=2.*M_PI*sample_time*freq;

	  for(int i=0;i<(int)pow(2.0,(double)N);i++)                            // calculate the phase for each symbol and save them into array phi;(remove -1 to solve the problem)
	    {
	      phi[i]=i*2*M_PI/(int)pow(2.0,(double)N);
	    }
	}

	void set_attributes(){
	};

	void processing()
	{
	  sc_uint<N> in_bv;
	  uint in_int=0;
	  sc_uint<N> in_pre;

	  in_bv=in.read();                                                      // read input port and convert the value to a integer value
	  if(in_int!=in_pre)                                                    // reset actval if a new symbol is received
	    actval=0;

	  in_pre=in.read();

	  in_int=(int)in_bv;
	  for(int i=0; i<out_rate; i++)  		                        // write <rate> data token for one symbol
	    {
	      out.write(sin( actval + phi[in_int]),i);
	      actval+=step_size;
	    }
	}

      public:
	// The constructor takes the value of carrier frequency and the data rate of the output port.
	psk(sc_core::sc_module_name n,double _freq,int _out_rate=1)
	  {
	    freq=_freq;
	    out_rate=_out_rate;
	    actval=0;
	  }
      };



/************************************Synchronizer***********************************/
template <typename T>
SCA_TDF_MODULE(buffer_synchronizer) {

   public:
   sca_tdf::sca_in<T> data_in;  	// incoming data token
   sca_tdf::sca_in<bool> valid_in;	// each data token is valid or invalid

   sca_tdf::sca_out<T> data_out;   	// the outgoing data token
   sca_tdf::sca_out<bool>* valid_outp;	// the valid now refers to the whole set of data token


   unsigned int rate;
   queue<T> q;			// a STL-queue for buffering
   T invalid_val;		// the value which replaces the invalid token
   bool val_out;

   void valid_out (sca_tdf::sca_signal<bool>* v_out);
   // the constructor takes the (input and output) data rate of the data token,
   // as well as the fill-in value for the invalid token. The latter defaults to 0
   buffer_synchronizer(sc_core::sc_module_name n, unsigned int _rate, bool vout=true, T _invalid_val  = static_cast<T>(0))
   {
     if(vout)
       {
	 valid_outp = new sca_tdf::sca_out<bool>;
	 val_out = true;
       }
     else
       val_out=false;
     rate = _rate;
     invalid_val = _invalid_val;
   };

   private:

   void set_attributes()
   {				// set the data rates
     data_in.set_rate(rate);
     valid_in.set_rate(rate);
     data_out.set_rate(rate);
     if(val_out)
       valid_outp->set_rate(1);
   }

   void processing()
   {
     for(unsigned int i = 0; i<rate; i++)		// read the data token
       if(valid_in.read(i)) q.push(data_in.read(i));	// and store them if they are valid

     if(q.size()>=rate)				// do we have enough token stored to write
     {						// to the output?
       for(unsigned int i = 0; i<rate; i++)	// => then do so!
       {
         data_out.write(q.front(),i);
         q.pop();
       }
       if(val_out)
	 valid_outp->write(true);		 // in this case, valid_out is true
     }
     else					// otherwise, we write rate copies of
     {						// invalid value to the output
       for(unsigned int i = 0; i<rate; i++)
         data_out.write(invalid_val,i);
       if(val_out)
	 valid_outp->write(false);			// in this case, valid_out is of course false
     }
   };
};

template<typename T>
void buffer_synchronizer<T>::valid_out (sca_tdf::sca_signal<bool>* v_out)
{
  (*valid_outp)(*v_out);
}


/**********************************************Network analyzer****************************************************/
/*sub module 1: measurement*/
SCA_TDF_MODULE(measure_sub) {

  sca_tdf::sca_in<double> in;			//In-Port
  sca_tdf::sca_in<double> ref;			//Out-Port

 private:
  ofstream output;			//Object for output file
  ostringstream ss1;			//String for amplitude response values
  ostringstream ss2;			//String for phase values

  ofstream measured_value;

  double scale_x;
  double scale_zero_x;			//Zero point for x-coordinate
  double scale_zero_y;			//Zero point for y-coordinate for amplitude response

  double time_buffer;			//Time buffer for period start
  double time;				//Time difference
  double period_number;			//Number of periods in one time window
  double frequenz;			//Current frequency values
  double start_f;				//Start frequency
  double stop_f;				//End frequency
  double step_f;				//Frequency steps
  double min_dB;				//Lower graph range for amplitude response
  double max_dB;				//Upper graph range for amplitude response
  double start_measurement_time;		//Time value at the beginning of measurement
  double stop_measurement_time;		//Time value at the end of measurement
  double start_measurement;		//Enable Measurement

  sca_vector <double> dut_buffer;		//Buffer for measurement values [0] = value; [1] = time
  sca_vector <double> ref_buffer;		//Buffer for reference values 	[0] = value; [1] = time

  double dB;				//Variable for amplitude response calculation
  double phase;				//Variable for phase calculation

  int find_beginning;			//Variable to find beginning of a period
  int x;					//Protection from incorrect measurement
  int start_ref;				//Start_ref == 1 --> the measurement of the ref signal begins
  int start_dut;				//Start_dut == 1 --> the measurement of the dut signal begins
  int start;				//If start == 0 --> the measurement of the current frequency stops

  void initialize();

  void processing();

  void measurement();

  void make_file();

  double sim_time();
 public:
  void finish();

  measure_sub(sc_core::sc_module_name n,  double _start_f, double _stop_f, double _step_f, double _min_dB, double _max_dB, double _period_number);

};

/*sub module 2: generator*/
SCA_TDF_MODULE(generator) {

  sca_tdf::sca_out<double> out;		//Out-Port to DUT
  sca_tdf::sca_out<double> ref;		//Out-Port to Measure-Module

 private:
  double amplitude;			//Amplitude of the signal
  double time_buffer;			//Time buffer for period start
  double time;				//Time difference
  double period_number;			//Number of periods in one time window
  double frequenz;			//Current frequency values
  double start_f;				//Start frequency
  double stop_f;				//End frequency
  double step_f;

  void initialize();
  void processing();

 public:
  generator(sc_core::sc_module_name n, double _amplitude, double _start_f, double _stop_f, double _step_f, double _period_number);
};

/*network analyzer*/

SC_MODULE(nwa) {

  sca_tdf::sca_in<double> in;		//In-Port to DUT
  sca_tdf::sca_out<double> out;	//Out-Port to DUT

 private:

  measure_sub* p_measure;
  generator* p_generator;

  sca_tdf::sca_signal<double> ref_signal;

 public:

  void finish();

  nwa(sc_core::sc_module_name n, double _amplitude, double _start_f, double _stop_f, double _step_f, double _min_dB, double
		_max_dB, double _period_number);
};

/**************************************** FFT/IFFT***************************************************************/

template <int N>

SCA_TDF_MODULE(fft_ifft) {

 public:
  sca_tdf::sca_in<double> in_real[N];
  sca_tdf::sca_in<double> in_imag[N];

  sca_tdf::sca_out<double> out_real[N];
  sca_tdf::sca_out<double> out_imag[N];

 private:
  string mode;
  int isign;
  double ZERO_THRESHOLD;
  bool debug;

  void set_attributes()
  {
  }

  void initialize()
  {
  }

  void processing() {

    vector<double> data;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;
    int istep;
    debug=false;

    for (int j=0;j<N;j++)
      {
	data.push_back(in_real[j].read());
	data.push_back(in_imag[j].read());
      }


    int data_size=data.size();
    if(debug==true)
      {
	cout<< "before rearrange" <<endl;
	for(int j=0;j< data_size;j++)
	  {
	    cout<< data.at(j) <<endl;
	  }
      }
    int j=1;                                                   //This is the bit-reversal section of the routine.

    for (int i=1; i < data_size; i+=2)
      {
	if (j > i) {
	  swap(data.at(j-1), data.at(i-1));
	  swap(data.at(j), data.at(i));
	}
	int m = data_size >> 1;
	while (m>=2 && j>m) {
	  j -= m;
	  m >>= 1;
	}
	j += m;
      }
    if(debug==true)
      {
	cout<< "after rearrange" <<endl;
	for(int j=0;j< data_size;j++)
	  {
	    cout<<data.at(j)<<endl;
	  }
      }
    int mmax=2;

    while (data_size>mmax)
      {
	istep = mmax << 1;
	//	theta = isign*(6.28318530717959/mmax);
	theta = isign*(2*M_PI/mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (int m=1; m<mmax; m+=2)
	  {
	    for (int i=m; i<=data_size; i+=istep)
	      {
		int j = i + mmax;
		tempr = wr*data[j-1]-wi*data[j];
		tempi = wr*data[j]+wi*data[j-1];
		data[j-1] = data[i-1]-tempr;
		data[j] = data[i]-tempi;
		data[i-1] += tempr;
		data[i] += tempi;
	      }
	    wr = (wtemp=wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	  }
	mmax = istep;
      }

    if(debug==true)
      {
	cout<< "after transform" <<endl;
	for(int j=0;j< data_size;j++)
	  {
	    cout<<data.at(j)<<endl;
	  }
      }

    if (mode=="IFFT")                           //devide by N if IFFT
      {
	for(int j=0;j< data_size;j++)
	  data[j]=data[j]/N;
      }
    if(debug==true)
      {
	cout<< "divided by N" <<endl;
	for(int j=0;j< data_size;j++)
	  {
	    cout<<data.at(j)<<endl;
	  }
      }

    int k=0;                                              //output
    for (int i=0;i<N;i++)
      {
	out_real[i].write(data.at(k));
	out_imag[i].write(data.at(k+1));
	k+=2;
      }
  }

 public:
  fft_ifft(sc_core::sc_module_name n, string _mode)
    {
      mode=_mode;
      if (mode=="FFT")
	isign=-1;
      else if(mode=="IFFT")
	isign=1;
      else
	cout<<"Error:mode can only be FFT or IFFT"<<endl;
    }
};

/*********************************************************Noise Generator****************************************************/

#ifndef TIMEADCDAC
#define TIMEADCDAC
#include "time.h"
#endif

class noiseGenerator {
    public:

    double nG_erww;	// Erwartungswert
    double nG_sigma;	// StdAbweichung
    double nG_factor;	// Gewichtung
    time_t ltime;	// Variable vom Typ time_t
            		// dient zum Merken der Systemzeit

    double uniform(double min, double max){

        double r, rmm;

        r = rand() / (RAND_MAX + 1.0);

        rmm = min + (max - min) * r;
        return rmm;
    }

    double uniform(double min, double max, double rnd_value){

        double r, rmm;

        // fuer rnd_value wird der Eingangswert des Wandlers verwendet:
        // da dieser Wert oft klein ist ==> Multiplikation mit 150
        srand(ltime * ( (unsigned int)(rnd_value * 150) ) );

        r = rand() / (RAND_MAX + 1.0);

        rmm = min + (max - min) * r;
        return rmm;
    }

    double gaussianValue(double value){

        double erg = 0;

        erg = (nG_factor / (nG_sigma * sqrt(2 * M_PI))) *
                exp(-pow(value - nG_erww, 2.0) / (2.0 * pow(nG_sigma, 2.0)));

        return erg;
    }


    double getSigma(double uref){

        int x = (int)uref;

        switch( x ){
            case 1:
                return uniform(0.2, 0.3);
                break;

            case 3:
                return uniform(0.5, 1.0);
                break;

            case 7:
                return uniform(2.0, 3.0);
                break;

            case 15:
                return uniform(4.0, 5.0);
                break;

            default:
                return uniform(4.0, 5.0);

        }
    }

    double getFactor(double uref){

        int x = (int)uref;

        switch( x ){
            case 1:
                return uniform(2.0, 3.0);
                break;

            case 3:
                return uniform(4.0, 10.0);
                break;

            case 7:
                return uniform(10.0, 25.0);
                break;

            case 15:
                return uniform(20.0, 60.0);
                break;

            default:
                return uniform(20.0, 60.0);
        }
    }

    //Konstruktor (siehe Dokumentation)
    noiseGenerator(unsigned int mul, double uref){

        time(&ltime);
        ltime = ltime * mul;
        srand(ltime);

        nG_erww = uref / 2.0;
        nG_sigma = getSigma(uref);
        nG_factor = getFactor(uref);

    }
};

/*********OFDM_Transmitter***********/
template <int N>
SC_MODULE(ofdm_se){

 public:
  //Ports:
  sca_tdf::sca_in<bool> in;           // input port
  sca_tdf::sca_out<double> out;        // output port


  //signals
  sca_tdf::sca_signal<bool> sig_pa[N];          //signals on s2p output ports
  sca_tdf::sca_signal<double> sig_real[N];      //signals on q_mapper i output ports
  sca_tdf::sca_signal<double> sig_imag[N];      //signals on q_mapper q output ports
  sca_tdf::sca_signal<double> sig_out_real[N];  //signals on fft real output ports
  sca_tdf::sca_signal<double> sig_out_imag[N];  //signals on fft imag output ports
  sca_tdf::sca_signal<double> sig_out_i;        //signals on p2s i output ports
  sca_tdf::sca_signal<double> sig_out_q;        //signals on p2s q output ports
  //  sca_tdf::sca_signal<double> sig_out;          //signals on transmitter output port

  // sca_tdf::sca_signal<double> sig_noise;
  // sca_tdf::sca_signal<double> sig_out_noise;
 private:
  //module instantiation
  s2p<bool,N>* s2p_sub;
  qam_map* qam_mapper_sub[N];
  fft_ifft<N>* ifft_sub;
  p2s<double,N>* p2s_r_sub;
  p2s<double,N>* p2s_i_sub;
  q_mixer_tr* q_mixer_tr_sub;

  // int qam_p_num;
  // noise_gauss* g_noise_sub;
  // add* add_sub;

 public:
  //Constructor
  ofdm_se(sc_core::sc_module_name n,double mixer_fc,int qam_p_num,double bit_f,int dout_rate,double _amp=1,int s2p_or=1,bool mixer_config=false,int p2s_rate=1)
  {
    //  qam_p_num=(int)pow(2.0,(double)(N/2));
    int mixer_rate=(int)floor(dout_rate*log2(qam_p_num)*mixer_fc/bit_f);
    s2p_sub = new s2p<bool,N>("s2p_sub",1);
    s2p_sub->in(in);
    for(int i=0;i<N;i++)
      s2p_sub->out[i](sig_pa[i]);

    for(int i=0;i<N;i++)
      {
        sc_core::sc_module_name nm(((string)"qam_mapper_sub_"+(char)(i+65)).c_str());
	qam_mapper_sub[i] = new qam_map(nm ,qam_p_num);
	qam_mapper_sub[i]->in(sig_pa[i]);
	qam_mapper_sub[i]->out_i(sig_real[i]);
	qam_mapper_sub[i]->out_q(sig_imag[i]);
      }

    ifft_sub = new fft_ifft<N>("ifft_sub","IFFT");
    for(int i=0;i<N;i++)
      {
	ifft_sub->in_real[i](sig_real[i]);
	ifft_sub->in_imag[i](sig_imag[i]);
	ifft_sub->out_real[i](sig_out_real[i]);
	ifft_sub->out_imag[i](sig_out_imag[i]);
      }

    p2s_r_sub = new p2s<double,N>("p2s_r_sub",p2s_rate);
    for(int i=0;i<N;i++)
      {
	p2s_r_sub->in[i](sig_out_real[i]);
      }
    p2s_r_sub->out(sig_out_i);

    p2s_i_sub = new p2s<double,N>("p2s_i_sub",p2s_rate);
    for(int i=0;i<N;i++)
      {
	p2s_i_sub->in[i](sig_out_imag[i]);
      }
    p2s_i_sub->out(sig_out_q);

    q_mixer_tr_sub = new q_mixer_tr("q_mixer_tr_sub",mixer_fc,_amp,mixer_rate,mixer_config,0.,0.);
    q_mixer_tr_sub->i_in(sig_out_i);
    q_mixer_tr_sub->q_in(sig_out_q);
    q_mixer_tr_sub->out(out); /***/

    /*
    g_noise_sub = new noise_gauss("noise_sub",10,0.0,mixer_rate);
    g_noise_sub->out(sig_noise);

    add_sub = new add("add_noise",mixer_rate);
    add_sub->in1(sig_noise);
    add_sub->in2(sig_out);
    add_sub->out(out);
    */
  }
};
//*********************************************************************************************
/*********OFDM_Receiver***********/

template <int N>
SC_MODULE(ofdm_re){

 public:
  //Ports
  sca_tdf::sca_in<double> in;           // input port
  sca_tdf::sca_out<bool> out;            // output port


  //Signals
  sca_tdf::sca_signal<double> sig_in_i;
  sca_tdf::sca_signal<double> sig_in_q;
  sca_tdf::sca_signal<double> sig_i_lp;
  sca_tdf::sca_signal<double> sig_q_lp;
  sca_tdf::sca_signal<double> sig_i_rounded;
  sca_tdf::sca_signal<double> sig_q_rounded;
  sca_tdf::sca_signal<double> sig_re_real[N];
  sca_tdf::sca_signal<double> sig_re_imag[N];
  sca_tdf::sca_signal<double> sig_re_out_real[N];
  sca_tdf::sca_signal<double> sig_re_out_imag[N];
  sca_tdf::sca_signal<bool> sig_demapper[N];
  sca_tdf::sca_signal<bool> sig_received;
 private:
  //module instantiation
  q_mixer_re* q_mixer_re_sub;
  lp* i_lp;
  lp* q_lp;
  downsample* i_round;
  downsample* q_round;
  s2p<double,N>* s2p_r_sub;
  s2p<double,N>* s2p_i_sub;
  fft_ifft<N>* fft_sub;
  qam_demap* qam_demapper_sub[N];
  p2s<bool,N>* p2s_sub;
  //  int demap_p;
 public:
  //Constructor
  ofdm_re(sc_core::sc_module_name n,double mixer_fc,int demap_p,double bit_f,int dout_rate,double _amp=1, int p2s_rate=1)
    {
      //   demap_p=(int)pow(2.0,(double)(N/2));
      int mixer_rate=(int)floor(dout_rate*log2(demap_p)*mixer_fc/bit_f);
      q_mixer_re_sub = new q_mixer_re("mixer_sub",mixer_fc,_amp,mixer_rate,false,0.,0.);
      q_mixer_re_sub->in(in);
      q_mixer_re_sub->i_out(sig_in_i);
      q_mixer_re_sub->q_out(sig_in_q);

      i_lp = new lp("i_lp",mixer_fc/1000.0);
      i_lp->in(sig_in_i);
      i_lp->out(sig_i_lp);

      q_lp = new lp("q_lp",mixer_fc/1000.0);
      q_lp->in(sig_in_q);
      q_lp->out(sig_q_lp);

      i_round = new downsample("i_round", mixer_rate, mixer_rate);
      i_round->in(sig_i_lp);
      i_round->out(sig_i_rounded);

      q_round = new downsample("q_round", mixer_rate, mixer_rate);
      q_round->in(sig_q_lp);
      q_round->out(sig_q_rounded);

      s2p_r_sub = new s2p<double,N>("s2p_r_sub",1);
      s2p_r_sub->in(sig_i_rounded);
      for(int i=0;i<N;i++)
	{
	  s2p_r_sub->out[i](sig_re_real[i]);
	}

      s2p_i_sub = new s2p<double,N>("s2p_i_sub",1);
      s2p_i_sub->in(sig_q_rounded);
      for(int i=0;i<N;i++)
	{
	  s2p_i_sub->out[i](sig_re_imag[i]);
	}

       fft_sub = new fft_ifft<N>("fft_sub","FFT");
       for(int i=0;i<N;i++)
	 {
	   fft_sub->in_real[i](sig_re_real[i]);
	   fft_sub->in_imag[i](sig_re_imag[i]);
	   fft_sub->out_real[i](sig_re_out_real[i]);
	   fft_sub->out_imag[i](sig_re_out_imag[i]);
	 }

       for(int i=0;i<N;i++)
	 {
           sc_core::sc_module_name nm(((string)"qam_demapper_sub_"+(char)(i+65)).c_str());
	   qam_demapper_sub[i] = new qam_demap(nm,demap_p);
	   qam_demapper_sub[i]->out(sig_demapper[i]);
	   qam_demapper_sub[i]->in_i(sig_re_out_real[i]);
	   qam_demapper_sub[i]->in_q(sig_re_out_imag[i]);
	 }

       p2s_sub = new p2s<bool,N>("p2s_sub",p2s_rate);
       for(int i=0;i<N;i++)
	 {
	   p2s_sub->in[i](sig_demapper[i]);
	 }
       p2s_sub->out(out);
    }
};

/*********************************************Eye Diagram******************************************************


TODO: separate it to header and .cpp
 Description of constructor:
 ______________________________
 eyediag(sc_core::sc_module_name n, double periodetime_, double sigamp_, int periodes_, int delay_)

 1. n 			<string> 	Modulname
 2. periodetime_	<double>	Symbol Period[s]
 3. sigamp_		<double>	Amplitude of signal[Vpp]
 4. periodes_		<int>		Number of periods to be plotted
 5. delay_		<int>		Number of periods to be ignored before plot
 6. in_rate             <integer>       data rate of input port, default value is 10
 ___________________________________________________________________________________________________

 INFO: 	In the end of your program(in main.cpp) the method finish() must be called to generate the SVG-File and finish
        the file handling.
***************************************************************************************************************/
SCA_TDF_MODULE(eyediag)
{
  sca_tdf::sca_in<double>   in;
 private:
  double periodetime;
  double sigamp;
  int wavemax;
  int delay;
  int lastx;
  int data_rate;

  ofstream output;
  ostringstream ss1;
 public:
  eyediag(sc_core::sc_module_name n, double periodetime_, double sigamp_, int periodes_, int delay_, int in_rate=10);

 private:
  void set_attributes();
  void processing();

public:
  void finish();

};


/**************************************************Scatter-Plot**************************************************
 Description of constructor:
 ______________________________
 eyediag(sc_core::sc_module_name n, double periodetime_, double sigamp_, int periodes_, int delay_)

 1. n 			<string> 	Modulname
 2. sigamp_		<double>	Signalamplitude [V]
 3. in_rate             <integer>       data rate of input ports(I and Q port), default value is 10
 ___________________________________________________________________________________________________

 INFO: 	In the end of your program(in main.cpp) the method finish() must be called to generate the SVG-File and finish
        the file handling.
*****************************************************************************************************************/

SCA_TDF_MODULE(scatter)
{
  sca_tdf::sca_in<double>   inI;
  sca_tdf::sca_in<double>   inQ;

 private:
  double sigamp;
  int data_rate;
  ofstream output;
  ostringstream ss1;

 public:
  scatter(sc_core::sc_module_name n, double sigamp_, int in_rate=10);

 private:
  void set_attributes();
  void processing();

 public:
  void finish()	;
};


/*************************************************ADC******************************************************/




#ifndef MATHADC
#define MATHADC
#include "math.h"
#endif

#ifndef NOISEADCDAC
#define NOISEADCDAC

#endif

#define VERBOSE 0
#define RNDA 20


template <int N>

SCA_TDF_MODULE(adc) {

  sca_tdf::sca_in<double> in;	         // Input port
  sca_tdf::sca_out<sc_bv<N> > out;            // Output port

  //variablen

  double gain_err;		// gain error
  double offset_err;		// offset error
  double u_ref;		        // converter's reference voltage
  double lsb;			// least significant bit
  long max;			// maximal output value
  int nl_mode;			// mode selection for liniarity error

  sc_bv<N> bv_erg;		// 16- bit Vektor for conv. value
  int first_inl;		// Merker fuer ersten Aufruf von INL


  //noise			// Objekt fuer Erzeugung von Zufallswerten
  // (see file: noise_ADCDAC.h)
  noiseGenerator *myNoise;

  void processing()
  {
    double analog;
    double temp;
    long erg;

    analog = in.read();	                //  read analog input

    temp = gainError(gain_err, analog); // temp with gain error
    temp = temp + offset_err;      	// temp with offset error
    temp = temp + NlError(nl_mode, analog);	// Nichtlinearitaetsfehler add.

    erg = roundValue(temp / lsb);	// calculate digital value and rounding
    erg = maximumValue(erg);    	// output limitation
    bv_erg = erg;	        	// save as bBitvektor

    out.write(bv_erg);          	// write result to output port

#if VERBOSE
    cout << " G_E: " << gain_err << " O_E: " << offset_err;
    cout << " analog: " << analog << " ADC: " << erg;
#endif
  }

  void initialize(){

    // objekt of class "noiseGenerator"
    myNoise = new noiseGenerator(RNDA, u_ref);
    // analog value according to the maximal digital output of adc
    max = (long)(pow(2.0, N - 1.0));
    // 1bit lost because of sign
    lsb = u_ref / (max);	// calculate lsb of the adc

    first_inl = 0;          // to note the first call of INL()


    //calculate a random number if parameter gain_err or offset_err is not 0;
    if(gain_err != 0)
      gain_err = myNoise->uniform(-(gain_err * lsb),(gain_err * lsb), gain_err);

    if(offset_err != 0)
      offset_err = myNoise->uniform(-(offset_err * lsb),(offset_err * lsb),offset_err);
  }

  void set_attributes(){}

  long roundValue(double val){
    return (long)lround(val);
  }

  double gainError(double err, double val){
    return(val * (1.0 + err));
  }

  long maximumValue(long value){

    if( abs(value) > max-1){

      if(value > 0)
	value = max-1;	//by positiv results => pos. maximum - 1 (according to the maximum analog value Uref-Uls)


      else
	value = - max;	//by negativ results => neg. Maximum
      // also -Uref
    }
    return value;
  }

  /*******************************************************************************
   !!!This function is still under test!!!
    Funktion: NlError
    zur Modellierung der Nichtlinearitaetsfehler eines ADC
    gibt den Wert des Gesamtfehlers zur��ck
    mode:
    0... keine Nichtlinearitaetsfehler
    1... integraler Nichtlinearitaetsfehler
    2... differentieller Nichtlinearitaetsfehler
    3... beide Fehler

    Paramter:
    int mode: Modus - siehe oben
    double val: Eingangspegel
  ************************************************************************************/
  double NlError(int mode, double val){

    double err = 0.0;

    switch(mode){

    case 0:		//keine Nichtlinearitaeten
      break;

    case 1:		//integrale Nichtlinearitaet
      err = err + Inl(val);
      break;

    case 2:		//differentielle Nichtlinearitaet
      err = err + Dnl(val);
      break;

    case 3:		// beide Nichtlinearitaeten
      err = err + Inl(val);
      err = err + Dnl(val);
      break;

    default:
      break;
    }

#if VERBOSE
    cout << " NL_E: " << err;
#endif

    return err;
  }

  double Inl(double value){

    double ret, sign;

    // wenn INL(value) zum ersten Mal aufgerufen wird:
    // bestimmen des Vorzeichens der Abweichung fuer den pos. Bereich
    if(first_inl == 0){
      // Zufallswert generieren
      sign = myNoise->uniform(-1.0, 1.0, RNDA);

      // je nach Vorzeichen des Zufallswertes folgt das
      // Vorzeichen der Abweichung im pos. Bereich
      if(sign >= 0.0)
	first_inl = 1;
      else
	first_inl = -1;
    }

    //liefert Wert von Gausskurve (Verteilung) in Abh. von value
    ret = myNoise->gaussianValue(abs(value));
    ret = ret * lsb * (double)first_inl;	// Vielfaches von Ulsb

    //Vorzeichen des Fehlers fuer jeweiligen Bereich der Kennlinie
    if(value < 0)
      ret = (-1.0) * ret;
    else
      ret = ret;

    return ret;
  }

  double Dnl(double value){

    double err, ret;
    double step_o, step_u, step_o_err, step_u_err;
    long round;

    round = roundValue(value / lsb);  // gerundeter digitalisierter wert

    // ermittelt den Fehler (err) fuer den Wert value
    // wenn Fehler = 0.0 => unveraendert zurueckgeben
    if(( err = myNoise->uniform(-1.10, 1.10, (double)round) ) == 0.0)
      return 0.0;

    if(abs(err) > 1.0 && err > 0)	// begrenze Fehler auf max. 1 Ulsb
      err = 1.0;

    if(abs(err) > 1.0 && err < 0)	// begrenze Fehler auf max. 1 Ulsb
      err = -1.0;

    err = err * lsb;		// Fehler in Ulsb schritten
    ret = 0.0;

    // entscheidung ueber neuen wert bei geaenderter stufenbreite
    // siehe dokumentation

    step_o = round * lsb + lsb / 2.0;
    step_u = (round - 1) * lsb + lsb / 2.0;
    step_o_err = step_o + err;
    step_u_err = step_u + err;

    if( step_o_err > step_u_err){

      if( (value >= step_u) && (value <= step_u_err) )
	ret = - lsb;

      if( (value >= step_o_err) && (value <= step_o) )
	ret = lsb;

      if( (value > step_u_err) && (value < step_o_err) )
	ret = 0.0;
    }
    else {

      if( value <= step_o_err )
	ret = - lsb;

      if( value >= step_u_err )
	ret = lsb;

      if( (value > step_o_err) && (value < step_u_err) )
	ret = 0.0;
    }

    // Rueckgabe des abgeaenderten Analogwertes
    return ret;
  }

  // Konstruktor (siehe Dokumentation)
  adc(sc_core::sc_module_name n, double uref=1.0, double gain_e=0.0, double offset_e=0.0, int nl_m=0)
    {
      gain_err = gain_e;
      offset_err = offset_e;
      nl_mode = nl_m;
      u_ref = uref;
    }
};



/********************************************DAC*************************************************/

#ifndef NOISEADCDAC
#define NOISEADCDAC
#endif

#define VERBOSE 0
#define RNDD 30

template <int N>
SCA_TDF_MODULE(dac){

    sca_tdf::sca_in<sc_bv<N> > in;	// input port
    sca_tdf::sca_out<double> out;    // output port

    //variablen


    double u_ref;	        // Referenzspannung
    double lsb;		        // Spannungswert des LSBs
    double gain_err;		// Verstaerkungsfehler
    double offset_err;		// Offset- Fehler
    long max;			// Maximal darstellbarer Wert (Aufloesungsbits)
    int nl_mode;		// Nichtlinearitaets- Modus

				// Objekt fuer Erzeugung von Zufallswerten
				// (see file: noise_ADCDAC.h)

    int first_inl;	        // Merker fuer ersten Aufruf von INL
    noiseGenerator *myNoise;

    void processing()
    {
        long digital;
        double analog;
        double erg;

        digital = in.read().to_int();	// "digitalen" Eingangswert lesen
        // ideal gewandelten Analogwert berechnen
        analog = (double)(digital * lsb);

        // Analogwert mit Verstaerkungsfehler manipulieren
        erg = gainError(gain_err, analog);
        // Analogwert mit Offsetfehler manipulieren
        erg = erg + offset_err;
        // Nichtlinearitaetsfehler hinzuaddieren (je nach Modus)
        erg = erg + NlError(nl_mode, analog);

        erg = roundValue(erg);       // Analogwert auf Ulsb- Schritte runden

        erg = maximumValue(erg);     // Wertbegrenzung

        out.write(erg);              // auf Ausgangsport schreiben

    #if VERBOSE
        cout << " G_E: " << gain_err << " O_E: " << offset_err;
        cout << " digital: " << digital << " DAC: " << erg << endl;
    #endif
    }

    void initialize(){

        myNoise = new noiseGenerator(RNDD, u_ref);  // Objekt "noiseGenerator"
        max = (long)(pow(2.0, N - 1.0));  // Maximaler Ergebniswert
        lsb = u_ref / (max);             // Spannungswert fuer LSB
        first_inl = 0;                   // merkt ersten Aufruf von INL()

        // berechnen eines Zufallswertes (fkt uniform())
        // wenn Parameter gain_err
        // bzw. Parameter offset_err nicht 0 sind.
        if(gain_err != 0)
            gain_err = myNoise->uniform(-(gain_err * lsb),
                (gain_err * lsb), gain_err);

        if(offset_err != 0)
            offset_err = myNoise->uniform(-(offset_err * lsb),
                (offset_err * lsb), offset_err);
    }

    void set_attributes(){


    }

    /*******************************************************************************
    Funktion: roundValue
        rundet einen gegebenen double Wert auf double mit vollem Ulsb - Schritt

    Paramter:
        double val: zu rundender double Wert
    ********************************************************************************/


    double roundValue(double val){

        return lround(val/lsb) * lsb;
    }


    /************************************************************************************
    Funktion: gainError
        zur Modellierung des Verstaerkungsfehlers eines DACs
        der double Wert "err" wird mit dem gegebenen
        Eingangspegel multipliziert
        ==> "err + 1" steht also fuer die fehlerhafte Steigung
        der Uebertragungskennlinie

    Paramter:
        double val: Eingangswert
        double err: abweichender Wert von der idealen Steigung 1.0
    ***************************************************************************************/

    double gainError(double err, double val){

        return (val * (err + 1.0));
    }


    /**********************************************************************************
    Funktion: maximumValue
        der Ausgabewert (analog) des DACs ist natuerlich durch
        die Aufloesung begrenzt bei Ueberschreitung des Maximalwertes
        wird der Maximalwert ausgegeben
        (sowohl positiver als auch negativer Bereich)

        Wertebereich von Uref - Ulsb bis -Uref
        MSB bezeichnet Vorzeichenbit

    Paramter:
        double val: "analoger" Ergebniswert
    **************************************************************************************/
    double maximumValue(double value){

        if( abs(value) > (u_ref - lsb)){

            if(value > 0)
                //bei positiven Ergebnis => pos. Maximum - 1
                // da Uref-Ulsb
                value = u_ref - lsb;

            else
            //bei negativen Ergebnis => neg. Maximum
    		// also - Uref
            value = - u_ref;
        }

        return value;		//Rueckgabe des begrenzten Wertes
    }

    /****************************************************************************************
    Funktion: NlError
        zur Modellierung der Nichtlinearitaetsfehler eines DAC
        gibt den Wert des Nichtlinearitaetsfehlers zurueck
        mode: 	0... keine Nichtlinearitaetsfehler
                1... integraler Nichtlinearitaetsfehler
                2... differentieller Nichtlinearitaetsfehler
                3... beide Fehler

    Paramter:
        int mode: Modus - siehe oben
        double val: Eingangspegel
    ************************************************************************************************/
    double NlError(int mode, double val){

        double err = 0.0;

        switch(mode){

            case 0:		//keine Nichtlinearitaeten
                break;

            case 1:		//integrale Nichtlinearitaet
                err = err + Inl(val);
                break;

            case 2:		//differentielle Nichtlinearitaet
                err = err + Dnl(val);
                break;

            case 3:		// beide Nichtlinearitaeten
                err = err + Inl(val);
                err = err + Dnl(val);
                break;

            default:
                break;
    }

    #if VERBOSE
        cout << " NL_E: " << err;
    #endif

        return err;
    }

    /*********************************************************************************************
    Funktion: Inl
        zur Modellierung des integralen Nichtlinearitaetsfehlers.
        in Abhaengigkeit des Eingangspegels wird ein groesserer oder kleinerer
        Fehler zum Eingangspegel addiert (Teile von Ulsb).

        in der Mitte der Charakteristik (jeweils pos, oder neg. Bereich)
        ist der I-Nl jeweils am groessten ==> daher aus Gausskurve.

    Paramter:
    double value: Eingangswert
    ***********************************************************************************************/
    double Inl(double value){

        double ret, sign;

        // wenn INL(value) zum ersten Mal aufgerufen wird:
        // bestimmen des Vorzeichens der Abweichung fuer den pos. Bereich
        if(first_inl == 0){
            sign = myNoise->uniform(-1.0, 1.0, RNDD);

            // je nach Vorzeichen des Zufallswertes folgt das
            // Vorzeichen der Abweichung im pos. Bereich
            if(sign >= 0.0)
                first_inl = 1;
            else
                first_inl = -1;
        }

        //liefert Wert von Gausskurve (Verteilung) in Abh. von value
        ret = myNoise->gaussianValue(abs(value));
        ret = ret * lsb * (double)first_inl;	// Vielfaches von Ulsb

        //Vorzeichen des Fehlers fuer jeweiligen Bereich der Kennlinie
        if(value < 0)
            ret = (-1.0) * ret;

        return ret;
    }

    /**********************************************************************************
    Funktion: Dnl
        zur Modellierung des differentiellen Nichtlinearitaetsfehlers.
        fuer jeden Eingangspegel wird ein Fehlerwert ermittelt.
        dieser Wert bleibt fuer die gesamte Simulationszeit derselbe.

        der diff. Nichtlinearitaetsfehler ist gleichverteilt
        ueber die gesamte Charakteristik und wirkt sich als zusaetzlicher
        Additiver Anteil auf die Kennlinie aus hat der Fehler die
        Groesse von Ulsb, so kann auch ein "Nichtmonotonie"- Fehler
        auftreten!

    Paramter:
        double value: Eingangswert
    *******************************************************************************************/


    double Dnl(double value){

    double err;

        // ermittelt den Fehler (err) fuer den Wert value
        // wenn Fehler = 0.0 => unveraendert zurueckgeben
        if(( err = myNoise->uniform(-1.50, 1.50, value) ) == 0.0)
            return 0.0;

        if(abs(err) > 1.50 && err > 0)	// begrenze Fehler auf max. 1.5 Ulsb
            err = 1.50;

        if(abs(err) > 1.50 && err < 0)	// begrenze Fehler auf max. 1.5 Ulsb
            err = -1.50;

        err = err * lsb;		// Fehler in Ulsb schritten

        return err;
    }

    //Konstruktor (siehe Dokumentation)
    dac(sc_core::sc_module_name n, double uref=1.0,
		double gain_e=0., double offset_e=0., int nl_m=0.)
    {
        gain_err = gain_e;
        offset_err = offset_e;
        nl_mode = nl_m;
        u_ref = uref;
    }
};

/********************************************Jitter*************************************************/

SCA_TDF_MODULE(jit) {

    sca_tdf::sca_in<double>  in;
    sca_tdf::sca_out<double>  out;

private:
    sca_tdf::sca_ltf_nd  ltf_1;
    sca_core::sca_time tdelay;

    sca_util::sca_vector<double> num, den;
    sca_util::sca_vector<double> state;
    double max_jit;

    void initialize() ;

    void processing();

    void set_attributes();

public:
    jit(sc_core::sc_module_name nm, double _max_jit);
};







  }//bb
}//TUV_ams_lib


#endif
