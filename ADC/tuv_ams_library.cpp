/*******************************************************************************************************
Copyright (c) 2008, Institute of Computer Technology (University of Technology Vienna)
Authors:
        Markus Damm <damm@ict.tuwien.ac.at>
        Jan Haase <haase@ict.tuwien.ac.at>
        Jiong Ou <ou@ict.tuwien.ac.at>
        Joseph Wenninger <wenninger@ict.tuwien.ac.at>
For the license please read LICENSE_BB.txt
**********************************************************************************************************/

#include "systemc-ams.h"
#include "tuv_ams_library.h"

namespace TUV_ams_lib {
  namespace bb{





 /**********************************Multiplier**********************************/




    void mul::set_attributes() { }

    void mul::initialize() { }

    void mul::processing() {

      out.write( in1.read() * in2.read() ) ;

    }

    mul::mul(sc_core::sc_module_name nm) { }

    /*******************************Divider*************************************/


    void divi::processing() {

      out.write(in1.read()/in2.read() ) ;

    }

    divi::divi(sc_core::sc_module_name nm) { }

    /********************************Adder************************************/

    void add::set_attributes() {in1.set_rate(data_rate);
      in2.set_rate(data_rate);
      out.set_rate(data_rate);
    }

    void add::initialize() { }

    void add::processing() {
      for(int i=0;i<data_rate;i++)
	{
	  out.write( in1.read(i) + in2.read(i),i ) ;
	}
    }

    add::add(sc_core::sc_module_name nm, int rate) {data_rate=rate; }

    /********************************Subtractor************************************/

    void sub::set_attributes() { }

    void sub::initialize() { }

    void sub::processing() {

      out.write( in1.read() - in2.read() ) ;

    }

    sub::sub(sc_core::sc_module_name nm) { }

/***********************************logarithm**********************************************************/

void logarithm::initialize(){};

void logarithm::processing()
  {
    double base_tmp;


    if(base=="10")
      out.write((BB_DBL)log10((BB_DBL)in.read()));
    else if(base=="e")
      out.write((BB_DBL)log((BB_DBL)in.read()));
    else if(base=="2")
      out.write((BB_DBL)log2((BB_DBL)in.read()));
    else
      {
	base_tmp=atof(base.data());
	if(base_tmp==0||base_tmp==1||base_tmp<0)
	  cout<<"base must be positive number and must not equal to 0 and 1"<<endl;
	else
	  out.write(log10(in.read())/log10(base_tmp));
      }
  }

logarithm::logarithm(sc_core::sc_module_name n, string _base="10")
{
  base =_base;
}


/********************************************tangnes**********************************************/

void tang::processing()
  {
    out.write(tanh(in.read()));
  }

tang::tang(sc_core::sc_module_name n){};

/********************************************* XNOR**************************************************/

    void xnor::processing()
    {
      if(in1.read()==in2.read())
	out.write(true);
      else
	out.write(false);
    }

    xnor:: xnor(sc_core::sc_module_name n)
    {
    }
    /************************** Compare ***************************/

    void compare::processing()
    {
      if(in.read()>=threshold)
	out.write(true);
      else
	out.write(false);
    }

    compare::compare(sc_core::sc_module_name,double _threshold){threshold=_threshold;}



/************************** Delay ***************************/

    void d_flipflop::set_attributes()
    {
      out.set_delay(delay);
    }

    void d_flipflop::initialize()
    {
      out.write(initial);
    }

    void d_flipflop::processing()
    {out.write(in.read());}

    d_flipflop::d_flipflop(sc_core::sc_module_name n,bool _init,int _delay)
    {
      initial=_init;
      delay=_delay;
    }


/****************************** NRZ *******************************/
    void nrz::processing()
    {
      if(in.read())
	out.write(vout);
      else
	out.write(-1.0*vout);
    }

nrz::nrz(sc_core::sc_module_name,double _vout){vout=_vout;}







/********************************Saturation************************************/

    void saturation::set_attributes() {}
    void saturation::initialize() {  }
    void saturation::processing() {

      if ( in.read() > OUT_MAX )
	out.write(OUT_MAX) ;

      else if ( in.read() < OUT_MIN )
	out.write(OUT_MIN);

      else
	out.write(in.read());

    }

    saturation::saturation(sc_core::sc_module_name nm, double _OUT_MAX, double _OUT_MIN) {
      OUT_MAX = _OUT_MAX;
      OUT_MIN = _OUT_MIN ;
    }

    /************************** Deadzone   ***************************/

    void deadzone::set_attributes() { }
    void deadzone::initialize() {  }
    void deadzone::processing() {

      internal = 0.0;

      if (( in.read() >= low) & (  in.read() <= high  )) {
	internal = 0.0 ;
      }

      else if (  in.read()  <= low  ) {
	internal = in.read() - low ;

      }

      else if (  in.read() >= high  ) {
	internal = in.read() - high ;

      }

      else internal = 0.0;


      out.write(internal);



    }

    deadzone::deadzone(sc_core::sc_module_name nm, double _high, double _low) {

      high = _high;
      low = _low ;

    }

    /*******************************Coulomb*************************************/



    void coulomb::set_attributes() { }

    void coulomb::initialize() {  }

    void coulomb::processing() {

      if (in.read() >= 0.0)
	out.write(gain * in.read() + y  );
      else
	out.write(gain * in.read() - y  );
    }

    coulomb::coulomb(sc_core::sc_module_name nm, double _gain, double _y ) {

      gain = _gain;
      y = _y ;

    }



/************************** INOUT**************************/
    void inout::processing()
    {
      out.write(in.read());
    }

    inout::inout(sc_core::sc_module_name n){}




    /*************************  Offset ***************************************/


    void offset::set_attributes() {}

    void offset::initialize() {

    }


    void offset::processing() {

      out.write( ( in.read() ) + ofs ) ;

    }

    offset::offset(sc_core::sc_module_name nm, double _offset)  {

      ofs= _offset ;

    }

    /************************** Integrator ***************************/

    void integ::processing()
    {
      out.write(ltf_1(B, A, in.read()+ initial));	// Perform the Laplace-transformation on the absolute value of the input signal and add the initial value to it
    }

    integ::integ(sc_core::sc_module_name n, double _initial)
    {
      B(0)= 1.0; A(0)=0.0; 		// values for the LTF
      A(1)= 1.0;	                // to describe laplace transfer function
      initial =  _initial;
    }

    /**********************************Up Sampler**********************************/


    void upsample::set_attributes() {

      out.set_rate(rate);

    }

    void upsample::initialize() {  }

    void upsample::processing() {

      for (int i = 0; i < rate; i++) {
	out.write(in.read(), i);
      }

    }

    upsample::upsample(sc_core::sc_module_name nm, int _rate) {

      rate= _rate ;

    }

/***************************Down Sampler******************************************************/


    void downsample::set_attributes() {

      in.set_rate(rate);

    }

    void downsample::initialize() {}

    void downsample::processing() {

      out.write(in.read(sel-1))	;

    }

    downsample::downsample(sc_core::sc_module_name nm, int _sel, int _rate_in) {

      rate= _rate_in ;
      sel=_sel;
    }

 /*****************************************************Fork with one output delayed one sample*********************************************************/
    void o2t::set_attributes(){
  out.set_rate(2);
}

void o2t::processing()
{
  inp=in.read();
  for(int i=0;i<2;i++)
    out.write(inp,i);
}

o2t::o2t(sc_core::sc_module_name){}

/************************************************Reverse Fork***************************************************************/

void t2o::set_attributes(){
  in.set_rate(2);
}

void t2o::processing()
{
  inp=in.read();

  out.write(inp);
}

t2o::t2o(sc_core::sc_module_name){}

//*******************


/*****************************************Nonlinearity****************************************************/
void nonlinearity::processing(){
  // check if input is out of the input-range

  if( abs(in.read()) <  in_max) {
    //      out.write(noise()+a*in.read() - c*pow(in.read(),3));  TODO:noise of module
    out.write(a*in.read() - c*pow(in.read(),3));
  }
  else if( in.read() < 0 ){
    out.write(-out_max);
  } else {
    out.write(out_max);
  }
}

nonlinearity::nonlinearity(sc_core::sc_module_name n, double gain, double ip3){
  // Calculation of polynomial faktors

  a = pow(10,gain/20.0);
  c = a * (4.0/3.0)/(pow(10,ip3/10)*2*50*0.001);

  // calculate the max input and output of LNA
  in_max = sqrt(a/(3*c));
  out_max = (2*a/3)*in_max;

  //     v_sq_noise = a*a*4.0 * (pow(10,(nf/10))-1) * 1.38 * pow(10.0,-28.0) * 290.0 * 50.0;
  //     v_noise = sqrt(v_sq_noise);

  //Display the calculated model parameters

  cout<<"Non-linear model parameters are "<<n<<": \n";
  cout<<"a: "<<a<<"\n";
  cout<<"c: "<<c<<"\n";
  cout<<"In_max: "<<in_max<<"\n";
  cout<<"Out_max: "<<out_max<<"\n";
  //     cout<<"Rausch-Spg: "<<v_noise<<"V rms"<<"\n\n";

}

/******************Uniformly distributed random Numbers****************/

      //void noise_uniform::initialize(){srand((unsigned)time(NULL));}                  // remove srand() so that the simulation results can be reproduced

    void noise_uniform::initialize(){};
    void noise_uniform::processing() {

      double ampl = max - min ;
      for(int i=0;i<rate;i++)
	{
	  rnd=ampl * ((double)rand())/((double)RAND_MAX)+min;                     // calculate random numbers;
	  out.write(rnd,i);                                                       // output rnd
	}
    }

    void noise_uniform::set_attributes()  {
      out.set_rate(rate);                                                         //set the data rate of the output port
    }

    noise_uniform::noise_uniform(sc_core::sc_module_name nm, double _min, double _max, int _rate) {
      min = _min;
      max = _max;
      rate = _rate;
    }



    /******************Uniformly distributed random Bits****************/

    //void rand_bool::initialize(){srand((unsigned)time(NULL));}                        //remove srand() so that the simulation results can be reproduced

    void rand_bool::initialize(){};

    void rand_bool::set_attributes(){                                                 //set the data rate of the output port
      out.set_rate(rate);

    }

    void rand_bool::processing() {

      for(int i=0;i<rate;i++)
	{
	  rnd =   (double)rand() / (double)RAND_MAX ;
	  if ( rnd >= 0.5 )
	    out.write(true,i);

	  else
	    out.write(false,i);


    }}

    rand_bool::rand_bool(sc_core::sc_module_name nm,int _rate)//:rate(_rate)
    {
      rate=_rate;
    }


 /******************Gausssian distributed random numbers****************/

    void noise_gauss::processing() {

      double rnd1;                                                                // variables for calculation of random numbers
      double rnd2;

      double G1;

      double Q;
      double Q1;
      double Q2;

      for(int i=0;i<rate;i++)
	{
	  do{

	    rnd1 = ((double)rand()) / ((double)RAND_MAX) ;                        // calculate the gausssian distributed random numbers
	    rnd2 = ((double)rand()) / ((double)RAND_MAX) ;                        // using polar method

	    Q1 = 2.0 * rnd1 - 1.0 ;
	    Q2 = 2.0 * rnd2 - 1.0 ;

	    Q = Q1 * Q1 + Q2 * Q2 ;

	  } while (Q > 1.0) ;

	  G1 = sqrt( - 2.0 * log(Q) / Q) * Q1 ;
	  out.write( mean + sqrt(variance) * G1 * Q2,i) ;
	}
    }

    void noise_gauss::set_attributes()  {                                             // set data rate of output port
      out.set_rate(rate);
    }

    void noise_gauss::initialize(){}

    noise_gauss::noise_gauss(sc_core::sc_module_name nm, double  _variance , double _mean, int _rate) {
      variance = _variance ;
      mean = _mean ;
      rate = _rate;
    }




    /************************** Sine Generator ***************************/
    void sine::set_attributes(){
      out.set_rate(rate);                                     // set the data rate of output port
    }

    void sine::initialize(){
      sample_time = out.get_timestep().to_seconds();                // read time step of output port
    }

    void sine::processing() {

      if(freq_recon==true)
	freq_reg=freq_in->read();                            // read the frequency from mode port,and save the value into parameter register

      if(ampl_recon==true) 	                             // read the amplitude from amp_in port,and save the value into parameter register
	amp_reg=amp_in->read();

      stepsize = sample_time*freq_reg*2.*M_PI;               // calculate stepsize according to the time step of output port

      for(int i=0; i<rate; i++)  		             // write <rate> data token for one bit
	{
	  out.write(amp_reg*sin( actval + phi )+offset,i);
	  actval+=stepsize;
	}
    }

    sine::sine(sc_core::sc_module_name n, double freq_def,double amp_def,double _phi, double _offset, bool amp_c, bool freq_c,int datarate)
    {
      if((amp_c)&&(!freq_c))
	{
	  amp_in = new sca_tdf::sca_in<double>;                   // create amp_in port if only amp_c is set to true
	  freq_recon =false;
	}
      else if((freq_c)&&(!amp_c))                            // create freq_in port if only freq_c is set to true
	{
	  freq_in = new sca_tdf::sca_in<double>;
	  ampl_recon =false;
	}
      else if(amp_c&&freq_c)                                 // create freq_in and amp_in ports if both freq_c and amp_c are set to true
	{
	  amp_in = new sca_tdf::sca_in<double>;
	  freq_in = new sca_tdf::sca_in<double>;
	}
      else                                                  // no creation if both freq_c and amp_c are set to false
	{
	  freq_recon =false;
	  ampl_recon =false;
	}
      amp_reg  = amp_def;
      freq_reg = freq_def;
      rate = datarate;
      actval = 0.;
      phi=_phi;
      offset=_offset;
    }

    void sine::freq_con (sca_tdf::sca_signal<double>* f_in)                         // function to refer value of signal "f_in" to the optional port "freq_in"
    {
      (*freq_in)(*f_in);
      freq_recon = true;
    }
    void sine::amp_con (sca_tdf::sca_signal<double>* a_in)                          // function to refer value of signal "a_in" to the optional port "ampl_in"
    {
      (*amp_in)(*a_in);
      ampl_recon = true;
    }


 /************************** Square Wave Generator ***************************/
    void sqr_gen::set_attributes()
    {
      out.set_rate(rate);                                                      // set the data rate of the output port
    }

    void sqr_gen::initialize()
    {
      stepsize = out.get_timestep().to_seconds()*2*M_PI*f;		               // compute stepsize of output signal according to the time step of output port
    }

    void sqr_gen::processing()
    {
      for(int i=0; i<rate; i++)  	                	               // write <rate> data token for one bit
	{
	  if(sin( actval )>0)
	    out.write(ofs+amp,i);
	  else
	    out.write(ofs-amp,i);
	  actval+=stepsize;
	}
    }

    sqr_gen::sqr_gen(sc_core::sc_module_name n, double freq, double _amp, double _ofs, int d_rate)
    {
      ofs=_ofs;
      amp=_amp;
      rate=d_rate;
      f=freq;
      actval=0.0;
    }




 /************************** Triangle Wave Generator ***************************/

    void tri_gen::set_attributes()
    {
      out.set_rate(rate);                                                       // set the data rate of the output port
    }

    void tri_gen::initialize()
    {
      stepsize = out.get_timestep().to_seconds()*2.0*M_PI*f;	                	// compute stepsize of output signal according to the time step of output port
    }

    void tri_gen:: processing()
    {
      for(int i=0; i<rate; i++)  	                                	// write <rate> data token for one bit
	{
	  out.write(ofs+2*amp/M_PI*asin(sin(actval)),i);
	  actval+=stepsize;
	}
    }

    tri_gen::tri_gen(sc_core::sc_module_name n, double freq, double _amp, double _ofs, int d_rate)
    {
      f=freq;
      ofs=_ofs;
      amp=_amp;
      rate=d_rate;
    }

    /************************** Sawtooth Wave Generator ***************************/

    void saw_gen::set_attributes()
    {
      out.set_rate(rate);                                                      // set the data rate of the output port
    }

    void saw_gen::initialize()
    {
      stepsize = out.get_timestep().to_seconds();	         	               // stepsize of output signal according to the time step of output port
    }

    void saw_gen::processing()
    {
      for(int i=0; i<rate; i++)  		                               // write <rate> data token for one bit
	{
#ifdef _MSC_VER
	  out.write(ofs+amp*2*(actval*f-floor((double)(actval*f+0.5))),i);
#else
		out.write(ofs+amp*2*(actval*f-floor(actval*f+0.5)),i);
#endif
	  actval+=stepsize;
	}
    }

    saw_gen::saw_gen(sc_core::sc_module_name n, double freq,double _amp, double _ofs, int d_rate)
    {
      ofs=_ofs;
      f=freq;
      amp=_amp;
      rate=d_rate;
      actval=0.0;
    }



/*************************************************wave generator********************************************/


wave::wave(sc_core::sc_module_name n, string type, double freq ,double amp, double ofs, int rate)
{
  w_type=type;
  if(type=="square")
    {
      sqr_sub = new sqr_gen("sqr_sub", freq, amp, ofs , rate);
      sqr_sub->out(out);
    }
  else if(type=="triangle")
    {
      tri_sub = new tri_gen("tri_sub", freq, amp, ofs , rate);
      tri_sub->out(out);
    }
  else if(type=="saw_tooth")
    {
      saw_sub = new saw_gen("tri_sub", freq, amp, ofs , rate);
      saw_sub->out(out);
    }
  else
    {
      cout<<"Error: Unsupported wave type!"<<endl;
      exit(0);
    }
}


    void wave::set_time_step(double sp ,char unit)
    {
      sc_time step(sp,unit);
      if(w_type=="square")
	{
	  sqr_sub->out.set_timestep(step);
	}
      else if(w_type=="triangle")
    {
      tri_sub->out.set_timestep(step);
    }
      else if(w_type=="saw_tooth")
	{
	  saw_sub->out.set_timestep(step);
	}
    }

  /*****************************************1. orderLowpass filter*************************/
void lp::initialize()
{
  out_tmp=0.0;
}

void lp::processing()
{
  // out_tmp=ltf_1(B, A, in.read());
  out.write(ltf_1(B, A, in.read()));	// Perform the Laplace-transformation on the absolute value of the carrier signal
}

lp::lp(sc_core::sc_module_name n, double freq_cut)
{
  B(0)= 1.0; A(0)=1.0; 		// values for the LTF
  A(1)= 1.0/(2.0*M_PI*freq_cut);	// to describe a lowpass-filter
}






 /**********************************Butter worth**********************************/

#if  ( (!defined(_MSC_VER)) && (!defined(__APPLE__)) && (!defined(__MACH__)) )

void btworth::initialize(){

  //Determination of the order
  int n = order();

  //Calculation of wc (half-power frequency)
  wc = wp/pow((exp10(-gp/10.)-1), (n/2));

  //Calculation of poles
  calc_poles(n);

  //Sort poles
  sort_poles(n);

  //Multiply all conjugated complex poles
  make_real(n);

  //Calculation of the characteristic polynomial
  calc_char_poly(n);
}

    void btworth::processing()
{
  out.write(ltf(B,A,in));
}

int btworth::order(){
  return (int)ceil(log10(( (exp10(-gs/10.)-1) / (exp10(-gp/10.)-1) )) / (2.*log10(ws/wp)));
}

void btworth::calc_poles(int n){

  //position in poles_unsorted
  int p = 0;

  for(int k = 1; k <= n; k+=1){

    poles_unsorted(p)= 1;					//s
    poles_unsorted(p+1) = -cos((M_PI/(2.*n))*(2.*k+n-1.));	//Re
    poles_unsorted(p+2) = -sin((M_PI/(2.*n))*(2.*k+n-1.));	//Im

    //if n odd --> middle pole is real
    if((n%2 == 1)&&(k == (n/2)+1)){

      poles_unsorted(p+2) = 0;			//Im = 0
    }

    p+=3;
  }
}

void btworth::sort_poles(int n){

  int dummy_odd = 0;
  int dummy_even = 0;
  int p = 0;

  for(int k = 0; k < n ; k+=1){

    if(k%2 == 1){
      poles_sorted(p) = poles_unsorted((n*3)-p+dummy_odd*3);
      poles_sorted(p+1) = poles_unsorted((n*3)-(p-1)+dummy_odd*3);
      poles_sorted(p+2) = poles_unsorted((n*3)-(p-2)+dummy_odd*3);

      dummy_odd+=1;

    }else{
      poles_sorted(p) = poles_unsorted(p-dummy_even*3);
      poles_sorted(p+1) = poles_unsorted(p+1-dummy_even*3);
      poles_sorted(p+2) = poles_unsorted(p+2-dummy_even*3);

      dummy_even+=1;
    }

    p += 3;
  }
}

void btworth::make_real(int n){

  //position in poles_sorted
  int p = 0;
  //position in poly_coef
  int q = 0;

  for(int k = 0; k < n; k+=1){

    sca_vector<double> result;
    sca_vector<double> poly_dummy_sorted;

    if(k%2 == 1){

      poly_dummy_sorted(0) = poles_sorted(p);
      poly_dummy_sorted(1) = poles_sorted(p+1);
      poly_dummy_sorted(2) = poles_sorted(p+2);

      poles_sorted(0) = poles_sorted(p+3);
      poles_sorted(1) = poles_sorted(p+4);
      poles_sorted(2) = poles_sorted(p+5);

      for(int i = 0; i < 3; i+=1){

	for(int j = 0; j < 3; j+=1){

	  result(i+j) = result(i+j) + poly_dummy_sorted(i)*poles_sorted(j);
	}

      }

      //Summation of all coefficients without s
      for(int x = 2; x <= 4; x+=1){

	result(2) = result(2) + abs(result(x+1));
      }

      //Store coefficients in factorized form
      for(int a = 0; a< 3; a+=1){

	poly_coef(a+q) = result(a);
      }

      p += 6;
      q += 3;
    }
  }
  if(n%2 == 1){

    poly_coef((n/2)*3) = poles_sorted(p);
    poly_coef((n/2)*3+1) = poles_sorted(p+1);
    poly_coef((n/2)*3+2) = poles_sorted(p+2);
  }
}

    void btworth::calc_char_poly(int n){

      int dummy_order = 3;
      int p = 3;

      poly_dummy(0) = poly_coef(0);
      poly_dummy(1) = poly_coef(1);
      poly_dummy(2) = poly_coef(2);

      while(poly_coef(p)){

	dummy_order = poly_multi(dummy_order, p);

	p += 3;
      }

      //Distinction between highpass and lowpass
      if(type == "lowpass"){

	for(int i = 0; i <= n; i+=1){

	  A(n-i) = poly_dummy(i)*pow((1/(wc)),(n-i));
	}

	B(0) = 1.;

      }else if(type == "highpass"){

	for(int i = 0; i <= n; i+=1){

	  A(i) = poly_dummy(i)*pow((1/(wc)),i);
	}

	B(n) = pow((1/wc),n);
      }else{

	cerr << "ERROR --> Filter must be highpass or lowpass!!!" << endl;
	exit(-1);
      }
    }

    int btworth::poly_multi(int m, int n){

    int i = 0;
    int j = 0;
    int k = 0;

    sca_vector<double> result;

    for(i = 0; i < 3; i+=1){

      for(j = 0; j < m; j+=1){

	result(i+j) = result(i+j) + poly_dummy(j)*poly_coef(i+n);

      }
    }

    for(k = 0; k <= (i+j-2); k+=1){

      poly_dummy(k) = result(k);

    }

    return i+j-1;
  }

    btworth::btworth(sc_core::sc_module_name nm, string _type, double _gp, double _gs, double _fp, double _fs)
    {
      type = _type;
      gp = _gp;
      gs = _gs;
      wp = 2*M_PI*_fp;
      ws = 2*M_PI*_fs;
    }


#endif


/*******************************chebyshev*************************************/


void chebyshev::processing(){

  if(type=="lowpass")                                                                 // determine the prototype lowpass filter
    {
      wp=wp_tmp;
      ws=ws_tmp;
    }
  else if(type=="highpass")                                                           // pre handel of characaters for frequency transformations
    {
      wp=1;
      ws=wp_tmp/ws_tmp;
    }

  double e;                                                 //height of ripples
  double x;
  double real;                                              //variable to save real part of a root
  double imag;                                              //variable to save imag part of a root

  back_door=true;

 vector<double> data;                                      //vector to save roots of filter
 vector<double> data2;                                     //tmp vector for calculating
  int data_size;                                            //size of vector data
  int data2_size;                                           //size of vector data2

vector<double> s;                                         //vector to save  coefficients of LPF
 vector<double> p;                                         //tmp vector for calculating
 vector<double> q;                                         //tmp vector for calculating

                                                            // calculate the order of normalized lowpass filter
  n=static_cast<int>(ceil(1.0/acosh(ws/wp)*acosh(pow((pow(10.,-gs/10)-1)/(pow(10.,r/10)-1),1./2.))));

  e=sqrt(pow(10,r/10)-1);                                                               // calculate the height of ripples

  x=1./n*asinh(1./e);

  if (back_door==true)
    {
      cout<< "Order of the filter is" << n <<endl;
      cout<< "E of the filter is" << e <<endl;
      cout<< "x of the filter is" << x <<endl;
    }

  for (int i=1; i<=n; i++)                                     //calculate roots of the "n" order normalized lowpass filter
    {
      real=-sin(((2.*i-1.)*M_PI)/(2*n))*sinh(x);
      imag=cos(((2.*i-1.)*M_PI)/(2*n))*cosh(x);

      if (back_door==true)
	{

	  cout<<"i is " << i <<endl;
	  cout<<"tmp_real is " << real <<endl;
	  cout<<"tmp_imag is " << imag <<endl;
	}
      data.push_back(real);
      data.push_back(imag);
    }

  if (back_door==true)
    {
      data_size=data.size();
      cout<< "data_size is " << data_size <<endl;
      for(int j=0;j<data_size;j++)
	cout<< data[j] <<endl;
    }

  data_size = data.size();

  if (data_size%4==0)                                          // check if the number of roots is even
    {
      for (int i=0;i<data_size/2;i+=2)
	{
	  data2.push_back(1.0);
	  data2.push_back(-2*data[i]);
	  data2.push_back(data[i]*data[i]+data[i+1]*data[i+1]);
	}
    }
  else
    {
      if(data_size==2)
	{
	  data2.push_back(0.);
	  data2.push_back(0.);
	  data2.push_back(-data[0]);
	}
      else
	{
	  for (int i=0;i<data_size/2;i+=2)
	    {
#ifdef _MSC_VER
		  if(i!=(int)ceil(((float)data_size)/2-1))	      
#else
		  if(i!=ceil(data_size/2-1))
#endif
		{
		  data2.push_back(1.0);
		  data2.push_back(-2*data[i]);
		  data2.push_back(data[i]*data[i]+data[i+1]*data[i+1]);
		}
	      else
		{
		  data2.push_back(0.);
		  data2.push_back(1.);
		  data2.push_back(-data[data_size/2-1]);
		}
	    }
	}
    }

  data2_size=data2.size();
  if (back_door==true)
    {
      data2_size=data2.size();
      cout<< "data2 is" <<endl;
      for(int j=0;j<data2_size;j++)
	cout<< data2[j] <<endl;
    }

  if(n==1)
    {
      for (int i=0; i<n+2; i++) s.push_back(0.0);
      s[2]=-data[0];
      s[1]=1;
      s[0]=0;
    }
  else if(n==2)
    {
      for (int i=0; i<n+2; i++) s.push_back(0.0);
      for(int i=0;i<3;i++)
	{
	  s[i]=data2[i];
	}
    }
  else if(n <1)
    cout<< "Error: Filter order must be a integer lager than 0. It is impossible to realize your specification using Chebyscev filter"<<endl;

  else
    {
      int k,m,l;
      m=3;
      l=3;
      k=6;

      for (int i=0; i<n+2; i++) s.push_back(0.0);       //initialize array s;

      if (back_door==true)
	{
	  cout<< "Size of s is " << s.size()<<endl;
	  cout<< "s is " <<endl;
	  for(int j=0;j<n+2;j++)
	    cout<< s[j] <<endl;
	}

      for (int i=0;i<m;i++)                            //initialize vector p and q;
	{
	  p.push_back(data2[i]);
	}

      for (int i=0;i<l;i++)
	{
	  q.push_back(data2[i+3]);
	}

      if (back_door==true)
	{
	  int p_size=p.size();
	  cout<< "p is " <<endl;
	  for(int j=0;j<p_size;j++)
	    cout <<p[j] <<endl;

	  int q_size=q.size();
	  cout<<"q is "<<endl;
	  for(int j=0;j<q_size;j++)
	    cout<< q[j] <<endl;
	}

      if (back_door==true)
	{
	  cout<<"data2_size is " <<endl;
	  cout<< data2.size() <<endl;
	}

      do                                                //calculate normalized coefficients of chebyshev polynomial
	{
	  for (int i=0; i<m+l-1; i++) s[i]=0.0;
	  for (int i=0; i<m; i++)
	    for (int j=0;j<l;j++)
	      s[i+j]=s[i+j]+p[i]*q[j];

	  p.clear();
	  q.clear();

	  if(k==data2_size) break;

	  for (int i=0;i<m+l-1;i++)
	    {
	      p.push_back(s[i]);
	    }

	  int p_size=p.size();
	  cout<< "p is now " <<endl;
	  for(int j=0;j<p_size;j++)
	    cout <<p[j] <<endl;

	  for (int i=0;i<l;i++)
	    {
	      q.push_back(data2[k]);
	      k++;
	    }

	  int q_size=q.size();
	  cout<<"q is now "<<endl;
	  for(int j=0;j<q_size;j++)
	    cout<< q[j] <<endl;
	  cout<< "k is now " << k  <<endl;

	  m=m+2;
	} while(k<=data2_size);
    }

  if (back_door==true)
    {
      cout<< "s at the end is " <<endl;
      for(int j=0;j<n+2;j++)
	cout<< s[j] <<endl;
    }	                                                                //prepair coefficients for LTF

  if (n%2==0)                                                                //when oder is even
    {
      for (int i=0; i<n+1;i++)
	A(n-i)=s[i]/pow(wp,n-i);


      B(0)=1.0*s[n]/(pow(10,r/20)*A(n));
      for (int i=0; i<n+1;i++)
	A(i)=A(i)/A(n);
    }
  else                                                                       //when oder is odd
    {
      for (int i=0; i<n+2;i++)
	{
	  A(n-i)=s[i+1]/pow(wp,n-i);
	}

      B(0)=s[n+1]/A(n);

      for (int i=0; i<n+1;i++)
	A(i)=A(i)/A(n);
    }

  if (back_door==true)
    {
      cout<< "B(0) at the end is " << B(0) <<endl;
      cout<< "A at the end is " <<endl;
      for(int j=0;j<n+1;j++)
	cout<< A(j) <<endl;
    }

  if(type=="highpass")                                                      //frequency transformations(lowpass filter -> highpass filter)
    {
      double tmp=A(0);
      B(n)=B(0)/tmp;
      for(int i=0;i<n;i++)
	B(i)=0;
      for(int i=0;i<n+1;i++)
	A(i)=A(i)*pow(wp_tmp,i)/tmp;
#ifdef _MSC_VER
	for(int i=0;i<floor((float)((n+1)/2));i++)    
#else
	for(int i=0;i<floor((n+1)/2);i++)  
#endif
	  swap(A(i),A(n-i));

    }

  if (back_door==true)
    {
      cout<< "A after swap is " <<endl;
      for(int j=0;j<n+1;j++)
	cout<< A(j) <<endl;

      cout<< "B after swap is " <<endl;
      for(int j=0;j<n+1;j++)
	cout<< B(j) <<endl;
    }

  out.write(ltf_1(B, A, in.read()));

}

chebyshev::chebyshev(sc_core::sc_module_name n, string _type, double _ratio, double _gs, double _wp, double _ws)
    {
      type=_type;
      r =_ratio;
      gs =_gs;
      wp_tmp=_wp;
      ws_tmp=_ws;
    }
/************************************Gaussian puls shaping**********************************************/
                      /*********in fact bessel lowpass filter***********/

void gauss_shaping::initialize()
{
  freq_cut=bts/ts*2*M_PI;
  B(0)= 2027025.0; A(0)=2027025.0; 		// values for the LTF
  A(1)= 2027025.0/(freq_cut);          	// to describe a lowpass-filter
  A(2)= 945945.0/pow(freq_cut,2);
  A(3)= 270270.0/pow(freq_cut,3);
  A(4)= 51975.0/pow(freq_cut,4);
  A(5)= 6930.0/pow(freq_cut,5);
  A(6)= 630.0/pow(freq_cut,6);
  A(7)= 36.0/pow(freq_cut,7);
  A(8)= 1.0/pow(freq_cut,8);
}

void gauss_shaping::processing()
{
  out.write(ltf_1(B, A, in.read()));	// Perform the Laplace-transformation on the absolute value of the carrier signal
}

gauss_shaping::gauss_shaping(sc_core::sc_module_name n, double _bts, double _ts)
{
  bts=_bts;
  ts=_ts;
}
 /********************************Gain Controller************************************/




    void gain::set_attributes() {in.set_rate(rate);out.set_rate(rate);}

    void gain::initialize() {  }


    void gain::processing() {
      for(int i=0;i<rate;i++)
	out.write( ( g *  in.read(i) ) + offset,i ) ;

    }

    gain::gain(sc_core::sc_module_name nm, double _gain, double _offset, int data_rate) {
      rate=data_rate;
      g = _gain ;
      offset = _offset ;

    }


/********************************Mixer**************************************/

void mixer::processing() {}
mixer::mixer(sc_core::sc_module_name n, double _gain, double ip3, bool _ideal)
{
  ideal_mixer = new mul("ideal_mixer_");

  ideal=_ideal;
  if(ideal)
    {
      i_gain = new gain("gain_sub",_gain);
      i_gain->in(sig_in);
      i_gain->out(sig_gout);
      ideal_mixer->in1(sig_gout);
    }
  else
    {
      mixer_nonlinearity = new nonlinearity("nonlinear_behavior", _gain, ip3);
      mixer_nonlinearity->in(sig_in);
      mixer_nonlinearity->out(sig_nonlinear);
      ideal_mixer->in1(sig_nonlinear);
    }

  ideal_mixer->in2(lo_in);
  ideal_mixer->out(out);
}

 /**********************************phase detector*******************************************************/

void phc::set_attributes()
{
  out.set_rate(rate);
}


void phc::processing()
{

  for(int i=0; i<rate; i++)  		// write <rate> data token for one bit
    {
      out.write(gain*in_ref.read()*in_vco.read(),i);
    }
}

phc::phc(sc_core::sc_module_name n, int datarate, double _gain)
  {
    rate=datarate;
    gain=_gain;
  }


/****************************************peak detector**************************************************/

void peak_dt::initialize()
  {
    peak=0.0;
  }

void peak_dt::processing()
  {
    p_dt(in.read());
    out.write(peak);                    //detect the peak of input signal
  }

peak_dt::peak_dt(sc_core::sc_module_name n){};

double peak_dt::p_dt(double sig){      //detect the anplitude of input signal
  if(sig>peak)
    {
      peak=sig;
    }
  return 0;
}

/********************************* -PI/4 Shifter****************************************/

nshift_45::nshift_45(sc_core::sc_module_name n, double w_if,double delta_R, double delta_C)
{
  r=1;
  c=1/w_if;

  vin = new sca_eln::sca_tdf::sca_vsource("vin");  // we pin the SDF-controlled
  vin->p(n1);		  // voltage source to n1
  vin->n(gnd);		  // and gnd
  vin->inp(in);

  r1 = new sca_eln::sca_r("r1");	// r1 connects n1 and n2
  r1->value=r+delta_R;
  r1->p(n2);
  r1->n(gnd);

  c1 = new sca_eln::sca_c("c1");	// and c1 connects n2 with gnd
  c1->value=c+delta_C;
  c1->p(n1);
  c1->n(n2);

  vout = new sca_eln::sca_tdf::sca_vsink("vout");  // this one transforms the
  vout->p(n2);  		    // between n2 and gnd
  vout->n(gnd);
  vout->outp(out);      // to a SDF-signal
}

/********************************* +PI/4 Shifter****************************************/

pshift_45::pshift_45(sc_core::sc_module_name n, double w_if,double delta_R, double delta_C)
{

  r=1;
  c=1/w_if;

  vin = new sca_eln::sca_tdf::sca_vsource("vin");  // we pin the SDF-controlled
  vin->p(n1);		  // voltage source to n1
  vin->n(gnd);		  // and gnd
  vin->inp(in);

  r1 = new sca_eln::sca_r("r1");	// r1 connects n1 and n2
  r1->value=r+delta_R;
  r1->p(n1);
  r1->n(n2);

  c1 = new sca_eln::sca_c("c1");	// and c1 connects n2 with gnd
  c1->value=c+delta_C;
  c1->p(n2);
  c1->n(gnd);

  vout = new sca_eln::sca_tdf::sca_vsink("vout");  // this one transforms the
  vout->p(n2);
  vout->n(gnd);	    // between n2 and gnd
  vout->outp(out);      // to a SDF-signal
}



/****** **********************************shifter ***********************************************/



shift_45::shift_45(sc_core::sc_module_name n,string mode, double w_if, double delta_r, double delta_c)
{
  i_gain = new gain("i_gain",sqrt(2.0));
  i_gain->in(in);
  i_gain->out(gain_out);
  if(mode=="pos")
    {
      i_pshift = new pshift_45("i_pshift",w_if,delta_r,delta_c);
      i_pshift->in(gain_out);
      i_pshift->out(out);
    }
  else if(mode=="neg")
    {
	i_nshift = new nshift_45("i_nshift",w_if,delta_r,delta_c);
	i_nshift->in(gain_out);
	i_nshift->out(out);
    }
  else
    cout<<"error: mode can either be pos or neg"<<endl;
}

    /**************************************digital VCO**************************************************/

void d_vco::initialize()
{
  sample_time = out.get_timestep().to_seconds();	// read sample time
  cout << "Sample-Time: " << sample_time << endl;
}

void d_vco::set_attributes()
{
  out.set_rate(rate);
}

void d_vco::processing() {

    double fr=freq+kvco*in.read();
    samples=(int)ceil(1/(fr*sample_time));
    for(int i=0; i<rate; i++)  		// write <rate> data token for one bit
    {
      out.write(actval,i);
      actsamples++;
      if(actsamples <= samples/2)
	{
	  actval=-1.0*gain;
	}
      else if((samples/2 < actsamples) & (actsamples < samples))
      	{
	  actval=1.0*gain;
      	}
      else
	{
	actsamples=0;
	actval=-1.0*gain;
	}
    }
  }

d_vco::d_vco(sc_core::sc_module_name n, double freq_data, int datarate, double _kvco, double _gain)
  {
    freq=freq_data;
    rate=datarate;
    actval = -1.0;
    actsamples=0;
    kvco=_kvco;
    gain=_gain;
  }

/******************************analog VCO*************************************************/
void a_vco::initialize()
{
  //sample_time = out.get_timestep().to_seconds();	// read sample time
  sample_time = out.get_timestep().to_seconds();
  cout << "Sample-Time: " << sample_time << endl;
}

void a_vco::set_attributes()
{
  out.set_rate(rate);
}

void a_vco::processing()
{

  double fr=freq+kvco*in.read();
  stepsize = sample_time*fr*2.*M_PI;

  for(int i=0; i<rate; i++)  		// write <rate> data token
    {
      out.write(gain*sin(actval),i);
      actval+=stepsize;
    }
}

a_vco::a_vco(sc_core::sc_module_name n, double freq_data, int datarate, double _kvco, double _gain)
{
  freq=freq_data;
  rate=datarate;
  actval = 0.0;
  kvco=_kvco;
  gain=_gain;
}
/*****************************PLL***************************************/

    pll::pll(sc_core::sc_module_name n,double phc_gain,double lp_fc,double vco_freq,double kvco,double vco_gain,int vco_out_rate,int phc_out_rate)
    {
      phc_sub = new phc("phc_sub",phc_out_rate,phc_gain);
      phc_sub->in_ref(ref);
      phc_sub->in_vco(vcoo);
      phc_sub->out(pco);

      lp_sub = new lp("lp_sub",lp_fc);
      lp_sub->in(pco);
      lp_sub->out(lpo);

      vco_sub = new a_vco("vco_sub",vco_freq,vco_out_rate,kvco,vco_gain);
      vco_sub->in(lpo);
      vco_sub->out(vcoo);
      vco_sub->out.set_delay(1);
    }

/***********************************************LNA**********************************************/


// Konstruktor
lna::lna(sc_core::sc_module_name n, double _gain, double _ip3, bool _ideal){

   // Parameterwerte zuweisen
   gain = _gain;
   ip3 = _ip3;
   ideal = _ideal;

   // Faktoren berechnen
   a = pow(10,gain/20);                             // Verstaerkungsfaktor
   b = 0.0;
   c = (a*4)/(3*(pow(10,ip3/10)*2*(50.0/1000.0)));

   // Kritische Werte berechnen
   input_max = sqrt(a/(3*c));
   output_max = ((2*a)/3)*input_max;

   #if _LNA_PROT
   char filename[25];
   strncpy(filename, n,14);
   strcat(filename, ".prot");
   prot.open(filename, ios_base::trunc);
   if(prot){
      prot << "\t\tProtokollierung fuer LNA" << endl;
      prot << "\t\t------------------------" << endl << endl;
      prot << "Paramter:" << endl << endl;
      prot << "\tgain:  " << gain << endl;
      prot << "\tip3:   " << ip3 << endl;
      prot << "\tideal: " << ideal << endl;
      prot << "\ta:     " << a << endl;
      prot << "\tb:     " << b << endl;
      prot << "\tc:     " << c << endl;
      prot << "\tinmax: " << input_max << endl;
      prot << "\toutmax:" << output_max << endl << endl;
      prot << "Daten: " << endl;
      prot << "-----------------------------------------------------" << endl << endl;
   }
   #endif
}

// Destruktor
lna::~lna(){
   #if _LNA_PROT
   // Protokollierungsfile schliessen wenn geoffnet
   if(prot){
      prot.close();
   }
   #endif
}

// Funktion wird jeden Simulationscyclus aufgerufen
void lna::processing(){
   double temp_out = 0.0;

   // Ueberpruefen ob idealer LNA
   if(ideal){
      // Ausgangswert berechnen
      temp_out = a*in.read();

   } else {  // realer LNA

      // Ueberpruefen ob Input zu gross
      if(abs(in.read()) < input_max){
	 // Input ok
	 temp_out = a*in.read()-b*pow(in.read(),2)-c*pow(in.read(),3);
      } else {
	 // Input zu gross
	 // Ueberpruefen ob negativer Eingang
	 if(in.read() < 0){
	    // Ausgang auf negativen Ausgangsmaximalwert setzen
	    temp_out = -output_max;
	 } else {
	    // Ausgang auf Ausgangsmaximalwert setzen
	    temp_out = output_max;
	 }
      }
   }

   // Am Ausgang schreiben
   out.write(temp_out);

   // Ausgabe
   //cout << "lna_in: " << in.read() << "    lna_out: " << temp_out << endl;

   #if _LNA_PROT
   if(prot){
      prot << "time: " << setw(10) << sca_get_time().to_seconds();
      prot << "  in: " << setw(10) << in.read();
      prot << "  out: " << setw(10) << temp_out << endl;
   }
   #endif
}

    /**********************************air**********************************************************/

void air::set_attributes()
{
  in.set_rate(rate);
  out.set_rate(rate);
}

void air::initialize()
{
  srand((unsigned)time(NULL));
}

void air::processing()
{
  if(mode=="uniform")
    {
      double ampl = max - min ;
      double rnd;
      for(int i=0;i<rate;i++)
	{
	  rnd=( ampl * ((double)rand()) / ((double)RAND_MAX) + min ) ;
	  out.write(gain*in.read(i)+rnd,i);
	}
    }
  else if(mode=="gauss_white")
    {
      double rnd1;
      double rnd2;

      double G1;

      double Q;
      double Q1;
      double Q2;

      for(int i=0;i<rate;i++)
	{
	  do{

	    rnd1 = ((double)rand()) / ((double)RAND_MAX) ;
	    rnd2 = ((double)rand()) / ((double)RAND_MAX) ;

	    Q1 = 2.0 * rnd1 - 1.0 ;
	    Q2 = 2.0 * rnd2 - 1.0 ;

	    Q = Q1 * Q1 + Q2 * Q2 ;

	  } while (Q > 1.0) ;

	  G1 = sqrt( - 2.0 * log(Q) / Q) * Q1 ;
	  out.write(gain*in.read(i) + mean + sqrt(variance) * G1,i) ;

	}
    }
}

air::air(sc_core::sc_module_name n, double atten, string n_art,double a,double b, int d_rate)
{
  mode=n_art;
  rate=d_rate;
  gain=atten;
  if(n_art=="uniform")
    {min=a;max=b;mean=0;variance=0;}
  else if(n_art=="gauss_white")
    {mean=b;variance=a;min=0;max=0;}
  else
    cout<<"Errot:Undefined noise type!This can only be uniform or gauss"<<endl;
}

/************************** Differential encoder ***************************/

    d_coding::d_coding(sc_core::sc_module_name n)
    {
      xnor_sub = new xnor("i_xnor");
      xnor_sub->in1(in);
      xnor_sub->in2(sig_out);
      xnor_sub->out(sig_xnor);

      d_fp_sub = new d_flipflop("i_d_fp");
      d_fp_sub->in(sig_xnor);
      d_fp_sub->out(sig_out);

      inout_sub = new inout("i_inout");
      inout_sub->in(sig_out);
      inout_sub->out(out);
    }
/************************** Differential decoder ***************************/
    d_decoding::d_decoding(sc_core::sc_module_name n)
    {
      xnor_sub = new xnor("i_xnor");
      xnor_sub->in1(in);
      xnor_sub->in2(sig_dff);
      xnor_sub->out(out);

      d_fp_sub = new d_flipflop("i_d_fp");
      d_fp_sub->in(in);
      d_fp_sub->out(sig_dff);
    }


/************************** AM(Amplitude Modulation) Modulator ***************************/
    void mod_am::set_attributes()
    {
      out.set_rate(rate);
    }

    void mod_am::initialize(){
      delta_T = 1.0/((double)rate * freq);                                                                        // calculate the time step of RF carrier
    }

    void mod_am::processing() {
      for (int i=0; i<rate;i++)
	{
	  out.write((ampl+(in.read()+offset))*cos(2.0*M_PI*freq*(get_time().to_seconds()+delta_T*i)),i);      // input value multi carrier and output it
	}
    }

    mod_am::mod_am(sc_core::sc_module_name nm, double _freq , double  _ampl, double _offset , int _rate) {
      freq = _freq ;
      ampl = _ampl;
      offset = _offset;
      rate = _rate;
    }


/************************** BPSK Modulator ***************************/

    void bpsk::initialize()
    {
      sample_time = out.get_timestep().to_seconds();    // read sample time
      step_size=2.*M_PI*sample_time*freq;        // set the time step of sine wave
    }

    void bpsk::set_attributes()
    {
      out.set_rate(out_rate);
    };

    void bpsk::processing()
    {
      inp=in.read();

      for(int i=0; i<out_rate; i++)  		               // write <rate> data token for one bit
	{
	  if(inp)
	    out.write(sin(actval),i);                          // output sin(wt) if input is true
	  else
	    out.write(-1.0*sin(actval),i);                     // else  -sin(wt)
	  actval+=step_size;
	}
    }

    bpsk::bpsk(sc_core::sc_module_name n, double _freq, int _out_rate)
    {
      freq=_freq;
      out_rate=_out_rate;
      actval=0;
    }

/************************** BPSK Demodulator ***************************/
    void bpsk_de::set_attributes(){in.set_rate(in_rate);}
void bpsk_de::processing(){}
    bpsk_de::bpsk_de(sc_core::sc_module_name n, double _freq,int _in_rate)
    {
      sine_sub = new sine("i_sine",_freq,1.0,0.0,0.0,false,false,1);
      sine_sub->out(sig_sine);

      mixer_sub = new mixer("i_mixer",1.0,0.0,true);
      mixer_sub->sig_in(in);
      mixer_sub->lo_in(sig_sine);
      mixer_sub->out(sig_mix);

      lp_sub = new lp("i_lp",_freq/10);
      lp_sub->in(sig_mix);
      lp_sub->out(sig_lp);

      decider = new downsample("i_decider",_in_rate-1,_in_rate);
      decider->in(sig_lp);
      decider->out(sig_dd);

      comp_sub = new compare("i_comp");
      comp_sub->in(sig_dd);
      comp_sub->out(out);

      in_rate=_in_rate;
    }

 /************************** BASK Modulator ***************************/

    void mod_bask::set_attributes() {
      out.set_rate(rate);
    }

    void mod_bask::initialize() {}
    void mod_bask::processing() {

      delta_T = out.get_time().to_seconds(); // / rate ;
      // besser waere in before_end of_elaboration(), gibt aber speicherzugriffsfehler :-(


      for (int i = 0; i < rate; i++)
	{
	  carrier = cos((2.0 * M_PI * freq * (get_time().to_seconds() + delta_T * i) ) + phi );
	  if (in.read()){
	    out.write(ampl1 * carrier, i);}                                             // output ampl1 * carrier if input==1
	  else {
	    out.write(ampl0 * carrier, i);}	                                        // else output ampl0 * carrier
	}
    }

    mod_bask::mod_bask(sc_core::sc_module_name nm, double _freq, double _ampl1, double _ampl0, double _phi, int _rate) {
      freq = _freq;
      ampl1 = _ampl1;
      ampl0 = _ampl0;
      phi = _phi;
      rate = _rate;
    }
/**************************BASK deModulator******************************************/

   void demod_bask::set_attributes() {

   	in.set_rate(rate);

    }

   void demod_bask::initialize() {

     	B(0) = 1.0 ;
     	A(0) = 1 ;
     	A(1)= 1.0/(2.0*M_PI*freq ); }

   void demod_bask::processing() {

	  internal2= 0.0;//?

	for (int i = 0; i < rate; i++) {

	  internal = ltf(B, A, S, sqrt(in.read(i)*in.read(i)));

	  internal2 = internal2 + internal;

	}

	  internal2 = internal2 / rate;

	if (internal2 >= level)
		out.write(true);
	else
		out.write(false);


	   }

   demod_bask::demod_bask(sc_core::sc_module_name nm, double _level, double _freq, int _rate) {

   	level = _level;
	freq = _freq ;
	rate = _rate;

   }

/******quadrature mixer for receiver ********/

void q_mixer_re::initialize()
{
  sample_time = i_out.get_timestep().to_seconds();            // read sample time
}

void q_mixer_re::set_attributes()
{
  in.set_rate(rate);
  i_out.set_rate(rate);
  q_out.set_rate(rate);
}

void q_mixer_re::processing()
{
  if(freq_recon==true)
    freq=freq_in->read();

  stepsize = sample_time*freq*2.*M_PI;
  for(int i=0;i<rate;i++)
    {
      i_out.write(in.read(i)*2*(1+e/2)*sin( actval +theta/2)*amp,i);
      q_out.write(in.read(i)*2*(1-e/2)*cos( actval -theta/2)*amp,i);
      actval+=stepsize;
    }
}
    q_mixer_re::q_mixer_re(sc_core::sc_module_name n,double _freq, double _amp,int data_rate,bool f_config,double _theta,double _e)
{
  freq_recon = false;
  if(f_config)
    {
      freq_in = new sca_tdf::sca_in<double>;
      freq_recon = true;
    }
  amp=_amp;
  actval = 0.;
  freq=_freq;
  rate=data_rate;
  e=_e;
  theta=_theta;
}

void q_mixer_re::freq_con (sca_tdf::sca_signal<double>* f_in)
{
  (*freq_in)(*f_in);
  freq_recon = true;
}



/******quadrature mixer for transmitter ********/
void q_mixer_tr::initialize()
{
  sample_time = out.get_timestep().to_seconds();            // read sample time
}

void q_mixer_tr::set_attributes()
{
  out.set_rate(rate);
}
void q_mixer_tr::processing()
{

  if(freq_recon==true)
    freq=freq_in->read();

  stepsize = sample_time*freq*2.*M_PI;

  for(int i=0; i<rate; i++)
    {
      out.write(i_in.read()*(1+e/2)*amp*sin( actval + theta/2)+q_in.read()*amp*(1-e/2)*cos( actval -theta/2),i);
      actval+=stepsize;
    }
}

q_mixer_tr::q_mixer_tr(sc_core::sc_module_name n, double _freq, double _amp, int data_rate, bool f_config,double _theta, double _e)
{
  freq_recon = false;
  if(f_config)
    {
      freq_in = new sca_tdf::sca_in<double>;
      freq_recon = true;
    }
  amp=_amp;
  actval = 0.;
  freq=_freq;
  rate=data_rate;
  theta=_theta;
  e=_e;
}

void q_mixer_tr::freq_con (sca_tdf::sca_signal<double>* f_in)
{
  (*freq_in)(*f_in);
  freq_recon = true;
}
/****** qam mapper ********/

    void qam_map_phi::set_attributes()
    {
      //  in.set_rate(rate_in);
    }

    void qam_map_phi::initialize()  {bit_num=0;wait=0;}

    void qam_map_phi::processing() {
      if(reset.read())                                      //?low active?
	{
	  cout<< bit_num<< endl;
	  symbol[bit_num] = in.read();
	  bit_num++;
	  if(bit_num==6)
	    {
	      switch( rate_in ) {

	      case 2  : i_o = (qam4(symbol.range((rate_in/2-1),0)))  ;
		q_o = (qam4(symbol.range(rate_in-1,rate_in/2))) ;
		break ;
	      case 4  : i_o = (qam16(symbol.range((rate_in/2-1),0))) ;
		q_o = (qam16(symbol.range(rate_in-1,rate_in/2))) ;
		break ;
	      case 6  : i_o = (qam64(symbol.range((rate_in/2-1),0))) ;
		q_o = (qam64(symbol.range(rate_in-1,rate_in/2))) ;
		break ;
	      case 8  : i_o = (qam256(symbol.range((rate_in/2-1),0)));
		q_o = (qam256(symbol.range(rate_in-1,rate_in/2)));
		break ;
	      default : cout << " data rate not supported  " << endl;
		sc_core::sc_stop();	// stop simulation
		break ;
	      }
	      bit_num=0;
	      valid.write(true);
	      out_i.write(i_o) ;
	      out_q.write(q_o) ;
	    }
	  else
	    {
	      valid.write(false);
	      out_i.write(0) ;
	      out_q.write(0) ;
	    }
	}
      else
	bit_num=0;
    }

    double qam_map_phi::qam4(sc_bv<1> data) {

      double val = 0.0 ;

      sc_uint<1> tmp = data;
      switch (tmp) {
      case 0: val = -1;break;
      case 1: val = 1;break;
      default : cout << "Strong noise or attenuation!" << data << "\n";
	break;
      }

      return val;
    }

    double qam_map_phi::qam16(sc_bv<2> data) {

      double val = 0.0 ;
      sc_uint<2> tmp = data;
      switch (tmp) {
      case 0: val = -3;
	break;
      case 1: val = -1;
	break;
      case 2: val = 3;
	break;
      case 3: val = 1;
	break;
      default : cout << "Strong noise or attenuation!" << data << "\n";
	break;
      }
      return val;
    }

    double qam_map_phi::qam64(sc_bv<3> data) {

      double val = 0.0 ;
      sc_uint<3> tmp = data;
      switch (tmp) {
      case 0: val = -7;
	break;
      case 1: val = -5;
	break;
      case 2: val = -1;
	break;
      case 3: val = -3;
	break;
      case 4: val = 7;
	break;
      case 5: val = 5;
	break;
      case 6: val = 1;
	break;
      case 7: val = 3;
	break;
      default : cout << "Strong noise or attenuation!" << data << "\n";
	break;
      }

      return val;
    }

    double qam_map_phi::qam256(sc_bv<4> data)  {

      double val = 0.0 ;

      sc_uint<4> tmp = data;
      switch (tmp)
	{
	case 0: val = -15; break;	// werte noch anpassen
	case 1: val = -13; break;
	case 2: val = -11; break;
	case 3: val = -9; break;
	case 4: val = -7;  break;
	case 5: val = -5;  break;
	case 6: val = -3;  break;
	case 7: val = -1;  break;
	case 8: val = 15; break;
	case 9: val = 13; break;
	case 10: val = 11; break;
	case 11: val = 9; break;
	case 12: val = 7;  break;
	case 13: val = 5;  break;
	case 14: val = 3;  break;
	case 15: val = 1;  break;
	default : cout << "Strong noise or attenuation!" << data << "\n";
	  break;
	}

      return val;
    }


    qam_map_phi::qam_map_phi(sc_core::sc_module_name nm, int dim_cons){
      if((dim_cons==4)||(dim_cons==16)||(dim_cons==64)||(dim_cons==256))
	rate_in = static_cast<int>(log2(dim_cons));
      else
	{
	  cout<<"Error:" << dim_cons <<" QAM is not supported yet!"<<endl;
	  cout<<"Possible value is 4,16,64 or 256"<<endl;
	  exit(0);
	}
    }
    //********************
    void qam_map::set_attributes()
    {
      in.set_rate(rate_in);
    }

    void qam_map::initialize()  {}

    void qam_map::processing() {
      for (int i = 0; i < rate_in ; i++) {
	symbol[i] = in.read(i);
      }

      switch( rate_in ) {

      case 2  : i_o = (qam4(symbol.range((rate_in/2-1),0)))  ;
	q_o = (qam4(symbol.range(rate_in-1,rate_in/2))) ;
	break ;
      case 4  : i_o = (qam16(symbol.range((rate_in/2-1),0))) ;
	q_o = (qam16(symbol.range(rate_in-1,rate_in/2))) ;
	break ;
      case 6  : i_o = (qam64(symbol.range((rate_in/2-1),0))) ;
	q_o = (qam64(symbol.range(rate_in-1,rate_in/2))) ;
	break ;
      case 8  : i_o = (qam256(symbol.range((rate_in/2-1),0)));
	q_o = (qam256(symbol.range(rate_in-1,rate_in/2)));
	break ;
      default : cout << " data rate not supported  " << endl;
	sc_core::sc_stop();	// stop simulation
	break ;
      }

      out_i.write(i_o) ;
      out_q.write(q_o) ;
    }

    double qam_map::qam4(sc_bv<1> data) {

      double val = 0.0 ;

      sc_uint<1> tmp = data;
      switch (tmp) {
      case 0: val = -1;
	break;
      case 1: val = 1;
	break;
      default : cout << "Strong noise or attenuation!" << data << "\n";
	break;
      }

      return val;
    }

    double qam_map::qam16(sc_bv<2> data) {

      double val = 0.0 ;
      sc_uint<2> tmp = data;
      switch (tmp) {
      case 0: val = -3;
	break;
      case 1: val = -1;
	break;
      case 2: val = 3;
	break;
      case 3: val = 1;
	break;
      default : cout << "Strong noise or attenuation!" << data << "\n";
	break;
      }

      return val;
    }

    double qam_map::qam64(sc_bv<3> data) {

      double val = 0.0 ;
      sc_uint<3> tmp = data;
      switch (tmp) {
      case 0: val = -7;
	break;
      case 1: val = -5;
	break;
      case 2: val = -3;
	break;
      case 3: val = -1;
	break;
      case 4: val = 7;
	break;
      case 5: val = 5;
	break;
      case 6: val = 3;
	break;
      case 7: val = 1;
	break;
      default : cout << "Strong noise or attenuation!" << data << "\n";
	break;
      }

      return val;
    }

    double qam_map::qam256(sc_bv<4> data)  {

	double val = 0.0 ;

	sc_uint<4> tmp = data;
	switch (tmp)
	  {
	  case 0: val = -15; break;	// werte noch anpassen
	  case 1: val = -13; break;
	  case 2: val = -11; break;
	  case 3: val = -9; break;
	  case 4: val = -7;  break;
	  case 5: val = -5;  break;
	  case 6: val = -3;  break;
	  case 7: val = -1;  break;
	  case 8: val = 15; break;
	  case 9: val = 13; break;
	  case 10: val = 11; break;
	  case 11: val = 9; break;
	  case 12: val = 7;  break;
	  case 13: val = 5;  break;
	  case 14: val = 3;  break;
	  case 15: val = 1;  break;
	  default : cout << "Strong noise or attenuation!" << data << "\n";
	    break;
	  }

	return val;
    }


    qam_map::qam_map(sc_core::sc_module_name nm, int _rate){
      if((_rate==4)||(_rate==16)||(_rate==64)||(_rate==256))
	rate_in = static_cast<int>(log2((BB_DBL)_rate));
      else
	{
	  cout<<"Error:" << _rate <<" QAM is not supported yet!"<<endl;
	  cout<<"Possible value is 4,16,64 or 256"<<endl;
	  exit(0);
	}

    }


//********************************************************************

/****** qam demapper ********/

void qam_demap::set_attributes()
{
  out.set_rate(rate_out);
}

void qam_demap::initialize()  {}

void qam_demap::processing() {

  vector<uint> i_tmp;
  vector<uint> q_tmp;

  switch( rate_out ) {

  case 2  : for(int i=0;i<1;i++)
	      {
		i_tmp.push_back(qam4(error_correction(in_i.read()))[i]);
	      }
	    for(int i=0;i<1;i++)
	      {
		q_tmp.push_back(qam4(error_correction(in_q.read()))[i]);
	      }
	    break ;
  case 4  : for(int i=0;i<2;i++)
	      {
		i_tmp.push_back(qam16(error_correction(in_i.read()))[i]);
	      }
	    for(int i=0;i<2;i++)
	      {
		q_tmp.push_back(qam16(error_correction(in_q.read()))[i]);
	      }
	    break ;
  case 6  : for(int i=0;i<3;i++)
	      {
		i_tmp.push_back(qam64(error_correction(in_i.read()))[i]);
	      }
	    for(int i=0;i<3;i++)
	      {
		q_tmp.push_back(qam64(error_correction(in_q.read()))[i]);
	      }
	    break ;
  case 8  : for(int i=0;i<4;i++)
	      {
		i_tmp.push_back(qam256(error_correction(in_i.read()))[i]);
	      }
	    for(int i=0;i<4;i++)
	      {
		q_tmp.push_back(qam256(error_correction(in_q.read()))[i]);
	      }
	    break ;
  default : cout << " data rate not supported  " << endl;
            sc_core::sc_stop();	// stop simulation
	    break ;
   }


  for (int n = 0; n < rate_out; n++)
    {
      if(n<rate_out/2)
	out.write(i_tmp[n] ,n);
      else
	out.write(q_tmp[n-rate_out/2],n);
    }
  i_tmp.clear();
  q_tmp.clear();

}

int qam_demap::error_correction(double x){

  int val=0;

  if(x>1-delta && x<1+delta)
    val=1;
  else if(x > -1-delta && x < -1+delta)
    val=-1;
  else if(x > 3-delta && x < 3+delta)
    val=3;
  else if(x > -3-delta && x < -3+delta)
    val=-3;
  else if(x > 5-delta && x < 5+delta)
    val=5;
  else if(x > -5-delta && x < -5+delta)
    val=-5;
  else if(x > 7-delta && x < 7+delta)
    val=7;
  else if(x > -7-delta && x < -7+delta)
    val=-7;
  else if(x > 9-delta && x < 9+delta)
    val=9;
  else if(x > -9-delta && x < -9+delta)
    val=-9;
  else if(x > 11-delta && x < 11+delta)
    val=11;
  else if(x > -11-delta && x < -11+delta)
    val=-11;
  else if(x > 13-delta && x < 13+delta)
    val=13;
  else if(x > -13-delta && x < -13+delta)
    val=-13;
  else if(x > 15-delta && x < 15+delta)
    val=15;
  else if(x > -15-delta && x <-15+delta)
    val=-15;
  else
    cout << "error input" << x << "\n";

  //  cout<< "hallo" << val <<endl;
  return val;
}

sc_uint<1> qam_demap::qam4(int data) {

	sc_uint<1> val = 0 ;

 	switch (data) {
		  case -1: val = 0;
		  	  break;
		  case 1:  val = 1;
		          break;
		  default : cout << "  " << data << "\n";
		          break;
  		}

	return val;
}

sc_uint<2> qam_demap::qam16(int data) {

	sc_uint<2> val = 0;

	switch (data) {
		  case -3: val = 00;
		          break;
		  case -1: val = 01;
		          break;
  		  case 3:  val = 10;
		          break;
  		  case 1:  val = 11;
		          break;
		  default : cout << "Inaviable input " << data << "\n";
			    break;
  		}

	return val;

}

sc_uint<3> qam_demap::qam64(int data) {

	sc_uint<3> val = 0 ;
	switch (data) {
		  case -7: val = 0;
		          break;
		  case -5: val = 1;
		          break;
  		  case -3: val = 2;
		          break;
  		  case -1: val = 3;
		          break;
		  case 7:  val = 4;
		          break;
		  case 5:  val = 5;
		          break;
		  case 3:  val = 6;
		          break;
		  case 1:  val = 7;
		          break;
		  default : cout << "  " << data << "\n";
			  break;
  		}
	return val;
}


sc_uint<4> qam_demap::qam256(int data)  {

        sc_uint<4> val = 0 ;
	switch (data)
	  {
	  case -15: val = 0; break;
	  case -13: val = 1; break;
	  case -11: val = 2; break;
	  case -9 : val = 3; break;
	  case -7:  val = 4; break;
	  case -5:  val = 5; break;
	  case -3:  val = 6; break;
	  case -1:  val = 7; break;
	  case 15:  val = 8; break;
	  case 13:  val = 9; break;
	  case 11:  val = 10; break;
	  case 9 :  val = 11; break;
	  case 7 :  val = 12; break;
	  case 5 :  val = 13; break;
	  case 3 :  val = 14; break;
	  case 1 :  val = 15; break;
	  default : cout << "  " << data << "\n";
	    break;
	  }
	return val;


}

qam_demap::qam_demap(sc_core::sc_module_name nm, int _rate, double _delta){
    if((_rate==4)||(_rate==16)||(_rate==64)||(_rate==256))
    {
	rate_out =static_cast<int>(log2((BB_DBL)_rate));
	delta=_delta;
    }
    else
    {
	cout<<"Error: " << _rate << "DeQAM is not supported yet!"<<endl;
	cout<<"Possible value is 4, 16, 64 or 256!";
	exit(0);
    }
}




/****************************** QPSK modulator*******************************/

qpsk::qpsk(sc_core::sc_module_name n, double _freq, int rate)
    {
      s2p_sub = new s2p<bool,2>("i_s2p",1);
      s2p_sub->in(in);
      s2p_sub->out[0](sig_i);
      s2p_sub->out[1](sig_q);

      nrz_i_sub = new nrz("i_sub",1.0);
      nrz_i_sub->in(sig_i);
      nrz_i_sub->out(sig_n_i);

      nrz_q_sub = new nrz("q_sub",1.0);
      nrz_q_sub->in(sig_q);
      nrz_q_sub->out(sig_n_q);

      mixer_sub = new q_mixer_tr("i_mix",_freq,1.0,rate,false);
      mixer_sub -> i_in(sig_n_i);
      mixer_sub -> q_in(sig_n_q);
      mixer_sub -> out(out);
    }

/****************************** QPSK demodulator*******************************/

qpsk_de::qpsk_de(sc_core::sc_module_name n, double _freq, int rate)
{
  p2s_sub = new p2s<bool,2>("i_s2p",1);
  p2s_sub->out(out);
  p2s_sub->in[0](sig_ci);
  p2s_sub->in[1](sig_cq);

  lp_i_sub = new lp("i_sub",_freq/10);
  lp_i_sub->in(sig_i);
  lp_i_sub->out(sig_lp_i);

  lp_q_sub = new lp("q_sub",_freq/10);
  lp_q_sub->in(sig_q);
  lp_q_sub->out(sig_lp_q);

  decider_i = new downsample("i_decider",rate-1,rate);
  decider_i->in(sig_lp_i);
  decider_i->out(sig_di);

  decider_q = new downsample("q_decider",rate-1,rate);
  decider_q->in(sig_lp_q);
  decider_q->out(sig_dq);

  comp_i_sub = new compare("i_comp");
  comp_i_sub->in(sig_di);
  comp_i_sub->out(sig_ci);

  comp_q_sub = new compare("q_comp");
  comp_q_sub->in(sig_dq);
  comp_q_sub->out(sig_cq);

  mixer_sub = new q_mixer_re("i_mix",_freq,1.0,rate,false);
  mixer_sub -> i_out(sig_i);
  mixer_sub -> q_out(sig_q);
  mixer_sub -> in(in);
}

/****************************** OQPSK modulator*******************************/

oqpsk:: oqpsk(sc_core::sc_module_name n, double _freq, int rate)
    {
      s2p_sub = new s2p<bool,2>("i_s2p",1);
      s2p_sub->in(in);
      s2p_sub->out[0](sig_i);
      s2p_sub->out[1](sig_q);

      delay_sub = new d_flipflop("i_d_fp",false);
      delay_sub->in(sig_q2);
      delay_sub->out(sig_dq);

      nrz_i_sub = new nrz("i_sub",1.0);
      nrz_i_sub->in(sig_i2);
      nrz_i_sub->out(sig_n_i);

      o2t_sub = new o2t("i_o2t");
      o2t_sub -> in(sig_q);
      o2t_sub -> out(sig_q2);

      o2t_i_sub = new o2t("i_io2t");
      o2t_i_sub -> in(sig_i);
      o2t_i_sub -> out(sig_i2);

      nrz_q_sub = new nrz("q_sub",1.0);
      nrz_q_sub->in(sig_dq);
      nrz_q_sub->out(sig_n_q);

      mixer_sub = new q_mixer_tr("i_mix",_freq,1.0,rate,false);
      mixer_sub -> i_in(sig_n_i);
      mixer_sub -> q_in(sig_n_q);
      mixer_sub -> out(out);
    }

/****************************** OQPSK demodulator*******************************/

oqpsk_de::oqpsk_de(sc_core::sc_module_name n, double _freq, int rate)
{
  p2s_sub = new p2s<bool,2>("i_s2p",1);
  p2s_sub->out(out);
  p2s_sub->in[0](sig_ito);
  p2s_sub->in[1](sig_qto);

  lp_i_sub = new lp("i_sub",_freq/10);
  lp_i_sub->in(sig_i);
  lp_i_sub->out(sig_lp_i);

  lp_q_sub = new lp("q_sub",_freq/10);
  lp_q_sub->in(sig_q);
  lp_q_sub->out(sig_lp_q);

  t2o_q_sub = new t2o("i_qo2t");
  t2o_q_sub -> in(sig_deq);
  t2o_q_sub -> out(sig_qto);

  t2o_i_sub = new t2o("i_io2t");
  t2o_i_sub -> in(sig_dei);
  t2o_i_sub -> out(sig_ito);

  delay_qsub = new d_flipflop("i_qd_fp",false);
  delay_qsub->in(sig_cq);
  delay_qsub->out(sig_deq);

  delay_isub = new d_flipflop("i_id_fp",false,2);
  delay_isub ->in(sig_ci);
  delay_isub ->out(sig_dei);

  decider_i = new downsample("i_decider",rate-1,rate);
  decider_i->in(sig_lp_i);
  decider_i->out(sig_di);

  decider_q = new downsample("q_decider",rate-1,rate);
  decider_q->in(sig_lp_q);
  decider_q->out(sig_dq);

  comp_i_sub = new compare("i_comp");
  comp_i_sub->in(sig_di);
  comp_i_sub->out(sig_ci);

  comp_q_sub = new compare("q_comp");
  comp_q_sub->in(sig_dq);
  comp_q_sub->out(sig_cq);

  mixer_sub = new q_mixer_re("i_mix",_freq,1.0,rate,false);
  mixer_sub -> i_out(sig_i);
  mixer_sub -> q_out(sig_q);
  mixer_sub -> in(in);
}


/****************************** DBPSK modulator *******************************/
    dbpsk::dbpsk(sc_core::sc_module_name n, double _freq, int rate)
    {
      d_code_sub  = new d_coding("i_dcode");
      d_code_sub->in(in);
      d_code_sub->out(sig_dcode);

      bpsk_sub = new bpsk("i_bpsk",_freq,rate);
      bpsk_sub->in(sig_dcode);
      bpsk_sub->out(out);
    }
/****************************** DBPSK demodulator *******************************/
    dbpsk_de::dbpsk_de(sc_core::sc_module_name n, double _freq, int rate)
    {
      d_decode_sub  = new d_decoding("i_decode");
      d_decode_sub->in(sig_dem);
      d_decode_sub->out(out);

      bpsk_de_sub = new bpsk_de("i_bpsk_de", _freq, rate);
      bpsk_de_sub->in(in);
      bpsk_de_sub->out(sig_dem);
    }




/****************************** DQPSK modulator*******************************/
dqpsk::dqpsk(sc_core::sc_module_name n, double _freq, int rate)
{
  s2p_sub = new s2p<bool,2>("i_s2p",1);
  s2p_sub->in(in);
  s2p_sub->out[0](sig_i);
  s2p_sub->out[1](sig_q);

  encode_i_sub = new d_coding("iencode");
  encode_i_sub -> in(sig_i);
  encode_i_sub -> out(sig_i_encode);

  encode_i_sub = new d_coding("qencode");
  encode_i_sub -> in(sig_q);
  encode_i_sub -> out(sig_q_encode);

  nrz_i_sub = new nrz("i_sub",1.0);
  nrz_i_sub->in(sig_i_encode);
  nrz_i_sub->out(sig_nrz_i);

  nrz_q_sub = new nrz("q_sub",1.0);
  nrz_q_sub->in(sig_q_encode);
  nrz_q_sub->out(sig_nrz_q);

  mixer_sub = new q_mixer_tr("i_mix",_freq,1.0,rate,false);
  mixer_sub -> i_in(sig_nrz_i);
  mixer_sub -> q_in(sig_nrz_q);
  mixer_sub -> out(out);

}
/****************************** DQPSK demodulator*******************************/
dqpsk_de::dqpsk_de(sc_core::sc_module_name n, double _freq, int rate)
{

  mixer_sub = new q_mixer_re("i_mix",_freq,1.0,rate,false);
  mixer_sub -> i_out(sig_i);
  mixer_sub -> q_out(sig_q);
  mixer_sub -> in(in);

  lp_i_sub = new lp("i_sub",_freq/10);
  lp_i_sub->in(sig_i);
  lp_i_sub->out(sig_lp_i);

  lp_q_sub = new lp("q_sub",_freq/10);
  lp_q_sub->in(sig_q);
  lp_q_sub->out(sig_lp_q);

  decider_i = new downsample("i_decider",rate-1,rate);
  decider_i->in(sig_lp_i);
  decider_i->out(sig_di);

  decider_q = new downsample("q_decider",rate-1,rate);
  decider_q->in(sig_lp_q);
  decider_q->out(sig_dq);

  comp_i_sub = new compare("i_comp");
  comp_i_sub->in(sig_di);
  comp_i_sub->out(sig_ci);

  comp_q_sub = new compare("q_comp");
  comp_q_sub->in(sig_dq);
  comp_q_sub->out(sig_cq);

  decode_isub  = new d_decoding("i_decode");
  decode_isub->in(sig_ci);
  decode_isub->out(sig_deci);

  decode_qsub  = new d_decoding("q_decode");
  decode_qsub->in(sig_cq);
  decode_qsub->out(sig_decq);

  p2s_sub = new p2s<bool,2>("i_p2s",1);
  p2s_sub->out(out);
  p2s_sub->in[0](sig_deci);
  p2s_sub->in[1](sig_decq);}

/*********Network analyzer***********/
/*sub module 1: measurement*/
void measure_sub::initialize()
{
  frequenz = start_f;
  time_buffer = 0;

  dut_buffer(0) = 0;
  dut_buffer(1) = 0;
  ref_buffer(0) = 0;
  ref_buffer(1) = 0;

  find_beginning = 0;
  x = 0;
  start_ref = 1;
  start_dut = 0;
  start = 1;
  start_measurement = 0;

  sim_time();			//Simulation time calculation
}

void measure_sub::processing()
{

  //Calculation of time differenc
  time = get_time().to_seconds() - time_buffer;

  //Measurement
  measurement();
}

void measure_sub::measurement()
{

  //Intercept when max. frequency is reached
  if(frequenz > stop_f){

  }//Measure current signal
  else if((time < stop_measurement_time)&&(start == 1)){

    if(start_measurement_time < time){

      //Seaching for the last periods
      if((find_beginning == 0)&&(ref.read() < 0)){

	find_beginning = 1;

      }else if(((find_beginning == 1) && (ref.read() > 0))||(start_measurement == 1)){

	//With start_measurement the Measurement
	start_measurement = 1;

	//Search for max. reference values in current period
	if((ref_buffer(0) < ref.read()) && start_ref == 1){

	  //Save reference value
	  ref_buffer(0) = ref.read();
	  //Save time of measurement
	  ref_buffer(1) = get_time().to_seconds();

	}else{

	  start_ref = 0;			//Stop with the measurement of the ref signal
	  start_dut = 1;			//Start with the measurement of the dut signal


	}

	//search for max. measurement values in current period
	if((dut_buffer(0) < in.read()) && (start_dut == 1)){

	  //Save measuremant value
	  dut_buffer(0) = in.read();
	  //Save time of measurement
	  dut_buffer(1) = get_time().to_seconds();
	  x = 1;
	}else  if((in.read() < dut_buffer(0)) && (x == 1)){

	  x = 0;
	  start = 0;			//Stop with the measurement
	  start_dut = 0;			//Stop with the measurement of the dut signal
	}
      }
    }


  }else{


    //When period has ended, the measurement values are processed and input to the svg file

    //Calculation of the amplitude response function in dB
    dB = 20 * log10(dut_buffer(0)/ref_buffer(0));

    //Calculation of the phase displacement in degrees
    phase = 360*frequenz*(ref_buffer(1) - dut_buffer(1));

    //Output
    cout << "Frequenz "<< "\t" << frequenz <<"\t"<< "H [dB]" << "\t\t" << dB <<endl;
    cout << "Frequenz "<< "\t" << frequenz <<"\t"<< "phase [Grad]" << "\t" << phase <<endl;
    cout << endl;

    //Storing values in the txt-file
    measured_value << "Frequenz [Hz] \t" << frequenz << "\t\t" << "H [dB] \t" << dB << "\t\t" << " phase [Grad] \t" << phase << endl;
    measured_value << endl;

    //svg output
    //Range limitation
    if(min_dB <= dB && dB <= max_dB){
      ss1 << 1280*(log10(frequenz)/scale_x) - scale_zero_x  << " " << -dB/(max_dB-min_dB)*400 + scale_zero_y << ", ";

      ss2 << 1280*(log10(frequenz)/scale_x) - scale_zero_x  << " " << 400-phase/360*200 << ", ";
    }

    //Buffer reinitialisation
    dut_buffer(0) = 0;
    dut_buffer(1) = 0;
    ref_buffer(0) = 0;
    ref_buffer(1) = 0;

    find_beginning = 0;
    x = 0;
    start_ref = 1;
    start_dut = 0;
    start = 1;
    start_measurement = 0;



    //Frequency increase
    frequenz = frequenz + step_f;
    time_buffer = get_time().to_seconds();

  }

}

void measure_sub::make_file(){

  //TXT-FILE----------------------------------------------------------------------------------------------
  string name_m=name();				//Use modul name as file name
  name_m+=".txt";					//Attach txt as suffix
  measured_value.open(name_m.c_str(), ios::out);	//Open file

  //SVG-FILE----------------------------------------------------------------------------------------------
  //Scaling of x-coordinate
  scale_x = log10(stop_f)-log10(start_f);
  //Determination of zero point for x-coordinate
  scale_zero_x = 1280*(log10(start_f)/scale_x);
  //Determination of zero point for x-coordinate for amplitude response
  scale_zero_y = max_dB/(max_dB-min_dB)*400;

  string filename=name();				//Use modul name as file name
  filename+=".svg";				//Attach svg as suffix
  output.open( filename.c_str(), ios::out);	//Open file

  //Add XML/SVG-Header
  output << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\" ?>" << endl;
  output << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20010904//EN\"" << endl;
  output << "  \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl;

  //SVG-Image Information
  output << "<svg width=\"1280\" height=\"600\" xmlns=\"http://www.w3.org/2000/svg\"" << endl;
  output << "  xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
  output << "  <title>Trace aus SystemC-AMS</title>" << endl;


  //Amplitude response graph ------------------------------------------------------------------------------
  //Frame for amplitude response graph
  output << "<rect x=\"0\" y=\"0\" width=\"1280\" height=\"400\" stroke=\"blue\" stroke-width=\"4px\" fill=\"none\"/>" << endl;
  //Creation of x-axis for amplitude response graph
  output << "<line stroke=\"red\" stroke-width=\"1px\" x1=\"0\" y1=\" "<< scale_zero_y << "\" x2=\"1280\" y2=\" " << scale_zero_y << "\"/>" << endl;

  // Y distance lines
  if((stop_f-start_f)/step_f <= 20 ){

    //Draw line at each measurement point(up to a maximum of 20 measurement point)

    int x = 0;		//Dummy variable for legend

    for(double i = start_f; i <= stop_f; i +=step_f){

      //Calculation of x-coordinate for distance lines
      int buffer = (int)((log10(i)/scale_x)*1280 - scale_zero_x);

      //Alternating legends (even measurement points below the x-axis;
      //Odd measurement points above the x-axis)
      if(x % 2){
	//Distance line
	output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\""<< buffer <<"\" y1=\"0\" x2=\""<< buffer <<"\" y2=\"600\"/>" << endl;
	//Legend
	output << "<text x=\""<< buffer+2 <<"\" y=\""<< scale_zero_y - 2 <<"\" fill=\"grey\" font-size=\"8px\">" << i <<"Hz</text>" << endl;
      }else{
	//Distance line
	output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\""<< buffer <<"\" y1=\"0\" x2=\""<< buffer <<"\" y2=\"600\"/>" << endl;
	//Legend
	output << "<text x=\""<< buffer+2 <<"\" y=\""<< scale_zero_y + 8 <<"\" fill=\"grey\" font-size=\"8px\">" << i <<"Hz</text>" << endl;
      }

      x += 1;
    }
  }else{

    //Draw line at each decade
    for(double i = start_f; i <= stop_f; i +=step_f){
      //Calculation of x-coordinate for distance lines
      int buffer = (int)((log10(i)/scale_x)*1280 - scale_zero_x);

      //Draw line at each decade
      if(log10(i) - (int)log10(i) == 0){
	//Distance line
	output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\""<< buffer <<"\" y1=\"0\" x2=\""<< buffer <<"\" y2=\"600\"/>" << endl;
	//Legend
	output << "<text x=\""<< buffer+2 <<"\" y=\""<< scale_zero_y + 8 <<"\" fill=\"grey\" font-size=\"8px\">" << i <<"Hz</text>" << endl;
      }
    }
  }


  // X distance lines
  if(max_dB-min_dB <= 30){

    //Draw line at each integer dB value (up to a maximum of 30 values)

    for(int i = (int)max_dB; i >= (int)min_dB; i -=1){

      //Calculation of y-coodinate
      int buffer = (int)(scale_zero_y - i/(max_dB - min_dB)*400);
      //Distance line
      output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\""<< buffer <<"\" x2=\"1280\" y2=\""<< buffer <<"\"/>" << endl;
      //Legend
      output << "<text x=\"10\" y=\""<< buffer-3 <<"\" fill=\"grey\" font-size=\"8px\">" << i <<"dB</text>" << endl;

    }
  }else{
    //Draw line in 5dB steps, if the range is more than 30 values
    for(int i = (int)max_dB; i >= (int)min_dB; i -=1){

      //Calculation of y-coodinate
      int buffer = (int)(scale_zero_y - i/(max_dB - min_dB)*400);

      //Draw only mod 5
      if(i % 5 == 0){
	//Distance line
	output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\""<< buffer <<"\" x2=\"1280\" y2=\""<< buffer <<"\"/>" << endl;
	//Legend
	output << "<text x=\"10\" y=\""<< buffer-3 <<"\" fill=\"grey\" font-size=\"8px\">" << i <<"dB</text>" << endl;
      }
    }
  }


  //Phase graph -------------------------------------------------------------------------------------
  //Frame for phase graph
  output << "<rect x=\"0\" y=\"400\" width=\"1280\" height=\"200\" stroke=\"blue\" stroke-width=\"4px\" fill=\"none\"/>" << endl;

  //Draw distance lines in 45-degree steps
  for(double i = 0; i < 360; i += 45){

    //Calculation of y-coodinate
    int buffer = (int)(i/360*200 + 400);
    //Distance line
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"" << buffer << "\" x2=\"1280\" y2=\"" << buffer << "\"/>" << endl;
    //Legend
    output << "<text x=\"10\" y=\""<< buffer+8 <<"\" fill=\"grey\" font-size=\"8px\"> -" << i << " Grad</text>" << endl;
  }
}

double measure_sub::sim_time(){

  double sum = 0;		//Dummy variable

  //Calculation of each period lag
  for(double i = start_f; i <= stop_f; i += step_f){

    //Addition of all period lags
    sum = sum + period_number*(1/start_f);

  }

  cout << "--> required SystemC-Time " << sum <<endl;
  measured_value << "--> required SystemC-Time " << sum << endl;
  measured_value << endl;

  return sum;
}

void measure_sub::finish(){
  string 	s1=ss1.str();		//Transformation of amplitude response values into a string
  //Creation of the Polyline
  output << "<polyline fill=\"none\" stroke=\"black\" stroke-width=\"1px\"" << endl;
  //Assignment of Polyline coordinates
  output << "points=\"" << s1 << "\" />" << endl;

  string 	s2=ss2.str();		//Transformation of phase values into a string
  //Creation of the Polyline
  output << "<polyline fill=\"none\" stroke=\"black\" stroke-width=\"1px\"" << endl;
  //Assignment of Polyline coordinates
  output << "points=\"" << s2 << "\" />" << endl;

  //Close svg file
  output << "</svg>" << endl;
}

measure_sub::measure_sub(sc_core::sc_module_name n,double _start_f,double _stop_f,double _step_f,double _min_dB,double _max_dB,double _period_number)
{

  //Initalisation of Parameter
  start_f = _start_f;
  stop_f = _stop_f;
  step_f = _step_f;
  min_dB = _min_dB;
  max_dB = _max_dB;
  period_number = _period_number;

  //The measurement starts with the two last periods of the start frequency
  start_measurement_time = period_number*(1/start_f)-2*(1/start_f);
  stop_measurement_time = period_number*(1/start_f);

  //Creating the svg-File
  make_file();
}

/*sub module 2: generator*/

void generator::initialize()
{
  frequenz = start_f;
  time_buffer = 0;
}

void generator::processing()
{

  //Calculation of time difference
  time = get_time().to_seconds() - time_buffer;

  //Intercept when max. frequency is reached
  if(frequenz > stop_f){

    out.write(0.0);

  }//generating current signal
  else if(time <= (period_number*(1/start_f))){

    out.write(amplitude*sin(2*M_PI*frequenz*get_time().to_seconds()));
    ref.write(amplitude*sin(2*M_PI*frequenz*get_time().to_seconds()));

  }else{

    //Frequency increase
    frequenz = frequenz + step_f;
    time_buffer = get_time().to_seconds();

  }
}

generator::generator(sc_core::sc_module_name n, double _amplitude, double _start_f, double _stop_f, double _step_f, double _period_number) {

  //Initalisation of Parameter
  amplitude = _amplitude;
  start_f = _start_f;
  stop_f = _stop_f;
  step_f = _step_f;
  period_number = _period_number;

}

/*network analyzer*/

void nwa::finish()
{
  p_measure->finish();

}
nwa::nwa(sc_core::sc_module_name n, double _amplitude, double _start_f, double _stop_f, double _step_f, double _min_dB, double
			     _max_dB, double _period_number){

  //Module
  p_measure = new measure_sub("measure", _start_f, _stop_f, _step_f, _min_dB, _max_dB, _period_number);
  p_generator = new generator("generator", _amplitude, _start_f, _stop_f, _step_f, _period_number);

  //--------------------------------------------------------------------------------------
  //Connection----------------------------------------------------------------------------

  //Reference signal
  p_generator->ref(ref_signal);
  p_measure->ref(ref_signal);

  //Out-Port
  p_generator->out(out);

  //In-Port
  p_measure->in(in);

}


/****************************************Eye-Diagram****************************************************/
eyediag::eyediag(sc_core::sc_module_name n, double periodetime_, double sigamp_, int periodes_, int delay_, int in_rate)
{
// Uebernehme Instanzvariablen
data_rate=in_rate;
periodetime=2*periodetime_;
sigamp=sigamp_;
wavemax=periodes_;
delay=delay_;
lastx=0;

string filename=name();
filename+=".svg";
output.open( filename.c_str(), ios::out); // öffnet Filehandle

// XML/SVG-Header
output << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\" ?>" << endl;
output << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20010904//EN\"" << endl;
output << "  \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl;

// SVG-Image Informationen
output << "<svg width=\"800\" height=\"600\" xmlns=\"http://www.w3.org/2000/svg\"" << endl;
output << "  xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
output << "  <title>Eye-Diagram by SystemC-AMS</title>" << endl;

// Zeichenebene Vorbereiten
// Blauer Rahmen um den Plot
output << "<rect x=\"0\" y=\"0\" width=\"800\" height=\"600\" stroke=\"blue\" stroke-width=\"4px\" fill=\"white\"/>" << endl;

// Rote Linien als x/y-Achsen
output << "<line stroke=\"red\" stroke-width=\"1px\" x1=\"0\" y1=\"300\" x2=\"800\" y2=\"300\"/>" << endl;
output << "<line stroke=\"red\" stroke-width=\"1px\" x1=\"400\" y1=\"0\" x2=\"400\" y2=\"600\"/>" << endl;

// Pfeile
output << "<polygon fill=\"red\" points=\"396 8 404 8 400 0\" />" << endl;
output << "<polygon fill=\"red\" points=\"792 296 800 300 792 304\" />" << endl;

// Raster
output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"200\" x2=\"800\" y2=\"200\"/>" << endl;
output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"100\" x2=\"800\" y2=\"100\"/>" << endl;
output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"400\" x2=\"800\" y2=\"400\"/>" << endl;
output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"500\" x2=\"800\" y2=\"500\"/>" << endl;

// Rasterlinien-Beschriftung
output << "<text x=\"5\" y=\"198\" fill=\"grey\" font-size=\"10px\">" << (sigamp/10*2) << "V</text>" << endl;
output << "<text x=\"5\" y=\"98\"  fill=\"grey\" font-size=\"10px\">" << (sigamp/10*4) << "V</text>" << endl;
output << "<text x=\"5\" y=\"298\" fill=\"grey\" font-size=\"10px\"> 0 V</text>" << endl;
output << "<text x=\"5\" y=\"398\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/10*2) << "V</text>" << endl;
output << "<text x=\"5\" y=\"498\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/10*4) << "V</text>" << endl;

//Zeitachsen-Beschriftung
output << "<text x=\"5\" y=\"312\" fill=\"red\" font-size=\"10px\">-" << (periodetime*1000/2) << " ms</text>" << endl;
output << "<text x=\"750\" y=\"312\" fill=\"red\" font-size=\"10px\">" << (periodetime*1000/2) << " ms</text>" << endl;
output << "<text x=\"402\" y=\"312\" fill=\"red\" font-size=\"10px\">0.0 ms</text>" << endl;

// Infotext ausgeben
output << "<text x=\"5\" y=\"595\" font-size=\"15px\">T=" << (periodetime*1000) << " ms | Vpp=" << sigamp << " V</text>" << endl;
output << "<text x=\"5\" y=\"17\" font-size=\"15px\">Eye-Diagram: \""<< name() <<"\" (" << wavemax << " Periodes | " << delay << " Periodes Delay)</text>" << endl;

// Polyline String vorbereiten
ss1 << "<polyline fill=\"none\" stroke=\"black\" stroke-width=\"1px\"" << endl;
ss1 << "points=\"" << endl;
}

void eyediag::set_attributes()
  {
    in.set_rate(data_rate);	// Liest 10 Token
  }

void eyediag::processing()
  {
    // Zentrierung des Auges
    double delta_t=periodetime/4;

    double time=get_time().to_seconds();

    time+=delta_t;

    int step=(int)(800*time/periodetime);
    int x=step%800;

    // Ausgabe nur der ersten X Perioden
    if( (time>(periodetime*delay)) && (time<(periodetime*(wavemax+delay))) )
    {
        if(x<lastx) // Beginnt neue Polylinie wenn das Ende der Zeichenebene erreicht wurde
	{
	  ss1 << "\" />" << endl;
   	  ss1 << "<polyline fill=\"none\" stroke=\"black\" stroke-width=\"1px\"" << endl;
	  ss1 << "points=\"" << endl;
	}

        lastx=x;

        // Normierung des Signals auf 1 - anschließende Pixelberechnung
        ss1 << x << " " << int(300 - in.read(0)/sigamp * 500) << ", ";
	// Durch read(0) wird nur jeder 10. Wert eingelesen

    }
  }

void eyediag::finish()
  {
    // Polyline abschließen
    ss1 << "\" />" << endl;

    // Gesamtes "Polyline-Informationen" in String umwandeln und in File schreiben
    string s=ss1.str();
    output << s << endl;

    // Abschluss des SVG-Dokuments und schließen des Filehandles
    output << "</svg>";
    output.close();
  }

/******************************************************Scatter-Plot************************************************/
scatter::scatter(sc_core::sc_module_name n, double sigamp_, int in_rate)
  {
    sigamp=sigamp_;
    data_rate=in_rate;
    string filename=name();
    filename+=".svg";

    output.open( filename.c_str(), ios::out);	// öffnet Filehandle

    // XML/SVG-Header
    output << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\" ?>" << endl;
    output << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20010904//EN\"" << endl;
    output << "  \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl;

    // SVG-Image Informationen
    output << "<svg width=\"600\" height=\"600\" xmlns=\"http://www.w3.org/2000/svg\"" << endl;
    output << "  xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
    output << "  <title>Scatter-Plot by SystemC-AMS</title>" << endl;

    // Zeichenebene Vorbereiten
    // Blauer Rahmen um den Plot
    output << "<rect x=\"0\" y=\"0\" width=\"600\" height=\"600\" stroke=\"blue\" stroke-width=\"4px\" fill=\"white\"/>" << endl;

    // Rote Linien als x/y-Achsen
    output << "<line stroke=\"red\" stroke-width=\"1px\" x1=\"0\" y1=\"300\" x2=\"600\" y2=\"300\"/>" << endl;
    output << "<line stroke=\"red\" stroke-width=\"1px\" x1=\"300\" y1=\"0\" x2=\"300\" y2=\"600\"/>" << endl;

    // Pfeile
    output << "<polygon fill=\"red\" points=\"296 8 304 8 300 0\" />" << endl;
    output << "<polygon fill=\"red\" points=\"592 296 600 300 592 304\" />" << endl;

    // Achsenbeschriftungen
    output << "<text x=\"529\" y=\"317\" fill=\"red\" font-size=\"15px\">In-Phase</text>" << endl;
    output << "<text x=\"303\" y=\"17\" fill=\"red\" font-size=\"15px\">Quadrature</text>" << endl;

    // Raster - horizontal
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"50\" x2=\"600\" y2=\"50\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"150\" x2=\"600\" y2=\"150\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"250\" x2=\"600\" y2=\"250\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"350\" x2=\"600\" y2=\"350\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"450\" x2=\"600\" y2=\"450\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"0\" y1=\"550\" x2=\"600\" y2=\"550\"/>" << endl;

    // Raster - vertikal
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"50\" y1=\"0\" x2=\"50\" y2=\"600\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"150\" y1=\"0\" x2=\"150\" y2=\"600\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"250\" y1=\"0\" x2=\"250\" y2=\"600\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"350\" y1=\"0\" x2=\"350\" y2=\"600\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"450\" y1=\"0\" x2=\"450\" y2=\"600\"/>" << endl;
    output << "<line stroke=\"grey\" stroke-width=\"1px\" x1=\"550\" y1=\"0\" x2=\"550\" y2=\"600\"/>" << endl;

    // Horizontal Rasterlinien-Beschriftung
    output << "<text x=\"5\" y=\"48\" fill=\"grey\" font-size=\"10px\">" << (sigamp/5*5) << "V</text>" << endl;
    output << "<text x=\"5\" y=\"148\" fill=\"grey\" font-size=\"10px\">" << (sigamp/5*3) << "V</text>" << endl;
    output << "<text x=\"5\" y=\"248\"  fill=\"grey\" font-size=\"10px\">" << (sigamp/5*1) << "V</text>" << endl;
    output << "<text x=\"5\" y=\"298\" fill=\"grey\" font-size=\"10px\"> 0 V</text>" << endl;
    output << "<text x=\"5\" y=\"348\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/5*1) << "V</text>" << endl;
    output << "<text x=\"5\" y=\"448\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/5*3) << "V</text>" << endl;
    output << "<text x=\"5\" y=\"548\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/5*5) << "V</text>" << endl;

    // Vertikale Rasterlinien-Beschriftung
    output << "<text x=\"552\" y=\"597\" fill=\"grey\" font-size=\"10px\">" << (sigamp/5*5) << "V</text>" << endl;
    output << "<text x=\"452\" y=\"597\" fill=\"grey\" font-size=\"10px\">" << (sigamp/5*3) << "V</text>" << endl;
    output << "<text x=\"352\" y=\"597\"  fill=\"grey\" font-size=\"10px\">" << (sigamp/5*1) << "V</text>" << endl;
    output << "<text x=\"302\" y=\"597\" fill=\"grey\" font-size=\"10px\"> 0 V</text>" << endl;
    output << "<text x=\"252\" y=\"597\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/5*1) << "V</text>" << endl;
    output << "<text x=\"152\" y=\"597\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/5*3) << "V</text>" << endl;
    output << "<text x=\"52\" y=\"597\" fill=\"grey\" font-size=\"10px\">-" << (sigamp/5*5) << "V</text>" << endl;

    // Infotext
    output << "<text x=\"5\" y=\"17\" font-size=\"15px\">Scatter-Plot: \""<< name() <<"\"</text>" << endl;
  }

void scatter::set_attributes()
{
  inI.set_rate(data_rate);	// Liest je 10 Token
  inQ.set_rate(data_rate);
}

void scatter::processing()
{
  // Normierung auf 1
  double sigI=inI.read(0)/sigamp; // Durch read(0) wird nur jeder 10. Wert eingelesen
  double sigQ=inQ.read(0)/sigamp;

  // Pixelberechnung
  int x=(int)(300 + sigI * 250);
  int y=(int)(300 + sigQ * 250);

  // Augabepuffer mit Koordinatenpaaren fuellen
  ss1 << x << " " << y << ", ";
}

void scatter:: finish()
  {
    string s=ss1.str(); // Erstellen des Polyline-Strings mit anschließender Ausgabe
    output << "<polyline fill=\"none\" stroke=\"black\" stroke-width=\"1px\"" << endl;
    output << "points=\"" << s << "\" />" << endl;
    output << "</svg>";

    // Abschluss des SVG-Dokuments und schließen des Filehandles
    output.close();
  }

/************************** Jitter ***************************/
void jit::set_attributes(){}

void jit::initialize()
     {
        num(0) = 1.0;
        den(0) = 1.0;
	 }

void jit:: processing() {
        
        double rnd1;
        double rnd2;
        double G1;
        double Q;
        double Q1;
        double Q2;
        double tmp;
        
	    do{
	        rnd1 = ((double)rand()) / ((double)RAND_MAX) ;
	        rnd2 = ((double)rand()) / ((double)RAND_MAX) ;

	        Q1 = 2.0 * rnd1 - 1.0 ;
	        Q2 = 2.0 * rnd2 - 1.0 ;

	        Q = Q1 * Q1 + Q2 * Q2 ;

	    } while (Q > 1.0) ;

	    G1 = sqrt( - 2.0 * log(Q) / Q) * Q1 ;
        tmp = max_jit + (max_jit/3) * G1;
        if(tmp > 2*max_jit)
            tdelay = sca_core::sca_time(2*max_jit, sc_core::SC_SEC);
        else if(tmp<0)
            tdelay =sca_core::sca_time(0, sc_core::SC_SEC);
        else
            tdelay =sca_core::sca_time(tmp, sc_core::SC_SEC);

        out.write( ltf_1( num,den,tdelay,state,in.read() ));
    }

jit::jit(sc_core::sc_module_name nm, double _max_jit) 
        {
            max_jit=_max_jit;
        }


  }//bb
}//TUV_ams_lib
