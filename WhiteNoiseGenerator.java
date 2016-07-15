package gb.esac.montecarlo;


import java.text.DecimalFormat;
import java.util.Arrays;

import cern.colt.list.DoubleArrayList;
import cern.jet.random.engine.RandomEngine;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.Exponential;
import cern.jet.random.Poisson;
import cern.jet.random.Normal;
import cern.jet.random.Uniform;

import org.apache.log4j.Logger;



/**
 * The class <code>WhiteNoiseGenerator</code> is used to simulate white noise.  The term "white noise" refers to a signal that has equal power at all frequencies, and thus the power spectrum of white noise is flat. White noise is uncorrelated or memory-less noise, which means that each event occurs independently, with no memory of the occurence time of the previous event: this is a Poisson process. As such, the distribution of inter-arrival times follows the single parameter exponential distribution with mean given by the inverse of the mean count rate (the decay constant is equal to the mean count rate). The methods of this class generate white noise in the form of arrival times, either with a constant mean or a sinusoidally modulated mean.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (May 2010, ESAC)
 */
public final class WhiteNoiseGenerator {


    static Logger logger = Logger.getLogger(WhiteNoiseGenerator.class);


    private static DecimalFormat number = new DecimalFormat("0.00000000");
    private static DecimalFormat freq = new DecimalFormat("0.000E0");

	

   /**
     * Generate white noise arrival times with constant mean.
     *
     * @param meanRate a <code>double</code> value
     * @param duration a <code>double</code> value
     * @return a <code>double[]</code> value
     */
    public static double[] generateArrivalTimes (double meanRate, double duration) {
		
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	return generateArrivalTimes(meanRate, duration, engine);
    }
	


    /**
     * Generate white noise arrival times with constant mean.
     *
     * @param meanRate a <code>double</code> value
     * @param duration a <code>double</code> value
     * @param engine a <code>RandomEngine</code> value
     * @return a <code>double[]</code> value
     */
    public static double[] generateArrivalTimes (double meanRate, double duration, RandomEngine engine) {

	logger.info("Generating white noise arrival times with constant mean rate");
	logger.info("Duration = "+duration);
	logger.info("Mean rate (specified) = "+meanRate);

	
	//  Generate arrival times
	Exponential randomExp = new Exponential(meanRate, engine);
	DoubleArrayList times = new DoubleArrayList();
	double dt = randomExp.nextDouble();
	double tzero = dt;
	double time = dt;
	times.add(time);
	while ( time < duration ) {

	    dt = randomExp.nextDouble();

	    time += dt;
	    times.add(time);
	}
	times.trimToSize();
	double[] arrivalTimes = times.elements();


	//  Adjust actual duration to specified duration
	arrivalTimes[arrivalTimes.length-1] = tzero + duration;
	int nevents = arrivalTimes.length;
	double mean = nevents/duration;

	logger.info("Arrival times generated");
	logger.info("nEvents = "+nevents);
	logger.info("Mean rate (actual) = "+mean);

	return arrivalTimes;
    }



	
    /**
     * Generate sinusoidally modullated white nosie arrival times. 
     *
     * @param meanRate a <code>double</code> value
     * @param duration a <code>double</code> value

     * @param period a <code>double</code> value
     * @param pulsedFrac a <code>double</code> value
     * @return a <code>double[]</code> value
     */
    public static double[] generateModulatedArrivalTimes(double meanRate, double duration, double period, double pulsedFrac) {
		
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	return generateModulatedArrivalTimes(meanRate, duration, period, pulsedFrac, engine);
    }

	
    /**
     * Generate sinusoidally modullated white nosie arrival times. 
     *
     * @param meanRate a <code>double</code> value
     * @param duration a <code>double</code> value
     * @param period a <code>double</code> value
     * @param pulsedFrac a <code>double</code> value
     * @param engine a <code>RandomEngine</code> value
     * @return a <code>double[]</code> value
     */
    public static double[] generateModulatedArrivalTimes(double meanRate, double duration, double period, double pulsedFrac, RandomEngine engine) {

	logger.info("Generating sinusoidally modulated white noise arrival times");
	logger.info("Mean rate (specified) = "+meanRate);
	logger.info("Duration = "+duration);
	logger.info("Period = "+period);
	logger.info("Pulsed Fraction = "+pulsedFrac);
		
	/**     Generating modulated arrival times
		 
	    Instantaneous modulated flux is given by:	
	    lambda(t) = a + b*(1 + sin(wt))
	    where  a = meanRate/( 1 + pulsedFrac/(wT)*(1 - cos(wT)) )
	           b = a*pulsedFrac
	    Probability density function is given by:
	    f(t) = [a + b*(1 - sin(w(t0 + t)))] * exp{ -at - (2b/w)(sin(w(t0 + t/2))*sin(wt/2)) }
		 
	**/

	//  Set up random generator distributions
	Uniform randUni = new Uniform(engine);
	Exponential randExp = new Exponential(1.0, engine);
	Poisson randPoiss = new Poisson(1.0, engine);
	//  Define Variables
	double pf = pulsedFrac;
	if ( pf == 1.0 ) pf = 0.999;
	double t_obs = duration;
	double w = 2*Math.PI/period;
	double phi_obs = w*t_obs;
	double lambda_0 = meanRate/( 1 + pf/phi_obs*(1 - Math.cos(phi_obs)) );
	double lambda_1 = pf*lambda_0;
	double thetaMax = lambda_0*t_obs + (lambda_1/w)*(1 - Math.cos(phi_obs));
	//  Construct set of theta[i]'s from which we will get the t[i]'s
	int nevents = randPoiss.nextInt(thetaMax);
	double[] theta = new double[nevents];
	for ( int i=0; i < nevents; i++ ) {
	    theta[i] = randUni.nextDouble()*thetaMax;
	}
	Arrays.sort(theta);
	//  Construct the t[i]'s using each theta[i] by solving
	//   theta[i] = lambda_0*t[i] - lambda_1/w * Math.cos(w*t[i])
	double[] times = new double[nevents];
	double t_new = 0;
	double f = 1;
	double fPrime = 0;
	//  Root finding by the False Position method
	for ( int i=0; i < nevents; i++ ) {
	    double a = -1;
	    double b = t_obs;
	    double t = 0;
	    int side = 0;
	    double fOfa = lambda_0*a + lambda_1/w*(1-Math.cos(w*a)) - theta[i];
	    double fOfb = lambda_0*b + lambda_1/w*(1-Math.cos(w*b)) - theta[i];
	    double fOft = lambda_0*t + lambda_1/w*(1-Math.cos(w*t)) - theta[i];
	    double delta = 1e-5;
	    double absDiff = Math.abs(b-a);
	    double limit = delta*Math.abs(b+a);
	    int iterations = 0;
	    while ( absDiff > limit && i < nevents ) {
		iterations++;
		t = (fOfa*b - fOfb*a) / (fOfa - fOfb);
		absDiff = Math.abs(b-a);
		limit = delta*Math.abs(b+a);
		fOft = lambda_0*t + lambda_1/w*(1-Math.cos(w*t)) - theta[i];
		if ( fOft*fOfb > 0 ) {
		    b = t;
		    fOfb = fOft;
		    if ( side == -1 ) {
			fOfa /= 2;
		    }
		    side = -1;
		}
		else if ( fOfa * fOft > 0 ) {
		    a = t;
		    fOfa = fOft;
		    if ( side == 1 ) {
			fOfb /= 2;
		    }
		    side = 1;
		}
		else break;
	    }
	    //System.out.println("t = "+number.format(t)+" (iterations = "+iterations+")");
	    times[i] = t;
	}
	//  Adjust actual duration to specified duration
	try {
	    times[nevents-1] = times[0] + duration;
	}
	catch (ArrayIndexOutOfBoundsException e) {
	    logger.warn("There are no elements in this array");
	}
	// 	/**  Root finding by Newton's method **/
	// 	for ( int i=0; i < nevents; i++ ) {
	// 	    double delta = 1;
	// 	    double t = theta[i]/lambda_0;
	// 	    //System.out.println("i = "+i+"(of "+nevents+"), \t First t = "+t);
	// 	    while ( delta > 1e-3 && i < nevents && f != 0 ) {
	// 		f = lambda_0*t + lambda_1/w*(1 - Math.cos(w*t)) - theta[i];
	// 		fPrime = lambda_0 + lambda_1*Math.sin(w*t);
	// 		t_new = t - f/fPrime;
	// 		//System.out.print("i = "+i+"(of "+nevents+"), \t (t-t_new) = "+number.format(t-t_new)+"\t");
	// 		delta = Math.abs((t - t_new)/t_new);
	// 		//System.out.println("delta = "+number.format(delta));
	// 		t = t_new;
	// 	    }
	// 	    //System.out.println("delta (upon exit) = "+delta);	    
	// 	    times[i] = t;
	// 	}
	double mean = nevents/duration;
	logger.info("Arrival times generated");
	logger.info("nEvents = "+nevents);
	logger.info("Mean rate (actual) = "+mean);
	return times;
    }



//     /**
//      * Generate white noise arrival times, each drawn from the distribution with the specified mean.
//      *
//      * @param meanRates a <code>double</code> value
//      * @param duration a <code>double</code> value
//      * @return a <code>double[]</code> value
//      */
//     public static double[] generateWhiteArrivalTimes (double[] meanRates, double duration) {
		
// 	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
// 	return generateWhiteArrivalTimes(meanRates, duration, engine);
//     }



//     /**
//      * Generate white noise arrival times, each drawn from the distribution with the specified mean.
//      *
//      * @param meanRates a <code>double</code> value
//      * @param duration a <code>double</code> value
//      * @param engine a <code>RandomEngine</code> value
//      * @return a <code>double[]</code> value
//      */
//     public static double[] generateWhiteArrivalTimes (double[] meanRates, double duration, RandomEngine engine) {


// 	logger.info("Generating white noise arrival times with variable mean rate");
// 	logger.info("Duration = "+duration);

// 	//  Generate arrival times
// 	Exponential randomExp = new Exponential(1, engine);
// 	DoubleArrayList times = new DoubleArrayList();
// 	int i=0;
// 	double dt = randomExp.nextDouble(meanRates[i]);
// 	double tzero = dt;
// 	double time = dt;
// 	times.add(time);
// 	while ( i < meanRates.length ) {

// 	    dt = randomExp.nextDouble(meanRates[i]);

// 	    time += dt;
// 	    times.add(time);
// 	    i++;
// 	}
// 	times.trimToSize();
// 	double[] arrivalTimes = times.elements();


// 	//  Adjust actual duration to specified duration
// 	arrivalTimes[arrivalTimes.length-1] = tzero + duration;
// 	int nevents = arrivalTimes.length;
// 	double meanRate = nevents/duration;

// 	logger.info("Arrival time have been generated");
// 	logger.info("nEvents = "+nevents);
// 	logger.info("Mean rate (actual) = "+meanRate);

// 	return arrivalTimes;
//     }
	

	
//     /**
//      * Generate white noise arrival times, each drawn from the distribution with the specified mean with specified gaussian standard deviation.
//      *
//      * @param meanRates a <code>double</code> value
//      * @param meanErrors a <code>double</code> value
//      * @param duration a <code>double</code> value
//      * @return a <code>double[]</code> value
//      */
//     public static double[] generateWhiteArrivalTimes(double[] meanRates, double[] meanErrors, double duration) {
		
// 	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
// 	return generateWhiteArrivalTimes(meanRates, meanErrors, duration, engine);
		
//     }

	

//     /**
//      * Generate white noise arrival times, each drawn from the distribution with the specified mean with specified gaussian standard deviation.
//      *
//      * @param meanRates a <code>double</code> value
//      * @param meanErrors a <code>double</code> value
//      * @param duration a <code>double</code> value
//      * @param engine a <code>RandomEngine</code> value
//      * @return a <code>double[]</code> value
//      */
//     public static double[] generateWhiteArrivalTimes(double[] meanRates, double[] meanErrors, double duration, RandomEngine engine) {
		

// 	logger.info("Generating white noise arrival times with variable mean rate");
// 	logger.info("Duration = "+duration);


// 	//  Generate arrival times
// 	Normal randomGauss = new Normal(0, 1, engine);
// 	Exponential randomExp = new Exponential(1, engine);
// 	DoubleArrayList times = new DoubleArrayList();
// 	int i=0;
// 	double mean = randomGauss.nextDouble(meanRates[i], meanErrors[i]);
// 	double dt = randomExp.nextDouble(mean);
// 	i++;
// 	double tzero = dt;
// 	double time = dt;
// 	times.add(time);
// 	while ( i < meanRates.length ) {

// 	    mean = randomGauss.nextDouble(meanRates[i], meanErrors[i]);
// 	    dt = randomExp.nextDouble(mean);

// 	    time += dt;
// 	    times.add(time);
// 	    i++;
// 	}
// 	times.trimToSize();
// 	double[] arrivalTimes = times.elements();


// 	//  Adjust actual duration to specified duration
// 	arrivalTimes[arrivalTimes.length-1] = tzero + duration;
// 	int nevents = arrivalTimes.length;
// 	double meanRate = nevents/duration;

// 	logger.info("Arrival time have been generated");
// 	logger.info("nEvents = "+nevents);
// 	logger.info("Mean rate (actual) = "+meanRate);

// 	return arrivalTimes;

//     }



}
