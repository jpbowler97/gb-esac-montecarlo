package gb.esac.montecarlo;


import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.binner.BinningException;
import gb.esac.binner.BinningUtils;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.periodogram.WindowFunctionException;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.BasicStats;
import gb.esac.tools.Complex;
import gb.esac.tools.FFT;
import java.text.DecimalFormat;
import java.util.Date;
import org.apache.log4j.Logger;



public final class TimmerKonig {

    private static Logger logger  = Logger.getLogger(TimmerKonig.class);
    private static DecimalFormat exp = new DecimalFormat("#.####E0#");
    private static DecimalFormat num = new DecimalFormat("#.####");

    /**
     * Method <code>generateComponentsFromSpectrum</code>  
     * is the core implementation of the Timmer-Konig algorithm.
     * It generates the randomized Fourier components from the input spectrum
     *
     * @param spec a <code>double[]</code> array that defines spectrum
     * @return a <code>Complex[]</code> array containing the (un-scaled) generated Fourier components
     */
    private static Complex[] generateComponentsFromSpectrum(double[] spec) {

	MersenneTwister64 engine = new MersenneTwister64(new Date());
	Normal gaussian = new Normal(0, 1, engine);
	double gauss1 = 0;
	double gauss2 = 0;
	double re, im = 0;
	Complex[] complex = new Complex[2*spec.length];
	for ( int i=0; i < spec.length; i++ ) {

	    //  Positive frequency
	    gauss1 = gaussian.nextDouble();
	    gauss2 = gaussian.nextDouble();
	    re = gauss1*Math.sqrt(0.5*spec[i]);
	    im = gauss2*Math.sqrt(0.5*spec[i]);
	    complex[i] = new Complex(re, im);
	    
	    // Timmer & Konig take the conjugate of the positive freq component
	    //complex[complex.length-1 -i] = complex[i].conjugate();

	    //  Negative frequency
	    gauss1 = gaussian.nextDouble();
	    gauss2 = gaussian.nextDouble();
	    re = gauss1*Math.sqrt(0.5*spec[i]);
	    im = gauss2*Math.sqrt(0.5*spec[i]);
	    complex[complex.length-1-i] = new Complex(re, im);
	}
	return complex;
    }

    private static Complex[] generateComponents(double[] frequencies, double alpha) {

	double[] spec = new double[frequencies.length];
	for ( int i=0; i < frequencies.length; i++ ) {
	    double omega = frequencies[i] * 2*Math.PI;
	    spec[i] = Math.pow(omega, -alpha);
	}
	return generateComponentsFromSpectrum(spec);
    }

    private static Complex[] generateComponents(double[] frequencies, double alpha1, double alpha2, double nuBreak) {

	//  Generate the two components of the spectrum
	double[] spec = new double[frequencies.length];
	int k=0;
	int l=0;
	for ( int i=0; i < frequencies.length; i++ ) {
	    double omega = frequencies[i] * 2*Math.PI;
	    if ( frequencies[i] < nuBreak ) {
		spec[i] = Math.pow(omega, -alpha1);
		k++;
	    }
	    else {
		spec[i] = Math.pow(omega, -alpha2);
		l++;
	    }
	}

	//  Adjust second component to join with first
	double omegaAtBreak = 2*Math.PI*frequencies[k];
 	double normFactor = Math.pow(omegaAtBreak, -alpha1)/Math.pow(omegaAtBreak, -alpha2);
 	for ( int i=0; i < l; i++ ) {
 	    spec[k+i] *= normFactor;
 	}

	return generateComponentsFromSpectrum(spec);
    }

    private static Complex[] generateComponents(double[] frequencies, double alpha1, double alpha2, double alpha3, double nuBreak1, double nuBreak2) {

	//  Generate the three components of the spectrum
	double[] spec = new double[frequencies.length];
	int k=0;
	int l=0;
	int m=0;
	double lastTerm = 0;
	for ( int i=0; i < frequencies.length; i++ ) {
	    double omega = frequencies[i] * 2*Math.PI;
	    if ( frequencies[i] < nuBreak1 ) {
		spec[i] = Math.pow(omega, -alpha1);
		k++;
	    }
	    else if ( frequencies[i] < nuBreak2 ) {
		spec[i] = Math.pow(omega, -alpha2);
		l++;
	    }
	    else {
		spec[i] = Math.pow(omega, -alpha3);
		m++;
	    }
	}

	//  Adjust second to join with first
	double omegaAtBreak1 = 2*Math.PI*frequencies[k];
 	double normFactor1 = Math.pow(omegaAtBreak1, -alpha1)/Math.pow(omegaAtBreak1, -alpha2);
 	for ( int i=0; i < l; i++ ) {
 	    spec[k+i] *= normFactor1;
 	}

	//  Adjust third to join with second
	double omegaAtBreak2 = 2*Math.PI*frequencies[k+l];
 	double normFactor2 = normFactor1*Math.pow(omegaAtBreak2, -alpha2)/Math.pow(omegaAtBreak2, -alpha3);
 	for ( int i=0; i < m; i++ ) {
 	    spec[k+l+i] *= normFactor2;
 	}

	//  Plot the spectrum
// 	try {
// 	    String[] h = new String[] {"DEV /XS", "LOG ON"};
// 	    AsciiDataFileWriter out = new AsciiDataFileWriter("spec.qdp");
// 	    out.writeData(h, frequencies, spec);
// 	}
// 	catch ( Exception e ) {};
	
	return generateComponentsFromSpectrum(spec);
    }
    
    /** 
    	IMPORTANT:

	nTimeBins = 2* nSpecBins
	
	A light curve defined on 64 bins allows for the testing of 32 independent frequencies.
	Since the FFT yields results for the 32 IFS and their negative, 
	we must suply 64 complex number to the FFT.
    **/


    /**  getFourierComponents  **/
    public static Complex[] getFourierComponents(double mean, double duration, double alpha) throws BinningException {
		
	return getFourierComponents(mean, duration, alpha, 1);
    }

    public static Complex[] getFourierComponents(double mean, double duration, double alpha, int nFreqsPerIFS) throws BinningException {

	double nuMin = 1d/duration;
	double nuMax = 2d*mean;
	return getFourierComponentsForFrequencies(alpha, nuMin, nuMax, nFreqsPerIFS);
    }


    /**  getFourierComponentsForFrequencies  **/
    public static Complex[] getFourierComponentsForFrequencies(double alpha, double nuMin, double nuMax) throws BinningException {

	return getFourierComponentsForFrequencies(alpha, nuMin, nuMax, 1);
    }

    public static Complex[] getFourierComponentsForFrequencies(double alpha, double nuMin, double nuMax, int nFreqsPerIFS) throws BinningException {

	double df = nuMin/nFreqsPerIFS;
	return getFourierComponentsForFrequencies(alpha, nuMin, nuMax, df);
    }

    public static Complex[] getFourierComponentsForFrequencies(double alpha, double nuMin, double nuMax, double df) throws BinningException {

	double[] frequencies = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, df);
	return generateComponents(frequencies, alpha);
    }

    public static Complex[] getFourierComponentsForFrequencies(double[] frequencies, double alpha) throws BinningException {

	return generateComponents(frequencies, alpha);
    }

    public static Complex[] getFourierComponentsForFrequencies(double[] frequencies, double alpha1, double alpha2, double nuBreak) {

	return generateComponents(frequencies, alpha1, alpha2, nuBreak);
    }

    public static Complex[] getFourierComponentsForFrequencies(double[] frequencies, double alpha1, double alpha2, double alpha3, double nuBreak1, double nuBreak2) {

	return generateComponents(frequencies, alpha1, alpha2, alpha3, nuBreak1, nuBreak2);
    }



    /**
     * The method <code>getRates</code> calls getFourierComponents and inverse Fourier transforms to get the associated rates that are then normalized to the value of meanRate. The normalization is done by first dividing the rates by their actual mean, and then multiplying by the specified mean. This ensures that the fractional RMS is preserved. This is important since the variance of the rates is a function of the spectral index of the red noise. Hence it must be preserved.
     *
     *
     * @param index a <code>double</code> value
     * @param duration a <code>double</code> value
     * @param nTimeBins an <code>int</code> value
     * @return a <code>double[]</code> value
     * @exception BinningException if an error occurs
     */
    public static double[] getRates(double mean, double duration, double alpha) throws BinningException {

	return getRates(mean, duration, alpha, 1);
    }

    public static double[] getRates(double mean, double duration, double alpha, int nFreqsPerIFS) throws BinningException {

	double nuMin = 1/duration;
	double nuMax = 2*mean;
	double df = nuMin/nFreqsPerIFS;
	double nFreqs = (nuMax - nuMin)/df;

	//  Adjust nuMax to have power of 2 number of frequencies
 	double exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
	int nNewBins = (int) Math.pow(2, exponent);
	if ( nFreqs != nNewBins ) {
	    logger.warn("Number of specified frequencies ("+nFreqs+") not a power of 2. Using "+nNewBins+" instead");
	    nFreqs = nNewBins;
	}
	nuMax = nuMin + df*nFreqs;
	
	//  Generate Fourier components, get the rates, and scale them to the specified mean rate
	Complex[] fourierComp = getFourierComponentsForFrequencies(alpha, nuMin, nuMax, df);
	double[] rates = getRatesFromFourierComponents(fourierComp);
	double scalingFactor = mean/BasicStats.getMean(rates);
	for ( int i=0; i < rates.length; i++ ) {
	    rates[i] *= scalingFactor;
	}
	return rates;
    }

    public static double[] getRatesFromFourierComponents(Complex[] fourierComp) throws BinningException {
	
	int n = fourierComp.length;
	if ( ! isPowerOfTwo(n) ) {
	    throw new BinningException("Number of Fourier components ("+n+") is not a power of 2.");
	}

	//  Inverse Fourier transform the components to get the light curve
	Complex[] ifft = FFT.ifft(fourierComp);

	//  Taking the norm of each complex number gets 
	//  all the info contained in the result of the inverse FFT
	//  and ensures that there are no negative numbers.
	double[] rates = Complex.norm(ifft);

	return rates;
    }

    private static boolean isPowerOfTwo(int n) {

	boolean isPowerOfTwo = true;
	double n1 = Math.floor(Math.log10(n)/Math.log10(2));
	double n2 = Math.log10(n)/Math.log10(2);
	if ( n1 != n2 ) {
	    isPowerOfTwo = false;
	}
	return isPowerOfTwo;
    }


    public static TimeSeries getTimeSeries(double mean, double duration, double alpha) throws BinningException {

	double[] rates = getRates(mean, duration, alpha, 1);
	return getTimeSeries(rates, duration);
    }

    public static TimeSeries getTimeSeries(double mean, double duration, double alpha, int nFreqsPerIFS) throws BinningException {

	double[] rates = getRates(mean, duration, alpha, nFreqsPerIFS);
	return getTimeSeries(rates, duration);
    }


    private static TimeSeries getTimeSeries(double[] rates, double duration) throws BinningException {
	
	int nTimeBins = rates.length;
	double[] binEdges  = BinningUtils.getBinEdges(0, duration, nTimeBins);
	double binWidth = duration/nTimeBins;
	double[] counts = new double[rates.length];
	for ( int i=0; i < rates.length; i++ ) {
	    counts[i] = rates[i]*binWidth;
	}
	return TimeSeriesMaker.makeTimeSeries(binEdges, counts);
    }


    /**
     * The method <code>getPowers</code> calls getFourierComponents and returns the corresponding powers.
     *
     * @param index a <code>double</code> value
     * @param duration a <code>double</code> value
     * @param nTimeBins an <code>int</code> value
     * @return a <code>double[]</code> value
     * @exception BinningException if an error occurs
     */
    public static double[] getPowers(double mean, double duration, double alpha, int nTimeBins) throws BinningException {
	
	Complex[] fourierComp = getFourierComponents(mean, duration, alpha, nTimeBins);
	return getPowers(fourierComp);
    }
	
    public static double[] getPowers(Complex[] fourierComp) {

	return Complex.normSquared(fourierComp);
    }


    private static Complex[] scale(Complex[] fourierComponents, double scalingFactor) {

	Complex[] result = new Complex[fourierComponents.length];
	for ( int i=0; i < result.length; i++ ) {
	    result[i] = fourierComponents[i].times(scalingFactor);
	}
	return result;
    }


}
