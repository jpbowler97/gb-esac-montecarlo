package gb.esac.montecarlo;


import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.engine.RandomEngine;
import gb.esac.binner.BinningException;
import gb.esac.binner.BinningUtils;
import gb.esac.binner.Resampler;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.timeseries.TimeSeriesException;
import gb.esac.tools.Complex;
import gb.esac.tools.Converter;
import gb.esac.tools.DistributionFunc;
import hep.aida.ref.histogram.Histogram1D;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Date;
import org.apache.log4j.Logger;



/**
 * The class <code>RedNoiseGenerator</code> is used to simulate red noise. The term "red noise" refers to a signal that has the most power at the lowest frequencies, i.e., the "red" part of the spectrum. More precisely, the power spectrum of red noise follows a power-law with a negative index. The index determines the "redness" of the noise: the higher the value the redder the noise. At the lower limit, red noise with a spectral index of 0 is, in fact, white noise.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0
 */
public final class RedNoiseGenerator {

    static Logger logger = Logger.getLogger(RedNoiseGenerator.class);
    static DecimalFormat sci = new DecimalFormat("0.0#E00");


    public static double[] generateArrivalTimes(double meanRate, double duration, double alpha) throws BinningException, TimeSeriesException {
		
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date()); // initialises random number generator called engine
	return generateArrivalTimes(meanRate, duration, alpha, engine);
    }

    public static double[] generateArrivalTimes(double meanRate, double duration, double alpha, RandomEngine engine) throws BinningException, TimeSeriesException {

	return generateArrivalTimes(meanRate, duration, alpha, 1, engine);
    }


    public static double[] generateArrivalTimes(double meanRate, double duration, double alpha, int nFreqsPerIFS) throws BinningException, TimeSeriesException {
		
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	return generateArrivalTimes(meanRate, duration, alpha, nFreqsPerIFS, engine);
    }
    
    public static double[] generateArrivalTimes(double meanRate, double duration, double alpha, int nFreqsPerIFS, RandomEngine engine) throws BinningException, TimeSeriesException {

	double nuMin = 1d/duration;
	double nuMax = 2d*meanRate; // 
	nuMax = meanRate; // 
	double df = nuMin/nFreqsPerIFS;
	double nFreqs = (nuMax - nuMin)/df;

	//  Adjust to have a power of 2 bins
 	double exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
	int nNewBins = (int) Math.pow(2, exponent);
	if ( nFreqs != nNewBins ) {
	    nNewBins = (int) Math.pow(2, exponent+1);
	    nFreqs = nNewBins;
	}
	int nTimeBins = (new Double(2*nFreqs)).intValue();

	return generateArrivalTimes(meanRate, duration, alpha, nFreqsPerIFS, nTimeBins, engine);
    }

    private static double[] generateArrivalTimes(double meanRate, double duration, double alpha, int nFreqsPerIFS, int nTimeBins, RandomEngine engine) throws BinningException, TimeSeriesException {

	if ( alpha == 0 ) {
	    return WhiteNoiseGenerator.generateArrivalTimes(meanRate, duration, engine);
	}

	DecimalFormat number = new DecimalFormat("0.000");
	DecimalFormat freq = new DecimalFormat("0.00#E00");

	logger.info("Generating red noise arrival times");
	logger.info("Mean rate (specified) = "+meanRate);
	logger.info("Duration (specified) = "+sci.format(duration));
	logger.info("Spectral index (specified) = "+alpha);

	//  Define effective Nyquist bin time
	double nuMin = 1/duration;
	double df = nuMin;
	double nuMax = 2*meanRate;
	double minBinTime = 1d/nuMax;
	int nIFS = (int) Math.floor( (nuMax-nuMin)/nuMin );
	int nBinsMin = (int) Math.ceil(duration/minBinTime);
	logger.info("Minimum frequency and IFS width (nuMin=df=1/duration) = "+freq.format(nuMin)+" Hz");
	logger.info("Effective Nyquist frequency (nuMax=2*meanRate) = "+freq.format(nuMax)+" Hz");
	logger.info("Effective Nyquist bintime (1/nuMax) = "+minBinTime+" s");
	logger.info("This results in "+nIFS+" IFS, and requires a minimum of "+nBinsMin+" bins for (effective) Nyquist resolution");
	
	//  Redefine nBinsMin to be the closest power of 2 larger than nBins
 	double exponent = Math.ceil(Math.log10(nBinsMin)/Math.log10(2));
	int nBinsMinClosestPowerOfTwo = (int) Math.pow(2, exponent);
	if ( nBinsMin != nBinsMinClosestPowerOfTwo ) {
	    logger.warn("Number of required bins ("+nBinsMin+") is not a power of 2. Using "+nBinsMinClosestPowerOfTwo+" instead");
	    nBinsMin = nBinsMinClosestPowerOfTwo;
	}

	//  Redefine the specified nTimeBins that is closest power of 2 
	logger.info("Specified number of Bins = "+nTimeBins);
 	exponent = Math.floor(Math.log10(nTimeBins)/Math.log10(2));
	int nNewBins = (int) Math.pow(2, exponent);
	if ( nTimeBins != nNewBins ) {
	    logger.warn("Number of specified bins ("+nTimeBins+") is not a power of 2. Using "+nNewBins+" time bins instead");
	    nTimeBins = nNewBins;
	}
	double actualBinTime = duration/nTimeBins;
	logger.info("Pure red noise time-domain signal will be defined on "+nTimeBins+" bins"); 
	logger.info("Thus the inherent time resolution of the signal is "+number.format(actualBinTime)+" s");

	//  Get the Timmer-Konig rates
 	double[] timmerRates = TimmerKonig.getRates(meanRate, duration, alpha, nFreqsPerIFS);

	//  Compare values of nBinsMin and nTimeBins
	if ( nBinsMin != nTimeBins ) {
	    logger.warn("Closest power of 2 greater than specified number of time bins ("+nTimeBins+") not equal to number of bins required for effective Nyquist resolution ("+nBinsMin+")");
	}
// 	if ( nTimeBins < nBinsMin ) {
// 	    logger.warn(nTimeBins+" < "+nBinsMin+": This would cause loss of resolution: Resampling.");
// 	    double[] oldBinEdges = BinningUtils.getBinEdges(0, duration, timmerRates.length);
// 	    double[] newBinEdges = BinningUtils.getBinEdges(0, duration, nTimeBins);
// 	    timmerRates = Resampler.resample(timmerRates, oldBinEdges, newBinEdges);	    
// 	}

	//  Define the number of events
	Poisson poisson = new Poisson(0, engine);
	int nevents = (new Double(meanRate*duration)).intValue();
	nevents = poisson.nextInt(nevents);

	//  Draw arrival times from the CDF of rates
	double tzero = 0;
	double tkBinTime = duration/nTimeBins;
	Histogram1D lcHisto = Converter.array2histo("light curve", tzero, tkBinTime, timmerRates);
	Histogram1D cdfHisto = DistributionFunc.getCDFHisto(lcHisto);
	double[] times = DistributionFunc.getRandom(cdfHisto, nevents);
	Arrays.sort(times);


	//  Adjust actual duration to specified duration
 	nevents=times.length;
 	times[nevents-1] = times[0] + duration;
	double actualMean = times.length/duration;
	
	logger.info("Arrival times generated");
	logger.info("nEvents = "+times.length);
	logger.info("Mean rate (actual) = "+actualMean);
	

	// START OF SECTION FOR TESTING

	//  	    /**  Write the full resolution Timmer LC and power spectrum  **/
	// 	    LightCurve timmerLC = TimmerKonig.getLightCurve(timmerRates, duration);
	// 	    timmerLC.scale(meanRate);
	// 	    timmerLC.writeAsQDPWithoutErrors("timmerLC.qdp");
	
	// 	    PeriodogramMaker psdMaker = new PeriodogramMaker();
	// 	    FFTPeriodogram timmerPSD = psdMaker.makeFFTPeriodogram(timmerLC, "leahy");
	// 	    timmerPSD.rebin(10, "papadakis");
	// 	    timmerPSD.addConstant(2.0);
	// 	    timmerPSD.writeAsQDP("timmerPSD.qdp");
	
	
	
	//              /**  Write CDF  **/
	// 	     LightCurve cdfLC = new LightCurve(lcTimes, cdf);
	// 	     cdfLC.writeAsQDP("cdf.qdp");
	
	
	// 	    private static DecimalFormat dec = new DecimalFormat("0.00");			
	//  	    /**  Rebin the Timmer LC using a bintime close to 40 s **/
	// 	    double expo = Math.floor(Math.log10(duration/40)/Math.log10(2));
	// 	    int nNewBins = (new Double(Math.pow(2, expo))).intValue();
	// 	    double nuMax = 0.5*nNewBins/duration;
	// 	    double lcBinTime = duration/nNewBins;
	// 	    System.out.println("Log  : New bintime for rebinning = "+number.format(lcBinTime));
	// 	    double[] binnedLC_TK = Binner.rebinRatesSimple(timmerRates, tkBinTime, lcBinTime);	    
	// 	    double binnedLC_TK_mean = BasicStats.getMean(binnedLC_TK);
	// 	    double binnedLC_TK_var = BasicStats.getVariance(binnedLC_TK);
	// 	    double binnedLC_TK_fracRMS = Math.sqrt(binnedLC_TK_var)/binnedLC_TK_mean;
	
	
	// 	    lcTimes = new double[binnedLC_TK.length];
	// 	    for ( int i=0; i < binnedLC_TK.length; i++ ) {
	// 		lcTimes[i] = (0.5 + i)*lcBinTime;
	// 	    }
	
	
	// 	    /**  Define output files and initialise output streams  **/
	// 	    String outfile_lc = "redNoise_lightcurves";
	// 	    String outfile_psd = "redNoise_powspecs";
	// 	    String outfile_lc_poisson = "redNoise_lightcurves_poisson";
	// 	    String outfile_psd_poisson = "redNoise_powspecs_poisson";
	// 	    PrintWriter pw_lc = null;
	// 	    PrintWriter pw_psd = null;
	// 	    PrintWriter pw_lc_poisson = null;
	// 	    PrintWriter pw_psd_poisson = null;
	// 	    try { 
	// 		pw_lc = new PrintWriter(new BufferedWriter(new FileWriter(outfile_lc+".qdp")));
	// 		pw_psd = new PrintWriter(new BufferedWriter(new FileWriter(outfile_psd+".qdp")));
	// 		pw_lc_poisson = new PrintWriter(new BufferedWriter(new FileWriter(outfile_lc_poisson+".qdp")));
	// 		pw_psd_poisson = new PrintWriter(new BufferedWriter(new FileWriter(outfile_psd_poisson+".qdp")));
	// 	    }
	// 	    catch ( IOException e ) { System.out.println("Error: Cannot create PrintWriter");}
	
	
	// 	    /**  Write header for lightcurves  **/
	// 	    head = new String[] {
	// 		"! QDP File",
	// 		"DEV /XS", "READ SERR 2",
	// 		"TIME OFF", "LINE ON",
	// 		"LAB T", "LAB F", "LW 3", "CS 0.9",
	// 		"LAB Y Arbitrary",
	// 		"LAB Y2 Cts/s",
	// 		"LAB Y3 Cts/s",
	// 		"LAB Y4 Cts/s",
	// 		"LAB Y5 Cts/s",
	// 		"VIEW 0.25 0.1 0.75 0.95",
	// 		"SKIP SINGLE",
	// 		"PLOT VERT",
	// 		"LAB X Time (s)",
	//  		"R X "+(-0.04*duration)+" "+(duration*1.04),
	// 		"!"
	// 	    };
	// 	    for ( int i=0; i < head.length; i++ ) {
	// 		pw_lc.println(head[i]);
	// 		pw_lc_poisson.println(head[i]);
	// 	    }
	
	// 	    /**  Write the rebinned Timmer LC  **/
	// 	    for ( int i=0; i < nNewBins; i++ ) {
	// 		pw_lc.println(lcTimes[i]+"\t"+binnedLC_TK[i]);
	// 		pw_lc_poisson.println(lcTimes[i]+"\t"+binnedLC_TK[i]);
	// 	    }
	// 	    pw_lc.println("NO NO NO");
	// 	    pw_lc_poisson.println("NO NO NO");
	
	
	// 	    /**  Make the PSD of the rebinned Timmer LC **/
	// 	    double[][] psd = Periodograms.makePSD_FFT(lcTimes, binnedLC_TK, "leahy-like");
	// 	    double[] f = psd[0];
	// 	    double[] p = psd[1];
	// 	    double[] logPow = Converter.lin2logSpace(p);
	// 	    double[] logFreq = Converter.lin2logSpace(f);
	
	
	// 	    /**  Define variable for loop on different rates **/
	// 	    double[] rates = new double[] {200, 20, 2, 0.2};
	// 	    double[] lc_err = new double[nNewBins];
	// 	    double[] lc_rate = new double[nNewBins];
	// 	    double[] scaledPoissonLC_rate = new double[nNewBins];
	// 	    double[] scaledPoissonLC_err = new double[nNewBins];
	// 	    double[] binnedLC = new double[nNewBins];
	// 	    double[] ymin = new double[rates.length];
	// 	    double[] ymax = new double[rates.length];
	// 	    double maxErr, max, min, mean, yPosDiff, yNegDiff = 0;
	// 	    double[] alphas = new double[rates.length+1];
	// 	    double[] alphasErr = new double[rates.length+1];
	// 	    double[] intercept = new double[rates.length+1];
	// 	    double[] alphas_poiss = new double[rates.length+1];
	// 	    double[] alphasErr_poiss = new double[rates.length+1];
	// 	    double[] intercept_poiss = new double[rates.length+1];
	
	
	// 	    /**  Fit the Timmer powspec in log-space  **/
	// 	    double[] fitResult = LeastSquaresFitter.leastSquaresFitLine(logFreq, logPow);
	// 	    alphas[0] = -fitResult[1];
	// 	    alphasErr[0] = fitResult[3];
	// 	    intercept[0] = fitResult[0];
	
	// 	    alphas_poiss[0] = -fitResult[1];
	// 	    alphasErr_poiss[0] = fitResult[3];
	// 	    intercept_poiss[0] = fitResult[0];
	
	
	// 	    /**  Define the header for the powspecs  **/
	// 	    head = new String[] {
	// 		"! QDP File",
	// 		"DEV /XS", "READ 2",
	// 		"TIME OFF", "LINE ON",
	// 		"LAB T", "LAB F", "LW 3", "CS 0.9",
	// 		"LAB Y Arbitrary",
	// 		"LAB X Frequency (Hz)",
	// 		"LOG X",
	// 		"VIEW 0.25 0.1 0.75 0.95",
	// 		"SKIP DOUBLE",
	// 		"PLOT VERT",
	//  		"R X "+(nuMin*0.8)+" "+(nuMax*1.2),
	// 		"!"
	// 	    };


	// 	    /**  Write the Timmer powspec **/
	// 	    for ( int i=0; i < head.length; i++ ) {
	// 		pw_psd.println(head[i]);
	// 		pw_psd_poisson.println(head[i]);
	// 	    }
	// 	    for ( int i=0; i < f.length; i++ ) {
	// 		pw_psd.println(f[i]+"\t"+logPow[i]);
	// 		pw_psd_poisson.println(f[i]+"\t"+logPow[i]);
	// 	    }
	// 	    pw_psd.println("NO NO");
	// 	    pw_psd_poisson.println("NO NO");
	// 	    for ( int i=0; i < f.length; i++ ) {
	// 		double fittedPow = intercept[0] -alphas[0]*logFreq[i];
	// 		pw_psd.println(f[i]+"\t"+fittedPow);
	// 		pw_psd_poisson.println(f[i]+"\t"+fittedPow);
	// 	    }
	// 	    pw_psd.println("NO NO"+"\n"+"NO NO");
	// 	    pw_psd_poisson.println("NO NO"+"\n"+"NO NO");


	// 	    /**  Loop on rates  **/
	// 	    for ( int r=0; r < rates.length; r++ ) {

	// 		/**  Pick arrival times from the Timmer LC  **/
	// 		int n = poisson.nextInt(rates[r]*duration);
	// 		double[] ts = Analysis.getRandom(cdfHisto, n);
	// 		Arrays.sort(ts);
	// 		double[] ts_adjusted = adjustTimes(ts, duration);

	// 		/**  bin them  **/
	// 		binnedLC = Binner.binOrderedData(ts_adjusted, nNewBins);

	// 		/**  Calculate count rate and errors in each bin  **/
	// 		for ( int m=0; m < binnedLC.length; m++  ) {
	// 		    lc_rate[m] = binnedLC[m]/lcBinTime;
	// 		    lc_err[m] = Math.sqrt(binnedLC[m])/lcBinTime;
	// 		}
	// 		double lc_evlist_mean = BasicStats.getMean(lc_rate);
	// 		double lc_evlist_var = BasicStats.getVariance(lc_rate);


	// 		/**  Scale the simulated light curve to have identical fractional rms as original **/
	// 		scaledPoissonLC_rate = new double[nNewBins];
	// 		for ( int k=0; k < binnedLC.length; k++ ) {
	// 		    scaledPoissonLC_rate[k] = binnedLC_TK[k]*rates[r];
	// 		    //scaledPoissonLC_rate[k] = (binnedLC_TK[k] - binnedLC_TK_mean)*Math.sqrt(lc_evlist_var/binnedLC_TK_var);
	// 		}


	// 		/**  Scale to the same mean as the event list light curve  **/
	// 		double scaledPoissonLC_mean = BasicStats.getMean(scaledPoissonLC_rate);
	// 		for ( int k=0; k < binnedLC.length; k++ ) {
	// 		    //scaledPoissonLC_rate[k] = (scaledPoissonLC_rate[k] - scaledPoissonLC_mean) + lc_evlist_mean;
	//    		    double eventsInBin = scaledPoissonLC_rate[k]*lcBinTime;
	//    		    int events = poisson.nextInt(eventsInBin);
	//    		    scaledPoissonLC_rate[k] = events/lcBinTime;
	// 		    //System.out.println(eventsInBin+"\t"+events+"\t"+scaledPoissonLC_rate[k]);
	// 		}

	// 		/**  Rescale the variance and mean to those of the event list light curve  **/
	// 		double scaledPoissonLC_var = BasicStats.getVariance(scaledPoissonLC_rate);
	// // 		for ( int k=0; k < binnedLC.length; k++ ) {
	// // 		    scaledPoissonLC_rate[k] *= Math.sqrt(lc_evlist_var/scaledPoissonLC_var);
	// // 		}
	// // 		scaledPoissonLC_mean = BasicStats.getMean(scaledPoissonLC_rate);
	// 		for ( int k=0; k < binnedLC.length; k++ ) {
	// // 		    scaledPoissonLC_rate[k] = (scaledPoissonLC_rate[k] - scaledPoissonLC_mean) + lc_evlist_mean;
	// 		    scaledPoissonLC_err[k] = Math.sqrt(scaledPoissonLC_rate[k]/lcBinTime);
	// 		}

	// 		System.out.print("V = "+number.format(lc_evlist_var)+"\t"+number.format(BasicStats.getVariance(scaledPoissonLC_rate)));
	// 		System.out.println("\t E = "+number.format(lc_evlist_mean)+"\t"+number.format(BasicStats.getMean(scaledPoissonLC_rate)));


	// 		/**  Write the LC  **/
	// 		maxErr = MinMax.getMax(lc_err);
	// 		max = MinMax.getMax(lc_rate) + maxErr;
	// 		min = MinMax.getMin(lc_rate) - maxErr;
	// 		mean = BasicStats.getMean(lc_rate);
	// 		yPosDiff = max - mean;
	// 		yNegDiff = mean - min;
	// 		ymax[r] =  mean + 1.1*(yPosDiff);
	// 		ymin[r] =  mean - 1.1*(yNegDiff);
	// 		for ( int m=0; m < binnedLC.length; m++ ) {
	// 		    pw_lc.println(lcTimes[m]+"\t"+lc_rate[m]+"\t"+lc_err[m]);
	// 		    pw_lc_poisson.println(lcTimes[m]+"\t"+scaledPoissonLC_rate[m]+"\t"+scaledPoissonLC_err[m]);
	// 		}
	// 		pw_lc.println("NO NO NO");
	// 		pw_lc_poisson.println("NO NO NO");


	// 		/**  Event list
	// 		     Make Powspec and fit in log-space  **/
	// 		psd = Periodograms.makePSD_FFT(lcTimes, lc_rate, "leahy-like");
	// 		f = psd[0];
	// 		p = psd[1];
	// 		logPow = Converter.lin2logSpace(p);
	// 		logFreq = Converter.lin2logSpace(f);
	// 		fitResult = LeastSquaresFitter.leastSquaresFitLine(logFreq, logPow);
	// 		alphas[r+1] = -fitResult[1];
	// 		alphasErr[r+1] = fitResult[3];
	// 		intercept[r+1] = fitResult[0];

	// 		/**  Write powspec in with power in log-space  **/
	// 		for ( int m=0; m < f.length; m++ ) {
	// 		    pw_psd.println(f[m]+"\t"+logPow[m]);
	// 		}
	// 		pw_psd.println("NO NO");
	// 		for ( int m=0; m < f.length; m++ ) {
	// 		    double fittedPow = intercept[r+1] -alphas[r+1]*logFreq[m];
	// 		    pw_psd.println(f[m]+"\t"+fittedPow);
	// 		}
	// 		pw_psd.println("NO NO"+"\n"+"NO NO");

	// 		/**  Scale poisson
	// 		     Make Powspec and fit in log-space  **/
	// 		psd = Periodograms.makePSD_FFT(lcTimes, scaledPoissonLC_rate, "leahy-like");
	// 		f = psd[0];
	// 		p = psd[1];
	// 		logPow = Converter.lin2logSpace(p);
	// 		logFreq = Converter.lin2logSpace(f);
	// 		fitResult = LeastSquaresFitter.leastSquaresFitLine(logFreq, logPow);
	// 		alphas_poiss[r+1] = -fitResult[1];
	// 		alphasErr_poiss[r+1] = fitResult[3];
	// 		intercept_poiss[r+1] = fitResult[0];

	// 		/**  Write powspec in with power in log-space  **/
	// 		for ( int m=0; m < f.length; m++ ) {
	// 		    pw_psd_poisson.println(f[m]+"\t"+logPow[m]);
	// 		}
	// 		pw_psd_poisson.println("NO NO");
	// 		for ( int m=0; m < f.length; m++ ) {
	// 		    double fittedPow = intercept_poiss[r+1] -alphas_poiss[r+1]*logFreq[m];
	// 		    pw_psd_poisson.println(f[m]+"\t"+fittedPow);
	// 		}
	// 		pw_psd_poisson.println("NO NO"+"\n"+"NO NO");

	// 	    }

	// 	    /**  Write the ranges for the lightcurves  **/
	// 	    pw_lc.println("R Y");
	// 	    pw_lc_poisson.println("R Y");
	// 	    for ( int r=0; r < rates.length; r++ ) {
	// 		pw_lc.println("R Y"+(r+2)+" "+number.format(ymin[r])+" "+number.format(ymax[r]));
	// 		pw_lc_poisson.println("R Y"+(r+2)+" "+number.format(ymin[r])+" "+number.format(ymax[r]));
	// 	    }
	// 	    pw_lc.println("!"+"\n"+"HARD "+outfile_lc+".ps/ps");
	// 	    pw_lc.flush();
	// 	    pw_lc.close();
	// 	    pw_lc_poisson.println("!"+"\n"+"HARD "+outfile_lc_poisson+".ps/ps");
	// 	    pw_lc_poisson.flush();
	// 	    pw_lc_poisson.close();


	// 	    /**  Write the results of the fits of alphas for the powspecs  **/
	// 	    String xPos = "0.28";
	// 	    String[] yPos = new String[] {"0.82", "0.65", "0.48", "0.31", "0.14"};
	// 	    for ( int r=0; r < alphas.length; r++ ) {
	// 		pw_psd.println("LAB "+(r+1)+" VPOS "+xPos+" "+yPos[r]+" \"\\ga = "+
	//        dec.format(alphas[r])+" +/- "+dec.format(alphasErr[r])+"\" JUST LEFT");
	// 		pw_psd_poisson.println("LAB "+(r+1)+" VPOS "+xPos+" "+yPos[r]+" \"\\ga = "+
	//        dec.format(alphas_poiss[r])+" +/- "+dec.format(alphasErr_poiss[r])+"\" JUST LEFT");
	// 	    }
	// 	    for ( int r=0; r < alphas.length; r++ ) {
	// 		pw_psd.println("LAB Y"+(r+2)+" \"Log Power\"");
	// 		pw_psd_poisson.println("LAB Y"+(r+2)+" \"Log Power\"");
	// 	    }
	// 	    pw_psd.println("!"+"\n"+"HARD "+outfile_psd+".ps/ps");
	// 	    pw_psd.flush();
	// 	    pw_psd.close();
	// 	    pw_psd_poisson.println("!"+"\n"+"HARD "+outfile_psd_poisson+".ps/ps");
	// 	    pw_psd_poisson.flush();
	// 	    pw_psd_poisson.close();


	// END    TEMPORARY FOR TESTING PURPOSES

	return times;

    }
	
	
    public static double[] generateModulatedArrivalTimes(double meanRate, double duration, double alpha, double period, double pulsedFrac) throws BinningException,  TimeSeriesException {

	int nFreqsPerIFS = 1;
	return generateModulatedArrivalTimes(meanRate, duration, alpha, nFreqsPerIFS, period, pulsedFrac);
    }

    public static double[] generateModulatedArrivalTimes(double meanRate, double duration, double alpha, int nFreqsPerIFS, double period, double pulsedFrac) throws BinningException, TimeSeriesException {
		
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	return generateModulatedArrivalTimes(meanRate, duration, alpha, nFreqsPerIFS, period, pulsedFrac, engine);
    }
	
    public static double[] generateModulatedArrivalTimes(double meanRate, double duration, double alpha, double period, double pulsedFrac, RandomEngine engine) throws BinningException, TimeSeriesException {
	
	int nFreqsPerIFS = 1;
	return generateModulatedArrivalTimes(meanRate, duration, alpha, nFreqsPerIFS, period, pulsedFrac, engine);
    }
	
    public static double[] generateModulatedArrivalTimes(double meanRate, double duration, double alpha, int nFreqsPerIFS, double period, double pulsedFrac, RandomEngine engine) throws BinningException, TimeSeriesException {
		
	if ( alpha == 0.0 ) {
	    return WhiteNoiseGenerator.generateModulatedArrivalTimes(meanRate, duration, period, pulsedFrac, engine);
	}

	logger.info("Generating sinusoidally modulated red noise arrival times");
	logger.info("Mean rate (specified) = "+meanRate);
	logger.info("Duration = "+duration);
	logger.info("Spectral index = "+alpha);
	logger.info("Period = "+period);
	logger.info("Pulsed Fraction = "+pulsedFrac);

	Poisson poisson = new Poisson(1, engine);
	int nevents = (new Double(meanRate*duration)).intValue();
	logger.info("Nominal number of events = "+nevents+" (duration*meanRate)");
	nevents = poisson.nextInt(nevents);
	logger.info("Random Poisson number of events = "+nevents);
	double nEvents = (new Double(nevents)).doubleValue();
	double nPulsedEvents = pulsedFrac*nEvents;
	double pulsedMeanRate = nPulsedEvents/duration;
	logger.info("Pulsed events = "+nPulsedEvents);
	logger.info("Pulsed mean rate = "+pulsedMeanRate);
	double nRedNoiseEvents = nEvents - nPulsedEvents;
	double redNoiseMeanRate = nRedNoiseEvents/duration;
	logger.info("Red noise events = "+nRedNoiseEvents);
	logger.info("Red noise mean rate = "+redNoiseMeanRate);
	double[] pulsedTimes = WhiteNoiseGenerator.generateModulatedArrivalTimes(pulsedMeanRate, duration, period, 0.999, engine);
	double[] redNoiseTimes = generateArrivalTimes(redNoiseMeanRate, duration, alpha, engine);

	logger.info("Combining red noise with sinusoidally modulated arrival times");

	if ( pulsedTimes.length <= 2 )
	    return redNoiseTimes;
	else {

	    int nTimes = pulsedTimes.length + redNoiseTimes.length -1;
	    double[] times = new double[nTimes];
	    int n = 0;
	    for ( int i=0; i < pulsedTimes.length-1; i++ ) {
		times[i] = pulsedTimes[i];
		n++;
	    }
	    for ( int i=0; i < redNoiseTimes.length; i++ ) {
		times[i+n] = redNoiseTimes[i];
	    }
	    Arrays.sort(times);

	    //  Adjust actual duration to specified duration
	    times[times.length-1] = times[0] + duration;
	    double actualMean = times.length/duration;

	    logger.info("Arrival times in list = "+times.length);
	    logger.info("Mean rate (actual) = "+actualMean);

	    return times;
	}
    }


    public static double[] generateTwoComponentArrivalTimes(double meanRate, double duration, double alpha1, double alpha2, double nuBreak) throws BinningException, TimeSeriesException {

	int nFreqsPerIFS = 1;
	return generateTwoComponentArrivalTimes(meanRate, duration, alpha1, alpha2, nuBreak, nFreqsPerIFS);
    }

    public static double[] generateTwoComponentArrivalTimes(double meanRate, double duration, double alpha1, double alpha2, double nuBreak, int nFreqsPerIFS) throws BinningException, TimeSeriesException {

	logger.info("Generating two-component red noise arrival times");
	logger.info("Mean rate (specified) = "+meanRate+" cps");
	logger.info("Duration = "+duration+" s");
	logger.info("Spectral index1 = "+alpha1);
	logger.info("Spectral index2 = "+alpha2);
	logger.info("Break frequency = "+nuBreak+" Hz");

	double nuMin = 1d/duration;
	double nuMax = 2d*meanRate;
	double df = nuMin/nFreqsPerIFS;
	double nFreqs = (nuMax - nuMin)/df;

	//  Adjust nuMax to have 2^x frequencies
 	double exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
	int nNewBins = (int) Math.pow(2, exponent);
	if ( nFreqs != nNewBins ) {
	    nNewBins = (int) Math.pow(2, exponent+1);
	    logger.warn("nFreqs ("+nFreqs+") was not a power of 2. Using "+nNewBins+" instead");
	    nFreqs = nNewBins;
	}
	nuMax = nuMin + df*(nFreqs);
	double[] frequencies = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, df);
	double powAtNuMin = meanRate*duration;
	Complex[] fourierComponents = TimmerKonig.getFourierComponentsForFrequencies(frequencies, alpha1, alpha2, nuBreak);
	double[] timmerRates = TimmerKonig.getRatesFromFourierComponents(fourierComponents);

// 	// TEMP
// 	int[] bins = new int[timmerRates.length];
// 	for ( int i=0; i < timmerRates.length; i++ ) bins[i] = i;
// 	String[] h = new String[] {"DEV /XS"};
// 	AsciiDataFileWriter out = null;
// 	try {
// 	    out = new AsciiDataFileWriter("timmerRates.qdp");
// 	    out.writeData(h, bins, timmerRates);
// 	}
// 	catch ( Exception e ) {};
// 	// 

	//  Define the number of events
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(0, engine);
	int nevents = (new Double(meanRate*duration)).intValue();
	logger.info("Nominal number of events = "+nevents+" (duration*meanRate)");
	nevents = poisson.nextInt(nevents);
	logger.info("Random Poisson number of events = "+nevents);

	//  Draw arrival times from the light curve's CDF
	double tzero = 0;
	double nTimeBins = 2*nFreqs;
	double tkBinTime = duration/nTimeBins;
	Histogram1D lcHisto = Converter.array2histo("light curve", tzero, tkBinTime, timmerRates);
	Histogram1D cdfHisto = DistributionFunc.getCDFHisto(lcHisto);
	double[] times = DistributionFunc.getRandom(cdfHisto, nevents);
	Arrays.sort(times);

// 	// TEMP
// 	try {
// 	    AsciiDataFileWriter lc = new AsciiDataFileWriter("lcHisto.qdp");
// 	    lc.writeHisto(lcHisto, "LC");
// 	}
// 	catch ( Exception e ) {}
// 	//

	//  Adjust actual duration to specified duration
	times[times.length-1] = times[0] + duration;
	double actualMean = times.length/duration;
	
	logger.info("Arrival times in list = "+times.length);
	logger.info("Mean rate (actual) = "+actualMean);

	return times;
    }

    public static double[] generateThreeComponentArrivalTimes(double meanRate, double duration, double alpha1, double alpha2, double alpha3, double nuBreak1, double nuBreak2) throws BinningException, TimeSeriesException {

	int nFreqsPerIFS = 1;

	logger.info("Generating two-component red noise arrival times");
	logger.info("Mean rate (specified) = "+meanRate);
	logger.info("Duration = "+sci.format(duration));
	logger.info("Spectral index1 = "+alpha1);
	logger.info("Spectral index2 = "+alpha2);
	logger.info("Spectral index3 = "+alpha3);
	logger.info("Break frequency1 = "+nuBreak1);
	logger.info("Break frequency2 = "+nuBreak2);

	double nuMin = 1d/duration;
	double nuMax = 2*meanRate;
	double df = nuMin/nFreqsPerIFS;
	double nFreqs = (nuMax - nuMin)/df;

	//  Adjust nuMax to have 2^x frequencies
 	double exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
	int nNewBins = (int) Math.pow(2, exponent);
	if ( nFreqs != nNewBins ) {
	    nNewBins = (int) Math.pow(2, exponent+1);
	    logger.warn("nFreqs ("+nFreqs+") was not a power of 2. Using "+nNewBins+" instead");
	    nFreqs = nNewBins;
	}
	nuMax = nuMin + df*nFreqs;
	double[] frequencies = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, df);
	double powAtNuMin = meanRate*duration;
	Complex[] fourierComponents = TimmerKonig.getFourierComponentsForFrequencies(frequencies, alpha1, alpha2, alpha3, nuBreak1, nuBreak2);
	double[] timmerRates = TimmerKonig.getRatesFromFourierComponents(fourierComponents);

	//  Define the number of events
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(0, engine);
	int nevents = (new Double(meanRate*duration)).intValue();
	logger.info("Nominal number of events = "+nevents+" (duration*meanRate)");
	nevents = poisson.nextInt(nevents);
	logger.info("Random Poisson number of events = "+nevents);

	//  Draw arrival times from the light curve's CDF
	double tzero = 0;
	double nTimeBins =  2*nFreqs;
	double tkBinTime = duration/nTimeBins;
	Histogram1D lcHisto = Converter.array2histo("light curve", tzero, tkBinTime, timmerRates);
	Histogram1D cdfHisto = DistributionFunc.getCDFHisto(lcHisto);
	double[] times = DistributionFunc.getRandom(cdfHisto, nevents);
	Arrays.sort(times);

	//  Adjust actual duration to specified duration
	times[nevents-1] = times[0] + duration;
	double actualMean = times.length/duration;
	
	logger.info("nEvents = "+times.length);
	logger.info("Mean rate (actual) = "+actualMean);

	return times;
    }

	
}
