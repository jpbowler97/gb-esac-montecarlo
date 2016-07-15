package gb.esac.montecarlo;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.binner.BinningUtils;
import gb.esac.binner.Resampler;;
import gb.esac.eventlist.EventList;
import gb.esac.periodogram.AggregatePeriodogram;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesResampler;
import gb.esac.tools.Converter;
import gb.esac.tools.DistributionFunc;
import hep.aida.ref.histogram.Histogram1D;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;


public class TimeMappingEffect {

    static DecimalFormat sci = new DecimalFormat("0.0#E0");

    public static void main(String[] args) throws Exception {

	//  Arguments
	double mean = 1;
	double duration = 1e3;
	double alpha = 0.5;
	int n = 100;
	if ( args.length == 4 ) {
	    mean = (Double.valueOf(args[0])).doubleValue();
	    duration = (Double.valueOf(args[1])).doubleValue();
	    alpha = (Double.valueOf(args[2])).doubleValue();
	    n = (Integer.valueOf(args[3])).intValue();
	}

	//  Define number of Timmer time bins
	double effectiveNyquistBinTime = 1/(2*mean);
	effectiveNyquistBinTime += Math.ulp(effectiveNyquistBinTime);
	double nbins = Math.ceil(duration/effectiveNyquistBinTime);
 	double exponent = Math.floor(Math.log10(nbins)/Math.log10(2));
	int nTimeBins = (int) Math.pow(2, exponent);
	double tkBinTime = duration/nTimeBins;

	//  Define the number of events to draw
	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(1, engine);
	int nevents = (new Double(mean*duration)).intValue();
	nevents = poisson.nextInt(nevents);
	
	//  Define periodogram aggregates
	AggregatePeriodogram avgShort = new AggregatePeriodogram();
	AggregatePeriodogram avgMed = new AggregatePeriodogram();
	AggregatePeriodogram avgLong = new AggregatePeriodogram();

	//  Loop on n
	for ( int i=0; i < n; i++ ) {

// 	    double[] times = RedNoiseGenerator.generateArrivalTimes(mean, duration, alpha);
// 	    EventList evlist = new EventList(times);
// 	    FFTPeriodogram psd = PeriodogramMaker.makePlainFFTPeriodogram(evlist);
// 	    avgShort.add(psd);

// 	    times = RedNoiseGenerator.generateArrivalTimes(mean, 10*duration, alpha);
// 	    evlist = new EventList(times);
// 	    psd = PeriodogramMaker.makePlainFFTPeriodogram(evlist);
// 	    avgMed.add(psd);

// 	    times = RedNoiseGenerator.generateArrivalTimes(mean, 100*duration, alpha);
// 	    evlist = new EventList(times);
// 	    psd = PeriodogramMaker.makePlainFFTPeriodogram(evlist);
// 	    avgLong.add(psd);


	    int [] resamplingFactors = new int[] {1, 16, 128};
	    double[] longDurations = new double[] {duration, 16*duration, 128*duration};
	    for ( int j=0; j < resamplingFactors.length; j++ ) {
		int factor = resamplingFactors[j];
		double longDuration = longDurations[j];

		//  Generate rates from freqs based on T=longDuration with matching time resolution
		double[] rates = TimmerKonig.getTimmerRates(alpha, longDuration, nTimeBins*factor);

		//  Resample the rates so that there are the same number of bins as if T=duration
// 		double[] oldBinEdges = BinningUtils.getBinEdges(0, longDuration, factor*nTimeBins);
// 		double[] newBinEdges = BinningUtils.getBinEdges(0, longDuration, nTimeBins);
// 		double[] timmerRates = Resampler.resample(rates, oldBinEdges, newBinEdges);
 		double[] timmerRates = rates;

		//  Map these rates onto a time line of T=duration
		double tzero = 0;
		tkBinTime = duration/timmerRates.length;
		Histogram1D lcHisto = Converter.array2histo("light curve", tzero, tkBinTime, timmerRates);
		Histogram1D cdfHisto = DistributionFunc.getCDFHisto(lcHisto);
		double[] times = DistributionFunc.getRandom(cdfHisto, nevents);
		Arrays.sort(times);
		times[nevents-1] = times[0] + duration;  //  Adjust actual duration to specified duration

		//  Make FFTPeriodogram and add to AggregatePeriodogram 
		EventList evlist = new EventList(times);
		FFTPeriodogram psd = PeriodogramMaker.makePlainFFTPeriodogram(evlist);
		if ( factor == resamplingFactors[0] ) {
		    avgLong.add(psd);
		}
		else if ( factor == resamplingFactors[1] ) {
		    avgMed.add(psd);
		}
		else {
		    avgShort.add(psd);
		}
	    }

	}
	AveragePeriodogram psdShort = avgShort.getPeriodogram();
	AveragePeriodogram psdMed = avgMed.getPeriodogram();
	AveragePeriodogram psdLong = avgLong.getPeriodogram();

	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("psd-mappingEffect3-mu-"+mean+"-T-"+sci.format(duration)+"-alpha-"+alpha+".qdp")));
	String[] header2 = new String[] {
		"DEV /XS",
		"READ SERR 2",
		"LAB T", "LAB F",
		"TIME OFF",
		"LINE STEP",
		"LOG ON",
		"LW 4", "CS 1.5",
		"LAB X Frequency (Hz)",
		"LAB Y Power",
		"VIEW 0.2 0.1 0.8 0.9",
		"SKIP SINGLE",
		"ERR OFF",
		"!"
	};
	for ( int i=0; i < header2.length; i++ ) 
	    pw.println(header2[i]);

	double[] freqs = psdLong.getFreqs();
	double[] powers = psdLong.getPowers();
	double[] errors = psdLong.getErrors();
	int i=0;
	while ( i < freqs.length && freqs[i] <= 0.25 ) {
	    pw.println(freqs[i]+"\t"+powers[i]+"\t"+errors[i]);
	    i++;
	}
	pw.println("NO NO NO");

	freqs = psdMed.getFreqs();
	powers = psdMed.getPowers();
	errors = psdMed.getErrors();
	i=0;
	while ( i < freqs.length && freqs[i] <= 0.25 ) {
	    pw.println(freqs[i]+"\t"+powers[i]+"\t"+errors[i]);
	    i++;
	}
	pw.println("NO NO NO");

	freqs = psdShort.getFreqs();
	powers = psdShort.getPowers();
	errors = psdShort.getErrors();
	i=0;
	while ( i < freqs.length && freqs[i] <= 0.25 ) {
	    pw.println(freqs[i]+"\t"+powers[i]+"\t"+errors[i]);
	    i++;
	}
	pw.close();

    }

}
