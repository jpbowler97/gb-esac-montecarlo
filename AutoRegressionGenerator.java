package gb.esac.montecarlo;


import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.tools.MinMax;
import java.text.DecimalFormat;
import java.util.Date;
import org.apache.log4j.Logger;

public final class AutoRegressionGenerator {

    static Logger logger = Logger.getLogger(RedNoiseGenerator.class);
    static DecimalFormat sci = new DecimalFormat("0.0#E0");


    
    /**
     * Method <code>generateAR1Process</code> 
     * The autocorrelation coefficient is labelled as rho
     *
     * @param order an <code>int</code> value of the order of the AR process
     * @param tau a <code>double</code> value of the characteristic timescale of the AR process
     * @param samplingTimes a <code>double[]</code> array specifying the arbitrarily spaced sampling times
     * @return a <code>double[]</code> value
     */
    public static double[] generateAR1Process(double[] samplingTimes, double tau) {

	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
 	Normal randomGauss = new Normal(0, 1, engine);
	int nTimes = samplingTimes.length;
	double[] amplitudes = new double[nTimes];
	amplitudes[0] = 0;
	for ( int i=1; i < nTimes; i++ ) {
	    double rho = Math.exp( -(samplingTimes[i] - samplingTimes[i-1])/tau );
	    double noiseVariance = -Math.exp( -2*(samplingTimes[i]-samplingTimes[i-1])/tau );
	    double noiseTerm = randomGauss.nextDouble(0, noiseVariance);
	    amplitudes[i] = rho*amplitudes[i-1] + noiseTerm;
	}
	double min = MinMax.getMin(amplitudes);
	min *= -1;
	for ( int i=0; i < nTimes; i++ ) {
	    amplitudes[i] += min;
	    amplitudes[i] /= 10;
	}
	return amplitudes;
    }


}
