package gb.esac.montecarlo;

import gb.esac.eventlist.EventList;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;

import java.util.Arrays;

public class TestTwoComponentRedNoise {

    public static void main(String[] args) throws Exception {

	double mean = 10;
	double duration = 1e4;
	double alpha1 = 1;
	double alpha2 = 2.;
	double alpha3 = 1;
	double nuBreak1 = 1e-3;
	double nuBreak2 = 1e-2;

	//double[] times = RedNoiseGenerator.generateTwoComponentArrivalTimes(mean, duration, alpha1, alpha2, nuBreak1);
	//double[] times = RedNoiseGenerator.generateThreeComponentArrivalTimes(mean, duration, alpha1, alpha2, alpha3, nuBreak1, nuBreak2);

	EventList evlist = new EventList(times);
	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(evlist, 10d);
	ts.writeCountsAsQDP("ts-multiCompRedNoise.qdp");
	FFTPeriodogram fft = PeriodogramMaker.makeWindowedFFTPeriodogram(ts, "Blackman", "Leahy");
	fft.writeAsQDP("fft-multiCompRedNoise.qdp");

    }

}
