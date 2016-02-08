package old;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

class Covariance{
	public static void main(String[] args){
		double[] one = {0,0,0,0,0};
		double[] two = {0,0,0,0,0};
		double res = new PearsonsCorrelation().correlation(one, two);
		
		System.out.println(res);
	}
}