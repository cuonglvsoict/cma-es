package functions;

public class Sphere extends Function {

	@Override
	public double getValue(double[] x) {
		double v[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i] - bias[i];
		}

		double sum = 0;
		for (double d : v) {
			sum += d * d;
		}
		return sum;

	}

}
