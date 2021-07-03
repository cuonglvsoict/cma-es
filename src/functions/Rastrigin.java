package functions;

public class Rastrigin extends Function {

	@Override
	public double getValue(double[] x) {
		double v[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i] - bias[i];
		}

		int dim = x.length;
		double obj = dim * 10;
		for (int i = 0; i < dim; i++) {
			obj = obj + (v[i] * v[i] - 10 * (Math.cos(2 * Math.PI * v[i])));
		}
		return obj;
	}

}
