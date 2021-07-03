import java.util.Random;

import basic.Individual;
import basic.Params;
import functions.Function;
import functions.Sphere;
import utils.Utils;

public class Main {
	public static void main(String args[]) {
		// setup benchmark
		Function func;
		
		func = new Sphere();
		func.dim = 50;
		func.LB = Utils.ones(func.dim, -100);
		func.UB = Utils.ones(func.dim, 100);
		func.bias = new double[func.dim];
		func.matrix = Utils.i_matrix(func.dim);
		for (int i = 0; i < func.dim; i++) {
			func.bias[i] = 0;
		}


//		func = new Rastrigin();
//		func.dim = 50;
//		func.LB = Utils.ones(func.dim, -50);
//		func.UB = Utils.ones(func.dim, 50);
//		func.bias = new double[func.dim];
//		func.matrix = Utils.i_matrix(func.dim);
//		int m = func.dim / 2;
//		for (int i = 0; i < m; i++) {
//			func.bias[i] = 40;
//		}
//		for (int i = m; i < func.dim; i++) {
//			func.bias[i] = -40;
//		}
		
		Params.rand = new Random(0);
		int lambda = (int) (4 + 3 * Math.floor(Math.log(func.dim)));
		CMAES solver = new CMAES(func.dim, func, lambda, false);
		solver.initialize();
		Individual best = solver.run();
		
		System.out.println("Best solution found: " + best.fitness);
	}
}
