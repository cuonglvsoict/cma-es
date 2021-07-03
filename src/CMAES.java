
import basic.Individual;
import basic.Params;
import basic.Population;
import functions.Function;
import utils.Utils;

public class CMAES {

	/**
	 * Number of offspring generated in each iteration.
	 */
	private int lambda;

	/**
	 * Number of individual selected for recombination.
	 */
	private int mu;

	/**
	 * overall standard deviation (step size).
	 */
	private double sigma;

	/**
	 * variance effectiveness
	 */
	private double mueff;

	/**
	 * learning rate
	 */
	private double ccov;

	/**
	 * learning rate when diag mode is active
	 */
	private double ccovsep;

	/**
	 * Expectation of ||N(0, I)||.
	 */
	private double chiN;

	/**
	 * Step size cumulation parameter.
	 */
	private double cs;

	/**
	 * Cumulation parameter.
	 */
	private double cc;

	/**
	 * Damping for step size.
	 */
	private double damps;

	/**
	 * Weights for recombination.
	 */
	private double[] weights;

	/**
	 * Scaling factors.
	 */
	private double[] diagD;

	/**
	 * Current centroid of the distribution.
	 */
	private double[] xmean;

	/**
	 * Evolution path.
	 */
	private double[] pc;

	/**
	 * Conjugate evolution path.
	 */
	private double[] ps;

	/**
	 * Coordinate system.
	 */
	private double[][] B;

	/**
	 * Current covariance matrix.
	 */
	private double[][] C;

	/**
	 * The current population.
	 */
	private Population population;

	/**
	 * Last iteration were the eigenvalue decomposition was calculated.
	 */
	private int lastEigenupdate;

	/**
	 * Best solution found so far
	 */
	private Individual best;

	/**
	 * search space dimension
	 */
	private int N;

	/**
	 * iteration counter
	 */
	private int iteration;

	/**
	 * if true, perform consistency checks to ensure the algorithm remains
	 * numerically stable.
	 */
	private boolean checkConsistency;

	/**
	 * The number of iterations in which only the covariance diagonal is used. This
	 * enhancement helps speed up the algorithm when there are many decision
	 * variables. Set to {@code 0} to always use the full covariance matrix.
	 */
	private int diagonalIterations;

	public CMAES(int dim, Function func, int lambda, boolean checkConsistency) {
		this.N = dim;
		this.lambda = lambda;
		this.checkConsistency = checkConsistency;
		this.population = new Population(dim, func);
	}

	public void initialize() {
		sigma = 0.5;
		diagonalIterations = 150 * N / lambda;
		diagD = new double[N];
		pc = new double[N];
		ps = new double[N];
		B = new double[N][N];
		C = new double[N][N];

		for (int i = 0; i < N; i++) {
			pc[i] = 0;
			ps[i] = 0;
			diagD[i] = 1;

			for (int j = 0; j < N; j++) {
				B[i][j] = 0;
			}

			for (int j = 0; j < i; j++) {
				C[i][j] = 0;
			}

			B[i][i] = 1;
			C[i][i] = diagD[i] * diagD[i];
		}

		xmean = new double[N];
		for (int i = 0; i < N; i++) {
			xmean[i] = Params.rand.nextDouble();
		}

		chiN = Math.sqrt(N) * (1.0 - 1.0 / (4.0 * N) + 1.0 / (21.0 * N * N));
		mu = (int) Math.floor(lambda / 2.0);
		weights = new double[mu];

		for (int i = 0; i < mu; i++) {
			weights[i] = Math.log(mu + 1) - Math.log(i + 1);
		}

		double sum = Utils.sum(weights);

		for (int i = 0; i < mu; i++) {
			weights[i] /= sum;
		}

		double sumSq = Utils.sumSq(weights);

		mueff = 1.0 / sumSq;

		cs = (mueff + 2) / (N + mueff + 3);
		damps = (1 + 2 * Math.max(0, Math.sqrt((mueff - 1.0) / (N + 1)) - 1)) + cs;
		cc = 4.0 / (N + 4.0);
		ccov = 2.0 / (N + 1.41) / (N + 1.41) / mueff
				+ (1 - (1.0 / mueff)) * Math.min(1, (2 * mueff - 1) / (mueff + (N + 2) * (N + 2)));
		ccovsep = Math.min(1, ccov * (N + 1.5) / 3.0);

	}

	/**
	 * Performs eigenvalue decomposition to update B and diagD.
	 */
	public void eigenDecomposition() {
		lastEigenupdate = iteration;

		if (diagonalIterations >= iteration) {
			for (int i = 0; i < N; i++) {
				diagD[i] = Math.sqrt(C[i][i]);
			}
		} else {
			// set B <- C
			for (int i = 0; i < N; i++) {
				for (int j = 0; j <= i; j++) {
					B[i][j] = B[j][i] = C[i][j];
				}
			}

			// eigenvalue decomposition
			double[] offdiag = new double[N];
			Utils.tred2(N, B, diagD, offdiag);
			Utils.tql2(N, diagD, offdiag, B);

//			if (checkConsistency) {
//				Utils.checkEigenSystem(N, C, diagD, B);
//			}

			// assign diagD to eigenvalue square roots
			for (int i = 0; i < N; i++) {
				if (diagD[i] < 0) { // numerical problem?
					System.err.println("an eigenvalue has become negative");
					diagD[i] = 0;
				}

				diagD[i] = Math.sqrt(diagD[i]);
			}
		}
	}

	/**
	 * Test and correct any numerical issues.
	 */
	public void testAndCorrectNumerics() {
		// flat fitness, test is function values are identical
		if (population.individuals.size() > 0) {
			population.sort();

			if (population.individuals.get(0).fitness == population.individuals
					.get(Math.min(lambda - 1, lambda / 2 + 1) - 1).fitness) {
				System.err.println("flat fitness landscape, consider reformulation of fitness, step size increased");
				sigma *= Math.exp(0.2 + cs / damps);
			}
		}

		// align (renormalize) scale C (and consequently sigma)
		double fac = 1.0;

		if (Utils.max(diagD) < 1e-6) {
			fac = 1.0 / Utils.max(diagD);
		} else if (Utils.min(diagD) > 1e4) {
			fac = 1.0 / Utils.min(diagD);
		}

		if (fac != 1.0) {
			sigma /= fac;

			for (int i = 0; i < N; i++) {
				pc[i] *= fac;
				diagD[i] *= fac;

				for (int j = 0; j <= i; j++) {
					C[i][j] *= fac * fac;
				}
			}
		}
	}

	public Individual samplePopulation() {
		iteration++;
		Individual best = null;

		if ((iteration - lastEigenupdate) > 1.0 / ccov / N / 5.0) {
			eigenDecomposition();
		}

		if (checkConsistency) {
			testAndCorrectNumerics();
		}

		population.clear();
		if (this.diagonalIterations >= iteration) {
			for (int i = 0; i < this.lambda; i++) {
				double[] x = new double[N];
				for (int j = 0; j < N; j++) {
					do {
						x[j] = xmean[j] + sigma * diagD[j] * Params.rand.nextGaussian();
					} while (x[j] > 1 || x[j] < 0);
				}

				Individual indiv = new Individual(x);
				indiv.fitness = indiv.calcFitness(population.func);
				population.addIndividual(indiv);

				if (best == null || best.fitness > indiv.fitness) {
					best = indiv;
				}
			}
		} else {
			for (int i = 0; i < this.lambda; i++) {
				boolean feasible;
				double[] x = new double[N];

				do {
					feasible = true;
					double[] artmp = new double[N];
					for (int j = 0; j < N; j++) {
						artmp[j] = diagD[j] * Params.rand.nextGaussian();
					}

					for (int j = 0; j < N; j++) {
						double sum = 0.0;
						for (int k = 0; k < N; k++) {
							sum += B[j][k] * artmp[k];
						}

						x[j] = xmean[j] + sigma * sum;

						if (x[j] > 1 || x[j] < 0) {
							feasible = false;
							break;
						}
					}
				} while (!feasible);

				Individual indiv = new Individual(x);
				indiv.fitness = indiv.calcFitness(population.func);
				population.addIndividual(indiv);

				if (best == null || best.fitness > indiv.fitness) {
					best = indiv;
				}
			}
		}

		return best;
	}

	public void updateDistribution() {
		double[] xold = xmean.clone();
		double[] BDz = new double[N];
		double[] artmp = new double[N];

		this.population.sort();

		// calculate xmean and BDz
		for (int i = 0; i < N; i++) {
			xmean[i] = 0;

			for (int j = 0; j < mu; j++) {
				xmean[i] += weights[j] * population.getIndividual(j).getGene(i);
			}

			BDz[i] = Math.sqrt(mueff) * (xmean[i] - xold[i]) / sigma;
		}

		// cumulation for sigma (ps) using B*z
		if (diagonalIterations >= iteration) {
			// given B=I we have B*z = z = D^-1 BDz
			for (int i = 0; i < N; i++) {
				ps[i] = (1.0 - cs) * ps[i] + Math.sqrt(cs * (2.0 - cs)) * BDz[i] / diagD[i];
			}
		} else {
			for (int i = 0; i < N; i++) {
				double sum = 0.0;

				for (int j = 0; j < N; j++) {
					sum += B[j][i] * BDz[j];
				}

				artmp[i] = sum / diagD[i];
			}

			for (int i = 0; i < N; i++) {
				double sum = 0.0;

				for (int j = 0; j < N; j++) {
					sum += B[i][j] * artmp[j];
				}

				ps[i] = (1.0 - cs) * ps[i] + Math.sqrt(cs * (2.0 - cs)) * sum;
			}
		}

		// calculate norm(ps)^2
		double psxps = 0;
		for (int i = 0; i < N; i++) {
			psxps += ps[i] * ps[i];
		}

		// cumulation for covariance matrix (pc) using B*D*z
		int hsig = 0;
		if (Math.sqrt(psxps) / Math.sqrt(1.0 - Math.pow(1.0 - cs, 2.0 * iteration)) / chiN < 1.4 + 2.0 / (N + 1)) {
			hsig = 1;
		}

		for (int i = 0; i < N; i++) {
			pc[i] = (1.0 - cc) * pc[i] + hsig * Math.sqrt(cc * (2.0 - cc)) * BDz[i];
		}

		// update of C
		for (int i = 0; i < N; i++) {
			for (int j = (diagonalIterations >= iteration ? i : 0); j <= i; j++) {
				C[i][j] = (1.0 - (diagonalIterations >= iteration ? ccovsep : ccov)) * C[i][j]
						+ ccov * (1.0 / mueff) * (pc[i] * pc[j] + (1 - hsig) * cc * (2.0 - cc) * C[i][j]);

				for (int k = 0; k < mu; k++) {
					C[i][j] += ccov * (1 - 1.0 / mueff) * weights[k]
							* (population.getIndividual(k).getGene(i) - xold[i])
							* (population.getIndividual(k).getGene(j) - xold[j]) / sigma / sigma;
				}
			}
		}

		// update of sigma
		sigma *= Math.exp(((Math.sqrt(psxps) / chiN) - 1) * cs / damps);
	}

	public Individual run() {
		Params.countEvals = 0;

		best = samplePopulation();
		System.out.println(iteration + "\t" + best.fitness);

		while (Params.countEvals < Params.MAX_EVALS) {
			updateDistribution();
			Individual indiv = samplePopulation();

			if (indiv.fitness < best.fitness) {
				best = indiv;
			}

			System.out.println(iteration + "\t" + best.fitness);
		}

		return best;
	}
}
