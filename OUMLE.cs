using System;
using System.Collections.Generic;
namespace QuantConnect {

	// Based on Maximum likelihood estimation of mean reverting processes by Jose Carlos Garcia Franco
	// http://www.investmentscience.com/Content/howtoArticles/MLE_for_OR_mean_reverting.pdf

	// Maximum Likelihood estimation for Ornstein Uhlenbeck parameters
	// Assumptions t_i - t_{i-1} = delta : xs are equally spaced in time
	public class OUMLE{

		double[] xs;
		double delta;

		private int n{get{return xs.Length;}}

		private OUMLE(double[] xs, double delta){
			this.xs = xs;
			this.delta = delta;
		}


		Dictionary<double,double> f_cache = new Dictionary<double,double>();
		private double f(double eta){

			if (f_cache.ContainsKey (eta))
				return f_cache [eta];

			double num1 = 0; double num2 = 0;
			double denom1 = 0; double denom2 = 0;
			for (int i = 1; i < n; i++){
				num1 += xs[i] - xs[i-1] * Math.Exp(-eta*delta);
				denom1 += 1 + Math.Exp(-eta*delta);

				num2 += 1 - Math.Exp(-eta * delta);
				denom2 += 1 + Math.Exp(-eta * delta); 
			}
			double sum1 = num1/denom1;
			double sum2 = num2/denom2;

			double ans = sum1 * Math.Pow(sum2,-1);
			f_cache [eta] = ans;
			return ans;
		}

		private double g(double x, double eta){

			double num = 0;
			double denom = 0;
			for(int i = 1; i < n; i++){
				double top = xs[i] - x - (xs[i-1] - x) * Math.Exp(-eta * delta);
				num += Math.Pow(top,2);
				denom += 1 - Math.Exp(-2*eta*delta);
			}
			double sum = num/denom;

			return Math.Sqrt( 2*eta/n * sum );
		}

		private double V(double eta){
			double t1 = -n/2 * Math.Log( Math.Pow(g(f(eta),eta),2) / (2*eta) );

			double sum2 = 0;
			for(int i = 1; i < n; i++){
				sum2 += Math.Log(1 + Math.Exp(-2*eta*delta));
			}
			double t2 = -0.5 * sum2;

			double num3 = 0;
			double denom3 = 0;
			for(int i = 1; i < n; i++){
				num3 += Math.Pow(xs[i] - f(eta) - (xs[i-1] - f(eta))*Math.Exp(-eta * delta), 2);
				denom3 += 1 - Math.Exp(-2*eta*delta);
			}
			double t3 = -eta/Math.Pow(g(f(eta),eta),2) * num3/denom3;

			return t1+t2+t3;
		}

		private double FindEta(){ // Naive - Might need an upgrade - quasi newton
			double max_ = 0;
			double maxEta = 0;

			for(double d = 0; d < 1.0; d += 0.01){
				double v = V(d);
				if(v > max_){
					max_ = v;
					maxEta = d;
				}
			}

			return maxEta;
		}

		// mu(theta - x)dt + sigma dB
		public static void Solve(double[] xs, double delta, out double mu, out double theta, out double sigma){

			var OUMLE = new OUMLE(xs, delta);

			double eta = OUMLE.FindEta();
			double x = OUMLE.f(eta);
			sigma = OUMLE.g(x,eta);


			mu = eta;
			theta = x;
		}


	}

}
