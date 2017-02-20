
import Jama.*;

public class BayesianCurveFitting {
	public static final double alpha = 0.005;
	public static final double beta = 11.1;

	public static final int x = 11;// rows # of input
	public static final int M = 15;// dimension #
	public static final int Mp = M + 1;
	public static int[] N = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	public static double datas[] = { 1222.52, 1222.51, 1222.29, 1222.80, 1222.37, 1222.272, 1222.356, 1222.53, 1222.69, 1222.90 };

	

	public static double[][] getPhiX(int a) {
		// ArrayList<Double> A = new ArrayList<Double>();
		double[][] PhiX = new double[Mp][Mp];
		int i = 0;
		while (i < Mp) {
			PhiX[i][i] = Math.pow(x, i);
			i++;
		}
		
		/*for (int j = 0; j < Mp; j++) {
			for (int k = 0; k < Mp; k++) {
				System.out.print(PhiX[j][k] + " ");
			}
			System.out.println();
		}*/

		return PhiX;
	}

	public static Matrix getInverseS() {
		double[][] B = new double[Mp][Mp];
		for (int i = 0; i < Mp; i++) {
			B[i][i] = 0.0;
		}
		Matrix I = Matrix.identity(Mp, Mp);// unit matrix I;
		Matrix temp = Matrix.constructWithCopy(B);
		Matrix sum = temp;//sigma(PhiX*PhiXT);
		for (int i = 0; i < N.length; i++) {
			Matrix PhiX = Matrix.constructWithCopy(getPhiX(N[i]));
			
			/*double[][] original = PhiX.getArray();
			double[][] col = new double[Mp][1];
			for(int j = 0; j < Mp; j++){
				col[j][0] = original[j][j];
			}*/
			
			//Matrix PhiX1Col = Matrix.constructWithCopy(col);
			Matrix PhiXT = PhiX.transpose();
			temp = PhiX.times(PhiXT);
			sum = sum.plus(temp);
		}
		
		Matrix BetaPhi = temp.times(beta);
		Matrix AlphaI = I.times(alpha);
		Matrix InverseS = AlphaI.plus(BetaPhi);
		
		return InverseS;
	}
	
	public static Matrix getMean(){
		Matrix PhiX = Matrix.constructWithCopy(getPhiX(x));
		
		/*double[][] original = PhiX.getArray();
		double[][] col = new double[Mp][1];
		for(int j = 0; j < Mp; j++){
			col[j][0] = original[j][j];
		}
		Matrix PhiX1Col = Matrix.constructWithCopy(col);*/
		
		Matrix PhiXT = PhiX.transpose();
		Matrix S = getInverseS().inverse();
		Matrix m1 = PhiXT.times(beta);
		Matrix m2 = m1.times(S);
		
		double[][] B = new double[Mp][Mp];
		for (int i = 0; i < Mp; i++) {
			B[i][i] = 0.0;
		}
		
		Matrix sum = Matrix.constructWithCopy(B);
		
		for(int i = 0; i < N.length; i++){
			Matrix PhiXn = Matrix.constructWithCopy(getPhiX(N[i]));
			Matrix temp = PhiXn.times(datas[i]);
			sum = sum.plus(temp);
		}
		
		return m2.times(sum);
	}
	
	public static Matrix getVariance(){
		Matrix PhiX = Matrix.constructWithCopy(getPhiX(x));
		
		/*double[][] original = PhiX.getArray();
		double[][] col = new double[Mp][1];
		for(int j = 0; j < Mp; j++){
			col[j][0] = original[j][j];
		}
		Matrix PhiX1Col = Matrix.constructWithCopy(col);*/
		
		Matrix PhiXT = PhiX.transpose();
		Matrix S = getInverseS().inverse();
		
		double p1 = 1/beta;
		Matrix m1 = PhiXT.times(S);
		Matrix m2 = m1.times(PhiX);
		Matrix unit = Matrix.identity(Mp,Mp);
		Matrix aux = unit.times(p1);
		
		return aux.plus(m2);
		
	}
	public static void main(String[] args) {
		Matrix mean = getMean();
		Matrix variance = getVariance();
		
		
		double[][] Mean = mean.getArray();
		double[] MeanX = new double[Mp];
		
		for(int i = 0; i < Mp; i++){
			for(int j = 0; j < Mp; j++){
				MeanX[i] = Mean[i][i]/10;
			}
		}
		
		double[][] Variance = variance.getArray();
		double[] VarianceX = new double[Mp];
		
		for(int i = 0; i < Mp; i++){
			for(int j = 0; j < Mp; j++){
				VarianceX[i] = Variance[i][i];
			}
		}
		
		for(int step = 1205; step < 1245;){
			double p1 = 1/Math.sqrt(2*Math.PI*VarianceX[0]);
			double top = (step - MeanX[0]);
			double p2 = Math.pow(Math.E, -(Math.pow(top, 2)/(2*VarianceX[0])));

			double result = p1*p2;
			step += 1;
			System.out.println(result + "    " + step);
		}
		
//		for(int i = 0; i < Mp; i++){
//			System.out.println(MeanX[i] + "     " + VarianceX[i]);
//		}
		
		
		
		/*for(int i = 0; i < Mp; i++){
			for(int j = 0; j < Mp; j++){
				System.out.print(qw[i][j] + " ");
			}
			System.out.println();
		}*/
	}
	
}