package npairs.utils;

import com.jmatio.types.MLDouble;
import com.jmatio.types.MLStructure;

import extern.NewMatFileReader;

import npairs.io.NiftiIO;
import npairs.shared.matlib.*;

/**<head>
 * <style>
 * <!--
 * .tab { margin-left: 40px; }
 * -->
* </style>
 * 
 * Implements the Procrustes rotation. Computes the optimal rotation matrix that maps 
 * a given set of points B onto a reference set of points A using the algorithm of Procrustes.
 * 
 *  <h4>Why is Procrustes useful?</h4>
 * <p>The data in B may need to be transformed before it can be compared directly 
 * 		with A.  For example, canonical variates from a canonical variate analysis (CVA) may live in
 * 		a different space for each subsample of given dataset, as variance measures
 * 		and other calculations involved in producing the CVA results will vary across these 
 * 		subsamples even though they are all drawn from the same underlying distribution.
 * 		CVA dimensions from a given set of results might be flipped in sign or permuted 
 * 		compared to another set of results.  They probably also live in a slightly different
 * 		space than the points of A, meaning that B might need to be rotated in order to live in A-space.
 * 
 * <h4>The Algorithm</h4>
 * <ul>
 * <li>Let A be a matrix containing the 'reference set' of data points.
 * <li>Let B be a matrix to be compared with A.
 * </ul>
 * <p>A and B must have the same dimensions (P rows X N columns, say). Data 
 * 	in each matrix must be arranged such that rows are variables ("dimensions") and columns are observations 
 * 	(N "points" in P-dimensional space).
 * 
 *<p>The Procrustes algorithm computes the linear transformation Q that minimizes the difference between 
 *A and B under the <i>Frobenius norm</i> || ||<sub>F</sub>:
 *
 * 			<p class = "tab">|| A - QB ||<sub>F</sub> :=	sqrt (&#8721;&#8721;|(A - QB)<sub><i>ij</i></sub>|<sup>2</sup>)
 * 					(i.e., the square root of the sum of the absolute squares of the elements of A - QB).</p>
 * 			
 * 			<p>It can be shown that || A - QB ||<sub>F</sub> is minimized when Q = UV<sup>t</sup>, where 
 *           U and V are the orthonormal eigenvector matrices from a singular value decomposition (SVD) of 
 *           AB<sup>t</sup>:
 * 			<p class = "tab">AB<sup>t</sup> = USV<sup>t</sup> (where S is a diagonal matrix of 
 *           singular values).</p>
 * 			
 * <p>To map B approximately into A-space, we can take a couple of approaches:  
 *<ol>
 *<li><b><i>Permutation and Reflection</i></b> of dimensions (rows) of B to match the order and sign of rows in A.
 *(<i>Reflection</i> means <i>negation</i>, i.e., multiplication of every element by -1). This is a simpler
 * approach than applying the full Procrustes transformation Q to B and ignores rotation. It is what is used 
 * in IDL NPAIRS.
 * 
 * <li><b><i>Rotation</i></b> of B into A-space (including also permutation and reflection at the same time). 
 * 
 * <p>To rotate B into A, we apply the transformation matrix Q: "new B" = QB. This is the full Procrustes 
 * transformation of B that minimizes the differences between A and B under the Frobenius norm. 
 * </ol>
 * 
 * <p>For more details, see Golub, G. H. and Van Loan, C. F. <i>Matrix Computations, 3rd ed. </i> Baltimore, MD: 
 * Johns Hopkins, 1996.
 *
 * <p>REQUIRED: Set npairs.shared.matlib.MatrixImpl.matlibType outside of this method.
 *	
 *</head>
 */
public class Procrustes {
	
	/**
	 * debug flag
	 */
	private static boolean debug = false;
	/**
	 * coordinates of reference data; this matrix has P rows and
     * N columns, i.e. fltarr(N,P) in idl notation, where P is the dimension of the
     * data set and N is the number of observations (samples).
	 */
	private Matrix A;
	/**
	 * coordinates of points to be compared with A; this data is transformed
	 * by Q to match A more closely.
	 */
	private Matrix B;
	/**
	 * optimal rotation matrix (such that || A - QB || --> min)
	 */
	private Matrix Q;
	/**
	 * vector specifying the distances between each 
     * point in A and its corresponding rotated point in QB
	 */
	private double[] res;
	/**
	 * Frobenius norm between A and QB: fnorm = || A - QB ||
	 */
	private double fnorm;
	/**
	 * refSign - an array indicating, for each input dim (row) of B, 
	 * whether it needs to be negated in order to map into A: -1 => negate;
	 * +1 => don't negate.
	 */
	private int[]refSign;
	/**
	 * refPermute - an array indicating, for each input dim (row) of B, 
	 * which dim (row index) it should occupy in order to map into A
	 * (in case dims are swapped in B vs A.)
	 */
	private int[]refPermute;
	/** 
	 * normedA, normedB are A and B, resp., normalized by dividing by max 
	 * magnitude value across BOTH A and B.  Hence max (normedA, normedB) == 1.
	 */
	private Matrix normedA;
	private Matrix normedB;
	/** 
	 * QB is the result after B has been rotated into A-space by Procrustes 
	 * transformation Q. 
	 */
	private Matrix QB;
	
	/** No. of columns (observations) in input data: */
	private int nPts;
	
	/** No. of rows (variables) in input data: */
	private int nDim;
	
	/** ABtMat = A * transpose(B); this is the Matrix to which SVD is applied: */
	private Matrix ABtMat;
	
	/** corrMat = corr (A,B) across dimensions (rows) (used to calculate permutations):  */
	private Matrix corrMat;
	/**
	 * SVD used in Procrustes rotation:
	 */
	private SingularValueDecomposition svd;
	
	/** Computes Procrustes transformation Q and associated
	 * reflection, permutation and residual error info.
	 * @param A matrix holding reference pattern 
	 * @param B matrix holding pattern to be mapped to A
	 */
	public Procrustes(Matrix A, Matrix B) {
                
		this.A = A;
        this.B = B;
        if ( A == null || B == null ) 
        	throw new IllegalArgumentException("Procrustes: Error - must specify " +
        			"input Matrices");
        
        nDim=A.numRows();
        nPts=A.numCols();
        if( nDim != B.numRows() || nPts != B.numCols() )
            throw new IllegalArgumentException("Procrustes: Error - input Matrices " +
            		"must be the same size");
        
        this.compute();
        
    } 
	
	/**
	 * Calculates maximum absolute value across both A and B.
	 * @return maximum value (absolute, i.e. magnitude) of both A and B points.
	 */
	private double max() {
		
		double mx = Math.abs(A.getQuick(0,0));
		
		// For each pt in each Matrix A and B, if current pt > current mx, 
		// set mx = current pt.
		for (int i=0; i < nDim; i++) {
			for (int j=0; j < nPts; j++) {
				if ( mx < Math.abs(A.getQuick(i,j)) ) 
					mx = Math.abs(A.getQuick(i,j));
				if ( mx < Math.abs(B.getQuick(i,j)) ) 
					mx = Math.abs(B.getQuick(i,j));
			}
		}	
		return mx;
		} //end method max()
		
	/**
	 * maxAbsLoc() find the location of maximum magnitude (i.e., absolute value) of 
	 * values in input array
	 * @param array input array
	 * @return index of maximum abs val in the input array 
	 */
	private int maxAbsLoc (double[] array) {
		
		int loc = 0;
		double temp = Math.abs(array[0]);
		for (int i=1; i < array.length; i++) {
			if ( Math.abs(array[i]) > temp) {
				loc = i;
				temp = Math.abs(array[i]);
			}
		}
		return loc;
	}
	
	/**
	 * Main computation method to compute Procrustes transformation Q and associated
	 * reflection, permutation and residual error info.
	 */
	private void compute() {
		
		refSign        = new int[nDim];
		refPermute     = new int[nDim];
		
		double mx = max();
		// Normalize data so absval (max(A,B)) == 1:
		normedA = A.mult(1/mx);
		normedB = B.mult(1/mx);
		
		if(nDim > 1) { // CASE 1 - Input data has more than one dim:
		
			/* Find Procrustes transformation matrix 
			 * 		Q = U * transpose(V), 
			 * where U and V come from SVD of A * transpose(B):
			 * 		A * transpose(B) =  U * S * transpose(V).
			 * (More precisely, we take SVD of ABtMat = (1/mx)^2 * A * transpose(B),
			 * 	which produces same U and V as SVD of A * transpose(B).)
			 */
			ABtMat = normedA.mult(normedB.transpose());
			svd    = ABtMat.svd();
			
//			if (debug) {
//				System.out.println("Finished running SVD(ABtMat)...");
//				System.out.println("U: ");
//				svd.getU().print();
//				System.out.println("V: ");
//				svd.getV().print();
//				System.out.println("S: ");
//				svd.getS().print();
//				System.out.println("S vals: ");
//				npairs.io.NpairsjIO.print(svd.getSVals());
//			}
			
			Q      = svd.getU().mult(svd.getV().transpose());
		
			/* Calculate permutation info:
			 */
			corrMat = corr(A, B);
			/* For each dimension (row) of input datasets, determine whether 
			 * 	maximum absolute value in corresponding row of transformation Matrix
			 * Q is positive or negative.  If negative, set corresponding element in 
			 * refSign array to -1; if positive, set corresp. refSign element to 1:
			 */
			int mxLoc;
			
			Matrix tempQ = Q.copy();
			if (debug) {
				System.out.println("Rotation mat: ");
				tempQ.print();
			}
			
			for (int i=0; i < nDim; i++)	{
				
				// Find index of max abs val in current row of 'tempQ':
				mxLoc = maxAbsLoc(tempQ.getRowQuick(i));
				
				/* Check sign of max abs val of current row and set corresp. element 
				 * of refSign to +1 (if positive tempQ[mxLoc]) or -1 (if neg. tempQ[mxLoc]):
				 */
				if (tempQ.getQuick(i,mxLoc) > 0) {
				   refSign[i] =  1;
				} else refSign[i] = -1;
				
				// Set refPermute:
				refPermute[i] = mxLoc;
				
				/* Following Jon Anderson's strategy in IDL NPAIRS code cva_reference.pro, 
				 * set mxLoc column of tempQ to zero (can't have more than 1 dim permuted into 
				 * the same row, anyway):
				 */
				for (int j = 0; j < nDim; j++) {
					tempQ.setQuick(j, mxLoc, 0);
				}
			}
			
		} // end nDim > 1
		
		else  { // CASE 2 - Input data has only 1 dim (row). Then Q = +/- [1](a 1X1 matrix).   
			
			int Q1D; // scalar value inside 1X1 matrix Q (== +/-1)
			
			Matrix diffMat, sumMat;
			double diffSumSqrs = 0;
			double sumSumSqrs  = 0;
			diffMat = normedA.minus(normedB);
			sumMat  = normedA.plus(normedB);
			
			// Calculate sums of squares of elements in diffMat and sumMat (for each Mat separately):
			for(int i=0; i < nDim; i++)
				for(int j=0; j < nPts; j++) {
					diffSumSqrs += diffMat.getQuick(i,j) * diffMat.getQuick(i,j);
					sumSumSqrs  +=  sumMat.getQuick(i,j) * sumMat.getQuick(i,j);
				}
			
			/*  If diffSumSqrs < sumSumSqrs, then corresp. points in A and B tend to have the same sign; 
			 * 	if diffSumSqrs > sumSumSqrs, then corresp. points in A and B tend to have opposite signs.
			 * 	Calculate Q:
			 */
			if (diffSumSqrs <= sumSumSqrs) Q1D =  1; 
			else 				   		   Q1D = -1;
			MatrixImpl mImpl = new MatrixImpl(new double[][] {{Q1D}});
			Q = mImpl.getMatrix();
			
			// For nDim == 1, refSign
			refSign = new int[1]; 
			refSign[0] = Q1D;
			refPermute = new int[1];
			refPermute[0] = 0;
		} // end if (nDim == 1)

		// Calculate QB:
		QB = Q.mult(B);
		// Compute residual error and Frobenius norm:
		computeRes(A, QB);
		
	} // end method compute()

	/**
	 * Computes the residual errors (Euclidean distance for each observation) and Frobenius norm of 
	 * A - QB.
	 */
	private void computeRes(Matrix A, Matrix QB) {
		
		res = new double[nPts];
		for(int j=0; j < nPts; j++) {
			for(int i=0; i < nDim; i++) {
				
				res[j] += (A.getQuick(i,j) - QB.getQuick(i,j)) * (A.getQuick(i,j) - QB.getQuick(i,j));
				fnorm  += (A.getQuick(i,j) - QB.getQuick(i,j)) * (A.getQuick(i,j) - QB.getQuick(i,j));
			}
			res[j] = Math.sqrt(res[j]);
		}
		fnorm = Math.sqrt(fnorm);
	}
	
	/** Computes correlation Matrix of ROWS OF input Matrices:
	 */
	private Matrix corr(Matrix A, Matrix B) {
		
		double   currAMean = 0;
		double   currBMean = 0;
		double[] stdDevA = new double[nDim];
		double[] stdDevB = new double[nDim];
		MatrixImpl mImpl = new MatrixImpl(nDim, nPts);
		Matrix aCentred = mImpl.getMatrix();
		Matrix bCentred = aCentred.copy();
		
		// correlate across rows:
		for (int j=0; j < nDim; j++) {
			for (int i=0; i < nPts; i++) {
				currAMean += A.getQuick(j,i)/nPts;
				currBMean += B.getQuick(j,i)/nPts;
			}
			for (int i=0; i < nPts; i++) {
				aCentred.setQuick(j, i, A.getQuick(j,i) - currAMean);
				bCentred.setQuick(j, i, B.getQuick(j,i) - currBMean);
				stdDevA[j] += aCentred.getQuick(j,i) * aCentred.getQuick(j,i);
				stdDevB[j] += bCentred.getQuick(j,i) * bCentred.getQuick(j,i);
			}
			stdDevA[j] = Math.sqrt(stdDevA[j]);
			stdDevB[j] = Math.sqrt(stdDevB[j]);
			for (int i=0; i < nPts; i++) {
				aCentred.setQuick(j, i, (aCentred.getQuick(j,i) / stdDevA[j]));
				bCentred.setQuick(j, i, (bCentred.getQuick(j,i) / stdDevB[j]));
			}
			currAMean = 0;
			currBMean = 0;
		}
		
		corrMat = aCentred.mult(bCentred.transpose());		
		return corrMat;
	}

	/** Returns the correlation matrix of rows (i.e., variables) of input data Matrices A, B.
	 * @return correlation matrix
	 */
	public Matrix getCorrMat() {
		return corrMat;
	}
	
	
	/**
	 * Returns the result of applying Procrustes transformation Q to B.
	 * @return rotated pattern QB 
	 */
	public Matrix getRotatedPatt()
	{
		return QB;
	}
	
	/**
	 * Returns the transformation (rotation) matrix Q that minimizes Frobenius norm 
	 * || A - QB ||<sub>F</sub>.
	 * @return optimized transformation matrix Q
	 */
	public Matrix getRot() {
		return Q;
	}
	
	/**
	 * Returns array of signed integers +/-1 indicating whether to negate each dimension in B.
	 * For each dimension (row) of input data sets, element with maximum absolute value in 
	 * corresponding row of transformation Matrix Q is either positive or negative.  If negative,
	 * corresponding element in reflection sign array = -1; if positive, corresponding reflection  
	 * sign array element = 1. 
	 * <p> Note that if dimensions are permuted as well as reflected, reflection is applied to 
	 * after permutation. 
	 * <p> E.g., if permuted B = {B<sub>1</sub>, 
	 *  B<sub>2</sub>, B<sub>0</sub>}' as in example from {@code #getRefPermute()} documentation
	 *  and reflection array = {1, 1, -1}, then B after permutation and reflection = 
	 *  {B<sub>1</sub>, B<sub>2</sub>, -B<sub>0</sub>}'. 
	 * 
	 * @return Array, length # rows in A and B, of signed integers indicating whether 
	 * corresponding rows of B should be negated (-1) or left alone (+1) when reflected to 
	 * match orientation of rows in A.
	 * 
	 * @see #getRefPermute()
	 */
	public int[] getRefSign() {
		return refSign;
	}
	
	/**
	 *  Returns array of permutation indices for dimensions (rows) in B. For each dimension, corresponding
	 *  permutation index indicates which dimension is to move to current location.  
	 *  E.g., if there are 3 dimensions in A and B and permute = {1, 2, 0}, permuted B = {B<sub>1</sub>, 
	 *  B<sub>2</sub>, B<sub>0</sub>}', where B<sub>i</sub> is ith row of original B (0-relative).
	 *  Permutation index for each dimension is set to index of element with maximum absolute 
	 *  value in corresponding row of transformation Matrix Q.
	 * @return Array, length # rows in A and B, containing indices of permuted rows of B.
	 * 
	 * @see #getRefSign()
	 */
	public int[] getRefPermute()
	{
		return refPermute;
	}
	
	/** Returns Frobenius norm || A - QB ||<sub>F</sub>.
	 * @return Minimized Frobenius norm value from Procrustes computation.
	 */
	public double getFnorm() {
		return fnorm;
	}
	
	/** Returns residual error values between observations in A and transformed 
	 *  Matrix QB. 
	 *  
	 *  @return Array, length # observations (columns) in A and B, containing
	 *          Euclidean distance between corresponding observations in A and QB
	 *          (i.e., the <i>l<sup>2</sup></i>-norm of each column of A - QB).
     */
	public double[] getRes() {
		return res;
	}
	
	/**
	 * main() function: for testing. 
	 */
	private static void main(String args[]) {

//		double[][] testdata1={{1.2,2.1,4.1}, {2.1,3.3,2.5}, {4.1,2.5,2.3}, {2.2,3.1,3.5}};
//		double[][] testdata2={{-1.2,-4.1,2.1}, {-2.1,-2.5,3.3}, {-4.1,-2.3,2.5}, {-2.2,-3.5,3.05}};
//		String[] niftiFilenames = {"/home/anita/Desktop/y246_5cond_aug4_2010_8_12pc_zs_cv1_NPAIRSJresult.nii",
//				"/home/anita/Desktop/y246_5cond_aug4_2010_8_12pc_zs_cv2_NPAIRSJresult.nii",
//				"/home/anita/Desktop/y246_5cond_aug4_2010_8_12pc_zs_cv3_NPAIRSJresult.nii"};
		String[] niftiFilenames = {"/home/anita/Desktop/y246_5cond_aug4_2010_8_12pc_zs_cv1_NPAIRSJresult.nii"};
		String resFile = "/home/anita/plsnpairs/grady/results/y246_5cond_aug4_2010/" +
			"y246_5cond_aug4_2010_008pc_NPAIRSJresult.mat";

		try {
		MLDouble stCoords = (MLDouble)new NewMatFileReader(resFile).
			getContent().get("st_coords");
		int[] maskCoords = stCoords.getIntFirstRowOfArray();
		int nDims = niftiFilenames.length;
		double[][] testdata1 = new double[nDims][maskCoords.length];
		double[][] testdata2 = new double[nDims][maskCoords.length];
		for (int i = 0; i < nDims; ++i) {
			testdata1[i] = NiftiIO.readNiftiData(niftiFilenames[i], 0, maskCoords);
			testdata2[i] = NiftiIO.readNiftiData(niftiFilenames[i], 1, maskCoords);
		}
		
		System.out.println("Size A: " + testdata1.length + " X " + testdata1[0].length);
		System.out.println("Size B: " + testdata2.length + " X " + testdata2[0].length);
		
		String matlibType = args[0];
		System.out.println("Using matlib type: " + matlibType);

			MatrixImpl mImpl1 = new MatrixImpl(testdata1, matlibType);
			Matrix test1 = mImpl1.getMatrix();
			MatrixImpl mImpl2 = new MatrixImpl(testdata2);
			Matrix test2 = mImpl2.getMatrix();
//			System.out.println("Input A: ");
//			test1.transpose().print();
//			System.out.println("Input B: ");
//			test2.transpose().print();

			Procrustes pTest=new Procrustes (test1, test2);
			//pTest.compute();
			System.out.println("Procrustes compute completed.");
//			System.out.println("RefPermute: ");
//			npairs.io.NpairsjIO.print(pTest.getRefPermute());
//		    System.out.println("RefSign: ");
//			npairs.io.NpairsjIO.print(pTest.getRefSign());
//			System.out.println("res = ");
//			npairs.io.NpairsjIO.print(pTest.res);
//		    System.out.println("fnorm = "+pTest.fnorm);
			System.out.println("Size Q: " + pTest.Q.numRows() + " X " + pTest.Q.numCols());
//			pTest.Q.print();
			System.out.println("Size QB: " + pTest.QB.numRows() + " X " + pTest.QB.numCols());
//			System.out.println("Transformed B ( == QB):");
//			pTest.getRotatedPatt().print();
			System.out.println("Correlation (A, B): "); 
			test1.transpose().correlate(test2.
				transpose()).print();
			System.out.println("Correlation (A, QB): ");
			test1.transpose().correlate(
					pTest.QB.transpose()).print();
			System.out.println("Correlation (B, QB): "); 
			test2.transpose().correlate(
					pTest.QB.transpose()).print();
			
			
			
		}
		catch (Exception e) {
			e.printStackTrace();
		}


	}

}
