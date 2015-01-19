package npairs.utils;

import java.io.IOException;
import java.io.PrintStream;

import npairs.io.NpairsDataLoader;
import npairs.shared.matlib.*;
import extern.niftijlib.Nifti1Dataset;

import pls.shared.MLFuncs;
//import npairs.Npairsj;
import npairs.NpairsException;
import npairs.NpairsSetupParams;

/** Performs Principal Component Analysis (PCA).
 * 
 * <p>PCA produces a representation of input data as a linear combination of orthonormal
 * components. Each component accounts for as much variance as possible uncorrelated with 
 * the other components: the 1st component accounts for the maximum amount of variance
 * that can be expressed in a single dimension; the 2nd component accounts for the maximum amount
 * of variance uncorrelated with the 1st component, etc. This means that for any given number 
 * <i>k</i> <= <i>d</i> (where <i>d</i> is total # of PC dimensions), the most information 
 * is preserved if the data is expressed via linear combinations of the first <i>k</i> PCs; i.e., 
 * the PCs serve as basis vectors for a minimally lossy representation of the data. 
 * 
 * <p>In NPAIRS, PCA is used as a "denoising" and dimension reduction step before data is 
 * passed into a CVA analysis. 
 * 
 * <h4>The Algorithm</h4>
 * <p>This class uses Eigenvalue Decomposition (EVD) to do PCA.
 * <ul>
 * <li>Let <b>M</b> be a matrix of data with <i>m</i> rows (observations) and <i>n</i> columns (variables).
 * <li> Let <i>d</i> = minimum of <i>m</i>, <i>n</i>.
 * </ul>
 * Then PCA of <b>M</b> is accomplished via EVD of the covariance matrix of <b>M</b>:
 * <ul>
 * <li>cov(<b>M</b>) = <b>V&#923;<sup>2</sup>V</b><sup>t</sup>, where 
 * <p>
 * <ul>
 * 	 <li>cov(<b>M</b>) = (<b>M</b> - <b>&#956;</b>)<sup>t</sup>(<b>M</b> - <b>&#956;</b>)/(<i>m</i> - 1)
 * 			(where <b>&#956;</b> is matrix containing sample (i.e., observation) mean vector in each row),
 * 	 <p>
 *	 <li><b>&#923;<sup>2</sup></b> is a diagonal <i>d</i> X <i>d</i> matrix of PC eigenvalues in descending order, and
 *   <p>
     <li><b>V</b> is an  <i>n</i> X <i>d</i> matrix of PC eigenvectors (in columns).
 * </ul>
 * </ul>
 * 
 * <p>Eigenvalue <i>&#955;<sub>i</sub></i> = amount of variance accounted for by <i>i<sup>th</sup></i> PC
 * component (eigenvector). If the <i>n</i> input variables are not linearly independent (always the case if 
 * <i>m</i> < <i>n</i>), then only the first <i>p</i> eigenvalues in <b>&#923;<sup>2</sup></b> 
 * will be non-zero, where <i>p</i> is the actual number of linearly independent input variables (i.e.,
 * <i>p</i> = rank of input data matrix).
 * 
 * <p>PC scores <b>P</b> are coefficients of mean-centred input data expressed as a linear combination
 * 		of <i>q</i> <= <i>d</i>  PC eigenvectors:
 * 			<ul>
 * 				<li>(<b>M</b> - <b>&#956;</b>)<sub><i>m X n</i></sub> &#8776; <b>P</b><sub><i>m X q</i></sub>  
 * 					<b>V</b><sup>t</sup><sub><i>q X n</i></sub> 
 * 				<p>If <i>q</i> = <i>p</i> then the equality is exact. 
 * 				</ul>
 * 
 * 			
 * In other words, <b>V</b><sup>t</sup> can be viewed as a linear transformation 
 * of mean-centred data into PC space.
 * <ul>
 * <li><b>P</b><sup>t</sup><sub><i>q X m</i></sub> = <b>V</b><sup>t</sup><sub><i>q X n</i></sub>
 * 		(<b>M</b> - <b>&#956;</b>)<sup>t</sup><sub><i>n X m</i></sub> 
 * </ul>
 * 
 * @author Anita Oder
 * 
 */
public class PCA {
	
//		private boolean normalizeBySD = false; // if true, PC Scores are normalized by
//											   // their standard deviations
//		                                       // to have variance 1.
		private boolean debug = false;
		private Matrix eigenvectors;  // contains PCA eigenvectors V in columns
		private double[] eigenvalues; // array containing PCA eigenvalues in descending
		                              // order (contents of main diagonal in S)
		private Matrix evalMat;       // diagonal Matrix containing eigenvalues down main diagonal
//		private Matrix covMat; // for debugging
		private Matrix pcScores; // projection of original data onto PC eigenvectors
		                         // (pcScores = MV)		
		private boolean pcaInOrigSpace; // Usually do PCA on feature-selected data, in which
		                                // case PCA will
		                                // not be in orig space unless explicitly projected back
		                                // into it.  
		private PrintStream logStream = null; // where to log info; null if no logging to be done
		private boolean getCpuTime = true; // if true, log cpu time
		
		
	/** Constructor for PCA object. Calculates PCA of input Matrix via eigenvalue decomposition.
	 * 
	 * @param M Matrix of data. 
	 * 			<p>Assumption: rows are observations, columns are variables.
	 * @param normalizePCsBySD
	 * 			If true, normalize PC scores to have variance 
	 * 			1 for each PC dimension by dividing them by their standard deviation over the 
	 * 			dimension.
	 * @param inputDataInOrigSpace
	 * 			If true, input data lives in original 'image' (voxel) space instead of having
	 * 			been transformed into a feature selection (i.e., reduced dimension) space.
	 * @param logStream
	 * 			Where to log info. Set to null if no info is to be logged.
	 */
	public PCA(Matrix M, boolean normalizePCsBySD, boolean inputDataInOrigSpace, 
			PrintStream logStream) {
		if (logStream != null) {
			this.logStream = logStream;
		}
//		double startTime = System.currentTimeMillis();
							
		pcaInOrigSpace = inputDataInOrigSpace;
		/** CPU time for PCA! **/
		timing.Time t = new timing.Time();
		long startCpuTime = t.getCpuTime();
		long startUserTime = t.getUserTime();
		long startSystemTime = t.getSystemTime();
		computePCA(M, normalizePCsBySD);
		double cpuTime = (t.getCpuTime() - startCpuTime) / Math.pow(10,  9);
		double userTime = (t.getUserTime() - startUserTime) / Math.pow(10,  9);
		double systemTime = (t.getSystemTime() - startSystemTime) / Math.pow(10,  9);
		
		logStream.println("PCA time in seconds (CPU): " + cpuTime );
		logStream.println("PCA time in seconds (User): " + userTime );
		logStream.println("PCA time in seconds (System): " + systemTime );
		
//		double totTime = (System.currentTimeMillis() - startTime) / 1000;
//		if (debug) {
//			logStream.println("Total time PCA: " + totTime + " s");
//		}
	}
	
	/** Constructor for PCA object. Initializes PCA results. 
	 * 
	 */
	private PCA(Matrix eigenvectors, double[] eigenvalues, Matrix pcScores, boolean pcaInOrigSpace) {
		this.eigenvectors = eigenvectors;
		this.eigenvalues = eigenvalues;
		this.evalMat = new MatrixImpl(eigenvalues.length, eigenvalues.length).getMatrix();
		try {
			evalMat.setDiag(eigenvalues);
		}
		catch (MatrixException me) { 
			// matrix is square by construction
		}
		this.pcScores = pcScores;	
		this.pcaInOrigSpace = pcaInOrigSpace;
	}
	
	
	/** Calculates PCA of input Matrix M via eigenvalue decomposition
	 * 
	 * @param M Matrix of doubles. 
	 *            Assumption: rows are observations, cols are variables
	 * @param normalizeBySD TODO
	 */
	private void computePCA(Matrix M, boolean normalizeBySD) {
		int nRows = M.numRows();
		int nCols = M.numCols();
		
		double sTime = System.currentTimeMillis();
		Matrix meanCentredM = M.meanCentreColumns();	// M is NOT modified
		double tTime = (System.currentTimeMillis() - sTime) / 1000;
		if (debug) {
			logStream.println("\tTime mean-centring PCA input data: " + tTime + " s");
		}
		
		if (nRows < nCols) {
			sTime = System.currentTimeMillis();
			Matrix sspMMt = meanCentredM.sspByRow();
			tTime = (System.currentTimeMillis() - sTime) / 1000;
			if (debug) {
				logStream.println("\tTime SSP of mean-centred input PCA data: " + tTime + " s");
			}
			
			sTime = System.currentTimeMillis();
			EigenvalueDecomposition evd = sspMMt.eigenvalueDecomposition();
			tTime = (System.currentTimeMillis() - sTime) / 1000;
			if (debug) {
				logStream.println("\tTime EVD of SSP in PCA: " + tTime + " s");
			}
			
			Matrix leftEigenvects = evd.getEvects(); // = U given MMt = US^2Ut
			Matrix invSqrtEvalMat = evd.getInvSqrtRealEvalMat();
			sTime = System.currentTimeMillis();
			eigenvectors = (meanCentredM.transpose()).mult(leftEigenvects.mult(invSqrtEvalMat));
			tTime = (System.currentTimeMillis() - sTime) / 1000;
			if (debug) {	
				logStream.println("\tTime calculating PCA eigenvectors: " + tTime + "s");
			}
			// eigenvalues are scaled because evd was done on ssp matrix, 
			// not (sample) covariance matrix
			eigenvalues = scale(evd.getRealEvals(), 1./(double)(nRows - 1));
		}
		
		else {
			sTime = System.currentTimeMillis();
			Matrix sspMtM = meanCentredM.sspByCol();
			tTime = (System.currentTimeMillis() - sTime) / 1000;
			if (debug) {
				logStream.println("\tTime SSP of mean-centred input PCA data: " + tTime + " s");
//				String saveName = "RRSD_tests/javaTest/pcaSSPMtMMat.debug." + count;
//				sspMtM.printToFile(saveName, "IDL");
			}
			sTime = System.currentTimeMillis();
			EigenvalueDecomposition evd = sspMtM.eigenvalueDecomposition();
			tTime = (System.currentTimeMillis() - sTime) / 1000;
			if (debug) {
				logStream.println("\tTime EVD of SSP in PCA: " + tTime + " s");
			}
			eigenvectors = evd.getEvects(); // == V given MtM = VS^2Vt
            // eigenvalues are scaled because evd was done on ssp matrix, 
			// not (sample) covariance matrix
			eigenvalues = scale(evd.getRealEvals(), 1./(double)(nRows - 1));
		}		
		
		// TODO: optimize calc. of pc scores.
		sTime = System.currentTimeMillis();
		pcScores = meanCentredM.mult(eigenvectors);
		tTime = (System.currentTimeMillis() - sTime) / 1000;
		if (debug) {
			logStream.println("\tTime calculating PC scores: " + tTime + " s");
		}
		if (normalizeBySD) {
			sTime = System.currentTimeMillis();
			// normalize each pcScore column (dim) to have variance 1
			// TODO: should we just divide by sqrt(eval) instead of calculating
			// stddev explicitly?  
			for (int d = 0; d < pcScores.numCols(); ++d) {
				double[] currPCs = pcScores.getColumnQuick(d);
				double stdDev = MLFuncs.std(currPCs);
				pcScores.setColumnQuick(d, MLFuncs.divide(currPCs, stdDev));
			}
			if (debug) {
				tTime = (System.currentTimeMillis() - sTime) / 1000;
				logStream.println("\tTotal time normalizing PC scores: " + tTime + " s");
			}
		}
		evalMat = new MatrixImpl(eigenvalues.length, eigenvalues.length).getMatrix();
		try {
			evalMat.setDiag(eigenvalues);
		}
		catch (MatrixException me) { 
			// matrix is square by construction
		}
	}
	
	/** Returns Matrix <b>V</b> of PC eigenvectors.
	 * 
	 * @return PC eigenvectors
	 * 		<ul> <li> Matrix dimensions: <i>n</i> X <i>d</i>, where 
	 * 			<p><i>n</i> = # columns (variables) in input data 
	 *          <p><i>d</i> = # PC dimensions (i.e., minimum of <i>m</i> ( = # rows in 
	 *          input data) and <i>n</i>).
	 *          <li><i>i<sup>th</sup></i> column contains eigenvector corresponding to 
	 *          <i>i<sup>th</sup></i> PC dimension.
	 * 		</ul>
	 * 
	 * @see getEvects(int[])
	 * @see getEvals()
	 * @see rotateEigimsToOrigSpace(int[], Matrix, Matrix)
	 */
	public Matrix getEvects() {
		return eigenvectors;
	}
	
	/** Returns Matrix of PC eigenvectors corresponding to 
	 * 	given PC dimensions only.
	 * 
	 * @param pcDims array containing PC dimensions (0-relative) to include in 
	 * 			returned eigenvector Matrix
	 * @return PC eigenvectors for given PC dimensions
	 * 		<ul> <li> Matrix dimensions: <i>n</i> X <i>k</i>, where 
	 * 				<p><i>n</i> = # of columns (variables) in input data
	 * 				<p><i>k</i> = # of elements in given pcDims array
	 * 
	 * 		</ul>
	 * @see getEvects()
	 * @see getEvals()
	 * @see rotateEigimsToOrigSpace(int[], Matrix, Matrix)
	 */
	public Matrix getEvects(int[] pcDims) {
		return eigenvectors.subMatrixCols(pcDims);
	}
	
	/** Returns PC eigenvalues.
	 * 
	 * @return Array containing eigenvalues in descending order.
	 * 		<ul>
	 * 		<li>Array length = # PC dimensions (eigenvectors).
	 *      <li><i>i<sup>th</sup></i> eigenvalue <i>&#955;<sub>i</sub></i>
	 *      	= variance of <i>i<sup>th</sup></i> PC eigenvector.
	 *      </ul>
	 * 
	 * @see getEvects()
	 * @see getEvalMat()
	 * @see rotateEigimsToOrigSpace(int[], Matrix, Matrix)
	 *
	 */
	public double[] getEvals() {
		return eigenvalues;
	}
	
	/** Returns diagonal Matrix <b>&#923;<sup>2</sup></b> containing PC eigenvalues.
	 * 
	 * @return Square Matrix containing PC eigenvalues in descending order down the main diagonal.
	 *         <ul><li>Matrix size: <i>d</i> X <i>d</i>, where 
	 *         <p><i>d</i> = # PC dimensions (eigenvectors).
	 *         </ul>
	 * @see getEvects()
	 * @see getEvals()
	 * @see rotateEigimsToOrigSpace(int[], Matrix, Matrix)
	 */
	public Matrix getEvalMat() {
		return evalMat;		
	}
	
	/** Returns representation of original data transformed into PC space, i.e.: 
	 * <ul><li>Matrix of coefficients corresponding to weights of PC components
	 * when input data is expressed as linear combination of PC components.
	 * </ul>  
	 * @return PC scores
	 * <ul> 
	 * 		<li>Matrix size: <i>m</i> X <i>d</i>, where      
	 *      <p><i>m</i> = # rows (observations) in input data
	 *      <p><i>d</i> = # PC dimensions (eigenvectors).
	 * </ul>
	 * @see getEvects()
	 * @see getPCScores(int[])
	 * @see rotateEigimsToOrigSpace(int[], Matrix, Matrix)
	 */
	public Matrix getPCScores() {
		return pcScores;		
	}
	
	
	/** Returns PC scores corresponding to PC dimensions specified in input 
	 * array.
	 * @param pcDims 
	 * 			Array of indices (0-relative) indicating which PC dimensions to include.
	 * 			<p> E.g., if input array is {0,1,4,5}, returns 
	 * 			scores for 1st, 2nd, 5th and 6th PC dimensions.
	 * 
	 * @return PC scores
	 * <ul> 
	 * 		<li>Matrix size: <i>m</i> X <i>k</i>, where
	 * 		<p><i>m</i> = # rows (observations) in input data
	 *      <p><i>k</i> = # of elements in given pcDims array
	 * </ul>
	 * 
	 * @see getPCScores()
	 * @see getEvects()
	 * @see rotateEigimsToOrigSpace(int[], Matrix, Matrix)
	 */
	
	public Matrix getPCScores(int[] pcDims) {
		return pcScores.subMatrixCols(pcDims);
	}
	
//	// Set method(s): used in Npairsj.rotateToOrigSpace(PCA)
//	public void setEvects(Matrix newEvects) {
////		if	((newEvects.numCols() != eigenvectors.numCols())) {
////			throw new IllegalArgumentException("Input Matrix and PCA eigenvectors "
////					+ "must have same number of columns (PC Dims)");
////		}
//		eigenvectors = newEvects;
//	}
	
	/** Indicates whether PCA results live in the same space as original data 
	 * before any transformations - i.e., whether PCA eigenvectors
	 * are represented in terms of original data variables (e.g., voxels). If original
	 * data consists of images, PC eigenvectors in original space may be considered
	 * 'eigenimages'.
	 * 
	 * @return 
	 * 		<p> <code>true</code> if PC eigenvectors are represented in original space
	 * 		<ol>This may happen in 2 cases:
	 * 			<li>Data used as PCA input lives in original space;
	 * 			<li>Data used in PCA lives in some other space but PC eigenvectors have 
	 * 				been transformed back to original space after PCA calculation.			
	 * 		</ol>
	 * 
	 * 		<p> <code>false</code> if PCA input data has been previously transformed into a 
	 * 				different (typically dimension-reduced) space and no rotation back into original
	 * 				space has been performed on PC eigenvectors after calculation.
	 */
	public boolean pcaInOrigSpace() {
		return pcaInOrigSpace;
		
	}
	
	
//	// Truncate pc scores and eigenvalues (and associated evalMat) so only input dims are
//	// included.  
//	private void truncateEvalsAndScores(int[] pcDims) {
//		pcScores = pcScores.subMatrixCols(pcDims);
//		eigenvalues = MLFuncs.getItemsAtIndices(eigenvalues, pcDims);
//		evalMat = evalMat.subMatrixCols(pcDims);
//	
//	}
	
		
	/** This method returns a PCA with only input PC eigenvector dimensions rotated back into 
	 * original space. The rest of the dimensions are excluded; otherwise transformation would 
	 * take a prohibitive length of time to calculate for large data.
	 * 
	 *  <p>If this PCA is already in original space, returned PCA contains subset of this 
	 *  PCA's eigenvectors.
	 * 
	 * <p>NOTE: Returned PC scores and eigenvalues are also included only for input PC dimensions.
	 * 
	 * <p>This PCA is not modified.
	 * 
	 * @param pcDims
	 * @param invXformedData <p>For regular EVD, this matrix is inverse of transformed data matrix 
	 * 						(i.e., inverse of M, the matrix used as PCA input)
	 *                      <p>TODO: update documentation explaining invXformedData for
	 *                        unweighted EVD  
	 * @param origData	data matrix in original space
	 * @return PCA containing transformed eigenvectors. Returned PCA contains only input 
	 * pcDims dimensions.
	 * 
	 */
	public PCA rotateEigimsToOrigSpace(int[] pcDims, Matrix invXformedData, Matrix origData) {

		if (pcaInOrigSpace) {
			// return PCA containing only input PC dimensions of this PCA
			return new PCA(getEvects(pcDims), MLFuncs.getItemsAtIndices(eigenvalues, pcDims),
					 getPCScores(pcDims), true);
		}
		
		// else return PCA (also containing only input PC dimensions) with rotated eigenvectors
		
		/** Assume that feature selection was done via SVD.
		 * Let original data D = WSYt. 
		 * 
		 * Then transformed data M = DY 
		 * and
		 * feature-selected (dimension-reduced) data M* = DY* (where Y* is a subset of columns of Y). 
		 * 
		 * PCA was done on M*. 
		 * 
		 * Transform PCA eigenvectors via inverted feature-selection 
		// transformation Matrix Y (i.e., via feature-selection 
		 * eigenvectors, a.k.a. "eigenimages", Y).
		 * <p>Note that Y^-1 = Yt.
		 * 

		// M = DY 
		 * ==> MYt> = D 
		 * ==> Yt = M^-1D.
		// Yt ( = SVD eigenimages, e.g.,) is very large to store, hence
		// we store M^-1 ( = inverse feature selection data )
		 * instead to reconstruct Yt on the fly:
		 * 
		// Yt = M^-1D 
		 * ==> MYt = MM-1D 
		 * ==> M1Vt = M1M^-1D = D1
		// where 
		 * 
		 * D1 = PCA eigenvectors in rows, in this case, and 
		 * M1 = PCA eigenvectors in rows transformed back into original voxel space.
		// 
		 * Given 
		 * 
		 * size D1 is reducDataDims X reducDataDims, 
		// size D^-1 is reducDataDims X data.numRows(), and
		// size M is data.numRows() X origDataDims,
		// 
		 * want M1t = origDataDims X reducDataDims().
		 * 
		// TODO: Want to retain only significant dimensions of PCs! 
		// (see IDL code for comparison)
		 */

		// size D1(invD) is reducDataDims X data.numRows(): 	

		Matrix D1invD = getEvects(pcDims).transpose().mult(invXformedData);

		// instead of multiplying D1 by invDM 
		// == (reducDataDims X reducDataDims) * (reducDataDims X origDataDims),
		// multiply D1(invD) by M
		// == (reducDataDims X data.numRows()) * (data.numRows() X origDataDims)
		// == fewer calculations

		Matrix eigenimages = D1invD.mult(origData); // eigenimages in rows
		
		return new PCA(eigenimages.transpose(), MLFuncs.getItemsAtIndices(eigenvalues, pcDims), 
				getPCScores(pcDims), true);		
	}
	

	/**
	 * Multiplies values in input double array by input scalar and returns new
	 * double array containing result.
	 * 
	 * @param values -
	 *            Array of values to be scaled by scalar
	 * @param scalar -
	 *            Factor by which to multiply elements in 'values'
	 * @return scaled array
	 */
	private static double[] scale(double[] values, double scalar) {
		double[] scaled_vals = new double[values.length];
		for (int i = 0; i < values.length; ++i) {
			scaled_vals[i] = values[i] * scalar;
		}
		return scaled_vals;
	}
}
