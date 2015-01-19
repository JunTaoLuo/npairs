package npairs.utils;

import npairs.shared.matlib.*;
import npairs.Npairs;
import npairs.NpairsException;
import npairs.io.*;
import pls.shared.MLFuncs;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

import java.util.*;
import java.io.IOException;
import java.io.PrintStream;

/**
 * <p>Performs Canonical Variates Analysis (CVA). 
 * 
 * <p> CVA is the multivariate version of Fisher's Linear Discriminant Analysis (LDA) technique 
 * whereby discrimination between groups or classes is optimized by finding the linear function 
 * maximizing the ratio of between-class covariances to within-class covariances. 
 * 
 * The CVA algorithm calculates eigenvectors and eigenvalues of 
 * <p>W<sup>-1</sup> * B</p> 
 * where
 * <ul> 
 * <li>W is the within-class sums-of-squares-and-products (SSP) 
 * matrix, and 
 * <li>B is the between-class SSP matrix. 
 * </ul>
 * <p>(An SSP matrix is an 'unscaled' 
 * covariance matrix, i.e., equivalent to a covariance matrix except for scaling 
 * by 1/n or 1/(n-1) where n is the number of samples.)
 * <p>See <i>Multivariate Analysis</i> (Mardia et al., 1979) for details.
 * 
 * <p>This class also produces the following metrics:
 * <ul>
 * <li>RMS (root mean-squared) error. See {@link #getRMSError() getRMSError}.</li>
 * <li>R<sup>2</sup> values. See {@link #getR2() getR2}. (optional)</li>
 * <li> Chi-squared values. See {@link #getChiSqrInfo() getChiSqrInfo}. (optional)</li>
 * </ul>
 * 
 * <p>REQUIRED: Set npairs.shared.matlib.Matrix_Impl.matlibType outside of this
 * class.
 * @see npairs.shared.matlib.Matrix
 * 
 * @author Anita Oder
 * 
 */
public class CVA {

	private boolean debug = false;

	protected double rmsError;

	protected Matrix cvaEvalMat;

	protected double[] cvaEvals;

	protected Matrix cvaEvectsSmall; // cva eigenimages before projection back into any other basis

	protected Matrix cvaScores;   // contains projection of input data onto cv
								// eigenimages,
								// i.e., representation of input data in cv space;
								// num rows = num rows in input data;
								// num cols = num CV dims
	protected Matrix cvaEigimsBig = null; // cva eigenimages projected back into
	                                 // 'original' space (if init. feature selection
	                                 // step was done before second-level basis
	                                 // decomposition step, this data will have
	                                 // been projected back through 2 transformations
	                                 // to get back into original space)
	
	protected Matrix cvaFeatSpEigims = null; // cva eigenimages projected back once, through basis 
	                          // decomposition space (e.g. PCA), but not through
	                          // init. feature selection step. 
	

	protected Matrix chiSqrInfo;

//	private double[][] confidReg;

	protected int nSigDim;
	
	private Matrix r2; // contains r^2 values calculated for each input data dim and CV dim
	                   // size(r2) = no. input data dims (e.g. no. PC dims) rows X 
	                   // no. CV dims cols

	protected int[] classLabels;

	protected int nCVDim;

	protected int nClasses;

	protected int[] clsSz; // contains class sizes for each class label (class labels
						 // in ascending order)

	protected int nVol;
	
	private Matrix W;
	private Matrix B;
	
	private boolean computeR2; // Recommended: set to false unless you want 
	                                   // to examine R2 output (since it's a non-trivial
	                                   // exercise to compute R2 values).
	private boolean computeChiSqr; // recommended: set to false unless
	                                        // you want to examine chi squared output 
	                                        // (since it's non-trivial to calculate)
//	protected boolean cvaInOrigSpace; // true if cva eigims live in original (image) space
	
	protected boolean logInfo = false; // true if info to be logged in logStream
	
	protected PrintStream logStream; // where to log info 
	
	/**
	 * Constructor for CVA object. Performs CVA. Input Matrix
	 * is NOT modified.
	 * 
	 * @param data 
	 *          Matrix with rows = observations, columns = variables
	 * @param labels 
	 *          Array (length # rows in data Matrix) containing class
	 *            labels for corresponding rows of data Matrix. Each class must
	 *            consist of at least 2 samples (rows).
	 * @param computeR2 
	 * 			Indicates whether to calculate R<sup>2</sup> stats. 
	 *             See {@link #getR2() getR2}.
	 * @param computeChiSqr 
	 *          Indicates whether to calculate Chi-squared stats.
	 * 				Threshold is set to 95%. See {@link #getChiSqrInfo 
	 *             getChiSqrInfo}.
	 * @param pcaEvects 
	 * 			If input data is in PCA space, i.e., if PCA was run on 
	 *             data before doing CVA, pcaEvects contains the 
	 *             eigenvectors from the PCA, i.e., the PCA basis vectors
	 *             of this input data. pcaEvects are used to generate CV eigenimages 
	 *             from CVA results by projecting CV eigenvectors back from PCA space.
	 *             <p> Set pcaEvects to null if no PCA was run.
	 * @param logStream 
	 * 			Where to print log information. Set to null to disable
	 *             logging info.
	 * 
	 */
	public CVA(Matrix data, int[] labels, boolean computeR2, boolean computeChiSqr, 
			Matrix pcaEvects, PrintStream logStream) throws NpairsException {
		if (logStream != null) {
			this.logInfo = true;
			this.logStream = logStream;
		}
		this.computeR2 = computeR2;
		this.computeChiSqr = computeChiSqr;
	
		if (debug) {
			System.out.println("CVA input labels (length = " + labels.length + "): ");
			NpairsIO.print(labels);
			System.out.println("Size input data: " + data.numRows() + " X " + data.numCols());
		}
			
		if (data.numRows() != labels.length) {
			throw new IllegalArgumentException(
					"No. rows in input Matrix does not match "
							+ "length of input array");
		}

		final double CHISQR_THRESH = 0.95;
		classLabels = labels;
		nVol = classLabels.length;
		
//		double sTime = System.currentTimeMillis();
		computeCVA(data);
//		this.cvaInOrigSpace = dataInOrigSpace;
//		if (debug) {
//			double tTime = (System.currentTimeMillis() - sTime) / 1000;
//			System.out.println("Total time CVA: " + tTime + " s");
//		}
		
		if (computeR2) {
			r2 = computeR2(data);
		}
		if (computeChiSqr) {

			try {
				if (logInfo) {
					logStream.print("\tCalculating chi sqr...");
				}
				double sTime = System.currentTimeMillis();
				chiSqrInfo = chiSqrInfo(CHISQR_THRESH);
				double tTime = (System.currentTimeMillis() - sTime) / 1000;
				if (logInfo) {
					logStream.println("[" + tTime + " s]");
				}

				//				System.out.println("Calculating confidence regions...");
				//				sTime = System.currentTimeMillis();
				//			confidReg = confidReg(chiSqrInfo.getRow(1), CHISQR_THRESH);
				//			if (debug) {
				//				double tTime = (System.currentTimeMillis() - sTime) / 1000;
				//				System.out.println("Finished calculating conf. reg.: " + tTime);
				//			}
			} 
			catch (MathException e) {
				throw new NpairsException("Problem encountered calculating chi squared info for CVA: " 
						+ e.getMessage());
			}
		}
		
		// Rotate eigenvectors back from pca space (if run) to create
		// CVA eigenimages. Eigenimages will generally be in initial 
		// feature-selection space after this step. If no feature-selection 
		// step was performed, eigenimages will be in original voxel space.
		long sTime = System.currentTimeMillis();
		if (logInfo) {
			logStream.print("Creating CVA eigenimages... ");
		}
		// If pca wasn't run, pcaEvects are null so CVA eigenimages
		// will be set equal to CVA eigenvectors.
		createFeatSpEigenimages(pcaEvects);
		if (logInfo) {
			long tTime = (System.currentTimeMillis() - sTime) / 1000;
			logStream.println("["+ tTime + " s]");
		}
	}

	private Matrix computeR2(Matrix data) {
		// for each CV dim, find correlation coefficients between CV scores and 
		// 'timecourse' for each input data dim (e.g. PC or voxel dim) 
		return data.correlate(cvaScores).squareElements();
	
	}

	// Super constructor for JavaQuadCVA, MatlabQD (for testing only?) 
	protected CVA() {};
	
	/**
	 * Performs CVA on input data. Input data is NOT modified.
	 * 
	 * @param data
	 * @throws NpairsException
	 */
	private void computeCVA(Matrix data) throws NpairsException {
	
		Matrix W = getW(data, classLabels);
//		if (debug) {

			EigenvalueDecomposition evd = W.eigenvalueDecomposition();
			double[] evals = evd.getRealEvals();
			
			double cond = MLFuncs.max(evals) / MLFuncs.min(evals);
			
			if (debug && logInfo) {
				logStream.print("(W cond. no. =  ");
				logStream.printf("%.3f) ", cond);
			}
			if (Math.abs(cond) > 1000 && logInfo) {
				logStream.printf("\nt*WARNING: W is nearly " 
						+ "singular! (W cond. no. = %.3f)* ", cond);
			}
//		}
		
		Matrix B = getB(data, W);


		// calculate W^(-1) * B

		Matrix invWTimesB = null;
		try { 
			invWTimesB = W.inverse().mult(B);
		} catch (IllegalArgumentException e) {
			throw new NpairsException("Within-class covariance matrix W is singular! " +
					"\nMaybe you included too many PC dimensions in the CVA analysis." +
					"\nSee analysis log file in results directory for more details.");
		}

		// calculate CVA eigenvalues/eigenvectors

		EigenvalueDecomposition evDecomp = invWTimesB.eigenvalueDecomposition();
		// TODO: verify that eigenvalues of W^(-1) * B are real and >= 0, 
		// so can ignore imaginary part
		
//		if (debug) {
//			tTime = (System.currentTimeMillis() - sTime) / 1000;
//			System.out.println("Time EVD of W^(-1)*B: " + tTime + " s");
//		}
		nCVDim = Math.min(nClasses - 1, data.numCols());
		cvaEvalMat = evDecomp.getRealEvalMat();
		cvaEvectsSmall = evDecomp.getEvects();

		// Calculate RMS error of |MV - VL|, where M is W^(-1) * B,
		// V is evects, L is evals
		rmsError = rmsError(invWTimesB, cvaEvectsSmall, cvaEvalMat);

		// trim evals/evects 
		cvaEvals = evDecomp.getRealEvals(nCVDim);
		cvaEvalMat = evDecomp.getRealEvalMat(nCVDim);
		cvaEvectsSmall = evDecomp.getEvects(nCVDim);

		// normalize eigenvectors so each group (class) has variance 1
		// and is uncorrelated to all other groups.
		// (compare this.normEvectsByLength())
		normEvectsByVar(W);

		// calculate canonical scores
		cvaScores = data.mult(cvaEvectsSmall);

	}

    /**	 Normalize eigenvectors so each group has variance 1 when represented
     *   in the eigenvector basis and is uncorrelated with all other groups, i.e.,
     *   normalize V so that it is the linear transformation that takes W 
     *   (normalized by degrees of freedom n-g) to the identity matrix I, i.e.:
     *   <p>given eigenvector matrix V, want V<sup>t</sup>(W/(n-g))V = I,
     *    where 
     *    <ul>
     *    <li>n = # observations (timepoints), and 
     *    <li>g = # classes;
     *    </ul>
     *    see e.g.(i) Strother et al. Neuroimage 15, 747-771 (2002) (p. 769),
     *            (ii) Mardia et al., 1979,p. 343
     *    
     * @param W within-class covariance matrix (actually sum of 
     *          individual class SSP Matrices) 
     * @see #normEvectsByLength()
     */
	private void normEvectsByVar(Matrix W) {
		
		int nObs = classLabels.length;
		Matrix scaledW = W.mult(1. / (nObs - nClasses));

		for (int i = 0; i < cvaEvectsSmall.numCols(); ++i) {
			Matrix currEvect = cvaEvectsSmall.subMatrixCols(new int[] { i });
			double d = currEvect.transpose().mult(scaledW).mult(currEvect).
				getQuick(0, 0);
			if (d > 0) { 
				// scale evect
				currEvect = currEvect.mult(1 / Math.sqrt(d));
				cvaEvectsSmall.setColumnQuick(i, currEvect.getColumnQuick(0));
			} 
			else { 
				// d = 0; set evect to 0.  
				for (int j = 0; j < cvaEvectsSmall.numRows(); ++j) {
					cvaEvectsSmall.setQuick(j, i, 0);
				}
			}
		}
	}

	/**	 Normalize eigenvectors so each has length 1.
	 * @see #normEvectsByVar(Matrix)
	 *
	 */
	@SuppressWarnings("unused")
	private void normEvectsByLength() {
		
		 for (int i = 0; i < cvaEvectsSmall.numCols(); ++i) {
		 double[] currEvect = cvaEvectsSmall.getColumn(i);
		 double sumOfSq = 0.0;
				
		 for (double val : currEvect) {
		 sumOfSq += (val * val);
								
		 }
							
		 for (int j = 0; j < currEvect.length; ++j) {
		 currEvect[j] = currEvect[j] / Math.sqrt(sumOfSq);
		 }
		 cvaEvectsSmall.setColumn(i, currEvect);
								
		 }
	}

	/**
	 * Returns RMS (root mean-squared) error of |MV - VL|, where 
	 * <ul>
	 * <li>M is W<sup>-1</sup> * B, 
	 * <li> V is eigenvectors, 
	 * <li> L is eigenvalues. 
	 * </ul>
	 * <p>REQUIRED: Matrix M must be square.
	 * 
	 * @param M
	 *            data matrix
	 * @param evects
	 *            V matrix
	 * @param evals
	 *            L matrix
	 * @return rms error
	 */
	private static double rmsError(Matrix M, Matrix evects, Matrix evals) {
		double rmsErr = 0;
		int numDims = M.numCols();

		// calculate MV - VL (result lives in diffMat)
		Matrix diffMat = M.mult(evects);
		diffMat = diffMat.plusEquals(evects.mult(evals).mult(-1));

		for (int dim = 0; dim < numDims; ++dim) {
			double sumOfSqrdDiffs = 0;
			for (int row = 0; row < numDims; ++row) {
				sumOfSqrdDiffs += Math.pow(diffMat.getQuick(row, dim), 2);
			}
			rmsErr += Math.sqrt(sumOfSqrdDiffs / numDims);
		}
		// get avg rmsErr over all dims
		rmsErr = rmsErr / numDims;
		return rmsErr;
	}
	
	/** Returns within-class SSP Matrix. (Actually SUM of individual class SSP 
	 * Matrices; note that this is equivalent in variance structure to the 
	 * *average* class SSP matrix.) Note individual class SSP matrices calculated
	 * on column-mean-centred individual class data.
	 * @return W - within-class SSP Matrix
	 */
	//NOTE: not used in PLSNPAIRS
	public Matrix getW() {
		return W;
	}

	/** Calculates and returns within-class sums-of-squares-and-products (SSP) Matrix 
	 *  (actually SUM of individual class SSP Matrices - NOTE that this 
	 *  is equivalent in variance structure to the *average* group SSP mat)).
	 *  Matrix of individual class data column-mean-centred before SSP calculated.
	 * @param data
	 * @param classLabels 
	 * @return W - within-class SSP Matrix
	 */
	protected Matrix getW(Matrix data, int[] classLabels) throws NpairsException {
		
		Hashtable<Integer, int[]> clsIndices = MLFuncs.getLabelIndices(classLabels);		
        //		 TODO: Using pls MLFuncs code here but not currently in
		//       getLabelIndices for same function 
		int[] sortedUniqClasses = MLFuncs.sortAscending(MLFuncs
				.unique(classLabels));
		nClasses = sortedUniqClasses.length;
		clsSz = new int[nClasses];
	    int nCols = data.numCols();
		
		// calculate W
//		Matrix W;
		try {
			W = (new MatrixImpl(nCols, nCols).getMatrix());

			for (int g = 0; g < nClasses; ++g) {
				int[] currClsIndices = clsIndices.get(sortedUniqClasses[g]);
				clsSz[g] = currClsIndices.length;
				if (clsSz[g] < 2) {
					throw new NpairsException(
					"Each class must contain at least 2 samples");
				}

				// compute SSP mat for current grp and add it to W
				Matrix currClsData = data.subMatrixRows(currClsIndices);
				currClsData = currClsData.meanCentreColumns();

				Matrix currSSPMat = currClsData.sspByCol();
                W = W.plusEquals(currSSPMat);
			}	
		} 
		catch (NullPointerException e) {
			throw new NpairsException(e.getMessage());
		}
		
		return W;
	}

	/** Returns between-class SSP Matrix B.
	 *  <p>B = T - W where 
	 *  <ul>
	 *  <li>T = total SSP Matrix (i.e., SSP Matrix 
	 *  for entire input data Matrix), and
	 *  <li>W = within-class SSP Matrix
	 *  <li> (Note: T and W SSP matrices calculated
	 *  on column-mean-centred data)
	 *  </ul>
	 * @return B - between-class SSP Matrix
	 */
	//NOTE: not used in PLSNPAIRS
	public Matrix getB() {
		return B;
	}
	
	/** Calculates and returns between-class SSP Matrix B
	 *  (B = T - W where T == total SSP Matrix and
	 *   W == within-class SSP Matrix; note that W is 
	 *   actually SUM of individual class SSP matrices.
	 *  This method of calculating B and W is exactly the
	 *  same as the method used in NPAIRS-2.5 cva_wb_mats.pro.
	 *  
	 * Input data Matrix is NOT modified.
	 * @param data
	 * @param W
	 * @return B
	 */
	protected Matrix getB(Matrix data, Matrix W) {
		
		data = data.meanCentreColumns(); 

//		double sTimeT = System.currentTimeMillis();
		Matrix T = data.sspByCol();
//		if (debug) {
//			double tTimeT = (System.currentTimeMillis() - sTimeT) / 1000;
//			System.out.println("Time calculating SSP of CVA input data: "
//					+ tTimeT + " s");
//		}

		// calculate B (between-groups SSP matrix) = T - W
//		Matrix B = T.copy();
		B = T.copy();
		B = B.plusEquals(W.mult(-1));
		return B;
	}

//	private Matrix getInvWTimesB() {
//		return invWTimesB;
//	}

	/**Returns RMS (root mean-squared) error, used to determine goodness of fit
	 * of eigenvalue/eigenvector estimates.
	 * @return RMS (root mean-squared) error of |MV - VL|, where 
	 * <ul>
	 * <li>M = W<sup>-1</sup> * B (see {@link #getW()}, {@link #getB()}),
	 * <li>V =  CV eigenvectors (see {@link #getEvects()}, but note that eigenvectors 
	 *  are only normalized *after* RMS error has been calculated, i.e., V here is not 
	 *  normalized), 
	 *  <li>L = CV eigenvalue matrix (see {@link #getEvalMat()}).
	 *  </ul> 
	 */
	// NOTE: not used in PLSNPAIRS
	public double getRMSError() {
		return rmsError;
	}

	/**Returns CV eigenvalue Matrix.
	 * 
	 * @return Diagonal Matrix containing CV eigenvalues down the main diagonal.
	 *         <p>Eigenvalues are in descending order unless current 
	 *         CVA results have been transformed to match reference CVA (see 
	 *         {@link #matchToRef(CVA, boolean, boolean)}).
	 *         <p>Matrix size: # CV dimensions X # CV dimensions
	 * @see #getEvals()
	 */
	// NOTE: not used in PLSNPAIRS
	public Matrix getEvalMat() {
		return cvaEvalMat;
	}

	/** Returns array of CV eigenvalues.
	 * 
	 * @return Array containing CV eigenvalues. 
	 *         <p>Eigenvalues are in descending order unless
	 *         current CVA results have been transformed to match
               reference CVA (see {@link #matchToRef(CVA, boolean, boolean)}).
	 *         <p>Array length: # CV dimensions
	 * @see #getEvalMat()
	 */
	public double[] getEvals() {
		return cvaEvals;
	}

	/** Returns CV scores, i.e.,
	 * <ul>
	 * <li>projection of input data onto CV
     *  eigenvectors;</li>
     *  <li>representation of input data in CV space.</li>
     *  </ul>
	 * Note that eigenvectors are normalized before CV scores are 
	 * calculated (see {@link #getEvects()} for details about normalization).
	 * @return Matrix containing CV scores, one column per CV dimension. 
	 *         <p>Matrix size: # rows (observations) in input data
	 *                             X # CV dimensions 
	 *         <p>Order of CV dimensions corresponds to order
	 *            of eigenvalues (see {@link #getEvals()}, {@link #getEvalMat()}).
	 * @see #avgCVScores()
	 */
	public Matrix getCVScores() {
		return cvaScores;
	}

	/** Returns CV eigenvectors, i.e., CV "eigenimages" before
	 *  projection back into any other basis.
	 * <p>Note -- Eigenvectors are normalized so each class has variance 1:
     *   <p>Given eigenvector matrix V and within-class SSP matrix W
     *   (see {@link #getW()}), V is normalized to satisfy
     *   <p>V<sup>t</sup> * W/(n-g) * V = I,
     *    where 
     *    <ul>
     *    <li>n = # observations (timepoints), and 
     *    <li>g = # classes;
     *    </ul>
     *    (see, e.g., Strother et al., Neuroimage 2001)
	 * @return Matrix containing normalized CV eigenvectors, one per column.
	 *        <p>Matrix size: # dimensions (columns) in input data
	 *                      X # CV dimensions
	 *         <p>Order of CV dimensions corresponds to order
	 *            of eigenvalues (see {@link #getEvals()}, {@link #getEvalMat()}).
	 */
	// NOTE: Not used in PLSNPAIRS
	public Matrix getEvects() {
		return cvaEvectsSmall;
	}

	/** Returns CV eigenimages after projection back into original 
	 *  (voxel) image basis, if data underwent initial feature selection
	 *  step. Must first create these eigenimages explicitly
	 *  by calling {@code rotateEigimsToOrigSpace(Matrix, Matrix)} before
	 *  calling this method. 
	 * 
	 * @return Matrix containing CV eigenimages, one per column.
	 *         <p>Matrix size: # dimensions (voxels) in original data
	 *                      before any basis changes 
	 *                      X # CV dimensions
	 *         <p>Order of CV dimensions corresponds to order
	 *            of eigenvalues (see {@link #getEvals()}, {@link #getEvalMat()}).
	 *         <p> Returns null if eigenimages in original space have not yet been 
	 *             created via {@code rotateEigimsToOrigSpace(Matrix, Matrix)}. 
	 *         <p>Note that this will be the case even if no initial feature selection step
	 *             was done, i.e., {@code rotateEigimsToOrigSpace(Matrix, Matrix)} must
	 *             always be called first to initialize return value for this method.
	 *             
	 * @see #rotateEigimsToOrigSpace(Matrix, Matrix)
	 * @see #getFeatSpEigims()
	 */
	public Matrix getEigimsBig() {
		return cvaEigimsBig;
	}
	

	/** Returns CV eigenimages after projection back through any 
	 *  basis changes (e.g., PCA) performed after initial feature selection step
	 *  and before CVA. If initial feature selection was performed,
	 *  eigenimages returned by this method will be in feature selection
	 *  space; if no initial feature selection was performed, eigenimages
	 *  returned by this method will be the same as those returned by 
	 *  {@code getEigimsBig()}, i.e., they will be in original 
	 *  data space.
	 * 
	 * @return Matrix containing CV eigenimages, one per column.
	 *         <p>Matrix size: # dimensions in data after optional 
	 *            initial feature-selection ('data-reduction')
	 *            step was performed on data X # CV dimensions
	 *         <p>Order of CV dimensions corresponds to order
	 *            of eigenvalues (see {@link #getEvals()}, {@link #getEvalMat()}).
	 *
	 * @see #rotateEigimsToOrigSpace(Matrix, Matrix)
	 * @see #getEigimsBig()
	 */
	public Matrix getFeatSpEigims() {
		return cvaFeatSpEigims;
	}

	/** Returns Chi-squared data for each CV dimension (+1) calculated using 
	 * algorithm of Jon Anderson's IDL NPAIRS code 'cva_sig_dims.pro'. I do 
	 * not know what information in the extra dimension represents.
	 * 
	 * @return Matrix chiSqrInfo, or null if info not calculated
	 * (i.e., if computeChiSqr set to False when CVA constructor called). 
	 * <p>Matrix size: 3 rows X (# CV dims + 1) columns.
	 * <li>chiSqrInfo(0, i) = Chi-squared value (Bartlett's statistic) for ith 
	 *        CV dimension
	 * <li>chiSqrInfo(1, i) = probability associated with corresponding Chi-squared
	 *         value 
	 * <li>chiSqrInfo(2, i) = degrees of freedom
	 * <p>Order of CV dimensions corresponds to order
	 *    of eigenvalues (see {@link #getEvals()}, {@link #getEvalMat()}).
	 * 

	 */
	//NOTE: not used in PLSNPAIRS
	public Matrix getChiSqrInfo() {
		return chiSqrInfo;
	}

//	public double[][] getConfidReg() {
//		return confidReg;
//	}

	// nSigDim only calculated if confidReg calculated
//	public int getNumSigCVDim() {
//		return nSigDim;
//	}

	/**
	 * Creates CVA eigenimages by projecting CVA eigenvectors onto input basis
	 * vectors (e.g., PCA eigenvectors, if PCA was run on data before passing it
	 * to CVA).
	 * <p>NOTE Dimensions of eigenimage Matrix will be:
	 * <p># data dims (i.e., # voxel or initial feature space variables) X # CVA dims
	 * 
	 * @param basisVects 
	 *            <p>Matrix of vectors onto which to project eigenvectors (columns =
	 *            vectors). 
	 *            <p>Set this argument to null if CVA already in initial feature selection 
	 *            or original voxel space, e.g., if no PCA was performed before CVA.
	 *            In this case, CVA eigenimages will be set to equal CVA 
	 *            eigenvectors.
	 *          
	 *            
	 * @see #getFeatSpEigims()
	 * @see #getEigimsBig()
	 * @see #rotateEigimsToOrigSpace(Matrix, Matrix)
	 */
	private void createFeatSpEigenimages(Matrix basisVects) {
		if (basisVects != null) {
			cvaFeatSpEigims = basisVects.mult(cvaEvectsSmall);
//			cvaInOrigSpace = basisVectsInOrigSpace;
		} else {
			cvaFeatSpEigims = cvaEvectsSmall;
		}
	}

	/** Not fully implemented! Applies given full Procrustes transformation to 
	 * CV scores and eigenvectors (and eigenimages in all available vector bases) 
	 * in this CVA, but still need to determine how the full transformation may be 
	 * applied to CV eigenvalues, Chi squared info &  R2 info.
	 * @param proc
	 * <p> Procrustes transformation to apply to this CVA data.
	 * @see npairs.utils.Procrustes
	 */
	private void applyFullProcrust(Procrustes proc) {
		
		// CV eigenimages: apply proc Xform
		cvaEvectsSmall = cvaEvectsSmall.mult(proc.getRot().transpose());
//		if (cvaFeatSpEigims != null) {
			cvaFeatSpEigims = cvaFeatSpEigims.mult(proc.getRot().transpose());
//		}
		if (cvaEigimsBig != null) {
			cvaEigimsBig = cvaEigimsBig.mult(proc.getRot().transpose());
		}
		
		// CV scores: apply proc Xform
		cvaScores = cvaScores.mult(proc.getRot().transpose());
		
		// CV evals, chi squared, R2: ??
		
	}

	/**
	 * Negates and permutes CV dims of this CVA result according to negation and
	 * permutation info supplied in input args. 
	 * <p>Fields affected: 
	 * <li>Eigenvalues - permuted; 
	 * <li>CV scores - negated and permuted; 
	 * <li>CV eigenimages - negated and permuted; 
	 * <li>Chi-squared info - permuted; 
	 * <li>R2 info - permuted
	 * 
	 * @param sign 
	 *            <p>Array consisting of 1's and -1's. 
	 *            <p>If ith permuted CV dim
	 *            is to be negated, sign[i] = -1; if not, sign[i] = 1.
	 * @param permutation 
	 *            <p>Array indicating new order of CV dims. 
	 * 
	 *            <p>E.g.: permutation = {0, 2, 1}, sign = {-1, 1, -1} means:
	 *            <p> newScores[0] = -1 *
	 *            oldScores[0] 
	 *            <p>newScores[1] = 1 * oldScores[2] 
	 *            <p>newScores[2] = -1 *
	 *            oldScores[1]
	 * 
	 * @throws NpairsException
	 *             if length of input args != # CV dims
	 */
	private void negateAndPermuteDims(int[] sign, int[] permutIndex)
			throws NpairsException {

		if (sign.length != nCVDim || permutIndex.length != nCVDim) {
			throw new NpairsException(
					"Input arrays must have length = no. of CV dims");
		}

		double[] tmpEvals = new double[nCVDim];
		for (int dim = 0; dim < nCVDim; ++dim) {
			tmpEvals[dim] = cvaEvals[permutIndex[dim]];
		}
		cvaEvals = tmpEvals;

		try {
			cvaEvalMat.setDiag(cvaEvals);
		}
		catch (MatrixException me) {
			// matrix is square by construction
		}
		
		cvaScores = cvaScores.permuteColumns(permutIndex);
		cvaEvectsSmall = cvaEvectsSmall.permuteColumns(permutIndex);
//		if (cvaFeatSpEigims != null) {
			cvaFeatSpEigims = cvaFeatSpEigims.permuteColumns(permutIndex);
//		}
		if (cvaEigimsBig != null) {
			cvaEigimsBig = cvaEigimsBig.permuteColumns(permutIndex);
		}
		
		if (computeR2) {
			r2 = r2.permuteColumns(permutIndex);
		}

		// TODO: Figure out how to handle extra chi-squared dim
		// (note: in idl cva_reference.pro, extra dim is ignored;
		// permutation is applied as if chi-squared info only has
		// no. CV dims dimensions.

		if (computeChiSqr) {
			int[] chiSqrPermutIndex = new int[nCVDim + 1];
			for (int i = 0; i < nCVDim; ++i) {
				chiSqrPermutIndex[i] = permutIndex[i];
			}
			chiSqrPermutIndex[nCVDim] = nCVDim;

			chiSqrInfo = chiSqrInfo.permuteColumns(chiSqrPermutIndex);
		}

		for (int dim = 0; dim < nCVDim; ++dim) {
			double[] newSignScores = MLFuncs.product(cvaScores.getColumn(dim),
					sign[dim]);
			cvaScores.setColumn(dim, newSignScores);
			double[] newSignEvects = MLFuncs.product(cvaEvectsSmall.getColumn(dim),
					sign[dim]);
			cvaEvectsSmall.setColumn(dim, newSignEvects);
//			if (cvaFeatSpEigims != null) {
				double [] newSignEigims = MLFuncs.product(cvaFeatSpEigims
						.getColumn(dim), sign[dim]);
				cvaFeatSpEigims.setColumn(dim, newSignEigims);
//			}
			if (cvaEigimsBig != null) {
				double[] newSignEigimsBig = MLFuncs.product(cvaEigimsBig
						.getColumn(dim), sign[dim]);
				cvaEigimsBig.setColumn(dim, newSignEigimsBig);
			}
		}
	}

	/**
	 * Calculates chi squared data for each CV dimension (+1?) using algorithm
	 * of Jon Anderson's IDL npairs code 'cva_sig_dims.pro'.
	 * 
	 * @param thresh
	 * @return Matrix chiSqrInfo 3 row X (num CV dims + 1) col. chiSqrInfo(0,i) =
	 *         chi squared value (Bartlett's statistic) for ith cv dimension
	 *         chiSqrInfo(1,i) = prob associated with corresponding chi squared
	 *         value chiSqrInfo(2,i) = degrees of freedom
	 */

	// TODO: figure out why there is one more dim than no. CV Dims and how to
	// handle
	// it - see CVA.negateAndPermuteDims(int[], int[])
	private Matrix chiSqrInfo(double thresh) throws MathException {

		int nPC = cvaEvals.length;

		// compute D statistic
		double[][] chiSqrInfo = new double[3][nCVDim + 1];
		for (int i = 0; i <= nCVDim; ++i) {
			for (int j = i + 1; j <= nPC; ++j) {
				chiSqrInfo[0][i] += Math.log(1 + cvaEvals[j - 1]);
			}
			chiSqrInfo[0][i] = (nVol - 1 - ((nPC + nClasses) / 2.0))
					* chiSqrInfo[0][i];
			chiSqrInfo[2][i] = (nPC - i) * (nClasses - i - 1);
			if (chiSqrInfo[2][i] > 0) {
				chiSqrInfo[1][i] = chiSqrProb(chiSqrInfo[0][i],
						chiSqrInfo[2][i], thresh);
			}
		}

		MatrixImpl mImpl = new MatrixImpl(chiSqrInfo);
		return mImpl.getMatrix();
	}

	/**
	 * Returns 2D array of confidence regions for each class ('group') and cv
	 * dimension
	 * 
	 * @param chiSqrProbs -
	 *            double array of prob vals for each dim
	 * @param thresh
	 * @return
	 * @throws NpairsException
	 */
	private double[][] confidReg(double[] chiSqrProbs, double thresh)
			throws NpairsException, MathException {

		double[][] confidReg = new double[nClasses][nCVDim];
		int numPAboveThresh = 0;
		for (double p : chiSqrProbs) {
			if (p >= thresh) {
				++numPAboveThresh;
			}
		}

		if (numPAboveThresh > 0) {
			double sumSqrdEv = 0;
			double sFact = 0;
			for (int i = 0; i < nClasses; ++i) {
				for (int j = 0; j < nCVDim; ++j) {
					// compute scaling factor for non-spherical variances
					for (int k = 0; k < cvaScores.numRows(); ++k) {
						sumSqrdEv += cvaScores.getQuick(k, j) * cvaScores.getQuick(k, j);
					}
					sFact = sumSqrdEv
							/ ((1 + cvaEvals[j]) * (nVol - nClasses));
					if (sFact == 0) {
						throw new NpairsException(
								"Error - cva scores cannot be all zeroes");
					}
					sFact = 1 / Math.sqrt(sFact);

					// compute radius for confidence circle
					confidReg[i][j] = (1 / sFact)
							* Math.sqrt((chiSqrCutoffVal(thresh,
									numPAboveThresh))
									/ clsSz[i]);
					sumSqrdEv = 0;
				}
			}
		}
		nSigDim = numPAboveThresh;
		
		return confidReg;
	}

	/**
	 * Returns probability of observing 'val' or something smaller from a
	 * chi-squared distribution with 'df' degrees of freedom: Probability(X <
	 * val) where X == chi-squared dist. with 'df' deg. of freedom
	 * 
	 * @param val
	 * @param df
	 * @param thresh TODO
	 * @return prob
	 */
	private static double chiSqrProb(double val, double df, double thresh) 
		throws MathException {
		
		ChiSquaredDistributionImpl chiSqrDist = new ChiSquaredDistributionImpl(df);
		double prob = -1;
		try {
			prob =  chiSqrDist.cumulativeProbability(val);
		}
		catch (MathException e) {
			double criticalVal = chiSqrDist.inverseCumulativeProbability(thresh);
//			System.out.println("Critical value for alpha = " + thresh + ": " + criticalVal);
//			System.out.println("Value to be evaluated: " + val);
			if (val > criticalVal) {
				// set probability of val > criticalVal to 1.0 to avoid convergence problems
				// in cumulativeProbability calculation of cdf for high vals.
				prob = 1.0;
			}
			else {
				throw new MathException("(CVA) Problem encountered calculating cumul. prob. of chi squared dist.: " +
					e.getMessage());					
			}
		}
		return prob;
	}

	/**
	 * Returns the cutoff value v such that
	 * 
	 * Probability(X < v) = a
	 * 
	 * where X is a random variable from the chi_sqr distribution with df degrees
	 * of freedom.
	 * 
	 * 
	 * @param prob
	 * @param df
	 * @return cutoff value 
	 */
	private static double chiSqrCutoffVal(double prob, double df)
			throws MathException {
		ChiSquaredDistributionImpl chiSqrDist = new ChiSquaredDistributionImpl(
				df);
		return chiSqrDist.inverseCumulativeProbability(prob);
	}

	/** Projects CV eigenvectors ('eigenimages') back into original space, 
	 * if initial feature selection step was performed before doing CVA.
	 * @param projFactorMat 
	 * <p>projection factor Matrix to be used to project data
	 *                      back into original space; see {@code 
	 *                      npairs.io.NpairsDataLoader.#getEVDProjFactorMat()}
	 *                      for details about projFactorMat.
	 *               <p> Set this Matrix to null if data is already in
	 *                   original space. 
	 * @param origData 
	 * <p>Data in original image space (before any data reduction or
	 *    other basis changes were performed) 
	 * <p>This argument is not read if projFactorMat set to null, i.e., if
	 *     data is already in original space.
	 * @see #getEigimsBig()
	 * @see #getFeatSpEigims()
	 * @see npairs.io.NpairsDataLoader.#getEVDProjFactorMat() 
	 */
	public void rotateEigimsToOrigSpace(Matrix projFactorMat, Matrix origData) {
		if (projFactorMat == null) {
			cvaEigimsBig = cvaFeatSpEigims;
			return;
		}
		// project cva results back onto projFactorMat (== inverted feature-selection 
		// projection Matrix for regular EVD; 
		// TODO: update documentation for unweighted EVD
		// (see PCA.rotateEigimsToOrigSpace(...) documentation for algebraic details)

		Matrix P1invP = cvaFeatSpEigims.transpose().mult(projFactorMat);

		Matrix voxSpaceCVAEigims = P1invP.mult(origData);

		cvaEigimsBig = voxSpaceCVAEigims.transpose();
	}
	
	/** Matches results in this CVA to reference CVA
	 * using Procrustes algorithm; see {@code npairs.utils.Procrustes(Matrix, Matrix)} for details 
	 * about Procrustes calculation and application of Procrustes to transform or negate/permute data.
	 
	 *
	 * @param refCVA 
	 * 			<p> CVA to use as reference. If useEigims = false, reference CVA and current CVA 
	 *              must have same # of CV dimensions and classes; if useEigims = true, reference 
	 *              CVA and current CVA must have same size and # of CV eigenimages.
	 * @param useEigims 
	 * 			<p> If true, use CV eigenimages to calculate 
	 * 	            Procrustes transformation.
	 *          <p> If false, use average CV scores for each class to calculate Procrustes transformation.
	 * @param useFullProcrustes 
	 *          <p> If true, apply Procrustes transformation [not currently implemented].
	 *          <p> If false, negate and permute CV dimensions instead. 
	 *          <p> Apply negation and permutation to:
	 *          <ul>
	 *          <li>CV eigenvalues
	 *          <li>CV eigenvectors (and all projected eigenimages)
	 *          <li>CV scores
	 *          <li>Chi-squared data
	 *          </ul>
	 *  
	 * <p>If both useEigims & useFullProcrustes == false:
	 *	Technique follows the one used in Jon Anderson's IDL NPAIRS code 'cva_reference.pro'.
	 * @see npairs.utils.Procrustes(Matrix, Matrix) 
	 * 
	*/
	
	public void matchToRef(CVA refCVA, boolean useEigims, boolean useFullProcrustes) 
		throws NpairsException {
		
		if (refCVA.getNumCVDims() < this.getNumCVDims()) {
			throw new NpairsException("Error - reference CVA must have at least " +
			"as many dimensions as current CVA");
		}

		double sTime = System.currentTimeMillis();
		Procrustes proc;
		if (!useEigims) {
			// Get average CV scores for each input condition
			Matrix avgRefCVScores = refCVA.avgCVScores();
			Matrix currAvgCVScores = this.avgCVScores();
			avgRefCVScores = avgRefCVScores.subMatrix(new int[] {0, avgRefCVScores.numRows() 
					- 1}, new int[] {0, currAvgCVScores.numCols() - 1});

			// apply Procrustes to the mean CV scores
			logStream.println("Calculating Procrustes using average CV scores... ");
			sTime = System.currentTimeMillis();
			proc = new Procrustes(avgRefCVScores.transpose(), 
					currAvgCVScores.transpose());
			double tTime = (System.currentTimeMillis() - sTime) / 1000;				
			logStream.println("[" + tTime + " s]");
		}
		else {
			// get CV eigenimages (in feat space if performed)
			Matrix refCVEigims = refCVA.getFeatSpEigims();
			Matrix currCVEigims = this.getFeatSpEigims();
			// apply Procrustes to eigenimages
			Npairs.getOutput().print("Calculating Procrustes using CV eigenimages.. ");
			sTime = System.currentTimeMillis();
			proc = new Procrustes(refCVEigims.transpose(), currCVEigims.transpose());
			double tTime = (System.currentTimeMillis() - sTime) / 1000;	
			Npairs.getOutput().println("[" + tTime + " s]");			
		}

		if (!useFullProcrustes) {
			int[] sign = proc.getRefSign();
			int[] permut = proc.getRefPermute();

			//			if (debug) {
			Npairs.getOutput().print("Signs: ");
			npairs.io.NpairsIO.print(sign, Npairs.getOutput());
			Npairs.getOutput().print("Permutations: ");
			npairs.io.NpairsIO.print(permut, Npairs.getOutput());
			//			}

			Npairs.getOutput().print("Negating and permuting dims... ");
			sTime = System.currentTimeMillis();
			this.negateAndPermuteDims(sign, permut);

			double ttTime = (System.currentTimeMillis() - sTime) / 1000;
			Npairs.getOutput().println("[" + ttTime + " s]");
		}
		else {
			// not fully implemented! 
			this.applyFullProcrust(proc);
		}
	}
	
	
	/**
	 * Saves CVA results to files. 
	 * <p>If original data has been transformed into new
	 * vector space (e.g., via initial feature selection and/or PCA), then this 
	 * method does NOT rotate data back into the original space first (must do 
	 * outside of this method).
	 * 
	 * <p>Output is saved in textfiles matching naming and content format used by 
	 * IDL NPAIRS. E.g., first row of files containing matrix data contains
	 * 2 integers corresponding to # columns and # rows in matrix; subsequent
	 * rows correspond to rows of matrix.
	 *   
	 * <p>Note that this method does NOT save R<sup>2</sup> data by default.  Must set
	 * saveR2 parameter to true for R<sup>2</sup> data to be saved here.
	 * 
	 * @param cvaSavePref
	 *            <p>Prefix (including path) of saved CVA files.
	 * @param saveAsSplit 
	 *            <p>Set to true if current CVA was run on split data
	 *                    (i.e., training subset of larger full data set); in
	 *                    this case, splitNum and splitHalf info will be 
	 *                    included in file suffix naming convention for all
	 *                    saved files. 
	 * @param splitNum 
	 *            <p>If saveAsSplit, which split number (in a given resampling 
	 *                 schema) corresponds to current split data.
	 * @param splitHalf 
	 *            <p>If saveAsSplit, which split half (in a given resampling
	 *                 schema) corresponds to current split data.
	 * @param saveR2 
	 *            <p>If true, save R<sup>2</sup> data for each CV dimension.
	 * @see {@link #npairs.Npairsj.saveR2() npairs.Npairsj.saveR2}
	 */
	// TODO: only save input number of CV dims.
	public void saveCVAResultsIDL(String cvaSavePref, boolean saveAsSplit,
			int splitNum, int splitHalf, boolean saveR2) throws IOException {
		String cvaEvalFile;
		String saveFormat = "IDL";
		String cvaEigimFile;
		String cvaScoreFile;
		String cvaRMSFile;
		String cvaClassInfoFile;
		String cvaChiSqrInfoFile;
		String cvaConfidRegFile;
		String cvaEvectFile;
		String cvaR2File = null;

		if (saveAsSplit) {
			cvaEvalFile = cvaSavePref + ".CVA." + splitNum + "." + splitHalf
					+ ".eigval";
			cvaEigimFile = cvaSavePref + ".CVA." + splitNum + "." + splitHalf
					+ ".eigim";
			cvaScoreFile = cvaSavePref + ".CVA." + splitNum + "." + splitHalf
					+ ".can";
			cvaRMSFile = cvaSavePref + ".CVA." + splitNum + "." + splitHalf
					+ ".rms";
			cvaClassInfoFile = cvaSavePref + ".CVA." + splitNum + "."
					+ splitHalf + ".group";
			cvaChiSqrInfoFile = cvaSavePref + ".CVA." + splitNum + "."
				+ splitHalf + ".chi";		
			cvaEvectFile = cvaSavePref + ".CVA." + splitNum + "." + splitHalf
					+ ".evect";
			cvaConfidRegFile = cvaSavePref + ".CVA." + splitNum + "."
					+ splitHalf + ".confidReg";
			if (saveR2) {
				cvaR2File = cvaSavePref + ".CVA." + splitNum + "."
					+ splitHalf + ".r2";
			}
		} else {
			cvaEvalFile = cvaSavePref + ".CVA.ALL.eigval";
			cvaEigimFile = cvaSavePref + ".CVA.ALL.eigim";
			cvaScoreFile = cvaSavePref + ".CVA.ALL.can";
			cvaRMSFile = cvaSavePref + ".CVA.ALL.rms";
			cvaClassInfoFile = cvaSavePref + ".CVA.ALL.group";
			cvaChiSqrInfoFile = cvaSavePref + ".CVA.ALL.chi";	
			cvaConfidRegFile = cvaSavePref + ".CVA.ALL.confidReg";
			cvaEvectFile = cvaSavePref + ".CVA.ALL.evect";
			if (saveR2) {
				cvaR2File = cvaSavePref + ".CVA.ALL.r2";
			}
		}

		NpairsIO.printToIDLFile(cvaEvals, cvaEvalFile);
		
		if (cvaEigimsBig != null) {
			cvaEigimsBig.printToFile(cvaEigimFile, saveFormat);
		}
		else {
			// save "feature space" eigenimages instead
			// (since these must actually be in orig space if no
			// initial selection was performed)
			cvaFeatSpEigims.printToFile(cvaEigimFile, saveFormat);
		}
		cvaScores.printToFile(cvaScoreFile, saveFormat);
		double[] rmsErrArray = { rmsError };
		NpairsIO.printToIDLFile(rmsErrArray, cvaRMSFile);
		NpairsIO.printRowToIDLFile(classLabels, cvaClassInfoFile);
		if (computeChiSqr) {
			chiSqrInfo.printToFile(cvaChiSqrInfoFile, "IDL");
		}
		// (note: confidreg info not normal idl npairs output)
		//MatrixImpl mImpl = new MatrixImpl(confidReg);
		//mImpl.getMatrix().printToFile(cvaConfidRegFile, "IDL");
		cvaEvectsSmall.printToFile(cvaEvectFile, saveFormat);
		if (saveR2) {
			r2.printToFile(cvaR2File, "IDL");
		}
	}

	 /** Returns number of CV dimensions in this analysis.
	 * @return # CV dimensions, i.e., the lesser of: 
	 *         <li> # input conditions (classes) - 1, or
	 *         <li> # dimensions (columns) in input data  
	 */
	public int getNumCVDims() {
		return nCVDim;
	}
	
	/** Calculates and returns average CV scores for each condition (class label)
	 *  specified for this CVA.
	 * @return Matrix containing average CV scores for each condition (class label)
	 *                    <li>Matrix size: # classes rows X # CV dimensions columns
	 *                    <li>rows correspond to class labels in 
	 *                      ascending order
	 * @see #getCVScores()
	 * @see pls.analysis.ResultSaver#avgCVScores(double[][], int[], int[], int[])
	 */
	public Matrix avgCVScores() {
			
		Hashtable<Integer, int[]> condIndices = MLFuncs.getLabelIndices(classLabels);
		int numCond = condIndices.size();
		int numCVDims = cvaScores.numCols();
		int[] sortedUniqCondLabels = new int[numCond];
		int i = 0;
		for (Enumeration<Integer> uniqCondLabels = condIndices.keys(); uniqCondLabels.hasMoreElements(); ) {
			sortedUniqCondLabels[i] = (Integer)uniqCondLabels.nextElement();
			++i;
		}
		sortedUniqCondLabels = MLFuncs.sortAscending(sortedUniqCondLabels);
		
		MatrixImpl avgScoresImpl = new MatrixImpl(numCond, numCVDims);
		Matrix avgScores = avgScoresImpl.getMatrix();
		for (int cond = 0; cond < numCond; ++cond) {
			int[] currCondIndices = condIndices.get(sortedUniqCondLabels[cond]);
			for (int dim = 0; dim < numCVDims; ++dim) {
				double[] currScores = MLFuncs.getItemsAtIndices(cvaScores.getColumn(dim),
						currCondIndices);
				avgScores.set(cond, dim, MLFuncs.avg(currScores));
			}
		}
		
		return avgScores;
	}
	
	/** Calculates and returns test CV scores, i.e.,
	 *  projection of test data onto this CVA's eigenimages. 
	 *  
	 * @param testData 
	 * 		<p>Matrix of data (rows = observations) for which
	 *        to calculate test CV scores using this CVA.
	 *        
	 *        <p>NOTE: If an initial feature selection (i.e., data-reduction) 
	 *        step such as initial eigenvalue decomposition
	 *        was performed on the data before calculating CVA, then 
	 *        testData must be in initial feature space as well, with
	 *        the same number of dimensions retained as in the CVA 
	 *        input data.  
	 *        
	 *        <p>This is because CV eigenimages used in this projection will 
	 *        then be in the initial feature space. (If PCA or some other
	 *        further projections were performed on data between initial 
	 *        feature selection and CVA, however, they should not also
	 *        be applied to test data, as CV eigenimages used in this
	 *        projection have already been projected back into either 
	 *        initial feature selection space, if relevant, or original
	 *        data space if no initial feature selection step was done.) 
	 *        
	 *           
	 *        	 
	 */
	public Matrix calcTestCVScores(Matrix testData) {
		return testData.mult(cvaFeatSpEigims);
	}

	/** Returns R<sup>2</sup> values calculated for each input data dimension and 
	 * CV dimension.
	 *  <p>R<sup>2</sup> is calculated as follows: 
	 *  <p>For each CV dimension, find correlation coefficients between CV scores and 
     *  'timecourse' for each input data dimension (e.g., PCA or voxel dimension). 
	 *  
	 * @return Matrix containing R<sup>2</sup> values. 
	 *         <p>Matrix size = # input data dimensions (e.g., # PCA dimensions) rows 
	 *                       X # CV dimensions columns.
	 * 
	 */
	public Matrix getR2() {
		return r2;
		
	}

	/** ******************************************************************************* */

	// Quick test of CVA.
	// (REQUIRED: input arg 'colt' or 'matlab' to indicate
	// matrix lib type)
//	public static void main(String[] args) {
//		String matlibType = "MATLAB";
//		MatrixImpl.setMatlibType(matlibType);
//		System.out.println("Matlib type: " + matlibType);
//
//		// using a 2-class, 2D dataset with 36 samples, 18 from each class
//		// (Mardia p. 329)
//		double[] x1 = { 191, 185, 200, 173, 171, 160, 188, 186, 174, 163, 190,
//				174, 201, 190, 182, 184, 177, 178, 186, 211, 201, 242, 184,
//				211, 217, 223, 208, 199, 211, 218, 203, 192, 195, 211, 187, 192 };
//		double[] x2 = { 131, 134, 137, 127, 118, 118, 134, 139, 131, 115, 143,
//				131, 130, 133, 130, 131, 127, 126, 107, 122, 114, 131, 108,
//				118, 122, 127, 125, 124, 129, 126, 122, 116, 123, 122, 123, 109 };
//		double[][] data = new double[2][36];
//		data[0] = x1;
//		data[1] = x2;
//		Matrix transM = new MatrixImpl(data).getMatrix();
//		Matrix M = null;
//		// JMatLink eng = null;
//		String idlFormatMatrixFname = "C:\\anita\\workspace\\PLSNPAIRSGoogleRepo"
//				+ "\\localTestData\\CVATest";
//
//		M = transM.transpose();
//
//		System.out.println("Input data: ");
//		M.print();
//		System.out.println("Saving matrix to file " + idlFormatMatrixFname);
//		M.printToFile(idlFormatMatrixFname, "IDL");
//		int[] labels = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
//				2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
//
//		try {
//			CVA cva = new CVA(M, labels);
//			System.out.println("W mat: ");
//			cva.getW().print();
//			System.out.println("B mat: ");
//			cva.getB().print();
//			System.out.println("W^(-1) * B:");
//			cva.getInvWTimesB().print();
//			System.out.println("labels: ");
//			for (int i = 0; i < labels.length; ++i) {
//				System.out.print(labels[i] + " ");
//
//			}
//			System.out.println();
//			System.out.println("Cva evals: ");
//			cva.getEvalMat().print();
//			System.out.println("Cva evects: ");
//			cva.getEvects().print();
//
//			System.out.println("RMS error: " + cva.getRMSError());
//
//			// String cvaSavePref = args[2];
//			// System.out.println("Saving CVA data to file " + cvaSavePref + "
//			// using saveCVAResults...");
//			//
//			// cva.saveCVAResultsIDL(cvaSavePref, false, 0, 1);
//		} catch (NpairsjException e) {
//			e.printStackTrace();
//		}
//		// catch (IOException e) {
//		// e.printStackTrace();
//		// }
//
//	}

	// testing chisqr/confidreg output:

//	private CVA(double[][] cvaScores, double[] cvaEvals, double[] classLabels) {
//		// set class variables required for chisqr/confidreg calculations
//		this.nVol = classLabels.length;
//		this.classLabels = new int[nVol];
//		for (int i = 0; i < nVol; ++i) {
//			this.classLabels[i] = (int) classLabels[i];
//		}
//
//		MatrixImpl mImpl = new MatrixImpl(cvaScores);
//		this.cvaScores = mImpl.getMatrix();
//		this.cvaEvals = cvaEvals;
//
//		this.nCVDim = this.cvaScores.numCols();
//		Hashtable<Integer, int[]> grpIndices = getLabelIndices(this.classLabels);
//		this.nGrp = grpIndices.size();
//		int[] sortedUniqGrpLabels = MLFuncs.sortAscending(MLFuncs
//				.unique(this.classLabels));
//		this.grpSz = new int[nGrp];
//		for (int g = 0; g < nGrp; ++g) {
//			grpSz[g] = grpIndices.get(sortedUniqGrpLabels[g]).length;
//		}
//	}

//	public static void main2(String[] args) {
//
//		MatrixImpl.setMatlibType("COLT");
//		String testNpairsPrefix = "C:\\anita\\workspace\\PLSwithNPAIRS\\localTestData\\testy2r4RasOrtho\\testy2r4RasOrtho";
//		try {
//			double[][] cvaScores = NpairsjIO.readFromIDLFile(testNpairsPrefix
//					+ ".CVA.ALL.can");
//			double[][] cvaEvals2D = NpairsjIO.readFromIDLFile(testNpairsPrefix
//					+ ".CVA.ALL.eigval");
//			double[][] classLabels2D = NpairsjIO
//					.readFromIDLFile(testNpairsPrefix + ".CVA.ALL.group");
//
//			double[] cvaEvals1D = cvaEvals2D[0];
//			double[] classLabels1D = classLabels2D[0];
//
//			// System.out.println("cvaEvals1D size: " + cvaEvals1D.length);
//			// System.out.println("classLabels1D size: " +
//			// classLabels1D.length);
//			System.out.println("Creating test CVA object...");
//			CVA cva = new CVA(cvaScores, cvaEvals1D, classLabels1D);
//			System.out.println("Testing chiSqrInfo...");
//			try {
//				Matrix chiSqrInfo = cva.chiSqrInfo(0.95);
//				System.out.println("Chi squared info: ");
//				chiSqrInfo.print();
//
//				double[][] confidReg = cva
//						.confidReg(chiSqrInfo.getRow(1), 0.95);
//				System.out.println("Confid regs: ");
//				NpairsjIO.print(confidReg);
//			} catch (MathException e) {
//				e.printStackTrace();
//			}
//
//			double val = 0.5;
//			double df = 2;
//
//			double chiSqrProb = chiSqrProb(val, df);
//			System.out.println("prob (X < " + val + "), df = " + df + ": "
//					+ chiSqrProb);
//
//			double chiSqrCutoffVal = chiSqrCutoffVal(chiSqrProb, df);
//			System.out.println("cutoff val s.t. prob(X < val) = " + chiSqrProb
//					+ "; df = " + df + ": " + chiSqrCutoffVal);
//		}
//
//		catch (MathException e) {
//			e.printStackTrace();
//		} catch (NpairsjException e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}

}

