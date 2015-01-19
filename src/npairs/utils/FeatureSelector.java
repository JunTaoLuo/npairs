package npairs.utils;

import java.io.IOException;

import npairs.Npairs;
import npairs.NpairsException;
import npairs.NpairsSetupParams;
import npairs.io.NpairsIO;
import npairs.shared.matlib.ColtMatrix;
import npairs.shared.matlib.EigenvalueDecomposition;
import npairs.shared.matlib.MatlabMatrix;
import npairs.shared.matlib.Matrix;
import npairs.shared.matlib.MatrixImpl;

public class FeatureSelector {
	
	boolean debug = false;
	NpairsSetupParams setupParams;
	Matrix origData;
	Matrix origSpaceXformFactor;
	String matlibType;
	Matrix featSelData;

	/** Constructs Matrix containing data transformed into reduced-dimension feature-selection format
	 * specified in input NpairsjSetupParams.
	 * @param setupParams contains parameters of given NPAIRS analysis, including type of initial
	 * 	feature selection and data reduction factor to be applied to data in feature selection space
	 * @param origData masked data in original (image) space; rows = scans, columns = voxels
	 * 				   </b>typically mean scans have been removed for each session in origData 
	 * 				   </b>(unless setupParams.doMSRBeforeInitEVD == false)
	 * @param matlibType type of Matrix library to use for feature selection
	 * @throws IOException
	 * @throws NpairsException
	 * @see #getFeatSelData()
	 */
	public FeatureSelector(NpairsSetupParams setupParams, Matrix origData, 
			String matlibType) throws IOException,
			NpairsException {
		this.setupParams = setupParams;
		this.origData = origData;
		this.matlibType = matlibType;
		featSelData = selectFeatures();
	}


	
	/** Returns 2D Matrix containing reduced-dimension feature-selected data using feature selection 
	 * technique specified in NpairsjSetupParams given to constructor. 
	 * 
	 */
	private Matrix selectFeatures() 
		throws IOException, NpairsException {
	
		Matrix featSelData = null;
		
		if (setupParams.doInitEVD()) {
			if (setupParams.loadEVD()) {
				Npairs.getOutput().print("Loading EVD from files... ");

				double sTime = System.currentTimeMillis();
				
				// load eigenvectors
				String svdEvectsFilename = setupParams.getEVDFilePref() + ".EVD.evects";
				Matrix svdEvects = null;
				double[] svdEvals1D = null;
				double[][] svdEvects2D = NpairsIO.readFromIDLFile(svdEvectsFilename);
				svdEvects = new MatrixImpl(svdEvects2D).getMatrix();
				
				// load eigenvalues
				String svdEvalsFilename = setupParams.getEVDFilePref() + ".EVD.evals";
				double[][] svdEvals2D = NpairsIO.readFromIDLFile(svdEvalsFilename);
				svdEvals1D = MatlabMatrix.trimTo1D(svdEvals2D);
				
				// calculate # of dims to retain using given data reduction factor
				double reducFactor = setupParams.getDataReductFactor();		
				int reducDataDims = 0;				
//				if (origData.numRows() < origData.numCols()) {
//					reducDataDims = (int)Math.round(origData.numRows() * reducFactor);
//				}
//				else {
//					reducDataDims = (int)Math.round(origData.numCols() * reducFactor);
////					throw new NpairsjException("WARNING: input data has more " +
////							"scans than masked voxels!");
//				}
				reducDataDims = (int)Math.round(svdEvects.numCols() * reducFactor);
			
				int[] rowRange = new int[] {0, svdEvects.numRows() - 1};
				int[] reducColRange = new int[] {0, reducDataDims - 1};
		
				// calculate reduced-dim eigenvectors 
				Matrix reducDimEvects = svdEvects.subMatrix(rowRange, reducColRange);
				
				// get singular values by finding sqrts of evd evals
				Matrix S = new MatrixImpl(reducDataDims, reducDataDims).getMatrix();
				for (int i = 0; i < reducDataDims; ++i) {
					S.setQuick(i,i, Math.sqrt(svdEvals1D[i]));
				}
				
				// set feature-selected data
				if (setupParams.normEVD()) { 
					// normed EVD ==> set variance of each eigenvector to 1 
					// so don't weight by eigenvalue
					featSelData = reducDimEvects;
				}
				else {
					featSelData = reducDimEvects.mult(S);	
				}

				// calculate Matrix origSpaceXformFactor - useful when transforming data 
				// back from feature-selection to original space.
				setOrigSpaceXformFactor(reducDimEvects, S);
				double tTime = (System.currentTimeMillis() - sTime) / 1000;
				
				Npairs.getOutput().println("[" + tTime + " s]");				
			}
			
			else {
				featSelData = selectFeaturesEVD(setupParams, origData, 
						setupParams.getDataReductFactor());
			}
		}
		
		return featSelData;
	}

	private void setOrigSpaceXformFactor(Matrix reducDimEvects, Matrix S) {
		int reducDataDims = reducDimEvects.numCols();
		Matrix invS = new MatrixImpl(reducDataDims, reducDataDims).getMatrix();
		for (int i = 0; i < reducDataDims; ++i) {
			invS.setQuick(i, i, 1 / S.getQuick(i, i));
		}

		origSpaceXformFactor = invS.mult(reducDimEvects.transpose());
	}

	
	/** Returns Matrix containing data represented in eigenvalue decomposition space.
	 * 
	 *  <p>Note that IDL uses EVD, too (actually eigenql of non-mc ssp mat) although it's
	 *  called 'SVD' in IDL NPAIRS. 
	 *  
	 *  <p>If using unweighted ("normed") EVD,
	then weighting matrix S is ignored in data representation P<sub>normed</sub>, i.e., 
	<ul>P<sub>normed</sub> = V  ( = MUS<sup>-1</sup> )  instead of </ul>
	<ul>P = VS  ( = MU ).   </ul>
 
<p>	This is equivalent to setting all (non-zero diagonal) elements of S to 1, 
	and  gets rid of differences in variance across basis vectors in U; hence 
	this can be thought of as a data "denoising" technique.
	In normed EVD, data is still represented in U-space, hence
	transformation back to original space is still a transformation back through U:
<ul>		P<sub>1</sub>U<sup>t</sup> = M<sub>1</sub> still holds.
</ul>	
Given P<sub>normed</sub> = V = MUS<sup>-1</sup>, U<sup>t</sup> = 
	S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup>M. 
	
<p>Hence we use S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup> in normed EVD where P<sup>-1</sup> is used
	in regular EVD. 
	
<p>	But note that 
<ul>	<li>P<sub>normed</sub> = V ==>  S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup> 
			=  S<sup>-1</sup>V<sup>t</sup>, 
and in regular EVD, 
	<li>P = VS ==> P<sup>-1</sup> = S<sup>-1</sup>V<sup>t</sup>. 
	</ul>	
<p>Therefore, we simply store  S<sup>-1</sup>V<sup>t</sup> in origSpaceXformFactor in both cases.
	 *   
	 * REQUIRED: 0.0 < dataReductionFactor <= 1.0
	 */
	private Matrix selectFeaturesEVD(NpairsSetupParams setupParams, Matrix data, 
			double dataReductionFactor) throws IOException, NpairsException {
		
		String taskMsg = "Running initial EVD with DRF = " + 
				dataReductionFactor + "... ";
		if (setupParams.normEVD()) {
			taskMsg = "Running initial EVD (normed) with DRF = " + 
					dataReductionFactor + "... ";
		}
		Npairs.getOutput().println(taskMsg);
		
		Matrix featSelData = null;
		double reducFactor = dataReductionFactor;		
		
		int reducDataDims = (int)Math.round(Math.min(data.numRows(), 
				data.numCols()) * reducFactor);
			
/*			 Given data = M, PCA(Mt) ==> MMt = V(S^2)Vt, where S^2 == diag Matrix containing
			 squared singular values from Matrix S in SVD(M) = VSUt along diagonal.
			 ==> featSelData = projection P of data M onto eigenimages U, 
			 i.e., featSelData P = MU
			                     = VSUtU 
			                     = VS, 
			 where U is data.numCols() X reducDataDims(rDD);
			      
			 		 S is rDD X rDD;
			 		 
			       V is data.numRows() X rDD;
			 			
			 ==> P = VS is data.numRows() X rDD
			       == dim(MU) as required */
			
			double sTime = System.currentTimeMillis();
			EigenvalueDecomposition evd = null;
			Matrix svdEvals = null;
			Matrix svdEvects = null;
			Matrix invSqrtEvals = null;
			if (!MatrixImpl.getMatlibType().equalsIgnoreCase(matlibType)) {
				//TODO: too many copies of data being made here!
				if (matlibType.equalsIgnoreCase("MATLAB")) {
					if (debug) {
						System.out.println("Doing EVD in Matlab...");
					}
					MatlabMatrix mlData = new MatlabMatrix(data.toArray());
					Npairs.getOutput().print("\tCreating SSP matrix from data matrix...");
					double sspSTime = System.currentTimeMillis();
					MatlabMatrix sspMat = mlData.mult(mlData.transpose());
					double sspTTime = (System.currentTimeMillis() - sspSTime) / 1000;
					Npairs.getOutput().println("[" + sspTTime + "s]");
					Npairs.getOutput().print("\tRunning EVD on SSP matrix...");
					double evdSTime = System.currentTimeMillis();
					evd = sspMat.eigenvalueDecomposition();
					double evdTTime = (System.currentTimeMillis() - evdSTime) / 1000;
					Npairs.getOutput().println("[" + evdTTime + "s]");
				}
				if (matlibType.equalsIgnoreCase("COLT")) {
					if (debug) {
						System.out.println("Doing EVD in Colt...");
					}
					ColtMatrix coltData = new ColtMatrix(data.toArray());
					Npairs.getOutput().print("\tCreating SSP matrix from data matrix...");
					double sspSTime = System.currentTimeMillis();
					ColtMatrix sspMat = coltData.mult(coltData.transpose());
					double sspTTime = (System.currentTimeMillis() - sspSTime) / 1000;
					Npairs.getOutput().println("[" + sspTTime + "s]");
					Npairs.getOutput().print("\tRunning EVD on SSP matrix...");
					double evdSTime = System.currentTimeMillis();
					evd = sspMat.eigenvalueDecomposition();
					double evdTTime = (System.currentTimeMillis() - evdSTime) / 1000;
					Npairs.getOutput().println("[" + evdTTime + "s]");
				}
				
				svdEvals = new MatrixImpl(evd.getRealEvalMat().toArray()).getMatrix();
				svdEvects = new MatrixImpl(evd.getEvects().toArray()).getMatrix();
				invSqrtEvals = new MatrixImpl(evd.getInvSqrtRealEvalMat().toArray()).getMatrix();
			}
			else {
				Npairs.getOutput().print("\tCreating SSP matrix from data matrix...");
				double sspSTime = System.currentTimeMillis();
				Matrix sspData = data.sspByRow();
				double sspTTime = (System.currentTimeMillis() - sspSTime) / 1000;
				Npairs.getOutput().println("[" + sspTTime + "s]");
				Npairs.getOutput().print("\tRunning EVD on SSP matrix...");
				double evdSTime = System.currentTimeMillis();
				evd = sspData.eigenvalueDecomposition();
				double evdTTime = (System.currentTimeMillis() - evdSTime) / 1000;
				Npairs.getOutput().println("[" + evdTTime + "s]");
				svdEvals = evd.getRealEvalMat();
				svdEvects = evd.getEvects();
				invSqrtEvals = evd.getInvSqrtRealEvalMat();
			}

//			if (debug) {
				double tTime = (System.currentTimeMillis() - sTime) / 1000;
				int hr = (int)(tTime / 3600);
				int min = (int)(tTime / 60) - (hr * 60);
				double s = tTime - (hr * 3600) - (min * 60) ;
				Npairs.getOutput().print("Total time doing EVD calculations: [" + hr + " h " +
						min + " min ");
				Npairs.getOutput().printf("%.3f", s);
				Npairs.getOutput().println(" s]");
//			}
			
			// save evals/evects, i.e. S^2 and V (TODO: save only if user specifies to do so?):
			String evalsFilename = setupParams.getResultsFilePrefix() + ".EVD.evals";
			String evectsFilename = setupParams.getResultsFilePrefix() + ".EVD.evects";
			double saveSTime = System.currentTimeMillis();
			Npairs.getOutput().print("Saving EVD info to file...");
			NpairsIO.printToIDLFile(evd.getRealEvals(), evalsFilename);
			svdEvects.printToFile(evectsFilename, "IDL");
			double saveTTime = (System.currentTimeMillis() - saveSTime) / 1000;
			Npairs.getOutput().println("[" + saveTTime + " s]");
			
			//***********************************************************
			// Save eigims too!
			//***********************************************************
//			System.out.print("Creating and saving svd eigims...");
//			double strtTime = System.currentTimeMillis();
//			String eigimsFilename = setupParams.resultsFilePrefix + ".EVD.eigims";
//			Matrix svdEigims = data.transpose().mult(svdEvects).mult(invSqrtEvals);
//			svdEigims.printToFile(eigimsFilename, "IDL");
//			double totTime = (System.currentTimeMillis() - strtTime) / 1000;
//			System.out.println("[" + totTime + " s]");
			int[] rowRange = new int[] {0, data.numRows() - 1};
			int[] reducColRange = new int[] {0, reducDataDims - 1};
			
			if (setupParams.normEVD()) {
				// ignore S so featSelData = V
				featSelData = svdEvects.subMatrix(rowRange, reducColRange);
			}
			else {
				// Given data M = VSUt, consider basis space to be
				// U, hence MU = VS, i.e., featSelData = VS
				Matrix S = svdEvals.subMatrix(reducColRange, reducColRange);
				for (int i = 0; i < reducDataDims; ++i) {
					double currEval = S.get(i, i);
					if (currEval > 0) {
						S.set(i, i, Math.sqrt(S.get(i, i)));
					}
					else S.set(i, i, 0);
				}

				featSelData = svdEvects.subMatrix(rowRange, reducColRange).
					mult(S);
			}
			
			double iSTime = System.currentTimeMillis();
			Npairs.getOutput().print("Calculating inverse EVD matrix ...");
			
			// see setOrigSpaceXformFactor()
			origSpaceXformFactor = invSqrtEvals.subMatrix(reducColRange, reducColRange).
					mult(svdEvects.subMatrix(rowRange, reducColRange).transpose());
			double iTTime = (System.currentTimeMillis() - iSTime) / 1000;
			Npairs.getOutput().println("[" + iTTime + "s]");
			        					
//		}
		
		return featSelData;
	}

/**	Returns a Matrix that is useful when transforming data back from dimension-reduced 
 * (e.g., EVD) space to original data (e.g., voxel) space.
 * 
 * <p> Let this matrix be F. 
 * 
 * <p> What is F?
 * 
 * <p> Assume data-reduction step is done via EVD, as this is what is currently implemented. Given original data M, 
 * where M is <i>m</i> observations X <i>n</i> variables (voxels), let M = VSU<sup>t</sup> be a singular
 * value decomposition (SVD) representation of M. 

<h4>Regular EVD</h4>
Let P = MU. Then P is M transformed into U-space. P usually represents M data
	in a much smaller matrix, since P is only <i>m</i> X <i>d</i>, where <i>d</i> = min(<i>m</i>, <i>n</i>),
	and usually <i>n</i> >> <i>m</i>.
	
<p> Data size may be further reduced by retaining only a subset of dimensions in P. This data reduction step
	may also be called 'feature selection', since EVD effectively separates the data into dimensions containing
	'features' composed of consistent patterns of variability, and the user may choose to exclude dimensions which
	reflect very little data variability or which consist of identifiably confounding features. 
	
<p> Then P* = MU* is a data-reduced representation of M, where U* is <i>n</i> X <i>k</i> and 
	contains only <i>k</i> of <i>d</i> dimensions.
	
<p> Hence M* = P*U*<sup>t</sup> (where M* = M exactly only if <i>k</i> = <i>d</i>), and to transform any other 
data  P<sub>1</sub> from U*-space to original space:

<p>M<sub>1</sub> = P<sub>1</sub>U*<sup>t</sup>

<p> But even with data reduction, U*<sup>t</sup> is very large to store. (Size of U* <= size of M.) 

<p> Hence we save F = P*<sup>-1</sup> = S*<sup>-1</sup>V*<sup>t</sup> instead, where S*, V* are the 
<i>k</i>-dimensional submatrices of S and V and dimensions of F are only <i>m</i> X <i>k</i>. 

<p> Then we could reconstruct U*<sup>t</sup> on the fly: U*<sup>t</sup> = FM. But we won't bother reconstructing
U*<sup>t</sup>. Instead we will calculate P<sub>1</sub>P*<sup>-1</sup>
before multiplying by M, since:

<p> M<sub>1</sub> = P<sub>1</sub>U*<sup>t</sup> = P<sub>1</sub>FM = P<sub>1</sub>P*<sup>-1</sup>M
<p> By calculating P<sub>1</sub>P*<sup>-1</sup> first, we only
   do matrix multiplication with huge dimension (<i>n</i>)
   once!
   
<h4>Normed EVD</h4> 
F = S*<sup>-1</sup>V*<sup>t</sup> for normed EVD as well. For regular EVD, P* = MU* = V*S*. 
But if normed EVD is used for feature selection instead of regular EVD, coefficients of data in U*-space are not 
weighted by S*:

	<ul>P* = V* = MU*S*<sup>-1</sup> </ul>
	<ul>==> U*<sup>t</sup> = S*<sup>-1</sup>P*<sup>-1</sup>M</ul>
	
	Therefore to transform any new data P<sub>1</sub> from U* to original space:
	<ul> M<sub>1</sub> = P<sub>1</sub>U*<sup>t</sup> = P<sub>1</sub>S*<sup>-1</sup>P*<sup>-1</sup>M</ul>
	
<p> Hence F = S*<sup>-1</sup>P*<sup>-1</sup> = S*<sup>-1</sup>V*<sup>t</sup>, the same as for Regular EVD. 

<p> F is <i>k</i> X <i>m</i>.

@return <p> Matrix containing S*<sup>-1</sup>V*<sup>t</sup> from data-reduced subset V*S*U*<sup>t</sup> of
	singular value decomposition M =  VSU<sup>t</sup>, given original data M.
			
	<p> Useful when transforming matrices calculated on transformed
			data back into original space (instead of calculating transform U*<sup>t</sup> every 
			time on the fly):
			<ul>U*<sup>t</sup> = S*<sup>-1</sup>V*<sup>t</sup>M
			<p>==></ul>
			<ul>P<sub>1</sub>U*<sup>t</sup> = P<sub>1</sub>(S*<sup>-1</sup>V*<sup>t</sup>)M
			</ul>
@see #getFeatSelData()
*/
public Matrix getOrigSpaceXformFactorMatrix() {
	//(NOTE  origSpaceXformFactor == ipcmat in IDL NPAIRS).
	return origSpaceXformFactor;
}


/** Returns Matrix containing data represented in feature selection space.
 * 
 *  Given original data M, let EVD(M<sup>t</sup>M) = US<sup>2</sup>U<sup>t</sup>.
 *  
 *  <p>If using unweighted ("normed") EVD: weighting matrix S is ignored in data 
 *  representation P<sub>normed</sub>, i.e., 
<ul>P<sub>normed</sub> = V  ( = MUS<sup>-1</sup> )  instead of </ul>
<ul>P = VS  ( = MU ).   </ul>

<p>	This is equivalent to setting all (non-zero diagonal) elements of S to 1, 
and  gets rid of differences in variance across basis vectors in U; hence 
this can be thought of as a data "denoising" technique.
In normed EVD, data is still represented in U-space, hence
transformation back to original space is still a transformation back through U:
<ul>		P<sub>1</sub>U<sup>t</sup> = M<sub>1</sub> still holds.
</ul>	
Given P<sub>normed</sub> = V = MUS<sup>-1</sup>, U<sup>t</sup> = 
S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup>M. 

<p>Hence we use S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup> in normed EVD where P<sup>-1</sup> is used
in regular EVD. 

<p>	But note that 
<ul>	<li>P<sub>normed</sub> = V ==>  S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup> 
		=  S<sup>-1</sup>V<sup>t</sup>, 
and in regular EVD, 
<li>P = VS ==> P<sup>-1</sup> = S<sup>-1</sup>V<sup>t</sup>. 
</ul>	
<p>Therefore, we simply store  S<sup>-1</sup>V<sup>t</sup> in origSpaceXformFactor in both cases.
*/
public Matrix getFeatSelData() {
	return featSelData;
}

}
