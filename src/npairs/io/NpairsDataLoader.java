package npairs.io;

import npairs.shared.matlib.*;
import npairs.utils.FeatureSelector;
import extern.niftijlib.Nifti1Dataset;
import npairs.Npairs;
import npairs.NpairsSetupParams;
import npairs.NpairsException;

import com.jmatio.types.MLDouble;

import extern.NewMatFileReader;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;

import pls.shared.MLFuncs;

/**<p> Loads data volumes listed in input NpairsjSetupParams file or object into 2D Matrix 
 * of given type. Performs feature selection (e.g.
 *  eigenvalue decomposition (EVD)) if NpairsjSetupParams includes an initial feature selection step
 *  and stores reduced-dimension feature-selected data in 2D Matrix.
 *  
 * <p> Also calculates and stores Matrix useful when
 *  transforming data back from feature-selection to original space.
 * 
 * <p> The currently implemented feature selection options are EVD and normed EVD.
 *  */

public class NpairsDataLoader {
	
	final boolean debug = true;
	
	private Matrix dataInOrigSpace;
	
//	private Matrix featSelData = null;
	
//	private Matrix origSpaceXformFactor;
	
	private String matlibType;
	private String matlibTypeForInitFeatSel;
	
	private float[] qOffset; 
	private double[] voxelSize;
	private int[] origin; // in 1-relative voxels a la Matlab PLS
	private int[] volDims3D; 
	
	private int[] maskCoords;
	
	private boolean loadDatamats = false;
	
	private FeatureSelector featSelector;
	
	private NpairsSetupParams setupParams;
	
	/** Loads NPAIRS masked data volumes into 2D Matrix of given type. Each row contains 
	 * a 3D volume in 1D format. Assuming 3D data is in [Z][Y][X] format, 1D coordinate 
	 * ordering is [x + y*xdim + z*ydim*xdim].
	 * 
	 * <p> Also transforms data into feature-selection space if initial feature selection
	 * option has been specified in setup parameters. Original masked data Matrix is
	 * retained. Only the number of dimensions indicated in setup parameters are 
	 * retained in feature-selected data - e.g., if data reduction factor = 0.3, then
	 * 30% of data dimensions are retained in feature selected data. 
	 * 
	 * @param npairsjSetupParamsMatFilename - name of NPAIRS setup parameters file
	 * 								(including full path and file extension)
	 * @param matlibType - which Matrix library to use (e.g. ParallelColt, Matlab)
	 * @param matlibTypeForInitFeatSel - which Matrix library to use for initial feature
	 * 					selection step (usually same as matlibType)
	 * @param loadDatamats - if true, load PLS-style datamat data instead of original image 
	 *                       volumes
	 *                     - datamats are 'unwound' so each data volume corresponding to a lag 
	 *                       within the given time window is contained in an individual row
	 *                       of the data Matrix - e.g., if window size = 8 then for each row
	 *                       of the PLS-style datamat there will be 8 consecutive rows in the 
	 *                       NPAIRS Matrix
	 * @throws NpairsException
	 * @throws IOException - if setup parameters file cannot be read
	 * @see NpairsDataLoader(NpairsjSetupParams, String, String, boolean)
	 * @see #getDataInOrigSpace()
	 * @see #getFeatSelData()
	 */
	public NpairsDataLoader(String npairsjSetupParamsMatFilename, 
			String matlibType, String matlibTypeForInitFeatSel, boolean loadDatamats) 
			throws NpairsException, IOException {
		
		this(new NpairsSetupParams(npairsjSetupParamsMatFilename, false),
				matlibType, matlibTypeForInitFeatSel, loadDatamats);
	}
	
	/** Loads NPAIRS masked data volumes into 2D Matrix of given type. Each row contains 
	 * a 3D volume in 1D format. Assuming 3D data is in [Z][Y][X] format, 1D coordinate 
	 * ordering is [x + y*xdim + z*ydim*xdim].
	 * 
	 * <p> Also transforms data into feature-selection space if initial feature selection
	 * option has been specified in setup parameters. Original masked data Matrix is still
	 * retained. Only the number of dimensions indicated in setup parameters are 
	 * retained in feature-selected data - e.g., if data reduction factor = 0.3, then
	 * 30% of data dimensions are retained in feature selected data. 
	 * 
	 * @param nsp - object containing setup parameters for current NPAIRS analysis
	 * @param matlibType - which Matrix library to use (e.g. ParallelColt, Matlab)
	 * @param matlibTypeForInitFeatSel - which Matrix library to use for initial feature
	 * 					selection step (usually same as matlibType)
	 * @param loadDatamats - if true, load PLS-style datamat data instead of original image 
	 *                       volumes
	 *                     - datamats are 'unwound' so each data volume corresponding to a lag 
	 *                       within the given time window is contained in an individual row
	 *                       of the data Matrix - e.g., if window size = 8 then for each row
	 *                       of the PLS-style datamat there will be 8 consecutive rows in the 
	 *                       NPAIRS Matrix
	 * @throws NpairsException
	 * @throws IOException - if setup parameters file cannot be read
	 * @see NpairsDataLoader(String, String, String, boolean)
	 * @see #getDataInOrigSpace()
	 * @see #getFeatSelData()
	 */
	public NpairsDataLoader(NpairsSetupParams nsp, String matlibType,
			String matlibTypeForInitFeatSel, boolean loadDatamats) throws NpairsException,
		IOException {
			this.matlibType = matlibType;
			this.matlibTypeForInitFeatSel = matlibTypeForInitFeatSel;
			this.loadDatamats = loadDatamats;
			this.setupParams = nsp;
			loadNpairsData();
	}
	
	/** Returns 2D Matrix containing masked data
	 * 
	 * @return 2D Matrix containing masked data
	 * @see NpairsDataLoader(String, String, String, boolean)
	 * @see #getFeatSelData()
	 */
	public Matrix getDataInOrigSpace() {
		return dataInOrigSpace;
	}

	/** Returns 2D Matrix containing masked reduced-dimension (i.e.,
	 * feature-selected) data. Only the number of dimensions indicated in 
	 * setup parameters are returned - e.g., if data reduction factor = 0.3, 
	 * then 30% of data dimensions are returned. Dimensions are ordered by
	 * amount of variability accounted for and retained dimensions account
	 * for as much variability as possible.
	 * 
	 * @return 2D Matrix containing feature-selected data
	 *		<p> Returns null if no feature selection was performed. 
	 * @see NpairsDataLoader(String, String, String, boolean)
	 * @see #getDataInOrigSpace()
	 * @see #getOrigSpaceXFormFactor()
	 */
	public Matrix getFeatSelData() {
		return featSelector.getFeatSelData();
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
		return featSelector.getOrigSpaceXformFactorMatrix();
	}
	
	/** 
	 * 
	 * @return True if PLS-style datamats are loaded instead of original data volumes; 
	 *         false if image volumes are loaded 
	 */
	public boolean loadDatamats() {
		return loadDatamats;
	}
	
	
	private void loadNpairsData() 
			throws NpairsException, IOException {
		

		// Load data (each row currently in IDL format: [Y][X][Z]):		
		try {
			double sTime = System.currentTimeMillis();
			double tTime = 0;
			if (!loadDatamats) {
				Npairs.getOutput().print("Loading data from image files...");
				dataInOrigSpace = loadOrigData(setupParams);
				tTime = (System.currentTimeMillis() - sTime) / 1000;
				setDataSpecs(setupParams.getDataFilenames()[0]);
			}
			else {
				Npairs.getOutput().print("Loading datamats...");
				dataInOrigSpace = loadDatamats(setupParams, matlibType);
				tTime = (System.currentTimeMillis() - sTime) / 1000;

			}
			Npairs.getOutput().println("[" + tTime + " s]");
			Npairs.getOutput().println("No. vols: " + dataInOrigSpace.numRows());
		}

		catch (OutOfMemoryError e) {
			throw new NpairsException("Ran out of memory loading data.  Try allocating more " +
					"memory when running plsnpairs.");
		}
		catch (MatrixException e) {
			throw new NpairsException("Error loading data volumes -  " + e.getMessage());
		}
	
		// Preprocessing (currently only mean session scan removal)
		if (setupParams.preProcess() && setupParams.doMSRBeforeInitEVD()) {
			preprocessData();
		}
		
		// Feature Selection:

		Npairs.getOutput().print("Do initial feature selection? ");
		if (setupParams.doInitFeatSelect()) {
			Npairs.getOutput().println("Yes");
		}
		else {
			Npairs.getOutput().println("No");
		}

		
		if (setupParams.doInitFeatSelect()) {
			featSelector = new FeatureSelector(setupParams, dataInOrigSpace, matlibTypeForInitFeatSel);
			//featSelData = featSelector.selectFeatures();
		}
	}
	
	private void preprocessData() throws NpairsException {
		// do MSR (mean session scan removal)? 
		if (setupParams.removeSessionMeans()) {
			long start = System.currentTimeMillis();
			if (debug) {
				System.out.println("Doing MSR before feature selection...");
			}
			Npairs.getOutput().print("Doing MSR before feature selection...");
			dataInOrigSpace = removeMeanSessionScans();
			long end = System.currentTimeMillis();
			long time = (end - start) / 1000;
			if (debug) {
				System.out.println("Done MSR [" + time + " s");
			}
			
			Npairs.getOutput().print("Done [" + time + " s]\n");
		}
		// more to be implemented
	}
	
	/** Remove mean scan from each session to get rid of session effects. Algorithm used is same
	 * as what is used in idl NPAIRS ssm_xform.pro (see 'transf' option 8: 'subtract out subject mean 
	 * profiles') 
	 */
	private Matrix removeMeanSessionScans() 
			throws NpairsException {
		int[] uniqSessLabs = MLFuncs.unique(setupParams.getSessLabels());
		int nSess = uniqSessLabs.length;
		int nDataDims = dataInOrigSpace.numCols();
		
		if (nSess > 1) {
			double[] grandMean = new double[nDataDims];
			// mean session scans are held in columns of sessMeans, not rows
			Matrix sessMeans = new MatrixImpl(nDataDims, nSess).getMatrix();
			int[] grandCount = new int[nDataDims];
		
			for (int s = 0; s < nSess; ++s) {
				int[] tmpCount = new int[nDataDims];
				int[] currScanLocs = MLFuncs.find(setupParams.getSessLabels(), uniqSessLabs[s]);
				int nCurrScans = currScanLocs.length;
				for (int i = 0; i < nCurrScans; ++i) {
					// add current scan data to running total for curr session in sessMeans
					double[] currScanData = dataInOrigSpace.getRow(currScanLocs[i]);
					int[] whereNZ = MLFuncs.findNonZero(currScanData);
				
					double[] currSessSumScans = sessMeans.getColumn(s);
					int nNZ = whereNZ.length;
					for (int n = 0; n < nNZ; ++n) {
						currSessSumScans[whereNZ[n]] += currScanData[whereNZ[n]];
						// increment element counter where curr scan has non-zero values	
						tmpCount[whereNZ[n]] += 1;
					}
					// set sessMeans to cumulative scan total for this session (for now)
					sessMeans.setColumn(s, currSessSumScans);
								
				}
				int[] currSessNZ = MLFuncs.findNonZero(tmpCount);
				if (currSessNZ.length > 0) {
					
					double[] currSessScanSum = sessMeans.getColumn(s);
					double[] currSessAvgScan = new double[nDataDims];
					for (int nzLoc : currSessNZ) {
						currSessAvgScan[nzLoc] = currSessScanSum[nzLoc] / tmpCount[nzLoc];	
						grandMean[nzLoc] += currSessScanSum[nzLoc];
						grandCount[nzLoc] += tmpCount[nzLoc];
					}
					sessMeans.setColumn(s, currSessAvgScan);		
				}
			}
			// grandMean == avg scan (actually avg of all non-zero values at each voxel 
			// or equivalent data basis component)
			int[] whereNZ = MLFuncs.findNonZero(grandCount);
			for (int nzLoc : whereNZ) {
				grandMean[nzLoc] /= grandCount[nzLoc];
			}
			
			// Remove session means from data (leaving in grand mean).  This
			// will not change the grand mean; it just means that each session 
			// mean will equal the grand mean.
			for (int s = 0; s < nSess; ++s) {
				// diffs is grandmean - sessmean, so e.g. if grandmean is 0, then 
				// diffs is -sessmean, hence *adding* diffs to curr scan data results
				// in removing current sessmean 
				double[] diffs = MLFuncs.subtract(grandMean, sessMeans.getColumn(s));
				int[] currScanLocs = MLFuncs.find(setupParams.getSessLabels(), uniqSessLabs[s]);
				int nCurrScans = currScanLocs.length;
				for (int i = 0; i < nCurrScans; ++i) {
					double[] currScanData = dataInOrigSpace.getRow(currScanLocs[i]);
					for (int j = 0; j < nDataDims; ++j) {
						if (currScanData[j] != 0) {
							currScanData[j] += diffs[j];
						}
					}
					dataInOrigSpace.setRow(currScanLocs[i], currScanData);
				}
			}
			
		} // end else
		
		// print loc of zeros for each data dim
//		System.out.println("Num scans: " + data.numRows());
//		for (int j = 0; j < nDataDims; ++j) {
//			int[] whereNZ = MLFuncs.findNonZero(data.getColumn(j));
//			int numZ = data.numRows() - whereNZ.length;
//			System.out.println("No. zeros (dim " + j + "): " + numZ);
//		}
		
//		// fill in zeros by regressing global mean against regional ('voxel') values
//		int nScans = data.numRows();
//		double[] scanMeans = new double[data.numRows()];
//		for (int i = 0; i < nScans; ++i) {
//			scanMeans[i] = MLFuncs.avg(data.getRow(i));
//		}
//		for (int j = 0; j < nDataDims; ++j) {
//			int[] whereNZ = MLFuncs.findNonZero(data.getColumn(j));
//			if (whereNZ.length > 0 && whereNZ.length < nScans) {	
//				// some but not all scans have 0 in current dim;
//				// extract the non-zero scan elements in curr dim and the
//				// corresponding scan means
//				double[] meansNZ = MLFuncs.getItemsAtIndices(scanMeans, whereNZ);
//				double[] dataNZ = data.subMatrix(whereNZ, new int[] {j}).getColumn(0);
//			}
//		}
		
		// now remove grand mean (i.e., column-centre the data).  * ASSUMES NO ZEROS *
		// TODO: consider whether there ever could be zeros in masked data; if not, don't
		// need non-zero-finding techniques above, either.
		// (NOTE: idl ssm_xform fills in zeros by regressing global mean against regional
		// ('voxel') values)
		dataInOrigSpace = dataInOrigSpace.meanCentreColumns();
		
		return dataInOrigSpace;
	}
	
	
	/** Sets voxel size, origin, qoffset and data dims using given dataset 
	 * 
	 * @param dataFilename - name of existing nifti data file from which to read header info
	 */
	private void setDataSpecs(String dataFilename) throws FileNotFoundException, IOException {

		Nifti1Dataset nDS = new Nifti1Dataset(dataFilename);
		nDS.readHeader();
		volDims3D = new int[] {(int)nDS.getXdim(), (int)nDS.getYdim(), (int)nDS.getZdim()};
		float[] voxSizeTmp = nDS.pixdim;
		// TODO: is this order of voxelSize always valid?
		voxelSize = new double[] {voxSizeTmp[1], voxSizeTmp[2], voxSizeTmp[3]}; 
		
		qOffset = nDS.qoffset;
		// take qoffset and calculate origin in voxels
		//  - for each dim, origin = -(qoffset) / (voxsize)
		origin = new int[3];
		for (int i = 0; i < 3; ++i) {
			origin[i] = (int) (-qOffset[i] / voxelSize[i]) + 1; // add one because 1-rel		                                                         // voxels are 1-rel.
		}	
	}
	
	private static Matrix loadDatamats(NpairsSetupParams setupParams, String matlibType) throws IOException,
		NpairsException, MatrixException {
		
		String[] datamatFilenames = setupParams.getDatamatFilenames();
		
		// see also ConcatenateDatamat.computeCommonCoordinates()
		int[] andMask = computeCommonCoordinates(datamatFilenames);
		int nRowsAll = 0;
		for (String d : datamatFilenames) {
			NewMatFileReader dmatReader = new NewMatFileReader(d);
			int currNConds = ((MLDouble)dmatReader.getMLArray("st_evt_list")).getN();
			int currWinSz = ((MLDouble)dmatReader.getMLArray("st_win_size")).get(0, 0).intValue();
			nRowsAll += currNConds * currWinSz;
		}
		
		int[] skipTmpts = setupParams.getSkipTmpts(); // which rows of concat. datamat to 
		                                         // exclude
		int nRows = nRowsAll - skipTmpts.length;
		
		Matrix allData = new MatrixImpl(nRows, andMask.length, matlibType).getMatrix();
		
		int beginRow = 0;
		for (String d : datamatFilenames) {
			NpairsReadDatamat nrd = new NpairsReadDatamat(d, setupParams.getCondSelection());
			Matrix currDatamat = nrd.getDatamat(andMask);
			allData.setSubMatrix(currDatamat, beginRow, 0);
			beginRow += currDatamat.numRows();
		}
		
		return allData;	
	}
	
	// see also ConcatenateDatamat.computeCommonCoordinates()
	private static int[] computeCommonCoordinates(String[] datamatFilenames) throws 
			IOException, NpairsException {
		int[] dims = ((MLDouble)new NewMatFileReader(datamatFilenames[0]).
				getMLArray("st_dims")).getIntFirstRowOfArray();
		int[] idxCount = new int[dims[0] * dims[1] * dims[2] * dims[3]]; 
		for (String dmfile : datamatFilenames) {
			// check dims
			NewMatFileReader fileReader = new NewMatFileReader(dmfile);
			int[] currDims = ((MLDouble)fileReader.
					getMLArray("st_dims")).getIntFirstRowOfArray();
			if (!Arrays.equals(dims, currDims)) {
				throw new NpairsException("All input data must be of same resolution.");
			}
			// get coords
			int[] stCoords = ((MLDouble)fileReader.getMLArray("st_coords")).
				getIntFirstRowOfArray();
			for (int c : stCoords) {
				++idxCount[c];
			}
		}
		
		int[] andMask = MLFuncs.find(idxCount, datamatFilenames.length);
		return andMask;
	}
	
	/** 
	 * 
	 * @param setupParams contains NPAIRS setup parameters
	 * @return Matrix containing masked input data
	 * 			<p>- rows = masked 3d volumes 
	 * 			<p>- assuming 3D data is in [Z][Y][X] format, 1D coordinate
 * 						ordering is [x + y*xdim + z*ydim*xdim]
	 * @throws MatrixException
	 * @throws IOException if files cannot be loaded
	 */
	private Matrix loadOrigData(NpairsSetupParams setupParams) throws MatrixException, IOException {
		double[] andMask = NiftiIO.getANDMask(setupParams.getMaskFilenames());
		maskCoords = MLFuncs.findNonZero(andMask);				
		Matrix maskedData =  NiftiIO.getMaskedDataMat(setupParams, matlibType, maskCoords);
	
		return maskedData;
	}
	
	
//	/** Returns Matrix containing feature-selected data using feature selection 
//	 * technique selected in setupParams.  
//	 * 
//	 * REQUIRED: 0.0 < nsp.dataReductionFactor <= 1.0
//	 */
//	private Matrix selectFeatures(NpairsjSetupParams setupParams) 
//		throws IOException, NpairsjException {
//	
//		Matrix featSelData = null;
//		
//		if (setupParams.doInitEVD()) {
//			if (setupParams.loadEVD()) {
//				Npairsj.output.print("Loading EVD from files... ");
//
//				double sTime = System.currentTimeMillis();
//				
//				// load eigenvectors
//				String svdEvectsFilename = setupParams.getEVDFilePref() + ".EVD.evects";
//				Matrix svdEvects = null;
//				double[] svdEvals1D = null;
//				double[][] svdEvects2D = NpairsjIO.readFromIDLFile(svdEvectsFilename);
//				svdEvects = new MatrixImpl(svdEvects2D).getMatrix();
//				
//				// load eigenvalues
//				String svdEvalsFilename = setupParams.getEVDFilePref() + ".EVD.evals";
//				double[][] svdEvals2D = NpairsjIO.readFromIDLFile(svdEvalsFilename);
//				svdEvals1D = MatlabMatrix.trimTo1D(svdEvals2D);
//				
//				// calculate # of dims to retain using given data reduction factor
//				double reducFactor = setupParams.getDataReductFactor();		
//				int reducDataDims = 0;				
////				if (origData.numRows() < origData.numCols()) {
////					reducDataDims = (int)Math.round(origData.numRows() * reducFactor);
////				}
////				else {
////					reducDataDims = (int)Math.round(origData.numCols() * reducFactor);
//////					throw new NpairsjException("WARNING: input data has more " +
//////							"scans than masked voxels!");
////				}
//				reducDataDims = (int)Math.round(svdEvects.numCols() * reducFactor);
//			
//				int[] rowRange = new int[] {0, svdEvects.numRows() - 1};
//				int[] reducColRange = new int[] {0, reducDataDims - 1};
//		
//				// calculate reduced-dim eigenvectors 
//				Matrix reducDimEvects = svdEvects.subMatrix(rowRange, reducColRange);
//				
//				// get singular values by finding sqrts of evd evals
//				Matrix S = new MatrixImpl(reducDataDims, reducDataDims).getMatrix();
//				for (int i = 0; i < reducDataDims; ++i) {
//					S.setQuick(i,i, Math.sqrt(svdEvals1D[i]));
//				}
//				
//				// set feature-selected data
//				if (setupParams.normEVD()) { 
//					// normed EVD ==> set variance of each eigenvector to 1 
//					// so don't weight by eigenvalue
//					featSelData = reducDimEvects;
//				}
//				else {
//					featSelData = reducDimEvects.mult(S);	
//				}
//
//				// calculate Matrix origSpaceXformFactor - useful when transforming data 
//				// back from feature-selection to original space.
//				setOrigSpaceXformFactor(reducDimEvects, S);
//				double tTime = (System.currentTimeMillis() - sTime) / 1000;
//				
//				Npairsj.output.println("[" + tTime + " s]");				
//			}
//			
//			else {
//				featSelData = selectFeaturesEVD(setupParams, origData, 
//						setupParams.getDataReductFactor());
//			}
//		}
//		
//		return featSelData;
//	}

//	private void setOrigSpaceXformFactor(Matrix reducDimEvects, Matrix S) {
//		int reducDataDims = reducDimEvects.numCols();
//		Matrix invS = new MatrixImpl(reducDataDims, reducDataDims).getMatrix();
//		for (int i = 0; i < reducDataDims; ++i) {
//			invS.setQuick(i, i, 1 / S.getQuick(i, i));
//		}
//
//		origSpaceXformFactor = invS.mult(reducDimEvects.transpose());
//	}

	
//	/** Returns Matrix containing feature-selected data using eigenvalue decomposition.
//	 * 
//	 *  <p>Note that IDL uses EVD, too (actually eigenql of non-mc ssp mat) although it's
//	 *  called 'SVD' in IDL NPAIRS. 
//	 *  
//	 *  <p>If using unweighted ("normed") EVD,
//	then weighting matrix S is ignored in data representation P<sub>normed</sub>, i.e., 
//	<ul>P<sub>normed</sub> = V  ( = MUS<sup>-1</sup> )  instead of </ul>
//	<ul>P = VS  ( = MU ).   </ul>
// 
//<p>	This is equivalent to setting all (non-zero diagonal) elements of S to 1, 
//	and  gets rid of differences in variance across basis vectors in U; hence 
//	this can be thought of as a data "denoising" technique.
//	In normed EVD, data is still represented in U-space, hence
//	transformation back to original space is still a transformation back through U:
//<ul>		P<sub>1</sub>U<sup>t</sup> = M<sub>1</sub> still holds.
//</ul>	
//Given P<sub>normed</sub> = V = MUS<sup>-1</sup>, U<sup>t</sup> = 
//	S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup>M. 
//	
//<p>Hence we use S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup> in normed EVD where P<sup>-1</sup> is used
//	in regular EVD. 
//	
//<p>	But note that 
//<ul>	<li>P<sub>normed</sub> = V ==>  S<sup>-1</sup>(P<sub>normed</sub>)<sup>-1</sup> 
//			=  S<sup>-1</sup>V<sup>t</sup>, 
//and in regular EVD, 
//	<li>P = VS ==> P<sup>-1</sup> = S<sup>-1</sup>V<sup>t</sup>. 
//	</ul>	
//<p>Therefore, we simply store  S<sup>-1</sup>V<sup>t</sup> in origSpaceXformFactor in both cases.
//	 *   
//	 * REQUIRED: 0.0 < dataReductionFactor <= 1.0
//	 */
//	private Matrix selectFeaturesEVD(NpairsjSetupParams setupParams, Matrix data, 
//			double dataReductionFactor) throws IOException, NpairsjException {
//		
//		String taskMsg = "Running initial EVD with DRF = " + 
//				dataReductionFactor + "... ";
//		if (setupParams.normEVD()) {
//			taskMsg = "Running initial EVD (normed) with DRF = " + 
//					dataReductionFactor + "... ";
//		}
//		Npairsj.output.println(taskMsg);
//		
//		Matrix featSelData = null;
//		double reducFactor = dataReductionFactor;		
//		
//		int reducDataDims = (int)Math.round(Math.min(data.numRows(), data.numCols()) * reducFactor);
//			
///*			 Given data = M, PCA(Mt) ==> MMt = V(S^2)Vt, where S^2 == diag Matrix containing
//			 squared singular values from Matrix S in SVD(M) = VSUt along diagonal.
//			 ==> featSelData = projection P of data M onto eigenimages U, 
//			 i.e., featSelData P = MU
//			                     = VSUtU 
//			                     = VS, 
//			 where U is data.numCols() X reducDataDims(rDD);
//			      
//			 		 S is rDD X rDD;
//			 		 
//			       V is data.numRows() X rDD;
//			 			
//			 ==> P = VS is data.numRows() X rDD
//			       == dim(MU) as required */
//			
//			double sTime = System.currentTimeMillis();
//			EigenvalueDecomposition evd = null;
//			Matrix svdEvals = null;
//			Matrix svdEvects = null;
//			Matrix invSqrtEvals = null;
//			if (!matlibType.equalsIgnoreCase(matlibTypeForInitFeatSel)) {
//				//TODO: too many copies of data being made here!
//				if (matlibTypeForInitFeatSel.equalsIgnoreCase("MATLAB")) {
//					if (debug) {
//						System.out.println("Doing EVD in Matlab...");
//					}
//					MatlabMatrix mlData = new MatlabMatrix(data.toArray());
//					Npairsj.output.print("\tCreating SSP matrix from data matrix...");
//					double sspSTime = System.currentTimeMillis();
//					MatlabMatrix sspMat = mlData.mult(mlData.transpose());
//					double sspTTime = (System.currentTimeMillis() - sspSTime) / 1000;
//					Npairsj.output.println("[" + sspTTime + "s]");
//					Npairsj.output.print("\tRunning EVD on SSP matrix...");
//					double evdSTime = System.currentTimeMillis();
//					evd = sspMat.eigenvalueDecomposition();
//					double evdTTime = (System.currentTimeMillis() - evdSTime) / 1000;
//					Npairsj.output.println("[" + evdTTime + "s]");
//				}
//				if (matlibTypeForInitFeatSel.equalsIgnoreCase("COLT")) {
//					if (debug) {
//						System.out.println("Doing EVD in Colt...");
//					}
//					ColtMatrix coltData = new ColtMatrix(data.toArray());
//					Npairsj.output.print("\tCreating SSP matrix from data matrix...");
//					double sspSTime = System.currentTimeMillis();
//					ColtMatrix sspMat = coltData.mult(coltData.transpose());
//					double sspTTime = (System.currentTimeMillis() - sspSTime) / 1000;
//					Npairsj.output.println("[" + sspTTime + "s]");
//					Npairsj.output.print("\tRunning EVD on SSP matrix...");
//					double evdSTime = System.currentTimeMillis();
//					evd = sspMat.eigenvalueDecomposition();
//					double evdTTime = (System.currentTimeMillis() - evdSTime) / 1000;
//					Npairsj.output.println("[" + evdTTime + "s]");
//				}
//				svdEvals = new MatrixImpl(evd.getRealEvalMat().toArray()).getMatrix();
//				svdEvects = new MatrixImpl(evd.getEvects().toArray()).getMatrix();
//				invSqrtEvals = new MatrixImpl(evd.getInvSqrtRealEvalMat().toArray()).getMatrix();
//			}
//			else {
//				Npairsj.output.print("\tCreating SSP matrix from data matrix...");
//				double sspSTime = System.currentTimeMillis();
//				Matrix sspData = data.sspByRow();
//				double sspTTime = (System.currentTimeMillis() - sspSTime) / 1000;
//				Npairsj.output.println("[" + sspTTime + "s]");
//				Npairsj.output.print("\tRunning EVD on SSP matrix...");
//				double evdSTime = System.currentTimeMillis();
//				evd = sspData.eigenvalueDecomposition();
//				double evdTTime = (System.currentTimeMillis() - evdSTime) / 1000;
//				Npairsj.output.println("[" + evdTTime + "s]");
//				svdEvals = evd.getRealEvalMat();
//				svdEvects = evd.getEvects();
//				invSqrtEvals = evd.getInvSqrtRealEvalMat();
//			}
//
////			if (debug) {
//				double tTime = (System.currentTimeMillis() - sTime) / 1000;
//				int hr = (int)(tTime / 3600);
//				int min = (int)(tTime / 60) - (hr * 60);
//				double s = tTime - (hr * 3600) - (min * 60) ;
//				Npairsj.output.print("Total time doing EVD calculations: [" + hr + " h " +
//						min + " min ");
//				Npairsj.output.printf("%.3f", s);
//				Npairsj.output.println(" s]");
////			}
//			
//			// save evals/evects, i.e. S^2 and V (TODO: save only if user specifies to do so?):
//			String evalsFilename = setupParams.getResultsFilePrefix() + ".EVD.evals";
//			String evectsFilename = setupParams.getResultsFilePrefix() + ".EVD.evects";
//			double saveSTime = System.currentTimeMillis();
//			Npairsj.output.print("Saving EVD info to file...");
//			NpairsjIO.printToIDLFile(evd.getRealEvals(), evalsFilename);
//			svdEvects.printToFile(evectsFilename, "IDL");
//			double saveTTime = (System.currentTimeMillis() - saveSTime) / 1000;
//			Npairsj.output.println("[" + saveTTime + " s]");
//			
//			//***********************************************************
//			// Save eigims too!
//			//***********************************************************
////			System.out.print("Creating and saving svd eigims...");
////			double strtTime = System.currentTimeMillis();
////			String eigimsFilename = setupParams.resultsFilePrefix + ".EVD.eigims";
////			Matrix svdEigims = data.transpose().mult(svdEvects).mult(invSqrtEvals);
////			svdEigims.printToFile(eigimsFilename, "IDL");
////			double totTime = (System.currentTimeMillis() - strtTime) / 1000;
////			System.out.println("[" + totTime + " s]");
//			int[] rowRange = new int[] {0, data.numRows() - 1};
//			int[] reducColRange = new int[] {0, reducDataDims - 1};
//			
//			if (setupParams.normEVD()) {
//				// ignore S so featSelData = V
//				featSelData = svdEvects.subMatrix(rowRange, reducColRange);
//			}
//			else {
//				// Given data M = VSUt, consider basis space to be
//				// U, hence MU = VS, i.e., featSelData = VS
//				Matrix S = svdEvals.subMatrix(reducColRange, reducColRange);
//				for (int i = 0; i < reducDataDims; ++i) {
//					double currEval = S.get(i, i);
//					if (currEval > 0) {
//						S.set(i, i, Math.sqrt(S.get(i, i)));
//					}
//					else S.set(i, i, 0);
//				}
//
//				featSelData = svdEvects.subMatrix(rowRange, reducColRange).
//					mult(S);
//			}
//			
//			double iSTime = System.currentTimeMillis();
//			Npairsj.output.print("Calculating inverse EVD matrix ...");
//			
//			// see setOrigSpaceXformFactor()
//			origSpaceXformFactor = invSqrtEvals.subMatrix(reducColRange, reducColRange).
//					mult(svdEvects.subMatrix(rowRange, reducColRange).transpose());
//			double iTTime = (System.currentTimeMillis() - iSTime) / 1000;
//			Npairsj.output.println("[" + iTTime + "s]");
//			        					
////		}
//		
//		return featSelData;
//	}
	
	public float[] getQOffset() {
		return qOffset;
	}
	public double[] getVoxSize() {
		return voxelSize;
	}
	
	public int[] getDims() {
		return volDims3D;
	}
	
	public int[] getOrigin() {
		return origin;
	}
	
	public int[] getMaskCoords() {
		return maskCoords;
	}
	
	// for testing 
	private static void main(String[] args) {
		String setupParamsFile = args[0];
		try {
			Matrix datamat = loadDatamats(new NpairsSetupParams(setupParamsFile, true), "Colt");
			System.out.println("Size of concatenated datamat: " + datamat.numRows() +
				" X " + datamat.numCols());
			datamat.printToFile("/home/anita/Desktop/concat_dmat.2D", "IDL");
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
