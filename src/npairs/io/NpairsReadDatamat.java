package npairs.io;

import java.util.ArrayList;
import java.util.Map;

import extern.NewMatFileReader;

import com.jmatio.io.MatFileFilter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLStructure;
import com.jmatio.types.MLChar;

import java.io.File;
import java.io.IOException;

import pls.shared.MLFuncs;
import npairs.NpairsException;
import npairs.shared.matlib.Matrix;
import npairs.shared.matlib.MatrixImpl;

/** Reads PLS-style datamat files into memory, arranging the data in NPAIRS instead
 *  of PLS format.
 *  
 *  PLS datamats are 2D arrays of data with dims
 *   	no. conditions X (no. voxels * window size)
 *  where window size is the number of consecutive timepoints included for each condition 
 *  NPAIRS datamats are rearranged so dims are
 *  	(no. conditions * window size) X no. voxels.
 *  So in NPAIRS, each row of data represents a brain volume.  
 *  (Note that 'conditions' in the above definition might actually mean 
 *  conditions*runs[*no. onsets] if runs weren't merged [and/or single subject analysis 
 *  was selected] when creating PLS datamat.)
 *  
 *  REQUIRED: MatrixImpl.matlibType is set outside of this method.
 *  
 * @author anita
 *
 */
public class NpairsReadDatamat {
	
	boolean debug = false;
	Matrix npairsDatamat;
	String dMatFilename;
	
	public NpairsReadDatamat (String dMatFilename, int[] conditionSelection) throws IOException, 
			NpairsException {
		this.dMatFilename = dMatFilename;
		System.out.println("Curr datamat file: " + dMatFilename);
		Map<String, MLArray> datamatMap = new NewMatFileReader(dMatFilename, 
				new MatFileFilter(new String[]{"st_datamat", "st_evt_list",
						"st_win_size", "create_datamat_info", "st_sessionFile",
						"session_info"})).getContent();
		double[][] dMatUnformatted = ((MLDouble)datamatMap.get("st_datamat")).getArray();
		int[] eventList = ((MLDouble)datamatMap.get("st_evt_list")).getIntArray()[0];
		int winSize = ((MLDouble)datamatMap.get("st_win_size")).get(0, 0).intValue();
		
		MLStructure dMatInfo = ((MLStructure)datamatMap.get("create_datamat_info"));
		int[] runIdx = null;
		int noRuns = 0;
		
		// If using sessiondata.mat combo file, num_runs info will be in 
		// current file
		if (dMatFilename.endsWith("sessiondata.mat")) {
			MLStructure sessInfo = (MLStructure)datamatMap.get("session_info");
			noRuns = ((MLDouble)sessInfo.getField("num_runs")).get(0,0).intValue();
		} 
		else {
			try {
				runIdx = ((MLDouble) dMatInfo.getField("run_idx"))
						.getIntFirstRowOfArray();
			} catch (NullPointerException npe) {
				// Java PLS calls it 'runIdx' and matlab PLS
				// calls it run_idx.
				try {
					runIdx = ((MLDouble) dMatInfo.getField("runIdx"))
							.getIntFirstRowOfArray();
				} catch (NullPointerException npe2) {
					// But matlab PLS from v 5.1105021 (May 2, 2011) ff. doesn't
					// include
					// this variable at all so must get noRuns info from session
					// file
					// instead. If run_idx/runIdx variable is missing, no runs
					// are excluded
					// from analysis.
					//

					String currSessFile = ((MLChar) datamatMap
							.get("st_sessionFile")).getString(0);
					// path for session file saved in datamat may be wrong so
					// use
					// datamat's current path instead (since session file and
					// datamat
					// are expected to live in same directory)
					String baseName = (new File(currSessFile)).getName();
					String currPathSessFile = (new File(dMatFilename))
							.getParent()
							+ System.getProperty("file.separator")
							+ baseName;
					System.out
							.println("Curr session file: " + currPathSessFile);
					MLStructure sessInfo = (MLStructure) new NewMatFileReader(
							currPathSessFile).getContent().get("session_info");
					noRuns = ((MLDouble) sessInfo.getField("num_runs")).get(0,
							0).intValue();
				}
			}
		}
		boolean mergeRuns = ((MLDouble)dMatInfo.getField(
				"merge_across_runs")).get(0,0) == 1;
		boolean singSubj = ((MLDouble)dMatInfo.getField(
				"single_subject_analysis")).get(0,0) == 1;
		
		if (runIdx != null) {
			noRuns = runIdx.length;
		}
		int noConds = MLFuncs.unique(eventList).length;
		if (!mergeRuns) {
			// each run treated as separate condition in eventList
			noConds /= noRuns; 
		}
		if (noConds != conditionSelection.length) {
			throw new IllegalArgumentException("Number of conditions in condition selection" +
					" array does not match number of conditions in datamat.");
		}
		
		int nDRows = eventList.length; // no. of elements in eventList == no. of rows in datamat
		int[] inclRows = new int[nDRows];
		int nInclRows = 0;
		for (int i = 0; i < nDRows; ++i) {
			// if runs aren't merged, each run gets separate event labels but we want same class
			// labels for each run
			if (!mergeRuns) {
			    for (int c = 0; c < noConds; ++c) {
			    	for (int r = 0; r < noRuns; ++r) {
			    		eventList[c*noRuns + r] = c + 1;
			    	}
			    }
			}
			int currIdx = eventList[i] - 1; // events are labelled consecutively beginning with 1
			
			if (conditionSelection[currIdx] == 1) {
				inclRows[i] = 1;
				++nInclRows;
			}
			else inclRows[i] = 0;
		}
		int nDCols = dMatUnformatted[0].length;
		double[][] dMatUnformattedExclConds = new double[nInclRows][nDCols];
		int currInclRow = 0;
		for (int i = 0; i < nDRows; ++i) {
			if (inclRows[i] == 1) {
				dMatUnformattedExclConds[currInclRow] = dMatUnformatted[i];
				++currInclRow;
			}
		}
		
		if (debug) {
			System.out.println("Dmat size (before excl. conds): " + nDRows + " X " + nDCols);
			System.out.println("Dmat size (after excl. conds): " + nInclRows + " X " + nDCols);
			System.out.println("Win size: " + winSize);
			System.out.println("Size of event list: " + eventList.length);
			System.out.println("Event list: ");
			npairs.io.NpairsIO.print(eventList);
			System.out.println("Run index list: ");
			npairs.io.NpairsIO.print(runIdx);
			System.out.println("Merge across runs? " + mergeRuns);
			System.out.println("Single subject analysis? " + singSubj);
			System.out.println("No. runs: " + noRuns);
			System.out.println("No. conditions: " + noConds);

			int noMaskedVox = nDCols / winSize;
			System.out.println("No. masked voxels: " + noMaskedVox);
		}
		// reshape data so lags are included in rows instead of columns
		
			npairsDatamat = new MatrixImpl(reshapePls2Npairs(dMatUnformattedExclConds, winSize)).getMatrix();
			if (debug) {
				System.out.println("Size of reshaped data: " + npairsDatamat.numRows() + 
					" X " + npairsDatamat.numCols());
			}
	}
	
	/** Returns new 2D double array containing input data reshaped so each lag
	 *  image is contained in a new row - i.e., instead of data being
	 *  R X (nVox * winSize), returned data is (R * winSize) X nVox.
	 *  E.g., if winSize is 8, then instead of 8 lag images all being
	 *  stored in a single data row, there will be 8 rows, each containing
	 *  the image for the corresponding lag.  
	 * @param origData: R X (num voxels * winSize).  Each row of data is ordered
	 *  as follows: vox1 lag1, vox1 lag2, ..., vox1 lagk, vox2 lag1, vox2 lag2, etc etc
	 *  where k = winSize
	 * @param winSize
	 * @return reordered data array
	 * @throws NpairsException
	 */
	private double[][] reshapePls2Npairs(double[][] data, int winSize) 
			throws NpairsException {
		int nRows = data.length;
		int nCols = data[0].length;
		int nVox = nCols / winSize;
		if (nCols % winSize != 0) {
			throw new NpairsException("Incompatible data and window sizes");
		}
		
		double[][] newData = new double[nRows * winSize][nVox];
		for (int r = 0; r < nRows; ++r) {
			for (int v = 0; v < nVox; ++v) {
				for (int w = 0; w < winSize; ++w) {
					newData[winSize*r + w][v] = data[r][winSize*v + w];
				}
			}
		}
		return newData;
	}
	
	protected Matrix getDatamat() {
		return npairsDatamat;
	}
	
	/** Returns subset of datamat located at input coordinates 
	 * 
	 * @param locs
	 * @return subset datamat
	 * @throws NpairsException if locs contains invalid index value; IOException if
	 * 			datamat could not be read 
	 */
	protected Matrix getDatamat(int[] locs) throws IOException, NpairsException {
		int[] stCoords = ((MLDouble)new NewMatFileReader(dMatFilename).getMLArray("st_coords")).
			getIntFirstRowOfArray();
		ArrayList<Integer> stCoordsAL = new ArrayList<Integer>(stCoords.length);
		for (int i : stCoords) { 
			stCoordsAL.add(new Integer(i));
		}
		ArrayList<Integer> locsAL = new ArrayList<Integer>(locs.length);
		for (int j : locs) {
			locsAL.add(new Integer(j));
		}
		
		int[] locCols = new int[locs.length];
		int c = 0;
		for (Integer loc : locsAL) {
			int locCol = stCoordsAL.indexOf(loc);
			if (locCol == -1) {
				throw new NpairsException("Index value in input argument array not valid masked data index");
			}
			locCols[c] = locCol;
			++c;
		}
		
		return npairsDatamat.subMatrixCols(locCols);
	}

//	private void checkValid(int[] locs) throws FileNotFoundException,
//			IOException, NpairsjException {
//		int[] stCoords = ((MLDouble)new NewMatFileReader(dMatFilename).getMLArray("st_coords")).
//			getIntFirstRowOfArray();
//		ArrayList<Integer> stCoordsAL = new ArrayList<Integer>(stCoords.length);
//		for (int i : stCoords) { 
//			stCoordsAL.add(new Integer(i));
//		}
//		ArrayList<Integer> locsAL = new ArrayList<Integer>(locs.length);
//		for (int j : locs) {
//			locsAL.add(new Integer(j));
//		}
//		if (!stCoordsAL.containsAll(locsAL)) {
//			throw new NpairsjException("Index value in input argument array not valid masked data index");
//		}
//	}

//	public static void main(String[] args) {
//		String filename = args[0];
//		System.out.println("Reading file: " + filename);
//		try {
//			NpairsReadDatamat nrd = new NpairsReadDatamat(filename);
//	
//		}
//		catch (IOException e) {
//			e.printStackTrace();
//		}
//		ColtMatrix cmat = new ColtMatrix(5,20);
//		cmat.setRandom();
//		System.out.println("Test matrix: ");
//		cmat.printf();
//		int winSize = 4;
//		System.out.println("Window size: " + winSize);
//		try {
//			double[][] reshapedData = reshapePls2Npairs(cmat.toArray(), winSize);
//			System.out.println("Reshaped matrix: ");
//			NpairsjIO.printf(reshapedData);
//		}
//		catch (NpairsjException ne) {
//			ne.printStackTrace();
//		}
//		
//	}
}

