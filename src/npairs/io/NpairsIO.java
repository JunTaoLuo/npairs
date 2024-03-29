package npairs.io;

import java.io.*;

import npairs.shared.matlib.*;
import npairs.utils.PCA;
import jmatlink.*;
import java.math.BigInteger;
import java.util.HashSet;
import java.util.Vector;

import pls.shared.MLFuncs;
import npairs.NpairsException;
import npairs.NpairsSetupParams;
import extern.ArrayFuncs;
import extern.niftijlib.Nifti1Dataset;

/** Contains convenience methods for reading and writing files associated with Npairsj.
 * @author anita
 *
 */
public class NpairsIO {
	
	private static boolean debug;
	
	/** Writes double array to textfile, one element per row.
	 * 
	 * @param data
	 * @param fileName (including full path)
	 * @throws IOException if error occurs writing to file
	 */
	//TODO: Possibly format data so that if length of array > some number, then
	// move into next column when printing to file.
	public static void printToIDLFile (double[] data, String fileName) throws IOException {
			
			PrintWriter pw = null;
			try {
				pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
				pw.println("1 " + data.length);
				for (int i = 0; i < data.length; ++i) {
					pw.println(data[i]);
				}
			}
//			catch (IOException e) {
//				System.err.println("Error opening file " + fileName);
//				e.printStackTrace();
//			}
			finally {
//				try {
					if (pw != null) {
						boolean error = pw.checkError();
						if (error) {
							throw new IOException("Error occurred writing to file "
									+ fileName);
						}
						pw.close();
					}
//				}
//				catch (IOException e) {
//					e.printStackTrace();
//				}
//				catch (Exception e) {
//					System.err.println("Error - could not close file " + fileName);
//					e.printStackTrace();
//				}
			}
	}
	
	
	/** Writes double array to textfile, one element per row.
	 * 
	 * @param data
	 * @param fileName (including full path)
	 * @throws IOException if error occurs writing to file
	 */
	//TODO: Possibly format data so that if length of array > some number, then
	// move into next column when printing to file.
	public static void printToFile (double[] data, String fileName) throws IOException {
			
			PrintWriter pw = null;
			try {
				pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
				
				for (int i = 0; i < data.length; ++i) {
					pw.println(data[i]);
				}
			}
//			catch (IOException e) {
//				System.err.println("Error opening file " + fileName);
//				e.printStackTrace();
//			}
			finally {
//				try {
					if (pw != null) {
						boolean error = pw.checkError();
						if (error) {
							throw new IOException("Error occurred writing to file "
									+ fileName);
						}
						pw.close();
					}
//				}
//				catch (IOException e) {
//					e.printStackTrace();
//				}
//				catch (Exception e) {
//					System.err.println("Error - could not close file " + fileName);
//					e.printStackTrace();
//				}
			}
	}
	
	/** Writes int array to textfile, one element per row.
	 * 
	 * @param data
	 * @param fileName (including full path)
	 * @throws IOException if error occurs writing to file
	 */
	//TODO: Possibly format data so that if length of array > some number, then
	// move into next column when printing to file.
	public static void printToFile (int[] data, String fileName) throws IOException {
			
			PrintWriter pw = null;
			try {
				pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
				
				for (int i = 0; i < data.length; ++i) {
					pw.println(data[i]);
				}
			}
//			catch (IOException e) {
//				System.err.println("Error opening file " + fileName);
//				e.printStackTrace();
//			}
			finally {
//				try {
					if (pw != null) {
						boolean error = pw.checkError();
						if (error) {
							throw new IOException("Error occurred writing to file "
									+ fileName);
						}
						pw.close();
					}
//				}
//				catch (IOException e) {
//					e.printStackTrace();
//				}
//				catch (Exception e) {
//					System.err.println("Error - could not close file " + fileName);
//					e.printStackTrace();
//				}
			}
	}
	
	public static void printToIDLFile(int[] data, String fileName) throws IOException {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
			
			pw.println("1 " + data.length);
			for (int i = 0; i < data.length; ++i) {
				pw.println(data[i]);
			}
		}
//		catch (IOException e) {
//			System.err.println("Error opening file " + fileName);
//			e.printStackTrace();
//		}
		finally {
//			try {
				if (pw != null) {
					boolean error = pw.checkError();
					if (error) {
						throw new IOException("Error occurred writing to file "
								+ fileName);
					}
					pw.close();
				}
//			}
//			catch (IOException e) {
//				e.printStackTrace();
//			}
//			catch (Exception e) {
//				System.err.println("Error - could not close file " + fileName);
//				e.printStackTrace();
//			}
		}
	}
	
	/** Prints input int array to file s.t. Jon Anderson's read_matrix.pro reads
	 *  data as row array.
	 * @param data
	 * @param fileName
	 */
	public static void printRowToIDLFile(int[] data, String fileName) throws IOException {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
			
			pw.println(data.length + " 1");
			for (int i = 0; i < data.length; ++i) {
				pw.println(data[i]);
			}
		}
//		catch (IOException e) {
//			System.err.println("Error opening file " + fileName);
//			e.printStackTrace();
//		}
		finally {
//			try {
				if (pw != null) {
					boolean error = pw.checkError();
					if (error) {
						throw new IOException("Error occurred writing to file "
								+ fileName);
					}
					pw.close();
				}
//			}
//			catch (IOException e) {
//				e.printStackTrace();
//			}
//			catch (Exception e) {
//				System.err.println("Error - could not close file " + fileName);
//				e.printStackTrace();
//			}
		}
	}
	
	/** Read 2D array from IDL-compatible 
	 *  textfile (data in textfile must be in read_matrix.pro 
	 *  format - see J. Anderson's IDL read_matrix.pro code).
	 *
	 */
	 public static int[][] readIntsFromIDLFile(String filename) throws
	 		IOException {
		 double[][] dataAsDbls = readFromIDLFile(filename);
		 int[][] dataAsInts = (int[][])ArrayFuncs.convertArray(dataAsDbls, int.class);
		 return dataAsInts;
	 }
	 
	/** Read 2D array from IDL-compatible 
	 *  textfile (data in textfile must be in read_matrix.pro 
	 *  format - see J. Anderson's IDL read_matrix.pro code).
	 *
	 */
	 public static double[][] readFromIDLFile(String filename) throws IOException {

		 double[][] data = null;
		 BufferedReader br = new BufferedReader(new FileReader(filename));
		 try {
			 String line = br.readLine();
			 String[] lineSplit = line.trim().split("\\s");
			 Vector<String> dataDims = new Vector<String>();
			 for (int i = 0; i < lineSplit.length; ++i) {
				 if (!lineSplit[i].equals("")) {
					 dataDims.add(lineSplit[i]);
				 }
			 }
			 if (dataDims.size() == 1) {
				 dataDims.add("1");
			 }
			 if (dataDims.size()!= 2) {
				 throw new IOException("Error - input file " + filename + " must be in "
						 + "read_matrix.pro format (see J.Anderson's NPAIRS"
						 + "read_matrix.pro code).");
			 }

			 int nCols = (new Integer(dataDims.get(0))).intValue();
			 int nRows = (new Integer(dataDims.get(1))).intValue();

			 // TODO: FIX THIS PART: Get data line-by-line into
			 //         2D data array.
			 data = new double[nRows][nCols];
			 int r = 0;
			 int c = 0;
			 while ((line = br.readLine()) != null){
				 String[] nextRowDataTemp = line.trim().split("\\s");
				 Vector<String> vNextRowData = new Vector<String>();
				 for (int i = 0; i < nextRowDataTemp.length; ++i) {
					 if (!nextRowDataTemp[i].equals("")) {
						 vNextRowData.add(nextRowDataTemp[i]);
					 }
				 }
				 for (int j = 0; j < vNextRowData.size(); ++j) {
					 data[r][c] = Double.parseDouble(vNextRowData.get(j));
					 if (c < nCols - 1) {
						 ++c;
					 }
					 else {
					 // move to next row
						 c = 0;
						 ++r;
					 }
				 }
			 }	 
		 }
		 catch (Exception e) {
			 throw new IOException("Error loading data from file " +
					 filename + ". " + e.getMessage());
		 }
		 finally {
			 br.close();
		 }
		 return data;
	 }
	 
	 /** Read array of doubles from 
	  *  textfile; reads data in order one row at a time
	  *  from beginning of file (does not matter how many elements
	  *  are given in each row)
	  * @param name of file containing data as doubles
	  */
	 public static double[] readDblsFromFile(String filename) throws IOException {

		 double[] data = null;
		 BufferedReader br = new BufferedReader(new FileReader(filename));
		 try {
			 Vector<Double> vNextRowDataDbls = new Vector<Double>();
			 String line;
			 while ((line = br.readLine()) != null){
				 String[] nextRowDataTemp = line.trim().split("\\s");
				 for (int i = 0; i < nextRowDataTemp.length; ++i) {
					 if (!nextRowDataTemp[i].equals("")) {
						 vNextRowDataDbls.add(new Double(nextRowDataTemp[i]));
					 }
				 }
			 }
			 data = new double[vNextRowDataDbls.size()];
			 for (int j = 0; j < vNextRowDataDbls.size(); ++j) {
				 data[j] = vNextRowDataDbls.get(j);
			 }
		 }
		 finally {
			 br.close();
		 }
		 return data;
	 }
	 
	 
	 /** Read array of ints from 
	  *  textfile; reads data in order one row at a time
	  *  from beginning of file (does not matter how many elements
	  *  are given in each row)
	  * @param name of file containing data as ints
	  */
	 public static int[] readIntsFromFile(String filename) throws IOException {

		 double[] dataAsDbls = readDblsFromFile(filename);
		 int[] dataAsInts = new int[dataAsDbls.length];
		 for (int i = 0; i < dataAsDbls.length; ++i) {
			 dataAsInts[i] = (int)dataAsDbls[i];
		 }
		 return dataAsInts;
	 }
	 
//	 /** Read 2d array of integers from textfile. Each column may contain a 
//	  *  different number of integers.
//	  * @param name of file containing 2d integer data
//	  * @return 2d int array int[numCols][] 
//	  */
//	 public static int[][] read2DIntsFromFile(String filename) throws IOException {
//		 int[][] data = null;
//		 BufferedReader br = new BufferedReader(new FileReader(filename));
//		 try {
//			 Vector<Integer> vNextRowDataInts = new Vector<Integer>();
//			 String[] nextRowTmp;
//			 String line = br.readLine();
//			 if (line != null) {
//				 nextRowTmp = line.trim().split("\\s");
//				 int nCols = nextRowTmp.length;
//				 data = new int[nCols][];
//				 
//			 }
//			 while (line != null) {
//				 nextRowTmp = line.trim().split("\\s");
//				 for (int i = 0; i < nextRowTmp.length; ++i) {
//					 if (!nextRowTmp[i].equals("")) {
//						 vNextRowDataInts.add(new Integer(nextRowTmp[i]));
//					 }
//				 }
//				 line = br.readLine();
//			 }
//		
//			 
//		 }
//		 finally {
//			 br.close();
//		 }
//		 return data;
//	 }

	 /** Prints 2D int array to IDL-compatible textfile
	 * 
	 * @param data
	 * @param filename
	 * 
	 */
	public static void printToIDLFile(int[][] data, String filename) throws IOException {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			
			pw.println(data[0].length + " " + data.length);
			for (int row = 0; row < data.length; ++row) {
				for (int col = 0; col < data[row].length; ++col) {
					pw.print(data[row][col] + " ");
				}
				pw.println();
			}
		}
//		catch (IOException e) {
//			System.err.println("Error opening file " + filename);
//			e.printStackTrace();
//		}
		finally {
//			try {
				if (pw != null) {
					boolean error = pw.checkError();
					if (error) {
						throw new IOException("Error occurred writing to file "
								+ filename);
					}
					pw.close();
				}
//			}
//			catch (IOException e) {
//				e.printStackTrace();
//			}
//			catch (Exception e) {
//				System.err.println("Error - could not close file " + filename);
//				e.printStackTrace();
//			}
		}
	}
	
	 /** Prints 2D double array to IDL-compatible textfile
	 * 
	 * @param data
	 * @param filename
	 * 
	 */
	public static void printToIDLFile(double[][] data, String filename) throws IOException {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
			
			pw.println(data[0].length + " " + data.length);
			for (int row = 0; row < data.length; ++row) {
				for (int col = 0; col < data[row].length; ++col) {
					pw.print(data[row][col] + " ");
				}
				pw.println();
			}
		}
//		catch (IOException e) {
//			System.err.println("Error opening file " + filename);
//			e.printStackTrace();
//		}
		finally {
//			try {
				if (pw != null) {
					boolean error = pw.checkError();
					if (error) {
						throw new IOException("Error occurred writing to file "
								+ filename);
					}
					pw.close();
				}
//			}
//			catch (IOException e) {
//				e.printStackTrace();
//			}
//			catch (Exception e) {
//				System.err.println("Error - could not close file " + filename);
//				e.printStackTrace();
//			}
		}
	}

	/** Prints 2D int array to stdout
	 * 
	 * @param array
	 */
	public static void print(int[][] data) {
		for (int row = 0; row < data.length; ++row) {
			for (int col = 0; col < data[row].length; ++col) {
				System.out.print(data[row][col] + " ");
			}
			System.out.println();
		}
	}
	
	/** Prints 2D BigInteger array to stdout
	 * 
	 * @param array
	 */
	public static void print(BigInteger[][] data) {
		for (int row = 0; row < data.length; ++row) {
			for (int col = 0; col < data[row].length; ++col) {
				System.out.print(data[row][col] + " ");
			}
			System.out.println();
		}
	}
	
	/** Prints 2D HashSet array to stdout
	 * 
	 * @param array
	 */
	public static void print(HashSet[][] data) {
		for (int row = 0; row < data.length; ++row) {
			for (int col = 0; col < data[row].length; ++col) {
				System.out.print(data[row][col] + " ");
			}
			System.out.println();
		}
	}
	
	/** Prints 2D double array to stdout
	 * 
	 * @param array
	 */
	public static void print(double[][] data) {
		for (int row = 0; row < data.length; ++row) {
			for (int col = 0; col < data[row].length; ++col) {
				System.out.print(data[row][col] + " ");
			}
			System.out.println();
		}
	}
	
	/** Prints 2D double array to stdout.  Prints 3 digits after 
	 *  decimal place.
	 * 
	 * @param array
	 */
	public static void printf(double[][] data) {
		for (int row = 0; row < data.length; ++row) {
			for (int col = 0; col < data[row].length; ++col) {
				System.out.printf("%.3f ", data[row][col]);
			}
			System.out.println();
		}
	}
	
	
	/** Prints 1D int array to stdout
	 * 
	 * @param array
	 */
	public static void print(int[] data) {
		print(data, System.out);
	}
	
	/** Prints 1D int array to given output stream
	 * 
	 * @param array
	 */
	public static void print(int[] data, PrintStream out) {
		for (int i = 0; i < data.length; ++i) {
			
			out.print(data[i] + " ");
		}
		out.println();
	}
	
	/** Prints 1D double array to stdout
	 * 
	 * @param array
	 */
	public static void print(double[] data) {
		for (int i = 0; i < data.length; ++i) {
			
			System.out.print(data[i] + " ");
		}
		System.out.println();
	}
	
	/** Prints 1D String array to stdout
	 * 
	 * @param array
	 */
	public static void print(String[] data) {
		for (int i = 0; i < data.length; ++i) {
			System.out.print(data[i] + " ");
		}
		System.out.println();
	}
	
	/** Prints 1D Boolean array to stdout
	 * 
	 * @param array
	 */
	public static void print(Boolean[] data) {
		for (int i = 0; i < data.length; ++i) {
			System.out.print(data[i] + " ");
		}
		System.out.println();
	}
	
	/** Prints 1D String array to given output stream
	 *
	 * @param array
	 */
	public static void print(String[] data, PrintStream outStream) {
		for (int i = 0; i < data.length; ++i) {
			outStream.print(data[i] + " ");
		}
		outStream.println();
	}
	
	/** Returns filename stripped of suffix (returns input String if it has no suffix)
	 * 
	 * @param String filename
	 * @return String filename stripped of suffix
	 */
	public static String stripSuff(String filename) {
//		System.out.println("Path separator: [" + System.getProperty("file.separator") + "]");
		int lastPathSepLoc = filename.lastIndexOf(System.getProperty("file.separator"));
//		System.out.println("Last path separator index: " + lastPathSepLoc);
		int suffSepLoc = filename.lastIndexOf(".");
//		System.out.println("Suffix separator index: " + suffSepLoc);
		String strippedFilename = filename;
		if (suffSepLoc > lastPathSepLoc) {
			strippedFilename = filename.substring(0, suffSepLoc);
		}
		return strippedFilename;
	}
	
	/** Copies input origFile to newFile (overwriting whatever already lives
	 *  in newFile).  If newFile and origFile are the same file, do nothing.
	 *  
	 * @param origFile - String containing name of file to be copied
	 * @param newFile  - String containing name of file to contain copy of origFile
	 */
	public static void copyFile (String origFile, String newFile) throws FileNotFoundException,
		IOException {
	
		File fOrig = new File(origFile);
		File fNew = new File(newFile);
		String pathOrig = fOrig.getAbsolutePath();
		String pathNew = fNew.getAbsolutePath();
		if (pathOrig.equals(pathNew)) {
			return;
		}
		
		InputStream in = new FileInputStream(fOrig);
		OutputStream out = new FileOutputStream(fNew);

		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0){
			out.write(buf, 0, len);
		}

		in.close();
		out.close();

		if (debug) {
			System.out.println("File copied.");
		}
	}
	
	
	// testing readDblsFromFile(String filename)...
	public static void main (String[] args) {
		// input arg is name of file containing doubles
		String filename = args[0];
		try {
			double[] data = readDblsFromFile(filename);
			System.out.println("Data in file '" + filename + "':");
			NpairsIO.print(data);
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	// testing stripSuff...
	public static void main0 (String[] args) {
		for (int i = 0; i < args.length; ++i) {
			String inFilename = args[i];
			String strippedFilename = stripSuff(inFilename);
			System.out.println("Input filename: " + inFilename);
			System.out.println("Filename with suffix stripped: " + strippedFilename);
		}
	}
	
	public static void main1 (String[] args) {
//		String filename = "/rri_disks/haier/anita/npairsj/npairsj_code_testing/javaPrintTest/dataArray.dat";
//		String filename2 = "/astraea/strother_lab/anita/workspace/newNpairsj2/dataArrayTest.dat";
//		int size = 10000;
//		double[] data = new double[size];
//		for (int i = 0; i < size; ++i) {
//			data[i] = (double)i;
//		}
		//System.out.println("Data array: ");
		//utils_tests.PCATest.printArray(data);
//		double startTime = System.currentTimeMillis();
//		printToFile(data, filename2);
//		double endTime = System.currentTimeMillis() - startTime;
//		System.out.println("Time to Print file: " + endTime);
//		startTime = System.currentTimeMillis();
//		writeToFile(data, filename2+".write");
//		endTime = System.currentTimeMillis() - startTime;
//		System.out.println("Time to Write file: " + endTime);
		
		String coltPrintFile = "/astraea/strother_lab/anita/workspace/newNpairsj2/coltPrintMat.mat";
		//String coltWriteFile = "/astraea/strother_lab/anita/workspace/newNpairsj2/coltWriteMat.mat";
		String matlabPrintFile = "/astraea/strother_lab/anita/workspace/newNpairsj2/matlabPrintMat.mat";
		//String matlabWriteFile = "/astraea/strother_lab/anita/workspace/newNpairsj2/matlabWriteMat.mat";
		
		int nRows = 10000;
		int nCols = 100;
		ColtMatrix cmat = new ColtMatrix(nRows, nCols);
		cmat.setRandom();
		JMatLink eng = new JMatLink();
		eng.engOpen();
		MatlabMatrix mmat = new MatlabMatrix(nRows, nCols);
		mmat.setRandom();
		
		double cStartTime = System.currentTimeMillis();
		cmat.printToFile(coltPrintFile, "DEFAULT");
		double cEndTime = System.currentTimeMillis() - cStartTime;
		System.out.println("Time to print colt file: " + cEndTime);
		cStartTime = System.currentTimeMillis();
		//cmat.writeToFile(coltWriteFile, "DEFAULT");
		cEndTime = System.currentTimeMillis() - cStartTime;
		System.out.println("Time to write colt file: " + cEndTime);
		double mlStartTime = System.currentTimeMillis();
		mmat.printToFile(matlabPrintFile, "DEFAULT");
		double mlEndTime = System.currentTimeMillis() - mlStartTime;
		System.out.println("Time to print matlab file: " + mlEndTime);
		mlStartTime = System.currentTimeMillis();
		//mmat.writeToFile(matlabWriteFile, "DEFAULT");
		mlEndTime = System.currentTimeMillis() - mlStartTime;
		System.out.println("Time to write matlab file: " + mlEndTime);
	}


	public static void print(short[] data)  {
		
		for (int i = 0; i < data.length; ++i) {
				
			System.out.print(data[i] + " ");
		}
		System.out.println();
	}


	public static void print(float[] data) {
		for (int i = 0; i < data.length; ++i) {
			System.out.print(data[i] + " ");
		}
		System.out.println();		
	}


	public static void print(Integer[] objIntArray) {
		for (int i = 0; i < objIntArray.length; ++i) {
			System.out.print(objIntArray[i] + " ");
		}
		System.out.println();
	}


	public static void print(boolean[] data) {
		for (int i = 0; i < data.length; ++i) {
			System.out.print(data[i] + " ");
		}
		System.out.println();	
	}
	
	/** Saves data reconstructed from PCA results. Data is mean-centred and 
	 *  transformed back to original image space (if PCA was done in some other space).  
	 *  Data reconstruction is typically done using dimension-reduced PCA results, i.e., 
	 *  data D = P<sup>*</sup>V<sup>*t</sup>), 
	 *  where 
	 *  <ul>
	 *  <li>* means dimension-reduced matrix, 
	 *  <li>P = PC scores, 
	 *  <li>V = PC eigenvectors.
	 *  <p>REQUIRED: PCA lives in original space (i.e., eigenvectors have been
	 *  rotated back into original space before calling this method if PCA was done in 
	 *  some other space).
	 * @param pca PCA from which to reconstruct and save data matrix
	 * @param nsp NpairsjSetupParams containing data/results filenames (data filename
	 * 		required for header info)
	 * @param ndl NpairsDataLoader containing mask coordinates and other data details
	 * @param logStream where to log info
	 *  
	 */
	public static void saveDataPostPCA(PCA pca, NpairsSetupParams nsp, NpairsDataLoader ndl, 
			PrintStream logStream) throws NpairsException, IOException {
		double sTime = 0;
		if (debug) {
			sTime = System.currentTimeMillis();
			logStream.print("\tCreating dim-reduced (denoised) input data...");
		}
		if (!pca.pcaInOrigSpace()) {
			throw new NpairsException("Must rotate PCA eigenvectors into original image " +
					"space before saving denoised (post-PCA) data.");
		}
		
		Matrix dimRedData = pca.getPCScores().mult(pca.getEvects().transpose());
	
		if (debug) {
			double tTime = (System.currentTimeMillis() - sTime) / 1000;
			logStream.println("[" + tTime + " s]");
			logStream.println("\tSize of denoised data matrix: " +
					dimRedData.numRows() + " X " + dimRedData.numCols());
		}

		int[] maskCoords = ndl.getMaskCoords(); // 0-relative
		
		// initialize nifti structure:
		Nifti1Dataset nifti = new Nifti1Dataset();
		// set header filename first and include .nii ext to set dataset file type to .nii
		// (note: also need to write header before writing data when saving .nii file)
		String saveName = nsp.getResultsFilePrefix() + ".DATA.POST-PCA.nii";
		nifti.setHeaderFilename(saveName);
		nifti.setDataFilename(saveName);
		
		// copy header info from first input data vol (and ndl info, which was also
		// read in from first input data vol) and write it to disk
		Nifti1Dataset vol1Nifti = new Nifti1Dataset(nsp.getDataFilenames()[0]);
		nifti.copyHeader(vol1Nifti);
		nifti.setHeaderFilename(saveName);
		nifti.setDataFilename(saveName);
		nifti.setDatatype((short)64);
		int[] volDims3D = ndl.getDims(); // might as well get 3D dims from ndl 
		int nVols = dimRedData.numRows();
		int xDim = volDims3D[0];
		int yDim = volDims3D[1];
		int zDim = volDims3D[2];
		nifti.setDims((short)4, (short)xDim, (short)yDim, (short)zDim,
				(short)nVols, (short)0, (short)0, (short)0);
		nifti.writeHeader();
		
		// set dimRedData rows into full 3d vols ([Z][Y][X] as per Nifti1Dataset specs)
		double[][][] currVol3D = new double[zDim][yDim][xDim];
		double[] currVol1D = new double[zDim * yDim * xDim];
		
		
		for (short v = 0; v < nVols; ++v) {
			currVol1D = MLFuncs.setVals(currVol1D, maskCoords, dimRedData.getRowQuick(v));
			// turn 1D array into 3D vol
			for (int z = 0; z < zDim; ++z) {
				for (int y = 0; y < yDim; ++y) {
					for (int x = 0; x < xDim; ++x) {
						currVol3D[z][y][x] = currVol1D[z*yDim*xDim + y*xDim + x];
					}
				}
			}
			
			nifti.writeVol(currVol3D, v);
		}
	}
	
	/** Saves PC results to IDL-format files.
	 * @param pcaSavePref prefix (including path) of saved PCA files.  
	 * @param pcaDims - dims to save.  If null, save all.
	 * @param saveAsSplit boolean
	 * @param splitNum
	 * @param splitHalf
	 * @param pca TODO
	 */
	
	public static void savePCAResultsIDL(String pcaSavePref, 
				 int[] pcaDims,
			     boolean saveAsSplit, 
			     int splitNum,
			     int splitHalf, PCA pca) {
//		String pcaEvalFile = "";
		String pcaEvectFile = "";
		String pcaScoreFile = "";
		if (saveAsSplit ) {					
			pcaEvectFile = pcaSavePref + ".PCA." + splitNum + "." + splitHalf
			+ ".evect";
			pcaScoreFile = pcaSavePref + ".PCA." + splitNum + "." + splitHalf 
			+ ".pcScore";

		}
		else {
			pcaEvectFile = pcaSavePref + ".PCA.ALL.evect";
			pcaScoreFile = pcaSavePref + ".PCA.ALL.pcScore";

		}

		if (pcaDims == null) {
			// save all pc dims 
			pca.getEvects().printToFile(pcaEvectFile, "IDL");
			pca.getPCScores().printToFile(pcaScoreFile, "IDL");
		}
		else {
			Matrix subsetEvects = pca.getEvects().subMatrixCols(pcaDims);
			subsetEvects.printToFile(pcaEvectFile, "IDL");
			Matrix subsetPCScores = pca.getPCScores().subMatrixCols(pcaDims);
			subsetPCScores.printToFile(pcaScoreFile, "IDL");
		}	
	}
}

	

