package npairs.io;

//import npairs.shared.matlib.*;
import java.io.*;
import java.util.*;

import com.jmatio.io.MatFileFilter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLStructure;

import npairs.NpairsException;
import npairs.NpairsSetupParams;
import npairs.shared.matlib.Matrix;
import npairs.shared.matlib.MatrixException;
import npairs.shared.matlib.MatrixImpl;
import extern.NewMatFileReader;
import extern.niftijlib.*;
//import npairs.NpairsjException;
import pls.shared.MLFuncs;
import npairs.io.NpairsIO;

/** Convenience class of static I/O methods for reading and writing Nifti data. (Also reads 
 * ANALYZE format files.)
 * 
 * <p>Uses {@link extern.niftijlib.Nifti1Dataset}.
 * 
 * @author Anita Oder
 *
 */
public class NiftiIO {

	final static boolean debug = false;
	// Valid DATA_TYPES are specified in Nifti-1 documentation; see
	// http://Nifti.nimh.nih.gov/pub/dist/src/Niftilib/Nifti1.h
	final static int[] DATA_TYPES = { 0, 1, 2, 4, 8, 16, 32, 64, 128, 255, 
								256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2304 }; 
	
	
	/** Writes 3D volume into Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}.  
	 *  File is of <i>.nii</i> format iff filename contains <i>.nii</i> extension; otherwise
	 *  data is saved into pair of files with <i>.img</i>/<i>.hdr</i> extensions.
	 *  
	 * @param vol3D - data to be saved as doubles; data must be in [Z][Y][X] order
	 * @param datatype - 
	 *             <p>Edited description of valid data types from 
	 *             <i>http://Nifti.nimh.nih.gov/pub/dist/src/Niftilib/Nifti1.h</i>:
	 *          
     *              <p>/*--- the original ANALYZE 7.5 type codes ---*
	                 				 
                  	<p>0     /* [unknown]         *
                  	<br>1     /* binary (1 bit/voxel)         *
          		 	<br>2     /* unsigned char (8 bits/voxel) *
          			<br>4     /* signed short (16 bits/voxel) *
           			<br>8     /* signed int (32 bits/voxel)   *
                 	<br>16    /* float (32 bits/voxel)        *
             		<br>32    /* complex (64 bits/voxel)      *
                 	<br>64    /* double (64 bits/voxel)       *
                   	<br>128   /* RGB triple (24 bits/voxel)   *
                	<br>255   /* not very useful (?)          *

 
                    <p>/*------------------- new codes for Nifti ---*
                 	<br>256   /* signed char (8 bits)         *
              		<br>512   /* unsigned short (16 bits)     *
              		<br>768   /* unsigned int (32 bits)       *
               		<br>1024  /* long long (64 bits)          *
             		<br>1280  /* unsigned long long (64 bits) *
            		<br>1536  /* long double (128 bits)       *
         			<br>1792  /* double pair (128 bits)       *
        			<br>2048  /* long double pair (256 bits)  *
              		<br>2304  /* 4 byte RGBA (32 bits/voxel)  *
     * @param voxSize 3-element double array containing voxel size of data {Xsize, Ysize, Zsize}
	 * @param filename
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @see #writeVol(double[][][], String)
	 * @see #writeVol(double[][][], int, String)
	 * @see #writeVol(double[][][], double[], String)
	 * @see #writeVol(double[], int[], double[], int[], String)
	 * @see #writeVol(double[], int[], int, double[], int[], String)
	 * @see #writeVol(double[], int[], boolean, int[], double[], int[], String)
	 * 
	 */
	public static void writeVol(double[][][] vol3D, int datatype, double[] voxSize, String filename) throws IOException, FileNotFoundException {
		Nifti1Dataset NiftiDS = new Nifti1Dataset();
		NiftiDS.setHeaderFilename(filename);
		NiftiDS.setDataFilename(filename);
		if (!MLFuncs.contains(DATA_TYPES, datatype)) {
			throw new IllegalArgumentException("Invalid datatype: " + datatype);
		}
		NiftiDS.setDatatype((short)datatype);
		NiftiDS.setVoxSize(voxSize);
		NiftiDS.setDims((short)3, (short)vol3D[0][0].length, (short)vol3D[0].length, (short)vol3D.length, 
				(short)0, (short)0, (short)0, (short)0);
		NiftiDS.writeHeader();
		NiftiDS.writeVol(vol3D, (short)0);
	}
	
	/** Writes 3D volume into Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}.  
	 *  File is of <i>.nii</i> format iff <code>filename</code> contains <i>.nii</i> extension; 
	 *  otherwise data is saved into pair of files with <i>.img</i>/<i>.hdr</i> extensions. 
	 *  Voxel size is set to default 1mm<sup>3</sup>.  
	 * @param vol3D - data to be saved as doubles; data must be in [Z][Y][X] order
	 * @param datatype - see {@link #writeVol(double[][][], int, double[], String)} for details
	 * @param filename name of file to be saved
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 * @see #writeVol(double[][][], String)
	 * @see #writeVol(double[][][], double[], String)
	 * @see #writeVol(double[][][], int, double[], String)
	 * @see #writeVol(double[], int[], double[], int[], String)
	 * @see #writeVol(double[], int[], int, double[], int[], String)
	 * @see #writeVol(double[], int[], boolean, int[], double[], int[], String)
	 * 
	 */
	public static void writeVol(double[][][] vol3D, int datatype, String filename) throws IOException, 
		FileNotFoundException {
		// use default voxel size {1, 1, 1};
		double[] defaultVoxSz = {1, 1, 1};
		writeVol(vol3D, datatype, defaultVoxSz, filename);
	}

	/** Writes 3D volume into Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}.  
	 *  File is of <i>.nii</i> format iff filename contains <i>.nii</i> extension; otherwise
	 *  data is saved into pair of files with <i>.img</i>/<i>.hdr</i> extensions.
	 *  Data is saved as floats (32 bits / voxel).
	 *  Voxel size is set to default 1mm<sup>3</sup>.
	 *  
	 * @param vol3D - data to be saved; data must be in [Z][Y][X] order
	 * @param filename
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 * @see #writeVol(double[][][], int, String)
	 * @see #writeVol(double[][][], double[], String)
	 * @see #writeVol(double[][][], int, double[], String)
	 * @see #writeVol(double[], int[], double[], int[], String)
	 * @see #writeVol(double[], int[], int, double[], int[], String)
	 * @see #writeVol(double[], int[], boolean, int[], double[], int[], String)
	 */
	public static void writeVol(double[][][] vol3D, String filename) throws IOException, FileNotFoundException {
		writeVol(vol3D, 16, filename);
	}
	
	/** Writes 3D volume into Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}.  
	 *  File is of <i>.nii</i> format iff filename contains <i>.nii</i> extension; otherwise
	 *  data is saved into pair of files with <i>.img</i>/<i>.hdr</i> extensions.
	 *  Data is saved as floats (32 bits / voxel). 
	 *  
	 * @param vol3D - data to be saved; data must be in [Z][Y][X] order
	 * @param voxSize - 3-element double array containing voxel size of data {Xsize, Ysize, Zsize}
	 * @param filename
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 * @see #writeVol(double[][][], String)	 
	 * @see #writeVol(double[][][], int, String)
	 * @see #writeVol(double[][][], int, double[], String)
	 * @see #writeVol(double[], int[], double[], int[], String)
	 * @see #writeVol(double[], int[], int, double[], int[], String)
	 * @see #writeVol(double[], int[], boolean, int[], double[], int[], String)
	 * 
	 */
	public static void writeVol(double[][][] vol3D, double[] voxSize, String filename) throws IOException, 
		FileNotFoundException {
		writeVol(vol3D, 16, voxSize, filename);
	}
	
	
	/** Writes 1D volume into 3D Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}.  
	 *  File is of <i>.nii</i> format iff filename contains <i>.nii</i> extension; otherwise
	 *  data is saved into pair of files with <i>.img</i>/<i>.hdr</i> extensions.
	 *  
	 * @param vol1D - data to be saved; data must be in [Z][Y][X] order, i.e.,
	 *              data[x + y*xdim + z*xdim*ydim] = data3D[z][y][x]
	 * @param volDims - 3-element int array containing dimensions of volume to be saved
	 *                  {xdim, ydim, zdim}
	 * @param datatype -
	 *          	see {@link #writeVol(double[][][], int, double[], String)} for details
	
     * @param voxSize - 3-element double array containing voxel size of data {Xsize, Ysize, Zsize}
     * @param origin - 3-element int array containing location of origin in voxel coordinates
	 * @param filename
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 * @see #writeVol(double[][][], String)	 
	 * @see #writeVol(double[][][], int, String)
	 * @see #writeVol(double[][][], double[], String)
	 * @see #writeVol(double[][][], int, double[], String)
	 * @see #writeVol(double[], int[], double[], int[], String)
	 * @see #writeVol(double[], int[], boolean, int[], double[], int[], String)
	 */
	public static void writeVol(double[] vol1D, int[] volDims, int datatype, double[] voxSize,
			int[] origin, String filename) throws IOException, FileNotFoundException {
		Nifti1Dataset NiftiDS = new Nifti1Dataset();
		NiftiDS.setHeaderFilename(filename);
		NiftiDS.setDataFilename(filename);
		if (!MLFuncs.contains(DATA_TYPES, datatype)) {
			throw new IllegalArgumentException("Invalid datatype: " + datatype);
		}
		NiftiDS.setDatatype((short)datatype);
		NiftiDS.setVoxSize(voxSize);
		NiftiDS.setDims((short)3, (short)volDims[0], (short)volDims[1], (short)volDims[2],
				(short)0, (short)0, (short)0, (short)0);
		float[] voxOffset = new float[] {(float)-22, (float)31, (float)18};
//		float mmOriginX = (float)((Math.abs(voxOffset[0]) + 0.5) * (Math.signum(voxOffset[0])*voxSize[0]));
//		float mmOriginY = (float)((Math.abs(voxOffset[1]) + 0.5) * (Math.signum(voxOffset[1])*voxSize[1])); 
//		float mmOriginZ = (float)((Math.abs(voxOffset[2]) + 0.5) * (Math.signum(voxOffset[2])*voxSize[2])); 
		float mmOriginX = (float)((Math.abs(voxOffset[0]) + 0.5) * voxSize[0]);
		float mmOriginY = (float)((Math.abs(voxOffset[1]) + 0.5) * voxSize[1]); 
    	float mmOriginZ = (float)((Math.abs(voxOffset[2]) + 0.5) * voxSize[2]); 
		NiftiDS.qoffset = new float[] {-mmOriginX, -mmOriginY, -mmOriginZ};
		NiftiDS.quatern = new float[] {(float)0.0, (float)1.0, (float)0.0};
		NiftiDS.srow_x = new float[] {(float)voxSize[0], (float)0.0, (float)0.0, -mmOriginX};
		NiftiDS.srow_y = new float[] {(float)0.0, (float)voxSize[1], (float)0.0, -mmOriginY};
		NiftiDS.srow_z = new float[] {(float)0.0, (float)0.0, (float)voxSize[2], -mmOriginZ};
		NiftiDS.sform_code = 2;
		NiftiDS.qform_code = 4;
		
		NiftiDS.qfac = -1;
		
		NiftiDS.writeHeader();
		NiftiDS.printHeader();
		double[][][] vol3D = MLFuncs.reshapeZYX(vol1D, volDims[0], volDims[1], volDims[2]);
		NiftiDS.writeVol(vol3D, (short)0);
	}
	
	/** Writes 1D volume into 3D Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}.  
	 *  File is of <i>.nii</i> format iff filename contains <i>.nii</i> extension; otherwise
	 *  data is saved into pair of files with <i>.img</i>/<i>.hdr</i> extensions.
	 *  Data is saved as floats (32 bits / voxel)
	 *  
	 * @param vol1D - data to be saved; data must be in [Z][Y][X] order, i.e.,
	 *              data[x + y*xdim + *xdim*ydim] = data3D[z][y][x]
	 * @param volDims - 3-element int array containing dimensions of volume to be saved
	 *                  {xdim, ydim, zdim}
	 * @param voxSize - 3-element double array containing voxel size of data {Xsize, Ysize, Zsize}
	 * @param origin - 3-element int array containing location of origin in voxel coordinates
	 * @param filename
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 * @see #writeVol(double[][][], String)	 
	 * @see #writeVol(double[][][], int, String)
	 * @see #writeVol(double[][][], double[], String)
	 * @see #writeVol(double[][][], int, double[], String)
	 * @see #writeVol(double[], int[], int, double[], int[], String)
	 * @see #writeVol(double[], int[], boolean, int[], double[], int[], String)
	 */
	public static void writeVol(double[] vol1D, int[] volDims, double[] voxSize, int[] origin,
			String filename) throws IOException, FileNotFoundException {
		writeVol(vol1D, volDims, 16, voxSize, origin, filename);
	}
	
	/** Writes masked 1D volume into 3D Nifti file <code>filename</code> using {@link extern.niftijlib.Nifti1Dataset}. 
	 *  File is of <i>.nii</i> format iff filename contains <i>.nii</i> extension; otherwise
	 *  data is saved into pair of Nifti files with <i>.img</i>/<i>.hdr</i> extensions.
	 *  Data is saved as floats (32-bits/ voxel)
	 *  
	 *  @param vol1D - data to be saved; data is in same order as in maskInds as
	 *               each element in data is mapped to corresponding maskInds 
	 *               location in 3D volume
	 *  @param maskInds - indices in 3D volume of elements stored in vol1D 
	 *                  - the 3D volume will be in [Z][Y][X] order
	 * @param oneRel - true if maskInds are 1-relative; false if they are 0-relative
	 * @param volDims - 3-element int array containing dimensions of volume to be saved
	 *                  {xdim, ydim, zdim}
	 * @param voxSize - 3-element double array containing voxel size of data {Xsize, Ysize, Zsize}
	 * @param origin - 3-element int array containing location of origin in voxel coordinates
	 * @param filename
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 * @see #writeVol(double[][][], String)	 
	 * @see #writeVol(double[][][], int, String)
	 * @see #writeVol(double[][][], double[], String)
	 * @see #writeVol(double[][][], int, double[], String)
	 * @see #writeVol(double[], int[], double[], int[], String)
	 * @see #writeVol(double[], int[], int, double[], int[], String)

	 */
	public static void writeVol(double[] vol1D, int[] maskInds, boolean oneRel, int[] volDims, 
			double[] voxSize, int[] origin, String filename) throws IOException, FileNotFoundException {
		if (vol1D.length != maskInds.length) {
			throw new IllegalArgumentException("Data and mask index arrays must " +
					"be of same length.");
		}
		if (volDims.length != 3) {
			throw new IllegalArgumentException("Volume dimension info must have 3 elements.");
	
		}
		if (voxSize.length != 3) {
			throw new IllegalArgumentException("Voxel size info must have 3 elements.");
		}
		if (origin.length != 3) {
			throw new IllegalArgumentException("Origin info must have 3 elements.");
		}

		int nVox = volDims[0] * volDims[1] * volDims[2];
		
		double[] fullVol1D = new double[nVox]; // embed input data into 1D vol 
		if (oneRel) {
			for (int i = 0; i < maskInds.length; ++i) {
				maskInds[i]--;
			}
		}
		for (int i = 0; i < maskInds.length; ++i) {
			int currIdx = maskInds[i];
			if (currIdx < 0 || currIdx >= nVox) {
				throw new IllegalArgumentException("Invalid mask index element: " + 
						currIdx + ".  Largest valid value: " + (nVox - 1));
			}
			fullVol1D[currIdx] = vol1D[i];
		}
		writeVol(fullVol1D, volDims, voxSize, origin, filename);
	}
	
	/**Writes given data volumes into a single 4D Nifti file.
	 * 
	 * @param cvs <p> 2D array of data volumes (e.g., CVs in NPAIRS), [# volumes][# masked voxels]
	 * 			<p>- each volume is a 1D array of voxel values corresponding to 
	 * 				masked data locations given in <code>maskInds</code>
	 * @param maskInds <p> array of indices giving location of masked data within 
	 * 			3D data volume of <code>volDims</code> size
	 * 			<p>- assumption: 3D data is arranged in [Z][Y][X] order and
	 * 				1D voxel ordering is [x + y*xdim + z*ydim*xdim]
	 * @param oneRel <p> true if given mask indices are 1-relative
	 * @param volDims <p> int array of 3D volume dimensions: {xDim, yDim, zDim}
	 * @param voxSize <p> double array containing voxel size in mm: {xSize, ySize, zSize}
	 * @param filename <p> name of file to be saved
	 * @throws IOException
	 * @throws FileNotFoundException
	 */	
	public static void writeVol4DNpairs(double[][] cvs, int[] maskInds,
			boolean oneRel, int[] volDims, double[] voxSize, String filename)
			throws IOException, FileNotFoundException{
		
		if (cvs[0].length != maskInds.length) {
			throw new IllegalArgumentException("Data and mask index arrays must " +
					"be of same length.");
		}
		if (volDims.length != 3) {
			throw new IllegalArgumentException("Volume dimension info must have 3 elements.");
	
		}
		if (voxSize.length != 3) {
			throw new IllegalArgumentException("Voxel size info must have 3 elements.");
		}

		if (oneRel) { // adjust maskInds to be 0-relative
			for (int i = 0; i < maskInds.length; ++i) {
				maskInds[i]--;
			}
		}
		
		int nVox = volDims[0] * volDims[1] * volDims[2];
		// embed input data into 2D volume, one 1D array for each CV volume:
		double[][] fullVol2D = new double[cvs.length][nVox]; 
	
		for (int i = 0; i < maskInds.length; ++i) {
			int currIdx = maskInds[i];
			if (currIdx < 0 || currIdx >= nVox) {
				throw new IllegalArgumentException("Invalid mask index element: " + 
						currIdx + ".  Largest valid value: " + (nVox - 1));
			}
			for(int cv = 0; cv < cvs.length; cv++){
				fullVol2D[cv][currIdx] = cvs[cv][i]; 
			}
		}
		writeVol4DNiftiDataset(volDims, voxSize, cvs.length, filename, fullVol2D);
	}
	
	//input lags should be zero relative.
	
	/** Writes given data volumes into a single 4D Nifti file.
	 * 
	 * @param singleLV <p> array containing data for multiple 3D volumes (<code>winSize</code> # of them), 
	 * 				interspersed in PLS style such that the first <code>winSize</code> values correspond
	 * 				to the first voxel in each of <code>winSize</code> volumes, the next 
	 * 				<code>winSize</code> values correspond to the second voxel, etc., i.e., 
	 * 				<p> {V1<sub>vol_1</sub>, V1<sub>vol_2</sub>,...,V1<sub>vol_winSize</sub>,
	 * 				V2<sub>vol_1</sub>,V2<sub>vol_2</sub>,...,V2<sub>vol_<code>winSize</code></sub>,...}
	 * @param maskInds <p>array of indices giving location of masked data within 
	 * 			3D data volume of <code>volDims</code> size
	 * 			<p>- assumption: 3D data is arranged in [Z][Y][X] order and
	 * 			1D voxel ordering is [x + y*xdim + z*ydim*xdim]
	 * @param oneRel <p> true if given mask indices are 1-relative
	 * @param winSize <p> number of volumes interspersed in <code>singleLV</code>
	 * @param lags <p> 0-relative array of 'lags' to include, i.e., which of the volumes
	 * 			contained in <code>singleLV</code> to save,
	 * 			<ul>e.g., 
	 * 			<li> <code>lags</code> = {0, 1,..., <code>winSize</code> - 1} means include all volumes;
	 * 			<li> <code>lags</code> = {2, 5} means include only the 3rd and 6th volumes
	 * 			</ul>
	 * @param volDims <p> int array of 3D volume dimensions: {xDim, yDim, zDim}
	 * @param voxSize <p> double array containing voxel size in mm: {xSize, ySize, zSize}
	 * @param filename <p> name of file to be saved
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 */
	public static void writeVol4DPLS(double[] singleLV, int[] maskInds, boolean oneRel, 
			int winSize, int[] lags, int[] volDims, double[] voxSize, String filename) 
			throws IOException, FileNotFoundException{
		
		if(singleLV.length / maskInds.length != winSize){
			throw new IllegalArgumentException("Data and mask index arrays must " +
			"be such that len(data) / len(mask) = number of lags");
		}

		if(lags.length == 0){
			throw new IllegalArgumentException("No lags specified to extract.");
		}
		
		if (volDims.length != 3) {
			throw new IllegalArgumentException("Volume dimension info must have 3 elements.");
	
		}
		if (voxSize.length != 3) {
			throw new IllegalArgumentException("Voxel size info must have 3 elements.");
		}
		
		int previous = -1;
		for(int lag : lags){
			if(previous < lag) previous = lag;
			else throw new IllegalArgumentException("Lags must be in ascending order");
		}
		
		if (oneRel) {
			for (int i = 0; i < maskInds.length; ++i) {
				maskInds[i]--;
			}
		}
		
		int nVox = volDims[0] * volDims[1] * volDims[2];
		double[][] lvMultiLags = new double[lags.length][nVox];
		int currentVoxel;
		int currIdx;
		
		for(int voxel = 0; voxel < singleLV.length; voxel+=winSize){
			currentVoxel = voxel / winSize;
			currIdx = maskInds[currentVoxel];
			
			if(currIdx < 0 || currIdx >= nVox){
				throw new IllegalArgumentException(
						"Invalid mask index element: " + 
						currIdx + ".  Largest valid value: " + (nVox - 1));
			}
			
			for(int i = 0; i < lags.length; i++){
				lvMultiLags[i][currIdx] =
					singleLV[voxel + lags[i]];
			}
		}
		
		writeVol4DNiftiDataset(volDims, voxSize, lags.length, filename,
				lvMultiLags);
	}

	/**
	 * @param volDims the voxel dimensions
	 * @param voxSize voxel size
	 * @param units the number of cvs / lags depending on the datatype.
	 * @param filename the filename to write out this Nifti information.
	 * @param vol2D either [lag#][lag] or [cv#][cv]
	 * @throws IOException 
	 * @throws FileNotFoundException
	 */
	private static void writeVol4DNiftiDataset(int[] volDims, double[] voxSize,
			int units, String filename, double[][] vol2D)
			throws IOException, FileNotFoundException {
		int datatype = 16;
		
		Nifti1Dataset NiftiDS = new Nifti1Dataset();
		NiftiDS.setHeaderFilename(filename);
		NiftiDS.setDataFilename(filename);
		if (!MLFuncs.contains(DATA_TYPES, datatype)) {
			throw new IllegalArgumentException("Invalid datatype: " + datatype);
		}
		NiftiDS.setDatatype((short)datatype);
		NiftiDS.setVoxSize(voxSize);
		NiftiDS.setDims((short)4, (short)volDims[0], (short)volDims[1], 
				(short)volDims[2],(short)units, (short)0, 
				(short)0, (short)0);
		
		NiftiDS.writeHeader();
	 
		for(int unit = 0; unit < units; unit++ ){
			double[][][] vol3D = MLFuncs.reshapeZYX(vol2D[unit], 
									volDims[0], volDims[1], volDims[2]);
			NiftiDS.writeVol(vol3D, (short)unit);
		}
	}
	
	/** Returns X, Y, Z dimensions of input Nifti volume.
	 * 
	 * @param niftiFilename Nifti data volume file name. Volume information is read from 
	 * information listed in Nifti file header for X, Y, Z dimensions. 
	 * 
	 * @return 3-element int array {xdim, ydim, zdim}
	 */
	public static int[] getVolDims3D(String niftiFilename) throws IOException {
		int[] volDims3D = null;		
		Nifti1Dataset NiftiDataset = new Nifti1Dataset(niftiFilename);
		// get dims
		try {
			NiftiDataset.readHeader();
		}
		catch (FileNotFoundException e) {
			throw new IOException(e.getMessage());
		}
		volDims3D = new int[3];
		volDims3D[0] = NiftiDataset.getXdim();
		volDims3D[1] = NiftiDataset.getYdim();
		volDims3D[2] = NiftiDataset.getZdim();
		
		return volDims3D;
	}

	/** Reads input data into a 2D array (returned) containing a 1D representation 
	 *  of data for each included timepoint t. 
	 *  
	 *  <p> Assumption: input data 3D volumes are arranged in [Z][Y][X] order, 
	 *  hence 2D array is arranged in [t][x + y*xdim + z*xdim*ydim] order. 
	 *  
	 *  <p> Note that input data file may be 3D or 4D.
	 *  
	 * @param NiftiFilename
	 * 		name of Nifti or ANALYZE file containing data
	 * @param skipTmpts 
	 * 		indices of timepoints to skip (0-relative)
	 *      <p> Note that <code>skipTmpts</code> can contain more points than the # of timepoints 
	 *      in input file, as the timepoints in the current input file may be only a subset of a larger 
	 *       collection of skipped timepoints. See <code>firstTmpt</code>.  
	 * @param firstTmpt which timepoint (in <code>skipTmpts</code>) corresponds to first timepoint of 
	 * 		current file. 
	 *  	E.g., if the full set of timepoints is contained in 3 input data files (file1, file2, file3) 
	 *  	and each file contains 10 timepoints, then to skip the first 3 timepoints in each file, 
	 *  	<p> <code>skipTmpts</code> = {0, 1, 2, 10, 11, 12, 20, 21, 22}.  When reading in file2 
	 *  	data, set <code>firstTmpt</code> = 10.
	 *  
	 * @return 2D double array containing data [# timepoints][# voxels]
	 * 			<p> if all timepoints in input file are skipped, return <code>null</code>
	 * @throws IOException
	 */
	public static double[][] readNiftiData(String NiftiFilename, int[] skipTmpts, int firstTmpt) 
		throws IOException {
		
		Nifti1Dataset NiftiDataset = new Nifti1Dataset(NiftiFilename);

		// get dims
		try {
			NiftiDataset.readHeader();
		}
		catch (FileNotFoundException e) {
			throw new IOException(e.getMessage());
		}

		int xDim = NiftiDataset.getXdim();
		int yDim = NiftiDataset.getYdim();
		int zDim = NiftiDataset.getZdim();
		int tDim = NiftiDataset.getTdim();

		if (debug) {
			System.out.println("Xdim: " + xDim);
			System.out.println("Ydim: " + yDim);
			System.out.println("Zdim: " + zDim);
			System.out.println("Tdim: " + tDim);
		}

		if (tDim == 0) tDim = 1;

		// load data 
		long xyz = (long)xDim * (long)yDim * (long)zDim;
		//System.out.println("xDim * yDim * zDim: " + xyz);
		
		
		int[] currSkipT = getCurrSkipT(firstTmpt, tDim, skipTmpts);
		int nSkip = currSkipT.length;
		if (debug) {
			System.out.println("Skipping " + nSkip + " tmpts...");
		}
		double[][] dataArray2D = null;
		if ((tDim - nSkip) > 0) {
			// not all tmpts are skipped
			dataArray2D = new double[tDim - nSkip][(int)xyz];
			int i = 0;
			for (short t = 0; t < tDim; ++t) {
				if (MLFuncs.contains(currSkipT, t)) {
					if (debug) {
						System.out.println("Skipping tmpt " + t);
					}
				}
				else {
					if (debug) {
						System.out.println("Reading tmpt " + t);
					}
					dataArray2D[i] = NiftiDataset.readDoubleVol1D(t);
					++i;
				}

//				double[] currVol = NiftiDataset.readDoubleVolReturnArray(t);
				// switch from row (x)  to col (y) dominant
//				for (int z = 0; z < zDim; ++z) {
//				for (int y = 0; y < yDim; ++y) {
//				for (int x = 0; x < xDim; ++x) {
//				dataArray2D[t][y + (yDim * x) + (xDim * yDim * z)] = 
//				currVol[x + (xDim * y) + (xDim * yDim *z)];
//				}
//				}
//				}

			}
		}
		return dataArray2D;
	}
	
	/** Reads input Nifti or ANALYZE data file into 2D byte array (returned). 
	 * <p> Doesn't check datatype!
	 * <p> 2D array is [# timepoints][# voxels] and contains a 1D representation 
	 * of a 3D image for each timepoint. 
	 * <p> Assuming that input data 3D volumes are arranged in [Z][Y][X] order,
	 *  each data in output 2D array is arranged in [t][x + y*xdim + z*xdim*ydim] order.
	 * <p> Note that input datafile can be 3D or 4D.
	 * 
	 * @param NiftiFilename
	 * 		name of file containing data
	 * @param skipTmpts 
	 * 		indices of timepoints to skip (0-relative). 
	 * See #readNiftiData(String NiftiFilename, int[] skipTmpts, int firstTmpt) for details.
	 * 
	 * @param firstTmpt which timepoint (in skipTmpts) corresponds to first timepoint of current file.  
	 * See #readNiftiData(String NiftiFilename, int[] skipTmpts, int firstTmpt) for details.
	 * 
	 * @return 2D byte array containing data [# timepoints][# voxels]; if all timepoints in input
	 * 	       file are skipped, return null
	 * @throws IOException
	 */
	public static byte[][] readNiftiDataBytes(String NiftiFilename, int[] skipTmpts, int firstTmpt) 
		throws IOException {
		
		Nifti1Dataset NiftiDataset = new Nifti1Dataset(NiftiFilename);

		// get dims
		try {
			NiftiDataset.readHeader();
		}
		catch (FileNotFoundException e) {
			throw new IOException(e.getMessage());
		}

		int xDim = NiftiDataset.getXdim();
		int yDim = NiftiDataset.getYdim();
		int zDim = NiftiDataset.getZdim();
		int tDim = NiftiDataset.getTdim();

		if (debug) {
			System.out.println("Xdim: " + xDim);
			System.out.println("Ydim: " + yDim);
			System.out.println("Zdim: " + zDim);
			System.out.println("Tdim: " + tDim);
		}

		if (tDim == 0) tDim = 1;

		// load data 
		long xyz = (long)xDim * (long)yDim * (long)zDim;
		//System.out.println("xDim * yDim * zDim: " + xyz);
		
		int[] currSkipT = getCurrSkipT(firstTmpt, tDim, skipTmpts);
		int nSkip = currSkipT.length;
		if (debug) {
			System.out.println("Skipping " + nSkip + " tmpts...");
		}
		byte[][] dataArray2D = null;
		if ((tDim - nSkip) > 0) {
			// not all tmpts are skipped
			dataArray2D = new byte[tDim - nSkip][(int)xyz];
			int i = 0;
			for (short t = 0; t < tDim; ++t) {
				if (MLFuncs.contains(currSkipT, t)) {
					if (debug) {
						System.out.println("Skipping tmpt " + t);
					}
				}
				else {
					if (debug) {
						System.out.println("Reading tmpt " + t);
					}
					dataArray2D[i] = NiftiDataset.readVolBlob(t);
					++i;
				}

//				double[] currVol = NiftiDataset.readDoubleVolReturnArray(t);
				// switch from row (x)  to col (y) dominant
//				for (int z = 0; z < zDim; ++z) {
//				for (int y = 0; y < yDim; ++y) {
//				for (int x = 0; x < xDim; ++x) {
//				dataArray2D[t][y + (yDim * x) + (xDim * yDim * z)] = 
//				currVol[x + (xDim * y) + (xDim * yDim *z)];
//				}
//				}
//				}

			}
		}
		return dataArray2D;
	}
	
	
	/** Returns Matrix containing masked data.  Rows are timepoints (scans); columns
	 *  are voxels.  Voxels <= 0 are NOT excluded.
	 * @param maskIndices indices of voxels to be included (if null, calculate 
	 *  using input data) - NOT checked for validity!  
	 *  
	 *  @return Matrix of masked data
	 */
	protected static Matrix getMaskedDataMat(NpairsSetupParams setupParams, String matlibType, int[] maskIndices) 
		throws MatrixException, IOException {
		
		int[] maskLocs = null;
		if (maskIndices == null) {
			// figure out mask indices from input data 
			double[] andMask = getANDMask(setupParams.getMaskFilenames());
			maskLocs = MLFuncs.findNonZero(andMask);
		}
		else maskLocs = maskIndices;
		
		int nMskVox = maskLocs.length;
		
		int nRowsIncl = setupParams.getNumVols();
		Matrix maskedData = new MatrixImpl(nRowsIncl, nMskVox, matlibType).getMatrix();
		
		// add data to Matrix one tmpt (row) at a time
		int row = 0; 
		int nFiles = setupParams.getNumDataFiles();
		int firstTmpt = 0; // first tmpt in each file
		for (int f = 0; f < nFiles; ++f) {
			int currNTmpts = setupParams.getNTmptsPerFile()[f];
			int[] currSkipTmpts = getCurrSkipT(firstTmpt, currNTmpts, 
					setupParams.getSkipTmpts());
			firstTmpt += currNTmpts;	
			if (currSkipTmpts.length < currNTmpts) {
				// not all tmpts skipped in curr file
//				System.out.println("Reading " + f + "th file: " + setupParams.getDataFilenames()[f]);
				Nifti1Dataset data = new Nifti1Dataset(setupParams.getDataFilenames()[f]);
				try {
					data.readHeader();
				}
				catch (FileNotFoundException e) {
					throw new IOException(e.getMessage());
				}
				
				for (int t = 0; t < currNTmpts; ++t) {
					if (!MLFuncs.contains(currSkipTmpts, t)) {
						double[] currData = data.readDoubleVol1D((short)t);
						double[] currMaskedData = new double[nMskVox];
						for (int i = 0; i < nMskVox; ++i) {	
							currMaskedData[i] = currData[maskLocs[i]];
						}
						maskedData.setRowQuick(row, currMaskedData);
						++row;
					}
				}
			}
//			else {
//				System.out.println("Skipping " + f + "th file: " + setupParams.getDataFilenames()[f]);
//			}
		}		
		return maskedData;
	}
		
		
	/** Returns array containing indices of subset of skipTmpts within input range
	 *  of tDim tmpts beginning with firstTmpt: [firstTmpt, firstTmpt + tDim - 1]
	 * 
	 * @param firstTmpt
	 * @param tDim
	 * @param skipTmpts
	 * @return
	 */
	private static int[] getCurrSkipT(int firstTmpt, int tDim, int[] skipTmpts) {
		int[] uniqSkipT = new int[0];
		if (skipTmpts != null) {
			uniqSkipT = MLFuncs.sortAscending(MLFuncs.unique(skipTmpts));
		}
		if (debug) {
			System.out.println("First tmpt: " + firstTmpt);
		}
		// find indices of skipped tmpts greater than or equal to lower bound (firstTmpt)
		int[] skipIdxBoundBelow = MLFuncs.findGreaterThanOrEqualTo(uniqSkipT, firstTmpt);
//		 find indices of skipped tmpts less than or equal to upper bound (lastTmpt)
		int lastTmpt = firstTmpt + tDim - 1;
		int[] skipIdxBoundAbove = MLFuncs.findLessThanOrEqualTo(uniqSkipT, lastTmpt);
		// find intersection of bounded tmpts to get current subset of skipped tmpts
		// (i.e., tmpt indices)
		Vector<Integer> vCurrSkipT = new Vector<Integer>();
		for (int i = 0; i < skipIdxBoundBelow.length; ++i) {
			if (MLFuncs.contains(skipIdxBoundAbove, skipIdxBoundBelow[i])) {
				vCurrSkipT.add(uniqSkipT[skipIdxBoundBelow[i]]);
			}
		}
		int[] currSkipT = new int[vCurrSkipT.size()];
		for (int i = 0; i < currSkipT.length; ++i) {
			currSkipT[i] = vCurrSkipT.get(i) - firstTmpt;
		}
		return currSkipT;
	}

	/** Reads input Nifti or ANALYZE data files into 2D array (returned). 
	 * <p> Returned array is [# timepoints][# voxels] and contains a 1D representation 
	 * of a 3D image for each timepoint. 
	 * <p> Assuming that input 3D data is arranged in [Z][Y][X] order,
	 *  each row of data in output 2D array is arranged in
	 *  [x + y*xdim + z*xdim*ydim] order.
	 *  
	 * <p> Note that input data can be 3D or 4D.
	 * 
	 * <p> Data in output matrix is stacked in order of input NiftiFilenames.
	 * <p> Required: all input 3D volumes have same number of voxels.
	 * 
	 * @param NiftiFilenames
	 * 			 String array of data filenames
	 * @return 2D double array containing data [# timepoints][# voxels]
	 */
	public static double[][] readNiftiData(String[] NiftiFilenames) throws IOException {

		double[][] dataArray2D = null;

		int numFiles = NiftiFilenames.length;
		if (numFiles == 0) {
			throw new IOException("Error - no NiftiFilenames in input arg");
		}

		int numTmpts = 0;
		int numVoxels = 0;
		Vector<double[]> data = new Vector<double[]>();

		for (int i = 0; i < numFiles; ++i) {
			if (debug) {
				System.out.println("Reading " + NiftiFilenames[i] + "...");
			}
			int[] skipTmpts = null;
			double[][] currData2D = readNiftiData(NiftiFilenames[i], skipTmpts, 0);
			numTmpts += currData2D.length;
			if (i == 0) {
				numVoxels = currData2D[0].length;
			}
			else if (numVoxels != currData2D[0].length) {
				throw new IOException("Error loading Nifti data: " + NiftiFilenames[i] + 
				" - no. of voxels in dataset does not match previously loaded datasets.");
			}
			for (int t = 0; t < currData2D.length; ++t) {
				data.add(currData2D[t]);
			}
		}
		if (debug) {
			System.out.println("Total number of tmpts: " + numTmpts);
			System.out.println("Number of voxels: " + numVoxels);
		}
		dataArray2D = new double[numTmpts][numVoxels];
		dataArray2D = data.toArray(dataArray2D);		

//		if (debug) {
//			System.out.println("Array to be returned by readData(String[]): ");
//			utils_tests.PCATest.printArray(dataArray2D);
//		}
		return dataArray2D;
	}
	
	
	/** Reads input Nifti or ANALYZE data files into 2D array (returned).
	 *  Input data files are masked with input masks, i.e.,
	 *  voxels set to 0 in input masks are also set to 0
	 *  in corresponding data.
	 *  Masks are combined into AND mask (see #getANDMask(String[])) which is then 
	 *  applied to all input 3D data volumes. Data is read in order presented in
	 *  NiftiFilenames.  
	 *  <p> Returned array is [# timepoints][# voxels] and contains a 1D representation 
	 *  of a 3D image for each timepoint. 
	 *  <p>Assuming that input data is arranged in [Z][Y][X] order,
	 *  each row of data in output 2D array is arranged in
	 *  [x + y*xdim + z*xdim*ydim] order.
	 *  
	 * @param NiftiFilenames 
	 * 			String array of Nifti or ANALYZE data filenames. Files may be 3D or 4D.
	 * @param maskFilenames
	 * 			String array of Nifti or ANALYZE mask filenames. Masks must be 3D.
	 * 
	 * @return 2D double array containing data [# timepoints][# voxels] 
	 */
	public static double[][] readNiftiData(String[] NiftiFilenames, String[] maskFilenames) 
			throws IOException {
		//TODO: save mask indices somewhere else, not in NiftiIO code.
		return readNiftiData(NiftiFilenames, maskFilenames, null, null, null);
	}

	/** Reads input Nifti or ANALYZE data files into 2D array (returned).
	 *  <p>Input data files are masked with input masks, i.e.
	 *  voxels set to 0 in input masks are also set to 0
	 *  in corresponding data.
	 *  <p>Masks are combined into ANDMask {see #getANDMask(String[])) which is then 
	 *  applied to all input 3D data volumes. Data is read in order presented in
	 *  NiftiFilenames, skipping the timepoints at indices in 'skipTmpts'.  
	 * <p> Returned array is [# timepoints][# voxels] and contains a 1D representation 
	 *  of a 3D image for each timepoint. 
	 *  <p> Assuming that input data is arranged in [Z][Y][X] order,
	 *  each row of data in output 2D array is arranged in
	 *  [x + y*xdim + z*xdim*ydim] order.
  	 *
	 * @param NiftiFilenames 
	 * 			String array of Nifti or ANALYZE data filenames. Data may be 3D or 4D.
	 * @param maskFilenames
	 * 			String array of Nifti or ANALYZE mask filenames. Masks must be 3D.
	 * @param saveMaskIndsFilename
	 *          String indicating filename under which to save indices of voxels included
	 *          in AND mask computed from input masks 
	 *          in an ASCII file (formatted to be read by idl read_matrix.pro)
	 *          - if null, then mask indices will not be saved
	 * @param skipTmpts 
	 * 			Timepoints to exclude when reading in data.  The indices 
	 *          run in order from 1st timepoint in NiftiFilenames[0] to
	 *          last timepoint in NiftiFilenames[n-1], where n = number of 
	 *          files in NiftiFilenames.  
	 *          - if null, include all timepoints
	 * @param nTmptsPerFile how many timepoints exist in each input file (in order)
	 * @return 2D double array containing data (timepoints X voxels) 
	 */
	public static double[][] readNiftiData(String[] NiftiFilenames, String[] maskFilenames, 
			String saveMaskIndsFilename, int[] skipTmpts, int[] nTmptsPerFile) throws IOException {
		int numDataFiles = NiftiFilenames.length;
		int numMaskFiles = maskFilenames.length;

		if (numDataFiles == 0) {
			throw new IllegalArgumentException("Error loading Nifti data-- input" +
					" NiftiFilenames arg is empty");
		}
		if (numMaskFiles == 0) {
			throw new IllegalArgumentException("Error loading mask files-- input" +
					" maskFilenames arg is empty");
		}
		
		double[][] maskedData = null; 
	//	ArrayList<Double[]> mData;
		
		byte[] andMask = getANDMaskBytes(maskFilenames);

		// get size and indices of mask
		int nVox = andMask.length;
		Vector<Integer> vMaskIndices = new Vector<Integer>();
		
		for (int i = 0; i < nVox; ++i) {
			if (andMask[i] != 0) {
				vMaskIndices.add(i);
			}
		}
		int nMskVox = vMaskIndices.size();	
		
		if (saveMaskIndsFilename != null) {
			int[] maskIndices = new int[nMskVox];
	
			for (int j = 0; j < nMskVox; ++j) {
				maskIndices[j] = vMaskIndices.get(j).intValue();		
			}
			NpairsIO.printToIDLFile(maskIndices, saveMaskIndsFilename);
		}

		//maskedData = new double[numDataFiles][nMskVox];
		int nScans = MLFuncs.sum(nTmptsPerFile); 
		int nSkipped = skipTmpts.length;
		
	//	mData = new ArrayList<Double[]>(nScans - nSkipped); 
		maskedData = new double[nScans - nSkipped][nMskVox];
		if (debug) {
			System.out.println("Size maskedData matrix: " + maskedData.length + " X " + 
					maskedData[0].length);
		}
		int currTmpt = 0;
		int currInclTmpt = 0;
		for (int i = 0; i < numDataFiles; ++i) {
			// load current dataset
			if (debug) {
				System.out.println("Reading datafile " + NiftiFilenames[i] + "... ");		
				System.out.println("Current tmpt: " + currTmpt);
			}
			double[][] currData = readNiftiData(NiftiFilenames[i], skipTmpts, currTmpt);
			if (currData == null) {
				// no tmpts included from current data file
				if (debug) {
					System.out.println("No tmpts incl from " + NiftiFilenames[i]);
				}
				currTmpt += nTmptsPerFile[i];
				continue;
			}
			else {
				if (debug) {
					System.out.println("Incl tmpts from " + NiftiFilenames[i]);
				}
				if (currData[0].length != nVox) {
					throw new IOException("Error-- data and mask images must have same "
							+ "dimensions");
				}
				
				int currNumScans = currData.length;
//				numScans += currNumScans - 1;
//				mData.ensureCapacity(numScans); // is this worth it?
				for (int k = 0; k < currNumScans; ++k) {
					int currRow = currInclTmpt + k;
//					Double[] currMaskedData = new Double[nMskVox];
					int maskedIndex = 0;
					for (int j = 0; j < nVox; ++j) {
						if (andMask[j] != 0) {
							maskedData[currRow][maskedIndex] = currData[k][j];
							++maskedIndex;
						}
					}
//					mData.add(currMaskedData);
				}
				currInclTmpt += currNumScans; // masked data row index incremented
			}
			currTmpt += nTmptsPerFile[i]; // total input data (incl and excl tmpts) index incremented
		}
//		maskedData = new double[nScans - nSkipped][nMskVox];
//		for (int n = 0; n < nScans - nSkipped; ++n) {
//			for (int v = 0; v < nMskVox; ++v) {
//				maskedData[n][v] = mData.get(n)[v];
//			}
//		}	
	
		if (maskedData == null) {
			throw new IOException("Could not create masked data");
		}
		return maskedData;
	}
	
	/** Returns data at given timepoint in input data file.
	 * 
	 * @param NiftiFilename
	 * 			name of file containing data in Nifti or ANALYZE format 
	 * @param tmpt
	 * 			timepoint to read
	 * @throws IllegalArgumentException
	 *          if timepoint is out of range
	 *          
	 * @return 1D double array containing input data at input timepoint
	 */
	public static double[] readNiftiData(String NiftiFilename, int tmpt) throws IOException {
		
		Nifti1Dataset NiftiDataset = new Nifti1Dataset(NiftiFilename);
		// get dims
		try {
			NiftiDataset.readHeader();
		}
		catch (FileNotFoundException e) {
			throw new IOException(e.getMessage());
		}
		
		int tDim = NiftiDataset.getTdim();
		if (tDim == 0) tDim = 1;
			
		if (tmpt >= tDim) {
			throw new IllegalArgumentException("Timepoint " + tmpt + " is out of range - " +
					"input file " + NiftiFilename + " has only " + tDim + " timepoints.");
		}		
		
		double[] data = NiftiDataset.readDoubleVol1D((short)tmpt);
		return data;
	}
	
	/** Returns masked data at given timepoint in input data file.
	 * 
	 * @param NiftiFilename
	 * 			name of file containing data in Nifti or ANALYZE format
	 * @param mask
	 * 			array containing 0's where data voxels are to be excluded
	 *  		- size of mask array must be same as 3d volume of data in input
	 *          Nifti file
	 * @param tmpt
	 * 			timepoint to read
	 * @throws IllegalArgumentException
	 *          if timepoint is out of range or mask is incorrect size
	 *          
	 * @return 1D double array containing masked input data at input timepoint
	 */
public static double[] readNiftiData(String NiftiFilename, double[] mask, int tmpt) 
		throws IOException {
		
		Nifti1Dataset NiftiDataset = new Nifti1Dataset(NiftiFilename);
		// get dims
		try {
			NiftiDataset.readHeader();
		}
		catch (FileNotFoundException e) {
			throw new IOException(e.getMessage());
		}
		int xDim = NiftiDataset.getXdim();
		int yDim = NiftiDataset.getYdim();
		int zDim = NiftiDataset.getZdim();
		int tDim = NiftiDataset.getTdim();
		if (tDim == 0) tDim = 1;
			
		if (tmpt >= tDim) {
			throw new IllegalArgumentException("Timepoint " + tmpt + " is out of range - " +
					"input file " + NiftiFilename + " has only " + tDim + " timepoints.");
		}
		long nVox = (long)xDim * (long)yDim * (long)zDim;		
		if (mask.length != nVox) {
			throw new IllegalArgumentException("Input mask and data are not the same size.");
		}
		
		double[] data = NiftiDataset.readDoubleVol1D((short)tmpt);
		int nMskVox = MLFuncs.findNonZero(mask).length;
		double[] maskedData = new double[nMskVox];
		int maskedIndex = 0;
		for (int j = 0; j < nVox; ++j) {
			if (mask[j] != 0) {
				maskedData[maskedIndex] = data[j];
				++maskedIndex;
			}
		}
		
		return maskedData;
	}

// Created to test Procrustes
/** Returns 1D array of masked data.
 * 
 * @param NiftiFilename name of data file in Nifti or ANALYZE format.
 * @param timepoint which timepoint in data file to read (put 0 if it's a 3D file)
 * @param maskCoords int array of 0-relative coordinates to be included in masked data
 * 					<p>- assuming 3D data is in [Z][Y][X] format, 1D coordinate
 * 						ordering is [x + y*xdim + z*ydim*xdim]
 * 
 * @return 1D double array containing masked data 
 */
public static double[] readNiftiData(String NiftiFilename, int tmpt, int[] maskCoords)
	throws IOException {
	Nifti1Dataset NiftiDataset = new Nifti1Dataset(NiftiFilename);
	try {
		NiftiDataset.readHeader();
	}
	catch (FileNotFoundException e) {
		throw new IOException(e.getMessage());
	}
	double[] data = NiftiDataset.readDoubleVol1D((short)tmpt);

	double[] maskedData = new double[maskCoords.length];
	int maskedIndex = 0;
	for (int j = 0; j < maskCoords.length; ++j) {
			maskedData[maskedIndex] = data[maskCoords[j] - 1];
			++maskedIndex;
	}
	return maskedData;
}

	
	
	/** Returns double 1D array containing AND mask, i.e., intersection of all 
	 * non-zero voxels, created from masks specified for each session file in input 
	 * setup parameter information.
	 * <p> Masks are created for each session file in one of 3 ways:
	 * <ol>
	 * <li>using input 3D mask volume file containing zeroes where voxels are to
	 * be excluded from mask
	 * <li>keeping all voxels (i.e. no data masked out)
	 * <li> via thresholding technique a la PLS, i.e.,
	 * 		data from each input scan is masked by applying threshold
	 * 	    = T * scan maximum voxel value (where T is input threshold 
	 *      value between 0 and 1); intersection of all scan masks
	 *      is calculated to produce AND mask across all scans
	 * </ol>
	 * <p> After a mask has been calculated for each session file via
	 * one of the 3 techniques listed above, this method takes the intersection across
	 * session file masks to produce a grand AND mask.
	 *   
	 * <p> All data and masks must have the same number of voxels.
	 * <p> Voxels to be included in mask are conventionally labelled with 1s, but there
	 * is no requirement that non-zero voxels in mask files be set to 1, and there is no 
	 * guarantee that non-zero voxels in returned AND mask will be set to 1 unless input
	 * mask files follow this convention.
	 * 
	 * @return AND mask created from masks specified in input setup parameter info, as
	 * 	double array (length number of voxels in mask and data volumes)
	 */
	private static double[] getANDMask(NpairsSetupParams nsp) throws NpairsException,
			IOException {
		String[] maskFilenames = nsp.getMaskFilenames();
		boolean[] useAllVox = nsp.inclAllVoxels();
		double[] threshVals = nsp.getMaskThreshVals();
		int nSess = maskFilenames.length;
		if (useAllVox.length != nSess
				||
			threshVals.length != nSess) {
			// Should never happen! If it does, there's a bug.
			throw new NpairsException("Error: mask, incl. voxel and threshold info arrays " +
						"should all be the same length!");
		}
		
		double[] andMask = null;
		int maskSize = 0;
		for (int i = 0; i < nSess; ++i) {
			if (!useAllVox[i]) {
				if (maskFilenames[i] != null) {
					double[] currMask = getANDMask(new String[] {maskFilenames[i]});
					if (andMask == null) {
						andMask = currMask;
						maskSize = currMask.length;
					}
					else {
						if (currMask.length != maskSize) {
							throw new IllegalArgumentException("All data and masks must have the " +
									"same voxel dimensions.");					
						}
						// zero out voxels in andMask that are zero in current mask
						for (int j = 0; j < maskSize; ++j) {
							if (currMask[0] == 0) {
								andMask[0] = 0;
							}
						}
					}
						
				}
				else if (threshVals[i] != -1) {
					// no mask file supplied for curr sess file; calculate mask 
					// using threshold value
//					double[] currMask = getThreshMask(nsp, threshVals[i],
//							false);
					
				}
				else {
					// no mask or thresh info supplied even though useAllVox = false
					throw new NpairsException("Error: data masking info for session file "
							+ "no. " + (nSess+1) + " is contradictory!");
				}
					
			}
			else {
				// using all voxels
			}
		}
			
		return null;
	}
	
	
	/** Returns mask generated by applying PLS-style thresholding to each
	 *  scan in current session and taking AND mask (intersection) of masks
	 *  generated for each scan.
	 *  PLS-style thresholding: multiply maximum voxel value in current scan
	 *  by input threshold to get lower bound on value of voxels to be included
	 *  in mask.  
	 * 
	 * @param NpairsSetupParams nsp - .mat file containing npairs setup parameter info
	 * @param double thresh - value to multiply scan maximum by to get calculated
	 *                threshold for each scan
	 * @param boolean considerAllVoxels - if false, only voxels strictly greater than
	 *                calculated thresh will be included in mask
	 *                                  - if true, voxels equal to calculated thresh 
	 *                will also be included 
	 * @see #pls.sessionprofile.RunGenerateDatamat.findOnsetCoords(Matrix, double, boolean)
	 * @return
	 */
	private static double[] getThreshMask(NpairsSetupParams nsp, int sessNum, double thresh, 
			boolean considerAllVoxels) {
		// which data files belong to this session?
		String[] sessDataFilenames = getSessDataFilenames(nsp, sessNum);
						
		return null;
	}

	
	/** Returns array (length total number of tmpts in given session file
	 * of data filenames by timepoint; if a given file is 4D
	 *  and contains multiple timepoints, said filename will be included once for each
	 *  timepoint it contains.
	 * @param nsp
	 * @param sessNum
	 * @return array of data filenames by timepoint
	 */
	private static String[] getSessDataFilenames(NpairsSetupParams nsp,
			int sessNum) {
		// if given vol not included,
	    // corresponding element in paddedSessLabels will be '-1'
		int[] paddedSessLabels = getPaddedSessLabels(nsp);
		int[] currSessTmpts = MLFuncs.find(paddedSessLabels, sessNum); // 0-relative tmpts
		Set<String> currSessDataFilenames = new HashSet<String>();
		int nDataFiles = nsp.getNumDataFiles();
		int currTmpt = 0;
		for (int i = 0; i < nDataFiles; ++i) {
			for (int j = 0; j < nsp.getNTmptsPerFile()[i]; ++j) {
				if (MLFuncs.contains(currSessTmpts, currTmpt)) {
					
				}
				++currTmpt;
			}
		}
		
		
		
		return null;
	}

	/** Returns session label array spanning all session
     file data volumes, not just the ones included
     in current analysis; if given vol not included,
     corresponding element in paddedSessLabels will be '-1'
     @param npairs setup parameters .mat file
     @return padded session label array
     */
	private static int[] getPaddedSessLabels(NpairsSetupParams nsp) {
		int nVolTotal = nsp.getNumVols() + nsp.getNSkipTmpts();
		int[] paddedSessLabels = new int[nVolTotal]; 
	    int skipIdx = 0;
	    int labIdx = 0;
		for (int i = 0; i < nVolTotal; ++i) {
	    	if (nsp.getSkipTmpts()[skipIdx] == i) {
	    		paddedSessLabels[i] = -1;
	    		++skipIdx;
	    	}
	    	else {
	    		paddedSessLabels[i] = nsp.getSessLabels()[labIdx];
	    		++labIdx;
	    	}	    	
	    }
		return paddedSessLabels;
	}

	/** Returns double 1D array containing AND mask created from input mask files, i.e.,
	 * intersection of all non-zero voxels in input mask files. 
	 * <p> Mask files must be 3D volumes in Nifti or ANALYZE format. All mask volumes must have 
	 * the same number of voxels. 
	 * <p> Mask voxels are conventionally set to either 1 or 0, but there is no requirement that 
	 * non-zero voxels in input files be set to 1, and there is no guarantee that non-zero voxels 
	 * in returned AND mask will be set to 1 unless input files follow this convention.
	 * 
	 * @param maskFileNames array of file names for masks
	 * @return AND mask created from input maskfiles, as double array 
	 * 			(length number of voxels in a single mask volume)
	 */
	public static double[] getANDMask(String[] maskFilenames) throws IOException {
		int numMaskFiles = maskFilenames.length;
		double[][] andMask = null;
		int maskSize = 0;
		Vector<String> uniqMaskFiles = new Vector<String>();
		for (int i = 0; i < numMaskFiles; ++i) {
			if (debug) {
				System.out.println("Reading maskfile " + maskFilenames[i] + "...");
			}
			if (maskFilenames[i] != null && !uniqMaskFiles.contains(maskFilenames[i])) {
				uniqMaskFiles.add(maskFilenames[i]);
			}
		}
		for (String mf : uniqMaskFiles) {
			try {
				int[] skipTmpts = null;
				double[][] currMaskData = readNiftiData(mf, skipTmpts, 0);
				if (currMaskData == null || currMaskData.length != 1) {
					throw new IllegalArgumentException("Input masks must be 3D");
				}
				if (andMask == null) {
					andMask = currMaskData;
					maskSize = currMaskData[0].length;				
				}
				else {
					if (currMaskData[0].length != maskSize) {
						throw new IllegalArgumentException("All input masks must have the same dimensions");					
					}
				
					// zero out voxels in andMask that are zero in current mask
					for (int j = 0; j < maskSize; ++j) {
						if (currMaskData[0][j] == 0) {
							andMask[0][j] = 0;
						}
					}
				}
			}
			catch (IOException e) {
				throw new IOException(e.getMessage());
			}
			catch (NullPointerException npe) {
				throw new IOException("Unable to create AND Mask");
			}
		}
		if (andMask == null) {
			throw new IOException("Unable to create AND Mask");
		}
		
		// TODO: set all non-zero voxels in AND mask to 1. 
		
		return andMask[0];
	}
	
	
	/** Returns double 2D array containing AND mask created from input maskfiles
	 * (Maskfiles must be Nifti or ANALYZE format)
	 *
	 * @param maskFileNames
	 * @return AND mask created from input maskfiles, as byte array 
	 * 			(length number of voxels in input mask files)
	 */
	private static byte[] getANDMaskBytes(String[] maskFilenames) 
		throws IOException {
		if (debug) {
			System.out.println("Getting byte mask...");
		}
		int numMaskFiles = maskFilenames.length;
		byte[][] andMask = null;
		int maskSize = 0;
		Vector<String> uniqMaskFiles = new Vector<String>();
		for (int i = 0; i < numMaskFiles; ++i) {
			if (debug) {
				System.out.println("Reading maskfile " + maskFilenames[i] + "...");
			}
			if (!uniqMaskFiles.contains(maskFilenames[i])) {
				uniqMaskFiles.add(maskFilenames[i]);
			}
		}
		for (String mf : uniqMaskFiles) {
			try {
				byte[][] currMaskData = readNiftiDataBytes(mf, null, 0);
				if (currMaskData == null || currMaskData.length != 1) {
					throw new IllegalArgumentException("Input masks must be 3D");
				}
				if (andMask == null) {
					andMask = currMaskData;
					maskSize = currMaskData[0].length;				
				}
				else {
					if (currMaskData[0].length != maskSize) {
						throw new IllegalArgumentException("All input masks must have the same dimensions");					
					}
				
					// zero out voxels in andMask that are zero in current mask
					for (int j = 0; j < maskSize; ++j) {
						if (currMaskData[0][j] == 0) {
							andMask[0][j] = 0;
						}
					}
				}
			}
			catch (IOException e) {
				throw new IOException(e.getMessage());
			}
			catch (NullPointerException npe) {
				throw new IOException("Unable to create AND Mask");
			}
		}
		if (andMask == null) {
			throw new IOException("Unable to create AND Mask");
		}
		// strip extraneous 2nd dim
		return andMask[0];
		
	}
	
	
//	For testing (args[1] == Nifti filename)
	
// NOTE: Careful!  x and y order is arranged so .img/hdr files look right when read in by 
// idl read_analyze; when saving as .nii, x and y dims were switched (so had [1][y][x] instead)
	private static void main2(String[] args) {
		Nifti1Dataset newNiftiDS;
		npairs.shared.matlib.ColtMatrix sampMat;
		npairs.shared.matlib.Matrix sampMat2;
		if (args[0].equals("mask")) {
			newNiftiDS = new Nifti1Dataset();
			
			newNiftiDS.setHeaderFilename(args[1] + ".nii");
		
			newNiftiDS.setDataFilename(args[1]);

			newNiftiDS.setDatatype((short)64);
			newNiftiDS.setDims((short)3, (short)4, (short)3, (short)1,
					(short)0, (short)0, (short)0, (short)0);
			sampMat = new npairs.shared.matlib.ColtMatrix(3,4);
			sampMat.setIdentity();
			System.out.println("Mask matrix: ");
			sampMat.print();
			sampMat2 = new npairs.shared.matlib.ColtMatrix(3,4);
			sampMat2.setIdentity();
			sampMat2.mult(2.0);
			System.out.println("Mask2 matrix: ");
			sampMat2.print();
		}

		else {
			newNiftiDS = new Nifti1Dataset();
			// set header filename first and include .nii ext to set dataset file type to .nii
			// (note: also need to write header before writing data when saving .nii file)
			// (note 2: doesn't matter whether .nii ext is included when setting data filename)
			newNiftiDS.setHeaderFilename(args[0] + ".nii");
			newNiftiDS.setDataFilename(args[0]);

			newNiftiDS.setDatatype((short)64);
			newNiftiDS.setDims((short)3, (short)4, (short)3, (short)2, 
					(short)0, (short)0, (short)0, (short)0);

			sampMat = new npairs.shared.matlib.ColtMatrix(3, 4);
			sampMat.setRandom();
			System.out.println("First random 2D colt matrix to be replicated in Nifti dataset: ");
			sampMat.print();
			sampMat2 = sampMat.mult(2.0);
			System.out.println("Second 2D matrix to be replicated in Nifti dataset: ");
			sampMat2.print();
		}
		// sampData3D is in [Z][Y][X] format, as per Nifti1Dataset.writeVol specs.
		double[][][] sampData3D = new double[2][3][4];
		sampData3D[0] = sampMat.toArray();
		sampData3D[1] = sampMat2.toArray();
		
		try {
			// remember: save header first when writing .nii files
			System.out.println("Writing header to disk... Filename: " + newNiftiDS.getHeaderFilename());
			newNiftiDS.writeHeader();
			System.out.println("Writing data to disk... Filename: " + newNiftiDS.getDataFilename());
			newNiftiDS.writeVol(sampData3D, (short)0);
			newNiftiDS.writeVol(sampData3D, (short)1);
			newNiftiDS.writeVol(sampData3D, (short)2);

		}
		catch (IOException e) {
			e.printStackTrace();
		}	
		
		// Now try loading info:
		try {
			int[] volDims = getVolDims3D(newNiftiDS.getDataFilename());
			System.out.println("Vol dims: ");
			NpairsIO.print(volDims);
		}
		catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
	}

	
	private static void main4(String[] args) {
		String datafile = "/home/anita/plsnpairs/spreng/data/RS4/001_bld008_rs4.nii";
		Nifti1Dataset dataNifti = new Nifti1Dataset(datafile);
		try {
			dataNifti.readHeader();
			dataNifti.printHeader();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	/** Test masked volume write to disk */
	private static void main5(String[] args) {
		
//		String dataFile = "/home/anita/plsnpairs/spreng/data/RS4/001_bld008_rs4.nii";
//		Nifti1Dataset dataNifti = new Nifti1Dataset(dataFile);
//		try {
//			dataNifti.readHeader();
//			dataNifti.printHeader();
//		}
//		catch (Exception e) {
//			e.printStackTrace();
//		}
		
		String resFile = "/home/anita/plsnpairs/spreng/results/10subj_30cls_allRuns_3splits/" +
			"10subj_30cls_allRuns_3splits_020pc_NPAIRSJresult.mat";
		// read in data and mask info
		ArrayList<String> fields = new ArrayList<String>();
		fields.add("st_coords");
		fields.add("npairs_result");
		fields.add("st_voxel_size");
		fields.add("st_origin");
		fields.add("st_dims");
		
		try {
			Map<String, MLArray> mResultInfo = new NewMatFileReader(
					resFile, new MatFileFilter(fields.toArray(new String[0]))).getContent();
			int[] st_coords = ((MLDouble) mResultInfo.get("st_coords")).getIntFirstRowOfArray();
			int[] volDims4D = ((MLDouble) mResultInfo.get("st_dims")).getIntFirstRowOfArray();
			int[] volDims = {volDims4D[0], volDims4D[1], volDims4D[3]};
			double[] voxSize = ((MLDouble) mResultInfo.get("st_voxel_size")).getFirstRowOfArray();
			voxSize[0] = -voxSize[0];
			int[] origin = ((MLDouble) mResultInfo.get("st_origin")).getIntFirstRowOfArray();
			
		    double[][] zs_eigim = ((MLDouble) ((MLStructure) mResultInfo.get("npairs_result")).
		    	getField("zscored_brainlv_avg")).getArray();
		    zs_eigim = MLFuncs.transpose(zs_eigim);
		    double[] zs_vol1D = zs_eigim[0];
		    System.out.println("Size st_coords: " + st_coords.length);
		    System.out.println("Size zs_vol1D: " + zs_vol1D.length);
		    
		    // save to disk
		    String saveName = "/home/anita/Desktop/10subj_30cls_allRuns_3splits_020pc.zs-eigim_voxqoffset.nii";
		    System.out.print("Saving " + saveName + "...");
		    writeVol(zs_vol1D, st_coords, true, volDims, voxSize, origin, saveName);
		    System.out.println("[DONE]");
		    
		 
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	

	 
}
