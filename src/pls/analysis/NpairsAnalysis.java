package pls.analysis;

import pls.chrome.shared.ProgressDialogWatcher;
import pls.shared.GlobalVariablesFunctions;
import pls.shared.NpairsAnalysisSetupERFileFilter;
import pls.shared.NpairsAnalysisSetupFileFilter;
import pls.shared.NpairsfMRIResultFileFilter;
import npairs.Npairs;
import npairs.io.NpairsDataLoader;
import npairs.io.NpairsIO;
import npairs.NpairsSetupParams;
import npairs.NpairsException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

import java.util.*;

/** This class is the executive manager of an entire NPAIRS analysis thread. It manages 
 * loading data, running the analysis and saving the results. 
 * 
 * @author Anita Oder
 *
 */
public class NpairsAnalysis extends ProgressDialogWatcher {
	
	String npairsSetupParamsMatFileName;
	String matlibType;
	String matlibTypeForInitFeatSel;
	
	int numNPAIRS = 1;
	
	boolean isBlocked;     // set to false for event-related fMRI
	boolean evdOnly;       // if true then just run initial EVD, save it and exit. Don't do analysis.
	int analysisNum = -1;  // if set to positive integer then run only that analysis number (given
	                       // multiple analyses that can be run from given setup file, e.g. when
	                       // pc range is set in setup file)
	boolean loadEVD;       // Overrides nsp.loadEVD.  Set to true if analysisNum is set ( > 0) 
	                       // because then we assume that EVD has already been run and saved,
	                       // even if nsp.loadEVD = false.
	
	public NpairsAnalysis(String npairsSetupParamsMatFileName, String matlibType, 
			String matlibTypeForInitFeatSel, boolean isBlocked, boolean evdOnly, int analysisNum) {
		this.npairsSetupParamsMatFileName = npairsSetupParamsMatFileName;
		this.matlibType = matlibType;
		this.matlibTypeForInitFeatSel = matlibTypeForInitFeatSel;
		this.isBlocked = isBlocked;
		this.evdOnly = evdOnly;
		this.analysisNum = analysisNum;
	}
	
	public NpairsAnalysis(String npairsSetupParamsMatFileName, String matlibType, 
			String matlibTypeForInitFeatSel, boolean isBlocked) {
		this(npairsSetupParamsMatFileName, matlibType, matlibTypeForInitFeatSel, isBlocked, false, -1);
	}
	
	public void doTask() throws Exception {
		String dataSource = "image files";
		String analysisType = "blocked";	
		if (!isBlocked) {
			dataSource = "datamats";
			analysisType = "event-related";
		}

	    String versionMessage = "Using plsnpairs version # " + 
			GlobalVariablesFunctions.getVersion() + ".\n";
		String libTypeMessage = "Using matrix library type: " + matlibType + ".\n";
		String analysisTypeMessage = "Doing " + analysisType + " fMRI analysis.\n";
		String evdOnlyMessage = "Doing initial EVD only. \nRunning and saving initial EVD.\n";
		String evdOnlyLoadEVDMessage = "Doing initial EVD only. \nEVD to be loaded from files.";
		String setupFileIDMessage = "Running NPAIRS setup file...";

		progress.postMessage(versionMessage);
		progress.postMessage(libTypeMessage);
	
		
		progress.startTask(setupFileIDMessage, "All tasks");
		
		NpairsSetupParams nsp = new NpairsSetupParams(npairsSetupParamsMatFileName, !isBlocked);
		
		File nspFileField = new File(nsp.getResultsFilePrefix().trim());
		File resultsDir = nspFileField.getParentFile();
		String resultFileName = nspFileField.getName();
		
		if (resultsDir == null) {
			// set results dir to current working directory if none specified
			resultsDir = new File(System.getProperty("user.dir")); 
		}
		
		// create results dir if it doesn't exist
		if (!resultsDir.exists()) {
			String createDirMess = "Creating results directory " + resultsDir.toString();
			progress.postMessage(createDirMess + "\n");
			Npairs.getOutput().println(createDirMess);
			boolean created = resultsDir.mkdir();
			if (!created) {
				String errorMessage = "Unable to create results directory " + resultsDir + ".";
				progress.postMessage(errorMessage + "\n");
				Npairs.getOutput().println(errorMessage);
				
				progress.endTask();
				progress.complete();
				return;
			}
		}
		
		if (!resultsDir.canWrite()) {
			String errorMessage = "Unable to save files in results directory " + resultsDir + ".";
			progress.postMessage(errorMessage + "\n");
			Npairs.getOutput().println(errorMessage);
		}
		
		else {

			loadEVD = nsp.loadEVD();

			if (evdOnly) {
				if (!loadEVD){
					progress.postMessage(evdOnlyMessage);
				}
				else {
					progress.postMessage(evdOnlyLoadEVDMessage);
				}
			}
			else {
				progress.postMessage(analysisTypeMessage);
			}

			nsp.initLogFile(evdOnly, analysisNum); // sets Npairsj.output to log file in results dir
			Npairs.getOutput().println(getDateTime());
			Npairs.getOutput().println(versionMessage);
			Npairs.getOutput().println(libTypeMessage);
			if (evdOnly) {
				if (!loadEVD) {
					Npairs.getOutput().println(evdOnlyMessage);
				}
				else Npairs.getOutput().println(evdOnlyLoadEVDMessage);
			}
			else {
				Npairs.getOutput().println(analysisTypeMessage);
			}
			Npairs.getOutput().println(setupFileIDMessage);

			
			try {
				// Copy the setup file into the result directory if it doesn't already exist there
				String ext;
				if (isBlocked) {
					ext = NpairsAnalysisSetupFileFilter.ext;
				} else {
					ext = NpairsAnalysisSetupERFileFilter.ext;
				}
				String setupFilename = nsp.getResultsFilePrefix() + ext;
				
				if (!new File(setupFilename).exists()) {
					String saveSetupMess = "Saving setup file " + setupFilename
							+ "... ";
					progress.postMessage(saveSetupMess + "\n");
					Npairs.getOutput().println(saveSetupMess);
					NpairsIO.copyFile(npairsSetupParamsMatFileName, 
							resultsDir.getAbsolutePath() 
							+ File.separator + resultFileName + ext);
				}
				else {
					String notSavingMsg = "Not saving setup file (already exists): " + setupFilename + "\n";
					progress.postMessage(notSavingMsg);
					Npairs.getOutput().println(notSavingMsg);
				}
			}
			catch (Exception e) {
				String errorMessage = "Unable to save copy of analysis file into results " +
				"directory.";
				progress.postMessage(errorMessage + "\n");
				Npairs.getOutput().println(errorMessage);

			}
			//TODO: warn user about permissions and let the user choose to stop or continue analysis
			//     String evdStatusMessage = "";
			NpairsDataLoader ndl = null;

			if (!(evdOnly && loadEVD)) {
				// if evdOnly and loadEVD, program will exit immediately; otherwise 
				// this part is executed
				String dataSourceMessage = "Loading data from " + dataSource + "...";
				progress.startTask(dataSourceMessage, "Data loading/EVD");
				Npairs.getOutput().println(dataSourceMessage);
				if (nsp.doInitFeatSelect()) {
					String ifsMsg = "Doing initial feature selection ";
					if (nsp.normEVD()) {
						ifsMsg = ifsMsg.concat("(normed EVD)...\n");
					}
					else if (nsp.doInitEVD()) {
						ifsMsg = ifsMsg.concat("(EVD)...\n");
					}
					progress.postMessage(ifsMsg);
					if (analysisNum > 0 || analysisNum == -1) {
						loadEVD = true; // assume EVD already run and saved if specifying which
						// analysis to run
						nsp.setLoadEVD(true); 
						nsp.setEVDFilePrefix(); // sets EVD file prefix to results file prefix (if not already set)
					}
					//        		else {
					//        			loadEVD = nsp.loadEVD;
					//        		}

					if (loadEVD) {
						progress.postMessage("Loading EVD from files.\n");
					}	
					else {
						String libTypeEVDMessage = "Using matrix library type: " + matlibTypeForInitFeatSel + 
						" to do initial feature selection.\n";
						progress.postMessage(libTypeEVDMessage);
						Npairs.getOutput().println(libTypeEVDMessage);
					}
				}
				ndl = new NpairsDataLoader(nsp, matlibType, 
						matlibTypeForInitFeatSel, analysisType.equals("event-related"));
				progress.endTask();
				//   	evdStatusMessage = "Finished running and saving initial EVD.";
			}


			if (evdOnly) {
				String evdDoneMessage = "Exiting program without running analysis.\n";
				progress.postMessage(evdDoneMessage);
				Npairs.getOutput().println(evdDoneMessage);
				return;
			}

			// Run NPAIRS
			checkMultNPAIRS(nsp);
			if (numNPAIRS == 1) {
				progress.startTask("Running NPAIRS analysis\n", "Analysis");
			} 
			else if (analysisNum == -1) { // run all N analyses (where N > 1)
				progress.postMessage("Starting FUUU loop.\n");
				
				
				


				double globalBest = 0.0;
				int bestQ = 0;
				int minQ = 5;
				int maxQ = numNPAIRS-1;
				int domain = maxQ - minQ;
				int lastAnalysisIdx = numNPAIRS-1;
				int i;
				double y = 0;
				
				
				
				ArrayList<ADCPoint>pointsSearched = new ArrayList<ADCPoint>();

				int sliceWidth = domain/5;
				
				
				for (int k = 0; k < domain/sliceWidth; k++){
					i = Math.min(maxQ, Math.max(0, k*sliceWidth + minQ));

					
					
					
		            if (numNPAIRS > 1) {
						String runMsg = "Running analysis # " + (i + 1) + "... ";
						progress.startTask(runMsg, 
								"Analysis # " + (i + 1));
						Npairs.getOutput().println(runMsg);
					}
					// fill in required setup parameters for current analysis
					boolean abortAnalyses = false;
					boolean pcsOutOfRange = false;
					try {
						pcsOutOfRange = nsp.setPCs(i, progress);
						if (pcsOutOfRange && i < lastAnalysisIdx) {
							//  pc range increases in each analysis so all subsequent analyses 
							// will also be out of range
							abortAnalyses = true;
							lastAnalysisIdx = i;
						}
					}
					catch (NpairsException npe) { // no pcs in range
						String errMsg = npe.getMessage() + " Stopping analysis.";
						Npairs.getOutput().println(errMsg);
						progress.postMessage("\n" + errMsg + "\n");
						progress.endTask();
						progress.complete();
						return;
					}
					if (nsp.setPCrange()) {
						nsp.setResultsFilePref(i);
					}
					try {
						Npairs npairsj = new Npairs(ndl, nsp, matlibType, i);

						
						y = npairsj.predictability;
						
						
						String saveResultMatMessage = "Saving NPAIRS results to file " + 
						nsp.getResultsFilePrefix() + NpairsfMRIResultFileFilter.EXTENSION + "...";
						progress.postMessage(saveResultMatMessage);
						Npairs.getOutput().print(saveResultMatMessage);
						double sTime = System.currentTimeMillis();
						new ResultSaver(npairsj, npairsSetupParamsMatFileName, nsp.getResultsFilePrefix());
						double tTime = (System.currentTimeMillis() - sTime) / 1000;
						progress.postMessage("Done [" + tTime + "s]\n");
						Npairs.getOutput().println("Done [" + tTime + "s]");
					}
					catch (NpairsException npe) {
						progress.printError(npe.getMessage());
						Npairs.getOutput().println(npe.getMessage());
					}
					progress.endTask(); // running current analysis
					if (abortAnalyses) {
						String stopAnalysesMsg = "Not running remaining analyses - PCs out of range.";
						progress.postMessage("\n" + stopAnalysesMsg + "\n");
						Npairs.getOutput().println(stopAnalysesMsg);
					}

					
					
					
					
					
					pointsSearched.add(new ADCPoint(i, y));
					if (y > globalBest){
						bestQ = i;
						globalBest = y;
					}
				}
				
				
				boolean done = false;
				while(!done){
				      ArrayList<ADCPoint> topPoints = new ArrayList<ADCPoint>(pointsSearched);
				      Collections.sort(topPoints, Collections.reverseOrder());

				      for (int k = 0; k < topPoints.size()/2; k++){
				        int q = topPoints.get(k).x;
				        for (int j = 1; j < pointsSearched.size(); j++){
				          int curX = pointsSearched.get(j).x;
				          int prevX = pointsSearched.get(j-1).x;
				          if (curX == q){
				            i = (curX+prevX)/2;
				            
				            
				            
				            
				            
				            if (numNPAIRS > 1) {
								String runMsg = "Running analysis # " + (i + 1) + "... ";
								progress.startTask(runMsg, 
										"Analysis # " + (i + 1));
								Npairs.getOutput().println(runMsg);
							}
							// fill in required setup parameters for current analysis
							boolean abortAnalyses = false;
							boolean pcsOutOfRange = false;
							try {
								pcsOutOfRange = nsp.setPCs(i, progress);
								if (pcsOutOfRange && i < lastAnalysisIdx) {
									//  pc range increases in each analysis so all subsequent analyses 
									// will also be out of range
									abortAnalyses = true;
									lastAnalysisIdx = i;
								}
							}
							catch (NpairsException npe) { // no pcs in range
								String errMsg = npe.getMessage() + " Stopping analysis.";
								Npairs.getOutput().println(errMsg);
								progress.postMessage("\n" + errMsg + "\n");
								progress.endTask();
								progress.complete();
								return;
							}
							if (nsp.setPCrange()) {
								nsp.setResultsFilePref(i);
							}
							try {
								Npairs npairsj = new Npairs(ndl, nsp, matlibType, i);

								
								y = npairsj.predictability;
								
								
								String saveResultMatMessage = "Saving NPAIRS results to file " + 
								nsp.getResultsFilePrefix() + NpairsfMRIResultFileFilter.EXTENSION + "...";
								progress.postMessage(saveResultMatMessage);
								Npairs.getOutput().print(saveResultMatMessage);
								double sTime = System.currentTimeMillis();
								new ResultSaver(npairsj, npairsSetupParamsMatFileName, nsp.getResultsFilePrefix());
								double tTime = (System.currentTimeMillis() - sTime) / 1000;
								progress.postMessage("Done [" + tTime + "s]\n");
								Npairs.getOutput().println("Done [" + tTime + "s]");
							}
							catch (NpairsException npe) {
								progress.printError(npe.getMessage());
								Npairs.getOutput().println(npe.getMessage());
							}
							progress.endTask(); // running current analysis
							if (abortAnalyses) {
								String stopAnalysesMsg = "Not running remaining analyses - PCs out of range.";
								progress.postMessage("\n" + stopAnalysesMsg + "\n");
								Npairs.getOutput().println(stopAnalysesMsg);
							}
				            
				            
				            
				            
				            pointsSearched.add(j, new ADCPoint(i, y));
				            
				            if (y > globalBest){
				              bestQ = i;
				              globalBest = y;
				            }
				            if (curX-prevX < 5)
				              done = true;
				            
				            break;
				          }
				        }
				      }
				    }
				
				
		        System.out.println("best q: " + bestQ);
		        System.out.println("best score: " + globalBest);
				
				
				
				
				
				
				
				
//		        SimulatedAnnealing SA = new SimulatedAnnealing(5, numNPAIRS - 1);
//		        int lastAnalysisIdx = numNPAIRS - 1;
//		        
//		        double bestScore = 0;
//		        int lastUpdate = 0;
//		        int iter = 0;
//
//		        int q = SA.getRandomNeighbour(SA.random.nextInt(SA.maxQ - SA.minQ) + SA.minQ, iter);
//		        double oldScore = 0;
//		        double newScore = oldScore;
//
//		        int bestQ = q;
//		        
//		        while (SA.temperature > SA.initialTemp/10000.0 && iter-lastUpdate < 10){
//		            SA.updateTemperature(iter);
//		            
//		            int newQ = SA.getRandomNeighbour(bestQ, iter);
//		            int i = newQ;
//		            
//		            
//		            
//		            if (numNPAIRS > 1) {
//						String runMsg = "Running analysis # " + (i + 1) + "... ";
//						progress.startTask(runMsg, 
//								"Analysis # " + (i + 1));
//						Npairs.getOutput().println(runMsg);
//					}
//
//		            
//					// fill in required setup parameters for current analysis
//					boolean abortAnalyses = false;
//					boolean pcsOutOfRange = false;
//					try {
//						pcsOutOfRange = nsp.setPCs(i, progress);
//						if (pcsOutOfRange && i < lastAnalysisIdx) {
//							//  pc range increases in each analysis so all subsequent analyses 
//							// will also be out of range
//							abortAnalyses = true;
//							lastAnalysisIdx = i;
//						}
//					}
//					catch (NpairsException npe) { // no pcs in range
//						String errMsg = npe.getMessage() + " Stopping analysis.";
//						Npairs.getOutput().println(errMsg);
//						progress.postMessage("\n" + errMsg + "\n");
//						progress.endTask();
//						progress.complete();
//						return;
//					}
//					
//					if (nsp.setPCrange()) {
//						nsp.setResultsFilePref(i);
//					}
//
//					try {
//						Npairs npairsj = new Npairs(ndl, nsp, matlibType, i);
//						newScore = npairsj.predictability;
//						
//						
//						String saveResultMatMessage = "Saving NPAIRS results to file " + 
//						nsp.getResultsFilePrefix() + NpairsfMRIResultFileFilter.EXTENSION + "...";
//						progress.postMessage(saveResultMatMessage);
//						Npairs.getOutput().print(saveResultMatMessage);
//						double sTime = System.currentTimeMillis();
//						new ResultSaver(npairsj, npairsSetupParamsMatFileName, nsp.getResultsFilePrefix());
//						double tTime = (System.currentTimeMillis() - sTime) / 1000;
//						progress.postMessage("Done [" + tTime + "s]\n");
//						Npairs.getOutput().println("Done [" + tTime + "s]");
//
//					}
//					catch (NpairsException npe) {
//						progress.printError(npe.getMessage());
//						Npairs.getOutput().println(npe.getMessage());
//					}
//					progress.endTask(); // running current analysis
//					
//					
//					if (abortAnalyses) {
//						String stopAnalysesMsg = "Not running remaining analyses - PCs out of range.";
//						progress.postMessage("\n" + stopAnalysesMsg + "\n");
//						Npairs.getOutput().println(stopAnalysesMsg);
//					}
//									
//		            
//		            
//					System.out.println("i," + iter + ",\tQ," + newQ + ",\tScore," + newScore + ";");
//					if (SA.accept(oldScore, newScore)){
//		                q = newQ;
//		                oldScore = newScore;
//		                if (newScore > bestScore){
//		                    bestScore = newScore;
//		                    bestQ = q;
//		                    lastUpdate = iter;
//		                    //if (oldScore > 0.233)   break;
//		                }
//		            }
//		            iter++;
//		        }
//
//				
//		        System.out.println(iter + " number of iterations!");
//		        System.out.println("best q: " + bestQ);
//		        System.out.println("best score: " + bestScore);
//				
			} 

			else{
				int firstAnalysisIdx = 0;
				int lastAnalysisIdx = numNPAIRS - 1;
				if (analysisNum > 0) { // just run the analysisNum'th analysis
					firstAnalysisIdx = analysisNum - 1;
					lastAnalysisIdx = analysisNum - 1;
				}
				for (int i = firstAnalysisIdx; i <= lastAnalysisIdx; ++i) {
					if (numNPAIRS > 1) {
						String runMsg = "Running analysis # " + (i + 1) + "... ";
						progress.startTask(runMsg, 
								"Analysis # " + (i + 1));
						Npairs.getOutput().println(runMsg);
					}
	
					// fill in required setup parameters for current analysis
					boolean abortAnalyses = false;
					boolean pcsOutOfRange = false;
					try {
						pcsOutOfRange = nsp.setPCs(i, progress);
						if (pcsOutOfRange && i < lastAnalysisIdx) {
							//  pc range increases in each analysis so all subsequent analyses 
							// will also be out of range
							abortAnalyses = true;
							lastAnalysisIdx = i;
						}
					}
					catch (NpairsException npe) { // no pcs in range
						String errMsg = npe.getMessage() + " Stopping analysis.";
						Npairs.getOutput().println(errMsg);
						progress.postMessage("\n" + errMsg + "\n");
						progress.endTask();
						progress.complete();
						return;
					}
					if (nsp.setPCrange()) {
						nsp.setResultsFilePref(i);
					}
	
					//				}
					//				else {
					//					if (!nsp.splitPCRangeValid()) {
					//						String splitPCRangeWarning = "Warning: PC Range for split data analyses" +
					//						"\nexceeds size of data.  Out-of-range PCs have been excluded.";
					//						progress.postMessage("\n" + splitPCRangeWarning + "\n");
					//						Npairsj.output.println(splitPCRangeWarning);
					//					}
					//					if (!nsp.fullPCRangeValid()) {
					//						String fullPCRangeWarning = "Warning: PC Range for full data analyses" +
					//						"\nexceeds size of data.  Out-of-range PCs have been excluded.";
					//						progress.postMessage("\n" + fullPCRangeWarning + "\n");
					//						Npairsj.output.println(fullPCRangeWarning);
					//					}
					//
					//				}
	
					try {
						Npairs npairsj = new Npairs(ndl, nsp, matlibType, i);
	
						String saveResultMatMessage = "Saving NPAIRS results to file " + 
						nsp.getResultsFilePrefix() + NpairsfMRIResultFileFilter.EXTENSION + "...";
						progress.postMessage(saveResultMatMessage);
						Npairs.getOutput().print(saveResultMatMessage);
						double sTime = System.currentTimeMillis();
						new ResultSaver(npairsj, npairsSetupParamsMatFileName, nsp.getResultsFilePrefix());
						double tTime = (System.currentTimeMillis() - sTime) / 1000;
						progress.postMessage("Done [" + tTime + "s]\n");
						Npairs.getOutput().println("Done [" + tTime + "s]");
	
					}
					catch (NpairsException npe) {
						progress.printError(npe.getMessage());
						Npairs.getOutput().println(npe.getMessage());
					}
					progress.endTask(); // running current analysis
					
					if (abortAnalyses) {
						String stopAnalysesMsg = "Not running remaining analyses - PCs out of range.";
						progress.postMessage("\n" + stopAnalysesMsg + "\n");
						Npairs.getOutput().println(stopAnalysesMsg);
					}				
				}				
			}
		}
		progress.endTask();
		progress.complete();
	}
	

	private void checkMultNPAIRS(NpairsSetupParams nsp) throws NpairsException, 
		IOException {
		
		numNPAIRS = nsp.getNumNPAIRS();
		if (numNPAIRS > 1) {
			progress.postMessage("Setup file contains " + numNPAIRS + " NPAIRS analyses.\n");
		}	
	}
			
    public static String getDateTime() {
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Date date = new Date();
        return dateFormat.format(date);
    }
	
//	public void finalize() throws Throwable {
//		System.out.println("NpairsjAnalysis killed.");
//		
//		super.finalize();
//	}
	
}

class SimulatedAnnealing {

    double initialTemp, temperature;
    Random random;
    int minQ, maxQ;
    
    ArrayList<Integer> Qs;

    public SimulatedAnnealing(int minQ, int maxQ){
        this.initialTemp = 10.0;
        this.temperature = this.initialTemp;
        this.random = new Random();
        
        this.minQ = minQ;
        this.maxQ = maxQ;

        this.Qs = new ArrayList<Integer>();
        for (int i = minQ; i < maxQ; i++)
            this.Qs.add(i);


    }
    

    public void updateTemperature(int iter){
        this.temperature /= 1.1;
    }


    public int getRandomNeighbour(int curQ, int iter){
        
        if (this.random.nextDouble() > 0.5){
            double range = this.maxQ - this.minQ;
            range = (range-5)*Math.exp(-iter/4.0) + 5;
            //System.out.println((int)range + ", ");

            if (this.random.nextDouble() > 0.5)
                return Math.min(Math.max(curQ + (this.random.nextInt((int)range)+1),this.minQ), this.maxQ);
            else
                return Math.min(Math.max(curQ - (this.random.nextInt((int)range)+1),this.minQ), this.maxQ);
        }
        else
            return this.random.nextInt(this.maxQ - this.minQ) + this.minQ;
    }

    public boolean accept(double oldScore, double newScore){
        if (newScore > oldScore)    return true;
        else{
            //System.out.println(newScore - oldScore);
            //System.out.println(Math.exp((newScore-oldScore)/this.temperature));
            return Math.exp((newScore-oldScore)/this.temperature) > random.nextDouble();
        }
    }


}


class ADCPoint implements Comparable{
    int x;
    double y;
    public ADCPoint(int x, double y){ this.x = x; this.y = y;}
    public int compareTo(Object o){
      ADCPoint p = (ADCPoint)o;
      if (this.y > p.y)  return 1;
      else if (this.y == p.y)  return 0;
      else  return -1;
    }
    public String toString(){
      return x + ", " + y;
    }
  }
