/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * NiftiExtractorBatch.java
 *
 * Created on 8-Nov-2010, 10:16:48 AM
 */

package pls.othertools.niftiextractor;

import java.awt.Component;
import java.awt.Container;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.filechooser.FileFilter;
import pls.shared.BfMRIResultFileFilter;
import pls.shared.NiftiFileFilter;
import pls.shared.NpairsfMRIResultFileFilter;
import pls.shared.PLSResultFileFilter;
import pls.shared.TextFileFilter;
import pls.shared.fMRIResultFileFilter;

/**
 * In order to edit this file properly you _need_ to use the swing gui builder
 * included with Netbeans. This is what I used to build this and this is why
 * you see so much ugly gui code. It is not something to be edited by hand.
 * 
 */
public class NiftiExtractorBatch extends javax.swing.JFrame {


    /** Creates new form Extractor */
    public NiftiExtractorBatch() {
        initComponents();
	}

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        dataRadioGroup = new javax.swing.ButtonGroup();
        filePathLabel = new javax.swing.JLabel();
        filepathField = new javax.swing.JTextField();
        cvRangeLabel = new javax.swing.JLabel();
        cvTextField = new javax.swing.JTextField();
        lagLabel = new javax.swing.JLabel();
        lagText = new javax.swing.JTextField();
        dataTypeLabel = new javax.swing.JLabel();
        datatype1 = new javax.swing.JRadioButton();
        datatype2 = new javax.swing.JRadioButton();
        datatype3 = new javax.swing.JRadioButton();
        pcRangeLabel = new javax.swing.JLabel();
        pcTextField = new javax.swing.JTextField();
        extractButton = new javax.swing.JButton();
        selectSetButton = new javax.swing.JButton();
        fileMenu = new javax.swing.JMenuBar();
        menuBar = new javax.swing.JMenu();
        loadMenuItem = new javax.swing.JMenuItem();
        extractMenuItem = new javax.swing.JMenuItem();
        clearMenuItem = new javax.swing.JMenuItem();
        exitMenuItem = new javax.swing.JMenuItem();
        switchViewMenu = new javax.swing.JMenu();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("NPAIRS Nifti Extractor");

        filePathLabel.setText("File path:");

        filepathField.setColumns(25);
        filepathField.setText("Enter filepath here complete with wildcard.");
        filepathField.setMaximumSize(new java.awt.Dimension(10, 27));

        cvRangeLabel.setText("CV:");

        cvTextField.setText("Enter the desired CV");

        lagLabel.setText("Lag:");

        lagText.setEnabled(false);

        dataTypeLabel.setText("Datatype:");

        dataRadioGroup.add(datatype1);
        datatype1.setSelected(true);
        datatype1.setText("Zscore");

        dataRadioGroup.add(datatype2);
        datatype2.setText("Avg Canon");

        dataRadioGroup.add(datatype3);
        datatype3.setText("Full Ref");

        pcRangeLabel.setText("PC Range:");

        pcTextField.setText("Enter the pc range");

        extractButton.setText("Extract");
        extractButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                extractButtonActionPerformed(evt);
            }
        });

        selectSetButton.setText("Select set");
        selectSetButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selectSetButtonActionPerformed(evt);
            }
        });

        menuBar.setText("File");

        loadMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_L, java.awt.event.InputEvent.CTRL_MASK));
        loadMenuItem.setText("Load");
        loadMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadMenuItemActionPerformed(evt);
            }
        });
        menuBar.add(loadMenuItem);

        extractMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_E, java.awt.event.InputEvent.CTRL_MASK));
        extractMenuItem.setText("Extract");
        extractMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                extractMenuItemActionPerformed(evt);
            }
        });
        menuBar.add(extractMenuItem);

        clearMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_C, java.awt.event.InputEvent.ALT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        clearMenuItem.setText("Clear");
        clearMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                clearMenuItemActionPerformed(evt);
            }
        });
        menuBar.add(clearMenuItem);

        exitMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_X, java.awt.event.InputEvent.ALT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        exitMenuItem.setText("Exit");
        exitMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitMenuItemActionPerformed(evt);
            }
        });
        menuBar.add(exitMenuItem);

        fileMenu.add(menuBar);

        switchViewMenu.setText("PLS extractor");
        switchViewMenu.addMenuListener(new javax.swing.event.MenuListener() {
            public void menuCanceled(javax.swing.event.MenuEvent evt) {
            }
            public void menuDeselected(javax.swing.event.MenuEvent evt) {
            }
            public void menuSelected(javax.swing.event.MenuEvent evt) {
                switchViewMenuMenuSelected(evt);
            }
        });
        fileMenu.add(switchViewMenu);

        setJMenuBar(fileMenu);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(cvRangeLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(filePathLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(pcRangeLabel, javax.swing.GroupLayout.Alignment.LEADING))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(lagLabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(dataTypeLabel, javax.swing.GroupLayout.Alignment.LEADING)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(lagText, javax.swing.GroupLayout.PREFERRED_SIZE, 332, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(datatype1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(datatype2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(datatype3, javax.swing.GroupLayout.DEFAULT_SIZE, 158, Short.MAX_VALUE))
                    .addComponent(filepathField, javax.swing.GroupLayout.DEFAULT_SIZE, 332, Short.MAX_VALUE)
                    .addComponent(cvTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 332, Short.MAX_VALUE)
                    .addComponent(pcTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 332, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(selectSetButton, javax.swing.GroupLayout.DEFAULT_SIZE, 111, Short.MAX_VALUE)
                    .addComponent(extractButton, javax.swing.GroupLayout.DEFAULT_SIZE, 111, Short.MAX_VALUE))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(filePathLabel)
                    .addComponent(filepathField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(selectSetButton, javax.swing.GroupLayout.PREFERRED_SIZE, 18, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pcRangeLabel)
                    .addComponent(pcTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cvRangeLabel)
                    .addComponent(cvTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(dataTypeLabel)
                    .addComponent(datatype1, javax.swing.GroupLayout.PREFERRED_SIZE, 22, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(datatype2)
                    .addComponent(datatype3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lagLabel)
                    .addComponent(lagText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(extractButton, javax.swing.GroupLayout.PREFERRED_SIZE, 18, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

	private void loadMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadMenuItemActionPerformed
		if(!currentlyNpairs){
			loader.removeChoosableFileFilter(txtFilter);
			loader.setFileFilter(plsResultFilter);
		}else{
			loader.setFileFilter(txtFilter);
			loader.removeChoosableFileFilter(plsResultFilter);
		}

		if(loader.showOpenDialog(this) == loader.APPROVE_OPTION){
			try{
				
				//PLS files do not have batch text files to load.
				if(!currentlyNpairs){
					String name = loader.getSelectedFile().getAbsolutePath();
					
					filepathField.setText(name);
					currentlyNpairs = false;
					
				}
				else{
					loadNpairsBatchFile(loader.getSelectedFile());
				}
			}catch(IOException e){
				JOptionPane.showMessageDialog(this, "IO error reading file",
						"Error", JOptionPane.ERROR_MESSAGE);
			}
		}
	}//GEN-LAST:event_loadMenuItemActionPerformed

	private void clearMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clearMenuItemActionPerformed
		filepathField.setText(null);
		pcTextField.setText(null);
		cvTextField.setText(null);
		lagText.setText(null);
		dataRadioGroup.clearSelection();
		datatype1.setSelected(true);
		controller.removeResultFiles();
	}//GEN-LAST:event_clearMenuItemActionPerformed
	
	private void extractButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_extractButtonActionPerformed
		File savedir;
		boolean ok;

		try {
			if(!currentlyNpairs){
				verifyFieldsPLS();
			}else{
				verifyFields();
			}
		} catch (NiftiParseError e) {
			return;
		}

		savedir = selectSaveDir();
		if (savedir == null) {
			return; //user cancelled.
		}
		
		//user has chosen the save directory.
		try{
			//remove previously loaded data before adding new data.
			controller.removeResultFiles();

			if(!currentlyNpairs){
				List<File> singleFile = new ArrayList<File>(1);

				singleFile.add(new File(filepathField.getText()));
				ok = loadFiles(singleFile);
			}else{
				ok = loadFiles(replaceWildCard(filepathField.getText(),
					                      parsePcRange(pcTextField.getText())));
			}
			
			if(ok){ //loading files ok.
				writeImages(savedir);
			}
		}catch(NiftiParseError e){
		 assert false : "this exception should never be generated (line 324)";
		}

	}//GEN-LAST:event_extractButtonActionPerformed

	private void extractMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_extractMenuItemActionPerformed
		extractButtonActionPerformed(null);
	}//GEN-LAST:event_extractMenuItemActionPerformed

	private void exitMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exitMenuItemActionPerformed
		this.dispose();
	}//GEN-LAST:event_exitMenuItemActionPerformed

	private void selectSetButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selectSetButtonActionPerformed

		if(currentlyNpairs){
			setChooser.removeChoosableFileFilter(plsResultFilter);
			setChooser.setFileFilter(npairsResultFilter);
		}else{
			setChooser.setFileFilter(plsResultFilter);
			setChooser.removeChoosableFileFilter(npairsResultFilter);
		}

		File selected;
		String name;
		setChooser.showOpenDialog(this);
		selected = setChooser.getSelectedFile();

		if(selected == null) return; //user cancelled
		
		name = selected.getName();
		//There are no pc sets for pls files so just let the user select an
		//individual file.
		if(!currentlyNpairs){
			filepathField.setText(selected.getAbsolutePath());
			return;
		}

		if(name.matches(".*\\d\\d\\dpc.*")){
			name = name.replaceFirst("\\d\\d\\dpc", "\\$pc");
			filepathField.setText(selected.getParent() + File.separator + name);
		}else{
			JOptionPane.showMessageDialog(this,"Please select a pc set");
		}
	
	}//GEN-LAST:event_selectSetButtonActionPerformed

	/**
	 * Switch the view so that the user can do either PLS or NPAIRS extraction.
	 * @param evt
	 */
	private void switchViewMenuMenuSelected(javax.swing.event.MenuEvent evt) {//GEN-FIRST:event_switchViewMenuMenuSelected
		
		clearMenuItemActionPerformed(null);
		currentlyNpairs = !currentlyNpairs;

		if(currentlyNpairs){
			this.setTitle("NPAIRS Nifti extractor");
			switchViewMenu.setText("PLS extractor");
			cvRangeLabel.setText("CV:");


			lagText.setEnabled(false);
			pcTextField.setEnabled(true);
			
			datatype1.setText("Zscore");
			datatype1.setSelected(true);
			datatype2.setText("Avg Canon");
			datatype3.setText("Fullref");
			datatype3.setVisible(true);
			selectSetButton.setText("Select set");
			currentlyNpairs = true;
		}else{
			this.setTitle("PLS Nifti extractor");
			switchViewMenu.setText("NPAIRS extractor");
			cvRangeLabel.setText("LV:");

			lagText.setEnabled(true);
			pcTextField.setEnabled(false);

			datatype1.setText("Bootstrap");
			datatype1.setSelected(true);
			datatype2.setText("Brain LV");
			datatype3.setVisible(false);
			selectSetButton.setText("Select file");
			currentlyNpairs = false;
		}
	}//GEN-LAST:event_switchViewMenuMenuSelected


	private void showError(String error){
		JOptionPane.showMessageDialog(this,
				error, "Error",JOptionPane.ERROR_MESSAGE);
	}

	private void throwError(String error) throws NiftiParseError{
		JOptionPane.showMessageDialog(this,
				error, "Error",JOptionPane.ERROR_MESSAGE);
		throw new NiftiParseError(error);
	}

	/**
	 * Bring up a dialog so the user can select a directory to save the
	 * .nii file into.
	 * @return Return the selected directory.
	 */
	private File selectSaveDir(){
		//get save directory.
		File savedir = null;
		String result;
		
		assert saverText != null;
		result = constructFileName();
		if(result == null) return null;
		
		saverText.setText(result);
		
		while(true){

			if (saver.showOpenDialog(this) == saver.APPROVE_OPTION) {
				savedir = saver.getCurrentDirectory();
				
				if (!savedir.canWrite()) {
					showError("Cannot write to chosen directory, try again.");
				}else{
					return savedir;
				}
			}else{
				return null; //user cancelled.
			}
		}
	}
	
	/**
	 * Construct the file name to be shown in the save file dialog and
	 * also to be used as the file name of the output .nii file. 
	 * @return The filename.
	 */
	private String constructFileName(){
		String resultFile = filepathField.getText();
		String suffix;
		String newprefix;
		String lagString = null;
		String datatype = null;
		
		int i;
		int[] lags = null;
		
		List<Integer> laglist = null;
		try{
			if(!currentlyNpairs){
				laglist = parsePcRange(lagText.getText());
				int laglen = laglist.size();
				lags = new int[laglen];
				for(i = 0; i < laglen; i++){
					lags[i] = laglist.get(i);
				}
			}
		}catch(NiftiParseError e){
			showError("Bad number range");
			return null; 
		}
		
		
		i = resultFile.lastIndexOf(File.separator);
		if(i != -1){
			resultFile = resultFile.substring(i+1,resultFile.length());
		}
		if(!currentlyNpairs){
			
			if (datatype1.isSelected()){
				datatype = BOOTSTRAP;
			}else if(datatype2.isSelected()){
				datatype = BRAINLV;
			}else{
				assert false : "Code should never have reached this point";
			}
			
			i = resultFile.indexOf(fMRIResultFileFilter.EXTENSION);
			
			if (i != -1) {
				suffix = fMRIResultFileFilter.EXTENSION;
			}
			else {
				i = resultFile.indexOf(BfMRIResultFileFilter.EXTENSION);
				suffix = BfMRIResultFileFilter.EXTENSION;
			}
	
			newprefix = resultFile.substring(0, i);
			newprefix += "_" + datatype;
			
			newprefix += "_lv" + cvTextField.getText().trim();
	
			lagString = "_lags" + controller.lagString(lags);
								
			newprefix += lagString + 
					suffix.substring(0,suffix.lastIndexOf(".mat"));
			resultFile = newprefix;
		
		}else{
			if (datatype1.isSelected()){
				datatype = ZSCORE;
			}else if(datatype2.isSelected()){
				datatype = AVGCANON;
			}else if(datatype3.isSelected()){
				datatype = FULLREF;
			}else{
				assert false : "Code should never have reached this point";
			}
			
			String pcString;
			pcString = pcTextField.getText().trim();
			
			i = resultFile.indexOf(NpairsfMRIResultFileFilter.EXTENSION);
			suffix = NpairsfMRIResultFileFilter.EXTENSION;

			newprefix = resultFile.substring(0,i);
			newprefix = newprefix.replace("$",pcString.replace(",", "_")); 
			newprefix += "_" + datatype;
			newprefix += "_cv" + (cvTextField.getText().trim());
			resultFile = newprefix + 
					suffix.substring(0,suffix.lastIndexOf(".mat"));
		}
		return resultFile +=".nii";
	}
	
	/**
	 * Verify the fields in the PLS GUI.
	 * @throws NiftiParseError
	 */
	private void verifyFieldsPLS() throws NiftiParseError{
		
		String path = filepathField.getText();
		String cv;
				
		if (path == null || path.trim().equals("")) {
			throwError("Enter a file to extract");
		}

		if (!path.endsWith(fMRIResultFileFilter.EXTENSION) &&
			!path.endsWith(BfMRIResultFileFilter.EXTENSION)){
			throwError("Specifed file must be a PLS file");
		}

		//Exception is thrown if pc's parsed are bad.
		String lagstring = lagText.getText();
		if(lagstring == null || lagstring.trim().equals("")){
			throwError("Lags must be specified");
		}
		parsePcRange(lagstring);

		if(!datatype1.isSelected() && !datatype2.isSelected()){
			throwError("Select a datatype");
		}

		cv = cvTextField.getText();
		if(cv == null || cv.trim().equals("")){
			throwError("Enter an LV");
		}else if(!cv.matches("\\s*\\d+\\s*")){
			throwError("Bad LV");
		}else if(Integer.parseInt(cv) <= 0){
			throwError("Bad LV");
		}
		
		cvTextField.setText(cv.trim());
		//All fields are ok
	}
	
	/**
	 * Verify the fields in the NPAIRS GUI view.
	 * @throws NiftiParseError
	 */
	private void verifyFields() throws NiftiParseError{
		List<String> fields = new ArrayList<String>();
		String path = filepathField.getText();
		String pcs = pcTextField.getText();
		String cv = cvTextField.getText();
		
		if(path == null || pcs == null || cv == null)
			throwError("Missing required fields");

		path = path.trim();
		pcs = pcs.trim();
		cv = cv.trim();
		
		fields.add(path);

		if(datatype1.isSelected()){
			fields.add(ZSCORE);
		}else if(datatype2.isSelected()){
			fields.add(AVGCANON);
		}else if(datatype3.isSelected()){
			fields.add(FULLREF);
		}else{
			assert false : "Extracting NPAIRS but no radio button selected";
		}

		fields.add("CV"+cv);
		for (String pc : pcs.split(","))
			fields.add(pc);

		parseNpairs(fields.toArray(new String[fields.size()]));
	}
	
	//All fields have been verified and data has been loaded once this function
	//is called. All that is left is to verify that the cv/lv and lags exist
	//for the chosen datatype.
	private void writeImages(File savelocation){
		
		int cv = Integer.parseInt(cvTextField.getText());
		List<Integer> lags = null;
		String datatype = null;
		String skipped = "";
		String retmessage = null;
		boolean validCV = false;

		if(!currentlyNpairs){
			try{
				lags = parsePcRange(lagText.getText());
			}catch(NiftiParseError e){
				assert false : "Parse exception thrown." +
						"\nShould never have occurred.";
			}

			if (datatype1.isSelected()){
				datatype = BOOTSTRAP;
			}else if(datatype2.isSelected()){
				datatype = BRAINLV;
			}else{
				assert false : "Code should never have reached this point";
			}
		}else{
			if (datatype1.isSelected()){
				datatype = ZSCORE;
			}else if(datatype2.isSelected()){
				datatype = AVGCANON;
			}else if(datatype3.isSelected()){
				datatype = FULLREF;
			}else{
				assert false : "Code should never have reached this point";
			}
		}

		if(currentlyNpairs){
			for (String resultFile : controller.getResultFiles()){
				validCV = controller.validCV(resultFile, cv, datatype,false);
				if(!validCV){
					showError("Aborting extraction.\nInvalid CV for file: " + resultFile);
					return;
				}
			}
			
			String filename = constructFileName();
			if(filename == null) 
				assert false : "PC field is incorrect";
			
			//the cv must be 0 based when passed to the controller.
			retmessage = controller.writeNiftiImage(
									savelocation.getAbsolutePath() 
									+ File.separator + filename, 
									controller.getResultFiles(), 
									cv-1, datatype);
			if(retmessage != null){
				showError(retmessage);
			}else{
				JOptionPane.showMessageDialog(
					this,"Extraction complete","Success",
					JOptionPane.INFORMATION_MESSAGE);
			}
			return;
		}
		
		int[] simpleLags = null;

		simpleLags = new int[lags.size()];
		for (int i = 0; i < lags.size(); i++) {
			//lags must be 0 based when passed to the controller.
			simpleLags[i] = lags.get(i);
		}
		
		for (String resultFile : controller.getResultFiles()){
			
			validCV = controller.validCV(resultFile, cv, datatype,true);

			int windowsize = controller.getWindowSize(resultFile);
			for(int lag : lags){
				if (lag > windowsize-1 || lag < 0) {
					showError("Invalid lag (" + lag + ") " 
							+ " size for file:\n"
							+ resultFile + 
							"\nAborting remaining extractions.");
					return;
				}
			}

			if(validCV){
				//the lv must be 0 based when passed to the controller.
				retmessage = controller.writeNiftiImage(
											savelocation, resultFile,
											cv-1, simpleLags, datatype);
				if(retmessage != null){
					showError("Aborting extraction due to error. " + retmessage);
					return;
				}
				
			}else{
				skipped +=  "\nNonexistent LV for file ";
				skipped += new File(resultFile).getName() + ". skipping";
			}
		}

		if(!skipped.equals("")){
			showError(skipped);
		}else{
			JOptionPane.showMessageDialog(this,"Extraction complete","Success",
					JOptionPane.INFORMATION_MESSAGE);
		}
	}
	
	

	/**
	 * Load an npairs batch text file.
	 * @param input
	 * @throws IOException
	 */
	private void loadNpairsBatchFile(File input) throws IOException{
		FileReader fr = new FileReader(input);
		BufferedReader bf = new BufferedReader(fr);

		String tokens[] = bf.readLine().split(",");
		String path;
		
		if(tokens.length == 0){
			try{
				throwError("Malformed file.");
			}catch(NiftiParseError e){}
			bf.close();
			return;
		}

		path = tokens[0].trim();

		if(path.endsWith(NpairsfMRIResultFileFilter.EXTENSION)){
			try{
				parseNpairs(tokens);
				currentlyNpairs = true;
			}catch(NiftiParseError e){}
		}else{
			try{
				throwError("Unknown file type specified in the filepath");
			}catch(NiftiParseError e){}
			bf.close();
			return;
		}
	}

	/**
	 * Parse a text file that contains information to fill in the GUI 
	 * fields. This function is also used when verifying the NPAIRS GUI 
	 * fields. See verifyFields() The fields are converted into a string 
	 * array, one that looks identical to the string array you would get
	 * if you loaded in a text file. The array is then parsed for errors.
	 * @param tokens Comma separated tokens.
	 * @throws NiftiParseError
	 */
	private void parseNpairs(String[] tokens) throws NiftiParseError{
		String path;
		String datatype;
		String CV;
		String pcRange = "";
		
		if(tokens.length < 4){
			throwError("Text file is missing fields");
		}

		path = tokens[0].trim();
		datatype = tokens[1].trim();
		CV = tokens[2].trim();

		for(int i=3; i < tokens.length-1; i++){
			pcRange += tokens[i] + ",";
		}pcRange += tokens[tokens.length-1];

		parsePcRange(pcRange);

		//Error check parsed fields.
		
		if(!datatype.equals(ZSCORE)
			&& !datatype.equals(AVGCANON)
			&& !datatype.equals(FULLREF)){
			throwError("Unrecognized NPAIRS datatype");
		}
		if(!CV.matches("CV\\d+")){
			throwError("Bad CV");
		}else{
			CV = CV.substring(2);
			if(Integer.parseInt(CV) <= 0) throwError("Bad CV");
		}

		int wildcard = path.indexOf("$");
		if (wildcard == -1){
			throwError("File path does not contain '$' wildcard");
		}
		if(!path.endsWith(NpairsfMRIResultFileFilter.EXTENSION)){
			throwError("Unrecognized NPAIRS file suffix.");
		}

		loadNpairsFields(path,datatype,CV,pcRange);
	}

	/**
	 * Set NPAIRS fields in the GUI.
	 * @param path Set the file name field to this string.
	 * @param datatype set the datatype radio button to this datatype.
	 * @param CV Set the CV to this CV.
	 * @param pcRange Set the pc range to this string.
	 */
	private void loadNpairsFields(String path, String datatype, String CV,
			String pcRange){

		filepathField.setText(path);
		cvTextField.setText(String.valueOf(Integer.parseInt(CV)));
		pcTextField.setText(pcRange);

		
		cvRangeLabel.setText("CV:");
		lagText.setEnabled(false);
		lagText.setText(null);
		datatype1.setText("Zscore");
		datatype2.setText("Avg Canon");
		datatype3.setText("Fullref");
		datatype3.setVisible(true);

		if(datatype.equals(ZSCORE)){
			datatype1.setSelected(true);
		}else if(datatype.equals(AVGCANON)){
			datatype2.setSelected(true);
		}else if(datatype.equals(FULLREF)){
			datatype3.setSelected(true);
		}else{
			assert false : "Unknown datatype " + datatype;
		}
	}

	/**
	 * Create a list of files by subtituting the pc number in place of the
	 * wild card.
	 * @param path The general file name containing the pc number.
	 * @param pclist The list of pc numbers.
	 * @return A list of files, each one belonging to a particular pc number.
	 */
	private List<File> replaceWildCard(String path, List<Integer> pclist){
		List<File> files = new ArrayList<File>(pclist.size());
		String pcstring;
		for(Integer pc : pclist){
			//Pad number with zeros so length = 3.
			if( (pc / 100 == 0) && (pc / 10 == 0)) pcstring ="00"+pc.toString();
			else if(pc / 100 == 0) pcstring = "0"+pc.toString();
			else pcstring = pc.toString();

			files.add(new File(path.replace("$", pcstring)));
		}
		return files;
	}

	/**
	 * Load data from a list of files and store it in memory. 
	 * @param files The files to extract data from.
	 * @return true if there was no error.
	 */
	private boolean loadFiles(List<File> files){
		String error = null;

		for(File file : files){
			if(!file.canRead()){
				showError(file + " \ndoes not exist or is not readable.");
				return false;
			}
		}
		
		for(File file : files){
			error = controller.extractData(file);

			if(error != null){
				showError("Error reading file: " + error);
				break;
			}
		}

		//Error on load. Remove all files loaded before error.
		if(error != null){
			for(File file : files){
				controller.removeResultFile(file.getAbsolutePath());
			}
			return false;
		}

		return true;
	}

	/**
	 * Parse a range of pc numbers given as a string and return the numbers
	 * in a list.
	 * @param pcRange The range given as a string. i.e "1,2-5,6,3"
	 * @return A list of integers parsed from the string sorted in 
	 * ascending order.
	 * @throws NiftiParseError
	 */
	private List<Integer> parsePcRange(String pcRange) 
			throws NiftiParseError{
		
		String[] tokens = pcRange.split(",");
		Set<Integer> pcs = new TreeSet<Integer>();

		for(String range : tokens){
			range = range.trim();
			if(range.matches("\\d+-\\d+s\\d+")){
				String [] parsed = range.split("s");
				int step = Integer.parseInt(parsed[1]);
				List<Integer> nrange = parsePcRange(parsed[0]);
				int nsize = nrange.size();
				
				for(int i = 0; i<nsize; i+=step){
					pcs.add(nrange.get(i));
				}
			}
			else if(range.matches("\\d+-\\d+")){
				String[] nums = range.split("-");
				
				if(nums.length != 2) throwError("Bad number range");
				
				Integer first = null;
				Integer second = null;

				try{
					first = Integer.parseInt(nums[0]);
					second = Integer.parseInt(nums[1]);
				}catch(NumberFormatException e){
					throwError("Bad number range");
				}
				
				if(currentlyNpairs){
					//parsing pcs, and a pc can't be 0.
					if(second <= 0 || first <= 0){
						throwError("Bad number range");
					}
				}else{
					//lags must be >=0.
					if(second < 0 || first < 0){
						throwError("Bad number range");
					}
				}
				
				if(second < first){
					throwError("Bad number range");
				}
				while(first < second+1){
					pcs.add(first++);
				}
			}else if(range.matches("\\d+")){
				try{
					Integer r = Integer.parseInt(range);
					if(currentlyNpairs){
						if(r <= 0) throwError("Bad number range");
					}else{
						if(r < 0)  throwError("Bad number range");
					}
					pcs.add(Integer.parseInt(range));
				}catch(NumberFormatException e){
					throwError("Bad number range");
				}
			}else{
				throwError("Bad number range");
			}
		}
		return new ArrayList<Integer>(pcs);
	}

	/**
	 * Return the textfield in the JFileChooser container;
	 * @param item a container (initially the JFileChooser extraction dialog).
	 */
	private static Container disableTextField(Container item){
		if(item instanceof JTextField){
			((JTextField)item).setEditable(false);
			
			return item; 
		}
		for(Component c : item.getComponents()){
			if(c instanceof Container){
				Container r = disableTextField((Container)c);
				if(r != null) return r;
			}
		}
		return null;
	}
    /**
    * @param args the command line arguments
    */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new NiftiExtractorBatch().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuItem clearMenuItem;
    private javax.swing.JLabel cvRangeLabel;
    private javax.swing.JTextField cvTextField;
    private javax.swing.ButtonGroup dataRadioGroup;
    private javax.swing.JLabel dataTypeLabel;
    private javax.swing.JRadioButton datatype1;
    private javax.swing.JRadioButton datatype2;
    private javax.swing.JRadioButton datatype3;
    private javax.swing.JMenuItem exitMenuItem;
    private javax.swing.JButton extractButton;
    private javax.swing.JMenuItem extractMenuItem;
    private javax.swing.JMenuBar fileMenu;
    private javax.swing.JLabel filePathLabel;
    private javax.swing.JTextField filepathField;
    private javax.swing.JLabel lagLabel;
    private javax.swing.JTextField lagText;
    private javax.swing.JMenuItem loadMenuItem;
    private javax.swing.JMenu menuBar;
    private javax.swing.JLabel pcRangeLabel;
    private javax.swing.JTextField pcTextField;
    private javax.swing.JButton selectSetButton;
    private javax.swing.JMenu switchViewMenu;
    // End of variables declaration//GEN-END:variables

	final static String AVGCANON = "ac";
	final static String ZSCORE = "zs";
	final static String FULLREF = "fr";
	final static String BOOTSTRAP = "bs";
	final static String BRAINLV = "bv";

	private NiftiExtractorController controller = new NiftiExtractorController();
	private JFileChooser loader = null;
	private JFileChooser saver = null;
	private JFileChooser setChooser = null;
	private FileFilter txtFilter = new TextFileFilter();
	private FileFilter plsResultFilter = new PLSResultFileFilter();
	private FileFilter npairsResultFilter = new NpairsfMRIResultFileFilter();
	private JTextField saverText = null;
	private boolean currentlyNpairs = true;

	{
		setTitle("NPAIRS Nifti Extractor");
		loader = new JFileChooser();
		loader.setMultiSelectionEnabled(false);
		loader.removeChoosableFileFilter(loader.getAcceptAllFileFilter());

		setChooser = new JFileChooser();
		setChooser.setMultiSelectionEnabled(false);
		setChooser.removeChoosableFileFilter(
				setChooser.getAcceptAllFileFilter());

		saver = new JFileChooser();
		saver.setFileHidingEnabled(true);
		saver.setMultiSelectionEnabled(false);
		/*saver.removeChoosableFileFilter(saver.getAcceptAllFileFilter());
		saver.addChoosableFileFilter(new FileFilter(){

			@Override
			public boolean accept(File arg0) {
				if(arg0.isDirectory()) return true;
				return false;
			}

			@Override
			public String getDescription() {
				return "Select a directory";
			}
			
		});*/
		saver.setDialogTitle("Select a directory");
		saverText = (JTextField)disableTextField(saver);
		
	}
	
	private class NiftiParseError extends Exception {

		String error = null;

		private NiftiParseError(){}
		public NiftiParseError(String err) {
			error = err;
		}

		@Override
		public String getMessage() {
			return error;
		}

	}
}