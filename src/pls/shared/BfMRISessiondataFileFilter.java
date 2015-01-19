package pls.shared;

import java.io.File;

public class BfMRISessiondataFileFilter extends ExtensionFileFilter {
	
	public static final String EXTENSION = "_BfMRIsessiondata.mat";
	
    public boolean accept(File f) {
        return f.isDirectory() || f.getName().endsWith(EXTENSION);
    }
    
    public String getDescription() {
        return "Block fMRI Sessiondata Files";
    }
    
    @SuppressWarnings("static-access")
	public String getExtension() {
    	return this.EXTENSION;
    }
}