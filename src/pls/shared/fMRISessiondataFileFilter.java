package pls.shared;

import java.io.File;

public class fMRISessiondataFileFilter extends ExtensionFileFilter {
	
	public static final String EXTENSION = "_fMRIsessiondata.mat";
	
    public boolean accept(File f) {
        return f.isDirectory() || f.getName().endsWith(EXTENSION);
    }
    
    public String getDescription() {
        return "Event-related fMRI Sessiondata Files";
    }
    
    @SuppressWarnings("static-access")
	public String getExtension() {
    	return this.EXTENSION;
    }
}