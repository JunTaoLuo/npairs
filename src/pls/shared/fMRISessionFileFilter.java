package pls.shared;

import java.io.File;

public class fMRISessionFileFilter extends ExtensionFileFilter {
	
	public static final String EXTENSION = "_fMRIsession.mat";
	
    public boolean accept(File f) {
        return f.isDirectory() || f.getName().endsWith(EXTENSION);
    }
    
    public String getDescription() {
        return "Event-related fMRI Session Files";
    }
    
    @SuppressWarnings("static-access")
	public String getExtension() {
    	return this.EXTENSION;
    }
}