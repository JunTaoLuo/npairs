package pls.shared;

import java.io.File;

public class BfMRISessionFileFilter extends ExtensionFileFilter {
	
	public static final String EXTENSION = "_BfMRIsession.mat";
	
    public boolean accept(File f) {
        return f.isDirectory() || f.getName().endsWith(EXTENSION);
    }
    
    public String getDescription() {
        return "Block fMRI Session Files";
    }
   

	@SuppressWarnings("static-access")
	public String getExtension() {
    	return this.EXTENSION;
    }
}