package pls.shared;

import java.io.File;

public class NpairsBlockSessionFileFilter extends ExtensionFileFilter {
	
	public static final String EXTENSION = "_NPAIRS_BfMRIsession.mat";
	
    public boolean accept(File f) {
        return f.isDirectory() || f.getName().endsWith(EXTENSION);
    }
    
    public String getDescription() {
        return "NPAIRS Block Session Files";
    }
    
    @SuppressWarnings("static-access")
	public String getExtension() {
    	return this.EXTENSION;
    }
}