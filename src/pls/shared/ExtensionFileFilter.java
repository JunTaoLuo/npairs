package pls.shared;

import java.io.File;

import javax.swing.filechooser.FileFilter;

abstract public class ExtensionFileFilter extends FileFilter {
	
//	public static final String EXTENSION = null;
	
    abstract public boolean accept(File f);
    
    abstract public String getDescription();
    
    abstract public String getExtension();
       
   
 
}