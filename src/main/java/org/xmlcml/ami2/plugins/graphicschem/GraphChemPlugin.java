package org.xmlcml.ami2.plugins.graphicschem;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xmlcml.ami2.plugins.AMIPlugin;
import org.xmlcml.ami2.plugins.sequence.SequenceArgProcessor;

/** 
 * Plug-in for extracting chemistry from diagrams.
 * 
 * @author pm286
 */
public class GraphChemPlugin extends AMIPlugin {

	private static final Logger LOG = Logger.getLogger(GraphChemPlugin.class);
	
	static {
		LOG.setLevel(Level.DEBUG);
	}
	
	public GraphChemPlugin(String[] args) {
		super();
		this.argProcessor = new GraphChemArgProcessor(args);
	}
	
	public static void main(String[] args) {
		GraphChemArgProcessor argProcessor = new GraphChemArgProcessor(args);
		argProcessor.runAndOutput();
	}

	public GraphChemPlugin(String args) {
		super();
		this.argProcessor = new GraphChemArgProcessor(args);
	}

}