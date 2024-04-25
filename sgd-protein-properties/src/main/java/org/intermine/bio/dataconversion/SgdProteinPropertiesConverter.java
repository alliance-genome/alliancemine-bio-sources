package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2011 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.Reader;
import java.util.Iterator;

import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.regex.Pattern;
import java.util.regex.Matcher;


/**
 * 
 * @author
 */
public class SgdProteinPropertiesConverter extends BioFileConverter
{
	//
	protected static final Logger LOG = Logger.getLogger(SgdProteinPropertiesConverter.class);
	private static final String DATASET_TITLE = "Protein Aggregation data";
	private static final String DATA_SOURCE_NAME = "SGD";
	private final Map<String, Item> proteinIdMap = new HashMap<String, Item>();
	private final Map<String, Item> pubmedIdMap = new HashMap<String, Item>();

	/**
	 * Constructor
	 * @param writer the ItemWriter used to handle the resultant items
	 * @param model the Model
	 */
	public SgdProteinPropertiesConverter(ItemWriter writer, Model model) {
		super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
	}

	/**
	 * 
	 *
	 * {@inheritDoc}
	 */
	public void process(Reader reader) throws Exception {
		processPropertiesDataFile(reader); 
		storeProteins();
	}


	/**
	 * 
	 * @param reader
	 * @throws Exception
	 * @throws ObjectStoreException
	 */
	private void processPropertiesDataFile(Reader preader) throws Exception, ObjectStoreException {

		/*
		 * Sample line 
		 * gene [tab] 
		 * ORFname [tab] 
         * Property [tab] 
         * Reversible? [tab] 
         * Condition [tab] 
         * location
         * Reference
		 */   	 
		System.out.println("Processing Protien Properties Data file....");    

		Iterator<?> tsvIter;
		try {
			tsvIter = FormattedTextParser.parseTabDelimitedReader(preader);
		} catch (Exception e) {
			throw new BuildException("cannot parse file: " + getCurrentFile(), e);
		}

		while (tsvIter.hasNext()) {

			String[] line = (String[]) tsvIter.next();

			if (line.length < 7) {
				LOG.error("Couldn't process line. Expected 7 cols, but was " + line.length);
				continue;
			}

			String protein =  line[1].trim();     
			String property = line[2].trim();
			String reversible =  line[3].trim(); 
			String condition = line[4].trim();
			String location = line[5].trim();
			String pmid = line[6].trim();

			newProduct(protein, property, reversible, condition, location, pmid);

		}

		preader.close();

	}

	/**
	 * 
	 * @param proteinId
	 * @param modSite
	 * @param modType
	 * @param source
	 * @param pmid
	 * @throws ObjectStoreException
	 * @throws Exception
	 */
	private void newProduct(String proteinId, String property, String reversible, String condition, String location, String pmid)
					throws ObjectStoreException, Exception {		

		Item protein = getProteinItem(proteinId);		
		Item pmods = getProteinProperty(property, reversible, condition, location, pmid);
		protein.addToCollection("proteinProperties", pmods.getIdentifier());

	}
	/**
	 * 
	 * @param geneId
	 * @return
	 * @throws ObjectStoreException
	 */

	private Item getProteinItem(String proteinId) throws ObjectStoreException{

		Item protein = proteinIdMap.get(proteinId);

		if (protein == null) {      	
			protein = createItem("Protein");
			proteinIdMap.put(proteinId, protein);
			protein.setAttribute("secondaryIdentifier", proteinId);                  		
		}

		return protein;

	}
	/**
	 * 
	 * @param modSite
	 * @param modType
	 * @param source
	 * @param pmid
	 * @return
	 */
	private Item getProteinProperty(String property, String reversible, String condition, String location,
			String pmid) throws ObjectStoreException {

		Item item = createItem("ProteinProperty");

		if(StringUtils.isNotEmpty(property)){  item.setAttribute("propertyName", property);}
		if(StringUtils.isNotEmpty(reversible)){  item.setAttribute("isReversible", reversible);}
		if(StringUtils.isNotEmpty(condition)){  item.setAttribute("condition", condition);}
		if(StringUtils.isNotEmpty(location)){ item.setAttribute("location", location);}
	
		item.setAttribute("source", "SGD");
		
		Item publication = pubmedIdMap.get(pmid);

		if (publication == null) {
			
			publication = createItem("Publication");			
			publication.setAttribute("pubMedId", pmid);			 
			pubmedIdMap.put(pmid, publication);
			item.setReference("publication", publication);  
			try {
				store(publication);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}				
		}else{
			item.setReference("publication", publication);   
		}
 
	

		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}	

		return item;
	}
	/**
	 * 
	 * @throws Exception
	 */
	private void storeProteins() throws Exception{
		for (Item protein : proteinIdMap.values()) {
			try {
				store(protein);
			} catch (ObjectStoreException e) {
				throw new Exception(e);
			}
		}

	}


}
