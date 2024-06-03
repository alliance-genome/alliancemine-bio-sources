package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2015 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.BufferedReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * 
 * @author
 */
public class SgdComplementationConverter extends BioFileConverter {
	//
	private static final String DATASET_TITLE = "Yeast Complementation";
	private static final String DATA_SOURCE_NAME = "SGD-BioGRID curated complementation";
	private final Map<String, Item> genes = new HashMap<String, Item>();
	private final Map<String, Item> publications = new HashMap<String, Item>();
	private final Map<String, String> homologs = new HashMap<String, String>();
	private static final String TAXON_ID = "559292"; //4932
	private static final String H_TAXON_ID = "9606";
	private Item yorganism;
	private Item horganism;

	/**
	 * Constructor
	 * 
	 * @param writer
	 *            the ItemWriter used to handle the resultant items
	 * @param model
	 *            the Model
	 */
	public SgdComplementationConverter(ItemWriter writer, Model model)
			throws ObjectStoreException {
		super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
		yorganism = createItem("Organism");
		yorganism.setAttribute("taxonId", TAXON_ID);
		store(yorganism);
		horganism = createItem("Organism");
		horganism.setAttribute("taxonId", H_TAXON_ID);
		store(horganism);
	}

	/**
	 * 
	 * 
	 * {@inheritDoc}
	 */
	public void process(Reader reader) throws Exception {
		processFile(reader);
		storeGenes();

	}

	/**
	 * 
	 * @param reader
	 * @throws Exception
	 * @throws ObjectStoreException
	 */
	private void processFile(Reader preader) throws Exception,
	ObjectStoreException {

		/*
		 * EntrezID - added by KK [0]
			Systematic name [1]
 			Gene name       
			HGNC Approved Symbol    
			HGNC ID 
			Direction of complementation  [5]  
			PMID    
			Source 
 			Note 
		 */
		System.out
		.println("Processing SGD-BioGRID complementation data file...."); 

		BufferedReader br = new BufferedReader(preader);
		String line = null;
		String notes = "";


		while ((line = br.readLine()) != null) {

			String[] array = line.split("\t", -1); //keep trailing empty
			if (array.length < 8) {
				throw new IllegalArgumentException(
						"Not enough elements (should be  8 not "
								+ array.length + ") in line: " + line);
			}
			String entrezId = array[0].trim();
			String yeastGene = array[1].trim();
			String complement = array[5].trim();
			String pmid = array[6].trim();
			String source = array[7].trim();
			if(array.length > 8) {
				 notes = array[8].trim();
			}
		    
			System.out.println("Processing line..." + yeastGene + "   "+ entrezId);

			Item ygene = getGeneItem(yeastGene, "secondaryIdentifier", yorganism);
			Item hgene = getGeneItem(entrezId, "primaryIdentifier", horganism);		
			
			getComplement(complement, notes, source, pmid, ygene, hgene);
			getComplement(complement, notes, source, pmid, hgene, ygene);

		}

		preader.close();

	}

	/**
	 * 
	 * @param c
	 * @param n
	 * @param source
	 * @param pmid
	 * @return
	 * @throws ObjectStoreException
	 */
	private void getComplement(String c, String n, String s, String pmid,
			Item yg, Item hg) throws ObjectStoreException {

		Item pub = publications.get(pmid);
		if (pub == null) {
			pub = createItem("Publication");
			pub.setAttribute("pubMedId", pmid);
			publications.put(pmid, pub);
			store(pub);
		}

		Item comp = createItem("Complement");
		
		comp.setAttribute("direction", c);
		if(StringUtils.isNotEmpty(n)) comp.setAttribute("notes", n);
		comp.setAttribute("source", s);
		comp.setReference("publication", pub);
		comp.setReference("gene", yg);
		comp.setReference("complement", hg);
		store(comp);
       
	}

    
	/**
	 * 
	 * @param geneId
	 * @return
	 * @throws ObjectStoreException
	 */
	private Item getGeneItem(String geneId, String identifier, Item org)
			throws ObjectStoreException {

		Item gene = genes.get(geneId);

		if (gene == null) {
			gene = createItem("Gene");
			genes.put(geneId, gene);
			gene.setAttribute(identifier, geneId);
			gene.setReference("organism", org);
		}

		return gene;

	}

	/**
	 * 
	 * @throws Exception
	 */
	private void storeGenes() throws Exception {
		for (Item gene : genes.values()) {
			try {
				store(gene);
			} catch (ObjectStoreException e) {
				throw new Exception(e);
			}
		}

	}

}
