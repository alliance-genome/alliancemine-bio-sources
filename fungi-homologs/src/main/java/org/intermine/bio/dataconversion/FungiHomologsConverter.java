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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class FungiHomologsConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "Orthologs from FungiDB using OrthoMCL";
    private static final String DATA_SOURCE_NAME = "FungiDB";
    private static HashMap taxonIds = new HashMap();
    private Map<String, String> genes = new HashMap<String, String>();
    
    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public FungiHomologsConverter(ItemWriter writer, Model model) throws Exception {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        loadTaxonIds();
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

    	Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
    	while (lineIter.hasNext()) {

    		String[] line = lineIter.next();
    		if (line.length < 3 && StringUtils.isNotEmpty(line.toString())) {
    			throw new RuntimeException("Invalid line, should be at least 3 columns but is '"
    					+ line.length + "' instead");
    		}

    		String gene2 = line[0];       //homolog ID
    		String gene1 = line[1];       //yeast gene - systematic name
    		String organism = line [2];   //homolog Organism
    		
    		String taxonName = (String) taxonIds.get(organism);
    		
    		if(taxonName != null) {
    		
    			String refId1 = null;
    			String refId2 = refId2 = parseGene(gene2, organism, "primaryIdentifier");

    			//System.out.println("gene1...: "+ gene1);
    			
    			String geneIds[] = gene1.split(",");
    			if(geneIds.length > 0){

    				for(int i=0; i < geneIds.length; i++) {

    					String gene = geneIds[i].trim();   	
    					String refId = parseGene(gene, "Saccharomyces cerevisiae", "secondaryIdentifier");
    					if (refId == null || refId2 == null) {
    						continue;
    					}
    					processHomologue(refId, refId2);
    					processHomologue(refId2, refId);

    				}

    			}else{
    				refId1 = parseGene(gene1, "Saccharomyces cerevisiae", "secondaryIdentifier");

    				if (refId1 == null || refId2 == null) {
    					continue;
    				}
    				processHomologue(refId1, refId2);
    				processHomologue(refId2, refId1);
    			}
    		}

    	}

    }
   
    // save homologue pair
    private void processHomologue(String gene1, String gene2)
        throws ObjectStoreException {
        Item homologue = createItem("Homologue");
        homologue.setReference("gene", gene1);
        homologue.setReference("homologue", gene2);
        homologue.setAttribute("type", "homologue");
        store(homologue);
    }

    private String parseGene(String identifier, String taxonName, String fieldName)
        throws ObjectStoreException {
        if (StringUtils.isBlank(identifier)) {
            return null;
        }
        String refId = genes.get(identifier);
        if (refId == null) {
            Item item = createItem("Gene");
            item.setAttribute(fieldName, identifier);
            String taxonId = (String) taxonIds.get(taxonName);
            item.setReference("organism", getOrganism(taxonId));
            store(item);
            refId = item.getIdentifier();
            genes.put(identifier, refId);
        }
        return refId;
    }
    
    
    private void loadTaxonIds() throws Exception {
    	
    	// if the input file changes, you have to make sure the same naming convention is used
    	//if there is no match here..nothing happens in the parser
    	
    	taxonIds.put("Saccharomyces cerevisiae","4932");
    	taxonIds.put("A. capsulatus G186AR","447093"); //3412 --ok
    	taxonIds.put("A. capsulatus NAm1", "339724"); //3472 -- ok
    	taxonIds.put("A. flavus NRRL3357", "332952"); //4004  -- ok
    	taxonIds.put( "A. fumigatus Af293", "330879"); //3679 -- ok 	
    	taxonIds.put("A. nidulans FGSC A4", "227321"); // 3734 	-- ok 
    	taxonIds.put( "A. niger ATCC 1015", "380704"); // 3817 -- ok 	
    	taxonIds.put("Aspergillus terreus NIH2624", "341663"); // 3722  	
    	taxonIds.put( "C. immitis H538.4", "396776"); // 3658 	--ok
    	taxonIds.put( "C. immitis RS", "246410"); // 3418 -- ok
     	taxonIds.put( "C. posadasii C735 delta SOWgp", "222929"); // 3383	-- ok  	
    	taxonIds.put( "N. fischeri NRRL 181", "331117"); //3820 -- ok   	    	
    	taxonIds.put("P. marneffei ATCC 18224", "441960");  // 3837 -- ok	
    	taxonIds.put( "Candida albicans SC5314", "237561"); //does not seem to appear in the file 	
    	taxonIds.put( "Schizosaccharomyces pombe 972h-", "284812"); //does not seem to appear in the file	
    	taxonIds.put( "M. oryzae 70-15", "242507"); //3675 -- ok 	
    	taxonIds.put( "N. crassa OR74A", "367110"); //3603  -- 3601 -- 2 less***	
    	taxonIds.put( "C. gattii R265", "294750");// 3034 -- ok
    	taxonIds.put("C. gattii WM276", "367775"); //3093 -- ok   	
    	taxonIds.put( "C. neoformans var. grubii H99", "235443"); //3137 - ok 
     	taxonIds.put( "Cryptococcus neoformans var. neoformans B-3501A", "283643"); //does not seem to appear in the file    	
     	taxonIds.put( "C. neoformans var. neoformans JEC21", "214684"); //3129 -- 3128 - 1 less***
     	taxonIds.put("U. maydis 521", "237631"); //3176 -- ok
    	
    }
}
