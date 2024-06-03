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
public class CgobHomologsConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "CGD/CGOB file download";
    private static final String DATA_SOURCE_NAME = "CGD";
    private Map<String, String> genes = new HashMap<String, String>();
    //order of columns in file to loaded should be the following
    private static final String TAXON_ID_1 = "4932"; //Saccharomyces cerevisiae S288C
    private static final String TAXON_ID_2 = "237561"; // Candida albicans SC5314
    private static final String TAXON_ID_3 = "294748"; //Candida albicans WO-1
    private static final String TAXON_ID_4 = "573826"; //Candida dubliniensis CD36 
    private static final String TAXON_ID_5 = "578454"; //Candida parapsilosis CDC317

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public CgobHomologsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
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
            if (line.length < 5 && StringUtils.isNotEmpty(line.toString())) {
                throw new RuntimeException("Invalid line, should be 5 columns but is '"
                        + line.length + "' instead");
            }

            String gene1 = line[0];
            String gene2 = line[1];
            String gene3 = line[2];
            String gene4 = line[3];
            String gene5 = line[4];

            if(gene1.equals("---") &&  gene2.equals("---") && gene3.equals("---") && gene4.equals("---") && gene5.equals("---")){
            	System.out.println("all columns empty");
            	continue;
            }
  
            String refId1 = null;
            if(!gene1.equals("---")){
            	refId1 = parseGene(gene1, TAXON_ID_1, "secondaryIdentifier"); //original yeast gene to tie up the homolog to
            }
            
            String refId2 = null;
            if(!gene2.equals("---")){
            	 refId2 = parseGene(gene2, TAXON_ID_2, "primaryIdentifier");
            }
            
            String refId3 = null;            
            if(!gene3.equals("---")){
            	 refId3 = parseGene(gene3, TAXON_ID_3, "primaryIdentifier");
            }
            
            String refId4 = null;
            if(!gene4.equals("---")){
            	 refId4 = parseGene(gene4, TAXON_ID_4, "primaryIdentifier");
            }
            
            String refId5 = null;
            if(!gene5.equals("---")){
            	 refId5 = parseGene(gene5, TAXON_ID_5, "primaryIdentifier");
            }
            
            
            if (refId1 == null || refId2 == null || refId3 == null || refId4 == null || refId5 == null) {
                continue;
            }

            processHomologue(refId1, refId2);
            processHomologue(refId2, refId1);
            
            processHomologue(refId1, refId3);
            processHomologue(refId3, refId1);
            
            processHomologue(refId1, refId4);
            processHomologue(refId4, refId1);
            
            processHomologue(refId1, refId5);
            processHomologue(refId5, refId1);

        }
    }
    	

    
    
    // save homologue pair
    private void processHomologue(String gene1, String gene2)
        throws ObjectStoreException {
        Item homologue = createItem("Homologue");
        homologue.setReference("gene", gene1);
        homologue.setReference("homologue", gene2);
        homologue.setAttribute("type", "homologue");
        homologue.setAttribute("source", "CGOB");
        store(homologue);
    }

    private String parseGene(String identifier, String taxonId, String fieldName)
        throws ObjectStoreException {
        if (StringUtils.isBlank(identifier)) {
            return null;
        }
        String refId = genes.get(identifier);
        if (refId == null) {
            Item item = createItem("Gene");
            item.setAttribute(fieldName, identifier);
            item.setReference("organism", getOrganism(taxonId));
            store(item);
            refId = item.getIdentifier();
            genes.put(identifier, refId);
        }
        return refId;
    }
    
}