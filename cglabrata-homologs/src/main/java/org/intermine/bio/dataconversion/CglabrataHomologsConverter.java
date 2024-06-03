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
public class CglabrataHomologsConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "C.glabrata homologs from CGOB and YGOB";
    private static final String DATA_SOURCE_NAME = "CGD";
    private Map<String, String> genes = new HashMap<String, String>();
    private static final String TAXON_ID_1 = "4932"; //Saccharomyces cerevisiae S288C
    private static final String TAXON_ID_2 = "284593"; // Candida glabrata  CBS 138


    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public CglabrataHomologsConverter(ItemWriter writer, Model model) {
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
            if (line.length < 3 && StringUtils.isNotEmpty(line.toString())) {
                throw new RuntimeException("Invalid line, should be 3 columns but is '"
                        + line.length + "' instead");
            }

            String gene1 = line[0];
            String gene2 = line[1];
            String source = line[2];


            String refId1 = null;
            if(!gene1.equals("---")){
            	refId1 = parseGene(gene1, TAXON_ID_1, "secondaryIdentifier"); //original yeast gene to tie up the homolog to
            }
            
            String refId2 = null;
            if(!gene2.equals("---")){
            	 refId2 = parseGene(gene2, TAXON_ID_2, "primaryIdentifier");
            }
 
            if (refId1 == null || refId2 == null) {
                continue;
            }

            processHomologue(refId1, refId2, source);
            processHomologue(refId2, refId1, source);

        }
    }
    
    // save homologue pair
    private void processHomologue(String gene1, String gene2, String source)
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
    