package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2014 YeastMine
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
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Converter that processes 2 column file Sc.(secondaryIdentifier) --> Pombe identifiers
 * YDR419W SPBC16A3.11
 * YFR027W SPBC16A3.11
 * YDR302W SPCC1450.15
 * YBR265W SPCC1450.15
 * @author
 */
public class PombeHomologsConverter extends BioFileConverter
{
    //
    private static final Logger LOG = Logger.getLogger(PombeHomologsConverter.class);

    private static final String DATASET_TITLE = "Manually Curated S.Pombe Homologs";
    private static final String DATA_SOURCE_NAME = "PomBase";
    private Map<String, String> genes = new HashMap<String, String>();
    private static final String SC_TAXON_ID = "4932";
    private static final String SP_TAXON_ID = "4896";
 

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public PombeHomologsConverter(ItemWriter writer, Model model) {
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
            if (line.length < 2 && StringUtils.isNotEmpty(line.toString())) {
                throw new RuntimeException("Invalid line, should be 2 columns but is '"
                        + line.length + "' instead");
            }

            String gene1 = line[0];
            String gene2 = line[1];

  
            String refId1 = parseGene(gene1, SC_TAXON_ID, "secondaryIdentifier");
            String refId2 = parseGene(gene2, SP_TAXON_ID, "primaryIdentifier");

            if (refId1 == null || refId2 == null) {
                continue;
            }

            processHomologue(refId1, refId2);
            processHomologue(refId2, refId1);

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
