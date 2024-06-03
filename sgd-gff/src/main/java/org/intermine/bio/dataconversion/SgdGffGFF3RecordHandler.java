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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.intermine.bio.dataconversion.IdResolver;
import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

/**
 * A converter/retriever for the SgdGff dataset via GFF files.
 */

public class SgdGffGFF3RecordHandler extends GFF3RecordHandler
{

    private final Map<String, Item> pubmedIdMap = new HashMap<String, Item>();
    private final Map<String, Item> geneIdMap = new HashMap<String, Item>();
    private static final String TAXON_ID = "559292"; //4932
    
    protected static final Logger LOG = Logger.getLogger(SgdGffGFF3RecordHandler.class);
    
    /**
     * Create a new SgdGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public SgdGffGFF3RecordHandler (Model model) {
        super(model);
        refsAndCollections.put("mRNA", "gene");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
    	// This method is called for every line of GFF3 file(s) being read.  Features and their
    	// locations are already created but not stored so you can make changes here.  Attributes
    	// are from the last column of the file are available in a map with the attribute name as
    	// the key.   For example:
    	//
    	//     Item feature = getFeature();
    	//     String symbol = record.getAttributes().get("symbol");
    	//     feature.setAttrinte("symbol", symbol);
    	//
    	// Any new Items created can be stored by calling addItem().  For example:
    	// 
    	//     String geneIdentifier = record.getAttributes().get("gene");
    	//     gene = createItem("Gene");
    	//     gene.setAttribute("primaryIdentifier", geneIdentifier);
    	//     addItem(gene);
    	//
    	// You should make sure that new Items you create are unique, i.e. by storing in a map by
    	// some identifier. 

    	//chrVIII	Pelechano_2013	mRNA	486202	486672	.	-	.	
    	//ID=Pelechano_chrVIII_486672_486202;Name=YHR191C.95;glu_count=0;gal_count=3;Note=Covering_one_intact_ORF;5pScore=1;3pScore=1;Parent=YHR191C
    	//ID=S000166779;Name=SC_Transcript_00009801;glu_count=0;gal_count=2;Note=Covering_one_intact_ORF;5pScore=2;3pScore=1;
    	//5pDataset=Miura_2006,Zhang_2005;3pDataset=Xu_2009_ORFs;Qualifier=Verified;Parent=YAL060

    	Item transcript = getFeature(); //class name is MRNA

    	String parentGene = record.getAttributes().get("Parent").get(0);

    	String glu_count = record.getAttributes().get("glu_count").get(0);
    	String gal_count = record.getAttributes().get("gal_count").get(0);
    	String note = record.getAttributes().get("Note").get(0);
    	String fivescore = record.getAttributes().get("5pScore").get(0);
    	String threescore = record.getAttributes().get("3pScore").get(0);
    	String fivesdataset = record.getAttributes().get("5pDataset").get(0);
    	String threedataset = record.getAttributes().get("3pDataset").get(0);
    	

    	Item gene = getGene(parentGene);
    	if (gene != null) { 	

    		if(StringUtils.isNotEmpty(glu_count)){
    			transcript.setAttribute("glucoseCount",glu_count);                    
    		}
    		if(StringUtils.isNotEmpty(gal_count)){
    			transcript.setAttribute("galactoseCount",gal_count);                    
    		}                
    		if(StringUtils.isNotEmpty(note)){
    			transcript.setAttribute("note",note);                    
    		}                
    		if(StringUtils.isNotEmpty(fivescore)){
    			transcript.setAttribute("fivePrimeScore",fivescore);                    
    		}                
    		if(StringUtils.isNotEmpty(threescore)){
    			transcript.setAttribute("threePrimeScore",threescore);                    
    		}
    		if(StringUtils.isNotEmpty(fivesdataset)){
    			transcript.setAttribute("fivePrimeDataSet",fivesdataset);                    
    		}
    		if(StringUtils.isNotEmpty(threedataset)){
    			transcript.setAttribute("threePrimeDataSet",threedataset);                    
    		}
    		
    		transcript.setReference("gene", gene.getIdentifier());

    	}

    }

    private Item getGene(String secondaryIdentifier) {
 
        Item gene = geneIdMap.get(secondaryIdentifier);
        if (gene == null) {
            gene = converter.createItem("Gene");
            geneIdMap.put(secondaryIdentifier, gene);
            gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
            gene.setReference("organism", getOrganism());
            addItem(gene);
        }
        return gene;
 
    }
}
