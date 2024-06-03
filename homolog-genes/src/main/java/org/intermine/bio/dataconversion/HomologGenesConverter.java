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

import org.intermine.metadata.Model;
import org.apache.commons.lang.StringUtils;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class HomologGenesConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "CGD and AspGD chromosomal_feature.tab files";
    private static final String DATA_SOURCE_NAME = "CGD and AspGD Download files";
    private Map<String, Item> genes = new HashMap();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public HomologGenesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     * Columns within the file:

	1.  Feature name; this is the primary name
	2.  Gene name, if available
	3.  Aliases (multiples separated by |)
	4.  Feature type
	5.  Chromosome or Contig name
	6.  Start Coordinate
	7.  Stop Coordinate
	8.  Strand
	9.  Primary AspGDID
	10. Secondary AspGDID (if any)
	11. Description
	12. Date Created
	13. Sequence Coordinate Version Date (if any)
	14. Blank
	15. Blank
	16. Date of gene name reservation (if any).
	17. Has the reserved gene name become the standard name? (Y/N)
	18. Blank

     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
    	
   
		Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);

		while (lineIter.hasNext()) {

			String[] line = (String[]) lineIter.next();

			String primaryIdentifier = line[0].trim();
			String name = line[1].trim();
			String alias = line[2].trim();
			String description = line[10].trim();
			
			System.out.println("Processing line.." + primaryIdentifier);
			
			Item gene = createItem("Gene");
			gene.setAttribute("primaryIdentifier", primaryIdentifier);
			if(StringUtils.isNotEmpty(name)) { gene.setAttribute("symbol", name); }
			if(StringUtils.isNotEmpty(description)) { gene.setAttribute("briefDescription", description);}
			if(StringUtils.isNotEmpty(alias)) { 
				String newalias = alias.replaceAll("\\|", " ");
				gene.setAttribute("sgdAlias", newalias);
				}
			
			
			store(gene);
			
			         
		}

    }
    
    
}
