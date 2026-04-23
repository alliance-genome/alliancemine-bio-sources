package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2026 AllianceMine
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
 * 
 * @author
 */
public class HomologGenesConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(HomologGenesConverter.class);

    private static final String DATASET_TITLE = "CGD and AspGD chromosomal_feature.tab files";
    private static final String DATA_SOURCE_NAME = "CGD and AspGD Download files";

    // Column layout of CGD/AspGD chromosomal_feature.tab
    private static final int COL_PRIMARY_ID = 0;
    private static final int COL_NAME = 1;
    private static final int COL_ALIAS = 2;
    private static final int COL_DESCRIPTION = 10;
    private static final int MIN_COLUMNS = COL_DESCRIPTION + 1;

    private Map<String, Item> genes = new HashMap<String, Item>();

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
			if (line.length < MIN_COLUMNS) {
				continue;
			}

			String primaryIdentifier = line[COL_PRIMARY_ID].trim();
			String name = line[COL_NAME].trim();
			String alias = line[COL_ALIAS].trim();
			String description = line[COL_DESCRIPTION].trim();

			LOG.debug("Processing line.." + primaryIdentifier);

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
