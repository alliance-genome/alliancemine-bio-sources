package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2016 FlyMine
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

/*
 *
 * @author
 */
public class AllianceDiseaseConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "Alliance Gene data set";
    private static final String DATA_SOURCE_NAME = "AGR";
    private String licence;
    private Map<String, String> chromosomes = new HashMap();
    private Map<String, String> plasmids = new HashMap();
    private Map<String, String> sequences = new HashMap();
    private Map<String, Item> genes = new HashMap();
    private Map<String, Item> locations = new HashMap();
    private Map<String, Item> organisms = new HashMap();
    private Map<String, Item> proteins = new HashMap();
    private Map<String, Item> genesName = new HashMap();
    private Map<String, String> genesAliases = new HashMap();
    private Map<String, String> synonyms = new HashMap();
    private Map<String, Item> publications = new HashMap();


    /**
     * Construct a new AllianceGenesConverter.
     * @param database the database to read from
     * @param model the Model used by the object store we will write to with the ItemWriter
     * @param writer an ItemWriter used to handle Items created
     */
    public AllianceDiseaseConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }


    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception, ObjectStoreException{
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        //Id      Name    Species Chromosome      Start   End     Strand  SoTerm
        int count = 0;
        while (lineIter.hasNext()) {

            Item organism = null;
            String[] line = (String[]) lineIter.next();
            System.out.println("Line count..." + count + "   "+ line.length);

            String origspecies = line[0].trim();
            String taxonId = origspecies.replace("NCBITaxon:","");
            if (count == 0) {
                //organism = newOrganism(taxonId);
                count++;
            }
            if (count > 0 && count <= 17){
                count++;
                System.out.println("count..." + count);
                continue;
            }
            String productId = line[3].trim();
            String doId = line[6].trim();
            String relation_type = line[5].trim();
            String withText = line[8].trim();
            String strEvidence = line[11].trim();
            String pub = line[13].trim();
            String date_assigned = line[14].trim();
            String dataSourceCode = line[15].trim();
            String annotType = "manually curated";
            String qualifier = "";
            String annotationExtension = "";
            System.out.println("Processing line.." + productId);

        }


        System.out.println("size of genes:  " + genes.size());


    }


}
