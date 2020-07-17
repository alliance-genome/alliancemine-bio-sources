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
public class AllianceGenesConverter extends BioFileConverter
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
    public AllianceGenesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }


    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception, ObjectStoreException{
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        //Id      Name    Species Chromosome      Start   End     Strand  SoTerm
        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            String primaryIdentifier = line[0].trim();
            String name = line[1].trim();
            String origspecies = line[2].trim();
            String species = origspecies.replace("NCBITaxon:","");
            String chr = line[3].trim();
            String start = line[4].trim();
            String end = line[5].trim();
            String strand = line[6].trim();
            String soTerm = line[7].trim();

            System.out.println("Processing line.." + primaryIdentifier);

            // ~~~ MOD and Chromosome ~~~
            String organism = getOrganism(species);
            String chromosome = getChromosome(chr, organism);
            System.out.println("Processing Chromosome.." + chromosome + "    "+organism);

            // ~~~ gene ~~~
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
            if(StringUtils.isNotEmpty(name)) { gene.setAttribute("symbol", name); }
            if(StringUtils.isNotEmpty(soTerm)) { gene.setAttribute("featureType", soTerm);}
            gene.setReference("organism", organism);
           // gene.setReference("chromosome", chromosome);

            // ~~~ location ~~~
//            if(!start.equals("null") || !end.equals("null")) {
//                String locationRefId = getLocation(gene, chromosome, start, end, strand);
//                gene.setReference("chromosomeLocation", locationRefId);
//            }

            genes.put(primaryIdentifier, gene);

        }
        System.out.println("size of genes:  " + genes.size());

        storeGenes();

    }

    private String getLocation(Item subject, String chromosomeRefId,
                               String startCoord, String stopCoord, String strand)
            throws ObjectStoreException {

        String start = startCoord;
        String end = stopCoord;

        if (!StringUtils.isEmpty(start) && !StringUtils.isEmpty(end)) {
            subject.setAttribute("length", getLength(start, end));
        }

        Item location = createItem("Location");

        if (!StringUtils.isEmpty(start))
            location.setAttribute("start", start);
        if (!StringUtils.isEmpty(end))
            location.setAttribute("end", end);
        if (!StringUtils.isEmpty(strand))
            location.setAttribute("strand", strand);

        location.setReference("feature", subject);
        location.setReference("locatedOn", chromosomeRefId);

        try {
            store(location);
        } catch (ObjectStoreException e) {
            throw new ObjectStoreException(e);
        }
        return location.getIdentifier();
    }


    private String getChromosome(String chr, String org) throws Exception {

        String unq = chr+":"+org;
        String chrId  = chromosomes.get(unq);
        if (StringUtils.isEmpty(chr)) {
            return null;
        }
        if(chrId == null) {
            Item chromosome = createItem("Chromosome");
            chromosome.setAttribute("primaryIdentifier", chr);
            chromosome.setReference("organism", org);
            chromosomes.put(unq, chromosome.getIdentifier());
            try {
                store(chromosome);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
        return chrId;
    }


    private String getLength(String start, String end)
            throws NumberFormatException {
        Integer a = new Integer(start);
        Integer b = new Integer(end);

        // if the coordinates are on the crick strand, they need to be reversed
        // or they result in a negative number
        if (a.compareTo(b) > 0) {
            a = new Integer(end);
            b = new Integer(start);
        }

        Integer length = new Integer(b.intValue() - a.intValue());
        return length.toString();
    }


    /**
     *
     * @throws ObjectStoreException
     */

    private void storeGenes() throws ObjectStoreException {
        for (Item gene : genes.values()) {
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }



}
