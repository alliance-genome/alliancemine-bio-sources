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
import org.intermine.xml.full.Item;

/*
 * 
 * @author
 */
public class AllianceGenesConverter extends BioFileConverter {

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

        //Id      Name    Description     Species   Chromosome      Start   End     Strand  SoTerm
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int count = 0;
        System.out.println("Processing Genes...");
        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            if(count == 0) { count++; continue;}
            String primaryIdentifier = line[0].trim();
            System.out.println("primaryId is ..." + primaryIdentifier);
            String secondaryIdentifier = line[1].trim();
            String name = line[2].trim();
            String description = line[3].trim();
            String autoDescription = line[4].trim();
            String origspecies = line[5].trim();
            String species = origspecies.replace("NCBITaxon:","");
            String chromosome = line[6].trim();
            String start = line[7].trim();
            String end = line[8].trim();
            String strand = line[9].trim();
            String feature_type = line[10].trim();

            String chr = "";
            if(species.equals("559292")){
                if(chromosome.equals("Mito")){
                    chr = "chrmt";
                }else {
                    chr = "chr" + chromosome;
                }
            }else{
                chr = chromosome;
            }
           Item g  = genes.get(primaryIdentifier);
           if (g != null){
               System.out.println("Is a duplicate line.." + primaryIdentifier);
               continue;
           }
            // ~~~ MOD and Chromosome ~~~
            String organism = getOrganism(species);
            String chrId = getChromosome(chr, organism);

            Item item = null;
            if (feature_type.equalsIgnoreCase("RNase_MRP_RNA_gene")) {
                item = createItem("RNaseMRPRNAGene");
            } else if (feature_type.equalsIgnoreCase("RNase_P_RNA_gene")) {
                item = createItem("RNasePRNAGene");
            } else if (feature_type.equalsIgnoreCase("SRP_RNA_gene")) {
                item = createItem("SRPRNAGene");
            } else if (feature_type.equalsIgnoreCase("antisense_lncRNA_gene")) {
                item = createItem("AntisenseLncRNAGene");
            } else if (feature_type.equalsIgnoreCase("bidirectional_promoter_lncRNA")) {
                item = createItem("BidirectionalPromoterLncRNA");
            } else if (feature_type.equalsIgnoreCase("biological_region")) {
                item = createItem("BiologicalRegion");
            } else if (feature_type.equalsIgnoreCase("blocked_reading_frame")) {
                item = createItem("BlockedReadingFrame");
            } else if (feature_type.equalsIgnoreCase("gene")) {
                item = createItem("Gene");
            } else if (feature_type.equalsIgnoreCase("gene_segment")) {
                item = createItem("GeneSegment");
            } else if (feature_type.equalsIgnoreCase("heritable_phenotypic_marker")) {
                item = createItem("HeritablePhenotypicMarker");
            } else if (feature_type.equalsIgnoreCase("lincRNA_gene")) {
                item = createItem("LincRNAGene");
            } else if (feature_type.equalsIgnoreCase("lncRNA_gene")) {
                item = createItem("LncRNAGene");
            } else if (feature_type.equalsIgnoreCase("miRNA_gene")) {
                item = createItem("MiRNAGene");
            } else if (feature_type.equalsIgnoreCase("ncRNA_gene")) {
                item = createItem("NcRNAGene");
            } else if (feature_type.equalsIgnoreCase("piRNA_gene")) {
                item = createItem("PiRNAGene");
            } else if (feature_type.equalsIgnoreCase("polymorphic_pseudogene")) {
                item = createItem("PolymorphicPseudogene");
            } else if (feature_type.equalsIgnoreCase("polypeptide")) {
                item = createItem("Polypeptide");
            } else if (feature_type.equalsIgnoreCase("protein_coding_gene")) {
                item = createItem("ProteinCodingGene");
            } else if (feature_type.equalsIgnoreCase("pseudogene")) {
                item = createItem("Pseudogene");
            } else if (feature_type.equalsIgnoreCase("pseudogenic_gene_segment")) {
                item = createItem("PseudogenicGeneSegment");
            } else if (feature_type.equalsIgnoreCase("rRNA_gene")) {
                item = createItem("RRNAGene");
            } else if (feature_type.equalsIgnoreCase("ribozyme_gene")) {
                item = createItem("RibozymeGene");
            } else if (feature_type.equalsIgnoreCase("scRNA_gene")) {
                item = createItem("ScRNAGene");
            } else if (feature_type.equalsIgnoreCase("sense_intronic_ncRNA_gene")) {
                item = createItem("SenseIntronicNcRNAGene");
            } else if (feature_type.equalsIgnoreCase("sense_overlap_ncRNA_gene")) {
                item = createItem("SenseOverlapNcRNAGene");
            } else if (feature_type.equalsIgnoreCase("snRNA_gene")) {
                item = createItem("SnRNAGene");
            } else if (feature_type.equalsIgnoreCase("snoRNA_gene")) {
                item = createItem("SnoRNAGene");
            } else if (feature_type.equalsIgnoreCase("tRNA_gene")) {
                item = createItem("TRNAGene");
            } else if (feature_type.equalsIgnoreCase("telomerase_RNA_gene")) {
                item = createItem("TelomeraseRNAGene");
            } else if (feature_type.equalsIgnoreCase("transposable_element_gene")) {
                item = createItem("TransposableElementGene");
            }

            item.setAttribute("primaryIdentifier", primaryIdentifier);
            item.setAttribute("secondaryIdentifier", secondaryIdentifier);
            if(StringUtils.isNotEmpty(name)) { item.setAttribute("symbol", name); }
            if(StringUtils.isNotEmpty(feature_type)) { item.setAttribute("featureType", feature_type);}
            if(StringUtils.isNotEmpty(description)) { item.setAttribute("briefDescription", description);}
            if(StringUtils.isNotEmpty(autoDescription)) { item.setAttribute("description", autoDescription);}
            item.setReference("organism", organism);
            item.setReference("chromosome", chrId);
            // ~~~ location ~~~
            if(!start.equals("null") || !end.equals("null")) {
                String locationRefId = getLocation(item, chrId, start, end, strand);
                item.setReference("chromosomeLocation", locationRefId);
            }
            genes.put(primaryIdentifier, item);
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

        if (StringUtils.isEmpty(chr)) {
            return null;
        }

        String unq = chr+":"+org;
        String chrId  = chromosomes.get(unq);

        if(chrId == null) {
            Item chromosome = createItem("Chromosome");
            chromosome.setAttribute("primaryIdentifier", chr);
            chromosome.setReference("organism", org);
            try {
                store(chromosome);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
            chrId = chromosome.getIdentifier();
            chromosomes.put(unq, chrId);
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
