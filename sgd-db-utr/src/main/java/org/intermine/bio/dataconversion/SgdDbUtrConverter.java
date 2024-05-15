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
import java.io.BufferedReader;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.intermine.objectstore.ObjectStoreException;


/**
 *
 * @author
 */
public class SgdDbUtrConverter extends BioFileConverter
{

    private static final String DATASET_TITLE = "SGD UTRs from DB";
    private static final String DATA_SOURCE_NAME = "SGD UTRs from DB";
    private final Map<String, Item> genes = new HashMap<String, Item>();
    private final Map<String, Item> transcripts = new HashMap<String, Item>();
    private Map<String, String> chromosomes = new HashMap();
    private static final String TAXON_ID = "559292"; //4932
    private Item organism;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public SgdDbUtrConverter(ItemWriter writer, Model model) throws ObjectStoreException{
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        organism = createItem("Organism");
        organism.setAttribute("taxonId", TAXON_ID);
        organism.setAttribute("genus", "Saccharomyces");
        organism.setAttribute("species", "cerevisiae");
        organism.setAttribute("name", "Saccharomyces cerevisiae");
        organism.setAttribute("shortName", "S. cerevisiae");
        store(organism);
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        processFile(reader);
        storeTranscripts();
        storeGenes();
    }

    private void processFile(Reader preader) throws Exception,
            ObjectStoreException {

		/*
		Gene.secondaryIdentifier	-- key for gene
		Gene.transcripts.primaryIdentifier	-- key for transcript
		Gene.transcripts.symbol
		Gene.transcripts.chromosome.primaryIdentifier
		Gene.chromosomeLocation.start
		Gene.transcripts.chromosomeLocation.start
		Gene.chromosomeLocation.end
		Gene.transcripts.chromosomeLocation.end
		Gene.chromosomeLocation.strand
		 */
        System.out.println("Processing SGD transcript info data file exported from YeastMine....");

        BufferedReader br = new BufferedReader(preader);
        String line = null;
        String notes = "";


        while ((line = br.readLine()) != null) {

            String[] array = line.split("\t", -1); //keep trailing empty
            if (array.length < 9) {
                throw new IllegalArgumentException(
                        "Not enough elements (should be  9 not "
                                + array.length + ") in line: " + line);
            }
            String geneId = array[0].trim();
            String transcriptId = array[1].trim();
            String chromosome = array[3].trim();
            String geneStart = array[4].trim();
            String transcriptStart = array[5].trim();
            String geneEnd = array[6].trim();
            String transcriptEnd = array[7].trim();
            String strand = array[8].trim();

            if(geneStart.equalsIgnoreCase(transcriptStart)) {
                System.out.println("TS and GS are same : "+ geneId + "   "+ transcriptId);
                continue;
            }

            if(geneEnd.equalsIgnoreCase(transcriptEnd)) {
                System.out.println("TE and GE are same : "+ geneId + "   "+ transcriptId);
                continue;
            }

            //System.out.println("Processing line..." + geneId + "   "+ transcriptId);
            getUTRs(geneId, transcriptId, chromosome, geneStart, transcriptStart, geneEnd, transcriptEnd, strand);

        }

        preader.close();

    }

    private void getUTRs(String geneId, String transcriptId, String chromosome, String geneStart,
                         String transcriptStart, String geneEnd, String transcriptEnd, String strand) throws ObjectStoreException {

        Item gene = getGeneItem(geneId);
        Item transcript = getTranscriptItem(transcriptId);

        String chromosomeId = getChromosome(chromosome);

        //add a five-prime-utr
        Item fiveutr = null;
        String label = "";
        if(strand.equals("-1")){
            fiveutr = createItem("ThreePrimeUTR");
            label = transcriptId+"-3prime-utr";
        }else{
            fiveutr = createItem("FivePrimeUTR");
            label = transcriptId+"-5prime-utr";
        }
        fiveutr.setAttribute("primaryIdentifier", label);
        Integer start = Integer.valueOf(geneStart) - 1;
        String fivePrimeLocationRefId = getLocation(fiveutr, chromosomeId, transcriptStart, start.toString(), strand);
        fiveutr.setReference("chromosome", chromosomeId);
        fiveutr.setReference("chromosomeLocation", fivePrimeLocationRefId);
        store(fiveutr);

        //add a three-prime-utr
        Item threeutr = null;
        String label2 = "";
        if(strand.equals("-1")){
            threeutr = createItem("FivePrimeUTR");
            label2 = transcriptId+"-5prime-utr";

        }else{
            threeutr = createItem("ThreePrimeUTR");
            label2 = transcriptId+"-3prime-utr";
        }
        threeutr.setAttribute("primaryIdentifier", label2);
        Integer geneend = Integer.valueOf(geneEnd) + 1;
        String threePrimeLocationRefId = getLocation(threeutr, chromosomeId, geneend.toString(), transcriptEnd, strand);
        threeutr.setReference("chromosome", chromosomeId);
        threeutr.setReference("chromosomeLocation", threePrimeLocationRefId);
        store(threeutr);

        transcript.addToCollection("UTRs", fiveutr);
        transcript.addToCollection("UTRs", threeutr);

        gene.addToCollection("transcripts", transcript);

    }


    private String getChromosome(String identifier) throws ObjectStoreException {
        if (StringUtils.isEmpty(identifier)) {
            return null;
        }
        String refId = chromosomes.get(identifier);
        if (refId == null) {
            Item item = createItem("Chromosome");
            item.setAttribute("primaryIdentifier", identifier);
            item.setReference("organism", organism);
            refId = item.getIdentifier();
            chromosomes.put(identifier, refId);
            try {
                store(item);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
        return refId;
    }

    private String getLocation(Item subject, String chromosomeRefId,
                               String startCoord, String stopCoord, String strand)
            throws ObjectStoreException {

        String start = startCoord;
        String end = stopCoord;

        if(new Integer(start) > new Integer(end)){
            start = stopCoord;
            end =  startCoord;
        }

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

    private Item getGeneItem(String geneId)
            throws ObjectStoreException {
        Item gene = genes.get(geneId);
        if (gene == null) {
            gene = createItem("Gene");
            genes.put(geneId, gene);
            gene.setAttribute("secondaryIdentifier", geneId);
        }
        return gene;
    }

    private Item getTranscriptItem(String transcriptId)
            throws ObjectStoreException {
        Item transcript = transcripts.get(transcriptId);
        if (transcript == null) {
            transcript = createItem("MRNA");
            transcripts.put(transcriptId, transcript);
            transcript.setAttribute("primaryIdentifier", transcriptId);
        }
        return transcript;
    }

    private void storeGenes() throws Exception {
        for (Item gene : genes.values()) {
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new Exception(e);
            }
        }
    }

    private void storeTranscripts() throws Exception {
        for (Item transcript : transcripts.values()) {
            try {
                store(transcript);
            } catch (ObjectStoreException e) {
                throw new Exception(e);
            }
        }

    }
}

