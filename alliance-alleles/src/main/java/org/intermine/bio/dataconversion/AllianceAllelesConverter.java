package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2021 AllianceMine
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
import org.intermine.dataconversion.ItemWriter;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;
import org.apache.commons.lang.StringUtils;

/*
 *
 * @author
 */
public class AllianceAllelesConverter extends BioFileConverter {

    private static final String DATASET_TITLE = "Alliance Alleles data set";
    private static final String DATA_SOURCE_NAME = "DiOPT";
    private String licence;
    private Map<String, Item> genes = new HashMap();
    private Map<String, Item> alleles = new HashMap();
    private Map<String, Item> variants = new HashMap();
    private Map<String, Item> variantdetails = new HashMap();
    private Map<String, Item> alleleNames = new HashMap();

    /**
     * Construct a new AllianceGenesConverter.
     *
     * @param database the database to read from
     * @param model    the Model used by the object store we will write to with the ItemWriter
     * @param writer   an ItemWriter used to handle Items created
     */
    public AllianceAllelesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception, ObjectStoreException {

        /*Taxon	SpeciesName	AlleleId	AlleleSymbol	AlleleSynonyms
        VariantId	VariantSymbol	VariantSynonyms	VariantCrossReferences
        AlleleAssociatedGeneId	AlleleAssociatedGeneSymbol	VariantAffectedGeneId
        VariantAffectedGeneSymbol	Category	VariantsTypeId	VariantsTypeName
        VariantsHgvsNames	Assembly	Chromosome	StartPosition
        EndPosition	SequenceOfReference	SequenceOfVariant	MostSevereConsequenceName
        VariantInformationReference	HasDiseaseAnnotations	HasPhenotypeAnnotations
         */
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int count = 0;
        System.out.println("Processing Variant-Alleles...");
        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            if (count < 15) {
                count++;
                continue;
            }
            //allele information
            String origspecies1 = line[0].trim().trim();
            String species = origspecies1.replace("NCBITaxon:", "").trim();
            String org = getOrganism(species).trim();
            String alleleId = line[2].trim();
            String alleleSymbol = line[3].trim();
            String alleleSynonym = line[4].trim();
            String variantId = line[5].trim();
            String variantSymbol = line[6].trim();
            String variantSynonym = line[7].trim();
            String variantCrossRefs = line[8].trim();
            String geneId = line[9].trim();
            String alleleType = line[13].trim();
            String variantType = line[15].trim();
            String variantHgvsName = line[16].trim();
            String assembly = line[17].trim();
            String chr = line[18].trim();
            String chrStart = line[19].trim();
            String chrEnd = line[20].trim();
            String seqRef = line[21].trim();
            String seqVariant = line[22].trim();
            String mostSevere = line[23].trim();
            String variantReference = line[24].trim();
            String hasDisease = line[25].trim();
            String hasPhenotype = line[26].trim();

            Item gene = null;
            if (StringUtils.isNotEmpty(geneId)) {
                gene = getGene(geneId, org);
            }
            Item allele = processAlleles(alleleId, alleleSymbol, alleleSynonym, alleleType, variantId, variantSymbol, variantSynonym, variantCrossRefs, variantHgvsName, variantType, assembly, chr,
                    chrStart, chrEnd, seqRef, seqVariant, mostSevere, variantReference, hasDisease, hasPhenotype);
            if (gene != null) {
                allele.setReference("gene", gene);
                gene.addToCollection("alleles", allele);
            }
        }
        System.out.println("size of alleles:  " + alleles.size());

        storeVariants();
        storeAlleles();
        storeGenes();
    }

    /**
     * @param g1
     * @param s1
     * @param g2
     * @param s2
     * @param algorithm
     * @param match
     * @param total
     * @param best
     * @param reverse
     * @throws ObjectStoreException
     */
    private Item processAlleles(String alleleId, String alleleSymbol, String alleleSynonym, String alleleType,
                                String variantId, String variantSymbol, String variantSynonym, String variantCrossRefs,
                                String variantHgvsName, String variantType, String assembly, String chr,
                                String chrStart, String chrEnd, String seqRef, String seqVariant, String mostSevere,
                                String variantReference, String hasDisease, String hasPhenotype)
            throws ObjectStoreException {

        Item allele = alleles.get(alleleId);

        if (allele == null) {

            allele = createItem("Allele");
            allele.setAttribute("featureType", "Allele");
            if (StringUtils.isNotEmpty(alleleId)) allele.setAttribute("alleleId", alleleId);
            if (StringUtils.isNotEmpty(alleleSymbol)) allele.setAttribute("alleleSymbol", alleleSymbol);
            //if (StringUtils.isNotEmpty(alleleSynonym)) allele.setAttribute("alleleSynonym", alleleSynonym);
            if (StringUtils.isNotEmpty(alleleType)) allele.setAttribute("alleleType", alleleType);
            //String refId = allele.getIdentifier();
                    /*if(pmrefNo != null ) {
                        Item publication = publications.get(pmrefNo);
                        if (publication == null) {
                            publication = createItem("Publication");
                            publication.setAttribute("pubMedId", pmid);
                            publications.put(pmrefNo, publication);
                        }
                        allele.addToCollection("publications", publication);
                    }
                    if(gene != null) {
                        allele.setReference("gene", gene);
                        gene.addToCollection("alleles", allele);
                    }*/
            alleles.put(alleleId, allele);
            alleleNames.put(alleleId, allele);
        } //allele

        if(alleleType.contains("with")) {

            Item variant = createItem("Variant");
            if (StringUtils.isNotEmpty(variantId)) variant.setAttribute("variantId", variantId);
            if (StringUtils.isNotEmpty(variantSymbol)) variant.setAttribute("variantSymbol", variantSymbol);
            //if (StringUtils.isNotEmpty(variantSynonym)) variant.setAttribute("variantSynonym", variantSynonym);
            //if (StringUtils.isNotEmpty(variantCrossRefs)) variant.setAttribute("variantCrossRefs", variantCrossRefs);
            if (StringUtils.isNotEmpty(variantHgvsName)) variant.setAttribute("VariantsHgvsNames", variantHgvsName);
            if (StringUtils.isNotEmpty(variantType)) variant.setAttribute("variantType", variantType);

            Item variantdetail = createItem("VariantDetails");
            if (StringUtils.isNotEmpty(assembly)) variantdetail.setAttribute("assembly", assembly);
            if (StringUtils.isNotEmpty(chr)) variantdetail.setAttribute("chr", chr);
            if (StringUtils.isNotEmpty(chrStart)) variantdetail.setAttribute("chrStartPosition", chrStart);
            if (StringUtils.isNotEmpty(chrEnd)) variantdetail.setAttribute("chrEndPosition", chrEnd);
            if (StringUtils.isNotEmpty(seqRef)) variantdetail.setAttribute("sequenceOfReference", seqRef);
            if (StringUtils.isNotEmpty(seqVariant)) variantdetail.setAttribute("sequenceOfVariant", seqVariant);
            if (StringUtils.isNotEmpty(mostSevere)) variantdetail.setAttribute("mostSevereConsequenceName", mostSevere);
            if (StringUtils.isNotEmpty(variantReference))
                variantdetail.setAttribute("variantInformationReference", variantReference);
            if (StringUtils.isNotEmpty(hasDisease)) variantdetail.setAttribute("hasDiseaseAnnotations", hasDisease);
            if (StringUtils.isNotEmpty(hasPhenotype))
                variantdetail.setAttribute("hasPhenotypeAnnotations", hasPhenotype);

            variant.addToCollection("variantdetails", variantdetail);
            variant.setReference("allele", allele);  //<---missed and wasted many hours!!??!!??!!

            try {
                store(variantdetail);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
            /*if(pmrefNo != null ) {
                Item publication = publications.get(pmrefNo);
                if (publication == null) {
                    publication = createItem("Publication");
                    publication.setAttribute("pubMedId", pmid);
                    publications.put(pmrefNo, publication);
                }
                allele.addToCollection("publications", publication);
            }*/

            allele.addToCollection("variants", variant);
            variants.put(variantId, variant);

        }

        return allele;
    }

    /**
     *
     * @param g
     * @param org
     * @return
     * @throws Exception
     */
    private Item getGene(String g, String org) throws ObjectStoreException {

        Item gene  = genes.get(g);
        if(gene == null) {
            //System.out.println("creating new  ..." + g);
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", g);
            gene.setReference("organism", org);
            genes.put(g, gene);
        }
        return gene; //.getIdentifier();
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
    /**
     *
     * @throws ObjectStoreException
     */

    private void storeAlleles() throws ObjectStoreException {
        for (Item allele : alleles.values()) {
            try {
                store(allele);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }

    /**
     *
     * @throws ObjectStoreException
     */

    private void storeVariants() throws ObjectStoreException {
        for (Item variant : variants.values()) {
            try {
                store(variant);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }

}
