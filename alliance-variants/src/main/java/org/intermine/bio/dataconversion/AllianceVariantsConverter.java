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
 * Reads variants.tsv (emitted by scripts/fetch_variants.py from the
 * Alliance VARIANT-ALLELE-JSON FMS bulk files) and produces Variant
 * items plus one VariantConsequence per consequence name. Allele +
 * Gene references use enrichment-merge against records produced by
 * alliance-alleles / alliance-genes via the standard primary-identifier
 * + organism keys.
 */
public class AllianceVariantsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceVariantsConverter.class);

    private static final String DATASET_TITLE = "Alliance Variants";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Must match scripts/fetch_variants.py COLUMNS
    private static final int COL_VARIANT_ID    = 0;
    private static final int COL_VARIANT_SYM   = 1;
    private static final int COL_VAR_TYPE      = 2;
    private static final int COL_VAR_TYPE_ID   = 3;
    private static final int COL_HGVS          = 4;
    private static final int COL_CHROMOSOME    = 5;
    private static final int COL_START         = 6;
    private static final int COL_END           = 7;
    private static final int COL_REF_SEQ       = 8;
    private static final int COL_VAR_SEQ       = 9;
    private static final int COL_CONSEQUENCES  = 10;
    private static final int COL_ALLELE_ID     = 11;
    private static final int COL_GENE_ID       = 12;
    private static final int COL_TAXON         = 13;
    private static final int COL_DATA_PROVIDER = 14;
    private static final int MIN_COLUMNS       = 15;

    private final Map<String, String> alleles = new HashMap<String, String>();
    private final Map<String, String> genes = new HashMap<String, String>();
    private final Map<String, String> soTerms = new HashMap<String, String>();

    public AllianceVariantsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance variants...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        int consequences = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("variantId".equals(line[COL_VARIANT_ID])) {
                continue;
            }
            consequences += processRow(line);
            rows++;
        }
        LOG.info("Variants: emitted " + rows + " Variant items + " + consequences + " VariantConsequence items");
    }

    private int processRow(String[] line) throws ObjectStoreException {
        String variantId = line[COL_VARIANT_ID].trim();
        if (variantId.isEmpty()) {
            return 0;
        }

        // Variant carries identifier + symbol + type; coordinate / sequence
        // facets live on a VariantDetails child to match the existing
        // genomic-additions schema. VariantConsequence captures the
        // most-severe-consequence rows produced by VEP.
        Item variant = createItem("Variant");
        variant.setAttribute("primaryIdentifier", variantId);
        variant.setAttribute("variantId", variantId);
        setIfPresent(variant, "variantSymbol", line[COL_VARIANT_SYM]);
        setIfPresent(variant, "variantType", line[COL_VAR_TYPE]);
        setIfPresent(variant, "variantsTypeId", line[COL_VAR_TYPE_ID]);
        setIfPresent(variant, "VariantsHgvsNames", line[COL_HGVS]);

        String alleleId = line[COL_ALLELE_ID].trim();
        String taxon = line[COL_TAXON].trim();
        if (!alleleId.isEmpty()) {
            variant.setReference("allele", getAllele(alleleId, taxon));
        }
        String geneCsv = line[COL_GENE_ID].trim();
        String firstGene = "";
        if (!geneCsv.isEmpty()) {
            firstGene = geneCsv.split("\\|", 2)[0].trim();
            if (!firstGene.isEmpty()) {
                variant.setReference("affectedGene", getGene(firstGene, taxon));
            }
        }
        if (StringUtils.isNotEmpty(taxon)) {
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            variant.setReference("organism", getOrganism(taxonId));
        }
        store(variant);

        // Coordinates / sequence go on a VariantDetails child if any are populated.
        boolean hasCoords = !line[COL_CHROMOSOME].trim().isEmpty()
            || !line[COL_START].trim().isEmpty()
            || !line[COL_REF_SEQ].trim().isEmpty();
        if (hasCoords) {
            Item details = createItem("VariantDetails");
            setIfPresent(details, "chr", line[COL_CHROMOSOME]);
            setIntIfPresent(details, "chrStartPosition", line[COL_START]);
            setIntIfPresent(details, "chrEndPosition", line[COL_END]);
            setIfPresent(details, "sequenceOfReference", line[COL_REF_SEQ]);
            setIfPresent(details, "sequenceOfVariant", line[COL_VAR_SEQ]);
            store(details);
            variant.addToCollection("variantdetails", details);
        }

        // One VariantConsequence per consequence-name token.
        int emitted = 0;
        String consequenceCsv = line[COL_CONSEQUENCES].trim();
        if (!consequenceCsv.isEmpty()) {
            for (String token : consequenceCsv.split("\\|")) {
                String code = token.trim();
                if (code.isEmpty()) {
                    continue;
                }
                Item cons = createItem("VariantConsequence");
                cons.setAttribute("consequenceCode", code);
                cons.setReference("variant", variant);
                if (!firstGene.isEmpty()) {
                    cons.setReference("gene", getGene(firstGene, taxon));
                }
                store(cons);
                emitted++;
            }
        }
        return emitted;
    }

    private String getAllele(String primaryId, String taxon) throws ObjectStoreException {
        String ref = alleles.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item item = createItem("Allele");
        item.setAttribute("primaryIdentifier", primaryId);
        if (StringUtils.isNotEmpty(taxon)) {
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            item.setReference("organism", getOrganism(taxonId));
        }
        store(item);
        alleles.put(primaryId, item.getIdentifier());
        return item.getIdentifier();
    }

    private String getGene(String primaryId, String taxon) throws ObjectStoreException {
        String ref = genes.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item item = createItem("Gene");
        item.setAttribute("primaryIdentifier", primaryId);
        if (StringUtils.isNotEmpty(taxon)) {
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            item.setReference("organism", getOrganism(taxonId));
        }
        store(item);
        genes.put(primaryId, item.getIdentifier());
        return item.getIdentifier();
    }

    private static void setIfPresent(Item item, String attrName, String value) {
        if (value != null) {
            String v = value.trim();
            if (!v.isEmpty() && !"-".equals(v)) {
                item.setAttribute(attrName, v);
            }
        }
    }

    private static void setIntIfPresent(Item item, String attrName, String value) {
        if (value == null) {
            return;
        }
        String v = value.trim();
        if (v.isEmpty()) {
            return;
        }
        try {
            Integer.parseInt(v);
            item.setAttribute(attrName, v);
        } catch (NumberFormatException nfe) {
            // skip non-numeric
        }
    }
}
