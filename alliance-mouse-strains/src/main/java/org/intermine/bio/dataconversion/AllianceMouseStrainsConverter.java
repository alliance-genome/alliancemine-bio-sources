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
 * Reads mouse-strains.tsv (emitted by scripts/fetch_mousemine_strains.py from
 * MouseMine via PathQuery REST) and produces Strain items plus their carried
 * Allele items. Each TSV row is one (strain, allele) pair; rows for the same
 * strain are folded into a single Strain Item with a populated alleles
 * collection.
 *
 * @author
 */
public class AllianceMouseStrainsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceMouseStrainsConverter.class);

    private static final String DATASET_TITLE = "MouseMine Strains";
    private static final String DATA_SOURCE_NAME = "MouseMine";
    private static final String MOUSE_TAXON = "10090";

    // Must match scripts/fetch_mousemine_strains.py COLUMNS
    private static final int COL_STRAIN_ID         = 0;
    private static final int COL_STRAIN_NAME       = 1;
    private static final int COL_ATTRIBUTE_STRING  = 2;
    private static final int COL_ALLELE_ID         = 3;
    private static final int MIN_COLUMNS           = 4;

    private final Map<String, Item> strains = new HashMap<String, Item>();
    private final Map<String, String> alleles = new HashMap<String, String>();
    private String organismRef = null;

    public AllianceMouseStrainsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        LOG.info("Processing MouseMine strains...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("strainId".equals(line[COL_STRAIN_ID])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        for (Item strain : strains.values()) {
            store(strain);
        }
        LOG.info("MouseStrains: processed " + rows + " (strain,allele) rows -> "
                + strains.size() + " distinct Strain items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String strainId = line[COL_STRAIN_ID].trim();
        if (strainId.isEmpty()) {
            return;
        }
        Item strain = strains.get(strainId);
        if (strain == null) {
            strain = createItem("Strain");
            strain.setAttribute("primaryIdentifier", strainId);
            setIfPresent(strain, "name", line[COL_STRAIN_NAME]);
            setIfPresent(strain, "attributeString", line[COL_ATTRIBUTE_STRING]);
            strain.setReference("organism", getMouseOrganism());
            strains.put(strainId, strain);
        }
        String alleleId = line[COL_ALLELE_ID].trim();
        if (!alleleId.isEmpty()) {
            strain.addToCollection("alleles", getAllele(alleleId));
        }
    }

    private String getMouseOrganism() throws ObjectStoreException {
        if (organismRef != null) {
            return organismRef;
        }
        organismRef = getOrganism(MOUSE_TAXON);
        return organismRef;
    }

    private String getAllele(String alleleId) throws ObjectStoreException {
        String ref = alleles.get(alleleId);
        if (ref != null) {
            return ref;
        }
        Item a = createItem("Allele");
        a.setAttribute("alleleId", alleleId);
        store(a);
        alleles.put(alleleId, a.getIdentifier());
        return a.getIdentifier();
    }

    private static void setIfPresent(Item item, String attrName, String value) {
        if (value != null) {
            String v = value.trim();
            if (!v.isEmpty() && !"-".equals(v)) {
                item.setAttribute(attrName, v);
            }
        }
    }
}
