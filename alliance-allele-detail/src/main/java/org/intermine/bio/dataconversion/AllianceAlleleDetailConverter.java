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
import java.util.Iterator;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads allele-detail.tsv (emitted by scripts/fetch_allele_detail.py from the
 * Alliance /allele/{id} endpoint) and emits partial Allele items carrying
 * three additional attributes: alterationType, apiCategory, apiCrossReference.
 *
 * The InterMine integration engine merges these partial Allele items into
 * those produced by alliance-alleles via the Allele.key_alleleid integration
 * key (alleleId). Order of source ingestion is decided by project.xml; both
 * sources work whether this one runs before or after alliance-alleles.
 *
 * @author
 */
public class AllianceAlleleDetailConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceAlleleDetailConverter.class);

    private static final String DATASET_TITLE = "Alliance Allele Detail";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Must match scripts/fetch_allele_detail.py COLUMNS
    private static final int COL_ALLELE_ID         = 0;
    private static final int COL_ALTERATION_TYPE   = 1;
    private static final int COL_API_CATEGORY      = 2;
    private static final int COL_API_CROSSREF      = 3;
    private static final int MIN_COLUMNS           = 4;

    public AllianceAlleleDetailConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance allele detail...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            String alleleId = line[COL_ALLELE_ID].trim();
            if (alleleId.isEmpty() || "alleleId".equals(alleleId)) {
                continue;
            }
            Item allele = createItem("Allele");
            allele.setAttribute("alleleId", alleleId);
            setIfPresent(allele, "alterationType", line[COL_ALTERATION_TYPE]);
            setIfPresent(allele, "apiCategory", line[COL_API_CATEGORY]);
            setIfPresent(allele, "apiCrossReference", line[COL_API_CROSSREF]);
            store(allele);
            rows++;
        }
        LOG.info("AlleleDetail: emitted " + rows + " partial Allele items");
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
