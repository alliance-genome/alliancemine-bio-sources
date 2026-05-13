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
 * Reads agms.tsv (emitted by scripts/fetch_agms.py from the Alliance FMS
 * AGM bulk files) and produces DiseaseModel items for each Affected
 * Genomic Model. Disease + Gene references are populated by the sibling
 * alliance-disease-models source (which carries the gene/disease links);
 * this module fills in the model identity (modelId, modelName,
 * modelSubtype, organism).
 */
public class AllianceAgmsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceAgmsConverter.class);

    private static final String DATASET_TITLE = "Alliance Affected Genomic Models";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Must match scripts/fetch_agms.py COLUMNS
    private static final int COL_MODEL_ID      = 0;
    private static final int COL_MODEL_NAME    = 1;
    private static final int COL_MODEL_SUBTYPE = 2;
    private static final int COL_TAXON         = 3;
    private static final int COL_DATA_PROVIDER = 4;
    private static final int MIN_COLUMNS       = 5;

    public AllianceAgmsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance AGMs...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("modelId".equals(line[COL_MODEL_ID])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("AGMs: emitted " + rows + " DiseaseModel items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String modelId = line[COL_MODEL_ID].trim();
        if (modelId.isEmpty()) {
            return;
        }
        Item model = createItem("DiseaseModel");
        model.setAttribute("modelId", modelId);
        setIfPresent(model, "modelName", line[COL_MODEL_NAME]);
        setIfPresent(model, "modelSubtype", line[COL_MODEL_SUBTYPE]);
        setIfPresent(model, "dataProvider", line[COL_DATA_PROVIDER]);
        store(model);
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
