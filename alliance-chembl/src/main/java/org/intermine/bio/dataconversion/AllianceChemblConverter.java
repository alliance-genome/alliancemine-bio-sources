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

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads chembl-targets.tsv (emitted by scripts/fetch_chembl.py) and
 * produces CrossReference items on the Protein class, source=ChEMBL.
 * Allows AllianceMine queries to pivot from any UniProt-keyed Protein
 * to its ChEMBL target page (mechanism/activity/drug indications).
 */
public class AllianceChemblConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceChemblConverter.class);

    private static final String DATASET_TITLE = "ChEMBL Target Mapping";
    private static final String DATA_SOURCE_NAME = "ChEMBL";

    private static final int COL_ACC      = 0;
    private static final int COL_CHEMBL   = 1;
    private static final int COL_NAME     = 2;
    private static final int COL_TYPE     = 3;
    private static final int COL_URL      = 4;
    private static final int MIN_COLUMNS  = 5;

    private final Map<String, String> proteinsByAcc = new HashMap<String, String>();
    private String chemblSourceRef;

    public AllianceChemblConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing ChEMBL target → UniProt cross-references...");
        if (chemblSourceRef == null) {
            Item ds = createItem("DataSource");
            ds.setAttribute("name", "ChEMBL");
            store(ds);
            chemblSourceRef = ds.getIdentifier();
        }
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("uniProtAccession".equals(line[COL_ACC])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("ChEMBL: emitted " + rows + " CrossReference items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String acc = line[COL_ACC].trim();
        String chembl = line[COL_CHEMBL].trim();
        if (acc.isEmpty() || chembl.isEmpty()) {
            return;
        }
        Item xr = createItem("CrossReference");
        xr.setAttribute("identifier", chembl);
        xr.setAttribute("url", line[COL_URL].trim());
        xr.setAttribute("dbxreftype", line[COL_TYPE].trim());
        xr.setReference("subject", getProtein(acc));
        xr.setReference("source", chemblSourceRef);
        store(xr);
    }

    private String getProtein(String acc) throws ObjectStoreException {
        String ref = proteinsByAcc.get(acc);
        if (ref != null) {
            return ref;
        }
        Item p = createItem("Protein");
        p.setAttribute("primaryAccession", acc);
        store(p);
        proteinsByAcc.put(acc, p.getIdentifier());
        return p.getIdentifier();
    }
}
