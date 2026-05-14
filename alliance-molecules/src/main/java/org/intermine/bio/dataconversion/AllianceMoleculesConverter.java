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

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads molecules.tsv (emitted by scripts/fetch_molecules.py from
 * Alliance FMS MOLECULE bulk) and produces Molecule items.
 */
public class AllianceMoleculesConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceMoleculesConverter.class);

    private static final String DATASET_TITLE = "Alliance Small Molecules";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    private static final int COL_ID       = 0;
    private static final int COL_CHEBI    = 1;
    private static final int COL_NAME     = 2;
    private static final int COL_IUPAC    = 3;
    private static final int COL_SMILES   = 4;
    private static final int COL_INCHI    = 5;
    private static final int COL_INCHIKEY = 6;
    private static final int COL_PROVIDER = 7;
    private static final int MIN_COLUMNS  = 8;

    public AllianceMoleculesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance molecules...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("primaryIdentifier".equals(line[COL_ID])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("Molecules: emitted " + rows + " Molecule items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String pid = line[COL_ID].trim();
        if (pid.isEmpty()) {
            return;
        }
        Item m = createItem("Molecule");
        m.setAttribute("primaryIdentifier", pid);
        setIfPresent(m, "chebiId", line[COL_CHEBI]);
        setIfPresent(m, "name", line[COL_NAME]);
        setIfPresent(m, "iupacName", line[COL_IUPAC]);
        setIfPresent(m, "smiles", line[COL_SMILES]);
        setIfPresent(m, "inchi", line[COL_INCHI]);
        setIfPresent(m, "inchiKey", line[COL_INCHIKEY]);
        setIfPresent(m, "dataProvider", line[COL_PROVIDER]);
        store(m);
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
