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
 * Reads htp-datasets.tsv + htp-samples.tsv (emitted by
 * scripts/fetch_htp.py from FMS HTPDATASET + HTPDATASAMPLE bulk files)
 * and produces DataSet items for parent studies, Sample items for the
 * individual sequencing/array samples, and SampleCharacteristic items
 * for the per-sample anatomy/stage/sex/assay facets.
 *
 * The converter dispatches on the first column header of each input file:
 *   "datasetId" -> processDatasetRow()
 *   "sampleId"  -> processSampleRow()
 */
public class AllianceHtpConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceHtpConverter.class);

    private static final String DATASET_TITLE = "Alliance High-Throughput Experiments";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Dataset TSV column indices
    private static final int DS_ID         = 0;
    private static final int DS_TITLE      = 1;
    private static final int DS_SUMMARY    = 2;
    private static final int DS_CATEGORIES = 3;
    private static final int DS_GEO        = 4;
    private static final int DS_PMIDS      = 5;
    private static final int DS_DATE       = 6;
    private static final int DS_PROVIDER   = 7;
    private static final int DS_MIN_COLS   = 8;

    // Sample TSV column indices
    private static final int S_ID            = 0;
    private static final int S_TITLE         = 1;
    private static final int S_ABUNDANCE     = 2;
    private static final int S_TYPE          = 3;
    private static final int S_STAGE_ID      = 4;
    private static final int S_STAGE_NAME    = 5;
    private static final int S_ANATOMY_ID    = 6;
    private static final int S_ANATOMY_STMT  = 7;
    private static final int S_BIOSAMPLE     = 8;
    private static final int S_SEX           = 9;
    private static final int S_ASSAY         = 10;
    private static final int S_ASSEMBLY      = 11;
    private static final int S_DATASET_IDS   = 12;
    private static final int S_TAXON         = 13;
    private static final int S_PROVIDER      = 14;
    private static final int S_MIN_COLS      = 15;

    private final Map<String, String> datasets = new HashMap<String, String>();
    private final Map<String, String> samples = new HashMap<String, String>();
    private final Map<String, String> publications = new HashMap<String, String>();
    private final Map<String, String> ontologyTerms = new HashMap<String, String>();

    public AllianceHtpConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        String[] header = null;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length == 0) {
                continue;
            }
            if (header == null && ("datasetId".equals(line[0]) || "sampleId".equals(line[0]))) {
                header = line;
                continue;
            }
            if (header == null) {
                // Skip pre-header comment rows.
                continue;
            }
            if ("datasetId".equals(header[0])) {
                if (line.length >= DS_MIN_COLS) {
                    processDatasetRow(line);
                }
            } else if ("sampleId".equals(header[0])) {
                if (line.length >= S_MIN_COLS) {
                    processSampleRow(line);
                }
            }
        }
    }

    private void processDatasetRow(String[] line) throws ObjectStoreException {
        String id = line[DS_ID].trim();
        if (id.isEmpty()) {
            return;
        }
        if (datasets.containsKey(id)) {
            return;
        }
        Item ds = createItem("DataSet");
        // Use the unique dataset CURIE as `name` — DataSet's primary key is
        // `name`, so collisions on title-based names abort the integrate. The
        // human-readable title becomes the description prefix.
        ds.setAttribute("name", id);
        ds.setAttribute("version", id);
        String title = line[DS_TITLE].trim();
        String summary = line[DS_SUMMARY].trim();
        if (!title.isEmpty() || !summary.isEmpty()) {
            String desc = title.isEmpty()
                ? summary
                : (summary.isEmpty() ? title : title + " — " + summary);
            ds.setAttribute("description", desc);
        }
        store(ds);
        datasets.put(id, ds.getIdentifier());

        // DataSet has a single `publication` reference (not a collection).
        // Use the first PMID; remaining PMIDs are reachable via Publication
        // primary-key dedup at the model level (Publication.key_pubmedid).
        String pmidCsv = line[DS_PMIDS].trim();
        if (!pmidCsv.isEmpty()) {
            String firstPmid = pmidCsv.split("\\|", 2)[0].trim();
            String pubRef = getPublication(firstPmid);
            if (pubRef != null) {
                ds.setReference("publication", pubRef);
            }
            // Materialise remaining PMIDs as standalone Publication items so
            // they appear in the warehouse even without a back-reference.
            for (int i = 1; i < pmidCsv.split("\\|").length; i++) {
                getPublication(pmidCsv.split("\\|")[i].trim());
            }
        }
    }

    private void processSampleRow(String[] line) throws ObjectStoreException {
        String id = line[S_ID].trim();
        if (id.isEmpty() || samples.containsKey(id)) {
            return;
        }
        Item sample = createItem("Sample");
        sample.setAttribute("primaryIdentifier", id);
        setIfPresent(sample, "primaryCharacteristic", line[S_TITLE]);
        setIfPresent(sample, "primaryCharacteristicType", line[S_TYPE]);
        setIfPresent(sample, "description", line[S_ABUNDANCE]);
        setIfPresent(sample, "materialType", line[S_SEX]);

        String taxon = line[S_TAXON].trim();
        if (StringUtils.isNotEmpty(taxon)) {
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            sample.setReference("organism", getOrganism(taxonId));
        }

        // Link sample → parent dataset(s).
        String dsCsv = line[S_DATASET_IDS].trim();
        if (!dsCsv.isEmpty()) {
            for (String token : dsCsv.split("\\|")) {
                String dsId = token.trim();
                String dsRef = datasets.get(dsId);
                if (dsRef == null) {
                    // Sample arrives before dataset row — create a shadow dataset
                    // that the dataset pass will merge into when keys collide.
                    Item shadow = createItem("DataSet");
                    shadow.setAttribute("name", dsId);
                    shadow.setAttribute("version", dsId);
                    store(shadow);
                    dsRef = shadow.getIdentifier();
                    datasets.put(dsId, dsRef);
                }
                sample.addToCollection("dataSets", dsRef);
            }
        }
        store(sample);
        samples.put(id, sample.getIdentifier());

        // Per-sample characteristics. Each is a SampleCharacteristic item with
        // type + value + optional ontologyTerm reference.
        addCharacteristic(sample, "stage", line[S_STAGE_NAME], line[S_STAGE_ID]);
        addCharacteristic(sample, "anatomy", line[S_ANATOMY_STMT], line[S_ANATOMY_ID]);
        addCharacteristic(sample, "sex", line[S_SEX], "");
        addCharacteristic(sample, "assay", line[S_ASSAY], line[S_ASSAY]);
        addCharacteristic(sample, "assembly", line[S_ASSEMBLY], "");
        addCharacteristic(sample, "biosample", line[S_BIOSAMPLE], "");
    }

    private void addCharacteristic(Item sample, String type, String value, String termId)
        throws ObjectStoreException {
        if (value == null) {
            return;
        }
        String v = value.trim();
        if (v.isEmpty() || "-".equals(v)) {
            return;
        }
        Item c = createItem("SampleCharacteristic");
        c.setAttribute("type", type);
        c.setAttribute("value", v);
        if (StringUtils.isNotEmpty(termId)) {
            c.setReference("ontologyTerm", getOntologyTerm(termId.trim()));
        }
        store(c);
        sample.addToCollection("characteristics", c);
    }

    private String getOntologyTerm(String identifier) throws ObjectStoreException {
        String ref = ontologyTerms.get(identifier);
        if (ref != null) {
            return ref;
        }
        Item term = createItem("OntologyTerm");
        term.setAttribute("identifier", identifier);
        store(term);
        ontologyTerms.put(identifier, term.getIdentifier());
        return term.getIdentifier();
    }

    private String getPublication(String pubToken) throws ObjectStoreException {
        if (pubToken == null) {
            return null;
        }
        String p = pubToken.trim();
        if (p.isEmpty()) {
            return null;
        }
        String ref = publications.get(p);
        if (ref != null) {
            return ref;
        }
        Item pub = createItem("Publication");
        if (p.startsWith("PMID:")) {
            pub.setAttribute("pubMedId", p.substring("PMID:".length()));
        } else {
            pub.setAttribute("pubXrefId", p);
        }
        store(pub);
        publications.put(p, pub.getIdentifier());
        return pub.getIdentifier();
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
