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
import java.util.*;

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
public class AllianceExpressionConverter extends BioFileConverter {

    private static final String DATASET_TITLE = "Alliance WT Expression data set";
    private static final String DATA_SOURCE_NAME = "Alliance WT Expression";
    private String licence;
    private Map<String, Item> genes = new HashMap();
    private Map<String, Item> expannotations = new HashMap();
    protected Map<String, String> ontoTerms = new LinkedHashMap<String, String>();

    /**
     * Construct a new AllianceGenesConverter.
     * @param database the database to read from
     * @param model the Model used by the object store we will write to with the ItemWriter
     * @param writer an ItemWriter used to handle Items created
     */
    public AllianceExpressionConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception, ObjectStoreException{

        /*Species SpeciesID  {1}     GeneID {2}  GeneSymbol      Location{4}        StageTerm{5}
        AssayID {6} AssayTermName
        CellularComponentID {8}     CellularComponentTerm   CellularComponentQualifierIDs CellularComponentQualifierTermNames
        SubStructureID {12}  SubStructureName        SubStructureQualifierIDs        SubStructureQualifierTermNames
         AnatomyTermID   {16} AnatomyTermName AnatomyTermQualifierIDs AnatomyTermQualifierTermNames
         SourceURL{20}   Source{21}  Reference {22}
         */
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int count = 0;
        System.out.println("Processing Expression...");
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
           // if(count < 17 ) { count++; continue;}
            String origspecies = line[1].trim();
            String species = origspecies.replace("NCBITaxon:","");
            String org = getOrganism(species);
            String geneId = line[2];
            String location = line[4];
            String stageterm = line[5];
            String assayID = line[6];
            String cellularcomponentID = line[8];
            String substructureID = line[12];
            String anatomyID = line[16];
            String sourceUrl = line[20];
            String source = line[21];
            String reference = line[22];
            processExpressionAnnotation(geneId, org, location, stageterm, assayID, cellularcomponentID, substructureID,
                    anatomyID, sourceUrl, source, reference);
        }
        System.out.println("size of expannots genes:  " + genes.size());
        storeGenes();
        storeExpAnnotations();

    }

    /**
     *
     * @param g
     * @param o
     * @param loc
     * @param stage
     * @param assay
     * @param cellular
     * @param structure
     * @param anatomy
     * @param sourceurl
     * @param ref
     * @throws ObjectStoreException
     */
    private void processExpressionAnnotation(String g, String o, String loc, String stage, String assay, String cellular, String structure,
                                   String anatomy, String sourceurl, String source, String ref)
            throws ObjectStoreException {

        String gene = getGene(g, o);

        if (gene == null) {
            return;
        }
        Item expAnnot = createItem("ExpressionAnnotation");
        expAnnot.setReference("gene", gene);
        expAnnot.setAttribute("location", loc);
        expAnnot.setAttribute("stageTerm", stage);
        expAnnot.setAttribute("source", source);
        expAnnot.setAttribute("sourceURL", sourceurl);
        expAnnot.setReference("assay", newOntologyTerm(assay));
        expAnnot.setReference("cellularcomponent", newOntologyTerm(cellular));
        expAnnot.setReference("substructure", newOntologyTerm(structure));
        expAnnot.setReference("anatomy", newOntologyTerm(anatomy));
        expAnnot.setAttribute("pubmedID", ref);
        expannotations.put(gene, expAnnot);
        expAnnot.setReference("gene", gene);
        gene.addToCollection("expressionAnnotations", expAnnot);

    }
    /**
     *
     * @param g
     * @param org
     * @return
     * @throws Exception
     */
    private String getGene(String g, String org) throws ObjectStoreException {

        Item gene  = genes.get(g);
        if(gene == null) {
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", g);
            gene.setReference("organism", org);
        }
        String geneId = gene.getIdentifier();
        genes.put(geneId, gene);
        return geneId;
    }

    private String newOntologyTerm(String identifier) throws ObjectStoreException {
        if (identifier == null) {
            return null;
        }
        String termIdentifier = ontoTerms.get(identifier);
        if (termIdentifier == null) {
            Item item = createItem(termClassName);
            item.setAttribute("identifier", identifier);
            store(item);
            termIdentifier = item.getIdentifier();
            ontoTerms.put(identifier, termIdentifier);
        }
        return termIdentifier;
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

    private void storeExpAnnotations() throws ObjectStoreException {
        for (Item exp : expannotations.values()) {
            try {
                store(exp);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }

}
