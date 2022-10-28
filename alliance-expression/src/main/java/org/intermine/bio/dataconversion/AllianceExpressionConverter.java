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
    private Map<Item, Item> expannotations = new HashMap();
    private Map<String, Item> ontoTerms = new HashMap();
    private Map<String, Item> publications = new HashMap();

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
            if(count < 15 ) { count++; continue;}
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
        storeOntologyTerms();
        storeExpAnnotations();
        storeGenes();

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

        Item gene = getGene(g, o);

        if (gene == null) {
            return;
        }
        Item expAnnot = createItem("ExpressionAnnotation");
        if(StringUtils.isNotEmpty(loc) ) { expAnnot.setAttribute("location", loc);}
        if(StringUtils.isNotEmpty(stage) ) { expAnnot.setAttribute("stageTerm", stage);}
        if(StringUtils.isNotEmpty(source) ) { expAnnot.setAttribute("source", source);}
        if(StringUtils.isNotEmpty(sourceurl) ) {expAnnot.setAttribute("sourceURL", sourceurl);}
        if(StringUtils.isNotEmpty(assay) ) { expAnnot.setReference("assay", newOntologyTerm(assay));}
        if(StringUtils.isNotEmpty(cellular) ) { expAnnot.setReference("cellularcomponent", newOntologyTerm(cellular));}
        if(StringUtils.isNotEmpty(structure) ) {expAnnot.setReference("substructure", newOntologyTerm(structure));}
        if(StringUtils.isNotEmpty(anatomy) ) { expAnnot.setReference("anatomy", newOntologyTerm(anatomy));}
        if(StringUtils.isNotEmpty(ref) ) {
            String[] s = ref.split(",");
            for(int i=0; i < s.length; i++) {
                String pub = getPublication(s[i]);
                expAnnot.addToCollection("publications", pub);
            }
        }
        expannotations.put(gene, expAnnot);
        //expAnnot.setReference("gene", gene);
        gene.addToCollection("expressionAnnotations", expAnnot);

    }


    private String getPublication(String pubMedId)
            throws ObjectStoreException {
        Item item = publications.get(pubMedId);
        if (item == null) {
            item = createItem("Publication");
            if (StringUtils.isNotEmpty(pubMedId)) {
                item.setAttribute("pubMedId", pubMedId);
            }
            try {
                store(item);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
            publications.put(pubMedId, item);
        }
        return item.getIdentifier();
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
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", g);
            gene.setReference("organism", org);
        }
        String geneId = gene.getIdentifier();
        genes.put(g, gene);
        return gene;
    }

    private String newOntologyTerm(String identifier) throws ObjectStoreException {

        Item term = ontoTerms.get(identifier);
        if (term == null) {
            String[] ont = identifier.split(":");
            if(ont[0].startsWith("ZFA")){
                term = createItem("ZFATerm");
            }else if(ont[0].startsWith("EMAPA")){
                term = createItem("EMAPATerm");
            }else if(ont[0].startsWith("WBbt")){
                term = createItem("WBBTTerm");
            }else if(ont[0].startsWith("FBbt")){
                term = createItem("FBBTTerm");
            }else{
                term = createItem("OntologyTerm");
            }
            System.out.println(" term.. " + identifier);
            term.setAttribute("identifier", identifier);
        }
        ontoTerms.put(identifier, term);
        return term.getIdentifier();
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


    /**
     *
     * @throws ObjectStoreException
     */

    private void storeOntologyTerms() throws ObjectStoreException {
        for (Item ont : ontoTerms.values()) {
            try {
                store(ont);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }

}
