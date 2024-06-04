package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2010 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.sql.Connection;
import java.util.ArrayList;
import java.sql.Array;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections.keyvalue.MultiKey;
import org.apache.commons.lang.StringUtils; 
//import org.biojava.bio.program.homologene.OrthoPairSet.Iterator;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.sql.Database;
import org.intermine.xml.full.Item;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Iterator;

/**
 * Converts results sets into intermine objects
 * 
 * @author Julie Sullivan
 */
public class SgdConverter extends BioDBConverter {

	// private static final Logger LOG = Logger.getLogger(SgdConverter.class);
	private static final String DATASET_TITLE = "SGD data set";
	private static final String DATA_SOURCE_NAME = "SGD";
	private Map<String, String> chromosomes = new HashMap();
	private Map<String, String> plasmids = new HashMap();
	private Map<String, String> sequences = new HashMap();
	private Map<String, Item> interactions = new HashMap();
	private Map<String, String> interactionitems = new HashMap();
	private Map<String, String> interactionterms = new HashMap<String, String>();
	private Map<MultiKey, Item> interactionsnew = new HashMap<MultiKey, Item>();
	private final Map<String, Item> ecoMap = new HashMap<String, Item>(); //regulation data
	private Map<String, String> literatureTopics = new HashMap();
	private Map<String, Item> genes = new HashMap();
	private Map<String, Item> hgncgenes = new HashMap();
	private Map<String, Item> alleles = new HashMap();
	private Map<String, Item> transcripts = new HashMap();
	private Map<String, Item> proteins = new HashMap();
	private Map<String, Item> genesName = new HashMap();
	private Map<String, String> genesChromosomes = new HashMap();
	private Map<String, Item> allelesName = new HashMap();
	private Map<String, String> genesAliases = new HashMap();
	private Map<String, String> alleleAliases = new HashMap();
	private Map<String, String> synonyms = new HashMap();
	private Map<String, Item> publications = new HashMap();
	private Map<String, Item> interactiontype = new HashMap();
	private Map<String, Item> interactiondetail = new HashMap();
	private Map<String, Item> experimenttype = new HashMap();
	private Map<String, Item> interactiondetectionmethods = new HashMap();
	private Map<String, Item> pathways = new HashMap();
	private Map<String, Item> pathwaysummarys = new HashMap();
	private Map<String, Item> phenotypes = new HashMap();
	private Map<String, HashMap<String, String>> phenotypeannots = new HashMap();
	private Map<String, String> datasources = new HashMap();
	private Map<String, Item> proteinDomains = new HashMap<String, Item>();
	private Map<String, Item> interproDomains = new HashMap<String, Item>();
	private Map<String, Item> proteinAbundance = new HashMap<String, Item>();
	private Map<String, Item> regulationSummary = new HashMap<String, Item>();

	private static final String TAXON_ID = "559292"; //4932
	private static final String H_TAXON_ID = "9606";
	private Item organism;
	private Map<String, String> featureMap = new HashMap();
	private static final boolean LOCAL = false;
	private String licence;


	private static final SgdProcessor PROCESSOR = new SgdProcessor();

	/**
	 * Construct a new SgdConverter.
	 * 
	 * @param database
	 *            the database to read from
	 * @param model
	 *            the Model used by the object store we will write to with the
	 *            ItemWriter
	 * @param writer
	 *            an ItemWriter used to handle Items created
	 * @throws ObjectStoreException
	 *             if organism can't be stored
	 */
	public SgdConverter(Database database, Model model, ItemWriter writer)
			throws ObjectStoreException {
		super(database, model, writer, DATA_SOURCE_NAME, DATASET_TITLE);
		organism = createItem("Organism");
		organism.setAttribute("taxonId", TAXON_ID);
		organism.setAttribute("genus", "Saccharomyces");
		organism.setAttribute("species", "cerevisiae");
		organism.setAttribute("name", "Saccharomyces cerevisiae");
		organism.setAttribute("shortName", "S. cerevisiae");
		store(organism);
	}

	/**
	 * {@inheritDoc}
	 */
	public void process() throws Exception {

		// a database has been initialized from properties starting with db.sgd
		Connection connection = getDatabase().getConnection();

		processChromosomeSequences(connection);
		processGenes(connection);
		processNISS(connection);
		processAliases(connection);
		processCrossReferences(connection);
		processGeneLocations(connection);
		processGeneChildrenLocations(connection);
		processProteins(connection);
		processAlleles(connection);
		processAlleleAliases(connection);
		processTranscripts(connection);
		processAllPubs(connection);
		processPubsWithFeatures(connection);
		
		processProteinAbundance(connection);
		processProteinHalfLife(connection);
		processProteinDomains(connection);
		processProteinModifications(connection);
		processProteinInfo(connection);
		
		processParalogs(connection);
		processFunctionSummary(connection);
		processRegulation(connection);
		processRegulationSummary(connection);

		processAllPathways(connection);
		processGenePathways(connection);
		storePathways();
		storeAlleles();
		storeTranscripts();

		if(!LOCAL) {

			 processPhysicalInteractions(connection);
			 processGeneticInteractions(connection);
			 processGeneticInteractionsWithAlleles(connection);
			 storeInteractionTypes();
			 storeInteractionExperiments();
			 //storeInteractionDetails(); <--keep commented
			 storeInteractions();

			 processPhenotypes(connection);
			 processPhenotypeSummary(connection);
			 storePhenotypes();

		}
		storePublications();
		storeRegulationSummaries();
		storeGenes();
		storeProteins();

	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processGenes(Connection connection) throws SQLException,
	ObjectStoreException {

		System.out.println("Processing Genes...");
		ResultSet res = PROCESSOR.getChromosomalFeatureResults(connection);

		while (res.next()) {

			String featureNo = res.getString("dbentity_id");
			if (genes.get(featureNo) == null) {
				//~~~ gene ~~~
				String primaryIdentifier = res.getString("sgdid");
				String secondaryIdentifier = res.getString("systematic_name");
				String symbol = res.getString("gene_name");
				String name = res.getString("name_description");
				String feature_type = res.getString("feature_type");
				String headline = res.getString("headline");
				String description = res.getString("description");
				String qualifier = res.getString("qualifier");
				//String feat_attribute = res.getString("feat_attribute");
				String status = res.getString("dbentity_status");

				Item item = null;
				if (feature_type.equalsIgnoreCase("ORF")) {
					item = createItem("ORF");
				} else if (feature_type.equalsIgnoreCase("pseudogene")) {
					item = createItem("Pseudogene");
				} else if (feature_type
						.equalsIgnoreCase("transposable element gene")) {
					item = createItem("TransposableElementGene");
				} else if (feature_type
						.equalsIgnoreCase("not physically mapped")) {
					item = createItem("NotPhysicallyMapped");
				} else if (feature_type
						.equalsIgnoreCase("long terminal repeat")) {
					item = createItem("LongTerminalRepeat");
				} else if (feature_type.equalsIgnoreCase("ARS")) {
					item = createItem("ARS");
				} else if (feature_type.equalsIgnoreCase("tRNA gene")) {
					item = createItem("TRNAGene");
				} else if (feature_type.equalsIgnoreCase("snoRNA gene")) {
					item = createItem("SnoRNAGene");
				} else if (feature_type.equalsIgnoreCase("not in systematic sequence of S288C")) {
					item = createItem("NotInSystematicSequenceOfS288C");
				} else if (feature_type.equalsIgnoreCase("LTR retrotransposon")) {
					item = createItem("Retrotransposon");
				} else if (feature_type.equalsIgnoreCase("X element combinatorial repeats")) {
					item = createItem("XElementCombinatorialRepeat");
				} else if (feature_type.equalsIgnoreCase("X element")) {
					item = createItem("XElement");
				} else if (feature_type.equalsIgnoreCase("telomere")) {
					item = createItem("Telomere");
				} else if (feature_type.equalsIgnoreCase("telomeric repeat")) {
					item = createItem("TelomericRepeat");
				} else if (feature_type.equalsIgnoreCase("rRNA gene")) {
					item = createItem("RRNAGene");
				} else if (feature_type.equalsIgnoreCase("Y prime element")) {
					item = createItem("YPrimeElement");
				} else if (feature_type.equalsIgnoreCase("centromere")) {
					item = createItem("Centromere");
				} else if (feature_type.equalsIgnoreCase("ncRNA gene")) {
					item = createItem("NcRNAGene");
				} else if (feature_type.equalsIgnoreCase("snRNA gene")) {
					item = createItem("SnRNAGene");
				}else if (feature_type.equalsIgnoreCase("blocked reading frame")) {
					item = createItem("BlockedReadingFrame");
				}else if (feature_type.equalsIgnoreCase("origin of replication")) {
					item = createItem("OriginOfReplication");
				}else if (feature_type.equalsIgnoreCase("matrix attachment site")) {
					item = createItem("MatrixAttachmentSite");
				}else if (feature_type.equalsIgnoreCase("telomerase RNA gene")) {
					item = createItem("TelomeraseRNAGene");
				}else if (feature_type.equalsIgnoreCase("gene group")) {
					item = createItem("GeneGroup");
				}else if (feature_type.equalsIgnoreCase("silent mating type cassette array")) {
					item = createItem("SilentMatingTypeCassetteArray");
				}else if (feature_type.equalsIgnoreCase("mating type region")) {
					item = createItem("MatingTypeRegion");
				}else if (feature_type.equalsIgnoreCase("intein encoding region")) {
					item = createItem("InteinEncodingRegion");
				} else if (feature_type.equalsIgnoreCase("recombination enhancer")) {
					item = createItem("RecombinationEnhancer");
				}else if (feature_type.equalsIgnoreCase("non transcribed region")) {
					item = createItem("NonTranscribedRegion");
				}

				// set for all types, so you can use LSF to query for these
				// different type of objects in a template.
				item.setAttribute("featureType", feature_type);
				item.setAttribute("primaryIdentifier", "SGD:"+primaryIdentifier);
				if (StringUtils.isNotEmpty(name)) item.setAttribute("name", name);
				item.setAttribute("secondaryIdentifier", secondaryIdentifier);
				item.setReference("organism", organism);
				if (StringUtils.isNotEmpty(symbol)) item.setAttribute("symbol", symbol);				
				if (StringUtils.isNotEmpty(description)) item.setAttribute("modDescription", description);
				//if (StringUtils.isNotEmpty(headline)) item.setAttribute("briefDescription", headline);
				if (qualifier != null) {
					if (StringUtils.isNotEmpty(qualifier)) {
						item.setAttribute("qualifier", qualifier);
					}
				}
				if (status != null) {
					if (StringUtils.isNotEmpty(status)) {
						item.setAttribute("status", status);
					}
				}
				String refId = item.getIdentifier();
				genes.put(featureNo, item);
				genesName.put(secondaryIdentifier, item);

				// ~~~ synonyms ~~~
				getSynonym(refId, "symbol", symbol);
				getSynonym(refId, "identifier", secondaryIdentifier);
			}
		}
		System.out.println("size of genes:  " + genes.size());
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processNISS(Connection connection) throws SQLException,
	ObjectStoreException {

		System.out.println("Processing NISS 55 features...");
		ResultSet res = PROCESSOR.getNISS(connection);

		while (res.next()) {

			String featureNo = res.getString("dbentity_id");
			if (genes.get(featureNo) == null) {
				//~~~ gene ~~~
				String primaryIdentifier = res.getString("sgdid");
				String secondaryIdentifier = res.getString("systematic_name");
				String symbol = res.getString("gene_name");
				String name = res.getString("name_description");
				String headline = res.getString("headline");
				String description = res.getString("description");
				String qualifier = res.getString("qualifier");
				String status = res.getString("dbentity_status");

				Item item = createItem("NotInSystematicSequenceOfS288C");

				// set for all types, so you can use LSF to query for these
				// different type of objects in a template.
				item.setAttribute("featureType", "not in systematic sequence of S288C");
				item.setAttribute("primaryIdentifier", "SGD:"+primaryIdentifier);
				if (StringUtils.isNotEmpty(name)) item.setAttribute("name", name);
				item.setAttribute("secondaryIdentifier", secondaryIdentifier);
				item.setReference("organism", organism);
				if (StringUtils.isNotEmpty(symbol)) item.setAttribute("symbol", symbol);				
				if (StringUtils.isNotEmpty(description)) item.setAttribute("description", description);
				//if (StringUtils.isNotEmpty(headline)) item.setAttribute("briefDescription", headline);
				if (qualifier != null) {
					if (StringUtils.isNotEmpty(qualifier)) {
						item.setAttribute("qualifier", qualifier);
					}
				}

				if (status != null) {
					if (StringUtils.isNotEmpty(status)) {
						item.setAttribute("status", status);
					}
				}
				String refId = item.getIdentifier();
				genes.put(featureNo, item);
				genesName.put(secondaryIdentifier, item);

				// ~~~ synonyms ~~~
				getSynonym(refId, "symbol", symbol);
				getSynonym(refId, "identifier", secondaryIdentifier);
			}
		}
	}

	/**
	 *
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processAlleles(Connection connection) throws SQLException,
			ObjectStoreException {

		System.out.println("Processing Alleles...");
		ResultSet res = PROCESSOR.getAlleleResults(connection);

		while (res.next()) {

			String alleleNo = res.getString("allele");
			String featureNo = res.getString("locus");

			System.out.println(alleleNo + "\t"+featureNo);

			Item gene = genes.get(featureNo);
			Item allele = alleles.get(alleleNo);

			if (gene != null) {

				String name = res.getString("allele_name");
				String description = res.getString("description");
				String aclass = res.getString("allele_class");
				String aliasName = res.getString("alias_name");
				String pmrefNo = res.getString("reference_id");
				String pmid = res.getString("pmid");
				String alleleSgdid = res.getString("allelesgdid");

				if (allele == null) {

					allele = createItem("Allele");
					allele.setAttribute("featureType", "Allele");
					if (StringUtils.isNotEmpty(name)) allele.setAttribute("name", name);
					if (StringUtils.isNotEmpty(description)) allele.setAttribute("description", description);
					if (StringUtils.isNotEmpty(aclass)) allele.setAttribute("alleleClass", aclass);
					if (StringUtils.isNotEmpty(alleleSgdid)) allele.setAttribute("alleleSgdid", alleleSgdid);
					if (StringUtils.isNotEmpty(aliasName)) allele.setAttribute("aliasName", aliasName);

					if(pmrefNo != null ) {
						Item publication = publications.get(pmrefNo);
						if (publication == null) {
							publication = createItem("Publication");
							if (StringUtils.isNotEmpty(pmid)) publication.setAttribute("pubMedId", pmid);
							publications.put(pmrefNo, publication);
						}
						allele.addToCollection("publications", publication);
					}
					alleles.put(alleleNo, allele);

				} else{

					if(pmrefNo != null ) {
						Item publication = publications.get(pmrefNo);
						if (publication == null) {
							publication = createItem("Publication");
							if (StringUtils.isNotEmpty(pmid)) publication.setAttribute("pubMedId", pmid);
							publications.put(pmrefNo, publication);
						}
						allele.addToCollection("publications", publication);
					}
				}
				//if (StringUtils.isNotEmpty(concatAliasName)) { allele.setAttribute("aliasName", concatAliasName); }
				allele.setReference("gene", gene.getIdentifier());
				gene.addToCollection("alleles", allele);

		  }//gene

		}//while
		System.out.println("size of alleles:  " + alleles.size());
	}

	/**
	 *
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processTranscripts(Connection connection) throws SQLException,
			ObjectStoreException {

		System.out.println("Processing Transcripts...");
		ResultSet res = PROCESSOR.getTranscriptResults(connection);

		while (res.next()) {

			String transcriptNo = res.getString("transcript");
			String name = res.getString("format_name");  //YAL009W_id006
			String [] n = name.split("_");
			String geneSecondaryId =  n[0];

			Item gene = genesName.get(geneSecondaryId);
			Item transcript = transcripts.get(transcriptNo);

			if (gene != null) {
				String start_index = res.getString("start_index");
				String end_index = res.getString("end_index");
				String strand = res.getString("strand");
				String newstrand = "";
				if (strand.equals("+")) {
					newstrand = "1";
				} else if (strand.equals("-")) {
					newstrand = "-1";
				} else if (strand.equals("0")) {
					newstrand = "0";
				}
				String in_gal = res.getString("in_gal");
				String in_ypd = res.getString("in_ypd");
				String in_ncbi = res.getString("in_ncbi");
				String ingal = "";
				if (in_gal.equals("t")) {
					ingal = "true";
				} else if (in_gal.equals("f")) {
					ingal = "false";
				}
				String inypd = "";
				if (in_ypd.equals("t")) {
					inypd = "true";
				} else if (in_ypd.equals("f")) {
					inypd = "false";
				}
				String incbi = "";
				if (in_ncbi.equals("t")) {
					incbi = "true";
				} else if (in_ncbi.equals("f")) {
					incbi = "false";
				}
				String pmrefNo = res.getString("reference_id");
				String pmid = res.getString("pmid");
				if (transcript == null) {
					transcript = createItem("MRNA");
					transcript.setAttribute("featureType", "mRNA");
					if (StringUtils.isNotEmpty(name)) transcript.setAttribute("primaryIdentifier", "SGD:"+name);
					if (StringUtils.isNotEmpty(in_gal)) transcript.setAttribute("in_gal", ingal);
					if (StringUtils.isNotEmpty(in_ypd)) transcript.setAttribute("in_ypd", inypd);
					if (StringUtils.isNotEmpty(in_ncbi)) transcript.setAttribute("in_ncbi", incbi);

					if(pmrefNo != null ) {
						Item publication = publications.get(pmrefNo);
						if (publication == null) {
							publication = createItem("Publication");
							publication.setAttribute("pubMedId", pmid);
							publications.put(pmrefNo, publication);
						}
						transcript.addToCollection("publications", publication);
					}
					transcript.setReference("gene", gene.getIdentifier());

					//location
					String refId = genesChromosomes.get(geneSecondaryId);
					transcript.setReference("chromosome", refId);
					String locationRefId = getLocation(transcript, refId, start_index, end_index, newstrand);
					transcript.setReference("chromosomeLocation", locationRefId);

					transcripts.put(transcriptNo, transcript);
				}
			} //gene
		}//while
		System.out.println("size of transcripts:  " + transcripts.size());
	}


	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processParalogs(Connection connection) throws SQLException,
	ObjectStoreException {

		System.out.println("Processing Paralog pairs...");
		ResultSet res = PROCESSOR.getParalogs(connection); // ordered by featureNo

		while (res.next()) {

			String parentFeatureNo = res.getString("parent_id");
			String childFeatureNo = res.getString("child_id");
			String refNo = res.getString("reference_id");		
			String source="";
			if(refNo.equalsIgnoreCase("526423")){ //hack for PMID 16169922
				source= "YGOB";
			}else{
				source = "SGD";
			}
			Item pmid = getExistingPub(refNo);
			Item parentGene = genes.get(parentFeatureNo);
			Item childGene = genes.get(childFeatureNo);

			if(parentGene!= null && childGene != null) {	
				processHomologues(parentGene.getIdentifier(), childGene.getIdentifier(), pmid.getIdentifier(), source);
				processHomologues(childGene.getIdentifier(), parentGene.getIdentifier(), pmid.getIdentifier(), source);

			}

		}
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processRegulation(Connection connection) throws SQLException,
	ObjectStoreException {

		System.out.println("Processing Regulation data...");
		ResultSet res = PROCESSOR.getRegulationData(connection); // ordered by featureNo

		while (res.next()) {

			String factorGene = res.getString("regulator_id");
			String targetGene = res.getString("target_id");
			String evidenceCode = res.getString("format_name");	
			String condition = res.getString("happens_during");	
			String regulationDirection = res.getString("direction");
			String regulationType = res.getString("regulation_type");	
			String regulatorType = res.getString("regulator_type");	
			String annotationType = res.getString("annotation_type");	
			String pmid = res.getString("pmid");	
			String refNo = res.getString("refNo");
			String strainBackground = res.getString("strain_background");	
			String source = res.getString("source");

			getRegulation(factorGene, targetGene, evidenceCode, condition,  regulationDirection,  pmid,
					refNo, source, strainBackground,regulationType,regulatorType, annotationType);

		}
	}

	/**
	 * 
	 * @param factorGene
	 * @param targetGene
	 * @param evidenceCode
	 * @param condition
	 * @param regulationDirection
	 * @param pmid
	 * @param source
	 * @param pvalue
	 * @param fdr
	 * @param strainBackground
	 * @param strain
	 * @throws ObjectStoreException
	 */
	private void getRegulation(String factorGene, String targetGene, 
			String evidenceCode, String condition, String regulationDirection, String pmid, String refNo, String source, 
			String strainBackground, String regulationType, String regulatorType, String annotationType) 
					throws ObjectStoreException {

		Item rGene = genes.get(factorGene);
		Item tGene = genes.get(targetGene);

		if (rGene != null && tGene != null) {       	               	           			

			Item bindingSite = createItem("TFBindingSite");  
			String name = factorGene + "_binding_site";

			bindingSite.setAttribute("name", name);
			bindingSite.setReference("regulator", rGene.getIdentifier());
			bindingSite.setReference("gene", tGene.getIdentifier());

			Item evidence = createItem("RegulationEvidence");

			if (StringUtils.isNotEmpty(evidenceCode)) {

				Item eco = ecoMap.get(evidenceCode);
				if(eco == null) {
					eco = createItem("ECOTerm");
					ecoMap.put(evidenceCode, eco);
					eco.setAttribute("identifier", evidenceCode); 
					try {
						store(eco);
					} catch (ObjectStoreException e) {
						throw new ObjectStoreException(e);
					}	
				}
				evidence.setReference("ontologyTerm", eco.getIdentifier());
			}


			try {
				store(evidence);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}			

			bindingSite.setReference("regEvidence", evidence.getIdentifier());                  	

			if (StringUtils.isNotEmpty(condition)) {
				bindingSite.setAttribute("experimentCondition", condition);
			}
			if (StringUtils.isNotEmpty(regulationDirection)) {
				bindingSite.setAttribute("regulationDirection", regulationDirection);
			}
			if (StringUtils.isNotEmpty(source)) { 
				bindingSite.setAttribute("datasource", source);
			}
			if (StringUtils.isNotEmpty(strainBackground)) {
				String q[] = strainBackground.split("_");
				int len = q.length;
				String strain;
				if(len == 3){
					 strain = q[len-1];
				}else{
					 strain = "Other";
				}
				bindingSite.setAttribute("strainBackground", strain);
			}
			if (StringUtils.isNotEmpty(regulationType)) {
				bindingSite.setAttribute("regulationType", regulationType);
			}
			if (StringUtils.isNotEmpty(regulatorType)) {
				bindingSite.setAttribute("regulatorType", regulatorType);
			}
			if (StringUtils.isNotEmpty(annotationType)) {
				bindingSite.setAttribute("annotationType", annotationType);
			}

			Item publication = publications.get(refNo);

			if(publication == null) {
				publication = createItem("Publication");
				publication.setAttribute("pubMedId", pmid); 
				publications.put(refNo, publication);	
			}

			bindingSite.addToCollection("publications", publication);      

			try {
				store(bindingSite);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}	

		}

	}

	/**
	 * 
	 * @param gene1
	 * @param gene2
	 * @param pmid
	 * @param source
	 * @throws ObjectStoreException
	 */
	private void processHomologues(String gene1, String gene2, String pmid, String source) throws ObjectStoreException {
		if (gene1 == null || gene2 == null) {
			return;
		}
		Item homologue = createItem("Homologue");
		homologue.setReference("gene", gene1);
		homologue.setReference("homologue", gene2);
		homologue.setAttribute("type", "paralogue"); 
		homologue.setAttribute("source", source);
		homologue.setReference("publication", pmid);
		store(homologue);
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processAliases(Connection connection) throws SQLException, ObjectStoreException {

		ResultSet res = PROCESSOR.getAliases(connection); // ordered by featureNo

		System.out.println("Processing ALiases...");
		while (res.next()) {

			String geneFeatureNo = res.getString("dbentity_id");
			String alias_type = res.getString("alias_type");
			String alias_name = res.getString("display_name");

			Item item = genes.get(geneFeatureNo);
			if (item != null) {
				getSynonym(item.getIdentifier(), alias_type, alias_name); // adding sgd_aliases as synonyms
				String name = genesAliases.get(geneFeatureNo);
				if (name == null) {
					genesAliases.put(geneFeatureNo, alias_name);
				} else {
					String newname = name + "|" + alias_name;
					genesAliases.put(geneFeatureNo, newname);
				}
			}

		}

		// Add the concatenated string as alias name
		Set<Map.Entry<String, String>> set = genesAliases.entrySet();
		java.util.Iterator<Map.Entry<String, String>> it = set.iterator();

		while (it.hasNext()) {
			Map.Entry<String, String> anEntry = it.next();
			String geneFeatureNo = anEntry.getKey();
			String alias = anEntry.getValue();
			Item item = genes.get(geneFeatureNo);
			if(StringUtils.isNotEmpty(alias)){ item.setAttribute("sgdAlias", alias); }
		}

	}


	/**
	 *
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processAlleleAliases(Connection connection) throws SQLException, ObjectStoreException {

		ResultSet res = PROCESSOR.getAlleleAliasesResults(connection);

		System.out.println("Processing Alleles ALiases...");
		while (res.next()) {
			String alleleFeatureNo = res.getString("allele");
			String allele_name = res.getString("allele_name");
			String alias_name = res.getString("alias_name");
			String s1 = alias_name.replace("\"", "");
			String s2 = s1.replace("{", "");
			String s3 = s2.replace("}", "");


			Item item = alleles.get(alleleFeatureNo);
			if (item != null) {
				if (StringUtils.isNotEmpty(s3)) {
					item.setAttribute("aliasName", s3);
				}
			}
		}

	}


	/**
	 *
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processAllPathways(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing ALL Pathways...");
		ResultSet res = PROCESSOR.getAllPathways(connection); // ordered by featureNo

		while (res.next()) {

			String dbxref_id = res.getString("biocyc_id"); //pathway name
			String dbxref_name = res.getString("display_name"); //pathway identifier i.e. short name
			String summary_type = res.getString("summary_type");
			String text = res.getString("text");
			String refNo = res.getString("reference_id"); // refs for pathways

			//Item pathway = pathways.get(dbxref_id);
			//if (pathway != null) {
			getPathway(dbxref_id, dbxref_name, summary_type, text, refNo);
			//}
		}
	}


	private String getPathway(String id, String name, String summaryType, String text, String refNo)
			throws ObjectStoreException {

		Item crf = pathways.get(id);
		Item ps = pathwaysummarys.get(id);

		if (crf == null) {

			crf = createItem("Pathway");
			crf.setAttribute("identifier", id);
			crf.setAttribute("name", name);

			if(ps == null) {
				ps = createItem("PathwaySummary");
				ps.setAttribute("summaryType", summaryType);
				if (!StringUtils.isEmpty(text)) { ps.setAttribute("text", text);}
			}else{
				ps.setAttribute("summaryType", summaryType);
				if (!StringUtils.isEmpty(text)) { ps.setAttribute("text", text);}
			}

			Item publication = publications.get(refNo);
			if(publication == null) {
				publication = createItem("Publication");
				publications.put(refNo, publication);
			}
			ps.addToCollection("publications", publication);
			pathwaysummarys.put(id,ps);
			pathways.put(id, crf);
		}else{

			if(ps == null) {
				ps = createItem("PathwaySummary");
				ps.setAttribute("summaryType", summaryType);
				if (!StringUtils.isEmpty(text)) { ps.setAttribute("text", text);}
			}else{
				ps.setAttribute("summaryType", summaryType);
				if (!StringUtils.isEmpty(text)) { ps.setAttribute("text", text);}
			}

			Item publication = publications.get(refNo);
			if(publication == null) {
				publication = createItem("Publication");
				publications.put(refNo, publication);
			}
			ps.addToCollection("publications", publication);
			pathwaysummarys.put(id,ps);
			pathways.put(id, crf);
		}

		crf.setReference("summary", ps);
		String refId = crf.getIdentifier();
		return refId;

	}


	/**
	 *
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processGenePathways(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing Pathways...");
		ResultSet res = PROCESSOR.getGenePathways(connection); // ordered by featureNo

		while (res.next()) {
			String geneFeatureNo = res.getString("dbentity_id");
			String dbxref_id = res.getString("biocyc_id"); //pathway name
			String dbxref_name = res.getString("display_name"); //pathway identifier i.e. short name
			//String summary_type = res.getString("summary_type");
			//String text = res.getString("text");
			String refNo = res.getString("reference_id"); // refs for pathways

			Item item = genes.get(geneFeatureNo);
			if (item != null) {
				getGenePathway(item.getIdentifier(), dbxref_id, dbxref_name, refNo); //summary_type, text,
			}
		}
	}

	private String getGenePathway(String geneIdentifier, String id, String name, String refNo)
			throws ObjectStoreException {

		Item crf = pathways.get(id);
		if (crf == null) {

			crf = createItem("Pathway");
			crf.setAttribute("identifier", id);
			crf.setAttribute("name", name);
			//crf.setAttribute("summaryType", summaryType);
			//if (!StringUtils.isEmpty(text)) { crf.setAttribute("text", text);}

			Item publication = publications.get(refNo);
			if(publication == null) {
				publication = createItem("Publication");
				publications.put(refNo, publication);
			}
			crf.addToCollection("publications", publication);

			pathways.put(id, crf);

			crf.addToCollection("genes", geneIdentifier);

		}else{
			Item publication = publications.get(refNo);
			if(publication == null) {
				publication = createItem("Publication");
				publications.put(refNo, publication);
			}
			crf.addToCollection("publications", publication);
			crf.addToCollection("genes", geneIdentifier);
		}

		String refId = crf.getIdentifier();
		return refId;

	}


	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processCrossReferences(Connection connection) throws SQLException, ObjectStoreException {

		System.out.println("Processing DbXRefs...");
		ResultSet res = PROCESSOR.getCrossReferences(connection); // ordered by featureNo

		while (res.next()) {
			String geneFeatureNo = res.getString(1);
			String dbxref_id = res.getString(2);
			String dbx_source = res.getString(3);
			String dbxref_type = res.getString(4);

			Item item = genes.get(geneFeatureNo);

			if (item != null) {
				String url = getCrossRefURL(dbx_source,dbxref_id, dbxref_type);
				getCrossReference(item.getIdentifier(), dbxref_id, dbx_source, dbxref_type, url);
			}
		}
	}

	/**
	 * 
	 * @param source
	 * @param dbxref_id
	 * @param dbxref_type
	 * @return
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private String getCrossRefURL(String source, String dbxref_id, String dbxref_type) throws SQLException, ObjectStoreException {

		String url = "";

		if(source.equalsIgnoreCase("AspGD")){
			url = "http://www.aspgd.org/cgi-bin/locus.pl?locus="+dbxref_id;
		}else if (source.equalsIgnoreCase("BioGRID")){
			url = "http://thebiogrid.org/"+dbxref_id;
		}else if (source.equalsIgnoreCase("CGD")){
			url = " http://www.candidagenome.org/cgi-bin/locus.pl?locus="+dbxref_id;
		}else if (source.equalsIgnoreCase("DIP")){
			url = "http://dip.doe-mbi.ucla.edu/dip/Browse.cgi?PK="+dbxref_id;
		}else if (source.equalsIgnoreCase("EBI")){
			url = "http://www.uniprot.org/uniprot/"+dbxref_id;
		}else if (source.equalsIgnoreCase("EUROSCARF")){
			url = "http://web.uni-frankfurt.de/fb15/mikro/euroscarf/data/"+dbxref_id;
		}else if (source.equalsIgnoreCase("IUBMB")){
			url = "http://www.expasy.org/enzyme/"+dbxref_id;
		}else if (source.equalsIgnoreCase("LoQate")){
			url = "http://www.weizmann.ac.il/molgen/loqate/gene/view/"+dbxref_id;
		}else if (source.equalsIgnoreCase("MetaCyc")){
			url = "http://pathway.yeastgenome.org/YEAST/new-image?type=PATHWAY&object="+dbxref_id;
		}else if (source.equalsIgnoreCase("PDB")){
			url = "http://www.rcsb.org/pdb/explore/explore.do?structureId="+dbxref_id;
		}else if (source.equalsIgnoreCase("NCBI") && (dbxref_type.equalsIgnoreCase("Gene ID"))){
			url = "http://www.ncbi.nlm.nih.gov/gene/"+dbxref_id;
		}else if (source.equalsIgnoreCase("NCBI") && (dbxref_type.equalsIgnoreCase("RefSeq protein version ID"))){
			url = "http://www.ncbi.nlm.nih.gov:80/entrez/viewer.fcgi?val="+dbxref_id;
		}else if (source.equalsIgnoreCase("PomBase")){
			url = "http://www.pombase.org/spombe/result/"+dbxref_id;
		}else if (source.equalsIgnoreCase("RNAcentral")){
			url = "http://rnacentral.org/rna/"+dbxref_id;
		}else if (source.equalsIgnoreCase("TCDB")){
			url = "http://www.tcdb.org/tcdb/index.php?tc="+dbxref_id;
		}

		return url;
	}

	/**
	 * 
	 * @param subjectId
	 * @param id
	 * @param source
	 * @param dbxref_type
	 * @param url
	 * @return
	 * @throws ObjectStoreException
	 */
	private String getCrossReference(String subjectId, String id, String source, String dbxref_type, String url)
			throws ObjectStoreException {

		String refId = "";
		Item crf = createItem("CrossReference");
		crf.setReference("subject", subjectId);
		crf.setAttribute("identifier", id);
		crf.setAttribute("dbxreftype", dbxref_type);
		if(! StringUtils.isEmpty(url)) {
			crf.setAttribute("url", url);	
		}

		String dsId = datasources.get(source);
		if (dsId == null) {
			Item ds = createItem("DataSource");
			//ds.setAttribute("name", source);
			String join = source+"/"+dbxref_type;
			ds.setAttribute("name", join);
			String sourceUrl = getSourceURL(source);
			if(! StringUtils.isEmpty(sourceUrl)) {
				ds.setAttribute("url", sourceUrl);
			}
			try {
				store(ds);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}

			crf.setReference("source", ds.getIdentifier());
			datasources.put(join, ds.getIdentifier());
		} else {
			crf.setReference("source", dsId);
		}

		try {
			store(crf);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}

		refId = crf.getIdentifier();
		return refId;

	}
	/**
	 * 
	 * @param source
	 * @return
	 * @throws ObjectStoreException
	 */
	private String getSourceURL(String source) throws ObjectStoreException {

		String url = "";

		if(source.equalsIgnoreCase("AspGD")){
			url = "http://www.aspgd.org/";
		}else if (source.equalsIgnoreCase("BioGRID")){
			url = "http://thebiogrid.org/";
		}else if (source.equalsIgnoreCase("CGD")){
			url = " http://www.candidagenome.org/";
		}else if (source.equalsIgnoreCase("DIP")){
			url = "http://dip.doe-mbi.ucla.edu/";
		}else if (source.equalsIgnoreCase("EBI")){
			url = "http://www.uniprot.org/uniprot/";
		}else if (source.equalsIgnoreCase("EUROSCARF")){
			url = "http://web.uni-frankfurt.de/";
		}else if (source.equalsIgnoreCase("IUBMB")){
			url = "http://www.expasy.org/";
		}else if (source.equalsIgnoreCase("LoQate")){
			url = "http://www.weizmann.ac.il/molgen/loqate/";
		}else if (source.equalsIgnoreCase("MetaCyc")){
			url = "http://pathway.yeastgenome.org/";
		}else if (source.equalsIgnoreCase("PDB")){
			url = "http://www.rcsb.org/pdb/";
		}else if (source.equalsIgnoreCase("NCBI")){
			url = "http://www.ncbi.nlm.nih.gov/";
		}else if (source.equalsIgnoreCase("PomBase")){
			url = "http://www.pombase.org/";
		}else if (source.equalsIgnoreCase("RNAcentral")){
			url = "http://rnacentral.org/rna/";
		}else if (source.equalsIgnoreCase("TCDB")){
			url = "http://www.tcdb.org/tcdb/";
		}

		return url;
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 * @throws Exception
	 */

	private void processGeneLocations(Connection connection)
			throws SQLException, ObjectStoreException, Exception {

		System.out.println("Processing GeneLocations...");
		ResultSet res = PROCESSOR.getChromosomalFeatureLocationResults(connection);

		while (res.next()) {
			String featureNo = res.getString("contig_id");
			String featureName = res.getString("format_name");
			String featureType = res.getString("feature_type");
			String geneFeatureNo = res.getString("dbentity_id");
			String geneFeatureName = res.getString("gene_name");
			String secondaryIdentifier = res.getString("systematic_name");
			String strand = res.getString("strand");
			String seq_length = res.getString(11);

			String newstrand = "";
			if (strand.equals("+")) {
				newstrand = "1";
			} else if (strand.equals("-")) {
				newstrand = "-1";
			} else if (strand.equals("0")) {
				newstrand = "0";
			}
			String fixed_chromosome_no = getFixedChrName(featureName);
			Item item = genes.get(geneFeatureNo);

			// ~~~ chromosome OR plasmid ~~~
			String refId = null;
			if (featureType.equalsIgnoreCase("plasmid")) {
				refId = getPlasmid(fixed_chromosome_no); 
				item.setReference("plasmid", refId);
			} else if (featureType.equalsIgnoreCase("chromosome")) {
				refId = getChromosome(fixed_chromosome_no);
				item.setReference("chromosome", refId);
			}

			// ~~~ location ~~~
			String locationRefId = getLocation(item, refId, res.getString("start_index"), res.getString("end_index"),newstrand); 

			if (featureType.equalsIgnoreCase("plasmid")) {
				item.setReference("plasmidLocation", locationRefId);
			} else if (featureType.equalsIgnoreCase("chromosome")) {
				item.setReference("chromosomeLocation", locationRefId);
			}
			// ~~ add sequence
			String seqRefId = getSequence(geneFeatureNo, res.getString("residues"), seq_length);
			item.setReference("sequence", seqRefId);
			item.setAttribute("length", seq_length);

			genesChromosomes.put(secondaryIdentifier, refId);

		}
		res.close();
		System.out.println("size of genes:  " + genes.size());
	}

	/**
	 * To add a chromosome number to the NISS and NPM features that do have a
	 * chromosome number as parent
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 *
	private void processChrLocations(Connection connection)
			throws SQLException, ObjectStoreException, Exception {

		ResultSet res = PROCESSOR.getChromosomeLocationResults(connection);
		System.out.println("Processing ChrLocations...");
		while (res.next()) {
			String featureNo = res.getString("feature_no");
			String geneFeatureNo = res.getString("gene_feature_no");
			String featureType = res.getString("feature_type");
			String featureName = res.getString("identifier");

			String fixed_chromosome_no = getFixedChrName(featureName);

			Item item = genes.get(geneFeatureNo);

			// ~~~ chromosome OR plasmid ~~~
			String refId = null;
			if (featureType.equalsIgnoreCase("plasmid")) {
				refId = getPlasmid(fixed_chromosome_no);
				item.setReference("plasmid", refId);
			} else if (featureType.equalsIgnoreCase("chromosome")) {
				refId = getChromosome(fixed_chromosome_no);
				item.setReference("chromosome", refId);
			}

			// ~~~ location ~~~
			String locationRefId = getLocation(item, refId, "1", "1", "0"); // was

			if (featureType.equalsIgnoreCase("plasmid")) {
				item.setReference("plasmidLocation", locationRefId);
			} else if (featureType.equalsIgnoreCase("chromosome")) {
				item.setReference("chromosomeLocation", locationRefId);
			}

		}
		res.close();
	}*/

	private void processGeneChildrenLocations(Connection connection) throws SQLException, ObjectStoreException, Exception {

		System.out.println("Processing GeneChildrenLocations...");
		ResultSet res = PROCESSOR.getChildrenFeatureLocationResults(connection);

		while (res.next()) {

			String geneFeatureNo = res.getString("parent_id");
			String parentFeatureType = res.getString("parent_type").trim();

			String geneChildFeatureNo = res.getString("child_id");
			String childFeatureType = res.getString("child_type").trim();

			String chromosome_no = res.getString("format_name"); //root chr.number
			//String secondaryIdentifier = res.getString("child_identifier"); //child identifier is wrong -- fix it 11/13
			String primaryIdentifier = res.getString("child_sgdid")+"_C"; // SXX

			String maxcoord = res.getString("child_end_coord");
			String mincoord = res.getString("child_start_coord");
			String strand = res.getString("strand");

			String seq = res.getString("residues");
			String seqLen = res.getString("seq_length");
			String child_status = res.getString("child_status");

			String newstrand = "";
			if (strand.equals("+")) {
				newstrand = "1";
			} else if (strand.equals("-")) {
				newstrand = "-1";
			} else if (strand.equals("0")) {
				newstrand = "0";
			}

			String fixed_chromosome_no = getFixedChrName(chromosome_no);

			// figure out why duplicates in the SQL..???..
			/*if (featureMap.get(geneChildFeatureNo) == null) {
				featureMap.put(geneChildFeatureNo, geneFeatureNo);
			} else {
				continue;
			}*/

			Item parent = genes.get(geneFeatureNo);
			// create the child Item
			Item childItem = getChildItem(childFeatureType);

			childItem.setAttribute("primaryIdentifier", primaryIdentifier);
			childItem.setReference("organism", organism);
			childItem.setAttribute("featureType", childFeatureType);
			childItem.setAttribute("status", child_status);

			if (childFeatureType.equalsIgnoreCase("intron")) {
				childItem.addToCollection("genes", parent.getIdentifier());
			} else {
				String refname = getReferenceName(childFeatureType,parentFeatureType);
				childItem.setReference(refname, parent.getIdentifier());
			}
			// ~~ add sequence
			String seqRefId = getSequence(geneChildFeatureNo, seq, seqLen);
			childItem.setReference("sequence", seqRefId);

			// ~~~ chromosome and location ~~~
			String refId = null;
			if (chromosome_no.equalsIgnoreCase("2-micron_plasmid")) {
				refId = getPlasmid(fixed_chromosome_no);
				childItem.setReference("plasmid", refId);
				String locationRefId = getLocation(childItem, refId, mincoord, maxcoord, newstrand);
				childItem.setReference("plasmidLocation", locationRefId);
			} else {
				refId = getChromosome(fixed_chromosome_no);
				childItem.setReference("chromosome", refId);
				String locationRefId = getLocation(childItem, refId, mincoord, maxcoord, newstrand);
				childItem.setReference("chromosomeLocation", locationRefId);
			}

			try {
				store(childItem);
			} catch (ObjectStoreException e) {
				e.printStackTrace();
				throw new ObjectStoreException(e);
			}

			// ~~~ store these last ~~~
			String childId = childItem.getIdentifier();
			getSynonym(childId, "identifier", primaryIdentifier);

		}
		res.close();

	}

	private String getReferenceName(String type, String ptype)
			throws ObjectStoreException {

		String name = "";

		if (type.equalsIgnoreCase("ARS_consensus_sequence") && ptype.equalsIgnoreCase("ARS")) {
			name = "ars";
		}else if (type.equalsIgnoreCase("CDS") && ptype.equalsIgnoreCase("blocked reading frame")) {
			name = "blockedreadingframe";
		}else if (type.equalsIgnoreCase("CDS") && ptype.equalsIgnoreCase("ORF")) {
			name = "orf";
		} else if (type.equalsIgnoreCase("uORF") && ptype.equalsIgnoreCase("ORF")) {
			name = "orf";
		}else if (type.equalsIgnoreCase("CDS") && ptype.equalsIgnoreCase("pseudogene")) {
			name = "pseudogene";
		} else if (type.equalsIgnoreCase("CDS") && ptype.equalsIgnoreCase("transposable element gene")) {
			name = "transposableelementgene";
		}else if (type.equalsIgnoreCase("centromere_DNA_Element_I")) {
			name = "centromere";
		} else if (type.equalsIgnoreCase("centromere_DNA_Element_II")) {
			name = "centromere";
		} else if (type.equalsIgnoreCase("centromere_DNA_Element_III")) {
			name = "centromere";
		} else if (type.equalsIgnoreCase("external_transcribed_spacer_region")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("five_prime_UTR_intron")) {
			name = "orf";
		}else if (type.equalsIgnoreCase("intein_encoding_region")) {
			name = "orf";
		}else if (type.equalsIgnoreCase("internal_transcribed_spacer_region")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("intron") && ptype.equalsIgnoreCase("ORF")) {
			name = "orf";
		}else if (type.equalsIgnoreCase("intron") && ptype.equalsIgnoreCase("rRNA gene")) {
			name = "rrna_gene";
		}else if (type.equalsIgnoreCase("intron") && ptype.equalsIgnoreCase("snoRNA gene")) {
			name = "snorna_gene";
		}else if (type.equalsIgnoreCase("intron") && ptype.equalsIgnoreCase("tRNA gene")) {
			name = "trna_gene";
		} else if (type.equalsIgnoreCase("non_transcribed_region")) {
			name = "ncrna_gene";
		}else if (type.equalsIgnoreCase("noncoding_exon") && ptype.equalsIgnoreCase("ncRNA gene")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("noncoding_exon") && ptype.equalsIgnoreCase("rRNA gene")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("noncoding_exon") && ptype.equalsIgnoreCase("snoRNA gene")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("noncoding_exon") && ptype.equalsIgnoreCase("snRNA gene")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("noncoding_exon") && ptype.equalsIgnoreCase("tRNA gene")) {
			name = "ncrna_gene";
		} else if (type.equalsIgnoreCase("plus_1_translational_frameshift") && ptype.equalsIgnoreCase("ORF")) {
			name = "orf";
		} else if (type.equalsIgnoreCase("plus_1_translational_frameshift") && ptype.equalsIgnoreCase("pseudogene")) {
			name = "pseudogene";
		} else if (type.equalsIgnoreCase("plus_1_translational_frameshift") && ptype.equalsIgnoreCase("transposable element gene")) {
			name = "transposableelementgene";
		} else if (type.equalsIgnoreCase("telomeric_repeat") && ptype.equalsIgnoreCase("telomere")) {
			name = "telomere";
		}  else if (type.equalsIgnoreCase("telomeric_repeat") && ptype.equalsIgnoreCase("ORF")) {
			name = "orf";
		} else if (type.equalsIgnoreCase("X_element") && ptype.equalsIgnoreCase("telomere")) {
			name = "telomere";
		} else if (type.equalsIgnoreCase("X_element_combinatorial_repeat") && ptype.equalsIgnoreCase("telomere")) {
			name = "telomere";
		}   else if (type.equalsIgnoreCase("Y_prime_element") && ptype.equalsIgnoreCase("telomere")) {
			name = "telomere";
		}  else if (type.equalsIgnoreCase("noncoding_exon") && ptype.equalsIgnoreCase("telomerase RNA gene")) {
			name = "ncrna_gene";
		}else if (type.equalsIgnoreCase("W_region") && ptype.equalsIgnoreCase("silent mating type cassette array")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("W_region") && ptype.equalsIgnoreCase("mating type region")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("Y_region") && ptype.equalsIgnoreCase("mating type region")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("X_region") && ptype.equalsIgnoreCase("silent mating type cassette array")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("Y_region") && ptype.equalsIgnoreCase("silent mating type cassette array")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("X_region") && ptype.equalsIgnoreCase("mating type region")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("Z1_region") && ptype.equalsIgnoreCase("silent mating type cassette array")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("Z1_region") && ptype.equalsIgnoreCase("mating type region")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("Z2_region") && ptype.equalsIgnoreCase("silent mating type cassette array")) {
			name = "matingtyperegion";
		}else if (type.equalsIgnoreCase("Z2_region") && ptype.equalsIgnoreCase("mating type region")) {
			name = "matingtyperegion";
		}

		return name;

	}

	private String getCasedName(String name) throws Exception {

		String newname = name;

		StringBuffer sb = new StringBuffer();
		Matcher m = Pattern
				.compile("([a-z])([a-z]*)", Pattern.CASE_INSENSITIVE).matcher(
						name);
		while (m.find()) {
			m.appendReplacement(sb, m.group(1).toUpperCase()
					+ m.group(2).toLowerCase());
		}
		newname = m.appendTail(sb) + "p".toString();

		return newname;

	}

	private void processProteins(Connection connection) throws SQLException,
	ObjectStoreException, Exception {

		System.out.println("Processing Proteins...");
		ResultSet res = PROCESSOR.getProteinResults(connection);
		
		while (res.next()) {
			String featureNo = res.getString("dbentity_id");
			String primaryIdentifier = res.getString("sgdid");
			String secondaryIdentifier = res.getString("format_name");
			String symbol = res.getString("display_name");
			String residues = res.getString("residues");
			String length = res.getString(6);

			Item item = genes.get(featureNo);

			// ~~~ sequence ~~~
			Item protein = createItem("Protein");
			protein.setAttribute("primaryIdentifier", primaryIdentifier);
			protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
			protein.setAttribute("length", length);
			protein.setReference("organism", organism);

			if (symbol != null) {
				String modSymbol = getCasedName(symbol);
				protein.setAttribute("symbol", modSymbol);
			}

			Item seq = createItem("Sequence");
			seq.setAttribute("residues", residues);
			seq.setAttribute("length", length);

			try {
				store(seq);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}

			protein.setReference("sequence", seq.getIdentifier());
			if(item != null) protein.addToCollection("genes", item.getIdentifier());
			proteins.put(featureNo, protein);

		}

	}

	private void processProteinHalfLife(Connection connection) throws SQLException,
	ObjectStoreException, Exception {

		System.out.println("Processing Proteins Half-life data...");
		ResultSet res = PROCESSOR.getProteinHalfLifeResults(connection);

		while (res.next()) {
			String featureNo = res.getString("dbentity_id");
			String units = res.getString("data_unit");
			String value = res.getString("data_value");
			if(value.equalsIgnoreCase(">= 100")) value = "10000";
			String refNo = res.getString("referencedbentity");
			String pmid = res.getString("pmid");
			String experiment = "protein half-life";	

			Item protein = proteins.get(featureNo);
			Item pmods = getProteinHalfLife(experiment, value, units, pmid, refNo);
			protein.addToCollection("proteinHalfLife", pmods.getIdentifier());

		}

	}


	private Item getProteinHalfLife(String experiment, String value, String units, String pmid, String refNo) throws ObjectStoreException {

		Item item = createItem("ProteinHalfLife");

		if(StringUtils.isNotEmpty(experiment)){  item.setAttribute("experiment", experiment);}
		if(StringUtils.isNotEmpty(value)){  item.setAttribute("value", value);}
		if(StringUtils.isNotEmpty(units)){  item.setAttribute("units", units);}
		item.setAttribute("source", "SGD");

		Item publication = publications.get(refNo);

		if (publication == null) {	
			publication = createItem("Publication");			
			publication.setAttribute("pubMedId", pmid);			 
			publications.put(refNo, publication);
		}
		item.setReference("publication", publication);   

		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}	
		return item;
	}


	private void processProteinDomains(Connection connection) throws SQLException, ObjectStoreException, Exception {

		System.out.println("Processing Proteins Domains data...");
		ResultSet res = PROCESSOR.getProteinDomainsResults(connection);

		while (res.next()) {
			String featureNo = res.getString("dbentity_id");
			String domainDesc = res.getString("description");
			String domainMatch = res.getString("display_name");
			String sstart = res.getString("start_index");
			String send = res.getString("end_index");
			String runDate = res.getString("date_of_run");
			String method = res.getString("url_type");
			String interpro = res.getString("interpro_id");

			Item protein = proteins.get(featureNo);
			Item pdomain = getDomain(domainMatch, domainDesc, sstart, send, runDate, method);
			protein.addToCollection("proteinDomains", pdomain.getIdentifier());

		}

	}

	private Item getDomain(String identifier, String description, String s_start, String s_end, String runDate, String method)
			throws ObjectStoreException{

		Item item = createItem("ProteinDomain");
		item.setAttribute("name", identifier);
		if(StringUtils.isNotEmpty(description)) item.setAttribute("description", description);
		item.setAttribute("start", s_start);
		item.setAttribute("end", s_end );
		item.setAttribute("runDate", runDate);            	            
		item.setAttribute("method", method);

		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}	
		return item;
	}

	private void processProteinModifications(Connection connection) throws SQLException, ObjectStoreException, Exception {
		
		System.out.println("Processing Proteins Modification data...");
		ResultSet res = PROCESSOR.getProteinModificationResults(connection);
		while (res.next()) {
			String featureNo = res.getString("dbentity_id");
			String siteindex = res.getString("site_index");
			String siteresidue = res.getString("site_residue");
			String modType = res.getString("display_name");
			String modifier = res.getString("modifier_name");
			String modSite = siteresidue+siteindex;
			String source =  "SGD";
			String pmid = res.getString("pmid");
			String refNo = res.getString("referencedbentity");

			Item protein = proteins.get(featureNo);						
			Item pmods = getProteinMod(modSite, modType, modifier, source, refNo, pmid);
			protein.addToCollection("proteinModificationSites", pmods.getIdentifier());				
		}

	}

	/**
	 * 
	 * @param modSite
	 * @param modType
	 * @param source
	 * @param pmid
	 * @return
	 */
	private Item getProteinMod(String modSite, String modType, String modifier, String source, String refNo, String pmid) throws ObjectStoreException {

		Item item = createItem("ProteinModificationSite");

		item.setAttribute("modificationType", modType);
		item.setAttribute("modificationSite", modSite);
		item.setAttribute("source", source);

		if (!StringUtils.isEmpty(modifier)) {
			//Item gene = genes.get(modifier);
			//item.addToCollection("modifier", gene.getIdentifier());
			item.setAttribute("modifier", modifier);
		}

		Item publication = publications.get(refNo);

		if(publication == null) {
			publication = createItem("Publication");
			publication.setAttribute("pubMedId", pmid);
			publications.put(refNo, publication);		
		}
		item.setReference("publication", publication);      

		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}	

		return item;
	}


	/**
	 * 
	 * @param modSite
	 * @param modType
	 * @param source
	 * @param pmid
	 * @return
	 */
	private Item getProteinTermini(String modSite, String modType, String pmid, String refNo, String category, String seq) throws ObjectStoreException {

		Item item = createItem("ProteinModificationSite");

		if(!StringUtils.isEmpty(modType)) { item.setAttribute("modificationType", modType); }
		if(!StringUtils.isEmpty(modSite)) {item.setAttribute("experimentalNterminalSite", modSite); }
		if(!StringUtils.isEmpty(category)) { item.setAttribute("category", category); }
		if(!StringUtils.isEmpty(category)) { item.setAttribute("experimentalNterminalSequence", seq); }


		Item publication = publications.get(refNo);

		if(publication == null) {
			publication = createItem("Publication");
			publication.setAttribute("pubMedId", pmid);       
			publications.put(refNo, publication);		
		}
		item.setReference("publication", publication);      

		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}	

		return item;
	}


	private void processProteinInfo(Connection connection) throws SQLException, ObjectStoreException, Exception {

		System.out.println("Processing Proteins Info data...");
		ResultSet res = PROCESSOR.getProteinInfoResults(connection);

		while (res.next()) {

			String featureNo = res.getString("dbentity_id");

			String molwt = res.getString("molecular_weight");
			String pi = res.getString("pi");

			String ntermseq = res.getString("n_term_seq");
			String ctermseq = res.getString("c_term_seq");
			String gravyScore = res.getString("gravy_score");
			String aromaticityScore = res.getString("aromaticity_score");
			String cai = res.getString("cai");
			String codonBias = res.getString("codon_bias");
			String fopScore = res.getString("fop_score");

			String ala = res.getString("ala");
			String cys = res.getString("cys");
			String asp = res.getString("asp");
			String glu = res.getString("glu");
			String phe = res.getString("phe");
			String gly = res.getString("gly");
			String his = res.getString("his");

			String ile = res.getString("ile");
			String lys = res.getString("lys");
			String leu = res.getString("leu");
			String met = res.getString("met");
			String asn = res.getString("asn");
			String pro = res.getString("pro");		
			String gln = res.getString("gln");
			String arg = res.getString("arg");
			String ser = res.getString("ser");		
			String thr = res.getString("thr");
			String val = res.getString("val");			
			String trp = res.getString("trp");
			String tyr = res.getString("tyr");

			String carbon = res.getString("carbon");
			String hydrogen = res.getString("hydrogen");
			String nitrogen = res.getString("nitrogen");
			String oxygen = res.getString("oxygen");
			String sulphur = res.getString("sulfur");
			String instabilityIndex = res.getString("instability_index");
			String allCysHalf = res.getString("all_cys_ext_coeff");
			String noCysHalf = res.getString("no_cys_ext_coeff");
			String aliphaticIndex = res.getString("aliphatic_index");

			Item protein = proteins.get(featureNo);		

			if (molwt != null) protein.setAttribute("molecularWeight", molwt);		       									
			if (pi != null) protein.setAttribute("pI", pi);					
			if (fopScore != null) protein.setAttribute("fopScore", fopScore);					
			if (gravyScore != null) protein.setAttribute("gravyScore", gravyScore);					
			if (aromaticityScore != null) protein.setAttribute("aromaticityScore", aromaticityScore);					
			if (cys != null) protein.setAttribute("cys", cys);					
			if (gln != null) protein.setAttribute("gln", gln);				
			if (glu != null) protein.setAttribute("glu", glu);					
			if (gly != null) protein.setAttribute("gly", gly);					
			if (his != null) protein.setAttribute("his", his);					
			if (ile != null) protein.setAttribute("ile", ile);					
			if (leu != null) protein.setAttribute("leu", leu);					
			if (lys != null) protein.setAttribute("lys", lys);
			if (met != null) protein.setAttribute("met", met);					
			if (phe != null) protein.setAttribute("phe", phe);					
			if (pro != null) protein.setAttribute("pro", pro);					
			if (ser != null) protein.setAttribute("ser", ser);					
			if (thr != null) protein.setAttribute("thr", thr);					
			if (trp != null) protein.setAttribute("trp", trp);
			if (tyr != null) protein.setAttribute("tyr", tyr);
			if (val != null) protein.setAttribute("val", val);					
			if (ala != null) protein.setAttribute("ala", ala);					
			if (arg != null) protein.setAttribute("arg", arg);					
			if (asn != null) protein.setAttribute("asn", asn);					
			if (asp != null) protein.setAttribute("asp", asp);					
			if (ntermseq != null) protein.setAttribute("ntermseq", ntermseq);					
			if (ctermseq != null) protein.setAttribute("ctermseq", ctermseq);					
			if (cai != null) protein.setAttribute("cai", cai);					
			if (codonBias != null) protein.setAttribute("codonBias", codonBias);							
			if (carbon != null) protein.setAttribute("carbon", carbon);				
			if (hydrogen != null) protein.setAttribute("hydrogen", hydrogen);					
			if (nitrogen != null) protein.setAttribute("nitrogen", nitrogen);					
			if (oxygen != null) protein.setAttribute("oxygen", oxygen);					
			if (sulphur != null) protein.setAttribute("sulphur", sulphur);					
			if (instabilityIndex != null) protein.setAttribute("instabilityIndex", instabilityIndex);					
			if (allCysHalf != null) protein.setAttribute("allCysHalf", allCysHalf);					
			if (noCysHalf != null) protein.setAttribute("noCysHalf", noCysHalf);					
			if (aliphaticIndex != null) protein.setAttribute("aliphaticIndex", aliphaticIndex);

		}

	}


	private void processProteinAbundance(Connection connection) throws SQLException,
	ObjectStoreException, Exception {

									
		System.out.println("Processing Proteins Abundance data...");
		ResultSet res = PROCESSOR.getProteinAbundanceResults(connection);
		while (res.next()) {

			String featureNo = res.getString("dbentity_id");
			String annotationId = res.getString("annotation_id");
			String units = res.getString("data_unit");
			String abundance = res.getString("data_value");
			String pmid = res.getString("pmid");
			String parentpmid = res.getString("parent_pmid");
			String refNo = res.getString("referencedbentity");
			String parentrefNo = res.getString("original_referencedbentity");
			String experiment = "abundance";
			String strainBackground = res.getString("strain_background");
			
			String strain = "";
			if (StringUtils.isNotEmpty(strainBackground)) {
				String q[] = strainBackground.split("_");
				int len = q.length;
				if(len == 3){
					 strain = q[len-1];
				}else{
					 strain = "Other";
				}
			}

			String assay = res.getString("assay");
			String media = res.getString("media");
			String process = res.getString("process");
			String foldchange = res.getString("fold_change");
			
			String chemical = res.getString("chemical");
			String conc_unit = res.getString("concentration_unit");
			String conc_value = res.getString("concentration_value");

			String treatment = "";
			if(conc_value != null && conc_unit != null && chemical !=null) {
				treatment = conc_value+" "+conc_unit+" "+chemical;
			}else if(chemical !=null) {
				treatment = chemical;
			}else if(process != null){
				treatment = process;
			}else {
				treatment = "untreated";
			}
			
			String time_unit = res.getString("time_unit");
			String time_value = res.getString("time_value");
			String treatment_time = "";
			if(time_unit != null & time_value != null) {
				treatment_time = time_value+" "+time_unit;
			}
			
			String median_value = res.getString("median_value");
			String median_abs_dev_value = res.getString("median_abs_dev_value");

			Item protein = proteins.get(featureNo);
			Item pmods = getProteinAbundance(abundance, pmid, refNo, units, annotationId, treatment, treatment_time,
					foldchange, process, media, assay, strain, parentpmid, parentrefNo);
			
			if(StringUtils.isNotEmpty(median_value)){ protein.setAttribute("median", median_value);}
			if(StringUtils.isNotEmpty(median_abs_dev_value)){ protein.setAttribute("MAD", median_abs_dev_value);}
			if(StringUtils.isNotEmpty(units)){ protein.setAttribute("units", units);}
			protein.addToCollection("proteinAbundance", pmods.getIdentifier());

		}

	}


	private Item getProteinAbundance(String abundance, String pmid, String refNo, 
			String units, String annotation_id, String treatment, String treatment_time,String fold_change, String process,
			String media, String assay, String strain, String parentpmid, String parentrefNo) throws ObjectStoreException {

		Item item = createItem("ProteinAbundance");

		if(StringUtils.isNotEmpty(abundance)){ item.setAttribute("abundance", abundance);}

		if(StringUtils.isNotEmpty(fold_change)){ item.setAttribute("foldChange", fold_change);}	
		if(StringUtils.isNotEmpty(treatment)){ item.setAttribute("treatment", treatment);}	
		if(StringUtils.isNotEmpty(treatment_time)){ item.setAttribute("treatmentTime", treatment_time);}	
		if(StringUtils.isNotEmpty(process)){ item.setAttribute("assay", process);}	
		if(StringUtils.isNotEmpty(media)){ item.setAttribute("media", media);}	
		if(StringUtils.isNotEmpty(assay)){ item.setAttribute("visualization", assay);}	
		if(StringUtils.isNotEmpty(strain)){ item.setAttribute("strainBackground", strain);}	
		if(StringUtils.isNotEmpty(units)){ item.setAttribute("units", units);}	
		
		item.setAttribute("source", "SGD");

		Item publication = publications.get(parentrefNo);
		if (publication == null) {
			publication = createItem("Publication");			
			publication.setAttribute("pubMedId", pmid);			 
			publications.put(parentrefNo, publication);			
		}
		item.setReference("publication", publication);   
		
		
		Item parentpublication = publications.get(refNo);
		if (parentpublication == null) {
			parentpublication = createItem("Publication");			
			parentpublication.setAttribute("pubMedId", pmid);			 
			publications.put(refNo, parentpublication);			
		}
		item.setReference("origPublication", parentpublication);   

		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}
		
	
		return item;
	}


	private String getFixedChrName(String name) throws Exception {

		String newname = "";

		if (name.equalsIgnoreCase("Chromosome_Mito")) {
			newname = "chrmt";
		}else if(name.equalsIgnoreCase("2-micron_plasmid")){
			newname = name;
		}else if(name.equalsIgnoreCase("Chromosome_I")){
			newname ="chrI";
		}else if(name.equalsIgnoreCase("Chromosome_II")){
			newname ="chrII";
		}else if(name.equalsIgnoreCase("Chromosome_III")){
			newname ="chrIII";
		}else if(name.equalsIgnoreCase("Chromosome_IV")){
			newname ="chrIV";
		}else if(name.equalsIgnoreCase("Chromosome_V")){
			newname ="chrV";
		}else if(name.equalsIgnoreCase("Chromosome_VI")){
			newname ="chrVI";
		}else if(name.equalsIgnoreCase("Chromosome_VII")){
			newname ="chrVII";
		}else if(name.equalsIgnoreCase("Chromosome_VIII")){
			newname ="chrVIII";
		}else if(name.equalsIgnoreCase("Chromosome_IX")){
			newname ="chrIX";
		}else if(name.equalsIgnoreCase("Chromosome_X")){
			newname ="chrX";
		}else if(name.equalsIgnoreCase("Chromosome_XI")){
			newname ="chrXI";
		}else if(name.equalsIgnoreCase("Chromosome_XII")){
			newname ="chrXII";
		}else if(name.equalsIgnoreCase("Chromosome_XIII")){
			newname ="chrXIII";
		}else if(name.equalsIgnoreCase("Chromosome_XIV")){
			newname ="chrXIV";
		}else if(name.equalsIgnoreCase("Chromosome_XV")){
			newname ="chrXV";
		}else if(name.equalsIgnoreCase("Chromosome_XVI")){
			newname ="chrXVI";
		}
		return newname;

	}

	private Item getChildItem(String childType) throws ObjectStoreException {

		Item item = null;

		if (childType.equalsIgnoreCase("CDS")) {
			item = createItem("CDS");
		} else if (childType.equalsIgnoreCase("intron")) {
			item = createItem("Intron");
		} else if (childType.equalsIgnoreCase("five_prime_UTR_intron")) {
			item = createItem("FivePrimeUTRIntron");
		} else if (childType.equalsIgnoreCase("plus_1_translational_frameshift")) {
			item = createItem("Plus1TranslationalFrameshift");
		} else if (childType.equalsIgnoreCase("ARS consensus sequence")) {
			item = createItem("ARSConsensusSequence");
		} else if (childType.equalsIgnoreCase("binding_site")) {
			item = createItem("BindingSite");
		} else if (childType.equalsIgnoreCase("insertion")) {
			item = createItem("Insertion");
		} else if (childType.equalsIgnoreCase("repeat_region")) {
			item = createItem("RepeatRegion");
		} else if (childType.equalsIgnoreCase("noncoding_exon")) {
			item = createItem("NoncodingExon");
		} else if (childType.equalsIgnoreCase("external_transcribed_spacer_region")) {
			item = createItem("ExternalTranscribedSpacerRegion");
		} else if (childType.equalsIgnoreCase("internal_transcribed_spacer_region")) {
			item = createItem("InternalTranscribedSpacerRegion");
		} else if (childType.equalsIgnoreCase("non_transcribed_region")) {
			item = createItem("NonTranscribedRegion");
		} else if (childType.equalsIgnoreCase("centromere_DNA_Element_I")) {
			item = createItem("CentromereDNAElementI");
		} else if (childType.equalsIgnoreCase("centromere_DNA_Element_II")) {
			item = createItem("CentromereDNAElementII");
		} else if (childType.equalsIgnoreCase("centromere_DNA_Element_III")) {
			item = createItem("CentromereDNAElementIII");
		} else if (childType.equalsIgnoreCase("intein_encoding_region")) {
			item = createItem("InteinEncodingRegion");
		} else if (childType.equalsIgnoreCase("ARS_consensus_sequence")) {
			item = createItem("ARSConsensusSequence");
		}else if (childType.equalsIgnoreCase("Y_prime_element")) {
			item = createItem("YPrimeElement");
		}else if (childType.equalsIgnoreCase("X_element_combinatorial_repeat")) {
			item = createItem("XElementCombinatorialRepeat");
		}else if (childType.equalsIgnoreCase("X_element")) {
			item = createItem("XElement");
		}else if (childType.equalsIgnoreCase("telomeric_repeat")) {
			item = createItem("TelomericRepeat");
		}else if (childType.equalsIgnoreCase("uORF")) {
			item = createItem("uORF");
		}else if (childType.equalsIgnoreCase("W_region")) {
			item = createItem("W_region");
		}else if (childType.equalsIgnoreCase("X_region")) {
			item = createItem("X_region");
		}else if (childType.equalsIgnoreCase("Y_region")) {
			item = createItem("Y_region");
		}else if (childType.equalsIgnoreCase("Z1_region")) {
			item = createItem("Z1_region");
		}else if (childType.equalsIgnoreCase("Z2_region")) {
			item = createItem("Z2_region");
		}

		return item;

	}

	private void processChromosomeSequences(Connection connection)
			throws SQLException, ObjectStoreException, Exception {

		System.out.println("Processing ChromosomeSequence...");
		ResultSet res = PROCESSOR.getChromosomeSequenceResults(connection);

		while (res.next()) {
			String featureNo = res.getString("contig_id");
			String chromosomeNo = res.getString("format_name");
			String feature_type = res.getString("display_name");
			String residues = res.getString("residues");
			String length = res.getString(5);

			String fixed_chromosome_no = getFixedChrName(chromosomeNo);

			if (feature_type.equalsIgnoreCase("chromosome")) {
				Item chr = createItem("Chromosome");
				chr.setAttribute("primaryIdentifier", fixed_chromosome_no);
				chr.setReference("organism", organism);
				chr.setAttribute("length", length);
				chr.setAttribute("featureType", feature_type);

				Item seq = createItem("Sequence");
				seq.setAttribute("residues", residues);
				seq.setAttribute("length", length);

				try {
					store(seq);
				} catch (ObjectStoreException e) {
					throw new ObjectStoreException(e);
				}

				chr.setReference("sequence", seq.getIdentifier());
				chromosomes.put(fixed_chromosome_no, chr.getIdentifier());

				try {
					store(chr);
				} catch (ObjectStoreException e) {
					throw new ObjectStoreException(e);
				}

			} else if (feature_type.equalsIgnoreCase("plasmid")) {

				Item item = createItem("Plasmid");
				item.setAttribute("primaryIdentifier", fixed_chromosome_no); 
				item.setReference("organism", organism);
				plasmids.put(fixed_chromosome_no, item.getIdentifier());
				try {
					store(item);
				} catch (ObjectStoreException e) {
					throw new ObjectStoreException(e);
				}
			}

		}

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

	private void storeTranscripts() throws ObjectStoreException {
		for (Item transcript : transcripts.values()) {
			try {
				store(transcript);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}

	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storeProteins() throws ObjectStoreException {
		for (Item protein : proteins.values()) {
			try {
				store(protein);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}
	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storeInteractionTypes() throws ObjectStoreException {
		for (Item type : interactiontype.values()) {
			try {
				store(type);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}


	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storeInteractionExperiments() throws ObjectStoreException {
		for (Item exp : experimenttype.values()) {
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

	private void storeInteractionDetails() throws ObjectStoreException {
		for (Item det : interactiondetail.values()) {
			try {
				store(det);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}


	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storePublications() throws ObjectStoreException {
		for (Item pub : publications.values()) {
			try {
				store(pub);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}

	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storePathways() throws ObjectStoreException {
		for (Item path : pathways.values()) {
			try {
				store(path);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
		for (Item ps : pathwaysummarys.values()) {
			try {
				store(ps);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}

	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storeInteractions() throws ObjectStoreException {
		for (Item intact : interactions.values()) {
			try {
				store(intact);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}

	/**
	 * 
	 * @throws ObjectStoreException
	 */

	private void storePhenotypes() throws ObjectStoreException {
		for (Item pheno : phenotypes.values()) {
			try {
				store(pheno);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}


	/*private void storePhenotypeAnnotations() throws ObjectStoreException {
		for (Item pheno : phenotypeannots.values()) {
			try {
				store(pheno);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}*/

	private void storeProteinAbundance() throws ObjectStoreException {
		for (Item pa : proteinAbundance.values()) {
			try {
				store(pa);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}
	
	private void storeRegulationSummaries() throws ObjectStoreException {
		for (Item rs : regulationSummary.values()) {
			try {
				store(rs);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processAllPubs(Connection connection) throws SQLException,
	ObjectStoreException {

		String prevReferenceNo = "";
		String prevPubMedId = "";
		String prevTitle = "";
		String prevCitation = "";
		String prevJournal = "";
		String prevVolume = "";
		String prevPages = "";
		String prevYear = "";
		String prevIssue = "";
		String prevAbst = "";
		String prevDbxRef = "";
		String prevDateCreated = "";

		ArrayList<String> hm = new ArrayList<String>();

		Item gene = null;
		boolean firstrow = true;

		System.out.println("Processing All Publications with Topics...");
		ResultSet res = PROCESSOR.getPubAllResults(connection);

		while (res.next()) {

			String referenceNo = res.getString("dbentity_id");
			String pubMedId = res.getString("pmid");
			String title = res.getString("title");
            String volume = res.getString("volume");
            String pages = res.getString("page");
            String year = res.getString("year");
            String issue = res.getString("issue");
			String citation = res.getString("citation");
			String topic = res.getString("topic");
			String journal = res.getString("med_abbr");
			String abstractText = res.getString("text");
			String dbxrefid = res.getString("sgdid");
			String date_created = res.getString("date_created");

			if (firstrow) {
				prevReferenceNo = referenceNo;
				firstrow = false;
			}

			if (!referenceNo.equalsIgnoreCase(prevReferenceNo)) {
				getPub(prevReferenceNo, prevTitle, prevPubMedId, prevCitation,
						hm, prevJournal, prevVolume, prevPages, prevYear,
						prevIssue, prevAbst, prevDbxRef, prevDateCreated);
				hm.clear();
			}

			if (topic != null) {
				hm.add(new String(topic));
			}

			prevReferenceNo = referenceNo;
			prevTitle = title;
			prevPubMedId = pubMedId;
			prevCitation = citation;
			prevJournal = journal;
			prevVolume = volume;
			prevYear = year;
			prevIssue = issue;
			prevPages = pages;
			prevAbst = abstractText;
			prevDbxRef = dbxrefid;
			prevDateCreated = date_created;

		}

		// process the very last reference group
		getPub(prevReferenceNo, prevTitle, prevPubMedId, prevCitation, hm,
				prevJournal, prevVolume, prevPages, prevYear, prevIssue,
				prevAbst, prevDbxRef, prevDateCreated);
		hm.clear();

	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processPubsWithFeatures(Connection connection)
			throws SQLException, ObjectStoreException {

		String prevGeneFeatureNo = "";
		String prevReferenceNo = "";
		String prevPubMedId = "";
		String prevTitle = "";
		String prevCitation = "";
		String prevJournal = "";
		String prevVolume = "";
		String prevPages = "";
		String prevYear = "";
		String prevIssue = "";
		String prevDbxRef = "";
		String prevDateCreated = "";

		ArrayList<String> hm = new ArrayList<String>();

		Item gene = null;
		boolean firstrow = true;
		System.out.println("Processing Publications With Chromosomal Features...");
		ResultSet res = PROCESSOR.getPubWithFeaturesResults(connection);

		while (res.next()) {

			String referenceNo = res.getString("referenceFeatureNo");
			String geneFeatureNo = res.getString("featureNo");
			String pubMedId = res.getString("pmid");
			String title = res.getString("title");
			String citation = res.getString("citation");
			String topic = res.getString("topic");
			String journal = res.getString("med_abbr");
			String volume = res.getString("volume");
			String pages = res.getString("page");
			String year = res.getString("year");
			String issue = res.getString("issue");
			String dbxrefid = res.getString("sgdid");
			String date_created = res.getString("date_created");

			if (!geneFeatureNo.equalsIgnoreCase(prevGeneFeatureNo)) {

				if (!firstrow) {

					gene = genes.get(geneFeatureNo);

					if(gene != null) {

						getPubAnnot(prevReferenceNo, prevTitle, prevPubMedId,
								prevCitation, hm, gene, prevJournal, prevVolume,
								prevPages, prevYear, prevIssue, prevDbxRef, prevDateCreated);
						hm.clear();
						prevReferenceNo = referenceNo;
					}
				}

			}

			if (firstrow) {
				prevReferenceNo = referenceNo;
				firstrow = false;
			}

			if (!referenceNo.equalsIgnoreCase(prevReferenceNo)) {

				if(gene != null) {
					getPubAnnot(prevReferenceNo, prevTitle, prevPubMedId,
							prevCitation, hm, gene, prevJournal, prevVolume,
							prevPages, prevYear, prevIssue, prevDbxRef, prevDateCreated);
					hm.clear();
				}
			}

			if (topic != null) {
				hm.add(new String(topic));
			}

			prevGeneFeatureNo = geneFeatureNo;
			prevReferenceNo = referenceNo;
			prevTitle = title;
			prevPubMedId = pubMedId;
			prevCitation = citation;
			prevJournal = journal;
			prevVolume = volume;
			prevYear = year;
			prevIssue = issue;
			prevPages = pages;
			prevDbxRef = dbxrefid;

		}

		// process the very last reference group
		if(gene != null) {
			getPubAnnot(prevReferenceNo, prevTitle, prevPubMedId, prevCitation, hm,
					gene, prevJournal, prevVolume, prevPages, prevYear, prevIssue, prevDbxRef, prevDateCreated);
			hm.clear();
		}

	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processPhenotypeSummary(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing Phenotype Summary....");
		ResultSet res = PROCESSOR.getPhenotypeSummary(connection);	

		while (res.next()) {

			String featureNo = res.getString("dbentity_id");
			String summary = res.getString("text");

			Item gene = genes.get(featureNo);		
			gene.setAttribute("phenotypeSummary", summary);

		}
	}
	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processFunctionSummary(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing Function Summary....");
		ResultSet res = PROCESSOR.getFunctionSummary(connection);	

		while (res.next()) {

			String featureNo = res.getString("dbentity_id");
			String summary = res.getString("text");

			Item gene = genes.get(featureNo);		
			gene.setAttribute("functionSummary", summary);

		}
	}
	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processGeneSummary(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing Gene Summary....");
		ResultSet res = PROCESSOR.getGeneSummary(connection);	

		while (res.next()) {

			String featureNo = res.getString("dbentity_id");
			String summary = res.getString("text");

			Item gene = genes.get(featureNo);		
			gene.setAttribute("geneSummary", summary);

		}
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */
	private void processRegulationSummary(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing Regulation Summary....");
		ResultSet res = PROCESSOR.getRegulationSummary(connection);	
		while (res.next()) {
			String featureNo = res.getString("dbentity_id");
			String summary = res.getString("text");
			String refNo = res.getString("reference_id");
			String pmid = res.getString("pmid");
			getRegulationSummary(featureNo, summary, refNo, pmid);
		}

	}
	
	private void getRegulationSummary(String featureNo, String text, String refNo, String pmid) 
			throws ObjectStoreException {

		Item gene = genes.get(featureNo);
		Item rs = regulationSummary.get(featureNo);

		if(rs == null){
			rs = createItem("RegulationSummary");
			rs.setAttribute("summaryParagraph", text);	
			Item publication = publications.get(refNo);
			if(publication == null) {
				publication = createItem("Publication");
				publication.setAttribute("pubMedId", pmid);     
				publications.put(refNo, publication);                		
			}
			rs.addToCollection("publications", publication);
			regulationSummary.put(featureNo, rs);
			gene.setReference("regulationSummary", rs.getIdentifier());
		}else{				
			Item publication = publications.get(refNo);
			if(publication == null) {
				publication = createItem("Publication");
				publication.setAttribute("pubMedId", pmid);     
				publications.put(refNo, publication);                		
			}
			rs.addToCollection("publications", publication);		
		}

	}
	
	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processPhysicalInteractions(Connection connection)
			throws SQLException, ObjectStoreException {

		String dsId = getBioGridDataSet();
		int count = 0;
		
		System.out.println("Processing Physical Interactions...");
		ResultSet res = PROCESSOR.getPhysicalInteractionResults(connection);

		while (res.next()) {
			count++;
			String geneFeatureName = res.getString("dbentity1_id");
			Item gene = genes.get(geneFeatureName);

			String interactionNo = res.getString("annotation_id");
			String referenceNo = res.getString("reference_id");
			String interactionType = "physical interactions";
			String experimentType = res.getString("biogrid_experimental_system");
			String annotationType = res.getString("annotation_type");
			String modification = res.getString("modification");

			String interactingGeneFeatureName = res.getString("dbentity2_id");
			Item interactingGene = genes.get(interactingGeneFeatureName);

			String action = res.getString("bait_hit");
			String[] a = action.split("-");
			String role1 = a[0];
			String role2 = a[1];
			String source = res.getString("display_name");
			String phenotype = ""; //res.getString("phenotype");
			String citation = res.getString("citation");
			String pubmed = res.getString("pmid");
			String title = res.getString("title");
			String volume = res.getString("volume");
			String page = res.getString("page");
			String year = res.getString("year");
			String issue = res.getString("issue");
			String abbreviation = res.getString("med_abbr");
			String firstAuthor = res.getString("first_author");
			String dbxrefid = res.getString("sgdid");
			String note = res.getString("note");

			String interactionRefId = getInteraction(interactionNo,
					referenceNo, interactionType, experimentType,
					annotationType, modification, interactingGene, role1, source,
					phenotype, citation, gene, pubmed, title, volume, page,
					year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);

			//store the reverse relationship so that template changes do not have to be made
			if(!geneFeatureName.equals(interactingGeneFeatureName)) {
				String interactionRefId2 = getInteraction(interactionNo,
						referenceNo, interactionType, experimentType,
						annotationType, modification, gene, role2, source,
						phenotype, citation, interactingGene, pubmed, title, volume, page,
						year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);
			}

			/*String interactionRefId = getInteractionNew(interactionNo,
					referenceNo, interactionType, experimentType,
					annotationType, modification, interactingGene, role1, role2, source,
					phenotype, citation, gene, pubmed, title, volume, page,
					year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);*/

		}
		System.out.println("physical interaction count is : " + count);
	}

	/**
	 * 
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processGeneticInteractions(Connection connection)
			throws SQLException, ObjectStoreException {

		String dsId = getBioGridDataSet();
		int count = 0;

		System.out.println("Processing Genetic Interactions...");
		ResultSet res = PROCESSOR.getGeneticInteractionResults(connection);

		while (res.next()) {
			count++;
			String geneFeatureName = res.getString("dbentity1_id");

			Item gene = genes.get(geneFeatureName); //can save on look-ups here

			String interactionNo = res.getString("annotation_id");
			String referenceNo = res.getString("reference_id");
			String interactionType = "genetic interactions"; //res.getString("interaction_type");
			String experimentType = res.getString("biogrid_experimental_system");
			String annotationType = res.getString("annotation_type");
			String modification = ""; //res.getString("modification");

			String interactingGeneFeatureName = res.getString("dbentity2_id");
			Item interactingGene = genes.get(interactingGeneFeatureName);

			String action = res.getString("bait_hit");
			String[] a = action.split("-");
			String role1 = a[0];
			String role2 = a[1];
			String source = res.getString("source");
			String phenotype = res.getString("phenotype");
			String citation = res.getString("citation");
			String pubmed = res.getString("pmid");
			String title = res.getString("title");
			String volume = res.getString("volume");
			String page = res.getString("page");
			String year = res.getString("year");
			String issue = res.getString("issue");
			String abbreviation = res.getString("med_abbr");
			String firstAuthor = res.getString("first_author");
			String dbxrefid = res.getString("sgdid");
			String note = res.getString("note");

			String interactionRefId = getInteraction(interactionNo,
					referenceNo, interactionType, experimentType,
					annotationType, modification, interactingGene, role1, source,
					phenotype, citation, gene, pubmed, title, volume, page,
					year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);
			String interactionRefId2 = "";
			if(!geneFeatureName.equals(interactingGeneFeatureName)) {//store the reverse relationship so that template changes do not have to be made; act1 in gene.X or participant.X

				 interactionRefId2 = getInteraction(interactionNo,
						referenceNo, interactionType, experimentType,
						annotationType, modification, gene, role2, source,
						phenotype, citation, interactingGene, pubmed, title, volume, page,
						year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);
			}
			String both = "";
			if(!interactionRefId2.isEmpty()) {
				both = interactionRefId+":"+interactionRefId2;
			}else {
				both = interactionRefId;
			}
			interactionitems.put(interactionNo, both);
			/*String interactionRefId = getInteractionNew(interactionNo,
					referenceNo, interactionType, experimentType,
					annotationType, modification, interactingGene, role1, role2, source,
					phenotype, citation, gene, pubmed, title, volume, page,
					year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);*/

		}
		System.out.println("genetic interaction count is : " + count);
	}

	/**
	 *
	 * @param connection
	 * @throws SQLException
	 * @throws ObjectStoreException
	 */

	private void processGeneticInteractionsWithAlleles(Connection connection)
			throws SQLException, ObjectStoreException {

		System.out.println("Processing Genetic Interactions with Alleles.....");
		ResultSet res = PROCESSOR.getGeneticInteractionWithAllelesResults(connection);
		int count = 0;

		while (res.next()) {

			String interactionNo = res.getString("interaction_id");
			String annotationNo = res.getString("annotation_id");
			String allele1_id =  res.getString("allele1_id");
			String allele2_id = res.getString("allele2_id");
			String sga_score = res.getString("sga_score");
			String pvalue =  res.getString("pvalue");

			Item interaction = null;
			Item interaction2 = null;

			String b = interactionitems.get(interactionNo);
			if(b.contains(":")){
				String[] bo = b.split(":");
				 interaction = interactions.get(bo[0]);
				 interaction2 = interactions.get(bo[1]);
			}else{
				interaction = interactions.get(b);
			}

			if(interaction != null){

				Item allint =  createItem("AlleleInteraction");

				if (StringUtils.isNotEmpty(sga_score)) allint.setAttribute("sgaScore", sga_score);
				if (StringUtils.isNotEmpty(pvalue)) allint.setAttribute("pvalue", pvalue);
				if (StringUtils.isNotEmpty(allele1_id)) {
					Item allele1 = alleles.get(allele1_id);
					allint.setReference("allele1", allele1);
				}
				if (StringUtils.isNotEmpty(allele2_id)) {
					Item allele2 = alleles.get(allele2_id);
					allint.setReference("allele2", allele2);
				}

				try {
					store(allint);
				} catch (ObjectStoreException e) {
					throw new ObjectStoreException(e);
				}

				interaction.addToCollection("alleleinteractions", allint.getIdentifier());
			}

			if(interaction2 != null){

				Item allint2 =  createItem("AlleleInteraction");

				if (StringUtils.isNotEmpty(sga_score)) allint2.setAttribute("sgaScore", sga_score);
				if (StringUtils.isNotEmpty(pvalue)) allint2.setAttribute("pvalue", pvalue);
				if (StringUtils.isNotEmpty(allele1_id)) {
					Item allele1 = alleles.get(allele1_id);
					allint2.setReference("allele1", allele1);
				}
				if (StringUtils.isNotEmpty(allele2_id)) {
					Item allele2 = alleles.get(allele2_id);
					allint2.setReference("allele2", allele2);
				}

				try {
					store(allint2);
				} catch (ObjectStoreException e) {
					throw new ObjectStoreException(e);
				}

				interaction2.addToCollection("alleleinteractions", allint2.getIdentifier());
			}

		}
		System.out.println("genetic allele interaction count is : " + count);
	}

	private void processPhenotypes(Connection connection) throws SQLException,
	ObjectStoreException {

		ResultSet res = PROCESSOR.getPhenotypeResults(connection);
		System.out.println("Processing Phenotypes...");
		while (res.next()) {
			
			String geneFeatureNo = res.getString("dbentity_id");
			String phenotypeAnnotNo = res.getString("annotation_id");	
			String groupNo = res.getString("group_id");
			String experimentType = res.getString("experiment");
			String experimentComment = res.getString("experiment_comment");
			String alleleComment = res.getString("allele_comment");
			String reporterComment = res.getString("reporter_comment");
			String mutantType = res.getString("mutant");
			String qualifier_observable = res.getString("phenotype");
			String strainBackground = res.getString("strain_name");
			String reporter = res.getString("reporter");
			String alleledbentity_id = res.getString("allele");
			String assay = res.getString("assay");
			String details = res.getString("details");
			String pmid = res.getString("pmid");
			String refNo = res.getString("refNo");
			
				
			Array cond_class = res.getArray("condclass");
			String[] str_cond_class = (String[])cond_class.getArray();
			
			Array cond_name = res.getArray("condname");
			String[] str_cond_name = (String[])cond_name.getArray();
			
			Array cond_value = res.getArray("condvalue");
			String[] str_cond_value = (String[])cond_value.getArray();
			
			Array cond_unit = res.getArray("condunit");
			String[] str_cond_unit = (String[])cond_unit.getArray();
			 
			String qualifier = "";
			String observable = "";
			if(qualifier_observable.contains(":")){
				String t[] = qualifier_observable.split(":");
				qualifier = t[1].trim();
				observable = t[0].trim();
			}else{
				observable = qualifier_observable.trim();
			}
		
			String strain = "";
			if (StringUtils.isNotEmpty(strainBackground)) {
				String q[] = strainBackground.split("_");
				int len = q.length;
				if(len == 3){
					 strain = q[len-1];
				}else{
					 strain = "Other";
				}
			}
			Item gene = genes.get(geneFeatureNo);
			
			String chemical = "";
			String condition = "";
			String chemcond = getPhenotypeCondition(str_cond_class, str_cond_name, str_cond_value, str_cond_unit);
			//System.out.println("chemcond is...." + chemcond);
			String cc[] = chemcond.split("_");
			chemical = cc[0];
			condition = cc[1];				
			
			getPhenotype(phenotypeAnnotNo, groupNo, qualifier, observable, experimentType, experimentComment,
					alleleComment, reporterComment, strain, mutantType, reporter, alleledbentity_id, assay, chemical, condition, details, pmid, refNo, gene);
			
		}


	}
	
	
	private String getPhenotypeCondition(String[] condClass, String[] condName, String[] condValue, String[] condUnits) {

		String chemical = " ";
		String condition = " ";
		for (int i=0; i< condClass.length; i++) {
			
			String cc = condClass[i];
			String cn = condName[i];
			String cv = condValue[i];
			String cu = condUnits[i];
			if(cc != null  && cn != null){	
				cc = Character.toUpperCase(cc.charAt(0)) + cc.substring(1);
				if(cc.equals("Chemical")){
				if(cv !=null) chemical += cv;
				if(cu !=null) {
					if(cu.equals("%")){
						chemical += cu;
					}else{
						chemical += " "+cu;
					}
				}
				if(cn !=null) chemical += " "+cn+", ";
				}else if(cc.equals("Temperature")){
					if(cc !=null) condition += cc+": ";
					if(cn !=null) condition += cn;
					if(cv !=null && cu !=null){
						condition += " ("+cv+cu+") ";
					}else{
					 if(cv !=null && !cv.isEmpty()) condition += " ("+cv+") ";
					}
				}else{
					if(cc !=null) condition += cc+": ";
					if(cv !=null) condition += cv;
					if(cu !=null) condition += cu;
					if(cn !=null) condition += " "+cn+", ";
				}
			}		
		}

		String new_chemical = " ";
		String new_condition = " ";
		
		if(condition.indexOf(",") > 0) {
			new_condition = condition.trim().replaceAll(",$","");
		}else{
			new_condition = condition;
		}
		if(chemical.indexOf(",") > 0)  {
			new_chemical = chemical.trim().replaceAll(",$","");
		}else{
			new_chemical = chemical;
		}
		
		return new_chemical+"_"+new_condition;

	}
	

	private void getPhenotype(String phenotypeAnnotNo, String groupNo, String qualifier, String observable, String experimentType, String experimentComment, 
			String alleleComment, String reporterComment, String strain_background, String mutantType, String reporter, String alleleId,
			String assay, String chemical, String condition, String details, String pmid, String refNo, Item gene ) throws ObjectStoreException {

		Item  pheno = createItem("Phenotype");
			
		if (qualifier != null && StringUtils.isNotEmpty(qualifier)) {
			pheno.setAttribute("qualifier", qualifier);
		}	
		if (observable != null && StringUtils.isNotEmpty(observable)) {
			pheno.setAttribute("observable", observable);
		}
		if (experimentType != null && StringUtils.isNotEmpty(experimentType)) {
			pheno.setAttribute("experimentType", experimentType);
		}
		if (experimentComment != null  && StringUtils.isNotEmpty(experimentComment)) {
			pheno.setAttribute("experimentComment", experimentComment);
		}
		if (alleleComment != null && StringUtils.isNotEmpty(alleleComment)) {
			pheno.setAttribute("alleleComment", alleleComment);
		}
		if (reporterComment != null && StringUtils.isNotEmpty(reporterComment)) {
			pheno.setAttribute("reporterComment", reporterComment);
		}
		if (strain_background != null && StringUtils.isNotEmpty(strain_background)) {
			pheno.setAttribute("strainBackground", strain_background);
		}
		if (mutantType != null && StringUtils.isNotEmpty(mutantType)) {
			pheno.setAttribute("mutantType", mutantType);
		}
		if (reporter != null && StringUtils.isNotEmpty(reporter)) {
			pheno.setAttribute("reporter", reporter);
		}
		if (alleleId != null && StringUtils.isNotEmpty(alleleId)) {
			Item allele = alleles.get(alleleId);
			if (allele != null) {
				pheno.setAttribute("allele", allele.getAttribute("name").getValue());
				if (allele.getAttribute("description") != null) {
					String description = allele.getAttribute("description").getValue();
					if( description != null && StringUtils.isNotEmpty(description)) {
						pheno.setAttribute("alleleDescription", description);
					}
				}
			}
		}
		if (assay != null && StringUtils.isNotEmpty(assay)) {
			pheno.setAttribute("assay", assay);
		}
		if (chemical != null && StringUtils.isNotEmpty(chemical)) {
			pheno.setAttribute("chemical", chemical);
		}
		if (condition != null && StringUtils.isNotEmpty(condition)) {
			pheno.setAttribute("condition", condition);
		}
		if (details != null && StringUtils.isNotEmpty(details)) {
			pheno.setAttribute("details", details);
		}
		
		Item pub = getPubPhenotype(refNo, pmid);

		pheno.addToCollection("genes", gene.getIdentifier());
		pheno.addToCollection("publications", pub.getIdentifier());
		pub.addToCollection("phenotypes", pheno.getIdentifier());
		
		String unq = phenotypeAnnotNo+":"+groupNo;
		phenotypes.put(unq, pheno);
	

	}

	
	private String getLocation(Item subject, String chromosomeRefId,
			String startCoord, String stopCoord, String strand)
					throws ObjectStoreException {

		String start = startCoord;
		String end = stopCoord;

		/*
		 * String start = ""; //(strand.equals("-") ? stopCoord : startCoord);
		 * //was C for crick String end = ""; // (strand.equals("-") ?
		 * startCoord : stopCoord); //was C for crick
		 * 
		 * if(strand.equals("-")){ start = stopCoord; end = startCoord;
		 * System.out.println("start and stop should be reveresed.. " + start +
		 * "  " + end); }else{ end = stopCoord; start = startCoord; }
		 */

		/*
		 * if (StringUtils.isEmpty(startCoord)) { start = "0"; } if
		 * (StringUtils.isEmpty(stopCoord)) { end = "0"; }
		 */

		if (!StringUtils.isEmpty(start) && !StringUtils.isEmpty(end)) {
			subject.setAttribute("length", getLength(start, end));
		}

		Item location = createItem("Location");

		if (!StringUtils.isEmpty(start))
			location.setAttribute("start", start);
		if (!StringUtils.isEmpty(end))
			location.setAttribute("end", end);
		if (!StringUtils.isEmpty(strand))
			location.setAttribute("strand", strand);

		location.setReference("feature", subject);
		location.setReference("locatedOn", chromosomeRefId);

		try {
			store(location);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}
		return location.getIdentifier();
	}

	private String getLength(String start, String end)
			throws NumberFormatException {
		Integer a = new Integer(start);
		Integer b = new Integer(end);

		// if the coordinates are on the crick strand, they need to be reversed
		// or they result in a negative number
		if (a.compareTo(b) > 0) {
			a = new Integer(end);
			b = new Integer(start);
		}

		Integer length = new Integer(b.intValue() - a.intValue());
		return length.toString();
	}

	private String getChromosome(String identifier) throws ObjectStoreException {
		if (StringUtils.isEmpty(identifier)) {
			return null;
		}
		String refId = chromosomes.get(identifier);
		if (refId == null) {
			Item item = createItem("Chromosome");
			item.setAttribute("primaryIdentifier", identifier);
			item.setReference("organism", organism);
			refId = item.getIdentifier();
			chromosomes.put(identifier, refId);
			try {
				store(item);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
		return refId;
	}

	private String getBioGridDataSet() throws ObjectStoreException {

		Item item = createItem("DataSet");
		item.setAttribute("name", "BioGRID interaction data set");


		Item ds = createItem("DataSource");
		ds.setAttribute("name", "BioGRID");
		ds.addToCollection("dataSets", item.getIdentifier());
		try {
			store(item);
			store(ds);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}

		return item.getIdentifier();

	}
	private String getPlasmid(String identifier) throws ObjectStoreException { // String
		// id,
		if (StringUtils.isEmpty(identifier)) {
			return null;
		}
		String refId = plasmids.get(identifier);
		if (refId == null) {
			Item item = createItem("Plasmid");
			item.setAttribute("primaryIdentifier", identifier);
			item.setReference("organism", organism);
			refId = item.getIdentifier();
			plasmids.put(identifier, refId);
			try {
				store(item);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
		return refId;
	}

	private String getSequence(String id, String residues, String length)
			throws ObjectStoreException {
		if (StringUtils.isEmpty(id)) {
			return null;
		}
		String refId = sequences.get(id);
		if (refId == null) {
			Item item = createItem("Sequence");
			item.setAttribute("residues", residues);
			item.setAttribute("length", length);
			// item.setReference("organism", organism);
			refId = item.getIdentifier();
			sequences.put(id, refId);
			try {
				store(item);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
		return refId;
	}

	private String getSynonym(String subjectId, String type, String value)
			throws ObjectStoreException {
		String key = subjectId + type + value;
		if (StringUtils.isEmpty(value)) {
			return null;
		}
		String refId = synonyms.get(key);
		if (refId == null) {
			Item syn = createItem("Synonym");
			syn.setReference("subject", subjectId);
			//syn.setAttribute("type", type);
			syn.setAttribute("value", value);
			refId = syn.getIdentifier();
			// synonyms.get(key);
			try {
				store(syn);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
		}
		return refId;
	}

	/**
	 *
	 * @param interactionNo
	 * @param referenceNo
	 * @param interactionType
	 * @param experimentType
	 * @param annotationType
	 * @param modification
	 * @param interactingGene
	 * @param action
	 * @param source
	 * @param phenotype
	 * @param citation
	 * @param gene
	 * @param pubMedId
	 * @param title
	 * @param volume
	 * @param page
	 * @param year
	 * @param issue
	 * @param journal
	 * @param dsetIdentifier
	 * @param firstAuthor
	 * @param dbxrefid
	 * @return
	 * @throws ObjectStoreException
	 */
	private String getInteractionNew(String interactionNo, String referenceNo,
								  String interactionType, String experimentType,
								  String annotationType, String modification, Item interactingGene, String role1, String role2,
								  String source, String phenotype, String citation, Item gene,
								  String pubMedId, String title, String volume, String page,
								  String year, String issue, String journal, String dsetIdentifier, String firstAuthor, String dbxrefid, String note)
			throws ObjectStoreException {

		Item item = getInteractionItem(gene.getIdentifier(), interactingGene.getIdentifier());
		Item detail = createItem("InteractionDetail");

		detail.setAttribute("type", interactionType);
		detail.setAttribute("annotationType", annotationType);
		if (StringUtils.isNotEmpty(modification)) detail.setAttribute("modification", modification);
		if (StringUtils.isNotEmpty(phenotype)) detail.setAttribute("phenotype", phenotype);
		detail.setAttribute("role1", role1);
		detail.setAttribute("role2", role2);
		detail.addToCollection("allInteractors", interactingGene.getIdentifier());
		detail.addToCollection("dataSets", dsetIdentifier);

		String shortType = interactionType.substring(0, interactionType.indexOf(' '));
		detail.setAttribute("relationshipType", shortType); //interactionType
		String unqName = firstAuthor+"-"+pubMedId+"-"+experimentType;
		if (StringUtils.isNotEmpty(note)) detail.setAttribute("note", note);

		//add publication as experiment type
		Item storedExperimentType = experimenttype.get(unqName);
		if(storedExperimentType == null) {

			storedExperimentType = createItem("InteractionExperiment");
			storedExperimentType.setAttribute("name", unqName);
			experimenttype.put(unqName, storedExperimentType);

			String storedTermId = interactionterms.get(experimentType);
			if (storedTermId != null) {
				storedExperimentType.addToCollection("interactionDetectionMethods", storedTermId );
			} else {
				storedTermId = getInteractionTerm(experimentType);
				storedExperimentType.addToCollection("interactionDetectionMethods", storedTermId );
			}
		}

		//add publication as reference on experiment
		Item storedRef = publications.get(referenceNo);

		if (storedRef != null) {
			storedExperimentType.setReference("publication", storedRef.getIdentifier());
		} else {

			Item pub = createItem("Publication");

			if (StringUtils.isNotEmpty(pubMedId)) {
				pub.setAttribute("pubMedId", pubMedId);
			}
			if (StringUtils.isNotEmpty(dbxrefid)) {
				pub.setAttribute("pubXrefId", dbxrefid);
			}
			if (StringUtils.isNotEmpty(title)) {
				pub.setAttribute("title", title);
			}
			if (StringUtils.isNotEmpty(citation)) {
				pub.setAttribute("citation", citation);
			}
			if (StringUtils.isNotEmpty(journal)) {
				pub.setAttribute("journal", journal);
			}
			if (StringUtils.isNotEmpty(volume)) {
				pub.setAttribute("volume", volume);
			}
			if (StringUtils.isNotEmpty(page)) {
				pub.setAttribute("pages", page);
			}
			if (StringUtils.isNotEmpty(year)) {
				pub.setAttribute("year", year);
			}
			if (StringUtils.isNotEmpty(issue)) {
				pub.setAttribute("issue", issue);
			}
			publications.put(referenceNo, pub);
			storedExperimentType.setReference("publication", pub.getIdentifier());
		}

		detail.setReference("experiment", storedExperimentType.getIdentifier());
		detail.setReference("interaction", item);
		interactiondetail.put(interactionNo, detail);
		/*try {
			store(detail);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}*/

		interactions.put(item.getIdentifier(), item);
		String refId = item.getIdentifier();
		return refId;
	}

	/**
	 * 
	 * @param interactionNo
	 * @param referenceNo
	 * @param interactionType
	 * @param experimentType
	 * @param annotationType
	 * @param modification
	 * @param interactingGene
	 * @param action
	 * @param source
	 * @param phenotype
	 * @param citation
	 * @param gene
	 * @param pubMedId
	 * @param title
	 * @param volume
	 * @param page
	 * @param year
	 * @param issue
	 * @param journal
	 * @param dsetIdentifier
	 * @param firstAuthor
	 * @param dbxrefid
	 * @return
	 * @throws ObjectStoreException
	 */
	private String getInteraction(String interactionNo, String referenceNo,
			String interactionType, String experimentType,
			String annotationType, String modification, Item interactingGene, String action,
			String source, String phenotype, String citation, Item gene,
			String pubMedId, String title, String volume, String page,
			String year, String issue, String journal, String dsetIdentifier, String firstAuthor, String dbxrefid, String note)
					throws ObjectStoreException {

		Item item = getInteractionItem(gene.getIdentifier(), interactingGene.getIdentifier());	
		Item detail = createItem("InteractionDetail");		   

		detail.setAttribute("type", interactionType);		
		detail.setAttribute("annotationType", annotationType);
		if (StringUtils.isNotEmpty(modification)) detail.setAttribute("modification", modification);
		if (StringUtils.isNotEmpty(phenotype)) detail.setAttribute("phenotype", phenotype);
		detail.setAttribute("role1", action);
		detail.addToCollection("allInteractors", interactingGene.getIdentifier());
		detail.addToCollection("dataSets", dsetIdentifier);		

		String shortType = interactionType.substring(0, interactionType.indexOf(' ')); 
		detail.setAttribute("relationshipType", shortType); //interactionType
		String unqName = firstAuthor+"-"+pubMedId+"-"+experimentType;
		if (StringUtils.isNotEmpty(note)) detail.setAttribute("note", note);

		//add publication as experiment type
		Item storedExperimentType = experimenttype.get(unqName);		
		if(storedExperimentType == null) {	

			storedExperimentType = createItem("InteractionExperiment");
			storedExperimentType.setAttribute("name", unqName);	
			experimenttype.put(unqName, storedExperimentType);			

			String storedTermId = interactionterms.get(experimentType);
			if (storedTermId != null) {
				storedExperimentType.addToCollection("interactionDetectionMethods", storedTermId );
			} else {
				storedTermId = getInteractionTerm(experimentType);
				storedExperimentType.addToCollection("interactionDetectionMethods", storedTermId );
			}
		}

		//add publication as reference on experiment
		Item storedRef = publications.get(referenceNo);

		if (storedRef != null) {
			storedExperimentType.setReference("publication", storedRef.getIdentifier());
		} else {

			Item pub = createItem("Publication");

			if (StringUtils.isNotEmpty(pubMedId)) {
				pub.setAttribute("pubMedId", pubMedId);
			}
			if (StringUtils.isNotEmpty(dbxrefid)) {
				pub.setAttribute("pubXrefId", dbxrefid);
			}
			if (StringUtils.isNotEmpty(title)) {
				pub.setAttribute("title", title);
			}
			if (StringUtils.isNotEmpty(citation)) {
				pub.setAttribute("citation", citation);
			}
			if (StringUtils.isNotEmpty(journal)) {
				pub.setAttribute("journal", journal);
			}
			if (StringUtils.isNotEmpty(volume)) {
				pub.setAttribute("volume", volume);
			}
			if (StringUtils.isNotEmpty(page)) {
				pub.setAttribute("pages", page);
			}
			if (StringUtils.isNotEmpty(year)) {
				pub.setAttribute("year", year);
			}
			if (StringUtils.isNotEmpty(issue)) {
				pub.setAttribute("issue", issue);
			}
			publications.put(referenceNo, pub);
			storedExperimentType.setReference("publication", pub.getIdentifier());
		}

		detail.setReference("experiment", storedExperimentType.getIdentifier());	
		detail.setReference("interaction", item);

		try {
			store(detail);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}

		interactions.put(item.getIdentifier(), item);
		//interactionitems.put(interactionNo, item);
		String refId = item.getIdentifier();
		return refId;

	}



	private Item getInteractionItem(String refId, String gene2RefId) throws ObjectStoreException {
		MultiKey key = new MultiKey(refId, gene2RefId);
		Item interaction = interactionsnew.get(key);
		if (interaction == null) {
			interaction = createItem("Interaction");
			interaction.setReference("participant1", refId); //gene1
			interaction.setReference("participant2", gene2RefId); //gene2
		}
		return interaction;
	}


	private String getLiteratureTopic(String topic) throws ObjectStoreException {
		String refId = literatureTopics.get(topic);
		if (refId == null) {
			Item item = createItem("LiteratureTopic");
			item.setAttribute("name", topic);
			try {
				store(item);
			} catch (ObjectStoreException e) {
				throw new ObjectStoreException(e);
			}
			refId = item.getIdentifier();
			literatureTopics.put(topic, refId);
		}
		return refId;
	}

	private void getPub(String referenceNo, String title, String pubMedId,
			String citation, ArrayList hm, String journal, String volume,
			String pages, String year, String issue, String abstractText, String dbxrefid, String datecreated)
					throws ObjectStoreException {

		Item storedRef = publications.get(referenceNo);

		if (storedRef == null) {

			Item item = createItem("Publication");

			if (StringUtils.isNotEmpty(pubMedId)) {
				item.setAttribute("pubMedId", pubMedId);
			}
			if (StringUtils.isNotEmpty(dbxrefid)) {
				item.setAttribute("pubXrefId", dbxrefid);
			}
			if (StringUtils.isNotEmpty(title)) {
				item.setAttribute("title", title);
			}
			if (StringUtils.isNotEmpty(citation)) {
				item.setAttribute("citation", citation);
			}
			if (StringUtils.isNotEmpty(journal)) {
				item.setAttribute("journal", journal);
			}
			if (StringUtils.isNotEmpty(volume)) {
				item.setAttribute("volume", volume);
			}
			if (StringUtils.isNotEmpty(pages)) {
				item.setAttribute("pages", pages);
			}
			if (StringUtils.isNotEmpty(year)) {
				item.setAttribute("year", year);
			}
			if (StringUtils.isNotEmpty(issue)) {
				item.setAttribute("issue", issue);
			}
			if (StringUtils.isNotEmpty(abstractText)) {
				item.setAttribute("abstractText", abstractText);
			}
			if (StringUtils.isNotEmpty(datecreated)) {
				item.setAttribute("dateCreated", datecreated);
			}

			Iterator iter = hm.iterator();
			while (iter.hasNext()) {
				String value = (String) iter.next();
				item.addToCollection("literatureTopics",getLiteratureTopic(value));
			}
			publications.put(referenceNo, item);
		}

	}

	private Item getExistingPub(String referenceNo)
			throws ObjectStoreException {
		Item storedRef = publications.get(referenceNo);
		return storedRef;
	}


	private Item getPubPhenotype(String prevReferenceNo, String pubMedId)
					throws ObjectStoreException {

		//Item storedPheno = phenotypepubs.get(phenoAnnotNo);
		//if (storedPheno != null) {

			Item storedRef = publications.get(prevReferenceNo);

			if (storedRef == null) {
				 storedRef = createItem("Publication");
				if (StringUtils.isNotEmpty(pubMedId)) {
					storedRef.setAttribute("pubMedId", pubMedId);
				}
				publications.put(prevReferenceNo, storedRef);
				//storedPheno.addToCollection("publications", item.getIdentifier());
				//item.addToCollection("phenotypes", storedPheno.getIdentifier());
			}
			return storedRef;

		//}

	}

	/**
	 * 
	 * @param featureNo
	 * @param issue
	 * @param volume
	 * @param pubMedId
	 * @param pages
	 * @param title
	 * @param year
	 * @param citation
	 * @return
	 * @throws ObjectStoreException
	 */
	private void getPubAnnot(String referenceNo, String title, String pubMedId,
			String citation, ArrayList hm, Item gene, String journal,
			String volume, String pages, String year, String issue, String dbxrefid, String datecreated)
					throws ObjectStoreException {

		Item pubAnnot = createItem("PublicationAnnotation");
		Item storedRef = publications.get(referenceNo);

		if (storedRef == null) {
			Item item = createItem("Publication");

			if (StringUtils.isNotEmpty(pubMedId)) {
				item.setAttribute("pubMedId", pubMedId);
			}
			if (StringUtils.isNotEmpty(dbxrefid)) {
				item.setAttribute("pubXrefId", dbxrefid);
			}
			if (StringUtils.isNotEmpty(title)) {
				item.setAttribute("title", title);
			}
			if (StringUtils.isNotEmpty(citation)) {
				item.setAttribute("citation", citation);
			}
			if (StringUtils.isNotEmpty(journal)) {
				item.setAttribute("journal", journal);
			}
			if (StringUtils.isNotEmpty(volume)) {
				item.setAttribute("volume", volume);
			}
			if (StringUtils.isNotEmpty(pages)) {
				item.setAttribute("pages", pages);
			}
			if (StringUtils.isNotEmpty(year)) {
				item.setAttribute("year", year);
			}
			if (StringUtils.isNotEmpty(issue)) {
				item.setAttribute("issue", issue);
			}
			if (StringUtils.isNotEmpty(datecreated)) {
				item.setAttribute("dateCreated", datecreated);
			}

			Iterator iter = hm.iterator();
			while (iter.hasNext()) {
				String value = (String) iter.next();
				pubAnnot.addToCollection("literatureTopics", getLiteratureTopic(value));
			}

			//item.addToCollection("genes", gene.getIdentifier());
			item.addToCollection("entities", gene.getIdentifier());

			publications.put(referenceNo, item);

			pubAnnot.setReference("publication", item.getIdentifier());
			pubAnnot.addToCollection("genes", gene.getIdentifier());

			// in order for the publication enrichment to work, genes needs to
			// have publications directly as a collection, I think.
			gene.addToCollection("publications", item.getIdentifier());

		} else {
			// this means this publication appears for a different gene, lets
			// add the literature topics
			Iterator iter = hm.iterator();
			while (iter.hasNext()) {
				String value = (String) iter.next();
				pubAnnot.addToCollection("literatureTopics", getLiteratureTopic(value));
			}

			pubAnnot.setReference("publication", storedRef.getIdentifier());
			pubAnnot.addToCollection("genes", gene.getIdentifier());

			//storedRef.addToCollection("genes", gene.getIdentifier());
			storedRef.addToCollection("entities", gene.getIdentifier());

			gene.addToCollection("publications", storedRef.getIdentifier());
		}

		try {
			store(pubAnnot);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}

	}

	private String getInteractionTerm(String identifier)  throws ObjectStoreException {
		String refId = interactionterms.get(identifier);
		if (refId != null) {
			return refId;
		}
		Item item = createItem("InteractionTerm");
		item.setAttribute("identifier", identifier);
		interactionterms.put(identifier, item.getIdentifier());
		try {
			store(item);
		} catch (ObjectStoreException e) {
			throw new ObjectStoreException(e);
		}
		return item.getIdentifier();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getDataSetTitle(String taxonId) {
		return DATASET_TITLE;
	}

	/**
	 * Set the data licence for these data.
	 *
	 * @param licence should be URI to data licence.
	 */
	public void setLicence(String licence) {
		this.licence = licence;
	}

	/**
	 * Get the data licence for these data.
	 *
	 * @return URI to data licence.
	 */
	public String getLicence() {
		return licence;
	}



}
