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
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import org.apache.log4j.Logger;

/**
 * Queries the sgd oracle database and returns result sets
 * @author Julie Sullivan
 */
public class SgdProcessor
{
	private static final Logger LOG = Logger.getLogger(SgdProcessor.class);  
	private static final String SCHEMA_OWNER = "nex.";

	/**
	 * Return the results of running a query for genes
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getChromosomalFeatureResults(Connection connection)
			throws SQLException {

		String query = "select l.dbentity_id, l.systematic_name, d.sgdid, l.gene_name, l.name_description, s.display_name as feature_type,  l.headline,"
				+ " l.description, l.qualifier, d.dbentity_status"
				+ " from nex.locusdbentity l, nex.contig c, nex.dnasequenceannotation a, nex.so s, nex.dbentity d "
				+ " where C.contig_id = A.contig_id "
				+ " and (C.format_name like 'Chromosome_%' OR  C.format_name like '2-micron%') "
				+ " and A.dbentity_id = L.dbentity_id "
				+ " and S.so_id = A.so_id "
				+ " and a.taxonomy_id = 274901 "
				+ " and a.dna_type = 'GENOMIC' "
				+ " and L.dbentity_id = D.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of running a query for alleles
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getAlleleResults(Connection connection)
			throws SQLException {

		String query = "select ad.dbentity_id as allele, db.sgdid as allelesgdid, la.locus_id as locus, ad.description, db.display_name as allele_name, "
			    + " aa.display_name as alias_name, ar.reference_id, rdb.pmid, s.display_name as allele_class "
				+ "from nex.alleledbentity ad "
				+ "inner join nex.dbentity db on ad.dbentity_id = db.dbentity_id "
				+ "inner join nex.so s on ad.so_id = s.so_id "
				+ "left join nex.locus_allele la on la.allele_id = ad.dbentity_id "
				+ "left join nex.literatureannotation ar on ad.dbentity_id = ar.dbentity_id "
				+ "left join nex.allele_alias aa on ad.dbentity_id = aa.allele_id "
				+ "left join nex.referencedbentity rdb on ar.reference_id = rdb.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of running a query for alleles
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getAlleleAliasesResults(Connection connection)
			throws SQLException {

		String query = "select ad.dbentity_id as allele, db.display_name as allele_name, "
				+ "array_agg (aa.display_name) as alias_name "
				+ "from nex.alleledbentity ad "
				+ "inner join nex.dbentity db on ad.dbentity_id = db.dbentity_id "
				+ "left join nex.allele_alias aa on ad.dbentity_id = aa.allele_id "
				+ "where aa.display_name is not null "
				+ "group by ad.dbentity_id , db.display_name "
				+ "order by ad.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}



	/**
	 * Return the results of running a query for alleles
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getTranscriptResults(Connection connection)
			throws SQLException {

		String query = "select ad.dbentity_id as transcript, db.format_name, start_index, end_index, strand, in_gal, in_ypd, in_ncbi, ar.reference_id, pmid "
				+ "from nex.transcriptdbentity ad "
				+ "inner join nex.dbentity db on ad.dbentity_id = db.dbentity_id "
				+ "inner join nex.dnasequenceannotation dsa on dsa.dbentity_id = ad.dbentity_id "
				+ "left join nex.transcript_reference ar on ad.dbentity_id = ar.transcript_id "
				+ "left join nex.referencedbentity rdb on ar.reference_id = rdb.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}



	/**
	 * Return the results of running a query for NISS 
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getNISS(Connection connection)
			throws SQLException {

		String query = "select l.dbentity_id, l.systematic_name, d.sgdid, l.gene_name, l.name_description, l.headline,"
				+ " l.description, l.qualifier, d.dbentity_status"
				+ " from nex.locusdbentity l, nex.dbentity d"
				+ " where not_in_s288c = true"
				+ " and L.dbentity_id = D.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of running a query for genes and chromosomes
	 * @param connection the connection
	 * @return the results
	 * modified to get Active only since the addition of seq_rel and release tables 9/7/11
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getChromosomalFeatureLocationResults(Connection connection)
			throws SQLException {
		String query = "select c.contig_id, c.format_name, s.display_name as feature_type, l.dbentity_id, l.gene_name, l.systematic_name, a.strand, a.end_index, a.start_index, a.residues, length(a.residues) "
				+ " from nex.locusdbentity l, nex.contig c, nex.dnasequenceannotation a, nex.so s "
				+ " where C.contig_id = A.contig_id "
				+ " and (C.format_name like 'Chromosome_%' OR  C.format_name like '2-micron%')"
				+ " and A.dbentity_id = L.dbentity_id "
				+ " and S.so_id = c.so_id "
				+ " and a.taxonomy_id = 274901 "
				+ " and a.dna_type = 'GENOMIC'";


		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 *
	 * @param connection
	 * @param parent_id
	 * @return
	 * @throws SQLException
	 */
	protected ResultSet getChildrenFeatureLocationResults(Connection connection)
			throws SQLException {

		String query = "select a.dbentity_id as parent_id, "
				+ "s.display_name as parent_type, a2.dbentity_id as child_id, a2.display_name as child_type, "
				+ " ld.systematic_name as child_identifier, d2.sgdid as child_sgdid, "
				+ " a2.contig_end_index as child_end_coord, a2.contig_start_index as child_start_coord, "
				+ " a.strand, c.format_name , a2.residues, length(a2.residues) as seq_length, d.dbentity_status as parent_status, d2.dbentity_status as child_status "
				+ " from   nex.contig c, nex.dnasequenceannotation a, nex.dnasubsequence a2, nex.dbentity d, nex.dbentity d2, nex.locusdbentity ld, nex.so s, nex.taxonomy t "
				+ " where c.contig_id = a.contig_id "
				+ " and (C.format_name like 'Chromosome_%' OR  C.format_name like '2-micron%') "
				+ " and    a.annotation_id = a2.annotation_id "
				+ " and    a.dbentity_id = d.dbentity_id "
				+ " and    a.so_id = s.so_id "
				+ " and    a2.dbentity_id = d2.dbentity_id "
				+ " and    d2.dbentity_id = ld.dbentity_id "
				+ " and    a.taxonomy_id = t.taxonomy_id "
				+ " and    t.taxid = 'TAX:559292'";
				//+ " and    t.display_name = 'Saccharomyces cerevisiae S288c'";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Returns the results of running a query for genes and their sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException  if there is a database problem
	 */
	protected ResultSet getChromosomeSequenceResults(Connection connection)
			throws SQLException {

		String query =  "SELECT c.contig_id, c.format_name, s.display_name, c.residues, length(residues) "
				+ "FROM "+ SCHEMA_OWNER + "contig c, "+ SCHEMA_OWNER + "so s  " 
				+ "WHERE s.so_id = c.so_id "
				+ "AND s.display_name in ('chromosome', 'plasmid') "
				+ "AND taxonomy_id = 274901";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of running a query for protein sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getProteinResults(Connection connection)
			throws SQLException {


		String query = "select db.dbentity_id, db.sgdid, db.format_name, db.display_name, residues, length(residues) - 1 "
				+ " from nex.dbentity db, nex.locusdbentity ldb, nex.proteinsequenceannotation ps "
				+ " where db.dbentity_id = ldb.dbentity_id "
				+ " and ps.dbentity_id = db.dbentity_id"
				+ " and ps.taxonomy_id = 274901";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of running a query for protein sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getProteinHalfLifeResults(Connection connection)
			throws SQLException {

		String query = "select db.dbentity_id, data_value, data_unit, pmid, rdb.dbentity_id as referencedbentity"
				+ " from nex.dbentity db, nex.proteinexptannotation pea, nex.referencedbentity rdb"
				+ " where pea.dbentity_id = db.dbentity_id"
				+ " and pea.reference_id = rdb.dbentity_id"
				+ " and experiment_type = 'half-life'";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of running a query for protein sequence info
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getProteinInfoResults(Connection connection)
			throws SQLException {

		String query = "select db.dbentity_id, molecular_weight, n_term_seq, c_term_seq, pi, cai, codon_bias, fop_score,"
				+ " gravy_score, aromaticity_score, aliphatic_index, instability_index, ala, arg, asn, asp, cys, gln,"
				+ " glu, gly, his, ile, leu, lys, met, phe, pro, ser, thr, trp, tyr, val, hydrogen, sulfur, nitrogen,"
				+ " oxygen, carbon, no_cys_ext_coeff, all_cys_ext_coeff"
				+ " from nex.dbentity db, nex.locusdbentity ldb, nex.proteinsequenceannotation ps, nex.proteinsequence_detail psd"
				+ " where db.dbentity_id = ldb.dbentity_id"
				+ " and ps.dbentity_id = db.dbentity_id"
				+ " and psd.annotation_id = ps.annotation_id"
				+ " and ps.taxonomy_id = 274901";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}




	/**
	 * Return the results of running a query for protein sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getProteinAbundanceResults(Connection connection)
			throws SQLException {
		
		String query = "select db.dbentity_id, pea.annotation_id, data_value, data_unit, rdb2.pmid, rdb1.pmid as parent_pmid, rdb1.dbentity_id as original_referencedbentity, rdb2.dbentity_id as referencedbentity , " 
				+ "median_value, median_abs_dev_value, time_unit, time_value, concentration_value, concentration_unit, " 
				+ "fold_change, g.display_name as process, ef.display_name as media, ec.display_name as assay, chb.display_name as chemical, t.format_name as strain_background "
				+ "from nex.dbentity db " 
				+ "inner join nex.proteinabundanceannotation pea on db.dbentity_id = pea.dbentity_id "
				+ "inner join  nex.referencedbentity rdb1 on pea.reference_id = rdb1.dbentity_id "
				+ "inner join  nex.referencedbentity rdb2 on pea.original_reference_id = rdb2.dbentity_id "
				+ "left join nex.taxonomy t on t.taxonomy_id = pea.taxonomy_id "
				+ "left join nex.go g on g.go_id = pea.process_id "
				+ "left join nex.efo ef on ef.efo_id = pea.media_id " 
				+ "left join nex.eco ec on ec.eco_id = pea.assay_id "
				+ "left join nex.chebi chb on chb.chebi_id = pea.chemical_id ";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of running a query for protein sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getProteinDomainsResults(Connection connection)
			throws SQLException {

		String query = "select db.dbentity_id, start_index, end_index, pd.display_name, interpro_id, description, url_type, date_of_run"
				+ " from nex.dbentity db, nex.proteindomainannotation pda, nex.proteindomain pd, nex.proteindomain_url pdu"
				+ " where pda.dbentity_id = db.dbentity_id"
				+ " and pd.proteindomain_id = pda.proteindomain_id"
				+ " and pdu.proteindomain_id = pd.proteindomain_id"
				+ " and pda.taxonomy_id = 274901";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of running a query for protein sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getProteinModificationResults(Connection connection)
			throws SQLException {

		String query = "select db.dbentity_id, site_index, site_residue, p.display_name, m.format_name as modifier_name, pmid, rdb.dbentity_id as referencedbentity"
				+ " from nex.dbentity db"
				+ " inner join nex.posttranslationannotation pda on pda.dbentity_id = db.dbentity_id"
				+ " inner join nex.referencedbentity rdb on rdb.dbentity_id = pda.reference_id"
				+ " inner join nex.psimod p on  p.psimod_id = pda.psimod_id"
				+ " left join nex.dbentity m on m.dbentity_id = pda.modifier_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}



	/**
	 * Return the results of running a query for all publications associated with chromosomal features  
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getPubWithFeaturesResults(Connection connection)
			throws SQLException {

		String query = "select ldb.dbentity_id as featureNo, r.dbentity_id as referenceFeatureNo, r.pmid, r.fulltext_status, r.title, r.volume,"
				+ " r.page, r.year, r.issue, r.citation, la.topic, j.med_abbr, db.sgdid, db.date_created"
				+ " from nex.literatureannotation la"
				+ " inner join nex.referencedbentity r on la.reference_id = r.dbentity_id"
				+ " inner join nex.locusdbentity ldb on la.dbentity_id  = ldb.dbentity_id"
				+ " inner join nex.dbentity db on ldb.dbentity_id = db.dbentity_id"
				+ " left join nex.journal j on j.journal_id = r.journal_id"
				+ " group by ldb.dbentity_id, r.dbentity_id, r.pmid, r.fulltext_status, r.title, r.volume, r.page, r.year, r.issue, r.citation, la.topic, j.med_abbr, db.sgdid, db.date_created"
				+ " order by ldb.dbentity_id, r.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of running a query for all publications.  
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getPubAllResults(Connection connection)
			throws SQLException {

		String query = "select r.dbentity_id, r.pmid, r.title, r.volume, r.page, r.year, r.issue,"
				+	" r.citation, la.topic, j.med_abbr, rd.text, db.sgdid, db.date_created"
				+ " from nex.dbentity db"
				+ " inner join nex.referencedbentity r on db.dbentity_id = r.dbentity_id"
				+ " left join nex.journal j on j.journal_id = r.journal_id"
				+ " left join nex.referencedocument rd on  r.dbentity_id = rd.reference_id"
				+ " left join nex.literatureannotation la on la.reference_id = r.dbentity_id"
				+ " where rd.document_type = 'Abstract'"
				+ " order by r.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res; 
	}


	/**
	 * Return the results of running a query for all publications.  
	 * TODO only retreive publications for phenotype_annot_no
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getPubForPhenotype(Connection connection)
			throws SQLException {

		String query = "select annotation_id, reference_id, pmid, fulltext_status, rdb.title, volume, page, year, issue, citation, med_abbr, db.sgdid, db.date_created"
				+ " from nex.phenotypeannotation pa, nex.referencedbentity rdb, nex.journal j, nex.dbentity db"
				+ " where pa.reference_id = rdb.dbentity_id"
				+" and rdb.journal_id = j.journal_id"
				+ " and rdb.dbentity_id = db.dbentity_id"
				+ " order by annotation_id, reference_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res; 
	}
	/**
	 * Return the results of running a query for phenotype summaries
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getPhenotypeSummary(Connection connection)
			throws SQLException {

		String query = "select distinct(text), db.dbentity_id"
				+ " from nex.phenotypeannotation pa, nex.dbentity db,  nex.locussummary ls"
				+ " where pa.dbentity_id = db.dbentity_id"
				+ " and db.dbentity_id  = ls.locus_id"
				+ " and summary_type = 'Phenotype'"
				+ " group by text, db.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res; 
	}

	/**
	 * Return the results of running a query for phenotype summaries
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getFunctionSummary(Connection connection)
			throws SQLException {

		String query = "select distinct(text), db.dbentity_id"
				+ " from nex.goannotation ga, nex.dbentity db,  nex.locussummary ls"
				+ " where ga.dbentity_id = db.dbentity_id"
				+ " and db.dbentity_id  = ls.locus_id"
				+ " and summary_type = 'Function'"
				+ " group by text, db.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res; 
	}

	/**
	 * Return the results of running a query for phenotype summaries
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getGeneSummary(Connection connection)
			throws SQLException {

		String query = "select distinct(text), db.dbentity_id"
				+ " from nex.locusdbentity ldb, nex.dbentity db,  nex.locussummary ls"
				+ " where ldb.dbentity_id = db.dbentity_id"
				+ " and db.dbentity_id  = ls.locus_id"
				+ " and summary_type = 'Gene'"
				+ " group by text, db.dbentity_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res; 
	}


	/**
	 * Return the results of running a query for phenotype summaries
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getRegulationSummary(Connection connection)
			throws SQLException {

		String query = "select distinct(text), db.dbentity_id, reference_id, pmid"
				+ " from nex.locusdbentity ldb, nex.dbentity db,  nex.locussummary ls, nex.locussummary_reference lsr, nex.referencedbentity rdb"
				+ " where ldb.dbentity_id = db.dbentity_id"
				+ " and db.dbentity_id  = ls.locus_id"
				+ " and lsr.summary_id = ls.summary_id"
				+ " and lsr.reference_id = rdb.dbentity_id"
				+ " and summary_type = 'Regulation'"
				+ " group by text, db.dbentity_id, reference_id, pmid";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res; 
	}

	/**
	 * Return the results of running a query for CDSs and their sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getPhysicalInteractionResults(Connection connection)
			throws SQLException {
		
		String query = "select annotation_id, dbentity1_id, dbentity2_id, biogrid_experimental_system, bait_hit, s.display_name, annotation_type, psi.display_name as modification,"
				+ " citation, pmid, rdb.title, volume, page, year, issue, med_abbr, reference_id, substring(citation, 0, position( ')' in citation)+1) as first_author, db.sgdid, "
				+ " pa.description as note"
				+ " from nex.physinteractionannotation pa"
				+ " left join nex.psimod psi on pa.psimod_id = psi.psimod_id" 
				+ " inner join nex.source s on s.source_id = pa.source_id"
				+ " inner join  nex.referencedbentity rdb on pa.reference_id = rdb.dbentity_id"
				+ " left join nex.journal j on rdb.journal_id = j.journal_id"
				+ " inner join nex.dbentity db on  db.dbentity_id = rdb.dbentity_id";
				//+ " and (pa.dbentity1_id = 1268334 or pa.dbentity2_id = 1268334)";
		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of running a query for CDSs and their sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getGeneticInteractionResults(Connection connection)
			throws SQLException {

		String query = "select annotation_id, dbentity1_id, dbentity2_id, biogrid_experimental_system, p.display_name as phenotype, bait_hit, s.display_name as source, annotation_type,"
				+ " citation, pmid, rdb.title, volume, page, year, issue, med_abbr, reference_id, substring(citation, 0, position( ')' in citation)+1) as first_author, db.sgdid, "
				+ " split_part(ga.description, '|', 1) as note"
				+ " from nex.geninteractionannotation ga"
				+ " left join nex.phenotype p on p.phenotype_id = ga.phenotype_id"
				+ " inner join nex.source s on s.source_id = ga.source_id"
				+ " inner join nex.referencedbentity rdb on rdb.dbentity_id = ga.reference_id"
				+ " left join nex.journal j on rdb.journal_id = j.journal_id"
				+ " inner join nex.dbentity db on  db.dbentity_id = rdb.dbentity_id";
				//+ " and (ga.dbentity1_id = 1268334 or ga.dbentity2_id = 1268334)";
		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of running a query for CDSs and their sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getGeneticInteractionWithAllelesResults(Connection connection)
			throws SQLException {

		String query = " select annotation_id, interaction_id, allele1_id, allele2_id, sga_score, pvalue"
			+ " from nex.geninteractionannotation ga"
			+ " inner join nex.allele_geninteraction ag  on  ag.interaction_id = ga.annotation_id"
		    //+ " and (ga.dbentity1_id = 1268334 or ga.dbentity2_id = 1268334)"
			+ " order by ga.annotation_id ";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}



	/**
	 * Return the results of running a query for CDSs and their sequences
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getPhenotypeResults(Connection connection)
			throws SQLException {

		String query = "  select db.dbentity_id, pa.annotation_id, pac.group_id, p.display_name as phenotype, t.format_name as strain_name, pa.details,"
				+ " pa.experiment_comment, pa.allele_comment, pa.reporter_comment, a1.display_name as experiment, a2.display_name as mutant,"
				+ " al.dbentity_id as allele, rp.display_name as reporter, o.display_name as assay, rdb.pmid, rdb.dbentity_id as refNo,"
				+ " array_agg(condition_class) as condclass,"
				+ " array_agg(condition_name) as condname , array_agg(condition_value) as condvalue,"
				+ " array_agg(condition_unit) as condunit"
				+ " from nex.phenotypeannotation pa"
				+ " inner join nex.dbentity db on  db.dbentity_id = pa.dbentity_id"
				+ " inner join nex.phenotype p on pa.phenotype_id = p.phenotype_id"
				+ " inner join nex.referencedbentity rdb on pa.reference_id = rdb.dbentity_id"
				+ " inner join nex.taxonomy t on  pa.taxonomy_id = t.taxonomy_id"
				+ " left join nex.phenotypeannotation_cond pac on pac.annotation_id = pa.annotation_id"
				+ " left join nex.apo a1 on pa.experiment_id = a1.apo_id"
				+ " left join nex.apo a2 on pa.mutant_id = a2.apo_id"
				+ " left join nex.alleledbentity al on al.dbentity_id = pa.allele_id"
				+ " left join nex.reporter rp on rp.reporter_id = pa.reporter_id"
				+ " left join nex.obi o on pa.assay_id = o.obi_id"
				+ " group by  db.dbentity_id, pa.annotation_id, pac.group_id, p.display_name, t.format_name, pa.details,"
				+ " pa.experiment_comment, pa.allele_comment, pa.reporter_comment, a1.display_name, a2.display_name,"
				+ " al.dbentity_id, rp.display_name, o.display_name, rdb.pmid, rdb.dbentity_id"
				+ " order by db.dbentity_id, pa.annotation_id, pac.group_id ";
	
		LOG.info("executing: " + query);        
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of getting aliases
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getAliases(Connection connection)
			throws SQLException {

		String query = "select  l.dbentity_id, la.display_name, alias_type "
				+ "from nex.locusdbentity l, nex.locus_alias la, nex.contig c, nex.dnasequenceannotation a, nex.so s, nex.dbentity d "
				+ "where l.dbentity_id = a.dbentity_id "
				+ "and l.dbentity_id = la.locus_id "
				+ "and c.contig_id = a.contig_id "
				+ "and c.so_id = s.so_id "
				+ "and s.display_name in ('chromosome', 'plasmid')"
				+ "and a.dna_type = 'GENOMIC' "
				+ "and a.taxonomy_id = 274901 "
				+ "and d.dbentity_id = l.dbentity_id "
				+ "and d.dbentity_status = 'Active' "
				+ "and alias_type in ('Uniform', 'Non-uniform', 'Retired name', 'NCBI protein name')";


		LOG.info("executing: " + query);        
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of getting paralogs
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getParalogs(Connection connection)
			throws SQLException {

		String query = "SELECT parent_id, child_id, reference_id "
				+ "from nex.locus_relation lr, nex.locusrelation_reference lrr "
				+ " where ro_id = 169738 "
				+ " and lr.relation_id = lrr.relation_id";

		LOG.info("executing: " + query);        
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of getting cross-references from dbxref
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getCrossReferences(Connection connection)
			throws SQLException {

		String query = "select l.dbentity_id, la.display_name, so.display_name, alias_type "
				+ "from nex.locusdbentity l, nex.locus_alias la, nex.contig c, nex.dnasequenceannotation a, nex.so s, nex.dbentity d, nex.source so "
				+ "where l.dbentity_id = a.dbentity_id "
				+ "and l.dbentity_id = la.locus_id "
				+ "and c.contig_id = a.contig_id "
				+ "and c.so_id = s.so_id "
				+ "and so.source_id = la.source_id "
				+ "and s.display_name in ('chromosome', 'plasmid') "
				+ "and a.dna_type = 'GENOMIC' "
				+ "and a.taxonomy_id = 274901 "
				+ "and d.dbentity_id = l.dbentity_id "
				+ "and d.dbentity_status = 'Active' "
				+ "and alias_type NOT in ('Uniform', 'Non-uniform', 'Retired name', 'NCBI protein name')";

		LOG.info("executing: " + query);        
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}

	/**
	 * Return the results of getting cross-references from dbxref
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getGenePathways(Connection connection)
			throws SQLException {

		String query = "select distinct db.dbentity_id, biocyc_id, db2.display_name, pa.reference_id" +
				" from nex.dbentity db" +
				" inner join nex.pathwayannotation pa on db.dbentity_id = pa.dbentity_id" +
				" inner join nex.pathwaydbentity pdf on pdf.dbentity_id = pa.pathway_id" +
				" inner join nex.dbentity db2 on db2.dbentity_id = pdf.dbentity_id" +
				" order by 1";

		LOG.info("executing: " + query);        
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results of getting cross-references from dbxref
	 * @param connection the connection
	 * @return the results
	 * @throws SQLException if there is a database problem
	 */
	protected ResultSet getAllPathways(Connection connection)
			throws SQLException {

		String query = "select biocyc_id, db.display_name, ps.summary_type, ps.text, pss.reference_id" +
				" from nex.pathwaydbentity pdb " +
				" inner join nex.dbentity db on db.dbentity_id = pdb.dbentity_id" +
				" inner join nex.pathwaysummary ps on ps.pathway_id = pdb.dbentity_id" +
				" left join nex.pathwaysummary_reference pss on pss.summary_id = ps.summary_id";

		LOG.info("executing: " + query);
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


	/**
	 * Return the results for regulation info
	 * @param connection the connectio
	 * @return the results
	 * @throws SQLException if there is a database proble
	 */
	protected ResultSet getRegulationData(Connection connection)
			throws SQLException {

		String query = "select target_id, regulator_id, regulator_type, regulation_type, direction, annotation_type, pmid,"
				+ " rdb.dbentity_id as refNo, e.format_name, g.display_name as happens_during, s.display_name as source, t.format_name as strain_background"
				+ " from nex.dbentity db"
				+ " inner join nex.regulationannotation ra on db.dbentity_id = ra.regulator_id"
				+ " inner join nex.referencedbentity rdb on rdb.dbentity_id = ra.reference_id"
				+ " inner join nex.dbentity db2 on ra.target_id = db2.dbentity_id"
				+ " left join nex.eco e on e.eco_id = ra.eco_id"
				+ " left join nex.go g on g.go_id = ra.happens_during"
				+ " left join nex.taxonomy t on t.taxonomy_id = ra.taxonomy_id"
				+ " left join nex.source s on s.source_id = ra.source_id"
				//+ " where db.dbentity_id = 1267652"
				+ " order by 1";

		LOG.info("executing: " + query);        
		Statement stmt = connection.createStatement();
		ResultSet res = stmt.executeQuery(query);
		return res;
	}


}
