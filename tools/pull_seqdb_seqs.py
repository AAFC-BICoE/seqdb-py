'''
Created on Feb 12, 2015

@author: korolo

Given a SeqDB URL and API key, uses SeqDB webservices to load sequences, filtered by parameters 
below, into a fasta file.

Usage (for most current usage run this script with no options: >python pull_seqdb_seqs.py):
pull_seqdb_seqs -c <Path to yaml config with SeqDB API info> 
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
with one of the following options:
 --its
 --consensus
 --raw
Other arguments:
   -h   help (prints this message)
'''

import argparse
import logging.config
import os
import requests.exceptions
import sys 

from api.BaseSeqdbApi import UnexpectedContent
from api.BaseSequenceEntity import BaseSequenceEntity
from api.ConsensusSequenceApi import ConsensusSequenceApi
from api.RawSequenceApi import RawSequenceApi
from config import config_root
import tools_helper


### Values below are used in Galaxy wrappers, so make sure you know what you're doing if you're changing any of them. 
# File name where the pulled sequences will be stored. 
output_file_name = "seqdb_sequences."

# File name where taxonomy for the sequences will be stored. Optional output. 
output_taxonomy_file_name = "seqdb_taxonomy_file.txt"

# Values for the types of sequences this script downloads. I.e. "its" loads 
# ITS sequences. Note that raw sequences are not implemented in SeqDB yet. "raw":"raw",
pull_types_dict = {"its":"its", "consensus":"consensus", "raw":"raw", "all":"all"}

return_types_dict = {"fasta":"fasta", "fastq":"fastq"}

# Taxonomy ranks that can be specified as a filter parameter with corresponding value.
# These are used as a drop-down value in the wrapper.
taxonomy_ranks = {"species", "genus", "family", "order", "class", "phylum"}


def set_up_logging():
    ''' Loads main configuration file and sets up logging for the script '''
    main_conf = tools_helper.load_config(config_root.path() + '/config.yaml')

    if not main_conf:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)
    
    logging.config.dictConfig(main_conf['logging'])
    
    logging.info("%s %s" % (tools_helper.log_msg_scriptExecutionWithParams, sys.argv))
    

def parse_input_args(argv):
    ''' Parses command line arguments 
    Args:
        argv: command line arguments 
    Returns: 
        seqdb api_key and base_url to use for web services requests.
    '''
    pull_types_set = frozenset(pull_types_dict.values())
    return_types_set = frozenset(return_types_dict.values())
    
    parser = argparse.ArgumentParser(description="Load sequences from SeqDB")
    parser.add_argument('seq_type', help="Type of sequences to load", type=str, choices=pull_types_set)
    parser.add_argument('-c', help="SeqDB config file", dest="config_file", required=False)
    parser.add_argument('-u', help="SeqDB API URL", dest="base_url", required=False)
    parser.add_argument('-k', help="SeqDB API key", dest="api_key", required=False)
    parser.add_argument('-r', help="Return file type: fasta (default) or fastq", dest="return_type", choices=return_types_set, required=False)    
    parser.add_argument('-t', help="Output taxonomy file as well as sequence file", dest="output_taxonomy_file", action='store_true', required=False)
    parser.add_argument('--specNums', help="Specimen number(s). If multiple, separate by comma.", dest="specimen_nums", required=False)    
    parser.add_argument('--seqName', help="Sequence name (keyword)", dest="sequence_name", required=False)
    parser.add_argument('--sampleName', help="Sample name (keyword)", dest="sample_name", required=False)   
    parser.add_argument('--geneRegion', help="Gene region name (keyword)", dest="gene_region_name", required=False)    
    parser.add_argument('--projectName', help="Project Name (keyword)", dest="project_name", required=False)    
    parser.add_argument('--collectionCode', help="Collection code (keyword)", dest="collection_code", required=False)    
    parser.add_argument('--pubRefSeqs', help="Public reference sequences", dest="pub_ref_seqs", action='store_true', required=False)    
    parser.add_argument('--taxRank', help="Taxonomic rank to filter sequences on (need to specify the value as well). Ex. --taxRank phylum --taxValue chordata ", dest="tax_rank", choices=taxonomy_ranks, required=False)
    parser.add_argument('--taxValue', help="Value for the taxonomic rank to be filtered on (need to specify the rank as well)", dest="tax_value", required=False)
    #parser.add_argument('-t', help="Type of sequences to load", dest="load_type", type=str, choices=set(("its","consensus")), required=True)
    
    args = parser.parse_args(argv)
    
    if not (args.return_type):
        args.return_type = "fasta"
    elif args.return_type not in return_types_set:
        parser.error('Return type (-r) should be one of the following: {}'.format(return_types_set))
    elif args.return_type == "fastq" and args.seq_type != pull_types_dict["raw"]:
        parser.error('Fastq file format is only possible for raw sequences.')
           
    if not (args.config_file or (args.base_url and args.api_key)):
        parser.error('Either -c <configuration file>, or -u <base_url> -k <api_key> have to be specified')
    
    if args.seq_type == pull_types_dict["its"] and (args.specimen_nums or args.sequence_name or args.pub_ref_seqs):
        parser.error('ITS sequences can not be restricted by filters at the moment. Please do not use --specNum, --seqName or any other additional options.')
    
    if bool(args.tax_rank) != bool(args.tax_value):
        parser.error('Either both --taxRank and --taxValue have to be specified, or none.')
         
    return args


def get_ITS_seq_ids(rawSequence):
    ''' Get all sequence ids, which are associated with the ITS regions
    Args:
        rawSequence: Raw Sequence Entity
    Returns:
        a list of ITS sequence IDs
    '''
    
    #set_up_logging()

    ### Get sequence IDs for the ITS regions
    its_seq_ids = set()
    its_region_names = ["18s", "its", "28s",  "ssu", "16s", "5.8s", "lsu", "23s", "25s", "internal transcribed spacer"]
    for its_region_keyword in its_region_names:
        #TODO: parallelize; use locking when appending to its_seq_ids
        try:
            
            BaseSequenceEntity.regionNameFilter = its_region_keyword
            curr_seq_ids = rawSequence.getIds()
            its_seq_ids.update(curr_seq_ids)
                
        except requests.exceptions.ConnectionError as e:
            logging.error(tools_helper.log_msg_noDbConnection)
            logging.error(e.message)
            if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                logging.error(tools_helper.log_msg_seqDBConnection)
            if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                logging.error(tools_helper.log_msg_invalidURL)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.ReadTimeout as e:
            logging.error(tools_helper.log_msg_slowConnection)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.HTTPError as e:
            logging.error(tools_helper.log_msg_httpError)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except UnexpectedContent as e:
            logging.error(tools_helper.log_msg_apiResponseFormat)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except Exception as e:
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
    
    msg_numITSseqs = "Number of ITS sequences retrieved:"
    logging.info("%s %i " % (msg_numITSseqs, len(its_seq_ids)))

    return list(its_seq_ids)

"""
def get_seq_ids(rawSequence, consensusSequence, pull_type, specimen_nums=None, sequence_name=None, sample_name=None, pub_ref_seqs=None, region_name=None, project_name=None, collection_code=None, taxonomy_rank=None, taxonomy_value=None):
    ''' Gets sequence ids based on specified parameters 
    Agrs:
        pull_type: string of pre-determined values. Values should correspond to the values of pull_types_dict
        specimen_nums: if specified, list of specimen numbers for which the sequence ids will be retrieved
        sequence_name: if specified, sequence name for which the sequence ids will be retrieved
        sample_name: if specified, sample name for which the sequence ids will be retrieved
        pub_ref_seqs: if specified, pub ref sequence for which the sequence ids will be retrieved
        region_name: if specified, region name for which the sequence ids will be retrieved
        project_name: if specified, project name for which the sequence ids will be retrieved
        collection_code: if specified, collection_code for which the sequence ids will be retrieved
        taxonomy_rank: if specified, taxonomy rank for which the sequence ids will be retrieved
        setaxonomy_value: if specified, taxonomy value for which the sequence ids will be retrieved
    Returns:
        a list of sequences IDs based on the filter that has been chosen
    '''

    BaseSequenceEntity.sequenceNameFilter = sequence_name
    BaseSequenceEntity.sampleNameFilter = sample_name
    BaseSequenceEntity.pubRefSeqFilter = pub_ref_seqs
    BaseSequenceEntity.regionNameFilter = region_name
    BaseSequenceEntity.projectNameFilter = project_name
    BaseSequenceEntity.collectionCodeFilter = collection_code
    BaseSequenceEntity.taxonomyRankFilter = taxonomy_rank
    BaseSequenceEntity.taxonomyValueFilter = taxonomy_value
        
    if pull_type not in pull_types_dict.values():
        msg = "Value for pull_type should be one of the following: %s" %pull_types_dict.values()
        logging.error(msg)
        sys.exit(tools_helper.log_msg_sysExit + msg)
        
    seq_ids = []
    
    if pull_type == pull_types_dict["its"]:
        seq_ids = get_ITS_seq_ids(rawSequence)  
    
    else:
        try:
             
            if pull_type == pull_types_dict["consensus"]:
                if not specimen_nums:
                    specimen_nums = [None]
                for specimen_num in specimen_nums:
                    BaseSequenceEntity.specimenNumFilter = specimen_num
                    curr_seq_ids = consensusSequence.getIds()
                    seq_ids.extend(curr_seq_ids)       
                log_msg = "Number of consensus sequences retrieved:"
            
            elif pull_type == pull_types_dict["all"]:
                if not specimen_nums:
                    specimen_nums=[None]
                for specimen_num in specimen_nums:
                    BaseSequenceEntity.specimenNumFilter = specimen_num
                    curr_seq_ids_raw = rawSequence.getIds()
                    curr_seq_ids_consensus = consensusSequence.getIds()
                    seq_ids = curr_seq_ids_raw + curr_seq_ids_consensus      
                log_msg = "Number of all sequences retrieved:"
            
            elif pull_type == pull_types_dict["raw"]:
                if not specimen_nums:
                    specimen_nums = [None]
                for specimen_num in specimen_nums:
                    BaseSequenceEntity.specimenNumFilter = specimen_num
                    curr_seq_ids = rawSequence.getIds()
                    seq_ids.extend(curr_seq_ids)
                log_msg = "Number of raw sequences retrieved:"
        
        except requests.exceptions.ConnectionError as e:
            logging.error(tools_helper.log_msg_noDbConnection)
            logging.error(e.message)
            if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                logging.error(tools_helper.log_msg_seqDBConnection)
            if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                logging.error(tools_helper.log_msg_invalidURL)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.ReadTimeout as e:
            logging.error(tools_helper.log_msg_slowConnection)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.HTTPError as e:
            logging.error(tools_helper.log_msg_httpError)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except UnexpectedContent as e:
            logging.error(tools_helper.log_msg_apiResponseFormat)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except Exception as e:
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        
        logging.info("%s %i " % (log_msg, len(seq_ids)))
        
    return seq_ids
"""

def write_to_file(result_fasta, file_name, file_type):
    ''' Writes the fasta sequences to a file
    Args:
        result_fasta: fasta sequences retrieved
        file_name: name of output file
        file_type: fasta or fastq file types
    '''
    
    with open(file_name + file_type, 'a+r') as output_file:   
        try:
            output_file.write(result_fasta)
        except requests.exceptions.ConnectionError as e:
            logging.error(tools_helper.log_msg_noDbConnection)
            logging.error(e.message)
            if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                logging.error(tools_helper.log_msg_seqDBConnection)
            if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                logging.error(tools_helper.log_msg_invalidURL)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.ReadTimeout as e:
            logging.error(tools_helper.log_msg_slowConnection)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.HTTPError as e:
            logging.error(tools_helper.log_msg_httpError)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except UnexpectedContent as e:
            logging.error(tools_helper.log_msg_apiResponseFormat)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except Exception as e:
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
     
    msg_fileName = "Sequences written to a file:"
    logging.info("%s %s" % (msg_fileName, os.path.abspath(output_file.name)))

    #TODO: find a way to return the number of sequences that have been written to file

    #msg_seqNum = "Number of sequences written:"
    #logging.info("%s %s" % (msg_seqNum, len(idList)) )
  
            
def write_sequence_file(rawSequence, its_seq_ids, file_name, file_type):
    ''' Gets fasta or fastq sequences for ITS only based on IDs and writes to a file 
    Args:
        rawSequence: Raw Sequence Entity
        its_seq_ids: sequence IDs
        file_name: name of output file
        file_type: fasta or fastq file types
    Returns
        a list of IDs that have been successfully written to file
    '''

    output_file = open(file_name + file_type, 'w')
    
    success_ids = []
    for seq_id in its_seq_ids:
        #TODO: threads?
        '''
        import threading
        
        t = threading.Thread(target=<this try code>, args(seq_id,)
        '''
        try:
            if file_type == 'fasta':
                sequence = rawSequence.getFastaSequence(seq_id)    
            
            elif file_type == 'fastq':
                sequence = rawSequence.getFastqSequence(seq_id)
            
            output_file.write(sequence)
            success_ids.append(seq_id)
            
        except requests.exceptions.ConnectionError as e:
            logging.error(tools_helper.log_msg_noDbConnection)
            logging.error(e.message)
            if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                logging.error(tools_helper.log_msg_seqDBConnection)
            if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                logging.error(tools_helper.log_msg_invalidURL)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.ReadTimeout as e:
            logging.error(tools_helper.log_msg_slowConnection)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.HTTPError as e:
            logging.error(tools_helper.log_msg_httpError)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except UnexpectedContent as e:
            logging.error(tools_helper.log_msg_apiResponseFormat)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except Exception as e:
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
     
    output_file.close()   

    msg_fileName = "Sequences written to a file:"
    logging.info("%s %s" % (msg_fileName, os.path.abspath(output_file.name)))

    msg_seqNum = "Number of sequences written:"
    logging.info("%s %s" % (msg_seqNum, len(success_ids)) )

    return success_ids
    
    
def write_taxonomy_file(rawSequence, seq_ids, output_file_name):
    ''' Gets fasta sequences based on IDs and writes to a file
    Args:
        rawSequence: Raw Sequence Entity
        seq_ids: a list of sequence IDs
        output_file_name: name of output file
    Returns:
        a list of IDs that have been successfully written to file
    '''     
    with open(output_file_name, 'w') as output_file:
        success_ids = []
        for seq_id in seq_ids:
            try:
                # Quiime/Mothur taxonomy file format as described in Unite:
                # https://unite.ut.ee/repository.php
                unclassified_keyword = u'unclassified'
                determ_jsn = rawSequence.getAcceptedSpecimenDetermination(seq_id)
                if determ_jsn:
                    t_kingdom = determ_jsn["taxonomy"]["kingdom"]
                    if t_kingdom:
                        t_kingdom.replace(" ", "_")
                    else:
                        t_kingdom = unclassified_keyword
                    
                    t_phylum = determ_jsn["taxonomy"]["phylum"]
                    if t_phylum:
                        t_phylum.replace(" ", "_")  
                    else:
                        t_phylum = unclassified_keyword
                    
                    t_class = determ_jsn["taxonomy"]["taxanomicClass"]
                    if t_class:
                        t_class.replace(" ", "_")
                    else: 
                        t_class = unclassified_keyword
                    
                    t_order = determ_jsn["taxonomy"]["taxanomicOrder"]
                    if t_order:
                        t_order.replace(" ", "_")
                    else:
                        t_order = unclassified_keyword
                    
                    t_family = determ_jsn["taxonomy"]["family"]
                    if t_family:
                        t_family.replace(" ", "_")
                    else:
                        t_family = unclassified_keyword
                    
                    t_genus = determ_jsn["taxonomy"]["genus"]
                    if t_genus:
                        t_genus.replace(" ", "_")
                    else:
                        t_genus = unclassified_keyword
                    
                    t_species = determ_jsn["taxonomy"]["species"]
                    if t_species or determ_jsn["taxonomy"]["genus"]:
                        if t_species:
                                t_species.replace(" ", "_")
                        else:
                            t_species = "sp."
                        t_species = "{}_{}".format(t_genus, t_species)      
                    else:
                        t_species = unclassified_keyword
                    
                    taxonomy_line = u'{}\tk__{};p__{};c__{};o__{};f__{};g__{};s__{}\n'.format(seq_id, 
                                                                                            t_kingdom,
                                                                                            t_phylum,
                                                                                            t_class,
                                                                                            t_order,
                                                                                            t_family,
                                                                                            t_genus,
                                                                                            t_species)
                    output_file.write(taxonomy_line)
                    success_ids.append(seq_id)
                    
            except requests.exceptions.ConnectionError as e:
                logging.error(tools_helper.log_msg_noDbConnection)
                logging.error(e.message)
                if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                    logging.error(tools_helper.log_msg_seqDBConnection)
                if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                    logging.error(tools_helper.log_msg_invalidURL)
                sys.exit(tools_helper.log_msg_sysExit)
            except requests.exceptions.ReadTimeout as e:
                logging.error(tools_helper.log_msg_slowConnection)
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
            except requests.exceptions.HTTPError as e:
                logging.error(tools_helper.log_msg_httpError)
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
            except UnexpectedContent as e:
                logging.error(tools_helper.log_msg_apiResponseFormat)
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
            except Exception as e:
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
     
    msg_fileName = "Taxonomy written to a file:"
    logging.info("%s %s" % (msg_fileName, os.path.abspath(output_file.name)))
    
    return success_ids

    
def execute_script(input_args, output_file_name=output_file_name, output_taxonomy_file_name=output_taxonomy_file_name):
    
    ### Load main configuration file and set up logging for the script
    set_up_logging()
    
    ### Parse sript's input arguments
    parsed_args = parse_input_args(input_args)
    
    if parsed_args.config_file:
        tool_config = tools_helper.load_config(parsed_args.config_file)
        base_url = tool_config['seqdb']['base_url'] 
        api_key = tool_config['seqdb']['api_key'] 
    else:
        base_url = parsed_args.base_url 
        api_key = parsed_args.api_key 
        
    logging.info("%s '%s'" % (tools_helper.log_msg_apiUrl, base_url))
    
    ### Script execution
    seq_ids = []
    
    if pull_types_dict["its"] == parsed_args.seq_type:
        logging.info(tools_helper.log_msg_ITSLoad)
        rawSequence = RawSequenceApi(api_key, base_url)
        seq_ids = get_ITS_seq_ids(rawSequence)
        success_seq_ids = write_sequence_file(rawSequence, seq_ids, output_file_name, parsed_args.return_type)
        print("Number of sequences retrieved from Sequence Database:  {}".format(len(success_seq_ids))) 

    #TODO: log number of consensus sequences retrieved 
    elif pull_types_dict["consensus"] == parsed_args.seq_type:
        consensusSequence = ConsensusSequenceApi(api_key, base_url)
        
        #TODO: Test when filtering for a list of specimen ids
        specimen_nums_list = None
        if not parsed_args.specimen_nums:
            parsed_args.specimen_nums = [None]
        elif parsed_args.specimen_nums:
            specimen_nums_list = parsed_args.specimen_nums.replace(" ", "").split(",")
            for specimen_num in specimen_nums_list:
                consensusSequence.specimenNumFilter = specimen_num
                    
        consensusSequence.sequenceNameFilter = parsed_args.sequence_name
        consensusSequence.sampleNameFilter = parsed_args.sample_name
        consensusSequence.pubRefSeqFilter = parsed_args.pub_ref_seqs
        consensusSequence.regionNameFilter = parsed_args.gene_region_name
        consensusSequence.projectNameFilter = parsed_args.project_name
        consensusSequence.collectionCodeFilter = parsed_args.collection_code
        consensusSequence.taxonomyRankFilter = parsed_args.tax_rank
        consensusSequence.taxonomyValueFilter = parsed_args.tax_value    

        if return_types_dict["fasta"] == parsed_args.return_type: 
            result_fasta, result_offset = consensusSequence.getFastaSequencesWithOffset(offset=0)
            write_to_file(result_fasta, output_file_name, parsed_args.return_type)
            while result_offset > 0:
                result_fasta, result_offset = consensusSequence.getFastaSequencesWithOffset(offset=result_offset)
                write_to_file(result_fasta, output_file_name, parsed_args.return_type)
        if (parsed_args.output_taxonomy_file):
            seq_ids = consensusSequence.getIds()
            
    #TODO: log number of raw sequences retrieved         
    elif pull_types_dict["raw"] == parsed_args.seq_type:
        rawSequence = RawSequenceApi(api_key, base_url)
            
        #TODO: Test when filtering for a list of specimen ids
        specimen_nums_list = None
        if not parsed_args.specimen_nums:
            parsed_args.specimen_nums = [None]
        elif parsed_args.specimen_nums:
            specimen_nums_list = parsed_args.specimen_nums.replace(" ", "").split(",")
            for specimen_num in specimen_nums_list:
                rawSequence.specimenNumFilter = specimen_num
                    
        rawSequence.sequenceNameFilter = parsed_args.sequence_name
        rawSequence.sampleNameFilter = parsed_args.sample_name
        rawSequence.pubRefSeqFilter = parsed_args.pub_ref_seqs
        rawSequence.regionNameFilter = parsed_args.gene_region_name
        rawSequence.projectNameFilter = parsed_args.project_name
        rawSequence.collectionCodeFilter = parsed_args.collection_code
        rawSequence.taxonomyRankFilter = parsed_args.tax_rank
        rawSequence.taxonomyValueFilter = parsed_args.tax_value    
        
        if return_types_dict["fasta"] == parsed_args.return_type:
            result_fasta, result_offset = rawSequence.getFastaSequencesWithOffset(offset=0)
            write_to_file(result_fasta, output_file_name, parsed_args.return_type)
            while result_offset > 0:
                result_fasta, result_offset = rawSequence.getFastaSequencesWithOffset(offset=result_offset)
                write_to_file(result_fasta, output_file_name, parsed_args.return_type)
        
        if return_types_dict["fastq"] == parsed_args.return_type:
            result_fastq, result_offset = rawSequence.getFastqSequencesWithOffset(offset=0, limit=20)
            write_to_file(result_fastq, output_file_name, parsed_args.return_type)
            while result_offset > 0:
                result_fastq, result_offset = rawSequence.getFastqSequencesWithOffset(offset=result_offset, limit=20)
                write_to_file(result_fastq, output_file_name, parsed_args.return_type)
        
        if (parsed_args.output_taxonomy_file):
            seq_ids = rawSequence.getIds()

    #TODO: log number of all sequences retrieved                 
    elif pull_types_dict["all"] == parsed_args.seq_type:
        consensusSequence = ConsensusSequenceApi(api_key, base_url)
        rawSequence = RawSequenceApi(api_key, base_url)
        
        #TODO: Test when filtering for a list of specimen ids
        specimen_nums_list = None
        if not parsed_args.specimen_nums:
            parsed_args.specimen_nums = [None]
        elif parsed_args.specimen_nums:
            specimen_nums_list = parsed_args.specimen_nums.replace(" ", "").split(",")
            for specimen_num in specimen_nums_list:
                consensusSequence.specimenNumFilter = specimen_num
                    
        consensusSequence.sequenceNameFilter = parsed_args.sequence_name
        consensusSequence.sampleNameFilter = parsed_args.sample_name
        consensusSequence.pubRefSeqFilter = parsed_args.pub_ref_seqs
        consensusSequence.regionNameFilter = parsed_args.gene_region_name
        consensusSequence.projectNameFilter = parsed_args.project_name
        consensusSequence.collectionCodeFilter = parsed_args.collection_code
        consensusSequence.taxonomyRankFilter = parsed_args.tax_rank
        consensusSequence.taxonomyValueFilter = parsed_args.tax_value    

        #TODO: Test when filtering for a list of specimen ids
        specimen_nums_list = None
        if not parsed_args.specimen_nums:
            parsed_args.specimen_nums = [None]
        elif parsed_args.specimen_nums:
            specimen_nums_list = parsed_args.specimen_nums.replace(" ", "").split(",")
            for specimen_num in specimen_nums_list:
                rawSequence.specimenNumFilter = specimen_num
                    
        rawSequence.sequenceNameFilter = parsed_args.sequence_name
        rawSequence.sampleNameFilter = parsed_args.sample_name
        rawSequence.pubRefSeqFilter = parsed_args.pub_ref_seqs
        rawSequence.regionNameFilter = parsed_args.gene_region_name
        rawSequence.projectNameFilter = parsed_args.project_name
        rawSequence.collectionCodeFilter = parsed_args.collection_code
        rawSequence.taxonomyRankFilter = parsed_args.tax_rank
        rawSequence.taxonomyValueFilter = parsed_args.tax_value    

        if return_types_dict["fasta"] == parsed_args.return_type:
            result_fasta_consensus, result_offset_consensus = consensusSequence.getFastaSequencesWithOffset(offset=0)
            write_to_file(result_fasta_consensus, output_file_name, parsed_args.return_type)
            while result_offset_consensus > 0:
                result_fasta_consensus, result_offset_consensus = consensusSequence.getFastaSequencesWithOffset(offset=result_offset_consensus, limit=20)
                write_to_file(result_fasta_consensus, output_file_name, parsed_args.return_type)
                
            result_fasta_raw, result_offset_raw = rawSequence.getFastaSequencesWithOffset(offset=0, limit=20)
            write_to_file(result_fasta_raw, output_file_name, parsed_args.return_type)
            while result_offset_raw > 0:
                result_fasta_raw, result_offset_raw = rawSequence.getFastaSequencesWithOffset(offset=result_offset_raw, limit=20)
                write_to_file(result_fasta_raw, output_file_name, parsed_args.return_type)
                
        if (parsed_args.output_taxonomy_file):
            seq_ids = consensusSequence.getIds().append(rawSequence.getIds())
          
    if (parsed_args.output_taxonomy_file) and seq_ids:
        rawSequence = RawSequenceApi(api_key, base_url)
        write_taxonomy_file(rawSequence, seq_ids, output_taxonomy_file_name)
        print("Taxonomy file is written to a file: '%s'" % output_taxonomy_file_name)

    ### Post-execution: messages and logging
    print("Sequences are written to a file: {}".format(output_file_name + parsed_args.return_type))
    print("Execution complete.")
    
    logging.info(tools_helper.log_msg_execEnded)


def main():
    ''' This method has to have no arguments to create entry point to the egg.
        execute_script method was extracted so that we can create the unit 
        tests for the whole script execution.
    '''
    execute_script(sys.argv[1:])


if __name__ == '__main__':
    main()