'''
Created on Jul 10, 2015

@author: korolo
'''
import time

class TaxonomyLineage(object):
    '''
    Sets up an in-memory NCBI taxonomy db to be queried for taxonomy lineage
    '''


    def __init__(self, taxonomy_location):
        '''
        Reads two taxonomy files: nodes.dmp and names.dmp in memory
        Initializes two datastructures:
            self.taxonomic_lineage_ids: dictionary of tax_id:[parent_tax_id, taxonomic_rank]. 
                Generated from nodes.dmp file
            self.taxonomy_names: dictionary of tax_id: taxonomic_name. 
                Generated from names.dmp. Note that for multiple records with the same tax_id,
                only scientific names are kept.
        '''
        start_time = time.time()
        # Read in (tax_id, parent tax_id, tax_rank) from the nodes.dmp file
        nodes_file_path = taxonomy_location + "nodes.dmp"
        self.taxonomic_lineage_ids = {}
        with open(nodes_file_path, "r") as nodes_file:
            for line in nodes_file:
                line_tokens = line.split('\t')
                tax_id = int(line_tokens[0])
                parent_tax_id = int(line_tokens[2])
                taxon_rank = line_tokens[4]
                self.taxonomic_lineage_ids[tax_id] = [parent_tax_id,taxon_rank]
                #if tax_id == 129355:
                #    print "Got 129355"
        
        print "Time for loading nodes.dmp: %s" %(time.time()-start_time)
        
        start_time = time.time()
        # Read in tax_id - Taxonomy name pairs from the names.dmp list
        # Note that the tax_id column (first column) in names file does not contain unique ids.
        #      But we only keep a one-to-one relationship of tax_id - tax_name, therefore for 
        #      multiple records, we only keep those records that say they are "scientific name"
        names_file_path = taxonomy_location + "names.dmp"
        self.taxonomy_names = {}
        
        # assume that names.dmp is sorted on tax_id: otherwise takes too much time for reverse lookup
        prev_tax_id = -1
        with open(names_file_path, "r") as names_file:
            for line in names_file:
                line_tokens = line.split('\t')
                tax_id = int(line_tokens[0])
                tax_name = line_tokens[2]
                identif_src = line_tokens[6]
                if identif_src != "misspelling":
                    if (tax_id == prev_tax_id and identif_src == "scientific name") \
                        or tax_id != prev_tax_id:
                        self.taxonomy_names[tax_id] = tax_name
                
                prev_tax_id = tax_id 
                      
        print "Time for loading names.dmp: %s" %(time.time()-start_time)


    def findLineage(self, tax_id):
        ''' Get taxonomic lineage for the provided taxonomy id
        Return:
            list of tuples [taxonomic_rank, taxonomic_name] in the order from the given tax_id up.
        '''
        lineage_ids = list()
        lineage_ranks = list()
        curr_id = tax_id
        while curr_id != 1:
            curr_rank = self.taxonomic_lineage_ids[curr_id][1]
            lineage_ids.append(curr_id)
            lineage_ranks.append(curr_rank)
            curr_id = self.taxonomic_lineage_ids[curr_id][0]
            
        lineage_names = list()
        for tax_id, rank in zip(lineage_ids,lineage_ranks):
            lineage_names.append([rank, self.taxonomy_names[tax_id]])
            
        return lineage_names  