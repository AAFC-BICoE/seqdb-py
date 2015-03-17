'''
Created on Feb 12, 2015

@author: korolo

Usage:
python pull_seqdb_its_seqs.py -k <SeqDB API key>
or
python pull_seqdb_its_seqs.py --seqdb_api_key <SeqDB API key>
-h is help
-p is to use production URL (uses ***REMOVED*** (Oksana's) machine by default, since we don't have UAT for SeqDB as of this time)
'''
import sys, getopt, requests
from api.seqdbWebService import seqdbWebService

usage_help_line = """Usage of the script: \npull_seqdb_its_seqs.py -k <SeqDB API key>
Other arguments:
   -h   help (prints this message)
   -p   use production url for web service requests
   -u <url> specify base url for web services requests"""

prod_url = "***REMOVED***/api/v1"
local_url = "***REMOVED***:8181/seqdb.web-2.5/api/v1"
# file name where the received sequences will be stores
output_file_name = 'seqdb_ITS_seqs.fasta'



# Parses command line arguments 
# Returns seqdb api_key and base url to use for web services requests
def parse_input_args(argv):
    seqdb_api_key = ''
    base_url=''
    prod=False
    user_url=''
    
    try:
        opts, args = getopt.getopt(argv,"hpk:u:",["seqdb_api_key=", "seqdb_ws_url="])
    except getopt.GetoptError:
        print usage_help_line
        sys.exit(2)
        
    if len(opts)==0:
        print usage_help_line
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print usage_help_line
            sys.exit()
        elif opt in ("-k", "--seqdb_api_key"):
            seqdb_api_key = arg
        elif opt == "-p":
            prod=True
        elif opt in ("-u", "--seqdb_ws_url"):
            user_url = arg
            
    if user_url:
        base_url = user_url
    elif prod:
        base_url = prod_url
    else:
        base_url = local_url
    
    return seqdb_api_key, base_url




def main(api_key,base_url):
    #print "Sending request to: " + base_url
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
    try:
        its_region_ids = seqdbWS.getItsRegionIds()
    except requests.exceptions.ConnectionError as e:
        print "Could not connect to Sequence DB. \n", e.message
        sys.exit(1)
    except requests.exceptions.ReadTimeout as e:
        print "Connection too slow for getting data from Sequence DB. \n"
        sys.exit(1)
    except requests.exceptions.HTTPError as e:
        print "HTTP error when getting region ids from Sequence DB. \n", e.message
        sys.exit(1)

   
   
    #Get sequence IDs for the ITS regions
    its_seq_ids = []
    for region_id in its_region_ids:
        try:
            curr_seq_ids = seqdbWS.getSeqIds(region_id)
            its_seq_ids.extend(curr_seq_ids)
        except requests.exceptions.ConnectionError as e:
            print "Could not connect to Sequence DB. \n", e.message
            sys.exit(1)
        except requests.exceptions.ReadTimeout as e:
            print "Connection too slow for getting data from Sequence DB. \n"
            sys.exit(1)
        except requests.exceptions.HTTPError as e:
            print "HTTP error when getting region ids from Sequence DB. \n", e.message
            sys.exit(1)

     
    # Get fasta sequences based on ids and write to a file 
    output_file = open(output_file_name, 'w')
    
    for seq_id in its_seq_ids:
        try:
            # Request sequence in fasto format from SeqDB:
            fastaSequence = seqdbWS.getFastaSeq(seq_id)
            output_file.write(fastaSequence)
        except requests.exceptions.ConnectionError as e:
            print "Could not connect to Sequence DB. \n", e.message
            sys.exit(1)
        except requests.exceptions.ReadTimeout as e:
            print "Connection too slow for getting data from Sequence DB. \n"
            sys.exit(1)
        except requests.exceptions.HTTPError as e:
            print "HTTP error when getting region ids from Sequence DB. \n", e.message
            sys.exit(1)

     
    output_file.close()   
    return its_seq_ids
    


if __name__ == '__main__':
     # Parse command line to get seqdb api key (necessary to request seqDB web services) and base url for web services requests
    api_key, base_url = parse_input_args(sys.argv[1:])
   
    seq_ids = main(api_key, base_url)
    
    print("Loaded %i sequences from Sequence Dababase. \n\n Base url for Sequence Database web services:\n   %s" % (len(seq_ids), base_url))