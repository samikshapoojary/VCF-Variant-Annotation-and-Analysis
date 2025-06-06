#importing the necessary modules
import logging  #for logging errors, info, etc
import argparse  #for creating and handling command-line arguments and take input files from user
import vcf  #for parsing VCF file
import gffutils  #for parsing GFF file
import os  #for interacting with the file system
from Bio import SeqIO  #for parsing FASTA sequence files
from Bio.Seq import Seq  #for sequence translation
import pandas as pd  #for creation and manipulation of tabular data 
import seaborn as sns  #for data visualization as plots
import matplotlib.pyplot as plt  #for creating statistical plots

#creating command-line interface via argparse to take file paths from the user
parser = argparse.ArgumentParser(description='Analyzing variants to determine type and proportion of mutations and their effect on protein if variants are in coding regions', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf', required=True, help='vcf file')  #takes in VCF file containing variant data 
parser.add_argument('--gff', required=True, help='gff file')  #takes in GFF file containing gene annotations and genomic features
parser.add_argument('--fasta', required=True, help='fasta file')  #takes in FASTA file consisting of genome sequence
parser.add_argument('--log', required=True, help='log file')  #takes in log text file to log events occuring during the script execution
parser.add_argument('--output', required=True, help='output file')  #takes in the output file where results will be save
args = parser.parse_args()  #parse arguments provided by the user

#setting up the logger to store information, errors and other events
log_file = args.log  #storing the user provided log file path
logger = logging.getLogger()  #initializing a logger object
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(levelname)s - %(asctime)s - %(message)s') #configuring the logging to log everything above info level
logger.info(f'Storing the log entries in file {args.log}\n')

#creating a function to connect to or create a GFF database for efficient feature lookup
def create_gff_db(gff_file):
    db = gff_file.replace('.gff', '.db')  #change the file extension from .gff to .db for the SQLite database
    if not os.path.isfile(db):  #if the database file does not exist yet, creating it
        logger.info(f'Creating gff database {db}...\n')  #logging the creation process
        try:
            db = gffutils.create_db(gff_file, dbfn=db, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        except Exception as error:
            logger.error(f'Cannot create database {db} due to {error}\n')  #log detailed error message if creation fails
    else:  #if the database already exists, connect to it
        logger.info(f'Connecting to existing database {db}...\n')  #logging the connection
        try:
            return gffutils.FeatureDB(db, keep_order=True)  #return the database object for querying the features
        except ValueError as error:
            logger.error(f'Database {db} could not be connected to due to {error}\n')  #logging error message if connecting to the database fails
    
#creating function to parse the VCF file and filter variants based on quality 
def parse_vcf(vcf_file, qual_threshold):
    hq_variants = []  #store high-quality variants based on quality score
    lq_count = 0  #count low-quality variants
    try: #try opening the VCF file provided by user
        with open(vcf_file, 'r') as file:
            vcf_reader = vcf.Reader(file)  #read the VCF file using the vcf.Reader
            for record in vcf_reader:  #iterate through all variant records in the file
                if record.QUAL > qual_threshold:  #only keep variants with a quality score above the threshold
                    hq_variants.append(record)  #add high-quality variants to the list
                else:
                    lq_count += 1  #count low-quality variants
    except Exception as error:
        logger.error(f'Error parsing vcf file due to {error}\n')  #logging error if the VCF file cannot be parsed or file path not found
    return hq_variants, lq_count  #return the high-quality variants and the low-quality count

#function to parse a FASTA file and store sequences in a dictionary
def parse_fasta(fasta_file):
    fasta_dict = {}  #dictionary to store sequences by their IDs
    try: #try opening the FASTA file provided by user
        with open(fasta_file, 'r') as file:
            for seq in SeqIO.parse(file, 'fasta'):  #using SeqIO to parse the FASTA file.
                fasta_dict[seq.id] = str(seq.seq)  #store the sequence using the sequence ID as the key
    except Exception as error:
        logger.error(f'Error loading fasta file due to {error}\n')  #logging error if the FASTA file cannot be parsed
    return fasta_dict  #return the dictionary

#function to process GFF data and store relevant feature information by chromosome
def parse_gff(gff_db):
    feature_info = {}
    for feature in gff_db.all_features():  #iterate over all features in the GFF database
        chrom, start, end, strand, feature_type = feature.chrom, feature.start, feature.end, feature.strand, feature.featuretype  #extracting feature details  
        if feature_type == 'CDS': #obtaining variants in coding region
            transcript_id = feature.attributes['Parent'][0] #get associated transcript IDs for CDS
        else:
            transcript_id = 'NA' #non-coding features dont have transcript ID, marking as NA
        if chrom not in feature_info: #if the chromosome is not in the dictionary, adding it
            feature_info[chrom] = []
        feature_info[chrom].append({'chrom': chrom, 'start': start, 'end': end, 'strand':strand, 'feature_type':feature_type, 'transcript_id':transcript_id})
    return feature_info #return the dictionary

#function to check if a position is within the boundaries of coding region
def in_coding_region(chrom, pos, feature_info):
    if chrom not in feature_info:
        return False, 'Non-coding'  #return None if chrom not found in dict  
    for feature in feature_info[chrom]: #loop through features
        if feature['feature_type'] == 'CDS' and feature['start'] <= pos <= feature['end']: #check if position is within a CDS region
            return True, feature['transcript_id']  #if position is in a CDS, return True and associated transcript ID
    return False, 'Non-coding'  #return False if no CDS feature is found i.e non-coding variant

#function to calculate the protein location after translation of variant
def get_pro_loc(chrom, pos, feature_info):
    for feature in feature_info[chrom]: #loop through features 
        if feature['feature_type'] == 'CDS' and feature['start'] <= pos <= feature['end']: #check if position is within a CDS region
            if feature['strand'] == '+': #for the forward strand, calculate protein position using the start position
                pro_loc = (pos - feature['start']) // 3  #calculate protein position (0-based in python)
                return pro_loc + 1  #return 1-based protein position
            elif feature['strand'] == '-': #for the reverse strand, calculate codon position using the start position
                pro_loc = (feature['end'] - pos) // 3  #calculate protein position (0-based in python)
                return pro_loc + 1 #return 1-based protein position
    return None  #return None if no valid protein location

#function to generate the mutated codon sequence by replacing the base in the codon at the given position with the alternate allele
def get_mut_codon(coding_seq, alt_allele, pos):
    codon_start = (pos - 1) // 3 * 3  #getting start position of the codon
    codon = coding_seq[codon_start:codon_start + 3]  #extracting the codon
    codon_pos = (pos - 1) % 3  #find the position of the variant within the codon
    mutated_codon = codon[:codon_pos] + alt_allele + codon[codon_pos + 1:] #replace the base at the mutation position with the alternate allele
    return codon, mutated_codon #return the reference codon and mutated codona

#function to create a list of information for each variant, to be added to the output file
def add_variant_info(chrom, pos, ref_allele, alt_allele, mutation_type, transcript_id, pro_loc, ref_aa, alt_aa):
    if mutation_type == 'Synonymous':
        alt_aa = 'NA'  #for synonymous mutations, alternate amino acid is the same, hence NA
    return [chrom, pos, ref_allele, alt_allele, mutation_type, transcript_id, pro_loc, ref_aa, alt_aa]

#main function to carry out analysis on the variants
try:
    #load the necessary files
    gff_db = create_gff_db(args.gff) #create or load the GFF database using the user inputed GFF file
    if gff_db is None:
        logger.error('Exiting as failed to create or load GFF database. \n')
        print('Database creation was not completed. Try running the script again.')
        exit(1)
    fasta_file = args.fasta #user inputed FASTA file
    hq_variants, lq_count = parse_vcf(args.vcf, qual_threshold=20) #parse the VCF file and get high-quality variants based on the threshold of 20

    #log the results to track file inputs and results
    logger.info(f'Processing the VCF file: {args.vcf}')
    logger.info(f'Processing the GFF file: {args.gff}')
    logger.info(f'Processing the FASTA file: {args.fasta}\n')
    logger.info(f'Number of variants with Qual <= 20 i.e low quality variants was {lq_count}\n')
    logger.info(f'Number of variants with Qual > 20 i.e high quality variants was {len(hq_variants)}\n')

    mutation_categories = {'Non-coding': 0, 'Synonymous': 0, 'Non-synonymous': 0} #count for mutation types
    variant_info = [] #list to store information about each variant
    fasta_dict = parse_fasta(fasta_file) #parse the FASTA file to get the genome sequence
    feature_info = parse_gff(gff_db) #parse GFF data to extract features

    for variant in hq_variants: #loop through each high-quality variant to classify and analyze it
        chrom = variant.CHROM #chromosome where the variant is located
        pos = variant.POS #position of the variant on the chromosome
        ref_base = variant.REF #reference base at the variant position
        alt_base = str(variant.ALT[0]) #alternate base at the variant position 

        #check if the variant is in a coding region
        is_coding, transcript_id = in_coding_region(chrom, pos, feature_info)
        
        if is_coding: #if the variant is in a coding region, analyze it further
            pro_loc = get_pro_loc(chrom, pos, feature_info) #get protein location
            coding_seq = fasta_dict[chrom] #get the genomic sequence for the chromosome
            
            #get reference and mutated codons
            ref_codon, mut_codon = get_mut_codon(coding_seq, alt_base, pos)

            #translate codons to amino acids
            ref_aa = str(Seq(ref_codon).translate())
            alt_aa = str(Seq(mut_codon).translate())

            #determine if the mutation is synonymous or non-synonymous on the basis of reference and alternate allelles
            mutation_type = 'Synonymous' if ref_aa == alt_aa else 'Non-synonymous'
            mutation_categories[mutation_type] += 1 #update mutation type count

        else:  #if the variant is not in a coding region, classify it as non-coding
            mutation_type = 'Non-coding'
            pro_loc, ref_aa, alt_aa = 'NA', 'NA', 'NA'  #non-coding variants dont have protein information
            mutation_categories[mutation_type] += 1 #update non-coding count
        
        #add variant info for both coding and non-coding variants
        variant_info.append(add_variant_info(chrom, pos, ref_base, alt_base, mutation_type, transcript_id, pro_loc, ref_aa, alt_aa))

    #write the processed variant information to an output file in tsv format
    with open(args.output, 'w') as output_file:
        output_file.write(f'Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein\tLocation\tRef AA\tAlt AA\n')
        for variant in variant_info:
            output_file.write(f'{variant[0]}\t{variant[1]}\t{variant[2]}\t{variant[3]}\t{variant[4]}\t{variant[5]}\t{variant[6]}\t{variant[7]}\t{variant[8]}\n')
    logger.info(f'Writing output to file {args.output}\n')

    #create a dataframe of mutation type and proportion obtained from mutation counts and total count  
    mutation_proportions = {}
    for mutation_type, count in mutation_categories.items():
        mutation_proportions[mutation_type] = count / len(hq_variants)
    mutation_proportions_df = pd.DataFrame(mutation_proportions.items(), columns=['Mutation Type', 'Proportion'])

    #create plot of Proportion v/s Mutation type
    sns.barplot(data=mutation_proportions_df, x='Mutation Type', y='Proportion')
    plt.title('Proportion of Mutation Types')
    plt.savefig('mut_prop_barplot.png')
    plt.close()
    logger.info(f"Saving the barplot of mutation proportions as 'mut_prop_barplot.png'")

except Exception as error:
    logger.error(f'Processing could not be completed due to {error}')  #log any error in the entire process
    raise
