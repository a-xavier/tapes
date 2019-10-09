'''
    Name: tapes.py
    Author: Alexandre Xavier
    Version; 0.1
    Contact: alexandre.xavier@uon.edu.au
    Python Version: 3.x
    from https://github.com/a-xavier/tapes
    Licence: MIT
'''

#!usr/bin/env python3
import sys
import argparse
import os
import json
from time import time
from numpy import array_split
from multiprocessing.pool import Pool
import src.vep_process as vp
import src.t_func as tf
import src.prs

parser = argparse.ArgumentParser(prog='tapes', usage=tf.help_message(), add_help=False)

parser.add_argument('option', type=str, nargs=1,
                    help="Choose the main function:\n"
                         "annotate: use annovar to annotate a vcf file\n"
                         "sort: use tapes to prioritise variants\n"
                         "db: see or download databases\n"
                         "analyse: use an already analysed file to produce a specific report",
                    choices=['annotate', 'sort', 'db', 'decompose', 'analyse', 'analyze'])

parser.add_argument('-i', '--input',
                    required=False,
                    help="Path of the input csv file annotated by annovar",
                    type=str)
parser.add_argument('-o', '--output',
                    required=False,
                    help="Path of the sorted output csv file",
                    type=str)
parser.add_argument('-b', '--build_db',
                    help="Download databases flagged for download : /humandb will be populated with relevant databases\n"
                         "Also downloads the necessary tapes databases if necessary",
                    action="store_true",
                    required=False)

parser.add_argument('-A', '--annovar',
                    nargs='?',
                    help="Path of the annovar FOLDER",
                    type=str,
                    required=False)
parser.add_argument('-s', '--see_db',
                    help="Path of the annovar FOLDER : /humandb will be searched for relevant databases and write the"
                         "result in db_config.json in tapes current directory",
                    action="store_true",
                    required=False)
parser.add_argument('-e', '--enrichr',
                    help="use --enrichr to analyse top mutations (>0.85 in Total score) with EnrichR API"
                         "Default is GO_Biological_Process_2018, see http://amp.pharm.mssm.edu/Enrichr/#stats to see "
                         "all available database",
                    nargs="?", const="GO_Biological_Process_2018",
                    required=False)

parser.add_argument("-d", "--disease",
                    help="excluse genes non involved in specified disease (based on Disease column description in a new"
                         " spreadsheet in xlsx file, use quotes '' to use a multiple words term (eg; 'multiple sclerosis' "
                         "Default = cancer",
                    nargs="?", const="cancer",
                    required=False)
parser.add_argument('--BIG',
                    help="If the vcf is too big, exported file remove probably benign variants before exporting",
                    nargs='?', const=0.35, default=None, type=float,
                    required=False)
parser.add_argument('--by_sample',
                    help='Add a spreadsheet to xlsx report with the 5 top pathogenic variants',
                    required=False,
                    action="store_true")
parser.add_argument('--list',
                    help="A txt file containing a list of refGene Genes symbols or a string of refGenes Genes symbol "
                         "in quotes",
                    required=False)
parser.add_argument("-a", "--assembly",
                    help="Name of the assembly, hg19 / hg38"
                         "Default is hg19", nargs="?", const="hg19", default="hg19", choices={"hg19", "hg38"},
                    required=False)
parser.add_argument("--ref_anno",
                    help="Genes format to use RefSeq 'refGene',  ENSEMBL 'ensGene' or USCS 'knownGene "
                         "default is refGene", nargs="?", const="refGene", default="refGene",
                    choices={"refGene", "ensGene", "knownGene"},
                    required=False)
parser.add_argument("--acmg",
                    help="Do the ACMG criteria computation",
                    action='store_true',
                    required=False)
parser.add_argument("--trio",
                    help="Specify pedigree data for trio analysis as .txt file' "
                         "default is None", nargs='?', default=None,
                    required=False)
parser.add_argument("--kegg",
                    help="Similar to --list."
                         " Indicate a kegg pathway and tapes will return a spreasheet with"
                         "genes only contained in that pathway",
                    required=False)
parser.add_argument('--cutoff',
                    help='Select cutoff for BS1.'
                         'Frequency for rare disease, default is 0.005', type=float,
                    nargs=1,
                    default=0.005,
                    required=False)
parser.add_argument('--acmg_db',
                    help='The location of the databases used to assign ACMG criteria\n'
                    'Default is the directory where tapes.py is located',
                    nargs='?',
                    default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'acmg_db'),
                    required=False)
parser.add_argument("--tab",
                    help="Export in tab-separated values format instead of csv format",
                    action='store_true',
                    required=False)
parser.add_argument('--by_gene',
                    help='Add a spreadsheet to xlsx report with gene-grouped variants',
                    required=False,
                    action="store_true")
parser.add_argument('-t','--threads',
                    help='numbers of threads / Unused for now', type=int,
                    default=1,
                    required=False)
parser.add_argument('--pp2_percent',
                    help='Threshold for PP2, considering PP2 positive if variant is missense in a gene where more than the threshold are pathogenic missense', type=int,
                    default=80,
                    required=False)
parser.add_argument('--pp2_min',
                    help='Number of minimum pathogenic missense variants in gene to consider PP2', type=int,
                    default=1,
                    required=False)
parser.add_argument('--bp1_percent',
                    help='Threshold for BP1, considering BP1 positive if variant is missense in a gene where less than the threshold are pathogenic missense', type=int,
                    default=15,
                    required=False)
parser.add_argument("--test",
                    help="test flag - Will override all options and run a test of the sort function",
                    action='store_true',
                    required=False)
parser.add_argument('--prs',
                    help='Disease or trait to calculate Polygenic Risk Score against', type=str,
                    required=False)
parser.add_argument('--prs_num',
                    help='Number of public samples to get genotype from', type=int,
                    default=100,
                    required=False)

args = parser.parse_args()


if args.test == True:
    args = parser.parse_args(['sort','-i./toy_dataset/test_input.csv', '-o./test_output/', '--tab', '--test', '--enrichr', '--by_gene' , '--by_sample'])
    args.list = "FH CDC73"
    args.kegg = "Pathways in cancer"
    args.disease = "autosomal dominant"
    args.trio = "./toy_dataset/trio.txt"
    args.prs = 'Schizophrenia'
    args.prs_num = 10
    for item in vars(args):
        print(item, '=', getattr(args, item))
else:
    args = parser.parse_args()
    
# Process paths to transform into absolute 
if args.annovar:
    args.annovar = os.path.abspath(args.annovar)
if args.input != None:
    args.input = os.path.abspath(args.input)
if args.output != None:
    if args.output[-4:] == ('.csv' or '.vcf' or '.txt' or '.tsv'):
        args.output = os.path.abspath(args.output)
    elif args.output[-1] == '/':
        args.output = (os.path.abspath(args.output)+'/')
    elif args.output[-1] =='\\':
        args.output = (os.path.abspath(args.output)+'\\')
    else:
        pass

# Convert inputs and outputs into absolute paths
script_absolute_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_absolute_path)


# PROCESS ACMG_DB PATH
try:
    with open('./src/db_config.json', 'r') as db_cfg:
        cfg_dict = json.load(db_cfg)
    try:
        acmg_db_path = cfg_dict["acmg_db_path"]['acmg_db_path']
    except KeyError:
        acmg_db_path = args.acmg_db
except FileNotFoundError:
    print("No acmg_db path given and no db_config.json found\n"
          "Default is: {}".format(args.acmg_db))
    acmg_db_path = args.acmg_db


def main():
    start_time = time()
    output_type = tf.output_type(args.output)  # Determine output type either directory or csv or txt/tsv
    # Do not discriminate by platform anymore so win users can also use normal slash
    if output_type == 'directory':
        if '/' in args.output:
            output_prefix = args.output.split('/')[-2]
        elif '\\' in args.output:
            output_prefix = args.output.split('\\')[-2]

    file_path = tf.process_path(args.input)  # process relative to absolute path

    full_stuff, soft_used = tf.open_csv_file(file_path, acmg_db_path)    # Load annotated csv in pandas
    print(tf.tmp_stmp()+'{} annotated file'.format(soft_used))

    ref_anno = tf.reference_used(full_stuff)      # Detect reference annotation used (refseq/UCSC/ENSEMBL)


    # IF NUMBER OF SAMPLE UNDER 15 PS4 will not be calculated using fisher's test but with the normal db

    number_of_samples = tf.counting_number_of_samples(full_stuff)

    tf.check_if_data_already_processed(full_stuff)  # Check if file was already processed

    full_stuff = tf.replace_dash_by_dots(full_stuff)  # Replace dashes by dots in ref/alt columns


    ################################ PROCESSING AND CLASSIFYING VARiANTS ################################
    # ACMG TAG CHECKS FOR ALL ANNOTATION REQUIREMENT
    if args.acmg is True:
        tf.check_acmg_tag_requirement(full_stuff, ref_anno)
    elif args.acmg is False:
        pass

    ### TESTS ###
    # FIRST LOADS ALL NECESSARY DATABASES
    global tup_databases
    if soft_used == 'annovar':
        tup_databases = tf.load_databases(acmg_db_path, ref_anno, args.assembly, args.trio, args.pp2_percent, args.pp2_min, args.bp1_percent)
    elif soft_used == 'vep':
        tup_databases = vp.load_vep_databases(acmg_db_path, args.assembly, ref_anno, args.trio, args.pp2_percent, args.pp2_min, args.bp1_percent)

    # IF TRIO DATA IS THERE EXCLUDE HEALTHY PARENTS
    if args.trio is not None:
        print(tf.tmp_stmp() + 'Excluding Healthy Parents from trios...', end='', flush=True)
        full_stuff_tup = tf.process_df_for_trio(tup_databases, full_stuff)
        full_stuff = full_stuff_tup[0]
        list_healthy = full_stuff_tup[1]
        print('Done')
    else:
        list_healthy = []
    # CHECK FOR MULTIPROCESSING
    if args.threads >= 2:
        # TEST FOR MULTIPROCESSING
        split_df = array_split(full_stuff, args.threads) # Split dataframe
        list_args = []
        for chunck in split_df:
            toup_argument = (chunck, ref_anno, number_of_samples, tup_databases, args.cutoff)
            list_args.append(toup_argument)

        #TODO MAYBE ADD MULTIPROCESSING TO GO FASTER
        # problem HIGH USAGE OF RAM FOR EACH PROCESS (Because of databases copied to child processes)

        pool = Pool(processes=args.threads)
        if soft_used == "annovar":
            result = pool.starmap(tf.process_data, list_args)
        elif soft_used == 'vep':
            result = pool.starmap(vp.process_data, list_args)
        pool.close()
        pool.join()
        final_stuff = tf.pd.concat(result)  # Merge all pieces for final dataframe
        if soft_used == 'vep':
            final_stuff = vp.rename_columns(final_stuff, ref_anno)
        final_stuff.drop(list_healthy, axis=1, inplace=True)
        ###END MUTLTIPROCESS###

    else:
        if soft_used == 'annovar':
            final_stuff = tf.process_data(full_stuff, ref_anno, number_of_samples, tup_databases, args.cutoff)
            final_stuff.drop(list_healthy, axis=1, inplace=True)
        elif soft_used == 'vep':
            final_stuff = vp.process_data(full_stuff, ref_anno, number_of_samples, tup_databases, args.cutoff)
            final_stuff = vp.rename_columns(final_stuff, ref_anno)
            final_stuff.drop(list_healthy, axis=1, inplace=True)

    final_stuff = tf.sort_by_pathogenicity(final_stuff)


    ##############################  OUTPUT RAW REPORT  ###############################################

    if args.BIG is not None:
        final_stuff = final_stuff[final_stuff['Probability_Path'] >= float(args.BIG)]

    if output_type == 'csv' or output_type == 'txt':
        tf.export_sorted_csv(final_stuff, tf.process_path(args.output), args.tab)
    elif output_type == 'directory':
        if args.tab is True:
            ext = '.txt'
        else:
            ext = '.csv'
        outfile = os.path.join(args.output, (output_prefix+ext))
        tf.export_sorted_csv(final_stuff, outfile, args.tab)

    ################################ OPTIONS ############################################################

    if args.kegg is not None:
        tf.kegg_inverted_patway(args.kegg, final_stuff, ref_anno, acmg_db_path)

    if args.by_gene is True:
        tf.gene_grouped_report(final_stuff, ref_anno, acmg_db_path, number_of_samples)

    if args.enrichr is not None:
        tf.pathway_analysis(final_stuff, args.enrichr, ref_anno)

    if args.by_sample is True:
        tf.individual_reports(final_stuff, ref_anno)

    if args.disease is not None:
        tf.create_cancer_dataframe(final_stuff, args.disease, ref_anno)

    if args.list is not None:
        if os.path.isfile(args.list) is True:
            with open(args.list, 'r') as list_file:
                gene_list = []
                for line in list_file.readlines():
                    gene_list.append(line.rstrip('\n'))
        else:
            gene_list = set(args.list.split(" "))
        tf.find_list_of_genes(final_stuff, gene_list, ref_anno)

    if args.prs is not None:
        output_folder = os.path.dirname(args.output)
        src.prs.prs_calculation(acmg_db_path, final_stuff, args.output, args.prs, args.prs_num, output_type)
        print(tf.tmp_stmp()+"PRS Done")

    ################################## OUTPUT EXTRA REPORT ################################################
    if output_type == 'csv':
        path_output = args.output.replace('.csv', '.xlsx')
        tf.export_to_formatted_xlsx(final_stuff, path_output, number_of_samples, ref_anno)
        print(tf.tmp_stmp() + "Output written : {}/xlsx".format(tf.process_path(args.output)))
    elif output_type == 'txt':
        no_ext = args.output[:-4]
        path_output = no_ext+'.xlsx'
        tf.export_to_formatted_xlsx(final_stuff, path_output, number_of_samples, ref_anno)
        print(tf.tmp_stmp() + "Output written : {}/xlsx".format(tf.process_path(args.output)))
    elif output_type == 'directory':
        tf.csv_report(final_stuff, args.output, output_prefix, ref_anno, args.tab)

    stop_time = time()
    print('|| Process finished in ' + str(tf.ceil(stop_time - start_time)) + " seconds")

# USED FOR TESTING PURPOSE IN TRAVIS
def test_main():
    print("ACTIVATING MAIN TEST")
    main()
    sys.exit(0)


if __name__ == "__main__":
    if args.test == True:
        test_main()

    # Process data from annovar annotated files # EVEN IF NOT FULLY FREESOME COMPLIANT without acmg flag
    elif (args.output and args.input and args.assembly) and args.option == ['sort']:
        print('''
        ***TAPES: SORT***
        ''')
        main()

    # BUILD annovar DATABASES
    elif args.option == ['db'] and args.build_db is True:
        print('''
        ***TAPES: DOWNLOAD DATABASE***
        ''')
        # PROCESS annovar PATH
        if args.annovar is None:  # If annovar folder not supplied, look at the config file
            try:
                with open('./src/db_config.json', 'r') as db_cfg:
                    cfg_dict = json.load(db_cfg)
                annovar_path = cfg_dict["annovar_path"]['annovar_path']
            except FileNotFoundError:
                print("No annovar path given and no db_config.json found")
        else:
            annovar_path = args.annovar
        tf.build_annovar_db(annovar_path, args.assembly, args.acmg)
        tf.check_for_acmg_databases(acmg_db_path)
        tf.check_online_annovar_dbs(annovar_path)
        tf.scan_for_db(annovar_path)
        tf.writing_db_to_file(annovar_path, args.acmg_db)
        tf.amcg_db_download(acmg_db_path)

    # ANNOTATE VCF FILE
    elif (args.input and args.output) is not None and args.option == ['annotate']:
        print('''
        ***TAPES: ANNOTATE***
        ''')
        # PROCESS annovar PATH
        if args.annovar is None:  # If annovar folder not supplied, look at the config file
            try:
                with open('./src/db_config.json', 'r') as db_cfg:
                    cfg_dict = json.load(db_cfg)
                annovar_path = cfg_dict["annovar_path"]['annovar_path']
            except FileNotFoundError:
                print("No annovar path given and no db_config.json found")
        else:
            annovar_path = args.annovar
        if args.output[-4:] == '.csv' or args.output[-4:] == '.txt':
            tf.process_annotate_vcf(args.input, args.output, annovar_path, args.assembly, args.ref_anno, args.acmg)
            tf.merge_sample_info_in_annotated_csv(args.output, args.assembly)
        elif args.output[-4:] == '.vcf':
            tf.process_annotate_vcf(args.input, args.output, annovar_path, args.assembly, args.ref_anno, args.acmg)
        else:
            print('|| Could not determine the desired output type, please use either ".csv" or ".vcf" extension for the'
                  ' output')
            tf.sys.exit(1)

    # SCAN annovar DATABASES ALREADY DOWNLOADED AND STORE THEM INSIDE OF A FILE

    elif args.option == ['db'] and args.see_db is True:
        print('''
        ***TAPES: SEE DATABASE***
        ''')
        # PROCESS annovar PATH
        if args.annovar is None:  # If annovar folder not supplied, look at the config file
            try:
                with open('./src/db_config.json', 'r') as db_cfg:
                    cfg_dict = json.load(db_cfg)
                annovar_path = cfg_dict["annovar_path"]['annovar_path']
            except FileNotFoundError:
                print("No annovar path given and no db_config.json found")
        else:
            args.annovar = os.path.abspath(args.annovar)
            annovar_path = args.annovar
        tf.check_online_annovar_dbs(annovar_path)
        tf.scan_for_db(annovar_path)
        tf.writing_db_to_file(annovar_path, args.acmg_db)

    #  Decompose vcf
    elif args.option == ['decompose'] and args.input is not None and args.output is not None :
        print('''
        ***TAPES: DECOMPOSE***
        ''')
        print(tf.tmp_stmp()+'Decomposing...', end=' ')
        tf.decompose_vcf(args.input, args.output)
        print(tf.tmp_stmp()+"Done")

    # Re-analyse TODO REFINE PROCESS for just the analysis wanted (don't do xlsx report)
    elif (args.output and args.input and args.assembly) and (args.option == ['analyse'] or args.option == ['analyze']):
        print('''
        ***TAPES: RE-ANALYSE***
        ''')
        final_stuff, soft_used = tf.open_csv_file(args.input, acmg_db_path)
        ref_anno = args.ref_anno
        number_of_samples = tf.counting_number_of_samples(final_stuff)
        output_type = tf.output_type(args.output)  # Determine output type either directory or csv or txt/tsv
        if output_type == 'directory':
            print('|| Please use output type : txt or csv')

        if args.BIG is not None:
            final_stuff = final_stuff[final_stuff['Probability_Path'] >= float(args.BIG)]

        if args.kegg is not None:
            tf.kegg_inverted_patway(args.kegg, final_stuff, ref_anno, acmg_db_path)
            
        if args.prs is not None:
            output_folder = os.path.dirname(args.output)
            src.prs.prs_calculation(acmg_db_path, final_stuff, args.output, args.prs, args.prs_num, output_type)

        if args.by_gene is True:
            tf.gene_grouped_report(final_stuff, ref_anno, acmg_db_path, number_of_samples)

        if args.enrichr is not None:
            tf.pathway_analysis(final_stuff, args.enrichr, ref_anno)

        if args.by_sample is True:
            tf.individual_reports(final_stuff, ref_anno)

        if args.disease is not None:
            tf.create_cancer_dataframe(final_stuff, args.disease, ref_anno)

        if args.list is not None:
            if os.path.isfile(args.list) is True:
                with open(args.list, 'r') as list_file:
                    gene_list = []
                    for line in list_file.readlines():
                        gene_list.append(line.rstrip('\n'))
            else:
                gene_list = set(args.list.split(" "))
            tf.find_list_of_genes(final_stuff, gene_list, ref_anno)

        ################################## OUTPUT EXTRA REPORT ################################################
        if output_type == 'csv':
            tf.csv_report_reanalyse(args.output, False)
            print(tf.tmp_stmp() + "Done")
        elif output_type == 'txt':
            tf.csv_report_reanalyse(args.output, True)
            print(tf.tmp_stmp() + "Done")

    else:
        print(tf.help_message())
        print("Wrong combination of arguments")
