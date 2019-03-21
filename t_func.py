#!usr/bin/env python3
import subprocess
import sys
import os
import requests
import json
from xlsxwriter import utility
from numpy import ceil
import gzip
from math import sqrt, log1p, exp
from datetime import datetime
from vcf_parser import VCFParser  # FROM https://github.com/moonso/vcf_parser
import pandas as pd
from pandas import errors
from scipy.stats import fisher_exact
from difflib import get_close_matches
if 'win' not in sys.platform:
    from pysam import VariantFile


def tmp_stmp():
    return str(datetime.now())[:19] + '.....'


def add_list_to_df(full_table, user_list, name_of_column):
    user_list_in_df = pd.Series(user_list)  # transform the vector in pandas series
    full_table[name_of_column] = user_list_in_df.values
    return full_table


def open_csv_file(csv_file):
    if os.path.isfile(csv_file) and csv_file[-4:] == '.csv':
        try:
            dataframe = pd.read_csv(csv_file, low_memory=False)
            dataframe[['Start', 'End']] = dataframe[['Start', 'End']].astype(int)
        except errors.ParserError:
            tmp_file = csv_file+'_tmp.csv'
            print(tmp_stmp()+'Parser error, processing input file...', end="")
            with open(csv_file, 'r') as infile:
                with open(tmp_file, 'w') as outfile:
                    for line in infile.readlines():
                        nline = line.replace(', ', ' ')
                        outfile.write(nline)
            dataframe = pd.read_csv(tmp_file, low_memory=False)
            dataframe[['Start', 'End']] = dataframe[['Start', 'End']].astype(int)
            os.remove(tmp_file)
            print("Done")
        return dataframe
    elif os.path.isfile(csv_file) and csv_file[-4:] == '.vcf':
        dataframe = process_annotated_vcf(csv_file)
        dataframe[['Start', 'End']] = dataframe[['Start', 'End']].astype(int)
        return dataframe
    elif os.path.isfile(csv_file) and (os.path.isfile(csv_file) and csv_file[-4:] == '.txt' or os.path.isfile(csv_file) and csv_file[-4:] == '.tsv'):
        dataframe = pd.read_csv(csv_file, sep='\t', low_memory=False) # TODO maybechange that, header=0, engine='python',index_col=False)
        dataframe[['Start', 'End']] = dataframe[['Start', 'End']].astype(int)
        return dataframe
    elif os.path.isfile(csv_file) is False:
        print('File does not exist')
        sys.exit(1)
    else:
        print('Could not open {}'.format(csv_file))
        sys.exit(1)


def check_if_data_already_processed(full_stuff):
    if 'Prediction_ACMG_tapes' in full_stuff.columns and 'Probability_Path' in full_stuff.columns:
        print('This file has already been processed by Tapes')
        sys.exit(0)
    else:
        pass
    return


def seg_dup_check(full_stuff):
    # THEN LET's FILTER OUT SEGMENTAL DUPLICATION

    superdup_string = full_stuff['genomicSuperDups']
    superdup_contrib = []

    for superdup_pred in superdup_string:
        if 'Score=' in superdup_pred:
            superdup_contrib.append(-250)
        else:
            superdup_contrib.append(0)

    superdup_in_df = pd.Series(superdup_contrib)
    full_stuff['Genomic_SuperDup_Coef'] = superdup_in_df.values

    return full_stuff


def sort_by_pathogenicity(export_dataframe):
    path_df = export_dataframe.loc[export_dataframe['Prediction_ACMG_tapes'] == 'Pathogenic'].sort_values(
        by=['Probability_Path'], ascending=False)

    like_path_df = export_dataframe.loc[
        export_dataframe['Prediction_ACMG_tapes'] == 'Likely Pathogenic'].sort_values(by=['Probability_Path'],
                                                                                         ascending=False)

    VUS_df = export_dataframe.loc[export_dataframe['Prediction_ACMG_tapes'] == 'VUS'].sort_values(
        by=['Probability_Path'], ascending=False)

    like_ben_df = export_dataframe.loc[export_dataframe['Prediction_ACMG_tapes'] == 'Likely Benign'].sort_values(
        by=['Probability_Path'], ascending=False)

    ben_df = export_dataframe.loc[
        export_dataframe['Prediction_ACMG_tapes'].isin(['Benign', 'Benign auto'])].sort_values(
        by=['Probability_Path'], ascending=False)

    frames = [path_df, like_path_df, VUS_df, like_ben_df, ben_df]
    final_df = pd.concat(frames)

    return final_df


def export_sorted_csv(export_dataframe, export_filepath, tab_flag):
    if tab_flag is True or export_filepath[-4:] == '.txt' or export_filepath[-4:] == '.tsv':
        separator = '\t'
    else:
        separator = ','
    export_dataframe.to_csv(path_or_buf=export_filepath, index=False, sep=separator)
    return


def process_path(entered_path):  # takes relatives or absolute paths and transforms it into absolute pathw
    processed_path = os.path.abspath(entered_path)
    return processed_path


def build_annovar_db(path_to_annovar_folder_rel, assembly, acmg_tag):
    current_wd = os.path.dirname(os.path.realpath(__file__))
    config_file = os.path.join(current_wd, 'db_config.json')
    path_to_annovar_folder = process_path(path_to_annovar_folder_rel)
    if os.path.isdir(path_to_annovar_folder):
        os.chdir(path_to_annovar_folder)
        print('This is where the annovar folder is : {}'.format(path_to_annovar_folder))
        if os.path.isdir('humandb') and os.path.isfile('annotate_variation.pl'):
            print(tmp_stmp()+'Found /humandb directory and annotate_variation.pl')
            if os.path.isfile("example/gene_fullxref_tapes.txt"):
                print(tmp_stmp()+"Found tapes xref file")
            else:
                print('Processing full_xref file to full_xref_tapes. Please wait...')
                subprocess.run(["""tr -d '"' <example/gene_fullxref.txt > example/gene_fullxref_tapes.txt"""], shell=True)
                print('Done')

            if os.path.isfile(config_file) and acmg_tag is False:
                with open(config_file) as cf:
                    conf_dict = json.load(cf)
                    list_of_db_to_dl_hg_19 = []
                    list_of_db_to_dl_hg_38 = []

                for key in conf_dict["hg19"]:
                    db = conf_dict["hg19"][key].strip().lower()
                    if "down" in db :
                        list_of_db_to_dl_hg_19.append(key)

                for key in conf_dict["hg38"]:
                    db = conf_dict["hg38"][key].strip().lower()
                    if "down" in db :
                        list_of_db_to_dl_hg_38.append(key)

                if len(list_of_db_to_dl_hg_19) == 0 and len(list_of_db_to_dl_hg_38) == 0:
                    print(tmp_stmp()+"Nothing is flagged for download, exiting...")
                    return

                print(tmp_stmp()+'Some databases flagged for download. They will be downloaded now '
                      ''
                      ' (THIS WILL TAKE SOME TIME):'
                      '')

                for telecharger in list_of_db_to_dl_hg_19:
                    if telecharger in reg_anno():
                        print(tmp_stmp()+"Using ANNOVAR to download {}".format(telecharger))
                        subprocess.run(['perl', 'annotate_variation.pl', '-build', 'hg19',
                                        '-downdb', '{}'.format(telecharger), 'humandb/'])
                    else:
                        print(tmp_stmp()+"Using ANNOVAR to download {}".format(telecharger))
                        subprocess.run(['perl', 'annotate_variation.pl', '-buildver', 'hg19', '-downdb', '-webfrom'
                                           ,'annovar', '{}'.format(telecharger), 'humandb/'])

                for telecharger in list_of_db_to_dl_hg_38:
                    if telecharger in reg_anno():
                        print(tmp_stmp()+"Using ANNOVAR to download {}".format(telecharger))
                        subprocess.run(['perl', 'annotate_variation.pl', '-build', 'hg38', '-downdb',
                                        '{}'.format(telecharger), 'humandb/'])
                    else:
                        print(tmp_stmp()+"Using ANNOVAR to download {}".format(telecharger))
                        subprocess.run(['perl', 'annotate_variation.pl', '-buildver', 'hg38', '-downdb',
                                        '-webfrom', 'annovar', '{}'.format(telecharger), 'humandb/'])
                os.chdir(current_wd)
                print("All done, you can now use tapes to prioritise your variants")
            elif os.path.isfile(config_file) and acmg_tag is True:

                with open(config_file) as cf:
                    conf_dict = json.load(cf)

                telecharger = []
                print(tmp_stmp()+'{} will be downloaded'.format(telecharger))

                for db in acmg_db()[0].split(','):
                    if conf_dict[assembly][db].lower() == 'missing':
                        telecharger.append(db)
                for item in telecharger:
                    print(tmp_stmp() + "Using ANNOVAR to download {}".format(item))
                    subprocess.run(['perl', 'annotate_variation.pl', '-build', '{}'.format(assembly),
                                    '-downdb', '-webfrom', 'annovar', '{}'.format(item), 'humandb/'])
            else:
                print(tmp_stmp()+"No config file found, please run 'python tapes.py --see-db -A /path/to/annovar/folder' first ")
                sys.exit(1)

        elif not os.path.isdir('humandb'):
            print(tmp_stmp()+'No /humandb directory found. Exiting...')
            sys.exit(1)
        elif not os.path.isfile('annotate_variation.pl'):
            print(tmp_stmp()+'annotate_variation.pl not found. Exiting...')
            sys.exit(1)
    elif not os.path.isdir(path_to_annovar_folder):
        print(tmp_stmp()+"This is not an existing directory. Exiting...")
        sys.exit(1)
    return


def acmg_db():
    return ('avsnp150,clinvar_20180603,dbnsfp35c,dbscsnv11,gnomad_genome,gnomad_exome', 'f,f,f,f,f,f')


def process_annotate_vcf(vcf_path, output_path, path_to_annovar_folder_rel, assembly, generef, acmg_tag):
    #avdb_list_file = os.path.join(path_to_annovar_folder_rel, '{}_avdblist.txt'.format(assembly))

    tapes_wd = os.path.dirname(os.path.realpath(__file__))

    global vcf_path_abs
    vcf_path_abs = process_path(vcf_path)

    if output_path[-4:] == '.csv' or output_path[-4:] == '.txt':
        vcf_out = False
    elif output_path[-4:] == '.vcf':
        vcf_out = True

    output_path_abs = process_path(output_path)

    output_path_abs_no_ext = os.path.splitext(output_path_abs)[0]

    path_to_annovar_folder = process_path(path_to_annovar_folder_rel)

    output_name = output_path.split('/')[-1]

    output_path_only = output_path_abs.replace(output_name, '')

    vcf_db_file = os.path.join(os.getcwd(), 'db_vcf.json')
    if os.path.isdir(path_to_annovar_folder):
        os.chdir(path_to_annovar_folder)
        if os.path.isfile('convert2annovar.pl') and os.path.isfile('table_annovar.pl') and os.path.isfile('example/gene_fullxref_tapes.txt'):
            print('Found convert2annovar.pl / table_annovar.pl and the tapes version of fullxref')

            vcf_path_abs = process_input_variant_file(vcf_path_abs)
            decomp_test = test_if_decomposed(vcf_path_abs)
            vcf_path_abs = decomp_test[0]  # Test if the vcf file is decomposed or not
            was_decomposed = decomp_test[1] # if not creates a decomposed vcf and then annotates that

            # Read the config file and create valid annovar command
            if acmg_tag is False:
                with open(vcf_db_file, 'r') as cf: # read config file
                    dict = json.load(cf)
                    list_to_use = []
                for key in dict[assembly]:
                    if dict[assembly][key] == 'YES':
                        list_to_use.append(key) # create the list of all databases present
                # Remove unusable databases or or other acessory databases
                remove_list = ['avdblist', 'kgXref', 'knownGeneMrna.fa', 'ensGeneMrna.fa', 'refGeneMrna.fa', 'refGeneVersion' ]
                list_to_use = [n for n in list_to_use if n not in remove_list]

                # Check version database
                '''version_df = pd.read_csv(avdb_list_file, sep='\t')
                version_df.columns = ['db', 'date', 'size']'''

                to_remove = ["knownGene", "refGene", 'ensGene']  # because it should be first and is included in
                #  the command line which ref to use
                final_list = [x for x in list_to_use if x not in to_remove]
                final_list_str = ','.join(final_list)
                string_type = []
                db_type = db_type_bank()
                # ASSING THE TYPE OF OPERATION AND BEST GUESS IS F IF THE TYPE IS UNKNOWN
                for item in final_list:
                    try:
                        string_type.append(db_type[item])
                    except KeyError:
                        string_type.append('f')
                string_type_str = ",".join(string_type)

            elif acmg_tag is True:
                print(tmp_stmp()+'Annotating using preset settings for ACMG classification')
                final_list_str = acmg_db()[0]
                string_type_str = acmg_db()[1]

            if vcf_out is False and output_path_abs[-4:] == ".csv":

                print(tmp_stmp()+'Converting to avinput')
                subprocess.run(['perl', 'convert2annovar.pl', '-format', 'vcf4', '{}'.format(vcf_path_abs), '-outfile',
                                '{}.avinput'.format(output_path_abs_no_ext), '-allsample', '-withfreq', '-includeinfo'])
                print(tmp_stmp()+'Converted to avinput at {}.avinput'.format(output_path_abs_no_ext))
                subprocess.run(['perl', 'table_annovar.pl', '{}.avinput'.format(output_path_abs_no_ext), 'humandb/',
                                '-buildver', '{}'.format(assembly), '-out', '{}'.format(output_path_abs_no_ext),
                                '-remove', '-protocol',
                                '{},{}'.format(generef, final_list_str),
                                '-operation', 'gx,{}'.format(string_type_str), '-nastring', '.', '-csvout', '-polish', '-xref',
                                'example/gene_fullxref_tapes.txt'])

                # '-thread', '{}'.format(multiprocessing.cpu_count()) # LEAVE MULTITHREAD OUT FOR NOW BECAUSE IT CRASHED MY OWN COMPUTER

                os.rename('{}.{}_multianno.csv'.format(output_path_abs_no_ext, assembly), '{}'.format(os.path.join(output_path_only, output_name)))

                #subprocess.run(['mv', '{}.{}_multianno.csv'.format(output_path_abs_no_ext, assembly), '{}'.format(os.path.join(output_path_only, output_name))])

            elif vcf_out is False and output_path_abs[-4:] == ".txt":

                print(tmp_stmp()+'Converting to avinput')
                subprocess.run(['perl', 'convert2annovar.pl', '-format', 'vcf4', '{}'.format(vcf_path_abs), '-outfile',
                                '{}.avinput'.format(output_path_abs_no_ext), '-allsample', '-withfreq', '-includeinfo'])
                print(tmp_stmp()+'Converted to avinput at {}.avinput'.format(output_path_abs_no_ext))
                subprocess.run(['perl', 'table_annovar.pl', '{}.avinput'.format(output_path_abs_no_ext), 'humandb/',
                                '-buildver', '{}'.format(assembly), '-out', '{}'.format(output_path_abs_no_ext),
                                '-remove', '-protocol',
                                '{},{}'.format(generef, final_list_str),
                                '-operation', 'gx,{}'.format(string_type_str), '-nastring', '.', '-polish', '-xref',
                                'example/gene_fullxref_tapes.txt'])

                # '-thread', '{}'.format(multiprocessing.cpu_count()) # LEAVE MULTITHREAD OUT FOR NOW BECAUSE IT CRASHED MY OWN COMPUTER

                #subprocess.run(['mv', '{}.{}_multianno.txt'.format(output_path_abs_no_ext, assembly), '{}'.format(os.path.join(output_path_only, output_name))])
                os.rename('{}.{}_multianno.txt'.format(output_path_abs_no_ext, assembly),
                          '{}'.format(os.path.join(output_path_only, output_name)))

            elif vcf_out is True:
                subprocess.run(['perl', 'table_annovar.pl', '{}'.format(vcf_path_abs), 'humandb/',
                                '-buildver', '{}'.format(assembly), '-out', '{}'.format(output_path_abs_no_ext),
                                '-remove', '-protocol',
                                '{},{}'.format(generef, final_list_str),
                                '-operation', 'gx,{}'.format(string_type_str), '-nastring', '.', '-vcfinput',
                                '-xref',
                                'example/gene_fullxref_tapes.txt'])

            os.chdir(tapes_wd) # Back to current working directory

        else:
            print('convert2annovar.pl and table_annovar.pl not found in {} or tapes_fullxref in /example'.format(path_to_annovar_folder))
            sys.exit(1)

    elif not os.path.isdir(path_to_annovar_folder):
        print("{} is not an existing directory. Exiting...".format(path_to_annovar_folder))
        sys.exit(1)
    return


def merge_sample_info_in_annotated_csv(csv_file, assembly):
    print(tmp_stmp()+'Merging Sample genotype and Annotated csv/txt...', end=" ", flush=True)
    abs_csv_file = process_path(csv_file)

    #real_csv_file = csv_file.replace('.csv', '.{}_multianno.csv'.format(assembly)) # using annovar the actual output has hg multianno appended
    if csv_file[-4:] == '.csv':
        full_stuff = pd.read_csv(csv_file, sep=',') # Full dataframe annotated
    elif csv_file[-4:] == '.txt':
        full_stuff = pd.read_csv(csv_file, sep='\t') # Full dataframe annotated

    if csv_file[-4:] == '.csv':
        export_filepath = csv_file.replace('.csv', '_with_samples.csv') # Path and name of exported file

        avinput_path = abs_csv_file.replace('.csv', '.avinput') # USE SAMPLE DATA FROM AVINPUT INSTEAD OF VCF
    elif csv_file[-4:] == '.txt':
        export_filepath = csv_file.replace('.txt', '_with_samples.txt')  # Path and name of exported file

        avinput_path = abs_csv_file.replace('.txt', '.avinput')  # USE SAMPLE DATA FROM AVINPUT INSTEAD OF VCF

    # extract sample names in vcf, skipping the header and capturing only the first line
    with open(vcf_path_abs, 'r') as vf:
        for line in vf.readlines():
            if "##" in line:
                pass
            else:
                headers_in_vcf = line
                break
    headers_in_vcf = headers_in_vcf.split()  # Split the vcf column line into a list

    headers_in_vcf = headers_in_vcf[8:]  # Removing first columns from the header (the one not related to the

    global number_of_samples

    number_of_samples = len(headers_in_vcf) - 1  # -1 because we keep the FORMAT column for future uses (mainly to
    # detect the presence of sample data or not
    if 900 > number_of_samples > 1:
        avinput_df = pd.read_csv(avinput_path, header=None, sep="\t")

        avinput_df.drop(avinput_df.columns[range(0, 16)], axis=1, inplace=True)
        avinput_df.columns = headers_in_vcf
        for column in list(avinput_df.columns):
            if column == 'FORMAT':
                pass
            else:
                avinput_df[column].replace('0', '0/0', inplace=True)

        # SPLIT THE GENOTYPE OF EACH SAMPLES
        sep = ':'
        for column in avinput_df.columns:
            value = avinput_df[column].str.split(sep, n=1, expand=True)
            avinput_df[column] = value[0]

        if full_stuff.shape[0] != avinput_df.shape[0]:
            print("\n"+tmp_stmp()+"Exiting because Avinput, vcf and annotated csv don't have the same number of rows. Try Decomposing your vcf")
            sys.exit(1)
        else:
            pass

        new_df = full_stuff.merge(avinput_df, left_index=True, right_index=True)

        hom_list = []
        het_list = []
        norm_list = []

        for line in range(0, new_df.shape[0]):
            values_of_line = new_df.iloc[line].values
            hom_count = 0
            het_count = 0
            normal_count = 0
            for i in values_of_line:
                if i in ["1/1", "2/2", "1|1", "2|2", '1/2', '1|2', '2|1']:
                    hom_count = hom_count + 1
                elif i in ['0/1', '0/2', '0|1', '0|2', '1|0', '2|0']:
                    het_count = het_count + 1
                elif i == '0/0' or i == '0|0':
                    normal_count = normal_count + 1
            hom_list.append(hom_count)
            het_list.append(het_count)
            norm_list.append(normal_count)

        norm_in_df = pd.Series(norm_list)
        het_in_df = pd.Series(het_list)
        hom_in_df = pd.Series(hom_list)

        new_df['WT count'] = norm_in_df.values
        new_df['Het count'] = het_in_df.values
        new_df['Hom count'] = hom_in_df.values

        if csv_file[-4:] == '.csv':
            new_df.to_csv(path_or_buf=export_filepath, index=False)
        elif csv_file[-4:] == '.txt':
            new_df.to_csv(path_or_buf=export_filepath, index=False, sep='\t')
        print("Done")

       ########" IF MORE THAN 900 SAMPLES #######

    elif number_of_samples >= 900:
        print('\n More than 900 samples detected, Keeping the numbers of WT/het/hom and discarding Genotype matrix')
        avinput_df = pd.read_csv(avinput_path, header=None, sep="\t")

        avinput_df.drop(avinput_df.columns[range(0, 16)], axis=1, inplace=True)
        avinput_df.columns = headers_in_vcf
        for column in list(avinput_df.columns):
            if column == 'FORMAT':
                pass
            else:
                avinput_df[column].replace('0', '0/0', inplace=True)

        # SPLIT THE GENOTYPE VALUE from other infos in SAMPLES column
        sep = ':'
        for column in avinput_df.columns:
            value = avinput_df[column].str.split(sep, n=1, expand=True)
            avinput_df[column] = value[0]

        hom_list = []
        het_list = []
        norm_list = []

        for line in range(0, avinput_df.shape[0]):
            values_of_line = avinput_df.iloc[line].values
            hom_count = 0
            het_count = 0
            normal_count = 0
            for i in values_of_line:
                if i in ["1/1", "2/2", "1|1", "2|2", '1/2', '1|2', '2|1']:
                    hom_count = hom_count + 1
                elif i in ['0/1', '0/2', '0|1', '0|2', '1|0', '2|0']:
                    het_count = het_count + 1
                elif i == '0/0' or i == '0|0':
                    normal_count = normal_count + 1
            hom_list.append(hom_count)
            het_list.append(het_count)
            norm_list.append(normal_count)

        norm_in_df = pd.Series(norm_list)
        het_in_df = pd.Series(het_list)
        hom_in_df = pd.Series(hom_list)

        full_stuff['WT count'] = norm_in_df.values
        full_stuff['Het count'] = het_in_df.values
        full_stuff['Hom count'] = hom_in_df.values

        if csv_file[-4:] == '.csv':
            full_stuff.to_csv(path_or_buf=export_filepath, index=False)
        elif csv_file[-4:] == '.txt':
            full_stuff.to_csv(path_or_buf=export_filepath, index=False, sep='\t')

        print("Done")

    else:
        print('\n'+tmp_stmp()+'Only one sample found in vcf')
        pass

    # REMOVE TEMPORARY VCF AND DECOMPOSED VCF IF CREATED
    if vcf_path_abs.endswith('_tmp_decomposed.vcf'):
        os.remove(vcf_path_abs)
        tmp2 = vcf_path_abs.replace('_tmp_decomposed.vcf', '_tmp.vcf')
        os.remove(tmp2)
        print(tmp_stmp()+'Removed temporary files')
    elif vcf_path_abs.endswith('_tmp.vcf'):
        os.remove(vcf_path_abs)
        print(tmp_stmp() + 'Removed temporary files')
    elif vcf_path_abs.endswith('_decomposed.vcf'):
        os.remove(vcf_path_abs)
        print(tmp_stmp() + 'Removed temporary files')
    else:
        pass
    return


def check_proportion_of_transcripts_affected(full_table, acmg_db_path, ref_anno):
    transcript_file = os.path.join(os.path.join(process_path(acmg_db_path), 'acmg_db'),
                                   'transcripts_per_genes.txt')

    # Create dictionary with number of transcripts per genes with gene names as key
    with open(transcript_file, 'r') as file:
        transcript_dict = {}
        for line in file.readlines():
            if 'number_of_transcripts' in line:
                pass  # Skip header
            else:
                gene = line.split('\t')[0]
                number_temp = line.split('\t')[1]
                number = number_temp.replace('\n', '')
                transcript_dict[gene] = int(number)
    ordered_list = []

    gene_series = full_table['Gene.{}'.format(ref_anno)]

    # Make list for number of trancripts actually affected using NM reference sequence
    transcript_affected_list = []
    transcripts_affected = full_table['AAChange.{}'.format(ref_anno)]
    for line in transcripts_affected:
        count = line.count('NM_')
        transcript_affected_list.append(count)

    # Make list of number of transcript per gene present
    for gene in gene_series:
        try:
            sep = ';'
            gene_name_split = gene.split(sep, 1)[0]
            ordered_list.append(transcript_dict[gene_name_split])
        except KeyError:
            ordered_list.append('Nan')  # gene not found in dictionary

    proportional_tuple = list(zip(transcript_affected_list, ordered_list))

    # now do the thang and calculate proportion

    transcripts_affected_contrib = []

    for affect, total in proportional_tuple:
        try:

            if affect / total > 1:
                transcripts_affected_contrib.append(5)
            elif affect / total >= 0.66:
                transcripts_affected_contrib.append(5)
            elif affect / total > 0:
                transcripts_affected_contrib.append(-5)
            elif affect / total == 0:
                transcripts_affected_contrib.append(0)
        except TypeError:
            transcripts_affected_contrib.append(0)

    transcript_proportion_in_df = pd.Series(transcripts_affected_contrib)
    print(transcript_proportion_in_df)
    full_table['Transcripts_proportion_Coef'] = transcript_proportion_in_df.values

    return full_table


def pathway_analysis(full_table, path_db, ref_anno):
    try:
        global pathway_db_used
        pathway_db_used = path_db

        print('|| Using EnrichR to do pathway analysis with {}'.format(pathway_db_used))

        full_table = full_table[full_table["Probability_Path"] >= 0.80]
        gene_name_list = set(list(full_table["Gene.{}".format(ref_anno)]))

        if len(gene_name_list) <= 0:
            print('|| No Potentially Pathogenic variants detected, skipping EnrichR analysis...')
            return

        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        genes_str = '\n'.join(gene_name_list)
        description = 'Tapes gene list'
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        data_enrichr = json.loads(response.text)

        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        user_list_id = data_enrichr['userListId']
        gene_set_library = pathway_db_used
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        data_res = json.loads(response.text)
        data_res = data_res[pathway_db_used][0:11]

        global pathways_results
        pathways_results = pd.DataFrame.from_dict(data_res)
        column_names = ["Rank", "Name", "P-value", "Z-score", "Combined score", "Genes", "Adjusted p-value", "7", "8"]
        pathways_results.columns = column_names
        pathways_results = pathways_results[["Name", "P-value", "Z-score", "Combined score", "Genes", "Adjusted p-value"]]
        pathways_results.sort_values(by='Adjusted p-value', ascending=True)
    except TimeoutError:
        print('|| EnrichR : TimeoutError')
    return


def create_cancer_dataframe(full_stuff, disease, ref_anno):
    global cancer_df
    global disease_name
    disease_name = disease
    print("|| Creating {}-genes only spreadsheet".format(disease_name))
    cancer_df = full_stuff[full_stuff["Disease_description.{}".format(ref_anno)].str.contains(disease)]
    return


def individual_reports(full_stuff, ref_anno):
    global df_list
    full_stuff = full_stuff[full_stuff['Probability_Path'] >= 0.85]
    column_sample_list = []
    df_list = []

    for column in full_stuff.columns:
        if ('0/0' and "0/1") in list(full_stuff[column]):
            column_sample_list.append(column)

    for sample in column_sample_list:
        hom_for_samp = full_stuff.loc[full_stuff[sample] == '1/1']
        het_for_samp = full_stuff.loc[full_stuff[sample] == '0/1']
        all_for_samp = pd.concat([hom_for_samp, het_for_samp])
        all_for_samp = all_for_samp.sort_values(by='Probability_Path', ascending=False).head(5)
        all_for_samp_tuple = (sample, all_for_samp)
        df_list.append(all_for_samp_tuple)
    return


def export_to_formatted_xlsx(full_stuff, output_file, number_of_samples, ref_anno):
    print(tmp_stmp()+"Writing xlsx file...")
    number_of_rows = str(int(full_stuff.shape[0]) + 1)
    number_of_columns = int(full_stuff.shape[1]) + 1

    # IF THERE IS SAMPLE DATA
    # Find starting column of sample data
    try:
        print("|| Found sample data, proceeding accordingly...")
        column_sample_start = utility.xl_col_to_name(full_stuff.columns.get_loc('FORMAT') + 1)
        column_sample_stop = utility.xl_col_to_name(full_stuff.columns.get_loc('FORMAT') + 1 + number_of_samples)
    except KeyError:
        print("|| No Sample genotyping data, continuing without...")

    # Name of the column where our own annotation stops (NAME OF THE TOTAL COLUMN)
    col_total = utility.xl_col_to_name(number_of_columns - 2)
    col_prob = utility.xl_col_to_name(full_stuff.columns.get_loc('Probability_Path'))

    total_range = col_prob+"2:"+col_prob+str(full_stuff.shape[0]+1)

    # Start working with pandas and xlsxwriter
    writer = pd.ExcelWriter(output_file, engine='xlsxwriter')  # create the writer pandas to excel
    full_stuff.to_excel(writer, sheet_name='Data', index=False) # pd.dataframe to excel without index

    workbook = writer.book    # Create Workbook
    worksheet = writer.sheets['Data']   # Create Worksheet

    # DO conditional formating of the coef data and total


    worksheet.conditional_format(total_range, {'type': '2_color_scale',
                                             'min_color': "#89FE14",
                                             'max_color': "#F26243",
                                             'mid_type': "num"})

    col_name_to_hide = utility.xl_col_to_name(full_stuff.columns.get_loc('BA_1_contrib'))

    try:
        column_to_hide_with_sample = utility.xl_col_to_name(full_stuff.columns.get_loc('FORMAT'))
        worksheet.set_column('I:' + column_sample_stop, None, None, {'hidden': True})  # HIDE MOST COLUMNS
    except (UnboundLocalError, KeyError):
        worksheet.set_column('I:' + col_name_to_hide, None, None, {'hidden': True})  # HIDE MOST COLUMNS



    # ADD A 2ND SHEET FOR GRAPHS (TODO ADD MORE GRAPHS)
    wot_wot = full_stuff['Func.{}'.format(ref_anno)].value_counts()
    number_of_types = str(int(wot_wot.shape[0])+1)
    wot_wot.to_excel(writer, sheet_name='Charts')
    worksheet_chart = writer.sheets['Charts']

    # CHART 1
    chart_location = workbook.add_chart({'type': 'column'})
    chart_location.add_series({'values': '=Charts!$B$2:$B$'+number_of_types,
                               'categories': '=Charts!$A$2:$A$'+number_of_types,
                               'name':       '=Charts!$B$1'
                               })
    chart_location.set_y_axis({'name': 'Number of Variants'})
    chart_location.set_x_axis({'name': 'Genomic location'})
    worksheet_chart.insert_chart('D2', chart_location)
    chart_location.set_size({'width': 720, 'height': 576})

    # TRY CYTOBAND GRAPH
    try:
        path_cyto = full_stuff.loc[full_stuff['Prediction_ACMG_tapes'].isin(['Pathogenic', 'Likely Pathogenic'])]
        cyto_series = path_cyto['cytoBand'].value_counts()
        cyto_series.to_excel(writer, sheet_name='Charts', startrow=31, startcol=0, index=True)
        number_of_cyto = str(int(cyto_series.shape[0]) + 1)

        chart_cyto= workbook.add_chart({'type': 'column'})
        chart_cyto.add_series({'values': '=Charts!$B$33:$B$' + number_of_cyto,
                                   'categories': '=Charts!$A$33:$A$' + number_of_cyto,
                                   'name': '=Charts!$B$32'
                                   })
        chart_cyto.set_y_axis({'name': 'Number of Path/L_Path Variants'})
        chart_cyto.set_x_axis({'name': 'Genomic location'})
        worksheet_chart.insert_chart('D31', chart_cyto)
        chart_cyto.set_size({'width': 720, 'height': 576})
    except KeyError:
        pass

    # ADD A 3rd spreadsheet fo pathway analysipathways_resultss

    try:
        pathways_results.to_excel(writer, sheet_name=pathway_db_used, index=False)
        worksheet_pathways = writer.sheets[pathway_db_used]
    except NameError:
        pass

    # ADD 4th spreadsheet if cancer flag
    try:
        cancer_df.to_excel(writer, sheet_name=disease_name, index=False)
        worksheet_cancer = writer.sheets[disease_name]

        worksheet_cancer.conditional_format(total_range, {'type': '2_color_scale',
                                                   'min_color': "#89FE14",
                                                   'max_color': "#F26243",
                                                   'mid_type': "num"})
    except NameError:
        pass
    # ADD 5th spreadsheet to do individual report for each sample
    try:
        df_list
        col_minus_one = utility.xl_col_to_name(full_stuff.columns.get_loc('Probability_Path') - 1)
        print("|| Writing per sample 5 most pathigenic genes spreadsheet...")
        worksheet_sample = workbook.add_worksheet('Samples')
        writer.sheets['Samples'] = worksheet_sample
        bold = workbook.add_format({'bold': True})

        x_row = 1
        y_eow = 1
        for sample, df in df_list:
            worksheet_sample.write('A'+str(y_eow), sample, bold)
            df.to_excel(writer, sheet_name='Samples', startrow=x_row, startcol=0, index=False)
            x_row = x_row + 7
            y_eow = y_eow + 7
        worksheet_sample.set_column('Q:' + col_minus_one, None, None, {'hidden': True})
    except NameError:
        pass

    # ADD 5.5th spreadsheet to do per-gene report

    try:
        whole
        print("|| Writing per_genes spreadsheet...")
        worksheet_genes = workbook.add_worksheet('Per-Gene')
        writer.sheets['Per-Gene'] = worksheet_genes
        bold = workbook.add_format({'bold': True})
        x_row = 1
        y_eow = 1
        for tup in whole:
            if tup[6] is True:
                worksheet_genes.write('A'+str(y_eow), tup[0], bold)
                worksheet_genes.write('B' + str(y_eow), tup[3], bold)
                worksheet_genes.write('C' + str(y_eow), "FLAGS GENE", bold)
            elif int(tup[4]) > 250000:
                worksheet_genes.write('A'+str(y_eow), tup[0], bold)
                worksheet_genes.write('B' + str(y_eow), tup[3], bold)
                worksheet_genes.write('C' + str(y_eow), "LONG GENE", bold)
            elif tup[5] is True:
                worksheet_genes.write('A'+str(y_eow), tup[0], bold)
                worksheet_genes.write('B' + str(y_eow), tup[3], bold)
                worksheet_genes.write('C' + str(y_eow), "SAMPLE NUMBER WARNING", bold)
            else:
                worksheet_genes.write('A' + str(y_eow), tup[0], bold)
                worksheet_genes.write('B' + str(y_eow), tup[3], bold)
            tup[1].to_excel(writer, sheet_name='Per-Gene', startrow=x_row, startcol=0, index=False)
            x_row = x_row + tup[1].shape[0]+ 2
            y_eow = y_eow + tup[1].shape[0] +2
    except NameError:
        pass

    # ADD 6th spreadsheet to show user list of genes
    try:
        list_of_genes_df
        print('|| Writing user gene-list spreadsheet...')
        worksheet_list = workbook.add_worksheet('Gene List')
        writer.sheets['Gene List'] = worksheet_list
        list_of_genes_df.to_excel(writer, sheet_name='Gene List')
    except NameError:
        pass

    # ADD 7th spreadsheet to show keeg pathway reverse search
    try:
        kegg_df
        print("|| Writing spreadsheet for {}...".format(name_of_pathway_to_report))
        worksheet_kegg = workbook.add_worksheet(name_of_pathway_to_report)
        writer.sheets[name_of_pathway_to_report] = worksheet_kegg
        kegg_df.to_excel(writer, sheet_name=name_of_pathway_to_report)
    except (NameError, UnboundLocalError) as e:
        pass
    # SAVES THE EXCEL FILE
    writer.save()
    return


def check_online_annovar_dbs(annovar_path):
    global tapes_wd
    tapes_wd = os.path.dirname(os.path.realpath(__file__))
    os.chdir(annovar_path)
    print(tmp_stmp()+"Fetching ANNOVAR Alldb file")
    subprocess.run(['perl', 'annotate_variation.pl', '-buildver', 'hg19', '-downdb', '-webfrom', 'annovar', 'avdblist', './'])
    outfile_hg19 = os.path.join(annovar_path, "hg19_avdblist.txt")

    subprocess.run(
        ['perl', 'annotate_variation.pl', '-buildver', 'hg38', '-downdb', '-webfrom', 'annovar', 'avdblist', './'])
    outfile_hg38 = os.path.join(annovar_path, "hg38_avdblist.txt")
    global all_db
    all_db = []
    with open(outfile_hg19, 'r') as file:
        for line in file.readlines():
            line = line.split("\t")[0]
            if ".idx.gz" in line:
                pass
            else:
                #line = line.replace("hg19_", "")
                line = line.replace(".gz", "")
                line = line.replace(".txt", "")
                line = line.replace(".zip", "")
                all_db.append(line)
    with open(outfile_hg38, 'r') as file:
        for line in file.readlines():
            line = line.split("\t")[0]
            if ".idx.gz" in line:
                pass
            else:
                #line = line.replace("hg38_", "")
                line = line.replace(".gz", "")
                line = line.replace(".txt", "")
                line = line.replace(".zip", "")
                all_db.append(line)

    all_db = set(all_db)

    return


def scan_for_db(annovar_path):
    annovar_path_absolute = process_path(annovar_path)
    human_db_path = os.path.join(annovar_path_absolute, "humandb")
    pre_loaded_db = ['hg19_refGeneWithVer.txt', 'genometrax-sample-files-gff', 'GRCh37_MT_ensGene.txt',
                     'GRCh37_MT_ensGeneMrna.fa', 'hg19_example_db_gff3.txt',
                     'hg19_refGeneWithVerMrna.fa', 'hg19_MT_ensGeneMrna.fa',
                     'hg19_example_db_generic.txt', 'hg19_MT_ensGene.txt']
    if os.path.isdir(human_db_path):
        all_files = os.listdir(human_db_path)
        all_added_files = list(set(all_files) - set(pre_loaded_db))
        global hg_19_list
        global hg_38_list
        hg_19_list = []
        hg_38_list = []
        for file_name in all_added_files:
            if 'hg19' in file_name:
                file_name = file_name.replace(".gz", "")
                file_name = file_name.replace("hg19_", "")
                file_name = file_name.replace(".txt", "")
                file_name = file_name.replace(".idx", "")
                file_name = file_name.replace(".zip", "")
                hg_19_list.append(file_name)

            elif 'hg38' in file_name:
                file_name = file_name.replace(".gz", "")
                file_name = file_name.replace("hg38_", "")
                file_name = file_name.replace(".txt", "")
                file_name = file_name.replace(".idx", "")
                file_name = file_name.replace(".zip", "")
                hg_38_list.append(file_name)
        hg_19_list = set(hg_19_list)
        hg_38_list = set(hg_38_list)

    else:
        print("Can't find /humandb folder")
    return


def writing_db_to_file(annovar_path, acmg_db_path):
    acmg_db_path = process_path(acmg_db_path)
    os.chdir(tapes_wd)
    all_db_hg19_x = [x for x in all_db if "hg19" in x]
    all_db_hg38_x = [x for x in all_db if "hg38" in x]
    # No idea why it doesnt want to do inplace replace
    all_db_hg19 = []
    for item in all_db_hg19_x:
        item = item.replace("hg19_", "")
        all_db_hg19.append(item)

    all_db_hg38 = []
    for item in all_db_hg38_x:
        item = item.replace("hg38_", "")
        all_db_hg38.append(item)

    try:
        # See if downloaded db is in the list of all db
        list_for_dict_hg19 = []
        for item in all_db_hg19:
            if item in hg_19_list:
                a = (item, 'OK')
                list_for_dict_hg19.append(a)
            elif item not in hg_19_list:
                a = (item, 'MISSING')
                list_for_dict_hg19.append(a)
        # Add downloaded db even if it is not in the list of all db (somme are missing like genomicSuperdup)
        for item in hg_19_list:
            if item not in all_db_hg19:
                a = (item, 'OK')
                list_for_dict_hg19.append(a)

        list_for_dict_hg38 = []
        for item in all_db_hg38:
            if item in hg_38_list:
                a = (item, 'OK')
                list_for_dict_hg38.append(a)
            elif item not in hg_38_list:
                a = (item, 'MISSING')
                list_for_dict_hg38.append(a)

        for item in hg_38_list:
            if item not in all_db_hg38:
                a = (item, 'OK')
                list_for_dict_hg38.append(a)

        # If genomic superdup and cytobands or other region anno are not there, add them as missing

        for reg in reg_anno():
            if reg not in list_for_dict_hg19[0]:
                a = (reg, 'MISSING')
                list_for_dict_hg19.append(a)
            elif reg in list_for_dict_hg19[0]:
                a = (reg, 'OK')
                list_for_dict_hg19.append(a)
            if reg not in list_for_dict_hg38[0]:
                a = (reg, 'MISSING')
                list_for_dict_hg38.append(a)
            elif reg in list_for_dict_hg38[0]:
                a = (reg, 'OK')
                list_for_dict_hg38.append(a)



        # Join both dictionaries plus annovar path
        dico_19 = dict(sorted(list_for_dict_hg19))
        dico_38 = dict(sorted(list_for_dict_hg38))

        final_dict = {"hg19": dico_19,
                      "hg38": dico_38,
                      "annovar_path": {"annovar_path": annovar_path},
                      "acmg_db_path": {"acmg_db_path": acmg_db_path}

        }
        # DATABASES TO USE FOR VCF ANNOTATION
        vcf_json = os.path.join(os.getcwd(), "db_vcf.json")
        if not os.path.isfile(vcf_json):
            open(vcf_json, "a").write("{}")
        with open(vcf_json, 'r') as vj_there:
            already_there_vcf_dict = json.load(vj_there)

            # Build new dict without taking care of the present dict
            db_to_use_vcf_19 = []
            db_to_use_vcf_38 = []

            for key in final_dict['hg19']:
                if final_dict['hg19'][key] == 'OK':
                    if ".fa" in key:
                        pass
                    elif key == 'kgXref' or key == 'refGeneVersion' or key == 'avdblist':
                        pass
                    else:
                        a = (key, 'NO')
                        db_to_use_vcf_19.append(a)

            for key in final_dict['hg38']:
                if final_dict['hg38'][key] == 'OK':
                    if ".fa" in key:
                        pass
                    elif key == 'kgXref' or key == 'refGeneVersion' or key == 'avdblist':
                        pass
                    else:
                        a = (key, 'NO')
                        db_to_use_vcf_38.append(a)

            vcf_annotation_dict = {"hg19": dict(db_to_use_vcf_19),
                                   "hg38": dict(db_to_use_vcf_38)}
            try:
                for key in vcf_annotation_dict['hg19']:
                    if key in already_there_vcf_dict['hg19']:
                        pass
                    else:
                        already_there_vcf_dict['hg19'][key] = 'NO'

                for key in vcf_annotation_dict['hg38']:
                    if key in already_there_vcf_dict['hg38']:
                        pass
                    else:
                        already_there_vcf_dict['hg38'][key] = 'NO'

                with open(vcf_json, "w") as vj:
                    json.dump(already_there_vcf_dict, vj, indent=1)
            except KeyError:
                with open(vcf_json, "w") as vj:
                    json.dump(vcf_annotation_dict, vj, indent=1)

        # DATABASE DL CONFIG FILE
        config_json = os.path.join(os.getcwd(), "db_config.json")
        # Create config file if it does not exist
        if not os.path.isfile(config_json):
            open(config_json, "a").write("{}")

        # Write (and overwrite) everything to file
        with open(config_json, "w") as dj:
            json.dump(final_dict, dj, indent=1)
        print("ANNOVAR humandb folder location : ")
        print(annovar_path)
        print('hg19 databases present in /humandb :')
        #  Nicer output than just print list
        foolist = [k for k,v in final_dict['hg19'].items() if v == 'OK']
        if len(foolist) <= 3:
            for x in foolist: print(x)
        else:
            for a, b, c in zip(foolist[::3], foolist[1::3], foolist[2::3]):
                print('{:<30}{:<30}{:<}'.format(a, b, c))

        print()
        print('hg38 databases present in /humandb :')
        foolist = [k for k, v in final_dict['hg38'].items() if v == 'OK']
        if len(foolist) <= 3:
            for x in foolist: print(x)
        else:
            for a, b, c in zip(foolist[::3], foolist[1::3], foolist[2::3]):
                print('{:<30}{:<30}{:<}'.format(a, b, c))
        print()
        print(tmp_stmp()+'Saved db_config.json in {}'.format(os.getcwd()))
    except NameError:
        pass
    return


def check_for_acmg_databases(acmg_db_path):
    acmg_db_path = process_path(acmg_db_path)
    list_of_db = [
        'BS2_hom_het.hg19',
        'BS2_hom_het.hg38',
        'BS2_rec_dom_ad.txt',
        'clinvar.PS1.hg19.txt',
        'clinvar.PS1.hg38.txt',
        'final_PVS1.txt',
        'PM1_domains_with_benigns.hg19',
        'PM1_domains_with_benigns.hg38',
        'PP2_BP1.txt',
        'PS4.variants.hg19',
        'PS4.variants.hg38',
        'REK_canon.hg19',
        'REK_canon.hg38',
        'repeat_dict.hg19.gz',
        'repeat_dict.hg38.gz',
        'transcripts_per_genes.txt'
    ]
    if os.path.isdir(acmg_db_path):
        for db_file in list_of_db:
            if os.path.isfile(os.path.join(acmg_db_path, db_file)):
                pass
            else:
                print('{} is missing, please download it from https://github.com/a-xavier/tapes'.format(db_file))
                sys.exit(1)
        print(tmp_stmp()+'All acmg_db files found')
    else:
        print("|| 'acmg_db' folder not found at {}, exiting...".format(acmg_db_path))
        sys.exit(1)


def find_list_of_genes(full_stuff, list_of_genes, ref_anno):
    global list_of_genes_df
    list_of_genes_df = full_stuff[full_stuff['Gene.{}'.format(ref_anno)].isin(list_of_genes)]
    return


def reg_anno():
    return ['genomicSuperDups', 'wgRna', 'cytoBand', 'tfbsConsSites', 'targetScanS', 'dgvMerged', 'gwasCatalog', 'phastConsElements46way']


def db_type_bank(): # as of 21-11-1018
    gene_base = [
        'refGene', 'knownGene', 'ensGene'
    ]

    filter_base = [
        'abraom',
        'AFR.sites.2012_04', 'ALL.sites.2010_11', 'ALL.sites.2011_05', 'ALL.sites.2012_02', 'ALL.sites.2012_04',
        'AMR.sites.2012_04', 'ASN.sites.2012_04', 'avgwas_20150121',
        'avsift',
        'avsnp138',
        'avsnp142',
        'avsnp144',
        'avsnp147',
        'avsnp150',
        'cadd13gt10',
        'cadd13gt20',
        'cadd13',
        'caddgt10',
        'caddgt20',
        'caddindel',
        'cadd',
        'cg46',
        'cg69',
        'clinvar_20130905',
        'clinvar_20131105',
        'clinvar_20140211',
        'clinvar_20140303',
        'clinvar_20140702',
        'clinvar_20140902',
        'clinvar_20140929',
        'clinvar_20150330',
        'clinvar_20150629',
        'clinvar_20151201',
        'clinvar_20160302',
        'clinvar_20161128',
        'clinvar_20170130',
        'clinvar_20170501',
        'clinvar_20170905',
        'clinvar_20180603',
        'cosmic64',
        'cosmic65',
        'cosmic67',
        'cosmic67wgs',
        'cosmic68',
        'cosmic68wgs',
        'cosmic70',
        'dann',
        'dbnsfp30a',
        'dbnsfp31a_interpro',
        'dbnsfp33a',
        'dbscsnv11',
        'eigen',
        'esp5400_aa',
        'esp5400_all',
        'esp5400_ea',
        'esp6500_aa',
        'esp6500_all',
        'esp6500_ea',
        'esp6500si_aa',
        'esp6500si_all',
        'esp6500si_ea',
        'esp6500siv2_aa',
        'esp6500siv2_all',
        'esp6500siv2_ea',
        'EUR.sites.2012_04',
        'exac01',
        'exac02',
        'exac03nonpsych',
        'exac03nontcga',
        'exac03',
        'fathmm',
        'gerp++elem',
        'gerp++gt2',
        'gerp++',
        'gme',
        'gnomad_exome',
        'gnomad_genome',
        'gwava',
        'hrcr1',
        'icgc21',
        'intervar_20170202',
        'intervar_20180118',
        'kaviar_20150923',
        'kgXref',
        'ljb23_fathmm',
        'ljb23_gerp++',
        'ljb23_lrt',
        'ljb23_ma',
        'ljb23_metalr',
        'ljb23_metasvm',
        'ljb23_mt',
        'ljb23_phylop',
        'ljb23_pp2hdiv',
        'ljb23_pp2hvar',
        'ljb23_sift',
        'ljb23_siphy',
        'ljb26_all',
        'ljb26_cadd',
        'ljb26_fathmm',
        'ljb26_gerp++',
        'ljb26_lrt',
        'ljb26_ma',
        'ljb26_metalr',
        'ljb26_metasvm',
        'ljb26_mt',
        'ljb26_phylop100way_vertebrate',
        'ljb26_phylop46way_placental',
        'ljb26_pp2hdiv',
        'ljb26_pp2hvar',
        'ljb26_sift',
        'ljb26_siphy',
        'ljb26_vest',
        'ljb2_fathmm',
        'ljb2_gerp++',
        'ljb2_lrt',
        'ljb2_ma',
        'ljb2_mt',
        'ljb2_phylop',
        'ljb2_pp2hdiv',
        'ljb2_pp2hvar',
        'ljb2_sift',
        'ljb2_siphy',
        'ljb_gerp++',
        'ljb_lrt',
        'ljb_mt',
        'ljb_phylop',
        'ljb_pp2',
        'ljb_sift',
        'mcap',
        'mitimpact24',
        'mitimpact2',
        'nci60',
        'popfreq_all_20150413',
        'popfreq_max_20150413',
        'refLink',
        'regsnpintron',
        'revel',
        'snp129NonFlagged',
        'snp129',
        'snp130NonFlagged',
        'snp130',
        'snp131NonFlagged',
        'snp131',
        'snp132NonFlagged',
        'snp132',
        'snp135NonFlagged',
        'snp135',
        'snp137NonFlagged',
        'snp137',
        'snp138NonFlagged',
        'snp138',
        'snp142',
        '1000g2010nov',
        '1000g2011may',
        '1000g2012apr',
        '1000g2012feb',
        '1000g2014aug',
        '1000g2014oct',
        '1000g2014sep',
        '1000g2015aug',
        'spidex'
    ]

    region_base = reg_anno()

    gene_based_tup = []
    for db in gene_base:
        gene_based_tup.append((db, 'gx'))

    filter_based_tup = []
    for db in filter_base:
        filter_based_tup.append((db, 'f'))

    region_based_tup = []
    for db in region_base:
        region_based_tup.append((db, 'r'))

    all_list_db_type = gene_based_tup + filter_based_tup + region_based_tup
    all_list_db_type_dict = dict(all_list_db_type)
    return all_list_db_type_dict


def replace_dash_by_dots(full_stuff):
    full_stuff['Ref'] = full_stuff['Ref'].replace("-", ".")
    full_stuff['Alt'] = full_stuff['Alt'].replace("-", ".")
    return full_stuff


def counting_number_of_samples(full_stuff):
    try:
        samps = []
        samp_df = full_stuff[['WT count', 'Het count', 'Hom count']]
        if samp_df.shape[0] > 20:
            for line in range(0, 20):
                values_of_line = samp_df.iloc[line].values
                values_of_line = list(map(int, values_of_line))
                samps.append(sum(values_of_line))
            num_of_samps = max(samps)
        else:
            for line in range(0, 1):
                values_of_line = samp_df.iloc[line].values
                samps.append(sum(values_of_line))
            num_of_samps = max(samps)

        return num_of_samps
    except KeyError:
        print('|| No sample data found')
        num_of_samps = 1
        return num_of_samps


def open_file_in_os(path):
    subprocess.call(('xdg-open', path))
    return


def kegg_inverted_patway(name_of_pathway, full_stuff, ref_anno, acmg_db_path):
    with open(os.path.join(acmg_db_path, 'kegg_dict.json'), "r") as dj:
        kegg_dict = json.load(dj)
    list_of_keys = [x for x in kegg_dict]
    global name_of_pathway_to_report
    name_of_pathway_to_report = name_of_pathway
    lower = name_of_pathway.lower()
    try:
        lower = get_close_matches(lower, list_of_keys, cutoff=0.7)[0]
        name_of_pathway_to_report = lower
        gene_kegg_list = kegg_dict[lower]
        global kegg_df
        kegg_df = full_stuff.loc[full_stuff['Gene.{}'.format(ref_anno)].isin(gene_kegg_list)]
    except KeyError:
        print(tmp_stmp()+'{} is not a usable KEGG pathway, Skipping...'.format(name_of_pathway_to_report))
    return


def chr_num():
    return ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13','14','15', '16', '17', '18', '19', '20',
            '21', '22', 'X', 'Y', 'M']


def calculate_probability(dataframe_acmg):
    PVS = dataframe_acmg['PVS1_contrib']
    PS = dataframe_acmg['PS1_contrib'] + dataframe_acmg['PS2_contrib_merged'] + dataframe_acmg['PS3_contrib'] + dataframe_acmg['PS4_contrib']
    PM = dataframe_acmg['PM1_contrib'] + dataframe_acmg['PM2_contrib'] + dataframe_acmg['PM4_contrib'] + dataframe_acmg['PM5_contrib']
    PP = dataframe_acmg['PP2_contrib'] + dataframe_acmg['PP3_contrib'] + dataframe_acmg['PP5_contrib']
    BS = dataframe_acmg['BS1_contrib'] + dataframe_acmg['BS2_contrib'] + dataframe_acmg['BS3_contrib']
    BP = dataframe_acmg['BP1_contrib'] + dataframe_acmg['BP3_contrib'] + dataframe_acmg['BP4_contrib'] + dataframe_acmg['BP6_contrib'] + dataframe_acmg['BP7_contrib']
    BA1 = dataframe_acmg['BA_1_contrib']

    pc = 0.10 # Prior C
    X = 2
    probability_list = []
    for pvs, ps, pm, pp, bs, bp, ba1 in zip(PVS, PS, PM, PP, BS, BP, BA1):
        if ba1 != 1:
            odds_path = 350**((pp/(X**3)) + (pm/X**2) + (ps/X) + (pvs/1) - (bp/X**3) - (bs/X)) # ODDs pathogenicity equation

            proba = (odds_path * pc) / (((odds_path - 1) * pc) +1)  # Probability equation

            probability_list.append(float("%.4f" % proba))  # 4 significant digits
        elif ba1 == 1:
            probability_list.append(0)

    prediction_list= []
    for pvs, ps, pm, pp, bs, bp, ba1 in zip(PVS, PS, PM, PP, BS, BP, BA1):

        if ba1 == 1:
            prediction_list.append('Benign auto')
        elif bs >= 2:
            prediction_list.append('Benign')
        elif (bs == 1 and bp >= 1) or (bp >=2):
            prediction_list.append('Likely Benign')
        elif (pvs == 1 and (ps >= 1 or pm >= 2 or (pm == 1 and pp == 1) or pp >= 2)) or (ps >= 2) or (ps == 1 and (pm >= 3 or (pm == 2 and pp >=2) or (pm == 1 and pp >= 4))):
            prediction_list.append('Pathogenic')
        elif (pvs == 1 and pm == 1) or (ps == 1 and pm >= 1) or (ps == 1 and pp >= 2) or (pm >= 3) or (pm >= 2 and pp >= 2) or (pm == 1 and pp >= 4):
            prediction_list.append('Likely Pathogenic')
        else:
            prediction_list.append('VUS')
    return (probability_list, prediction_list)


def list_to_dataframe(list_of_lists, list_name_of_columns):
    scale = len(list_of_lists[0])
    for item in list_of_lists:
        if len(item) != scale:
            print('Inconsistencies in contribution list, exiting')
            sys.exit(1)
    data = {}
    for list_user, name in zip(list_of_lists, list_name_of_columns):
        data[name] = list_user
    new_dataframe = pd.DataFrame.from_dict(data)
    return new_dataframe


def check_PVS1_criteria(full_stuff, PVS1_list, splice_list, ref_anno, rek_dict):
    exonic_type = full_stuff['ExonicFunc.{}'.format(ref_anno)]
    func_type = full_stuff['Func.{}'.format(ref_anno)]
    gene_series = full_stuff['Gene.{}'.format(ref_anno)]
    aa_change_series = full_stuff['AAChange.{}'.format(ref_anno)]
    start_series = full_stuff['Start']

    PVS1_contrib = []

    for exonic, func, gene, splice, aa, start in zip(exonic_type, func_type, gene_series, splice_list, aa_change_series, start_series):
        if gene.split(';')[0] in PVS1_list or gene in PVS1_list:  #  In case there is a ; in the gene symbol, check for first element only
            if exonic == 'stopgain' or exonic == 'frameshift deletion' or exonic == 'frameshift insertion' or exonic == \
                    'frameshift_deletion' or exonic == 'frameshift_insertion':

                #TODO INCLUDE OR REMOVE VARIANTS CLOSE TO END OF THE GENE
                if ref_anno == 'refGene':
                    #exon_nums = [x.replace('exon', '') for x in aa.split(':') if 'exon' in x]
                    trans_num = [x for x in aa.split(':') if 'NM_' in x]
                elif ref_anno == 'ensGene':
                    #exon_nums = [x.replace('exon', '') for x in aa.split(':') if 'exon' in x]
                    trans_num = [x for x in aa.split(':') if 'ENST' in x]
                elif ref_anno == 'knownGene':
                    #exon_nums = [x.replace('exon', '') for x in aa.split(':') if 'exon' in x]
                    trans_num = [x for x in aa.split(':') if 'uc0' in x]
                for transcript in trans_num:
                    try:
                        stop = rek_dict[transcript]
                        break
                    except KeyError:
                        pass
                try:
                    if int(stop) - int(start) <= 110:  # To offset untranslated regions
                        PVS1_contrib.append(0)
                    else:
                        PVS1_contrib.append(1)
                except UnboundLocalError:
                    PVS1_contrib.append(0)
                #PVS1_contrib.append(1)

            elif func == 'splicing' and splice >= 0.6:
                PVS1_contrib.append(1)
            else:
                PVS1_contrib.append(0)
        else:
            PVS1_contrib.append(0)
    return PVS1_contrib


# TODO ADD -50bp of last exon to not do anything if null variant
def check_PS1_criteria(full_stuff, PS1_dict,ref_anno):
    chr_series = full_stuff['Chr']
    ref_series = full_stuff['Ref']
    alt_series = full_stuff['Alt']
    start_series = full_stuff['Start']
    end_series = full_stuff['End']
    aa_change_series = full_stuff['AAChange.{}'.format(ref_anno)]
    exonic_func_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]
    PS1_contrib = []

    for chr, ref, alt, start, stop, aa_change, exo_func in zip(chr_series, ref_series, alt_series, start_series, end_series, aa_change_series, exonic_func_series):
        if "frameshift" not in exo_func:
            if len(aa_change) >= 11 and (len(ref) and len(alt)) == 1 and (ref and alt) != ('-' or '.'):# To avoid catching indels with aa changes
                isolated = aa_change.split(":p.")[-1] # Take the last segment of aachange just after 'p.'
                first = isolated[:1]
                last = isolated[-1:]
                try:
                    if first == PS1_dict[(chr, start, stop)]["AA_ref"] and last == PS1_dict[(chr, start, stop)]["AA_alt"]:
                        PS1_contrib.append(1)
                    else:
                        PS1_contrib.append(0)
                except KeyError:
                    PS1_contrib.append(0)
            else:   #  IF not frameshift and no aa change data available
                PS1_contrib.append(0)
        elif "frameshift" in exo_func:  # Exclude indels and nonframeshift substitution
            PS1_contrib.append(0)
    return PS1_contrib


def check_PS2_criteria(full_stuff, whole_dict):
    PS2_total_contrib = []
    try:
        for key in whole_dict: # for each family in the dictionary
            mother_series = full_stuff[whole_dict[key]['m']]
            father_series = full_stuff[whole_dict[key]['f']]
            offspring_series = full_stuff[whole_dict[key]['o']]
            PS2_contrib_per_fam = []
            for mom, dad, child in zip(mother_series, father_series, offspring_series):
                if mom == '0/0' and dad == '0/0' and child == ('0/1' or '1/1' or '0/2' or '2/2'):
                    PS2_contrib_per_fam.append(1)
                else:
                    PS2_contrib_per_fam.append(0)
            b = (key, PS2_contrib_per_fam)
            PS2_total_contrib.append(b)
    except KeyError:
        print(tmp_stmp()+
              "Could not find specified sample. Make sure that the input file includes genotyping data"
              " and that sample IDs match in the trio file. Exiting...")
        sys.exit(1)
    return PS2_total_contrib


def check_PS3(full_stuff):
    review_series = full_stuff['CLNREVSTAT']
    sig_series = full_stuff['CLNSIG']

    PS3_contrib = []
    for rev, sig in zip(review_series, sig_series):
        if ('pathogenic' in sig.lower() or 'drug_response' in sig.lower()) and\
                ('practice_guideline' in rev or 'reviewed_by_expert_panel' in rev):
            PS3_contrib.append(1)
        else:
            PS3_contrib.append(0)

    return PS3_contrib


def or_check(n_1, n_2, n_3, n_4):
    try:
        #OR =(n_1 * n_4) / (n_3 * n_2)
        fish_test = fisher_exact([[n_1, n_2], [n_3, n_4]])
        OR = fish_test[0]
        OR_sig = float("%.2f" % OR)
        SE_lnOR = sqrt((1 / n_1) + (1 / n_2) + (1 / n_3) + (1 / n_4))
        confidence_interval = (exp(log1p(OR) - (1.96 * SE_lnOR)), exp(log1p(OR) + (1.96 * SE_lnOR)))
        # Down to 2 decimal digits
        confidence_interval_sig = (float("%.2f" % confidence_interval[0]), float("%.2f" % confidence_interval[1]))
        p_value = fish_test[1]
        p_value = float("%.4f" % p_value)
    except ZeroDivisionError:
        OR_sig = float('nan')
        confidence_interval_sig = float('nan')
        p_value = float('nan')
    return [OR_sig, confidence_interval_sig, p_value]


def check_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df, number_of_samples):
    if ('WT count' and 'Het count' and 'Hom count') in full_stuff.columns and number_of_samples >= 15:
        print(tmp_stmp()+"Calculating Odd ratios for mutation enrichment", end=' ', flush=True)
        PS4_contrib = []
        o_ratio_list = []
        conf_inter_list = []
        p_value_list = []

        hom_col = pd.to_numeric(full_stuff['Hom count'], errors='coerce')
        het_col = pd.to_numeric(full_stuff['Het count'], errors='coerce')
        wt_col = pd.to_numeric(full_stuff['WT count'], errors='coerce')

        freq_exo_series.fillna(0, inplace=True)
        freq_geno_series.fillna(0, inplace=True)

        fish = pd.concat([hom_col, het_col, wt_col , freq_geno_series, freq_exo_series], axis=1, sort=False)

        tuple_for_fish = (list(fish.itertuples(index=False, name=None)))

        # print(full_stuff['Het count'], full_stuff['Hom count'], full_stuff['WT count'], freq_col)
        for het, hom, wt, freq_geno, freq_exo in tuple_for_fish:

            # print('het: ', het, 'hom :', hom, 'wt :', wt, 'freq : ', freq)

            #########################################################
            #                   # With variant    #  Without variant
            #########################################################
            # cancer(our data)   #   n_1           #       n_3
            ##########################################################
            # no cancer(gnomad)  #   n_2           #      n_4
            ##########################################################

            # OR = BAD OUTCOME (Our data) group = n_1 control = n_2 GOOD OUTCOME (gnomad) group = n_3 control = n_4
            # Calculation OR = (n_1/n_3)/(n_2/n_4) = (n_1 * n_4)/(n_3 * n_2)
            # SE{ln(OR)} = sqrt( (1/n_1) + (1/n_2) + (1/n_3) + (1/n_4)
            #  95% CI = range( log1p(OR) - (1.96 * SE{ln(OR)}), log1p(OR) + (1.96 * SE{ln(OR)}))
            n_1 = int(het + hom)

            n_2 = int(wt)

            freq = max(freq_exo, freq_geno)

            # Frequency is MAF so allele frequency rather than individuals so total allele divided
            # By 2 represents total individuals and since most are heterzygours n_3 is kept as is
            # ONLY WORKS IF WE ASSUME VARIANT IS HETEROZYGOUS
            if freq <= 0.001 and freq != 0 and n_1 >= 2:
                sci = '{:.1e}'.format(freq)

                power = int(str(sci.split("e-")[1]))

                n_3 = ceil(float(sci.split("e-")[0]))

                n_4 = ((10 ** power) / 2) - n_3
                fish_test = or_check(n_1, n_2, n_3, n_4)
                odd = fish_test[0]
                o_ratio_list.append(odd)
                ci = fish_test[1]
                conf_inter_list.append(ci)
                p_val = fish_test[2]
                p_value_list.append(p_val)

                if odd > 20 and not (ci[0] <= 1 <= ci[1]) and p_val <= 0.01:
                    PS4_contrib.append(1)
                else:
                    PS4_contrib.append(0)

            elif 0.001 < freq <= 0.05 and n_1 >= 2:
                sci = '{:.1e}'.format(freq)

                power = int(str(sci.split("e-")[1]))
                n_3 = (ceil(float(sci.split("e-")[0]))) * 100
                n_4 = (((10 ** power) / 2) * 100) - n_3

                fish_test = or_check(n_1, n_2, n_3, n_4)
                odd = fish_test[0]
                o_ratio_list.append(odd)
                ci = fish_test[1]
                conf_inter_list.append(ci)
                p_val = fish_test[2]
                p_value_list.append(p_val)
                if odd > 20 and not (ci[0] <= 1 <= ci[1]) and p_val <= 0.01:
                    PS4_contrib.append(1)
                else:
                    PS4_contrib.append(0)
            elif freq > 0.01:
                o_ratio_list.append('.')
                conf_inter_list.append('.')
                p_value_list.append('.')
                PS4_contrib.append(0)
            else:
                o_ratio_list.append('.')
                conf_inter_list.append('.')
                p_value_list.append('.')
                PS4_contrib.append(0)
        print("Done")
        return PS4_contrib, o_ratio_list, conf_inter_list, p_value_list

    elif ('WT count' and 'Het count' and 'Hom count') not in full_stuff.columns or number_of_samples < 15:
        PS4_contrib = []
        snp_col = [s for s in full_stuff.columns if 'snp' in s][0]
        try:
            rs_number_series = full_stuff[snp_col]
            list_of_rsid = PS4_df['SNPID'].tolist()
            for rsid in rs_number_series:
                if rsid in list_of_rsid:
                    PS4_contrib.append(1)
                else:
                    PS4_contrib.append(0)
        except KeyError:
            print('|| Could not find SNP id column to assign PS4')
            sys.exit(1)
        return PS4_contrib


def check_PM1_criteria(full_stuff, benign_domain_dict, ref_anno):
    PM1_list = []
    domain_series = full_stuff['Interpro_domain']
    exonic_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]

    for domain, exonic in zip(domain_series, exonic_series):
        if domain in benign_domain_dict or domain == '.':  # If the domain is not present or domain is in the benign dic
            PM1_list.append(0)
        elif domain not in benign_domain_dict and 'nonsynonymous' in exonic:
        #elif domain not in benig_domain_dict and ('stop' in exonic or 'frameshift' in exonic or 'nonsynonymous' in exonic):
        # Not sure why but intervar excludes null variants ? look for rationale
            PM1_list.append(1)
        else:
            PM1_list.append(0)
    return PM1_list


def check_PM2(full_stuff, freq_exo_series, freq_geno_series, ref_anno):
    pli_series = pd.to_numeric(full_stuff['pLi.{}'.format(ref_anno)], errors='coerce')
    prec_series = pd.to_numeric(full_stuff['pRec.{}'.format(ref_anno)], errors='coerce')
    PM2_contrib = []

    #  USE PLI SCORE FOR RECESSIVE DOMINANT CALLING
    #  Dominant / Haploinsufficient if pli score > 0.85
    #  Recessive if prec score > 0.85
    #TODO INCLUDE IN PM2 DOMINANT VARIANT WITH MAF < 10e-6 (Probably means only one het individual in the cohort)
    for pli, prec, exo, geno in zip(pli_series, prec_series, freq_exo_series, freq_geno_series):
        #print('exo ' + str(exo) + " geno " + str(geno) + " pli "+ str(pli))
        if (pd.isna(exo) and pd.isna(geno) and pli >= 0.85) or (exo <= 0.000001 and geno <=0.000001 and pli >= 0.85) or \
                (exo <= 0.000001 and pd.isna(geno) and pli >= 0.85) or (pd.isna(exo) and geno <=0.000001 and pli >= 0.85):  # Dominant/Haploinsufficient
            PM2_contrib.append(1)
        elif ((exo and geno <= 0.005 and pli) or (pd.isna(exo) and pd.isna(geno)) or (exo <= 0.005 and pd.isna(geno)) or\
                (geno <= 0.005 and pd.isna(exo))) and (pli <= 0.85 or pd.isna(pli)):  # Recessive (no pLi requirement)

            PM2_contrib.append(1)
        else:

            PM2_contrib.append(0)
    return PM2_contrib


def check_PM4_and_BP3_repeat_region(repeat_dict, full_stuff, ref_anno):
    PM4_contrib = []
    BP3_contrib = []
    chr_series = full_stuff['Chr']
    snp_start_series = full_stuff['Start']
    in_frame_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]
    interpro_series = full_stuff['Interpro_domain']
    for chrom, snp_start, in_frame, domain in zip(chr_series, snp_start_series, in_frame_series, interpro_series):
        if 'nonframeshift' in in_frame:
            list_of_all_first_elements = [item[0] for item in repeat_dict[
                'chr' + str(chrom)]]  # Create a list of all the starting point of repeat regions
            list_of_all_second_elements = [item[1] for item in repeat_dict[
                'chr' + str(chrom)]]  # Create a list of all the ending point of repeat regions

            closest1 = min(list_of_all_first_elements,
                           key=lambda x: abs(x - snp_start))  # Find the closest starting point
            closest2 = min(list_of_all_second_elements,
                           key=lambda x: abs(x - snp_start))  # Find the closest ending point

            # Find the index of the closest starting and ending point (might be the same or different)
            index_of_first = list_of_all_first_elements.index(closest1)
            index_of_second = list_of_all_second_elements.index(closest2)
            # Create range of closest repeat regions
            range_1 = range(repeat_dict['chr' + str(chrom)][index_of_first][0],
                            repeat_dict['chr' + str(chrom)][index_of_first][1])
            range_2 = range(repeat_dict['chr' + str(chrom)][index_of_second][0],
                            repeat_dict['chr' + str(chrom)][index_of_second][1])
            if snp_start in range_1 or snp_start in range_2:  # If in repeat region PM4 = 0 CHeck BP3 domain condition
                PM4_contrib.append(0)
                if domain == '.' or len(domain) <= 1:  # Check if known domain (either dot or any one character NA string
                    BP3_contrib.append(1)
                else:
                    BP3_contrib.append(0)
            else:                                             # If not in repeat region PM4 = 1 and BP3 = 0
                PM4_contrib.append(1)
                BP3_contrib.append(0)

        elif 'stoploss' == in_frame:
            list_of_all_first_elements = [item[0] for item in repeat_dict[
                'chr' + str(chrom)]]  # Create a list of all the starting point of repeat regions
            list_of_all_second_elements = [item[1] for item in repeat_dict[
                'chr' + str(chrom)]]  # Create a list of all the ending point of repeat regions

            closest1 = min(list_of_all_first_elements,
                           key=lambda x: abs(x - snp_start))  # Find the closest starting point
            closest2 = min(list_of_all_second_elements,
                           key=lambda x: abs(x - snp_start))  # Find the closest ending point

            # Find the index of the closest starting and ending point (might be the same or different)
            index_of_first = list_of_all_first_elements.index(closest1)
            index_of_second = list_of_all_second_elements.index(closest2)
            # Create range of closest repeat regions
            range_1 = range(repeat_dict['chr' + str(chrom)][index_of_first][0],
                            repeat_dict['chr' + str(chrom)][index_of_first][1])
            range_2 = range(repeat_dict['chr' + str(chrom)][index_of_second][0],
                            repeat_dict['chr' + str(chrom)][index_of_second][1])
            if snp_start in range_1 or snp_start in range_2:  # If in repeat region PM4 = 0 and BP3 = 0 because stoploss not nonframeshift
                PM4_contrib.append(0)
                BP3_contrib.append(0)
            else:                                             # Else PM4 = 1 and BP3 = 0 because stoploss not nonframeshift
                PM4_contrib.append(1)
                BP3_contrib.append(0)
            # IF NOT IN REPEAT REGION and is a nonframeshift variant or a stoploss PM4=1
        else: # If not a nonframeshift or stoploss
            PM4_contrib.append(0)
            BP3_contrib.append(0)

    return PM4_contrib, BP3_contrib


def check_PM5_criteria(full_stuff, PS1_dict, ref_anno):
    chr_series = full_stuff['Chr']
    start_series = full_stuff['Start']
    end_series = full_stuff['End']
    aa_change_series = full_stuff['AAChange.{}'.format(ref_anno)]
    exonic_func_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]

    PM5_contrib = []

    '''if len(aa_change) >= 11 and (len(ref) and len(alt)) == 1 and (ref and alt) != (
            '-' or '.'):  # To avoid catching indels with aa changes
        isolated = aa_change.split(":p.")[-1]  # Take the last segment of aachange just after 'p.' '''

    for chr, start, stop, aa_change, exo_func in zip(chr_series, start_series, end_series, aa_change_series, exonic_func_series):
        if 'frameshift' not in exo_func:
            if len(aa_change_series) >= 10 and ':p.' in aa_change_series:
                isolated = aa_change.split(":p.")[-1]  # Take the last segment of aachange just after 'p.'
                first = isolated[:1]
                last = isolated[-1:]
                try:
                    # Check if last is different so it doesn't add up with PS1
                    if first == PS1_dict[(chr, start, stop)]["AA_ref"] and last != PS1_dict[(chr, start, stop)]["AA_alt"]:
                        PM5_contrib.append(1)
                    else:
                        PM5_contrib.append(0)
                except KeyError:
                    PM5_contrib.append(0)
            else:  # IF not frameshift but no aa change data available
                PM5_contrib.append(0)
        elif 'frameshift' in exo_func:  # Exclude indels and non/frameshift substitution
            PM5_contrib.append(0)
    return PM5_contrib


def check_PP2(full_stuff, PP2_list, ref_anno):
    gene_series = full_stuff['Gene.{}'.format(ref_anno)]
    type_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]

    PP2_contrib = []
    for gene, type_snv in zip(gene_series, type_series):
        if gene in PP2_list and type_snv == 'nonsynonymous SNV':
            PP2_contrib.append(1)
        else:
            PP2_contrib.append(0)
    return PP2_contrib


def check_PP3_and_BP4_dbnsfp(full_stuff):

    # STARTING TO USE THE IN SILICO PREDICTION TOOLS

    # LET'S START WITH SIFT_score

    sift_string = full_stuff["SIFT_score"]
    sift_mixed = pd.to_numeric(sift_string, errors='coerce')
    sift_contrib = []

    for s_contrib in sift_mixed:
        if s_contrib <= 0.05:
            sift_contrib.append(1)
        elif s_contrib > 0.05:
            sift_contrib.append(-1)
        else:
            sift_contrib.append(0)

    # THEN LRT_pred

    lrt_string = full_stuff['LRT_pred']
    lrt_contrib = []

    for lrt_pred in lrt_string:
        if lrt_pred == "D":
            lrt_contrib.append(1)
        elif lrt_pred == "N":
            lrt_contrib.append(-1)
        else:
            lrt_contrib.append(0)

    #  THEN MutationTaster_pred

    mut_taster_string = full_stuff['MutationTaster_pred']
    mut_taster_contrib = []

    for mut_taster_pred in mut_taster_string:
        if mut_taster_pred == "A":
            mut_taster_contrib.append(1)
        elif mut_taster_pred == "D":
            mut_taster_contrib.append(1)
        elif mut_taster_pred == "N":
            mut_taster_contrib.append(-1)
        elif mut_taster_pred == "P":
            mut_taster_contrib.append(-1)
        else:
            mut_taster_contrib.append(0)

    # THEN MutationAssessor_pred

    mut_assess_string = full_stuff['MutationAssessor_pred']
    mut_assess_contrib = []

    for mut_assess_pred in mut_assess_string:
        if mut_assess_pred == 'H':
            mut_assess_contrib.append(1)
        elif mut_assess_pred == 'M':
            mut_assess_contrib.append(1)
        elif mut_assess_pred == 'N':
            mut_assess_contrib.append(-1)
        elif mut_assess_pred == 'L':
            mut_assess_contrib.append(-1)
        else:
            mut_assess_contrib.append(0)

    # THEN FATHMM_pred

    fathmm_string = full_stuff['FATHMM_pred']
    fathmm_contrib = []

    for fathmm_pred in fathmm_string:
        if fathmm_pred == 'D':
            fathmm_contrib.append(1)
        elif fathmm_pred == 'T':
            fathmm_contrib.append(-1)
        else:
            fathmm_contrib.append(0)

    # THEN PROVEAN_score
    provean_string = full_stuff['PROVEAN_score']
    provean_mixed = pd.to_numeric(provean_string, errors='coerce')
    provean_contrib = []

    for provean_pred in provean_mixed:
        if provean_pred <= -2.5:
            provean_contrib.append(1)
        elif provean_pred > -2.5:
            provean_contrib.append(-1)
        else:
            provean_contrib.append(0)

    # Then META SVM
    meta_svm_string = full_stuff['MetaSVM_pred']
    meta_svm_contrib = []
    for pred in meta_svm_string:
        if pred == 'D':
            meta_svm_contrib.append(1)
        elif pred == 'T':
            meta_svm_contrib.append(-1)
        else:
            meta_svm_contrib.append(0)

    # Then META LR
    meta_lr_string = full_stuff['MetaLR_pred']
    meta_lr_contrib = []

    for pred in meta_lr_string:
        if pred == 'D':
            meta_lr_contrib.append(1)
        elif pred == 'T':
            meta_lr_contrib.append(-1)
        else:
            meta_lr_contrib.append(0)

    # Then M_CAP

    mcap_string = full_stuff['M-CAP_pred']
    mcap_contrib = []

    for pred in mcap_string:
        if pred == 'D':
            mcap_contrib.append(1)
        elif pred == 'T':
            mcap_contrib.append(-1)
        else:
            mcap_contrib.append(0)

    # Then FATHMM MKL
    mkl_string = full_stuff['fathmm-MKL_coding_pred']
    mkl_contrib = []
    for pred in mkl_string:
        if pred == 'D':
            mkl_contrib.append(1)
        elif pred == 'N':
            mkl_contrib.append(-1)
        else:
            mkl_contrib.append(0)

    # Then geno_canyon with cutoff of 0.5 as used in their paper
    canyon_int = pd.to_numeric(full_stuff['GenoCanyon_score'], errors='coerce')
    canyon_contrib = []
    for pred in canyon_int:
        if pred <= 0.5:
            canyon_contrib.append(-1)
        elif pred > 0.5:
            canyon_contrib.append(1)
        else:
            canyon_contrib.append(0)

    # Then GERP++ use cutoff at 2 or 2.5
    gerp_float = pd.to_numeric(full_stuff['GERP++_RS'], errors='coerce')
    gerp_contrib = []
    for pred in gerp_float:
        if pred >= 2.25:
            gerp_contrib.append(1)
        elif pred < 2.25:
            gerp_contrib.append(-1)
        else:
            gerp_contrib.append(0)

    # Make a matrix with all with shape check
    PP3_contrib = [] # Pathogenic
    BP4_contrib = []     # Benign
    all_contrib = [sift_contrib, lrt_contrib, mut_taster_contrib,
                   mut_assess_contrib, fathmm_contrib, provean_contrib, meta_svm_contrib,
                   meta_lr_contrib, mcap_contrib,
                   mkl_contrib, canyon_contrib, gerp_contrib]
    non_empty = [x for x in all_contrib if x != []]
    number_of_lists = len(non_empty)
    number_of_variants = len(sift_contrib)

    # Scoring is pathogenic if > 3 (not included and benign if inferior to 0
    for x in range(0, number_of_variants):
        total = 0
        for i in range(0, number_of_lists):
            total = total + non_empty[i][x]
        if total > 3:
            PP3_contrib.append(1)
            BP4_contrib.append(0)
        elif 0 <= total <= 3:
            PP3_contrib.append(0)
            BP4_contrib.append(0)
        elif total < 0:
            PP3_contrib.append(0)
            BP4_contrib.append(1)
        else:
            PP3_contrib.append(0)
            BP4_contrib.append(0)
    return [PP3_contrib, BP4_contrib]


def check_PP5(full_stuff):
    path_series = full_stuff['CLNSIG']
    PP5_contrib = []
    for pred in path_series:
        if pred.lower().find('pathogenic') >= 0 or pred.lower().find('conflict') >= 0: # TODO INCLUDE OR NOT CONFLICTING
            PP5_contrib.append(1)
        else:
            PP5_contrib.append(0)
    return PP5_contrib


def check_BA1(freq_geno_series, freq_exo_series):
    BA1_contrib = []
    freq_geno_series.fillna(0, inplace=True)
    freq_exo_series.fillna(0, inplace=True)
    for geno, exo in zip(freq_geno_series, freq_exo_series):
        freq = max(geno, exo)
        if freq >= 0.05:
            BA1_contrib.append(1)
        else:
            BA1_contrib.append(0)
    return BA1_contrib


def check_BS1(full_stuff, cutoff):
    # Select all frequency columns
    filter_columns = [col for col in full_stuff if col.startswith('gnomAD')] + [col for col in full_stuff if col.startswith('ExAC')]
    # Check the max for all populations
    max_freq_series = full_stuff[filter_columns].replace('.', 0)
    max_freq_series = max_freq_series.astype("float")
    max_freq_series = max_freq_series.max(axis=1)

    #max_freq_series = pd.to_numeric(max_freq_series, errors='coerce')
    BS1_contrib = []
    for freq in max_freq_series:
        if type(freq) is not float:
            BS1_contrib.append(0)
        elif freq > cutoff:
            BS1_contrib.append(1)
        else:
            BS1_contrib.append(0)
    return BS1_contrib


def check_BS2(full_stuff, bs2_het_list, bs2_hom_list, rec_list, dom_list, adult_list, ref_anno):
    start_series = full_stuff['Start']
    chr_series = full_stuff['Chr']
    gene_series = full_stuff['Gene.{}'.format(ref_anno)]
    ref_snv = full_stuff['Ref']
    alt_snv = full_stuff['Alt']

    all_tuple = zip(chr_series, start_series, ref_snv, alt_snv, gene_series)
    BS2_contrib =[]
    for chr, start, ref, alt, gene in all_tuple:
        if gene in adult_list:  # Remove all in adult onset
            BS2_contrib.append(0)
        elif gene in rec_list:  # Check for recessive first
            key_snv = str(chr) + '_' + str(start) + '_' + ref + '_' + alt
            if key_snv in bs2_hom_list:
                BS2_contrib.append(1)
            else:
                BS2_contrib.append(0)
        elif gene in dom_list:  # Then check dominant # What if both dominant and recessive
            key_snv = str(chr) + '_' + str(start) + '_' + ref + '_' + alt
            if key_snv in bs2_het_list:
                BS2_contrib.append(1)
            else:
                BS2_contrib.append(0)
        else:
            BS2_contrib.append(0)
    return BS2_contrib


def ckeck_BS3(full_stuff):
    review_series = full_stuff['CLNREVSTAT']
    sig_series = full_stuff['CLNSIG']
    BS3_contrib = []
    for rev, sig in zip(review_series, sig_series):
        if 'benign' in sig.lower() and ('practice_guideline' in rev or 'reviewed_by_expert_panel' in rev):
            BS3_contrib.append(1)
        else:
            BS3_contrib.append(0)
    return BS3_contrib


def check_BP1(full_stuff, BP1_list, ref_anno):
    gene_series = full_stuff['Gene.{}'.format(ref_anno)]
    type_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]
    BP1_contrib = []

    for gene, type_snp in zip(gene_series, type_series):
        if gene in BP1_list and type_snp == 'nonsynonymous SNV':
            BP1_contrib.append(1)
        else:
            BP1_contrib.append(0)
    return BP1_contrib


def check_BP6(full_stuff):
    path_series = full_stuff['CLNSIG']
    BP6_contrib = []
    for pred in path_series:
        if pred.lower().find('benign') > 0:
            BP6_contrib.append(1)
        else:
            BP6_contrib.append(0)
    return BP6_contrib


def check_BP7(full_stuff, ref_anno):
    type_series = full_stuff['ExonicFunc.{}'.format(ref_anno)]
    scsnv_ada_string = pd.to_numeric(full_stuff['dbscSNV_ADA_SCORE'], errors='coerce')
    scsnv_rf_string = pd.to_numeric(full_stuff['dbscSNV_RF_SCORE'], errors='coerce')

    scsnv_contrib =[]
    for func, ada, rf in zip(type_series, scsnv_ada_string, scsnv_rf_string):
        if func == 'synonymous SNV':
                if max(ada, rf) < 0.6:
                    scsnv_contrib.append(1)
                else:
                    scsnv_contrib.append(0)
        else:
            scsnv_contrib.append(0)
    return scsnv_contrib


def decompose_vcf(vcf_input, decomposed_vcf_output):
    try:
        my_parser = VCFParser(infile=vcf_input, split_variants=True, check_info=True)

        with open(decomposed_vcf_output, 'w') as outfile:
            for line in my_parser.metadata.print_header():
                outfile.write(line + '\n')
            for variant in my_parser:
                outfile.write('\t'.join([variant[head] for head in my_parser.header]))
                outfile.write('\n')
    except ModuleNotFoundError:
        print('|| Cannot decompose vcf, please install vcf_parser'
              'try : pip install vcf_parser')
    return


def process_data(full_stuff, ref_anno, number_of_samples, tup_database, cutoff):
    # Load all databases
    print(tmp_stmp()+'Starting...')
    PVS1_list, rek_dict, PS1_dict, bs2_het_list, bs2_hom_list, rec_list, dom_list, adult_list, PS4_df, repeat_dict,\
    PP2_list, BP1_list, data_pm1, whole_dict = tup_database

    # TODO REST INDEXES ?
    full_stuff.reset_index(drop=True, inplace=True)

    if 'chr' in str(full_stuff['Chr'].iloc[0]).lower():
        full_stuff['Chr'].replace(['chr' + s for s in chr_num()], chr_num(), inplace=True)
    df_length = full_stuff.shape[0]
    # TRY PVS1
    try:
    # PVS1 Null variant in a gene where LOF is known mechanism of disease
        # Make a list for splicing scores from dbscSNV
        ada_series = full_stuff['dbscSNV_ADA_SCORE']
        rf_series = full_stuff['dbscSNV_RF_SCORE']
        splice_list = []
        for ada, rf in zip(ada_series, rf_series):
            try:
                spl = max(float(ada), float(rf))
                splice_list.append(spl)
            except ValueError:
                splice_list.append(float('nan'))
                pass

        PVS1_contrib = check_PVS1_criteria(full_stuff, PVS1_list, splice_list, ref_anno, rek_dict)
        print(tmp_stmp()+'PVS1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        print('|| Cannnot calculate PVS1, no splicing annotation. Please annotate with dbscSNV')
        PVS1_contrib = [0] * df_length
        pass

    try:
        # PS1  AA change is the same as a known pathogenic mutation
        PS1_contrib = check_PS1_criteria(full_stuff, PS1_dict, ref_anno)
        print(tmp_stmp() + 'PS1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        print('|| Cannot calculare PS1')
        PS1_contrib = [0] * df_length


    # PS2 De novo (test against parents and no family histories)
    # Process PS2 data Because each family has their own one
    try:
        PS2_contrib = check_PS2_criteria(full_stuff, whole_dict)
        PS2_contrib_merged = []
        ps2_length = len(PS2_contrib)
        ps2_df = list_to_dataframe([x[1] for x in PS2_contrib], ['Trio_'+s+'_PS2' for s in [x[0] for x in PS2_contrib]])
        for index, row in ps2_df.iterrows():
            if 1 in list(row):
                PS2_contrib_merged.append(1)  #  IF at least one is de novo, put 1 in PS2 for this variant
            elif 1 not in list(row):
                PS2_contrib_merged.append(0)
        print(tmp_stmp() + 'PS2 done')
    except (UnboundLocalError, NameError, TypeError) as e:
        PS2_contrib_merged = [0] * df_length
        print("|| No trio data, skipping PS2")

    # PS3 Well established in vitro in vivo functional studies *
    # Considering Clinvar status "reviewed by panel of expert" and "practice guidelines"

    if 'CLNSIG' in full_stuff.columns and 'CLNREVSTAT' in full_stuff.columns:
        PS3_contrib = check_PS3(full_stuff)
        print(tmp_stmp() + 'PS3 done')
    else:
        print('|| No "CLNSIG" column found. Please annotate your data with a recent clinvar database')
        PS3_contrib = [0] * df_length

    # PS4 ODD ratios or database

    if 'gnomAD_exome_ALL' in full_stuff.columns and 'gnomAD_genome_ALL' in full_stuff.columns:
        freq_geno_series = pd.to_numeric(full_stuff['gnomAD_genome_ALL'], errors='coerce')
        freq_exo_series = pd.to_numeric(full_stuff['gnomAD_exome_ALL'], errors='coerce')
        PS4_contrib = check_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df, number_of_samples)
        print(tmp_stmp() + 'PS4 done')
    elif 'gnomAD_genome_ALL' in full_stuff.columns and 'ExAC_ALL' in full_stuff.columns:
        freq_geno_series = pd.to_numeric(full_stuff['gnomAD_genome_ALL'], errors='coerce')
        freq_exo_series = pd.to_numeric(full_stuff['ExAC_ALL'], errors='coerce')
        PS4_contrib = check_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df,number_of_samples)
        print(tmp_stmp() + 'PS4 done')
    else:
        print('|| All frequency data not found')
        PS4_contrib = [0] * df_length

    if len(PS4_contrib) == 4:  # Differentiate if PS4 was aquired by direct calculation or just using db
        odd_ratios_list = PS4_contrib[1]

        confidence_interval_list = PS4_contrib[2]

        p_value_list = PS4_contrib[3]

        PS4_contrib = PS4_contrib[0]

    # PM1 List of begning domains (to check against interpro domain affected)
    try:
        full_stuff['Interpro_domain']
        PM1_contrib = check_PM1_criteria(full_stuff, data_pm1['Domain'], ref_anno)
        print(tmp_stmp() + 'PM1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        print("|| No domain data found")
        PM1_contrib = [0] * df_length


    # PM2 Absent from controls (eg gnomad = 0) be careful about indels IN TABLE
    # Use gnomad_genome and either gnomad exome or exac for exome data
    if 'gnomAD_exome_ALL' in full_stuff.columns and 'gnomAD_genome_ALL' in full_stuff.columns:
        freq_geno_series = full_stuff['gnomAD_genome_ALL']
        freq_exo_series = full_stuff['gnomAD_exome_ALL']
        freq_exo_series = pd.to_numeric(freq_exo_series, errors='coerce')
        freq_geno_series = pd.to_numeric(freq_geno_series, errors='coerce')
        PM2_contrib = check_PM2(full_stuff, freq_exo_series, freq_geno_series, ref_anno)
        print(tmp_stmp() + 'PM2 done')
    elif 'gnomAD_genome_ALL' in full_stuff.columns and 'ExAC_ALL' in full_stuff.columns:
        freq_geno_series = full_stuff['gnomAD_genome_ALL']
        freq_exo_series = full_stuff['ExAC_ALL']
        freq_exo_series = pd.to_numeric(freq_exo_series, errors='coerce')
        freq_geno_series = pd.to_numeric(freq_geno_series, errors='coerce')
        PM2_contrib = check_PM2(full_stuff, freq_exo_series, freq_geno_series, ref_anno)
        print(tmp_stmp() + 'PM2 done')
    else:
        PM2_contrib = [0] * df_length

    # PM4 change of length (stoploss or nonframeshift mutation) in a non repeat region (use repeatmask file)
        # Load repeat mask dictionary
    try :
        PM4_contrib, BP3_contrib = check_PM4_and_BP3_repeat_region(repeat_dict, full_stuff, ref_anno)

        print(tmp_stmp() + 'PM4 and BP3 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PM4_contrib = [0] * df_length
        BP3_contrib = [0] * df_length

    # PM5 is basically PS1 with different amino acid
    try:
        PM5_contrib = check_PM5_criteria(full_stuff, PS1_dict, ref_anno)
        print(tmp_stmp() + 'PM5 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PM5_contrib = [0] * df_length
    # Check : IF R89H and R89L are both pathogenic and we have R89H does it still apply since R89L exist and it is
    #  already caught by PS1

    # PP2 Low rate of missense variants (a list of genes) Using Own DB
    try:
        PP2_contrib = check_PP2(full_stuff, PP2_list, ref_anno)
        print(tmp_stmp() + 'PP2 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PP2_contrib = [0] * df_length

    # PP5 is clinvar data

    if 'CLNSIG' in full_stuff.columns:
        PP5_contrib = check_PP5(full_stuff)
        print(tmp_stmp() + 'PP5 done')
    else:
        PP5_contrib = [0] * df_length
    try:
        PP3_BP4_contribs = check_PP3_and_BP4_dbnsfp(full_stuff)
        PP3_contrib = PP3_BP4_contribs[0]
        BP4_contrib = PP3_BP4_contribs[1]
        print(tmp_stmp() + 'PP3 and BP4 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PP3_contrib = [0] * df_length
        BP4_contrib = [0] * df_length



    try:
        BA_1_contrib = check_BA1(freq_exo_series, freq_geno_series)
        print(tmp_stmp() + 'BA1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BA_1_contrib = [0] * df_length

    try:
        BS1_contrib = check_BS1(full_stuff, cutoff)
        print(tmp_stmp() + 'BS1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BS1_contrib = [0] * df_length

    # BS2 Benign variants observed in dominant (heterozygote) and recessive (homozygote) diseases (or X linked)
    try:
        BS2_contrib = check_BS2(full_stuff, bs2_het_list, bs2_hom_list, rec_list, dom_list, adult_list, ref_anno)
        print(tmp_stmp() + 'BS2 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BS2_contrib = [0] * df_length

    # BS3 ALSO CLINVAR with reviewed by expert panel
    try:
        BS3_contrib = ckeck_BS3(full_stuff)
        print(tmp_stmp() + 'BS3 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BS3_contrib = [0] * df_length

    try:
        # BP1 Missense variant in gene where missense variants are not pathogenic
        BP1_contrib = check_BP1(full_stuff, BP1_list, ref_anno)
        print(tmp_stmp() + 'BP1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BP1_contrib = [0] * df_length

    # BP6 Clinvar benign maybe we can merge BP6 and PP5 to save time ?
    # Intervar auto db is not too reliable with this one
    try:
        BP6_contrib = check_BP6(full_stuff)
        print(tmp_stmp() + 'BP6 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BP6_contrib = [0] * df_length

    try:
        BP7_contrib = check_BP7(full_stuff, ref_anno)
        print(tmp_stmp() + 'BP7 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BP7_contrib = [0] * df_length

    print(tmp_stmp()+'ACMG Criteria assigned')

    all_criterias = [
    PVS1_contrib,
    PS1_contrib,
    PS2_contrib_merged,
    PS3_contrib,
    PS4_contrib,
    PM1_contrib ,
    PM2_contrib,
    PM4_contrib,
    PM5_contrib,
    PP2_contrib,
    PP3_contrib,
    PP5_contrib,
    BS1_contrib,
    BS2_contrib,
    BS3_contrib ,
    BP1_contrib,
    BP3_contrib,
    BP4_contrib,
    BP6_contrib,
    BP7_contrib,
    BA_1_contrib
    ]

    names_list = ['PVS1_contrib', 'PS1_contrib', 'PS2_contrib_merged', 'PS3_contrib', 'PS4_contrib', 'PM1_contrib', 'PM2_contrib',
    'PM4_contrib', 'PM5_contrib', 'PP2_contrib', 'PP3_contrib','PP5_contrib', 'BS1_contrib', 'BS2_contrib',
                  'BS3_contrib', 'BP1_contrib', 'BP3_contrib', 'BP4_contrib','BP6_contrib', 'BP7_contrib', 'BA_1_contrib']

    df_all = list_to_dataframe(all_criterias, names_list)
    prob_pred_list = calculate_probability(df_all)
    prob_list = prob_pred_list[0]
    pred_list = prob_pred_list[1]

    # Merging with original dataframe

    # THINGS TO ADD FIRST => ODD ratio and confidence interval / PS2 for each family
    try:
        full_stuff = pd.concat([full_stuff, ps2_df], axis=1, sort=False)
    except UnboundLocalError:
        pass

    try:
        odd_series = pd.Series(odd_ratios_list)
        ci_series = pd.Series(confidence_interval_list)
        p_value_series = pd.Series(p_value_list)
        full_stuff['Odd_ratio'] = odd_series.values
        full_stuff['CI'] = ci_series.values
        full_stuff['p_value'] = p_value_series
    except UnboundLocalError:
        pass

    final_df = pd.concat([full_stuff, df_all], axis=1, sort=False) # TODO FIX THIS SHIT

    prob_series = pd.Series(prob_list)
    pred_series = pd.Series(pred_list)
    final_df['Probability_Path'] = prob_series.values
    final_df['Prediction_ACMG_tapes'] = pred_series.values

    return final_df


def test_if_decomposed(vcf_file):
    process_path(vcf_file)
    decomposed = True
    # Test the first 2000 entries in vcf to see if there is multi allelic sites
    with open(vcf_file, 'r') as vcf:
        for i in range(1, 2000):
            try:
                line = vcf.readline().split('\t')
                if len(line[4].split(',')) > 1:
                    decomposed = False
                    break
                else:
                    pass
            except IndexError:
                pass
    if not decomposed:
        print(tmp_stmp()+'Decomposing vcf', end=' ', flush=True)
        outfile = vcf_file.replace('.vcf', '_decomposed.vcf')
        decompose_vcf(vcf_file, outfile)
        print("Done")
    elif decomposed:
        outfile = vcf_file
        print(tmp_stmp()+'Vcf already decomposed')
    # Return both the name of the vcf output and the state (if the input was decomposed or not)
    return outfile, decomposed


def reference_used(full_stuff):
    filter_col = [col for col in full_stuff if col.startswith(('Gene.', 'GeneDetail.'))]
    if filter_col[0].split('.')[1] == filter_col[1].split('.')[1]:
        ref_used = filter_col[0].split('.')[1]
        print(tmp_stmp()+'Gene annotation : '+ref_used)
    else:
        print('Could not determine annotation used, exiting...')
        sys.exit(1)
    return ref_used


def check_acmg_tag_requirement(full_table, ref_anno):
    required_annotations = [
        'Func.{}'.format(ref_anno),
        'Gene.{}'.format(ref_anno),
        'GeneDetail.{}'.format(ref_anno),
        'ExonicFunc.{}'.format(ref_anno),
        'AAChange.{}'.format(ref_anno),
        'pLi.{}'.format(ref_anno),
        'pRec.{}'.format(ref_anno),
        'pNull.{}'.format(ref_anno),
        'CLNSIG',
        'Interpro_domain',
        'dbscSNV_ADA_SCORE',
        'dbscSNV_RF_SCORE',
        'gnomAD_genome_ALL',
        "SIFT_score",
        'LRT_pred',
        'MutationTaster_pred',
        'MutationAssessor_pred',
        'FATHMM_pred',
        'PROVEAN_score',
        'MetaSVM_pred',
        'MetaLR_pred',
        'M-CAP_pred',
        'fathmm-MKL_coding_pred',
        'GenoCanyon_score',
        'GERP++_RS'
    ]

    if set(required_annotations) <= set(list(full_table.columns)) and ('ExAC_ALL' in list(full_table.columns) or
                                                                       'gnomAD_exome_ALL' in list(full_table.columns)):
        print(tmp_stmp()+'All required annotations found...')

    else:
        print(tmp_stmp()+'All required annotations not found. '
                         'Try to remove the --acmg tag\n'
                         +tmp_stmp()+'EXITING')
        sys.exit()


def output_type(output_arg):
    if output_arg[-1] == '/':
        if os.path.isdir(output_arg) and len(os.listdir(output_arg)) > 0:
            print(tmp_stmp()+'This directory already exists and is not empty, exiting...')
            sys.exit(1)
        elif os.path.isdir(output_arg) and len(os.listdir(output_arg)) == 0:
            print(tmp_stmp()+'Output type: FOLDER')
            return 'directory'
        elif os.path.isdir(output_arg) is False:
            print(tmp_stmp()+'Output type: FOLDER')
            os.makedirs(output_arg)
            return 'directory'
    elif output_arg[-1] == '\\':
        if os.path.isdir(output_arg) and len(os.listdir(output_arg)) > 0:
            print(tmp_stmp()+'This directory already exists and is not empty, exiting...')
            sys.exit(1)
        elif os.path.isdir(output_arg) and len(os.listdir(output_arg)) == 0:
            print(tmp_stmp()+'Output type: FOLDER')
            return 'directory'
        elif os.path.isdir(output_arg) is False:
            print(tmp_stmp()+'Output type: FOLDER')
            os.makedirs(output_arg)
            return 'directory'
    elif output_arg[-4:] == '.csv':
        print(tmp_stmp()+'Output type: CSV + XLSX')
        return 'csv'
    elif output_arg[-4:] == '.txt' or output_arg[-4:] == '.tsv':
        print(tmp_stmp()+'Output type: TXT/TSV + XLSX')
        return 'txt'

    else:
        print(tmp_stmp()+'Could not determine output type. Exiting...')
        sys.exit(1)


def csv_report(full_stuff, output_path, suffix, ref_anno, tab_output):

    if tab_output is True:
        separator = '\t'
        extension = '.txt'
    else:
        separator = ','
        extension = '.csv'

    wot_wot = full_stuff['Func.{}'.format(ref_anno)].value_counts()
    plot_loc = wot_wot.plot(kind='pie', figsize=(15, 12), legend=False, fontsize=12, rot=15,
                            title='Genomic Location of all variants')
    fig_loc = plot_loc.get_figure()
    fig_loc.savefig(os.path.join(output_path, (suffix + '_location')))
    fig_loc.clear()

    type_series = full_stuff['ExonicFunc.{}'.format(ref_anno)].value_counts()
    plot_type = type_series.plot(kind='pie', figsize=(15, 12), legend=False, fontsize=12, rot=15,
                                 title='Type of variants')
    fig_type = plot_type.get_figure()
    fig_type.savefig(os.path.join(output_path, (suffix + '_type')))
    fig_type.clear()

    pred_series = full_stuff['Prediction_ACMG_tapes'].value_counts()
    plot_pred = pred_series.plot(kind='pie', figsize=(15, 12), legend=False, fontsize=12, rot=15,
                                 title='Tapes prediction')
    fig_pred = plot_pred.get_figure()
    fig_pred.savefig(os.path.join(output_path, (suffix + '_prediction')))
    fig_pred.clear()

    try:
        path_cyto = full_stuff.loc[full_stuff['Prediction_ACMG_tapes'].isin(['Pathogenic', 'Likely Pathogenic'])]
        cyto_series = path_cyto['cytoBand'].value_counts()
        plot_cyto = cyto_series.plot(kind='bar', figsize=(15, 12), legend=False, fontsize=10, rot=35,
                                     title='Genomic Location of pathogenic variants by cytoband')
        fig_cyto = plot_cyto.get_figure()
        fig_cyto.savefig(os.path.join(output_path, (suffix + '_cytoband')))
        fig_cyto.clear()
    except (KeyError, TypeError) as e:
        pass

    try:
        pathways_results.to_csv((os.path.join(output_path, (suffix + '_' + pathway_db_used+extension))), index=False, sep=separator)
    except NameError:
        pass

    try:
        cancer_df.to_csv((os.path.join(output_path, (suffix + '_' + disease_name+extension))), index=False, sep=separator)
    except NameError:
        pass

    try:
        df_list
        print("|| Exporting most pathogenic variants per sample")
        for sample, df in df_list:
            with open((os.path.join(output_path, (suffix + '_by_sample'+extension))), 'a') as out_csv:
                out_csv.write(sample+'\n')
            df.to_csv((os.path.join(output_path, (suffix + '_by_sample'+extension))), mode='a', header=True, index=False, sep=separator)
            with open((os.path.join(output_path, (suffix + '_by_sample'+extension))), 'a') as out_csv:
                out_csv.write('\n')
    except NameError:
        pass

    try:
        list_of_genes_df.to_csv((os.path.join(output_path, (suffix + '_gene_list'+extension))), index=False, sep=separator)
        print("|| Exporting file for user supplied gene list")
    except NameError:
        pass

        # ADD 7th spreadsheet to show keeg pathway reverse search
    try:
        kegg_df
        print("|| Exporting file for {}...".format(name_of_pathway_to_report))
        kegg_df.to_csv((os.path.join(output_path, (suffix + '_' + name_of_pathway_to_report + extension))), index=False, sep=separator)
    except (NameError, UnboundLocalError) as e:
        pass

    try:
        whole
        print('|| Exporting per gene report...')
        for tup in whole:
            with open((os.path.join(output_path, (suffix + '_' + 'per_gene' + extension))), 'a') as outfile:
                if tup[6] is True:
                    outfile.write(tup[0] + separator + str(tup[3]) + separator + "FLAGS GENE" + '\n')
                elif int(tup[4]) > 250000:
                    outfile.write(tup[0] + separator + str(tup[3]) + separator + "LONG GENE" + '\n')
                elif tup[5] is True:
                    outfile.write(tup[0] + separator + str(tup[3]) + separator + "SAMPLE NUMBER WARNING" + '\n')
                else:
                    outfile.write(tup[0] + separator + str(tup[3]) + '\n')
            tup[1].to_csv((os.path.join(output_path, (suffix + '_' + 'per_gene' + extension))), mode='a', index=False, sep=separator)
            with open((os.path.join(output_path, (suffix + '_' + 'per_gene' + extension))), 'a') as outfile:
                outfile.write('\n')
    except (NameError, UnboundLocalError) as e:
        pass

    return


def process_input_variant_file(variant_file_path):
    path = process_path(variant_file_path)
    if not path.endswith('.vcf'):
        try:
            bcf_in = VariantFile(path)  # auto-detect input format
            print(tmp_stmp()+'Converting to vcf...')
            with VariantFile(path+'_tmp.vcf', 'w', header=bcf_in.header) as bcf_out:
                for rec in bcf_in.fetch():
                    bcf_out.write(rec)
            return path+'_tmp.vcf'

        except NotImplementedError:
            my_parser = VCFParser(infile=path, split_variants=True, check_info=True)
            print(tmp_stmp()+'Converting to vcf...')
            with open(path+'_tmp.vcf', 'w') as outfile:
                for line in my_parser.metadata.print_header():
                    outfile.write(line + '\n')
                for variant in my_parser:
                    outfile.write('\t'.join([variant[head] for head in my_parser.header]))
                    outfile.write('\n')
            return path + '_tmp.vcf'
    elif path.endswith('.vcf'):
        return path
    else:
        print(tmp_stmp()+"Input does not seem to be a vcf file. Exiting...")
        sys.exit(1)


def process_annotated_vcf(input_path):
    '''outfile_path = input_path + '_tmp.vcf'

    # REMOVE HEADER AND WRITE TEMP FILE
    with open(outfile_path, 'w') as outfile:
        with open(input_path, 'r') as input_vcf:
            for line in input_vcf.readlines():
                if line.startswith('##'):
                    pass
                else:
                    outfile.write(line)'''
    print(tmp_stmp()+'Processing vcf...', end=' ', flush=True)
    header_row = 0
    with open(input_path, 'r') as input_vcf:
        for line in input_vcf.readlines():
            if line.startswith('##'):
                header_row = header_row + 1
            else:
                break

    df = pd.read_csv(input_path, delimiter='\t',header=header_row)
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)  # REPLACE FIRST COLUMN NAME

    # Check if vcf is actually annotated
    if 'ANNOVAR_DATE=' in str(df['INFO'].iloc[0]):
        pass
    else:
        print('|| This vcf does not seem to be annotated, please use :\n || python tapes.py annotate')
        sys.exit(1)

    # Extract genotyping data

    geno_df = df.drop(df.columns[[0,1,2,3,4,5,6,7]], axis=1)
    number_of_samples_anno = geno_df.shape[1] - 1   # -1 To account for format column
    if number_of_samples_anno > 1:
        # SPLIT THE GENOTYPE OF EACH SAMPLES
        sep = ':'
        for column in geno_df.columns:
            value = geno_df[column].str.split(sep, n=1, expand=True)
            geno_df[column] = value[0]

        hom_list = []
        het_list = []
        norm_list = []

        for line in range(0, geno_df.shape[0]):
            values_of_line = geno_df.iloc[line].values
            hom_count = 0
            het_count = 0
            normal_count = 0
            for i in values_of_line:
                if i in ["1/1", "2/2", "1|1", "2|2", '1/2', '1|2', '2|1']:
                    hom_count = hom_count + 1
                elif i in ['0/1', '0/2', ' 0|1', '0|2', '1|0', '2|0']:
                    het_count = het_count + 1
                elif i == '0/0' or i == '0|0':
                    normal_count = normal_count + 1
            hom_list.append(hom_count)
            het_list.append(het_count)
            norm_list.append(normal_count)

        norm_in_df = pd.Series(norm_list)
        het_in_df = pd.Series(het_list)
        hom_in_df = pd.Series(hom_list)

        geno_df['WT count'] = norm_in_df.values
        geno_df['Het count'] = het_in_df.values
        geno_df['Hom count'] = hom_in_df.values

    else:
        print(tmp_stmp() + 'Only one sample found in vcf')

    # PROCESS VCF TO COLUMNS LOOKING LIKE  ANNO=VALUE
    info_df = df['INFO']
    info_df = info_df.str.split('ALLELE_END', n=1, expand=True)
    info_df = info_df[0].str.split('ANNOVAR_DATE=', n=1, expand=True)
    info_df = info_df[1].str.split(';', n=1, expand=True)
    info_df = info_df[1].str.split(';', n=-1, expand=True)
    info_df.drop(info_df.columns[len(info_df.columns) - 1], axis=1, inplace=True)  # TO REMOVE EMPTY LAST COLUMN

    # GET NAMES OF COLUMNS USING THE FIRST LINE (EVERYTHING BEFORE THE = SIGN)
    columns_tmp = list(info_df.iloc[0])
    column_list = []
    for thingy in columns_tmp:
        item = thingy.split('=')[0]
        column_list.append(item)
    #  CREATE A DICTIONARY TO EASILY TRANSFORM INTO DATAFRAME
    df_dict = {}
    for index, value in enumerate(column_list):
        df_dict[value] = list(info_df[index])
    # REPLACE ALL ITEMS IN LIST BY ONLY THE VALUE
    for key in df_dict:
        new_list = []
        for item in df_dict[key]:
            item = item.split('=')[1]
            new_list.append(item)
        df_dict[key] = new_list

    new_df = pd.DataFrame.from_dict(df_dict)
    basic_info_df = df[['CHROM', 'POS', 'REF', 'ALT']].copy()  #  Use .copy() to remove SettingWithCopyWarning

    pos_series = basic_info_df['POS']
    ref_series = basic_info_df['REF']
    end_list = []

    for pos, ref in zip(pos_series, ref_series):
        end = int(pos) + len(ref) - 1
        end_list.append(end)

    end_in_df = pd.Series(end_list)
    basic_info_df['End'] = end_in_df.values
    basic_info_df.columns = ['Chr', 'Start', 'Ref', 'Alt', 'End']
    basic_info_df = basic_info_df[['Chr', 'Start', 'End', 'Ref', 'Alt']]

    if basic_info_df.shape[0] == new_df.shape[0]:
        full_stuff_tmp = pd.concat([basic_info_df, new_df], axis=1, sort=False)

        if full_stuff_tmp.shape[0] == geno_df.shape[0]:
            full_stuff = pd.concat([full_stuff_tmp, geno_df], axis=1, sort=False)
        else:
            print('|| Dataframes shapes not matching during vcf processing, exiting...')
            sys.exit(1)
        print('Done')
        return full_stuff

    else:
        print('|| Dataframes shapes not matching during vcf processing, exiting...')
        sys.exit(1)


def gene_grouped_report(full_stuff, ref_anno, acmg_db_path, number_of_samples):
    global whole
    length_dict = {}
    with open(os.path.join(acmg_db_path, "gene_length.txt"), 'r') as infile:
        for line in infile.readlines():
            gene = line.split('\t')[2].strip()
            length = line.split('\t')[1].strip()
            length_dict[gene] = length

    full_stuff = full_stuff[full_stuff['Probability_Path'] >= 0.80]  # Remove non pathogenic variants
    grouped = full_stuff.groupby('Gene.{}'.format(ref_anno))  # GROUP BY GENES
    whole = []
    #  LOOP OVER ALL GENES TO remove columns/samples that only contains 0/0

    for truc in grouped:
        gene_df = truc[1]  # iterate over individual df for each gene
        list_cols = list(gene_df.columns)  # get list of columns
        format_index = list_cols.index("FORMAT")  # get index of format column
        wt_index = list_cols.index("WT count")  # get index of wt count column
        samp_list = list_cols[format_index + 1:wt_index]  # make the list of sample columns based on indexes
        samp_df = gene_df[samp_list]
        for column in samp_df.columns:  # make df of only genotyping matrix
            #  Remove all columns with no variants
            if ("0/1" or '0/2' or "1/1" or "1/2" or "1|0" or "0|2" or "0|1" or "2|0" or "1|1" or "2|2" or "1|2" or "2|1") \
                    not in list(samp_df[column]):
                samp_df.pop(column)
        list_samp_per_var = []
        for line in range(0, samp_df.shape[0]):  # for each variant of each gene
            serie = samp_df.iloc[line]
            for column in serie.index:
                if serie[column] == '0/0' or serie[column] == './.':
                    serie.drop(index=column, inplace=True)
            list_samp_per_var.append(", ".join(serie.axes[0].tolist()))
        del list_cols[7:-2]  # Remove all columns unless required
        gene_df = gene_df[list_cols].copy()
        gene_df['Samples'] = pd.Series(list_samp_per_var).values
        longueur = gene_df.shape[0]
        proba_series = gene_df['Probability_Path']
        samp_series = gene_df['Samples']
        score_top = 0
        for proba, samp in zip(proba_series, samp_series):
            score_top = score_top + (proba * (len(samp.split(','))))
        try:
            gene_length = length_dict[truc[0].strip()]
        except KeyError:
            gene_length = 0

        # Creates a warning if :
        # -at least half of the variants in any gene are affecting half or more of the total cohort
        # Suspicious because too abundant

        warning_sample_count = 0 # To raise warning about too many samples affected by variants in gene
        for line in range(0, gene_df.shape[0]):  # for each variant of each gene
            samp_affected = len(gene_df.iloc[line]['Samples'].split(','))
            if samp_affected >= number_of_samples/2:
                warning_sample_count = warning_sample_count + 1
        if warning_sample_count >= gene_df.shape[0]: # if more than half of the variant in this gene are suspicious
            suspicious_number_of_samples = True
        else:
            suspicious_number_of_samples = False

        #  truc[0] is the name of the gene
        # Setting up FLAGS warning
        if truc[0] in flags():
            flag_warning = True
        else:
            flag_warning = False

        tup = (truc[0], gene_df, longueur, score_top, int(gene_length), suspicious_number_of_samples, flag_warning)
        whole.append(tup)
    whole = sorted(whole, key=lambda key: key[3], reverse=True)
    return


def csv_report_reanalyse(output_path, tab_output):

    if tab_output is True:
        separator = '\t'
    else:
        separator = ','

    try:
        pathways_results.to_csv(output_path, index=False, sep=separator)
    except NameError:
        pass

    try:
        cancer_df.to_csv(output_path, index=False, sep=separator)
    except NameError:
        pass

    try:
        df_list
        print("|| Exporting most pathogenic variants per sample")
        for sample, df in df_list:
            with open(output_path, 'a') as out_csv:
                out_csv.write(sample+'\n')
            df.to_csv(output_path, mode='a', header=True, index=False, sep=separator)
            with open(output_path, 'a') as out_csv:
                out_csv.write('\n')
    except NameError:
        pass

    try:
        list_of_genes_df.to_csv(output_path, index=False, sep=separator)
        print("|| Exporting file for user supplied gene list")
    except NameError:
        pass

        # ADD 7th spreadsheet to show keeg pathway reverse search
    try:
        kegg_df
        print("|| Exporting file for {}...".format(name_of_pathway_to_report))
        kegg_df.to_csv(output_path, index=False, sep=separator)
    except (NameError, UnboundLocalError) as e:
        pass

    try:
        whole
        print('|| Exporting per gene report...')
        if os.path.isfile(output_path):
            os.remove(output_path)
        for tup in whole:
            with open(output_path, 'a') as outfile:
                if tup[6] is True:
                    outfile.write(tup[0] + separator + str(tup[3]) + separator + "FLAGS GENE" + '\n')
                elif int(tup[4]) > 250000:
                    outfile.write(tup[0] + separator + str(tup[3]) + separator + "LONG GENE" + '\n')
                elif tup[5] is True:
                    outfile.write(tup[0] + separator + str(tup[3]) + separator + "SAMPLE NUMBER WARNING" + '\n')
                else:
                    outfile.write(tup[0] + separator + str(tup[3]) + '\n')
            tup[1].to_csv(output_path, mode='a', index=False, sep=separator)
            with open(output_path, 'a') as outfile:
                outfile.write('\n')
    except (NameError, UnboundLocalError) as e:
        pass
    return


# Load all databases and then do the analysis (takes around 30 sec to load) TODO optimise libaries format/loagind time
def load_databases(acmg_db_path, ref_anno, assembly, trio_data, pp2_percent, pp2_min, bp1_percent):
    acmg_db_path = process_path(acmg_db_path)
    print(tmp_stmp()+'Loading databases...', end=" ", flush=True)

    ##### PVS1 LIST ########
    PVS1_list = []
    with open(os.path.join(acmg_db_path, "final_PVS1.txt"), 'r') as pvs1_file:
        for line in pvs1_file:
            new_line = line.replace('\n', '')
            PVS1_list.append(new_line)
    #####  rek dict ########
    rek_df = pd.read_csv(os.path.join(acmg_db_path, "REK_canon.{}".format(assembly)), sep='\t')
    if ref_anno == 'refGene':
        rek_custom = rek_df[['refgene', 'end']]
        rek_custom.columns = ['ref_anno', 'end']
        rek_custom = rek_custom[rek_custom['ref_anno'] != '.']
    elif ref_anno == 'knownGene':
        rek_custom = rek_df[['Transcript(ucsc/known)', 'end']]
        rek_custom.columns = ['ref_anno', 'end']
        rek_custom = rek_custom[rek_custom['ref_anno'] != '.']
    elif ref_anno == 'ensGene':
        rek_custom = rek_df[['ensgene', 'end']]
        rek_custom.columns = ['ref_anno', 'end']
        rek_custom = rek_custom[rek_custom['ref_anno'] != '.']
    rek_dict = dict(zip(rek_custom['ref_anno'], rek_custom['end']))

    ##### PS1 dict ########
    PS1_dict = {}
    PS1_df = pd.read_csv(open(os.path.join(acmg_db_path, 'clinvar.PS1.{}.txt'.format(assembly)), 'r'),
                         sep='\t',
                         header=0,
                         names=['Chr', 'Start', 'Stop', 'Ref', 'Alt', 'Pred', 'AA_change', 'AA_ref', 'AA_alt'],
                         dtype={'Chr': str, 'Start': int, 'Stop': int, 'Ref': str, 'Alt': str, 'AA_ref': str,
                                'AA_alt': str, 'AA_change': str, 'Pred': str})
    for chr, start, stop, aref, aalt in zip(PS1_df['Chr'], PS1_df['Start'], PS1_df['Stop'], PS1_df['AA_ref'], PS1_df['AA_alt']):
        PS1_dict[chr, start, stop] = {'AA_ref': aref, 'AA_alt': aalt}

    ##### BS2 whole thing ########

    bs2_het_list = []
    bs2_hom_list = []
    with gzip.open(os.path.join(acmg_db_path, "BS2_hom_het.{}".format(assembly)), 'rt') as bs2_file:
        for line in bs2_file.readlines():
            chr = line.split(' ')[0].strip()
            pos = line.split(' ')[1].strip()
            ref = line.split(' ')[2].strip()
            alt = line.split(' ')[3].strip()
            hom = line.split(' ')[4].strip()
            het = line.split(' ')[5].strip()
            if str(het) == '1':
                bs2_het_list.append(str(chr) + '_' + str(pos) + '_' + ref + '_' + alt)
            if str(hom) == '1':
                bs2_hom_list.append(str(chr) + '_' + str(pos) + '_' + ref + '_' + alt)

    bs2_hom_het_ad_df = pd.read_csv(os.path.join(acmg_db_path, 'BS2_rec_dom_ad.txt'), sep='\t')

    rec_list = []
    dom_list = []
    adult_list = []
    for gene, rec, dom, adult in zip(bs2_hom_het_ad_df['gene'], bs2_hom_het_ad_df['recessive'],
                                     bs2_hom_het_ad_df['dominant'], bs2_hom_het_ad_df['adult_onset']):
        if rec == 1:
            rec_list.append(gene)
        if dom == 1:
            dom_list.append(gene)
        if adult == 1:
            adult_list.append(gene)

    ##### wholde dict for trio ########
    if trio_data is not None: # IF there is trio data provided
        fam_df = pd.read_csv(trio_data, sep='\t', header=None, names=['role', 'fam_id', 'samp_id'])
        # Import trio data as dataframe
        fam_id_series = fam_df['fam_id']
        family_list = fam_id_series.value_counts().index.values.tolist() # Detect all the unique family IDs
        whole_dict = {}
        for famid_from_list in family_list: # For each family ID do a dictionary
            tmp_df = fam_df[fam_df['fam_id'] == famid_from_list]
            tmp_df = tmp_df[['role', 'samp_id']]
            famid_from_list_dict = {
                'm' : tmp_df[tmp_df['role'] == 'm']['samp_id'].to_string(index=False),
                'f': tmp_df[tmp_df['role'] == 'f']['samp_id'].to_string(index=False),
                'o': tmp_df[tmp_df['role'] == 'o']['samp_id'].to_string(index=False)
            }
            whole_dict[famid_from_list] = famid_from_list_dict
    else:
        whole_dict = None

    ##### PS4 df  ########
    PS4_df = pd.read_csv(open(os.path.join(acmg_db_path, "PS4.variants.{}".format(assembly)), 'r'), sep='\t',
                         dtype={"CHR": str, 'POS_hg19': int, 'SNPID': str, 'REF': str, 'ALT': str})
    ##### repeat dict ########

    with gzip.open(os.path.join(acmg_db_path, 'repeat_dict.{}.gz'.format(assembly)), "r") as dj:
        repeat_dict = json.load(dj)

    ##### PP2 BP1 ########

    pp2_bp1_path = os.path.join(acmg_db_path, "PP2_BP1.txt")
    pp2_bp1_df = pd.read_csv(pp2_bp1_path, sep='\t')
    pathogenic_80 = pp2_bp1_df.loc[pp2_bp1_df['Percent_Missense_Path'] >= pp2_percent]
    pathogenic_80_two = pathogenic_80.loc[pathogenic_80['Missense_Path'] >= pp2_min]

    PP2_list = list(pathogenic_80_two['Gene'])
    benign_10 = pp2_bp1_df.loc[pp2_bp1_df['Percent_Missense_Path'] <= bp1_percent]
    BP1_list = list(benign_10['Gene'])

    ###### PM1 BENIGN DICT ######

    with open(os.path.join(acmg_db_path, "PM1_domains_with_benigns.{}".format(assembly)), 'r') as pm1_file:
        list_chr = []
        list_gene = []
        list_domain = []
        for line in pm1_file:
            chr = line.split(" ")[0]
            gene = line.split(" ")[1]
            domain = " ".join(list(line.split(" ")[2:])).replace('\n', '')

            list_chr.append(chr)
            list_gene.append(gene)
            list_domain.append(domain)

        data_pm1 = {'Chr': list_chr, 'Gene': list_gene, 'Domain': list_domain}

    print("Done")

    tup_database = (PVS1_list
    , rek_dict
    , PS1_dict
    , bs2_het_list, bs2_hom_list, rec_list, dom_list, adult_list
    , PS4_df
    , repeat_dict
    , PP2_list
    , BP1_list
    , data_pm1
    , whole_dict)

    return tup_database


def amcg_db_download(folder_path):
    if os.path.isdir(folder_path):
        pass
    else:
        os.mkdir(folder_path)

    list_db = ['BS2_hom_het.hg19',
               'BS2_hom_het.hg38',
               'BS2_rec_dom_ad.txt',
               'clinvar.PS1.hg19.txt',
               'clinvar.PS1.hg38.txt',
               'final_PVS1.txt',
               'gene_length.txt',
               'kegg_dict.json',
               'PM1_domains_with_benigns.hg19',
               'PM1_domains_with_benigns.hg38',
               'PP2_BP1.txt',
               'PS4.variants.hg19',
               'PS4.variants.hg38',
               'REK_canon.hg19',
               'REK_canon.hg38',
               'repeat_dict.hg19.gz',
               'repeat_dict.hg38.gz',
               'transcripts_per_genes.txt']

    sizes = [21648,
             21642,
             63,
             298,
             295,
             14,
             382,
             59,
             32,
             32,
             25,
             13,
             13,
             621,
             1005,
             40759,
             40489,
             74
             ]

    for db_name, size in zip(list_db, sizes):

        if os.path.isfile(os.path.join(folder_path, db_name)):
            pass
        else:
            url = 'https://github.com/a-xavier/tapes/blob/master/acmg_db/{}?raw=true'.format(db_name)
            r = requests.get(url, allow_redirects=True, stream=True)
            if r.status_code == 200:
                with open(os.path.join(folder_path, db_name), 'wb') as f:
                    x = 0
                    for chunk in r.iter_content(chunk_size=1024):
                        x = x + 1
                        per = str(x)
                        print('{} downloading, {}/{}kb'.format(db_name, per, size), end='\r')
                        if chunk:  # filter out keep-alive new chunks
                            f.write(chunk)
                    print('')
            else:
                print('Failed to download {}'.format(db_name))
    return


def process_df_for_trio(tup_databases, full_stuff):
    list_healthy = []
    list_affected = []

    for key, value in tup_databases[13].items():
        list_healthy.append(tup_databases[13][key]['m'])
        list_healthy.append(tup_databases[13][key]['f'])
        list_affected.append(tup_databases[13][key]['o'])

    list_healthy = list(set(list_healthy))
    list_affected = list(set(list_affected))

    for affected in list_affected:
        if affected in list_healthy:
            list_healthy.remove(affected)

    healthy_df = full_stuff[list_healthy]
    full_stuff.drop(list_healthy, axis=1, inplace=True)
    full_stuff.drop(['WT count', 'Het count', 'Hom count'], axis=1, inplace=True)

    hom_list = []
    het_list = []
    norm_list = []

    for line in range(0, full_stuff.shape[0]):
        values_of_line = full_stuff.iloc[line].values
        hom_count = 0
        het_count = 0
        normal_count = 0
        for i in values_of_line:
            if i in ["1/1", "2/2", "1|1", "2|2", '1/2', '1|2', '2|1']:
                hom_count = hom_count + 1
            elif i in ['0/1', '0/2', '0|1', '0|2', '1|0', '2|0']:
                het_count = het_count + 1
            elif i == '0/0' or i == '0|0':
                normal_count = normal_count + 1
        hom_list.append(hom_count)
        het_list.append(het_count)
        norm_list.append(normal_count)

    norm_in_df = pd.Series(norm_list)
    het_in_df = pd.Series(het_list)
    hom_in_df = pd.Series(hom_list)

    full_stuff['WT count'] = norm_in_df.values
    full_stuff['Het count'] = het_in_df.values
    full_stuff['Hom count'] = hom_in_df.values

    full_stuff = pd.concat([full_stuff, healthy_df], axis=1, sort=False)
    return full_stuff, list_healthy


def help_message():
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        HIGHLIGHT = '\033[7m'
    if sys.stdout.isatty() and 'win' not in sys.platform:
        message = """
{high}---------------------------------------------------
TAPES : a Tool for Assessment and Prioritisation in  
                    Exome Studies                  
---------------------------------------------------{end}

python tapes.py function --input --output --options

{blue}---EXAMPLES---{end}                              
python tapes.py {red}sort{end} -i /path/to/input_file.txt -o /path/to/output_folder/ --options

python tapes.py {red}analyse{end} -i /path/to/input_file.csv -o /path/to/output_file.csv --options

python tapes.py {red}annotate{end} -i /path/to/input_file.vcf -o /path/to/output_file.vcf --options

python tapes.py {red}decompose{end} -i /path/to/input_file.vcf -o /path/to/output_file.vcf --options

{blue}---ANNOVAR DB MANAGEMENT---{end}
python tapes.py {red}db{end} -s -A /path/to/annovar  | See annovar databases
python tapes.py {red}db{end} -b -A /path/to/annovar  | Download annovar databases


{blue}---Main Function---{end}
{red}annotate{end}, {red}sort{end}, {red}db{end}, {red}decompose{end}, {red}analyse{end}
                    Choose the main function: 
                    annotate: use ANNOVAR to annotate a vcf file 
                    sort: use tapes to prioritise variants 
                    db: see or download databases 
                    analyse: use an already analysed file to produce a specific report

{blue}---optional arguments---{end}
-h , --help           Show this help message and exit

-i , --input          Path of the input file annotated by ANNOVAR

-o , --output         Path of the sorted output csv file

-b , --build_db       Download databases flagged for download : /humandb
                      will be populated with relevant databases Also
                      downloads the necessary tapes databases if
                      necessary

-A , --annovar        Path of the ANNOVAR FOLDER

-s , --see_db         annovar_path/humandb will be searched
                      for relevant databases and write the result in
                      db_config.json in tapes current directory

-e , --enrichr        Use --enrichr to analyse top mutations (>0.85 in Total
                      score) with EnrichR APIDefault is
                      GO_Biological_Process_2018, see
                      http://amp.pharm.mssm.edu/Enrichr/#stats to see all
                      available database

-d , --disease        Use the "Disease" Column to search for a specific term and
                      create a report, use quotes ' ' to use a multiple words term
                      (eg: "multiple sclerosis") Default = cancer

--BIG [BIG]           If the vcf is too big, exported file remove probably
                      benign variants before exporting

--by_sample           Create a report with the 5 top pathogenic variants per
                      sample 

--list                A txt file containing a list of refGene Genes symbols
                      or a string of refGenes Genes symbol in quotes

-a , --assembly [hg38 or hg19]
                      Name of the assembly, hg19 / hg38
                      Default is hg19

--ref_anno [refGene, knownGene or ensGene]
                      Genes format to use RefSeq 'refGene', ENSEMBL
                     'ensGene' or USCS 'knownGene default is refGene

--acmg                Use acmg preset databases / check presence of all
                      necessary annotations

--trio                Specify pedigree data for trio analysis as .txt file'
                      default is None

--kegg                Similar to --list. Indicate a kegg pathway and
                      tapes will return a spreasheet withgenes only
                      contained in that pathway

--cutoff              Select cutoff for BS1.Frequency for rare disease,
                      default is 0.005

--acmg_db             The location of the databases used to assign ACMG
                      criteria Default is the directory where tapes.py is
                      located

--tab                 Export in tab-separated values format instead of csv
                      format

--by_gene             Create a report with variants grouped by genes with
                      sample informations and gene burden.

-t , --threads        Numbers of threads / NOT OPTIMISED YET   

--pp2_percent         Threshold for PP2, considering PP2 positive if variant
                      is missense in a gene where more than the threshold
                      are pathogenic missense

--pp2_min             Number of minimum pathogenic missense variants in gene
                      to consider PP2

--bp1_percent         Threshold for BP1, considering BP1 positive if variant
                      is missense in a gene where less than the threshold
                      are pathogenic missense
    """.format(blue=bcolors.OKBLUE,
               green=bcolors.OKGREEN,
               red=bcolors.FAIL,
               end=bcolors.ENDC,
               head=bcolors.HEADER,
               high=bcolors.HIGHLIGHT)

    else:
        message = """
        ---------------------------------------------------
        TAPES : a Tool for Assessment and Prioritisation in
                          Exome Studies
        ---------------------------------------------------

        python tapes.py function --input --output --options

        ---EXAMPLES---                              
        python tapes.py sort -i /path/to/input_file.txt -o /path/to/output_folder/ --options

        python tapes.py analyse -i /path/to/input_file.csv -o /path/to/output_file.csv --options

        python tapes.py annotate -i /path/to/input_file.vcf -o /path/to/output_file.vcf --options

        python tapes.py decompose -i /path/to/input_file.vcf -o /path/to/output_file.vcf --options

        ---ANNOVAR DB MANAGEMENT---
        python tapes.py db -s -A /path/to/annovar  | See annovar databases
        python tapes.py db -b -A /path/to/annovar  | Download annovar databases


        ---Main Function---
        annotate, sort, db, decompose, analyse
                            Choose the main function: 
                            annotate: use ANNOVAR to annotate a vcf file 
                            sort: use tapes to prioritise variants 
                            db: see or download databases 
                            analyse: use an already analysed file to produce a specific report

        ---optional arguments---
        -h , --help           Show this help message and exit

        -i , --input          Path of the input file annotated by ANNOVAR

        -o , --output         Path of the sorted output csv file

        -b , --build_db       Download databases flagged for download : /humandb
                              will be populated with relevant databases Also
                              downloads the necessary tapes databases if
                              necessary

        -A , --annovar        Path of the ANNOVAR FOLDER

        -s , --see_db         annovar_path/humandb will be searched
                              for relevant databases and write the result in
                              db_config.json in tapes current directory

        -e , --enrichr        Use --enrichr to analyse top mutations (>0.85 in Total
                              score) with EnrichR APIDefault is
                              GO_Biological_Process_2018, see
                              http://amp.pharm.mssm.edu/Enrichr/#stats to see all
                              available database

        -d , --disease        Use the "Disease" Column to search for a specific term and
                              create a report, use quotes ' ' to use a multiple words term
                              (eg: "multiple sclerosis") Default = cancer

        --BIG [BIG]           If the vcf is too big, exported file remove probably
                              benign variants before exporting

        --by_sample           Create a report with the 5 top pathogenic variants per
                              sample 

        --list                A txt file containing a list of refGene Genes symbols
                              or a string of refGenes Genes symbol in quotes

        -a , --assembly [hg38 or hg19]
                              Name of the assembly, hg19 / hg38
                              Default is hg19

        --ref_anno [refGene, knownGene or ensGene]
                              Genes format to use RefSeq 'refGene', ENSEMBL
                             'ensGene' or USCS 'knownGene default is refGene

        --acmg                Use acmg preset databases / check presence of all
                              necessary annotations

        --trio                Specify pedigree data for trio analysis as .txt file'
                              default is None

        --kegg                Similar to --list. Indicate a kegg pathway and
                              tapes will return a spreasheet withgenes only
                              contained in that pathway

        --cutoff              Select cutoff for BS1.Frequency for rare disease,
                              default is 0.005

        --acmg_db             The location of the databases used to assign ACMG
                              criteria Default is the directory where tapes.py is
                              located

        --tab                 Export in tab-separated values format instead of csv
                              format

        --by_gene             Create a report with variants grouped by genes with
                              sample informations and gene burden.

        -t , --threads        Numbers of threads / NOT OPTIMISED YET   

        --pp2_percent         Threshold for PP2, considering PP2 positive if variant
                              is missense in a gene where more than the threshold
                              are pathogenic missense

        --pp2_min             Number of minimum pathogenic missense variants in gene
                              to consider PP2

        --bp1_percent         Threshold for BP1, considering BP1 positive if variant
                              is missense in a gene where less than the threshold
                              are pathogenic missense
            """
    return message


def flags():
    everything = [
        "TTN",
        "MUC16",
        "OBSCN",
        "AHNAK2",
        "SYNE1",
        "FLG",
        "MUC5B",
        "DNAH17",
        "PLEC",
        "DST",
        "SYNE2",
        "NEB",
        "HSPG2",
        "LAMA5",
        "AHNAK",
        "HMCN1",
        "USH2A",
        "DNAH11",
        "MACF1",
        "MUC17",
        "DNAH5",
        "GPR98",
        "FAT1",
        "PKD1",
        "MDN1",
        "RNF213",
        "RYR1",
        "DNAH2",
        "DNAH3",
        "DNAH8",
        "DNAH1",
        "DNAH9",
        "ABCA13",
        "APOB",
        "SRRM2",
        "CUBN",
        "SPTBN5",
        "PKHD1",
        "LRP2",
        "FBN3",
        "CDH23",
        "DNAH10",
        "FAT4",
        "RYR3",
        "PKHD1L1",
        "FAT2",
        "CSMD1",
        "PCNT",
        "COL6A3",
        "FRAS1",
        "FCGBP",
        "DNAH7",
        "RP1L1",
        "PCLO",
        "ZFHX3",
        "COL7A1",
        "LRP1B",
        "FAT3",
        "EPPK1",
        "VPS13C",
        "HRNR",
        "MKI67",
        "MYO15A",
        "STAB1",
        "ZAN",
        "UBR4",
        "VPS13B",
        "LAMA1",
        "XIRP2",
        "BSN",
        "KMT2C",
        "ALMS1",
        "CELSR1",
        "TG",
        "LAMA3",
        "DYNC2H1",
        "KMT2D",
        "BRCA2",
        "CMYA5",
        "SACS",
        "STAB2",
        "AKAP13",
        "UTRN",
        "VWF",
        "VPS13D",
        "ANK3",
        "FREM2",
        "PKD1L1",
        "LAMA2",
        "ABCA7",
        "LRP1",
        "ASPM",
        "MYOM2",
        "PDE4DIP",
        "TACC2",
        "MUC2",
        "TEP1",
        "HELZ2",
        "HERC2",
        "ABCA4"]
    return everything

