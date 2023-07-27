#!usr/bin/env python3
'''
    Name: tapes.py
    Author: Alexandre Xavier
    Version; 0.1
    Contact: alexandre.xavier@uon.edu.au
    Python Version: 3.x
    from https://github.com/a-xavier/tapes
    Licence: MIT
'''

import pandas as pd
import src.t_func as tf
import json
import sys
import os
import gzip


def vep_process_vcf(vcf_path, acmg_db_path):

    head_line = 0
    print(tf.tmp_stmp()+"VEP vcf processing", end=' ', flush=True)

    with open(vcf_path, 'r') as in_vcf:
        for line in in_vcf.readlines():
            if line.startswith('##') and "ID=CSQ" not in line and "ID=ANN" not in line:
                head_line = head_line + 1
            elif line.startswith('##') and "ID=CSQ" in line:
                csq_line = line
                csq_line = csq_line.replace('">', '')
                head_line = head_line + 1
            elif line.startswith('##') and "ID=ANN" in line:
                csq_line = line
                csq_line = csq_line.replace('">', '')
                head_line = head_line + 1
            else:
                break
    tp = pd.read_csv(vcf_path, sep='\t', header=head_line, iterator=True, chunksize=1000)
    df = pd.concat(tp, ignore_index=True)
    # df = pd.read_csv(vcf_path, sep='\t', header=head_line)
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
    df['CHROM'] = df['CHROM'].astype(str)
    df['CHROM'] = df['CHROM'].str.replace('chr', '')

    try:
        csq_list = csq_line.split('Format: ')[1].split('|')
        csq_list = [x.strip() for x in csq_list]
    except UnboundLocalError:
        print(' WARNING: No annotation found from VCF header. Exiting...')
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
        print(tf.tmp_stmp() + 'Only one sample found in vcf')

    ## Process info column

    info_df = df['INFO']

    info_df = info_df.str.split('CSQ=', n=1, expand=True)

    info_df = info_df[1].str.split(',', n=1, expand=True)

    info_df = info_df[0]

    info_df = info_df.str.split('|', n=-1, expand=True)

    info_df.columns = csq_list

    info_df = info_df.replace('', '.')

    ## ADD PLI IF NOT THERE
    pli_dict_path = os.path.join(acmg_db_path, "pli_dict.json")
    with open(pli_dict_path, 'r') as cf:
        pli_dict = json.load(cf)

    pli_list = []
    prec_list = []
    pnull_list = []

    for gene in info_df['SYMBOL']:
        try:
            pli_list.append(pli_dict[gene]['pli'])
            prec_list.append(pli_dict[gene]['prec'])
            pnull_list.append(pli_dict[gene]['pnull'])
        except KeyError:
            pli_list.append('.')
            prec_list.append('.')
            pnull_list.append('.')

    if str(info_df['Feature'][0]).startswith('ENST'):
        ref_anno = 'ensGene'
    elif str(info_df['Feature'][0]).startswith('NM_'):
        ref_anno = 'refGene'
    elif str(info_df['Feature'][0]).startswith('uc0'):
        ref_anno = 'knownGene'
    else:
        print('|| Cannot determine reference annotation used...')
        sys.exit(1)

    info_df['pLi'] = pli_list
    info_df['pRec'] = prec_list
    info_df['pNull'] = pnull_list

    ## ADD DISEASE Column

    dis_dict_path = os.path.join(acmg_db_path, 'Dis_dict.json')

    with open(dis_dict_path, "r") as dj:
        dis_dict = json.load(dj)

    dis_list = []
    for gene in info_df['SYMBOL']:
        try:
            dis_list.append(dis_dict[gene])
        except KeyError:
            dis_list.append('.')
    info_df['Disease_description.{}'.format(ref_anno)] = dis_list

    #info_df.to_csv('/home/alex/Desktop/info_df.txt', sep='\t', index=False)

    basic_info_df = df[['CHROM', 'POS', 'REF', 'ALT']].copy()

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


    if basic_info_df.shape[0] == info_df.shape[0]:
        full_stuff_tmp = pd.concat([basic_info_df, info_df], axis=1, sort=False)
        if full_stuff_tmp.shape[0] == geno_df.shape[0]:
            full_stuff = pd.concat([full_stuff_tmp, geno_df], axis=1, sort=False)
            print('Done')
            # FOR TESTING PURPOSE full_stuff.to_csv('/home/alex/Desktop/processedVEP.txt', sep='\t', index=False)
            return full_stuff
        else:
            print('|| Dataframes shapes not matching during vcf processing, exiting...')
            sys.exit(1)


def load_vep_databases(acmg_db_path, assembly, ref_anno, trio_data, pp2_percent, pp2_min, bp1_percent):
    acmg_db_path = tf.process_path(acmg_db_path)
    print(tf.tmp_stmp()+'Loading databases...', end=" ", flush=True)

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
    for chr, start, stop, aref, aalt in zip(PS1_df['Chr'], PS1_df['Start'], PS1_df['Stop'], PS1_df['AA_ref'],
                                            PS1_df['AA_alt']):
        PS1_dict[chr, start, stop] = {'AA_ref': aref, 'AA_alt': aalt}

    #####  PS2 wholde dict for trio ########
    if trio_data is not None:  # IF there is trio data provided
        fam_df = pd.read_csv(trio_data, sep='\t', header=None, names=['role', 'fam_id', 'samp_id'])
        # Import trio data as dataframe
        fam_id_series = fam_df['fam_id']
        family_list = fam_id_series.value_counts().index.values.tolist()  # Detect all the unique family IDs
        whole_dict = {}
        for famid_from_list in family_list:  # For each family ID do a dictionary
            tmp_df = fam_df[fam_df['fam_id'] == famid_from_list]
            tmp_df = tmp_df[['role', 'samp_id']]
            famid_from_list_dict = {
                'm': tmp_df[tmp_df['role'] == 'm']['samp_id'].to_string(index=False).strip(),
                'f': tmp_df[tmp_df['role'] == 'f']['samp_id'].to_string(index=False).strip(),
                'o': tmp_df[tmp_df['role'] == 'o']['samp_id'].to_string(index=False).strip()
            }
            whole_dict[famid_from_list] = famid_from_list_dict
    else:
        whole_dict = None

    ##### PS4 df  ########
    PS4_df = pd.read_csv(open(os.path.join(acmg_db_path, "PS4.variants.{}".format(assembly)), 'r'), sep='\t',
                         dtype={"CHR": str, 'POS_hg19': int, 'SNPID': str, 'REF': str, 'ALT': str})

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

    ##### BS2 whole thing ########

    BS2_hom_het_dict = {"het": {'1': {},
                                '2': {},
                                '3': {},
                                '4': {},
                                '5': {},
                                '6': {},
                                '7': {},
                                '8': {},
                                '9': {},
                                '10': {},
                                '11': {},
                                '12': {},
                                '13': {},
                                '14': {},
                                '15': {},
                                '16': {},
                                '17': {},
                                '18': {},
                                '19': {},
                                '20': {},
                                '21': {},
                                '22': {},
                                'X': {},
                                'Y': {},
                                'M': {},
                                'MT': {}},
                        "hom": {'1': {},
                                '2': {},
                                '3': {},
                                '4': {},
                                '5': {},
                                '6': {},
                                '7': {},
                                '8': {},
                                '9': {},
                                '10': {},
                                '11': {},
                                '12': {},
                                '13': {},
                                '14': {},
                                '15': {},
                                '16': {},
                                '17': {},
                                '18': {},
                                '19': {},
                                '20': {},
                                '21': {},
                                '22': {},
                                'X': {},
                                'Y': {},
                                'M': {},
                                'MT': {}}
                        }

    with gzip.open(os.path.join(acmg_db_path, "BS2_hom_het.{}".format(assembly)), 'rt') as bs2_file:
        for line in bs2_file.readlines():
            chr = line.split(' ')[0].strip()
            pos = line.split(' ')[1].strip()
            ref = line.split(' ')[2].strip()
            alt = line.split(' ')[3].strip()
            hom = line.split(' ')[4].strip()
            het = line.split(' ')[5].strip()
            if str(het) == '1':
                # bs2_het_list.append(str(chr) + '_' + str(pos) + '_' + ref + '_' + alt)
                BS2_hom_het_dict["het"][chr][pos] = ref + '_' + alt
            if str(hom) == '1':
                # bs2_hom_list.append(str(chr) + '_' + str(pos) + '_' + ref + '_' + alt)
                BS2_hom_het_dict["hom"][chr][pos] = ref + '_' + alt

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

    print("Done")
    
    tup_databases = (PVS1_list, 
    rek_dict, 
    PS1_dict, 
    BS2_hom_het_dict ,
    rec_list, 
    dom_list, 
    adult_list,
    PS4_df, 
    repeat_dict, 
    PP2_list, 
    BP1_list, 
    data_pm1, 
    whole_dict)
    
    return tup_databases


def check_vep_PVS1_criteria(full_stuff, PVS1_list, ref_anno, rek_dict):
    conseqence_series = full_stuff['Consequence']
    gene_series = full_stuff["SYMBOL"]
    feature_series = full_stuff['Feature']
    start_series = full_stuff['Start']

    ada_series = full_stuff['ada_score']
    rf_series = full_stuff['rf_score']
    splice_list = []
    for ada, rf in zip(ada_series, rf_series):
        try:
            spl = max(float(ada), float(rf))
            splice_list.append(spl)
        except ValueError:
            splice_list.append(float('nan'))
            pass

    PVS1_contrib = []

    for conseq, gene, splice, feat, start in zip(conseqence_series, gene_series, splice_list, feature_series, start_series):
        if gene.split(';')[0] in PVS1_list or gene in PVS1_list:  #  In case there is a ; in the gene symbol, check for first element only
            if conseq.split('&')[0] == 'stop_gained' or conseq.split('&')[0] == 'frameshift_variant':
                #TODO INCLUDE OR REMOVE VARIANTS CLOSE TO END OF THE GENE
                if ref_anno == 'refGene':
                    #exon_nums = [x.replace('exon', '') for x in aa.split(':') if 'exon' in x]
                    trans_num = feat.split('.')[0]
                # TODO Check for other reference annotation
                elif ref_anno == 'ensGene':
                    #exon_nums = [x.replace('exon', '') for x in aa.split(':') if 'exon' in x]
                    trans_num = feat
                elif ref_anno == 'knownGene':
                    #exon_nums = [x.replace('exon', '') for x in aa.split(':') if 'exon' in x]
                    trans_num = feat
                try:
                    stop = rek_dict[trans_num]
                except KeyError:
                    pass
                try:
                    if int(stop) - int(start) <= 110:  # To offset untranslated regions
                        PVS1_contrib.append(0)
                    else:
                        PVS1_contrib.append(1)
                except UnboundLocalError:
                    PVS1_contrib.append(1)
                #PVS1_contrib.append(1)

            elif (conseq.split('&')[0] == 'splice_acceptor_variant'
                  or conseq.split('&')[0] == 'splice_donor_variant') and splice >= 0.6:
                PVS1_contrib.append(1)
            else:
                PVS1_contrib.append(0)
        else:
            PVS1_contrib.append(0)
    return PVS1_contrib


def check_vep_PS1_criteria(full_stuff, PS1_dict, ref_anno):
    chr_series = full_stuff['Chr']
    ref_series = full_stuff['Ref']
    alt_series = full_stuff['Alt']
    start_series = full_stuff['Start']
    end_series = full_stuff['End']
    aa_change_series = full_stuff['Amino_acids']
    exonic_func_series = full_stuff['Consequence']
    PS1_contrib = []

    for chr, ref, alt, start, stop, aa_change, exo_func in zip(chr_series, ref_series, alt_series, start_series, end_series, aa_change_series, exonic_func_series):
        if exo_func == "missense_variant":
            if len(aa_change) == 3:# To avoid catching indels with aa changes
                first = aa_change.split('/')[0]
                last = aa_change.split('/')[1]
                try:
                    if first == PS1_dict[(chr, start, stop)]["AA_ref"] and last == PS1_dict[(chr, start, stop)]["AA_alt"]:
                        PS1_contrib.append(1)
                    else:
                        PS1_contrib.append(0)
                except KeyError:
                    PS1_contrib.append(0)
            else:   #  IF not frameshift and no aa change data available
                PS1_contrib.append(0)
        else:
            PS1_contrib.append(0)
    return PS1_contrib


def check_vep_PS2_criteria(full_stuff, whole_dict):
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
        print(tf.tmp_stmp() +
              "Could not find specified sample. Make sure that the input file includes genotyping data"
              " and that sample IDs match in the trio file. Exiting...")
        sys.exit(1)
    return PS2_total_contrib


def check_vep_PS3(full_stuff):
    review_series = full_stuff['clinvar_golden_stars']
    sig_series = full_stuff['clinvar_clnsig']

    PS3_contrib = []
    for rev, sig in zip(review_series, sig_series):
        if (rev.split('&')[0] == '3' or rev.split('&')[0] == '4') and \
                (sig.split('&')[0] == '4' or sig.split('&')[0] == '5' or sig.split('&')[0] == '6'):
            PS3_contrib.append(1)
        else:
            PS3_contrib.append(0)

    return PS3_contrib


def check_vep_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df, number_of_samples):
    return tf.check_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df, number_of_samples)


def check_vep_PM1_criteria(full_stuff, benign_domain_dict):
    PM1_list = []
    domain_series = full_stuff['Interpro_domain']
    exonic_series = full_stuff['Consequence']

    for domain, exonic in zip(domain_series, exonic_series):
        if domain.split('&')[0] in benign_domain_dict or domain == '.':  # If the domain is not present or domain is in the benign dic
            PM1_list.append(0)
        elif domain.split('&')[0] not in benign_domain_dict and 'missense_variant' in exonic:
        #elif domain not in benig_domain_dict and ('stop' in exonic or 'frameshift' in exonic or 'nonsynonymous' in exonic):
        # Not sure why but intervar excludes null variants ? look for rationale
            PM1_list.append(1)
        else:
            PM1_list.append(0)
    return PM1_list


def check_vep_PM2(full_stuff, freq_exo_series, freq_geno_series, ref_anno):
    pli_series = pd.to_numeric(full_stuff['pLi'], errors='coerce')
    prec_series = pd.to_numeric(full_stuff['pRec'], errors='coerce')
    PM2_contrib = []

    #  USE PLI SCORE FOR RECESSIVE DOMINANT CALLING
    #  Dominant / Haploinsufficient if pli score > 0.85
    #  Recessive if prec score > 0.85
    #TODO INCLUDE IN PM2 DOMINANT VARIANT WITH MAF < 10e-6 (Probably means only one het individual in the cohort)
    for pli, prec, exo, geno in zip(pli_series, prec_series, freq_exo_series, freq_geno_series):
        if (pd.isna(exo) and pd.isna(geno) and pli >= 0.85) or (exo <= 0.000001 and geno <=0.000001 and pli >= 0.85) or \
                (exo <= 0.000001 and pd.isna(geno) and pli >= 0.85) or (pd.isna(exo) and geno <=0.000001 and pli >= 0.85):  # Dominant/Haploinsufficient
            PM2_contrib.append(1)
        elif ((exo and geno <= 0.005) or (pd.isna(exo) and pd.isna(geno)) or (exo <= 0.005 and pd.isna(geno)) or
                (geno <= 0.005 and pd.isna(exo))) and (pli <= 0.85 or pd.isna(pli)):  # Recessive (no pLi requirement)

            PM2_contrib.append(1)
        else:
            PM2_contrib.append(0)
    return PM2_contrib


def check_vep_PM4_and_BP3_repeat_region(repeat_dict, full_stuff, ref_anno):
    PM4_contrib = []
    BP3_contrib = []
    chr_series = full_stuff['Chr']
    snp_start_series = full_stuff['Start']
    in_frame_series = full_stuff['Consequence']
    interpro_series = full_stuff['Interpro_domain']
    for chrom, snp_start, in_frame, domain in zip(chr_series, snp_start_series, in_frame_series, interpro_series):
        if 'inframe' in in_frame:
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
                if domain == '.' or len(
                        domain) <= 1:  # Check if known domain (either dot or any one character NA string
                    BP3_contrib.append(1)
                else:
                    BP3_contrib.append(0)
            else:  # If not in repeat region PM4 = 1 and BP3 = 0
                PM4_contrib.append(1)
                BP3_contrib.append(0)

        elif 'stop_lost' == in_frame:
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
            else:  # Else PM4 = 1 and BP3 = 0 because stoploss not nonframeshift
                PM4_contrib.append(1)
                BP3_contrib.append(0)
            # IF NOT IN REPEAT REGION and is a nonframeshift variant or a stoploss PM4=1
        else:  # If not a nonframeshift or stoploss
            PM4_contrib.append(0)
            BP3_contrib.append(0)

    return PM4_contrib, BP3_contrib


def check_vep_PM5_criteria(full_stuff, PS1_dict, ref_anno):
    chr_series = full_stuff['Chr']
    start_series = full_stuff['Start']
    end_series = full_stuff['End']
    aa_change_series = full_stuff['Amino_acids']
    exonic_func_series = full_stuff['Consequence']

    PM5_contrib = []

    for chr, start, stop, aa_change, exo_func in zip(chr_series, start_series, end_series, aa_change_series, exonic_func_series):
        if exo_func == "missense_variant":
            if len(aa_change) == 3:
                first = aa_change.split('/')[0]
                last = aa_change.split('/')[1]
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
        else:  # Exclude indels and non/frameshift substitution
            PM5_contrib.append(0)
    return PM5_contrib


def check_vep_PP2(full_stuff, PP2_list, ref_anno):
    gene_series = full_stuff['SYMBOL']
    type_series = full_stuff['Consequence']

    PP2_contrib = []
    for gene, type_snv in zip(gene_series, type_series):
        if gene in PP2_list and type_snv == 'missense_variant':
            PP2_contrib.append(1)
        else:
            PP2_contrib.append(0)
    return PP2_contrib


def check_vep_PP5(full_stuff):
    path_series = full_stuff['clinvar_clnsig']
    PP5_contrib = []
    for pred in path_series:
        if '4' in pred or '5' in pred or '6' in pred:  # TODO INCLUDE OR NOT CONFLICTING
            PP5_contrib.append(1)
        else:
            PP5_contrib.append(0)
    return PP5_contrib


def check_vep_BA1(freq_geno_series, freq_exo_series):
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


def check_vep_BS1(full_stuff, cutoff):
    # Select all frequency columns
    filter_columns = [col for col in full_stuff if col.endswith('_AF')]
    # Check the max for all populations
    max_freq_series = full_stuff[filter_columns].replace('.', 0)
    for column in max_freq_series.columns:
        max_freq_series[column] = pd.to_numeric(max_freq_series[column], errors='coerce')  # Convert each column to floats
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


def check_vep_BS2(full_stuff, BS2_hom_het_dict, rec_list, dom_list, adult_list, ref_anno):
    start_series = full_stuff['Start']
    chr_series = full_stuff['Chr']
    gene_series = full_stuff['SYMBOL']
    ref_snv = full_stuff['Ref']
    alt_snv = full_stuff['Alt']

    all_tuple = zip(chr_series, start_series, ref_snv, alt_snv, gene_series)
    BS2_contrib =[]

    for chr, start, ref, alt, gene in all_tuple:

        if gene in adult_list:  # Remove all in adult onset
            BS2_contrib.append(0)
        elif gene in rec_list:  # Check for recessive first
            key_snv = ref + '_' + alt
            try:
                if key_snv == BS2_hom_het_dict["hom"][str(chr)][str(start)]:
                    BS2_contrib.append(1)
                else:
                    BS2_contrib.append(0)
            except KeyError:
                BS2_contrib.append(0)
        elif gene in dom_list:  # Then check dominant # What if both dominant and recessive
            key_snv = ref + '_' + alt
            try:
                if key_snv == BS2_hom_het_dict["het"][str(chr)][str(start)]:
                    BS2_contrib.append(1)
                else:
                    BS2_contrib.append(0)
            except KeyError:
                BS2_contrib.append(0)
        else:
            BS2_contrib.append(0)
    return BS2_contrib


def check_vep_BS3(full_stuff):
    review_series = full_stuff['clinvar_golden_stars']
    sig_series = full_stuff['clinvar_clnsig']
    BS3_contrib = []
    for rev, sig in zip(review_series, sig_series):
        if (rev.split('&')[0] == '3' or rev.split('&')[0] == '4') and \
                (sig.split('&')[0] == '2' or sig.split('&')[0] == '3'):
            BS3_contrib.append(1)
        else:
            BS3_contrib.append(0)
    return BS3_contrib


def check_vep_BP1(full_stuff, BP1_list, ref_anno):
    gene_series = full_stuff['SYMBOL']
    type_series = full_stuff['Consequence']
    BP1_contrib = []

    for gene, type_snp in zip(gene_series, type_series):
        if gene in BP1_list and type_snp == 'missense_variant':
            BP1_contrib.append(1)
        else:
            BP1_contrib.append(0)
    return BP1_contrib


def check_vep_BP6(full_stuff):
    path_series = full_stuff['clinvar_clnsig']
    BP6_contrib = []
    for pred in path_series:
        if pred.split("&")[0] == "2" or pred.split("&")[0] == "3":
            BP6_contrib.append(1)
        else:
            BP6_contrib.append(0)
    return BP6_contrib


def check_vep_BP7(full_stuff, ref_anno):
    type_series = full_stuff['Consequence']
    scsnv_ada_string = pd.to_numeric(full_stuff['ada_score'], errors='coerce')
    scsnv_rf_string = pd.to_numeric(full_stuff['rf_score'], errors='coerce')

    scsnv_contrib =[]
    for func, ada, rf in zip(type_series, scsnv_ada_string, scsnv_rf_string):
        if func == 'synonymous_variant':
                if max(ada, rf) < 0.6:
                    scsnv_contrib.append(1)
                else:
                    scsnv_contrib.append(0)
        else:
            scsnv_contrib.append(0)
    return scsnv_contrib


def check_vep_PP3_and_BP4_dbnsfp(full_stuff):

    # STARTING TO USE THE IN SILICO PREDICTION TOOLS

    # LET'S START WITH SIFT_score

    sift_string = full_stuff["SIFT_score"]
    sift_string = sift_string.str.split("&", n=1, expand=False).str[0]
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
    mut_taster_string = mut_taster_string.str.split("&", n=1, expand=False).str[0]
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
    mut_assess_string = mut_assess_string.str.split("&", n=1, expand=False).str[0]
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
    fathmm_string = fathmm_string.str.split("&", n=1, expand=False).str[0]
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
    provean_string = provean_string.str.split("&", n=1, expand=False).str[0]
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


def process_data(full_stuff, ref_anno, number_of_samples, tup_database, cutoff):
    # Load all databases
    print(tf.tmp_stmp()+'Starting...')
    PVS1_list, rek_dict, PS1_dict, BS2_hom_het_dict, rec_list, dom_list,\
     adult_list, PS4_df, repeat_dict, PP2_list, BP1_list, data_pm1, whole_dict = tup_database

    # TODO REST INDEXES ?
    full_stuff.reset_index(drop=True, inplace=True)

    if 'chr' in str(full_stuff['Chr'].iloc[0]).lower():
        full_stuff['Chr'].replace(['chr' + s for s in tf.chr_num()], tf.chr_num(), inplace=True)
    df_length = full_stuff.shape[0]
    # TRY PVS1
    try:
    # PVS1 Null variant in a gene where LOF is known mechanism of disease
        # Make a list for splicing scores from dbscSNV

        PVS1_contrib = check_vep_PVS1_criteria(full_stuff, PVS1_list, ref_anno, rek_dict)
        print(tf.tmp_stmp()+'PVS1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        print('|| Cannnot calculate PVS1, no splicing annotation. Please annotate with dbscSNV')
        PVS1_contrib = [0] * df_length
        pass

    try:
        # PS1  AA change is the same as a known pathogenic mutation
        PS1_contrib = check_vep_PS1_criteria(full_stuff, PS1_dict, ref_anno)
        print(tf.tmp_stmp() + 'PS1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        print('|| Cannot calculare PS1')
        PS1_contrib = [0] * df_length


    # PS2 De novo (test against parents and no family histories)
    # Process PS2 data Because each family has their own one
    try:
        PS2_contrib = check_vep_PS2_criteria(full_stuff, whole_dict)
        PS2_contrib_merged = []
        ps2_length = len(PS2_contrib)
        ps2_df = tf.list_to_dataframe([x[1] for x in PS2_contrib], ['Trio_'+s+'_PS2' for s in [x[0] for x in PS2_contrib]])
        for index, row in ps2_df.iterrows():
            if 1 in list(row):
                PS2_contrib_merged.append(1)  #  IF at least one is de novo, put 1 in PS2 for this variant
            elif 1 not in list(row):
                PS2_contrib_merged.append(0)
        print(tf.tmp_stmp() + 'PS2 done')
    except (UnboundLocalError, NameError, TypeError) as e:
        PS2_contrib_merged = [0] * df_length
        print("|| No trio data, skipping PS2")

    # PS3 Well established in vitro in vivo functional studies *
    # Considering Clinvar status "reviewed by panel of expert" and "practice guidelines"

    if 'clinvar_golden_stars' in full_stuff.columns and 'clinvar_clnsig' in full_stuff.columns:
        PS3_contrib = check_vep_PS3(full_stuff)
        print(tf.tmp_stmp() + 'PS3 done')
    else:
        print('|| No "clinvar_golden_stars" or "clinvar_clnsig"  column found. Please annotate your data with a recent clinvar database')
        PS3_contrib = [0] * df_length

    # PS4 ODD ratios or database

    if 'gnomAD_exomes_AF' in full_stuff.columns and 'gnomAD_genomes_AF' in full_stuff.columns:
        freq_geno_series = pd.to_numeric(full_stuff['gnomAD_genomes_AF'], errors='coerce')
        freq_exo_series = pd.to_numeric(full_stuff['gnomAD_exomes_AF'], errors='coerce')
        PS4_contrib = check_vep_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df, number_of_samples)
        print(tf.tmp_stmp() + 'PS4 done')
    elif 'gnomAD_genomes_AF' in full_stuff.columns and 'ExAC_AF' in full_stuff.columns:
        freq_geno_series = pd.to_numeric(full_stuff['gnomAD_genomes_AF'], errors='coerce')
        freq_exo_series = pd.to_numeric(full_stuff['ExAC_AF'], errors='coerce')
        PS4_contrib = check_vep_PS4(full_stuff, freq_geno_series, freq_exo_series, PS4_df,number_of_samples)
        print(tf.tmp_stmp() + 'PS4 done')
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
        PM1_contrib = check_vep_PM1_criteria(full_stuff, data_pm1['Domain'])
        print(tf.tmp_stmp() + 'PM1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        print("|| No domain data found")
        PM1_contrib = [0] * df_length


    # PM2 Absent from controls (eg gnomad = 0) be careful about indels IN TABLE
    # Use gnomad_genome and either gnomad exome or exac for exome data
    if 'gnomAD_exomes_AF' in full_stuff.columns and 'gnomAD_genomes_AF' in full_stuff.columns:
        freq_geno_series = pd.to_numeric(full_stuff['gnomAD_genomes_AF'], errors='coerce')
        freq_exo_series = pd.to_numeric(full_stuff['gnomAD_exomes_AF'], errors='coerce')
        PM2_contrib = check_vep_PM2(full_stuff, freq_exo_series, freq_geno_series, ref_anno)
        print(tf.tmp_stmp() + 'PM2 done')
    elif 'gnomAD_genomes_AF' in full_stuff.columns and 'ExAC_AF' in full_stuff.columns:
        freq_geno_series = pd.to_numeric(full_stuff['gnomAD_genomes_AF'], errors='coerce')
        freq_exo_series = pd.to_numeric(full_stuff['ExAC_AF'], errors='coerce')
        PM2_contrib = check_vep_PM2(full_stuff, freq_exo_series, freq_geno_series, ref_anno)
        print(tf.tmp_stmp() + 'PM2 done')
    else:
        PM2_contrib = [0] * df_length

    # PM4 change of length (stoploss or nonframeshift mutation) in a non repeat region (use repeatmask file)
        # Load repeat mask dictionary
    try :
        PM4_contrib, BP3_contrib = check_vep_PM4_and_BP3_repeat_region(repeat_dict, full_stuff, ref_anno)

        print(tf.tmp_stmp() + 'PM4 and BP3 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PM4_contrib = [0] * df_length
        BP3_contrib = [0] * df_length

    # PM5 is basically PS1 with different amino acid
    try:
        PM5_contrib = check_vep_PM5_criteria(full_stuff, PS1_dict, ref_anno)
        print(tf.tmp_stmp() + 'PM5 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PM5_contrib = [0] * df_length
    # Check : IF R89H and R89L are both pathogenic and we have R89H does it still apply since R89L exist and it is
    #  already caught by PS1

    # PP2 Low rate of missense variants (a list of genes) Using Own DB
    try:
        PP2_contrib = check_vep_PP2(full_stuff, PP2_list, ref_anno)
        print(tf.tmp_stmp() + 'PP2 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PP2_contrib = [0] * df_length

    # PP5 is clinvar data

    if 'clinvar_clnsig' in full_stuff.columns:
        PP5_contrib = check_vep_PP5(full_stuff)
        print(tf.tmp_stmp() + 'PP5 done')
    else:
        PP5_contrib = [0] * df_length
    try:
        PP3_BP4_contribs = check_vep_PP3_and_BP4_dbnsfp(full_stuff)
        PP3_contrib = PP3_BP4_contribs[0]
        BP4_contrib = PP3_BP4_contribs[1]
        print(tf.tmp_stmp() + 'PP3 and BP4 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        PP3_contrib = [0] * df_length
        BP4_contrib = [0] * df_length



    try:
        BA_1_contrib = check_vep_BA1(freq_exo_series, freq_geno_series)
        print(tf.tmp_stmp() + 'BA1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BA_1_contrib = [0] * df_length

    try:
        BS1_contrib = check_vep_BS1(full_stuff, cutoff)
        print(tf.tmp_stmp() + 'BS1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BS1_contrib = [0] * df_length

    # BS2 Benign variants observed in dominant (heterozygote) and recessive (homozygote) diseases (or X linked)
    try:
        BS2_contrib = check_vep_BS2(full_stuff, BS2_hom_het_dict, rec_list, dom_list, adult_list, ref_anno)
        print(tf.tmp_stmp() + 'BS2 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BS2_contrib = [0] * df_length

    # BS3 ALSO CLINVAR with reviewed by expert panel
    try:
        BS3_contrib = check_vep_BS3(full_stuff)
        print(tf.tmp_stmp() + 'BS3 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BS3_contrib = [0] * df_length

    try:
        # BP1 Missense variant in gene where missense variants are not pathogenic
        BP1_contrib = check_vep_BP1(full_stuff, BP1_list, ref_anno)
        print(tf.tmp_stmp() + 'BP1 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BP1_contrib = [0] * df_length

    # BP6 Clinvar benign maybe we can merge BP6 and PP5 to save time ?
    # Intervar auto db is not too reliable with this one
    try:
        BP6_contrib = check_vep_BP6(full_stuff)
        print(tf.tmp_stmp() + 'BP6 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BP6_contrib = [0] * df_length

    try:
        BP7_contrib = check_vep_BP7(full_stuff, ref_anno)
        print(tf.tmp_stmp() + 'BP7 done')
    except (KeyError, UnboundLocalError, NameError) as e:
        BP7_contrib = [0] * df_length

    print(tf.tmp_stmp()+'ACMG Criteria assigned')

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

    df_all = tf.list_to_dataframe(all_criterias, names_list)
    prob_pred_list = tf.calculate_probability(df_all)
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


def rename_columns(full_stuff, ref_anno):
    dict_rename = {
        'SYMBOL': "Gene.{}".format(ref_anno),
        "BIOTYPE": "Func.{}".format(ref_anno),
        "Consequence": "ExonicFunc.{}".format(ref_anno),
    }

    full_stuff.rename(columns= dict_rename, inplace=True)

    return full_stuff
