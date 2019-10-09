#### TEST POLYGENIC RISK SCORE ####
import requests, sys, json, os
import pandas as pd


def get_public_genotype(rs_num):
    # Request the genotype for that snp
    url = "https://grch37.rest.ensembl.org/variation/homo_sapiens/{}?genotypes=1?phenotypes=1".format(rs_num)
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.get(url, headers=headers)

    if not r.ok:
        print('|| Unable to get Genotyping information for {}'.format(rs_num))
        decoded = {}
    pass

    decoded = r.json()
    # print("Got Genotype for {}".format(rs_num))
    return decoded


def make_df_from_public_genotype(list_of_rs_num, number_of_public_samp_to_get):

    final_geno_df = pd.DataFrame()

    y = len(list_of_rs_num)
    x=0
    print('|| Retrieving public genotypes from ENSEMBL')
    for rs_num in list_of_rs_num:

        x = x + 1
        print("{}/{}".format(x,y), end='\r')

        snp_serie = pd.Series()
        whole_result = get_public_genotype(rs_num)
        genotype_list = whole_result["genotypes"]
        minor_allele = whole_result["minor_allele"]
        # Only the first 100 individuals
        #for key, _ in enumerate(genotype_list):
        if minor_allele:
            for key in range(0, number_of_public_samp_to_get):
                if genotype_list[key]['genotype'].count(minor_allele) == 0:
                    genotype_list[key]['genotype'] = '0/0'
                elif genotype_list[key]['genotype'].count(minor_allele) == 1:
                    genotype_list[key]['genotype'] = '0/1'
                elif genotype_list[key]['genotype'].count(minor_allele) == 2:
                    genotype_list[key]['genotype'] = '1/1'
                snp_serie['snpid'] = rs_num
                snp_serie[genotype_list[key]['sample']] = genotype_list[key]['genotype']
        else:
            print('|| NO GENOTPYPE found for {}'.format(rs_num))
        if len(snp_serie) > 0:
            final_geno_df = final_geno_df.append(snp_serie, ignore_index=True)

    # final_geno_df.to_csv('/home/alex/Desktop/final_geno_df.txt', sep='\t', index=False)
    # DROP NAs and Empties
    final_geno_df.replace("", None, inplace=True)
    final_geno_df.dropna(axis=1, inplace=True)

    return final_geno_df


def prs_calculation(acmg_db_path, full_stuff, output_folder, disease, number_of_public_samp, output_type):
    if output_type == "directory":
        output_folder = os.path.dirname(output_folder)
        prefix = os.path.basename(output_folder)
    else:
        prefix = os.path.basename(output_folder[:-4])
        output_folder = os.path.dirname(output_folder)
    
    df = pd.read_csv(os.path.join(acmg_db_path, 'gwas_catalog_v1.0.2-associations_e96_r2019-07-30_ORonly.tsv'), sep='\t', low_memory=False)
    # LOAD ethnicity related snps
    ethno_list = []
    with open(os.path.join(acmg_db_path, "SNP_ethnicity.txt"), 'r') as ethno_file:
        for line in ethno_file.readlines():
            line = line.strip()
            ethno_list.append(line)

    # Remove all lines containing unit (to simplify)
    df = df[~df["95% CI (TEXT)"].str.contains('unit', na=False)]
    # REMOVE ETHNO SNPS
    df = df[~df['SNPS'].isin(ethno_list)]
    # Load TAPES dataframe
    tapes_df = full_stuff
    # Remove all unwanted columns from GWAS catalogue dataframe UNNECESSARY
    only_req = df[["SNPS","DISEASE/TRAIT","BETA", ]]
    # Groupe by column disease trait
    grouped = df.groupby(["DISEASE/TRAIT"])
    # Optional Write all gwas terms to a file with the number of snps associated
    '''with open('/home/alex/Desktop/gwas_terms2.txt', 'w') as outfile:
        for item in grouped:
            if item[1].shape[0] > 5:
                outfile.write(item[0]+'\t'+str(item[1].shape[0])+'\n')'''

    # Select a disease/trait
    #disease = "Crohn's disease"
    #disease = "Neuroticism"
    #disease = "Colorectal cancer"
    #disease = "Systemic lupus erythematosus"
    #disease = "Atrial fibrillation"

    # Initialise the main dict for the disease
    dict_disease = {}
    # Create lists of snps with associated beta, gene, CI and frequency
    list_of_rs_in_disease = list(df[df['DISEASE/TRAIT'] == disease]['SNPS'])
    list_of_beta_in_disease = list(df[df['DISEASE/TRAIT'] == disease]['BETA'])
    list_of_genes_in_disease =list(df[df['DISEASE/TRAIT'] == disease]['REPORTED GENE(S)'])
    list_of_CI = list(df[df['DISEASE/TRAIT'] == disease]['95% CI (TEXT)'])
    list_of_freq = list(df[df['DISEASE/TRAIT'] == disease]['RISK ALLELE FREQUENCY'])


    # Populate main dict with the loop with snp number as key
    for rs, beta, gene, ci, freq in zip(list_of_rs_in_disease, list_of_beta_in_disease, list_of_genes_in_disease, list_of_CI, list_of_freq):
        dict_disease[rs] = {"beta": float(beta), "gene": gene, "ci": ci, "freq": freq}
        
    # GET THE SNP COLUMN 
    list_col_snp = list(tapes_df.columns)
    list_col_snp = [i for i in list_col_snp if "snp" in i]
    snp_col = list_col_snp[0]

    # INTERSECT THE TAPES DF WITH SNPS FROM DISEASE (to remove snps not sequences/genotypes in cases)
    tapes_df = tapes_df[tapes_df[snp_col].isin(list_of_rs_in_disease)]
    tapes_df.reset_index(inplace=True)
    print(len(tapes_df[snp_col]), "SNPs in common with GWAS Catalog {} disease/trait".format(disease))

    if tapes_df.shape[0] == 0:
        print('|| NO SNPS IN COMMON BETWEEN DISEASE SNP and CASES SNP')
        return

    # CREATE DF FOR PUBLIC GENOTYPES # GET SNP COLUMN
    
    public_df = make_df_from_public_genotype(list(tapes_df[snp_col]), number_of_public_samp)

    ### FIND SAMPLE COLUMNS ###
    list__samp_columns = []
    # Find sample based on at least one heterozygous sample in column
    for column in tapes_df.columns:
        if '0/1' in list(tapes_df[column]) or '0|1' in list(tapes_df[column]):
            list__samp_columns.append(column)

    # Initialise list that will store the risk score
    list_of_risks = []
    # Initialise the dict to later create the graph
    graph_dict = {}

    # Iterate over all own samples to calculate PRS
    for sample in list__samp_columns:  # FOr all samples
        samp_df = tapes_df[tapes_df[sample].isin(["0/1","0/2","1/1","1/2","2/2"])]  # Only select rows where this sample has a variant
        rs_series = samp_df[snp_col].str.split("&", n=1, expand=False).str[0]   # Get the rs numbers of snps in those variants
        zygosy_series = samp_df[sample] # create series with zygosy for the sample associated with the rs_series
        # Initialise the polygenic risk factor and the lists that will store informations
        risk_factor = 0
        list_of_genes_involved = []
        list_of_rs_involved = []
        list_of_ci_involved = []
        list_of_freq_involved = []
        # Iterate over snps to calculate PRS
        for snp_samp, zyg in zip(rs_series, zygosy_series):
            # Get the genotype from public data
            if snp_samp in list_of_rs_in_disease and (zyg == "1/1" or zyg == "1/2" or zyg == "2/2"):
                risk_factor = risk_factor + (dict_disease[snp_samp]["beta"]*2)
                list_of_rs_involved.append(snp_samp)
                list_of_genes_involved.append(dict_disease[snp_samp]["gene"])
                list_of_ci_involved.append(dict_disease[snp_samp]["ci"])
                list_of_freq_involved.append(dict_disease[snp_samp]["freq"])
            elif snp_samp in list_of_rs_in_disease and (zyg == "0/1" or zyg == "0/2" or zyg == "./1"):
                risk_factor = risk_factor + (dict_disease[snp_samp]["beta"])
                list_of_rs_involved.append(snp_samp)
                list_of_genes_involved.append(dict_disease[snp_samp]["gene"])
                list_of_ci_involved.append(dict_disease[snp_samp]["ci"])
                list_of_freq_involved.append(dict_disease[snp_samp]["freq"])
            else:
                pass
        # Populate the graph dict with risk factor for each sample
        graph_dict[sample] = risk_factor
        list_of_risks.append(risk_factor)

        # print('Polygenic Risk Factor in sample {} for {} is {} with {} SNPs involved'.format(sample, disease, risk_factor,len(list_of_genes_involved)))
        # print(list_of_ci_involved)
        # print(list_of_freq_involved)

    # PRS FOR PUBLIC GENOTYPE

    list_public_samples = list(public_df.columns[:-1])

    for sample in list_public_samples:
        rs_series = public_df['snpid']
        zygosy_series = public_df[sample]
        risk_factor = 0
        list_of_genes_involved = []
        list_of_rs_involved = []
        list_of_ci_involved = []
        list_of_freq_involved = []
        # Iterate over snps to calculate PRS
        for snp_samp, zyg in zip(rs_series, zygosy_series):
            # Get the genotype from public data
            if snp_samp in list_of_rs_in_disease and (zyg == "1/1"):
                risk_factor = risk_factor + (dict_disease[snp_samp]["beta"] * 2)
                list_of_rs_involved.append(snp_samp)
                list_of_genes_involved.append(dict_disease[snp_samp]["gene"])
                list_of_ci_involved.append(dict_disease[snp_samp]["ci"])
                list_of_freq_involved.append(dict_disease[snp_samp]["freq"])
            elif snp_samp in list_of_rs_in_disease and (zyg == "0/1"):
                risk_factor = risk_factor + (dict_disease[snp_samp]["beta"])
                list_of_rs_involved.append(snp_samp)
                list_of_genes_involved.append(dict_disease[snp_samp]["gene"])
                list_of_ci_involved.append(dict_disease[snp_samp]["ci"])
                list_of_freq_involved.append(dict_disease[snp_samp]["freq"])
            else:
                pass
        # Populate the graph dict with risk factor for each sample
        graph_dict[sample] = risk_factor
        list_of_risks.append(risk_factor)


    # print('Max is: {}'.format(max(list_of_risks)))

    # MAKE COLOR LIST

    ser = pd.Series(graph_dict)
    ser.sort_values(inplace=True)

    color_list = []
    for sample in ser.keys():
        if sample in list__samp_columns:
            color_list.append('r')
        else:
            color_list.append('b')

    # STATS
    import scipy

    mean_cases = 0
    total_cases = len(list__samp_columns)
    cases_prs_list = []
    mean_controls = 0
    total_controls = len(list_public_samples)
    controls_prs_list = []

    plot_labels = []

    for sample in ser.keys():
        if sample in list__samp_columns:
            mean_cases = mean_cases + float(ser[sample])
            cases_prs_list.append(float(ser[sample]))
            plot_labels.append(sample)
        elif sample in list_public_samples:
            mean_controls = mean_controls + float(ser[sample])
            controls_prs_list.append(float(ser[sample]))
            plot_labels.append('')
        else:
            print('THIS SHOULD NOT HAPPEN')

    mean_cases = mean_cases/total_cases
    mean_controls = mean_controls/total_controls

    stat, pval = scipy.stats.ttest_ind(cases_prs_list, controls_prs_list, equal_var=False)
    man_stat, man_pval = scipy.stats.mannwhitneyu(cases_prs_list, controls_prs_list, alternative='two-sided')


    '''print('MEAN CASES')
    print(mean_cases)
    print("MEAN CONTROLS")
    print(mean_controls)

    print('STAT T TEST')
    print(stat)
    print('P VALUE')
    print(pval)

    print('STAT MAN TEST')
    print(man_stat)
    print('MAN P VALUE')
    print(man_pval)'''

    plot_prs = ser.plot(kind='bar', figsize=(15, 12), fontsize=12, rot=55,
                                title='Polygenic Risk Score for {} (Using {} variants \n Means Case = {:.3f}; '
                                      'Mean controls = {:.5f}; p-non-parametric = {:.5f})'.format(disease,
                                                                                                  len(tapes_df[snp_col]),
                                                                                                  mean_cases, mean_controls,
                                                                                                  man_pval), color=color_list)
    plot_prs.set_xticklabels(plot_labels, fontdict={'fontsize': 6})

    #plot_prs.legend(loc='best', bbox_to_anchor=(0, 0., 0.5, 0.5), labels=ser.index)
    fig_prs = plot_prs.get_figure()
    fig_prs.savefig((os.path.join(output_folder, (prefix + '_PRS_'+ disease + '.png'))))


    # MAKE FINAL DF

    final_prs_df = pd.concat([tapes_df[list__samp_columns], public_df.reindex(tapes_df[list__samp_columns].index)], axis=1)
    final_prs_df['INFOS'] = final_prs_df['snpid'].map(dict_disease)
    final_prs_df.to_csv(os.path.join(output_folder, (prefix + '_PRS_'+ disease + '.txt')), sep='\t', index=False)
    
    with open(os.path.join(output_folder, (prefix + '_PRS_'+ disease + '.txt')), 'a') as prsfile:
        prsfile.write('Mean PRS of cases\t{}\n'.format(mean_cases))
        prsfile.write('Mean PRS of controls\t{}\n'.format(mean_controls))
        prsfile.write('p-value\t{}\n'.format(man_pval))


    return
