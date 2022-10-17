#!/usr/bin/python

import os
import re
import subprocess
import pdb
import xml.etree.ElementTree as ET
import gzip
import wget
import pandas as pd
import itertools
import zipfile
from io import StringIO
import progressbar
import json
import tarfile
import scipy
from scipy import stats
import multiprocessing
from tqdm import tqdm
import numpy as np
import glob

#Set to the established download directory
dwn_dir = '/home/epineiro/Analysis/PanDrugs/2.0/downloads/'

#Set to the desired processed files directory
pro_dir = '/home/epineiro/Analysis/PanDrugs/2.0/processed/'

def process_DGIdb():

    print('Formating DGIdb file...')

    dgidb = 'DGIdb_interactions.tsv'
    dgidb_dwn = 'DGIdb_interactions_dwn.tsv'

    #Clean in the output file those lines that are matches not found and header lines introduced in each access
    outputh = open (pro_dir+dgidb, 'w')
    inputh = open (dwn_dir+dgidb_dwn, 'r')
    outputh.write('\t'.join(['gene_name','drug_name','interaction_type','source','gene_categories'])+'\n')
    for line in inputh:
        line = line.rstrip('\n')
        line_a = line.split('\t')
        if re.search('Possible suggestions:',line_a[0]) == None and re.search('Unmatched search term:',line_a[0]) == None and re.search('gene_name',line_a[0]) == None:
            if line_a[3] not in ['CIViC', 'OncoKB']:
                outputh.write(line+'\n')
    inputh.close()
    outputh.close()

def process_sabdab():

    print('Formating SAbDab file...')

    sabdab_file = 'sabdab.tsv'
    sabdab_dwn = 'TheraSAbDab_SeqStruc_OnlineDownload.csv'

    sabdab_targets = {'Anthrax Protective Antigen': '', 'APP Abeta 1-40/42;TFRC': 'APP;TFRC', 'AR-C124910XX': '', 'CALCA&CALCB': 'CALCA;CALCB', 'Canine NGFB': '', 'CD41_CD61': 'CD41;CD61', 'CEACAM5&CD3E;CD3E': 'CEACAM5;CD3E', 'Dabigatran': '', 'DNA/H1 Complex': '', 'ERBB2 (Domain II);ERBB2 (Domain IV)': 'ERBB2', 'FN extra domain B': 'FN1', 'FZD Family': 'FZD1;FZD2;FZD3;FZD4;FZD5;FZD6;FZD7;FZD8;FZD9;FZD10', 'Ganglioside GD2': 'B4GALNT1', 'Heat Shock Protein 90 Homolog': 'HSP90', 'HHV gB AD-1': '', 'HIV-1 gp120': '', 'HIV-1 gp120 CD4bs': '', 'HIV-1 gp120 V3': '', 'Idiotope of anti-NeuGc-ganliosides': '', 'IGF1&IGF2': 'IGF1;IGF2', 'IL31 (Canine)': '', 'IL31 (Feline)': '', 'Influenza A HA': '', 'Influenza A HA2': '', 'Influenza B Virus': '', 'ITGA2B_ITGB3': 'ITGA2B;ITGB3', 'ITGA4_ITGB7': 'ITGA4;ITGB7', 'ITGA4_ITGB7&ITGAE_ITGB7': 'ITGA4;ITGAE;ITGB7', 'ITGAV_ITGB3': 'ITGAV;ITGB3', 'MRSV envelope protein': '', 'NGFB (Canine)': '', 'Non-Binding': '', 'NOTCH2&3': 'NOTCH2;NOTCH3', 'PcrV type III secretion system': '', 'PcrV type III secretion system;Polysaccharide synthesis locus (Pseudomonas)': '', 'PDCD1 (Canine)': '', 'Phosphatidylserine': 'PTDSS1', 'pro-GDF-8': 'MSTN', 'Rabies Virus Spike Glycoprotein G': '', 'Rabies Virus Strain ERA GP Ectodomain Epitope G-II': '', 'Rabies Virus Strain ERA GP Ectodomain Epitope G-III': '', 'Rabies Virus Surface Glycoprotein 4 (gp4) Epitope 1': '', 'RhD': 'RHD', 'RSV gpF': '', 'RSV gpF;RSV gpF': '', 'RV Antigenic Site III': '', 'SARS-CoV-2 Spike': '', 'SARS-CoV-2 Spike RBD': '', 'Serotype IATS O11': '', 'Shiga Toxin Type 1': '', 'Shiga Toxin Type 2': '', 'SIRPα': 'SIRPA', 'TGFB1 (Canine)': '', 'Toxin A': '', 'VP2 (Canine)': '', 'Zaire Ebolavirus GP': '', 'Zaire Ebolavirus GP1': '', 'α4β7': 'ITGB7', 'RSV': '', 'SpA': ''}

    sabdab = pd.read_csv(dwn_dir+sabdab_dwn, low_memory=False)
    sabdab.fillna('', inplace=True)

    for index, row in sabdab.iterrows():
        if row['Target'] in sabdab_targets.keys():
            sabdab.loc[index, 'Target'] = sabdab_targets[row['Target']]

    sabdab = sabdab.rename(columns={'Therapeutic': 'drug_name', 'Target': 'gene_name'})
    sabdab = sabdab.drop(sabdab.loc[sabdab['gene_name'] == ''].index)
    sabdab.to_csv(pro_dir+sabdab_file, index=False, sep='\t', header=True)

def process_moalmanac():

    print('Formating MOAlmanac file...')

    moalmanac_file = 'moalmanac.tsv'
    moalmanac_dwn = 'moalmanac_dwn.tsv'

    def process_alterations(selection):
        alterations = []
        for ele in list(set(selection['feature_type'].tolist())):
            alterations.append(ele+' ('+','.join(list(set(selection[selection['feature_type'] == ele]['protein_change'].dropna().tolist())))+')')
        return(';'.join(alterations))

    inputf = pd.read_csv(dwn_dir+moalmanac_dwn, sep ='\t', header=None, low_memory=False)
    inputf.columns = ['assertion_id', 'feature_type', 'gene', 'gene1', 'gene2', 'drug', 'resistance', 'sensitivity', 'therapy_type', 'protein_change', 'predictive_implication', 'validated']
    inputf = inputf[inputf['validated'] == True]
    inputf = inputf[((inputf['drug'].str.contains('\+') == False) & (inputf['predictive_implication'].isin(['FDA-Approved', 'Guideline', 'Clinical trial', 'Clinical evidence']))) | ((inputf['drug'].str.contains('\+')) & (inputf['predictive_implication'].isin(['FDA-Approved', 'Guideline'])))]

    subset_add = pd.DataFrame()
    for i in ('gene', 'gene1', 'gene2'):
        subset_add1 = inputf[[i,'drug']]
        subset_add1.columns = ['gene', 'drug']
        subset_add = subset_add.append(subset_add1[~subset_add1.isnull().any(axis=1)])

    subset_add = subset_add.drop_duplicates()

    moalmanac = pd.DataFrame(columns=['gene_name','drug_name','source','response','alteration','therapy_type'])
    source = 'MOAlmanac'
    for index, row in subset_add.iterrows():
        selection = inputf[((inputf['gene'] == row['gene']) | (inputf['gene1'] == row['gene']) | (inputf['gene2'] == row['gene'])) & (inputf['drug'] == row['drug'])]

        if len(selection) > 0:
            drug = row['drug']
            if re.search('\+', row['drug']):
               drug = ' + '.join(sorted(drug.split(' + ')))

            response = []
            if True in list(set(selection['sensitivity'].tolist())):
                response.append('Sensitivity')
                alteration_sen = process_alterations(selection[selection['sensitivity'] == True])

            if True in list(set(selection['resistance'].tolist())):
                response.append('Resistance')
                alteration_res = process_alterations(selection[selection['resistance'] == True])

            if len(response) > 0:
                if len(response) == 2:
                   alteration = alteration_sen+'/'+alteration_res
                elif response[0] == 'Sensitivity':
                   alteration = alteration_sen
                elif response[0] == 'Resistance':
                    alteration = alteration_res

                therapy = list(set(selection['therapy_type'].tolist()))[0]
                moalmanac = moalmanac.append({'gene_name' : row['gene'], 'drug_name' : drug, 'source' : source, 'response' : '/'.join(response), 'alteration' : alteration, 'therapy_type' : therapy} , ignore_index=True)

    moalmanac.to_csv(pro_dir+moalmanac_file, index=False, sep='\t', header=True)

def process_GDSC():

    print('Formating GDSC...')

    GDSC = 'GDSC.tsv'
    GDSC_dwn = 'GDSC_features.csv'
    GDSC_tsv = 'GDSC_features.tsv'

    #Obtain the name of the ANOVA downloaded file from GDSC
    ANOVA_dwn = ''
    listdir = os.listdir(dwn_dir)
    for efile in listdir:
        if re.search('^ANOVA', efile) != None: ANOVA_dwn = efile

    filei = open(dwn_dir+GDSC_dwn,'r')
    fileo = open(pro_dir+GDSC_tsv,'w')

    for line in filei:
        line = line.strip('\n')
        line_a = line.split(',')
        if line_a[0] == 'cell_line_name':
            fileo.write('\t'.join(line_a)+'\n')
        else:
            fileo.write('\t'.join(line_a[0:8])+'\t'+','.join(line_a[8:])+'\n')

    filei.close()
    filei.close()

    def label_sens(row):
        if row['feature_delta_mean_ic50'] < 0: return 'sensitivity'
        else: return 'resistance'

    def obtain_features(row):
        if re.search('_mut',row['feature_name']) != None:
            return row['feature_name'].replace('_mut','').replace('.',',')
        elif re.search('loss.|gain.',row['feature_name']) != None and re.search('\.\.',row['feature_name']) != None:
            field_a = row['feature_name'].split('.')

            if re.search('loss|gain',field_a[3]) != None:
                pan_list = []
                for ele in row['feature_name'].split('.'):
                    if re.search('cnaPANCAN',ele): pan_list.append(ele)
                field_a = list(set([item for sublist in [y.split(',') for y in [str(x) for x in gfeatures[(gfeatures['genetic_feature'].isin(pan_list)) & (gfeatures['is_mutated'] == 1)]['genes_in_segment']]] for item in sublist]))
                return ','.join(field_a)
            else:
                #Correction of NKX2.1 value
                if row['feature_name'] == "gain.cnaPANCAN393..FOXA1.NKX2.1.PSMA6.":
                    field_a[field_a.index('1')-1:field_a.index('1')+1] = ['-'.join(field_a[field_a.index('1')-1:field_a.index('1')+1])]
                return ','.join(field_a[3:len(field_a)-1])

        elif re.search('loss.|gain.',row['feature_name']) != None and re.search('\.\.',row['feature_name']) == None:
            field_a = list(set([item for sublist in [y.split(',') for y in [str(x) for x in gfeatures[(gfeatures['genetic_feature'] == row['feature_name'].split('.')[1]) & (gfeatures['is_mutated'] == 1)]['genes_in_segment']]] for item in sublist]))

            while 'UNCLASSIFIED' in field_a: field_a.remove('UNCLASSIFIED')
            return ','.join(field_a)

    def obtain_alt(row):
        if re.search('loss.',row['feature_name']) != None: return 'deletion'
        elif re.search('gain.',row['feature_name']) != None: return 'amplification'
        elif re.search('_mut',row['feature_name']) != None: return 'mutation'

    gfeatures = pd.read_csv(pro_dir+'GDSC_features.tsv', sep='\t', header = 0, dtype={'cell_line_name': str, 'cosmic_sample_id': int, 'gdsc_desc1': str, 'gdsc_desc2': str, 'tcga_desc': str, 'genetic_feature': str, 'is_mutated': int, 'recurrent_gain_loss': str, 'genes_in_segment': str})

    filei = pd.read_excel(dwn_dir+ANOVA_dwn, sheet_name = 'PANCANCER_ANOVA', header = 0, dtype={'drug_name': str, 'drug_id': int, 'target': str, 'target_pathway': str, 'feature_name': str, 'n_feature_pos': int, 'n_feature_neg': int, 'ic50_effect_size': float, 'log_ic50_mean_pos': float, 'log_ic50_mean_neg': float, 'log_max_conc_tested': float, 'feature_ic50_t_pval': float, 'feature_delta_mean_ic50': float, 'feature_pos_ic50_var': float, 'feature_neg_ic50_var': float, 'feature_pval': float, 'tissue_pval': float, 'msi_pval': float, 'fdr': float, 'tissue_type': str, 'dataset_version': float})

    filei = filei[(filei['feature_pval'] < 0.001) & (filei['fdr'] < 25)]

    filei['response'] = filei.apply(lambda row: label_sens(row), axis = 1)
    filei['gene'] = filei.apply(lambda row: obtain_features(row), axis = 1)
    filei['alteration'] = filei.apply(lambda row: obtain_alt(row), axis = 1)
    filei['source'] = 'GDSC'

    filei.rename(columns={'gene': 'gene_name'}, inplace=True)

    filei.to_csv(pro_dir+GDSC, sep='\t', index=False, header=True)

def process_KEGG_ATC():

    print('Formating KEGG ATC...')

    KEGG_ATC = 'TargetBasedClassificationKEGG_formated.tsv'
    KEGG_ATC_dwn = 'br08310.keg'

    filei = open(dwn_dir+KEGG_ATC_dwn, 'r')
    fileo = open(pro_dir+KEGG_ATC, 'w')

    fileo.write('\t'.join(['drug', 'familiy1', 'family2', 'family3'])+'\n')

    elementA = ''
    elementB = ''
    elementC = ''
    elementD = ''

    for line in filei:
        line = line.rstrip()
        line = line.replace('<b>', ' ')
        line = line.replace('</b>', ' ')
        line_a = line.split('\t')

        if re.search('^[ABCDEF] *D[\d]{5}', line_a[0]):
            pattern = re.compile('^[ABCDEF] *D[\d]{5} +(.+?) \(')
            name_s = pattern.search(line_a[0])
            if name_s == None:
                pattern = re.compile('^[ABCDEF] *D[\d]{5} +(.+?)$')
                name_s = pattern.search(line_a[0])
            name = name_s.group(1)
            tipo = ''
            if len(line_a) > 1:
                elementA_w = elementA + '_'+line_a[1]
                elementB_w = elementB + '_'+line_a[1]
                elementC_w = elementC + '_'+line_a[1]
                elementD_w = elementD + '_'+line_a[1]
            else:
                (elementA_w, elementB_w, elementC_w, elementD_w) = (elementA, elementB, elementC, elementD)

            if line_a[0][0] == 'C': elementB_w = ''
            if line_a[0][0] == 'D': elementC_w = ''
            if line_a[0][0] == 'E': elementD_w = ''


            fileo.write('\t'.join([name,elementA_w,elementB_w,elementC_w,elementD_w])+'\n')
        else:
            if re.search('^A', line_a[0]):

                pattern = re.compile('^A (.+)')
                elementA_s = pattern.search(line_a[0])
                elementA = elementA_s.group(1)
                elementB = ''
                elementC = ''
                elementD = ''
            elif re.search('^B', line_a[0]):

                pattern = re.compile('^B  (.+) \[')
                elementB_s = pattern.search(line_a[0])
                if elementB_s == None:
                    pattern = re.compile('^B  (.+)$')
                    elementB_s = pattern.search(line_a[0])
                elementB = elementB_s.group(1)
                elementC = ''
                elementD = ''
            elif re.search('^C', line_a[0]):

                pattern = re.compile('^C    (.+?) \[')
                elementC_s = pattern.search(line_a[0])
                if elementC_s == None:
                    pattern = re.compile('^C    (.+?)$')
                    elementC_s = pattern.search(line_a[0])
                elementC = elementC_s.group(1)
                elementD = ''
            elif re.search('^D', line_a[0]):

                pattern = re.compile('^D      (.+?) \[')
                elementD_s = pattern.search(line_a[0])
                if elementD_s == None:
                    pattern = re.compile('^D      (.+?)$')
                    elementD_s = pattern.search(line_a[0])
                elementD = elementD_s.group(1)

    filei.close()
    fileo.close()

def process_cmap():

    print('Formating CMAP MOA...')

    CMAP_MOA = 'cmap_moa.tsv'

    fileo = open(pro_dir+CMAP_MOA ,'w')
    fileo.write('\t'.join(['drug', 'moa'])+'\n')

    with open(dwn_dir+CMAP_MOA) as f:
        for line in f:
            data = json.loads(line)
            for e in data:
                fileo.write('\t'.join([e['pert_iname'], e['name']])+'\n')
    fileo.close()

def process_FDA():

    print('Formating FDA file...')

    fda_status_file = 'fda_status.tsv'
    #Obtain the name of the FDA zip downloaded
    fda_dwn = ''
    listdir = os.listdir(dwn_dir)
    for efile in listdir:
        if re.search('^drugsatfda', efile) != None: fda_dwn = efile

    df_dict = dict()
    with zipfile.ZipFile(dwn_dir+fda_dwn) as z:
        for filename in z.namelist():
            if filename == 'MarketingStatus.txt':
                # read the file
                with z.open(filename) as f:
                    data = pd.read_csv(f, sep ='\t', header=0, low_memory=False)
                    df_dict['markstat'] = data
            elif filename == 'Products.txt':
                with z.open(filename) as f:
                    content = f.read()
                    content = str(content, 'utf-8').replace('\t\r\n', '\r\n')
                    data = pd.read_csv(StringIO(content), sep ='\t', header=0, low_memory=False)
                    df_dict['products'] = data
            elif filename == 'MarketingStatus_Lookup.txt':
                with z.open(filename) as f:
                    data = pd.read_csv(f, sep ='\t', header=0, low_memory=False)
                    df_dict['markstatlook'] = data

    #create the output file
    fda_drug_list = {}
    for index, row in df_dict['products'].iterrows():
        compound = row['ActiveIngredient'].split('; ')
        applno = row['ApplNo']
        for comp in compound:
            status = df_dict['markstatlook'][df_dict['markstatlook']['MarketingStatusID'] == df_dict['markstat'][df_dict['markstat']['ApplNo'] == applno]['MarketingStatusID'].tolist()[0]]['MarketingStatusDescription'].tolist()
            if comp in fda_drug_list.keys():
                fda_drug_list[comp][0] = fda_drug_list[comp][0]+status
                fda_drug_list[comp][1] = fda_drug_list[comp][1]+[applno]
            else:
                fda_drug_list[comp] = [status, [applno]]

    for comp in fda_drug_list:
        if 'Prescription' in fda_drug_list[comp][0] or 'Over-the-counter' in fda_drug_list[comp][0]: fda_drug_list[comp][0] = 'Approved'
        elif 'Discontinued' in fda_drug_list[comp][0]: fda_drug_list[comp][0] = 'Withdrawn'
        else: fda_drug_list[comp][0] = ''
        fda_drug_list[comp][1] = ';'.join([str(x) for x in list(set(fda_drug_list[comp][1]))])

    fda_status = pd.DataFrame(list(fda_drug_list.items()), columns=['drug','data'])
    fda_status = pd.concat([fda_status, pd.DataFrame(fda_status['data'].to_list(), columns = ['status', 'ApplNo'])], axis=1)
    fda_status = fda_status.drop('data', axis=1)
    fda_status.to_csv(pro_dir+fda_status_file, index=False, sep='\t')

def process_FDA_label():

    print('Formating FDA label file...')

    FDA_label = 'fda_labels.tsv'

    fileo = open(pro_dir+FDA_label,'w')
    fileo.write("\t".join(['indication', 'ApplNo'])+'\n')

    #Obtain the name of the FDA zip downloaded
    listdir = os.listdir(dwn_dir)

    for efile in listdir:

        if re.search('^drug\-label', efile) != None:

            with zipfile.ZipFile(dwn_dir+efile) as z:
                for filename in z.namelist():
                    # read the file
                    with z.open(filename) as f:
                        data = json.load(f)
                        for e in data['results']:
                            if 'indications_and_usage' in e.keys() and 'application_number' in e['openfda'].keys():
                                fileo.write("\t".join([e['indications_and_usage'][0], e['openfda']['application_number'][0]])+'\n')

    fileo.close()

def process_EMA():

    print('Formating EMA file...')

    ema_file = 'ema_status.tsv'
    EMA_dwn = 'Medicines_output_european_public_assessment_reports.xlsx'

    filei = pd.read_excel (dwn_dir+EMA_dwn)
    filei.columns = filei.loc[filei.iloc[:,0] == 'Category'].values.tolist()[0]
    selection = filei.loc[filei['Category'] == 'Human'][['Medicine name', 'Therapeutic area', 'International non-proprietary name (INN) / common name', 'Active substance', 'Authorisation status', 'Condition / indication']]

    selection.to_csv(pro_dir+ema_file, index=False, sep='\t')

def process_ct():

    print('Formating Clinical Trial file...')

    ct = 'clinicaltrials.tsv'
    ct_dwn = 'AllPublicXML.zip'

    fileo = open(pro_dir+ct,'w')
    fileo.write("\t".join(['condition', 'title', 'status', 'drug'])+'\n')

    with zipfile.ZipFile(dwn_dir+ct_dwn) as z:
        for i in progressbar.progressbar(range(len(z.infolist()))):
            if re.search('.xml$', z.infolist()[i].filename):
                with z.open(z.infolist()[i].filename) as f:
                    [title,status,condition,drug] = ['','','','']
                    tree = ET.ElementTree(file=f)
                    for elem in tree.iterfind('official_title'):
                        title = elem.text
                    for elem in tree.iterfind('overall_status'):
                        status = elem.text
                    for elem in tree.iterfind('condition'):
                        condition = elem.text
                    for elem in tree.iterfind('intervention/intervention_type'):
                        if elem.text == 'Drug':
                            for elem in tree.iterfind('intervention/intervention_name'):
                                drug = elem.text

                    if drug != '': fileo.write('\t'.join([condition,title,status,drug])+'\n')

    fileo.close()

def process_cgc_oncovar():

    print('Creating file oncovar scores...')

    oncovar_file = 'oncovar_scores.tsv'
    oncovar_dwn = 'TCGA.PanCancer.all.genes.OncoVar.tsv.gz'

    oncovar_tumor_dir = 'All_genes_OncoVar_TCGA.tar.gz'
    oncovar = ''

    tar = tarfile.open(name=dwn_dir+oncovar_tumor_dir)

    output = open(pro_dir+oncovar_file, 'w')
    output.write('\t'.join(['gene', 'cancer_type', 'P_value', 'FDR', 'OncoVar_Score', 'Consensus_Score', 'Driver_Level'])+'\n')

    for member in tar.getnames():
        if re.search('.gz$', member):
            f=tar.extractfile(member)
            data = pd.read_csv(f, compression='gzip', sep='\t', low_memory=False, header=0, encoding = 'ISO-8859-1')
            for index, row in data.iterrows():
                if row['Cancer'] == 'PanCancer':
                    cancer = row['Cancer']
                else:
                    cancer = row['Cancer'].split('[')[1]
                    cancer = cancer.replace(']', '')
                output.write('\t'.join([row['Gene_symbol'], cancer, str(row['P_value']), str(row['FDR']), str(row['OncoVar_Score']), str(row['Consensus_Score']), str(row['Driver Level'])])+'\n')
            if re.search(oncovar_dwn, member):
                oncovar = data

    print('Creating file with ONC and TSG classification...')

    cgc_file = 'cancer_gene_census.csv'
    cgc = pd.read_csv(dwn_dir+cgc_file, low_memory=False, header=0, encoding = 'ISO-8859-1')
    cgc = cgc.fillna('')

    fileo = open(pro_dir+'generole.tsv','w')
    fileo.write("\t".join(['gene', 'CGC', 'OncoKB', '2020Rule', 'CTAT', 'Oncogene', 'TSgene']) + '\n')

    genes = list(set(cgc['Gene Symbol'].tolist() + oncovar['Gene_symbol'].tolist()))

    for gene in genes:

        categories = [cgc.loc[cgc['Gene Symbol'] == gene]['Role in Cancer'], oncovar.loc[oncovar['Gene_symbol'] == gene]['OncoKB'], oncovar.loc[oncovar['Gene_symbol'] == gene]['2020Rule'], oncovar.loc[oncovar['Gene_symbol'] == gene]['CTAT'], oncovar.loc[oncovar['Gene_symbol'] == gene]['Oncogene'], oncovar.loc[oncovar['Gene_symbol'] == gene]['TSgene']]

        for idx, cat in enumerate(categories):
            list_cat = []
            if not cat.empty:
                value = cat.tolist()[0]
                if isinstance(value, str):
                    subcats = value.split(', ')
                    subcats = [x.split('/') for x in subcats]
                    subcats = [item for sublist in subcats for item in sublist]

                    for subcat in subcats:
                        if subcat not in list_cat: list_cat.append(subcat.upper())

                    if idx == 4:
                        if list_cat[0] == 'Y': list_cat[0] = 'ONC'
                    elif idx == 5:
                        if list_cat[0] == 'Y': list_cat[0] = 'TSG'

            list_cat = ['ONC' if x in ['FUSION', 'ONCOGENE'] else x for x in list_cat]
            list_cat = ['' if x in ['N', 'Y'] else x for x in list_cat]

            categories[idx] = '/'.join(list(set(list_cat)))

        fileo.write("\t".join([gene]+categories) + '\n')

    fileo.close()

    print('Creating file with CGC info for scores...')

    cgc = 'cgc_scores.tsv'

    fileo = open(pro_dir+cgc,'w')
    fileo.write("\t".join(['gene', 'cancer_type_source', 'Tier']) + '\n')

    for index, row in cgc.iterrows():
        tumors = row['Tumour Types(Somatic)'].split(', ')+row['Tumour Types(Germline)'].split(', ')
        tumors = [x for x in tumors if x != '']
        for t in list(set(tumors)):
            fileo.write("\t".join([row['Gene Symbol'], t, str(row['Tier'])]) + '\n')

    fileo.close()

def process_KEGG_ind():

    print('Creating pathway member information...')

    paths = ['hsa03320', 'hsa04010', 'hsa04012', 'hsa04014', 'hsa04015', 'hsa04020', 'hsa04022', 'hsa04024', 'hsa04062', 'hsa04064', 'hsa04066', 'hsa04068', 'hsa04071', 'hsa04110', 'hsa04114', 'hsa04115', 'hsa04150', 'hsa04151', 'hsa04152', 'hsa04210', 'hsa04261', 'hsa04270', 'hsa04310', 'hsa04330', 'hsa04340', 'hsa04350', 'hsa04370', 'hsa04390', 'hsa04510', 'hsa04520', 'hsa04530', 'hsa04540', 'hsa04611', 'hsa04620', 'hsa04621', 'hsa04622', 'hsa04630', 'hsa04650', 'hsa04660', 'hsa04662', 'hsa04664', 'hsa04666', 'hsa04668', 'hsa04670', 'hsa04722', 'hsa04910', 'hsa04912', 'hsa04914', 'hsa04915', 'hsa04916', 'hsa04919', 'hsa04920', 'hsa04921', 'hsa04922', 'hsa04971', 'hsa05010', 'hsa05012', 'hsa05160', 'hsa05200', 'hsa05205', 'hsa05212', 'hsa05214', 'hsa05218', 'hsa05231']

    global dwn_dir
    global pro_dir

    dwn_dir_KEGG = dwn_dir+'KEGGmodeled/'
    pro_dir_KEGG = pro_dir+'KEGGmodeled'

    if not os.path.exists(pro_dir_KEGG): os.makedirs(pro_dir_KEGG)

    if not os.path.exists(pro_dir_KEGG+'/'+'components'): os.makedirs(pro_dir_KEGG+'/'+'components')
    if not os.path.exists(pro_dir_KEGG+'/'+'relations'): os.makedirs(pro_dir_KEGG+'/'+'relations')
    if not os.path.exists(pro_dir_KEGG+'/'+'subpaths'): os.makedirs(pro_dir_KEGG+'/'+'subpaths')
    if not os.path.exists(pro_dir_KEGG+'/'+'dictionaries'): os.makedirs(pro_dir_KEGG+'/'+'dictionaries')

    pro_dir_KEGG = pro_dir_KEGG+'/'

    gene_symbol_dic = {}

    def retrieve_symbol(gene):

        if gene not in gene_symbol_dic.keys():

            import sys
            from Bio import Entrez

            Entrez.email = "" #place email

            request = Entrez.epost("gene", id=gene)

            result = Entrez.read(request)

            webEnv = result["WebEnv"]
            queryKey = result["QueryKey"]
            data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
            annotations = Entrez.read(data)

            gene_symbol = annotations["DocumentSummarySet"]["DocumentSummary"][0]["NomenclatureSymbol"]

            gene_symbol_dic[gene] = gene_symbol

        else: gene_symbol = gene_symbol_dic[gene]

        return gene_symbol

    def parse_kegg_xml():

        for path in paths:

            output_components = open(pro_dir_KEGG+'components/'+path+'.components.tsv','w')
            output_relations = open(pro_dir_KEGG+'relations/'+path+'.relations.tsv','w')

            tree = ET.parse(dwn_dir_KEGG+path+'.kgml.xml')
            root = tree.getroot()
#            print (root.tag)
#            print (root.attrib)
#            for child in root:
#                print(child.tag, child.attrib)
#                pdb.set_trace()
            group = {}
            for entry in root.iter('entry'):
                if entry.attrib['type'] not in ['ortholog', 'map']:
                    ID = entry.attrib['id']
                    if entry.attrib['type'] == 'group':
                        group[ID] = []
                        for component in entry.iter('component'):
                            group[ID].append(component.attrib['id'])
                    else:
                        components = entry.attrib['name'].replace('hsa:','')
                        components = components.replace('cpd:','')
                        components = components.replace('dr:','')
                        components = components.replace('gl:','').split(' ')
                        output_components.write('\t'.join([ID,':'.join(components)])+'\n')

            for relation in root.iter('relation'):
                entry1 = [relation.attrib['entry1']]
                entry2 = [relation.attrib['entry2']]
                TYPE = relation.attrib['type']

                #if ID belongs to a group the relationship is repited for each component and the name will be the combined one
                if entry1[0] in group.keys():
                    entry1 = group[entry1[0]]

                if entry2[0] in group.keys():
                    entry2 = group[entry2[0]]

                r_product = list(itertools.product(entry1, entry2))

                (name, value) = ([], [])
                for subtype in relation.iter('subtype'):
                    name.append(subtype.attrib['name'])
                    value.append(subtype.attrib['value'])

                    for r in r_product:
                        output_relations.write('\t'.join([r[0],r[1],TYPE,','.join(name),','.join(value)])+'\n')

            output_components.close()
            output_relations.close()

    def nodes_ids_and_symbols():
        for i in progressbar.progressbar(range(len(paths))):
            relations = pd.read_csv(pro_dir_KEGG+'relations/'+paths[i]+'.relations.tsv', sep='\t', low_memory=False, header=None)
            relations.columns = ['entry1', 'entry2', 'type', 'name', 'value']
            components = pd.read_csv(pro_dir_KEGG+'components/'+paths[i]+'.components.tsv', sep='\t', low_memory=False, header=None)
            components.columns = ['ID', 'components']

            fileo = open(pro_dir_KEGG+'dictionaries/'+paths[i]+'.codes.tsv', 'w')
            fileo.write('\t'.join(['node', 'id', 'symbol'])+'\n')

            nodes = list(set(relations['entry1'].tolist() +relations['entry2'].tolist()))

            for node in nodes:
                components_id = components.loc[components['ID'] == node]['components'].tolist()[0]
                if re.search('\.', components_id):
                    ids = components.loc[components['ID'] == entry]['components'].tolist()[0].split('.')
                    for i, t in enumerate(transl):
                        transl[i] = components.loc[components['ID'] == int(t)]['components'].tolist()[0]
                else:
                    ids = components.loc[components['ID'] == node]['components'].tolist()

                ids = [x for ele in ids for x in ele.split(':')]

                symbols = []
                for el in ids:
                    if re.search('(C|G)', el) == None:
                        symbol = retrieve_symbol(el)
                        if symbol != '': symbols.append(symbol)

                fileo.write('\t'.join([str(node), ':'.join(ids),':'.join(symbols)])+'\n')

            fileo.close()

    def translate_codes_to_symbols():
        for i in progressbar.progressbar(range(len(paths))):
            relations = pd.read_csv(pro_dir_KEGG+'relations/'+paths[i]+'.elemcode.tsv', sep='\t', low_memory=False, header=None)
            relations.columns = ['entry1', 'entry2', 'type', 'name', 'value']
            fileo = open(pro_dir_KEGG+'relations/'+paths[i]+'.symbol.tsv', 'w')

            for index, row in relations.iterrows():
                entries = [row['entry1'], row['entry2']]
                for idx, entry in enumerate(entries):
                    transl = []
                    for el in entry.split(':'):
                        if re.search('(C|G)', el) == None:
                            symbol = retrieve_symbol(el)
                            transl.append(symbol)
                        else:
                            transl.append('')

                    entries[idx] = ':'.join(transl)

                fileo.write('\t'.join(entries+row[2:].tolist())+'\n')

            fileo.close()

    def obtain_subroutes():
        import networkx as nx
        for i in progressbar.progressbar(range(len(paths))):
            fileo = open(pro_dir_KEGG+'subpaths/'+paths[i]+'.symbol.tsv','w')
            filei = pd.read_csv(pro_dir_KEGG+'relations/'+paths[i]+'.relations.tsv', sep='\t', low_memory=False, header=None)
            filei.columns = ['entry1', 'entry2', 'type', 'name', 'value']
            dictionary = pd.read_csv(pro_dir_KEGG+'dictionaries/'+paths[i]+'.codes.tsv', sep='\t', low_memory=False, header=0)
            generole = pd.read_csv(pro_dir+'generole.tsv', sep='\t', low_memory=False, header=0)

            entry_list = [list(a) for a in zip(filei.entry1, filei.entry2)]

            # Create graph
            G = nx.Graph()
            # Fill graph with data
            G.add_edges_from(entry_list)

            all_paths = []
            #obtain all combinations of paths of 4 components from node to node
            for start in G.nodes:
                for end in G.nodes:
                    if start != end:
                        all_paths = all_paths + list(nx.all_simple_paths(G, start, end, 3))

            #remove duplicates
            import itertools
            all_paths.sort()
            all_paths_dedup = [list(a) for a in (all_paths for all_paths,_ in itertools.groupby(all_paths))]

            #discard elements after a non-valid association, discarded next valid elements will appear in the begining of another path retrieved by all_simple_paths
            #this allows also to remove the nodes linked in backward direction from all_simple_paths
            generole = pd.read_csv(pro_dir+'generole.tsv', sep='\t', low_memory=False, header=0)
            all_paths_checked = []
            import more_itertools
            for path in all_paths_dedup:
                path_it = []
                for pair in list(more_itertools.windowed(path,n=2, step=1)):
                    associations = [item for sublist in [x.split(',') for x in filei.loc[(filei['entry1'] == pair[0]) & (filei['entry2'] == pair[1])]['name'].tolist()] for item in sublist]

                    if len(associations) > 0: #valid: association is in relations file
                        if len(set(associations) & set(['indirect effect'])) == 0: #valid: asociation is not indirect
                            if str(dictionary.loc[dictionary['node'] == pair[0]]['symbol'].tolist()[0]) != 'nan' and str(dictionary.loc[dictionary['node'] == pair[1]]['symbol'].tolist()[0]) != 'nan': #valid: node has gene symbol attached
                                #update source node
                                genes = []
                                for gene in dictionary.loc[dictionary['node'] == pair[0]]['symbol'].tolist()[0].split(':'):
                                    #decide the role of the gene
                                    roles = ':'.join(list(set([x for x in generole.loc[generole['gene'] == gene].values.flatten().tolist()[1:] if str(x) != 'nan'])))
                                    #check if there is a restrictive association for the pair, in that case the pair is discarded
                                    if (roles in ['ONC', ''] and len(set(associations) & set(['activation', 'expression'])) > 0) or (roles == 'TSG' and len(set(associations) & set(['inhibition', 'repression'])) > 0) or (roles not in ['ONC', '', 'TSG']) or (len(set(['activation', 'expression', 'inhibition', 'repression']) & set(associations)) == 0):
                                        genes.append(gene)

                                gene_names = ':'.join(genes)
                                if gene_names != '': #association only saved if genes in node
                                    if len(path_it) == 0: path_it = [gene_names]
                                    else: path_it[-1] = gene_names

                                    #report destination node
                                    genes = []
                                    for gene in dictionary.loc[dictionary['node'] == pair[1]]['symbol'].tolist()[0].split(':'):
                                       genes.append(gene)
                                    gene_names = ':'.join(genes)
                                    path_it.append(gene_names)

                            else:
                                break
                        else:
                            break
                    else:
                        break

                all_paths_checked.append(path_it)

            #remove duplicates
            all_paths_checked = [x for x in all_paths_checked if len(x) > 1]
            all_paths_checked.sort()
            all_paths_checked_dedup = [list(a) for a in (all_paths_checked for all_paths_checked,_ in itertools.groupby(all_paths_checked))]

            #remove pathways included in other ones
            all_paths_unified = []
            for path in all_paths_checked_dedup:
                found = False
                for path_search in all_paths_checked_dedup:
                    #if the length of both the path being searched and the path for the search is the same in the match, the search is on itself because redundances have been removed
                    if any(path == path_search[i:i+len(path)] for i in range(len(path_search))) and len(path) != len(path_search):
                        found = True
                        break
                if not found:        
                    all_paths_unified.append(path)

            for path in all_paths_unified:
                fileo.write('\t'.join([str(x) for x in path])+'\n')
            fileo.close()

    def obtain_upstream_genes():
        upstream_genes = {}
        for i in progressbar.progressbar(range(len(paths))):
            filei = open(pro_dir_KEGG+'subpaths/'+paths[i]+'.symbol.tsv', 'r')
            for line in filei:
                line = line.rstrip('\n')
                line_a = line.split('\t')
                for idx, col in enumerate(line_a[1:]):
                    for gene in col.split(':'):
                        if gene == 'MAPK1' and 'CDK5' in [x for ele in line_a[:idx+1] for x in ele.split(':')]:
                            print(paths[i])
                            print(line_a)
                        if gene in upstream_genes.keys():
                            upstream_genes[gene] = upstream_genes[gene] + [x for ele in line_a[:idx+1] for x in ele.split(':')]
                        else:
                            upstream_genes[gene] = [x for ele in line_a[:idx+1] for x in ele.split(':')]
        filei.close()
        fileo = open(pro_dir_KEGG+'upstream_genes.tsv', 'w')
        fileo.write('\t'.join(['gene', 'upstream_genes'])+'\n')
        for gene in upstream_genes:
            if gene != '':
                upstream_genes_list = list(set(upstream_genes[gene]))
                upstream_genes_list = list(filter(None, upstream_genes_list))
                if len(upstream_genes_list) > 0:
                    fileo.write('\t'.join([gene, '|'.join(upstream_genes_list)])+'\n')
        fileo.close()

    parse_kegg_xml()
    nodes_ids_and_symbols()
    obtain_subroutes()
    obtain_upstream_genes()

def process_gene_names():

    print('Creating file with checked gene symbols...')

    genes_chk = 'genes_checked.tsv'

    genes = []

    files = [dgidb, sabdab_file, moalmanac_file, GDSC]

    for f in files:
        inputf=open(pro_dir+f,'r')
        for line in inputf:
            line = line.rstrip('\n')
            line_a = line.split('\t')
            if 'gene_name' in line_a: gene_index = line_a.index('gene_name')
            if f == GDSC: genes = genes+line_a[gene_index].split(',')
            elif f == sabdab_file: genes = genes+line_a[gene_index].split(';')
            else: genes.append(line_a[gene_index])
        inputf.close()

    outputf = open(pro_dir+genes_chk, 'w')
    outputf.write('\t'.join(['gene_name', 'checked_gene_symbol'])+'\n')

    import httplib2 as http
    import json

    try:
        from urlparse import urlparse
    except ImportError:
        from urllib.parse import urlparse

    headers = {
        'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    method = 'GET'
    body = ''
    h = http.Http()

    genes = list(set(genes))

    for i in progressbar.progressbar(range(len(genes))):
        found = False
        symbol = ''

        path = '/search/symbol/'+genes[i]

        target = urlparse(uri+path)

        response, content = h.request(target.geturl(), method, body, headers)

        if response['status'] == '200':
        # assume that content is a json reply
        # parse content with the json module 
            data = json.loads(content)
            if len(data['response']['docs']) > 0:
                outputf.write('\t'.join([genes[i], data['response']['docs'][0]['symbol']])+'\n')
                found = True
        else:
            print ('Error detected: ' + response['status'])

        if not found:
            path = '/search/alias_symbol/'+genes[i]
            target = urlparse(uri+path)

            response, content = h.request(target.geturl(), method, body, headers)

            if response['status'] == '200':
                data = json.loads(content)
                if len(data['response']['docs']) > 0:
                    outputf.write('\t'.join([genes[i], data['response']['docs'][0]['symbol']])+'\n')
                    found = True
            else:
                print ('Error detected: ' + response['status'])

        if not found:
            path = '/search/prev_symbol/'+genes[i]
            target = urlparse(uri+path)

            response, content = h.request(target.geturl(), method, body, headers)

            if response['status'] == '200':
                data = json.loads(content)
                if len(data['response']['docs']) > 0:
                    outputf.write('\t'.join([genes[i], data['response']['docs'][0]['symbol']])+'\n')
                    found = True
            else:
                print ('Error detected: ' + response['status'])

        if not found:
            print('Not found:'+genes[i])
            outputf.write('\t'.join([genes[i], genes[i]])+'\n')

    outputf.close()

def process_civic():

    print('Processing civic file...')

    civic_file = 'civic.tsv'
    civic_dwn = 'civic_evidence.tsv'

    output = open(pro_dir+civic_file, 'w')
    output.write('\t'.join(['drug_name', 'gene_name', 'variation', 'response', 'source'])+'\n')
    civic = pd.read_csv(dwn_dir+civic_dwn, sep ='\t', low_memory=False)

    civic_select = civic.loc[(civic['status'] == 'ACCEPTED') & (civic['evidence_direction'] == 'SUPPORTS') & (civic['clinical_significance'].isin(['RESISTANCE', 'SENSITIVITYRESPONSE'])) & (civic['evidenceLevel'].isin(['A']))] # A: Validated association

    drug_gene = civic_select.drop_duplicates(subset=['drug', 'gene'])

    for index, row in drug_gene.iterrows():
        alteration_sen = []
        alteration_res = []
        civic_subset = civic_select.loc[(civic_select['drug'] == row['drug']) & (civic_select['gene'] == row['gene'])]

        for index2, row2 in civic_subset.iterrows():
            if row2['clinical_significance'] == 'SENSITIVITYRESPONSE':
                alteration_sen.append(row2['variant'])
            else:
                alteration_res.append(row2['variant'])

        if len(alteration_sen) > 0 and len(alteration_res) > 0:
            response = 'sensitivity / resistance'
            variation = '; '.join(list(set(alteration_sen)))+' / '+'; '.join(list(set(alteration_res)))
        elif len(alteration_sen) > 0:
            response = 'sensitivity'
            variation = '; '.join(list(set(alteration_sen)))
        else:
            response = 'resistance'
            variation = '; '.join(list(set(alteration_res)))


        output.write('\t'.join([row['drug'], row['gene'], variation, response, 'CIViC'])+'\n')
    output.close()

def process_oncoKB():

    print('Processing oncoKB file...')

    oncokb_file = 'oncokb.tsv'
    oncokb_dwn = 'oncokb_biomarker_drug_associations.tsv'

    output = open(pro_dir+oncokb_file, 'w')
    output.write('\t'.join(['drug_name', 'gene_name', 'variation', 'response', 'source'])+'\n')
    oncokb = pd.read_csv(dwn_dir+oncokb_dwn, sep ='\t', low_memory=False)

    oncokb_select = oncokb.loc[(oncokb['Level'].isin(['1', '2', 'R1'])) & (oncokb['Gene'] != 'Other Biomarkers')] # 1, 2 and R1 (top evidence)
    drug_gene = oncokb_select.drop_duplicates(subset=['Drugs (for therapeutic implications only)', 'Gene'])

    for index, row in drug_gene.iterrows():
        alteration_sen = []
        alteration_res = []
        oncokb_subset = oncokb_select.loc[(oncokb_select['Drugs (for therapeutic implications only)'] == row['Drugs (for therapeutic implications only)']) & (oncokb_select['Gene'] == row['Gene'])]

        for index2, row2 in oncokb_subset.iterrows():
            if row2['Level'] in ['1', '2', '3A']:
                alteration_sen.append(row2['Alterations'])
            else:
                alteration_res.append(row2['Alterations'])

        if len(alteration_sen) > 0 and len(alteration_res) > 0:
            response = 'sensitivity / resistance'
            variation = '; '.join(list(set(alteration_sen)))+' / '+'; '.join(list(set(alteration_res)))
        elif len(alteration_sen) > 0:
            response = 'sensitivity'
            variation = '; '.join(list(set(alteration_sen)))
        else:
            response = 'resistance'
            variation = '; '.join(list(set(alteration_res)))

        for d in row['Drugs (for therapeutic implications only)'].split(', '):
            if re.search('\+', d):
                d = ' + '.join(sorted(d.split(' + ')))
            output.write('\t'.join([d, row['Gene'], variation, response, 'OncoKB'])+'\n')
    output.close()

def process_intogen():

    print('Processing intogen file...')

    intogen = 'intogen.tsv'
    intogen_dwn = 'download'

    output = open(pro_dir+intogen, 'w')
    output.write('\t'.join(['gene_name', 'cohort', 'cancer_type', 'qvalue_combination'])+'\n')

    with zipfile.ZipFile(dwn_dir+intogen_dwn) as z:
        for i in z.infolist():
            if re.search('Compendium_Cancer_Genes.tsv', i.filename):
                with z.open(i.filename) as f:
                    data = pd.read_csv(f, sep ='\t', header=0, low_memory=False)
                    for index, row in data.iterrows():
                        output.write('\t'.join([row['SYMBOL'], row['COHORT'], row['CANCER_TYPE'], str(row['QVALUE_COMBINATION'])])+'\n')

    output.close()

def process_depmap():

    print('Processing DepMap public score + Chronos...')

    depmap_dwn = 'CRISPR_gene_effect.csv'
    depmap_pro = 'chronos_skew.tsv'

    matrix = pd.read_csv(dwn_dir+depmap_dwn, low_memory=False)
    matrix = matrix.set_index('DepMap_ID')
    matrix.columns = [x.split(' ')[0] for x in matrix.columns.tolist()]
    skewness = matrix.apply(lambda x : scipy.stats.skew(x, nan_policy='omit'))

    skewness.to_csv(pro_dir+depmap_pro, index=True, sep='\t', header=False)

def process_KEGG_pathways():

    print('Processing KEGG gene pathway file...')

    genepathway_file = 'gene_pathway.tsv'

    inputf = pd.read_csv(dwn_dir+genepathway_file, sep ="\t", header=None, low_memory=False)

    inputf.columns =['KEGG Gene ID', 'KEGG Pathway ID']
    inputf['KEGG Gene ID'] = inputf['KEGG Gene ID'].str.replace('hsa:','')
    inputf['KEGG Pathway ID'] = inputf['KEGG Pathway ID'].str.replace('path:','')

    array_agg = lambda x: '|'.join(x.astype(str))
    inputf = inputf.groupby('KEGG Gene ID').agg({'KEGG Pathway ID': array_agg})

    inputf.to_csv(pro_dir+genepathway_file, sep = "\t")

    print('Processing KEGG pathway descriptions file...')

    pathwaydesc_file = 'pathway_desc.tsv'

    inputf = pd.read_csv(dwn_dir+pathwaydesc_file, sep ="\t", header=None, low_memory=False)

    inputf.columns =['KEGG Pathway ID', 'KEGG Pathway desc']
    inputf['KEGG Pathway ID'] = inputf['KEGG Pathway ID'].str.replace('path:','')
    inputf['KEGG Pathway desc'] = inputf['KEGG Pathway desc'].str.replace(' \- Homo sapiens \(human\)','')

    inputf.to_csv(pro_dir+pathwaydesc_file, sep = "\t", index=False)

def process_SL():

    print('Processing SL dependencies file...')

    sl_file = 'all_genetic_dependencies.tsv'

    inputf = pd.read_csv(dwn_dir+sl_file, sep ="\t", low_memory=False)
    outputf = open(pro_dir+'genetic_dependencies.tsv','w')
    outputf.write('\t'.join(['gene','genetic_dependency'])+'\n')

    sel = inputf.loc[inputf['padj'] < 0.05].groupby("dependency")

    for dependency, frame in sel:
        dep_genes = []
        for index, row in sel.get_group(dependency).iterrows():
            if not row['gene'] == dependency:
                dep_genes.append(row['gene']+'('+row['alteration']+')')
        if len(dep_genes) > 0:
            outputf.write('\t'.join([dependency,'|'.join(dep_genes)])+'\n')

    outputf.close()

def create_cosmic_temp_files(idx):
    output_file = pro_dir+'temp/'+'COSMIC_'+str(idx)+'.tsv'
    select_column = ['Gene name', 'ID_sample', 'GENOMIC_MUTATION_ID', 'GRCh', 'Mutation genome position', 'FATHMM prediction', 'Mutation somatic status', 'HGVSC']
    outputf = open(output_file,'w')
    outputf.write('\t'.join(select_column)+'\n')
    outputf.close()

    cleaned_genes_group = [x for x in genes_group[idx] if str(x) != 'None']
    gene_list = ['^'+x+'_' for x in cleaned_genes_group]

    command = 'zcat '+pro_dir+outputn+' | egrep \''+'|'.join(gene_list)+'\' >> '+ output_file
    subprocess.call(command,shell=True)

def process_file(filename):
    print(filename)
    cosmic = {}
    gene_freq = {}
    mut_freq = {}

    inputf = pd.read_csv(pro_dir+'temp/'+filename, sep ="\t", header=0, low_memory=False)

    if not inputf.empty:
        inputf[['Gene name']] = inputf['Gene name'].str.split('_',expand=True)[[0]]
        inputf[['HGVSC_Transcript','HGVSC']] = inputf['HGVSC'].str.split(':',expand=True)
        inputf[['HGVSC_Transcript']] = inputf['HGVSC_Transcript'].str.split('.',expand=True)[[0]]
        outputf = open(pro_dir+'temp/'+filename.replace('.tsv', '_prc.tsv'), 'w')

        gene_list = inputf['Gene name'].unique().tolist()
        for gene in gene_list:
            samples_gene = len(inputf.loc[inputf['Gene name'] == gene]['ID_sample'].unique().tolist())
            mut_list = inputf.loc[inputf['Gene name'] == gene]['GENOMIC_MUTATION_ID'].unique().tolist()
            for mut in mut_list:
                samples_mut = len(inputf.loc[inputf['GENOMIC_MUTATION_ID'] == mut]['ID_sample'].unique().tolist())
                transcript_list = inputf.loc[(inputf['Gene name'] == gene) & (inputf['GENOMIC_MUTATION_ID'] == mut)]['HGVSC_Transcript'].unique().tolist()
                for trans in transcript_list:
                    FATHMM = ['' if x is np.nan else x for x in inputf.loc[(inputf['Gene name'] == gene) & (inputf['GENOMIC_MUTATION_ID'] == mut) & (inputf['HGVSC_Transcript'] == trans)]['FATHMM prediction'].tolist()][0]
                    HGVSc = inputf.loc[(inputf['Gene name'] == gene) & (inputf['GENOMIC_MUTATION_ID'] == mut) & (inputf['HGVSC_Transcript'] == trans)]['HGVSC'].tolist()[0]

                    outputf.write('\t'.join([':'.join([gene, trans, HGVSc]), mut, FATHMM, str(samples_gene), str(samples_mut), str(total_cosmic_rec)])+'\n')
        outputf.close()

def create_cosmic_temp_files(idx):
    outputn = cosmic_file.replace('.tsv.gz','.tsv.filtered.gz')
    output_file = pro_dir+'temp/'+'COSMIC_'+str(idx)+'.tsv'
    select_column = ['Gene name', 'ID_sample', 'GENOMIC_MUTATION_ID', 'GRCh', 'Mutation genome position', 'FATHMM prediction', 'Mutation somatic status', 'HGVSC']
    outputf = open(output_file,'w')
    outputf.write('\t'.join(select_column)+'\n')
    outputf.close()

    cleaned_genes_group = [x for x in genes_group[idx] if str(x) != 'None']
    gene_list = ['^'+x+'_' for x in cleaned_genes_group]

    command = 'zcat '+pro_dir+outputn+' | egrep \''+'|'.join(gene_list)+'\' >> '+ output_file
    subprocess.call(command,shell=True)

def process_cosmic():

    global cosmic_file
    cosmic_file = 'CosmicMutantExport.tsv.gz'

    print('Processing COSMIC data...')

    print('-->Filtering COSMIC...')
    #filtering columns
    inputf = gzip.open(dwn_dir+cosmic_file,'r')
    first_line = inputf.readline().rstrip().decode().split('\t')
    global outputn
    outputn = cosmic_file.replace('.tsv.gz','.tsv.filtered.gz')
    select_column = ['Gene name', 'ID_sample', 'GENOMIC_MUTATION_ID', 'GRCh', 'Mutation genome position', 'FATHMM prediction', 'Mutation somatic status', 'HGVSC']
    inputf.close()

    #filtering records
    idx = [first_line.index(i) for i in select_column]
    cols = [str(i + 1) for i in idx]
     #algunos registros tienen el identificador vacio y no tienen por hgvsc (no los cuento)
    command = 'zcat '+dwn_dir+cosmic_file+' | cut -f '+','.join(cols)+' | awk -F"\t" \'($3 != "" && $4 == "38") || $1 == "Gene name" {print}\' | gzip > '+pro_dir+outputn
    subprocess.call(command,shell=True)

    #obtaining list of genes
    gfile = pro_dir+'gene_list.txt'
    command = 'zcat '+pro_dir+outputn+' | cut -f 1 | cut -d\'_\' -f 1 | grep -v Gene | sort | uniq > '+gfile
    subprocess.call(command,shell=True)

    print('-->Creating COSMIC file...')
    genes_file = pd.read_csv(gfile, sep ="\t", header=None, low_memory=False)
    gene_list = genes_file[0].tolist()

    def grouper(n, iterable, fillvalue=None):
        import itertools
        args = [iter(iterable)] * n
        return itertools.zip_longest(*args, fillvalue=fillvalue)

    global genes_group
    genes_group = list(grouper(10, gene_list))

    os.mkdir(pro_dir+'temp/')
    a_pool = multiprocessing.Pool(processes=9)

    result_list_tqdm = []
    for result in tqdm(a_pool.imap(create_cosmic_temp_files, range(0,len(genes_group))), total=len(genes_group)):
        result_list_tqdm.append(result)

    with gzip.open(pro_dir+outputn,'r') as f:
        global total_cosmic_rec
        total_cosmic_rec = len(f.readlines()) - 1 #removing the header for the count

    a_pool = multiprocessing.Pool(processes=9)

    result_list_tqdm = []
    for result in tqdm(a_pool.imap(process_file, os.listdir(pro_dir+'temp/')), total=len(os.listdir(pro_dir+'temp/'))):
        result_list_tqdm.append(result)

    read_files = glob.glob(pro_dir+'temp/'+"*_prc.tsv")
    pdb.set_trace()
    with open(pro_dir+'COSMIC.tsv', 'ab') as outfile:
        outfile.write(bytes('\t'.join(['Cosmic_key', 'cosmic_id', 'FATHMM', 'Gene_freq', 'Mut_freq', 'Total'])+'\n','utf-8'))
        for f in read_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())

    shutil. rmtree(pro_dir+'temp/')

def process_xml_clinvar(ofile):

    outputf = open(ofile,'a')

    tree = ET.parse(xml_section)
    root = tree.getroot()

    #print (root.tag)
    #print (root.attrib)
    #for child in root:
    #    print(child.tag, child.attrib)
    for ClinVarSet in root.iter('ClinVarSet'):
        (ACC, ASSEMBLY, CHR, VCF_POS, VCF_REF, VCF_ALT, SYMBOL, TRAIT, SIG) = ([], [], [], [], [], [], [], [], [])
        for ReferenceClinVarAssertion in ClinVarSet.iter('ReferenceClinVarAssertion'):
            for ClinVarAccession in ReferenceClinVarAssertion.iter('ClinVarAccession'):
                ACC.append(ClinVarAccession.attrib['Acc'])
            for MeasureSet in ReferenceClinVarAssertion.iter('MeasureSet'):
                for Measure in MeasureSet.iter('Measure'):
                    for SequenceLocation in Measure.iter('SequenceLocation'):
                        if 'positionVCF' in SequenceLocation.attrib.keys():
                            ASSEMBLY.append(SequenceLocation.attrib['Assembly'])
                            CHR.append(SequenceLocation.attrib['Chr'])
                            VCF_POS.append(SequenceLocation.attrib['positionVCF'])
                            VCF_REF.append(SequenceLocation.attrib['referenceAlleleVCF'])
                            VCF_ALT.append(SequenceLocation.attrib['alternateAlleleVCF'])
                    for Symbol in Measure.iter('Symbol'):
                        for ElementValue in Symbol.iter('ElementValue'):
                            if ElementValue.attrib['Type'] == 'Preferred': SYMBOL.append(ElementValue.text)
            for TraitSet in ReferenceClinVarAssertion.iter('TraitSet'):
                for Trait in TraitSet:
                    for Name in Trait.iter('Name'):
                        for ElementValue in Name.iter('ElementValue'):
                            if ElementValue.attrib['Type'] == 'Preferred': TRAIT.append(ElementValue.text)
            for ClinicalSignificance in ReferenceClinVarAssertion.iter('ClinicalSignificance'):
                for Description in ClinicalSignificance.iter('Description'):
                    SIG.append(Description.text)

        for idx in range(len(ASSEMBLY)):
            outputf.write('\t'.join(['::'.join(list(set(ACC))),ASSEMBLY[idx],CHR[idx],VCF_POS[idx],VCF_REF[idx],VCF_ALT[idx],'::'.join(list(set(SYMBOL))),'::'.join(list(set(TRAIT))),'::'.join(list(set(SIG)))])+'\n')

    outputf.close()

def process_clinvar():

    print('Processing ClinVar data...')

    clinvar_file = 'ClinVarFullRelease_00-latest.xml.gz'

    print('-->Filtering ClinVar...')

    #filter records
    inputf = gzip.open(dwn_dir+clinvar_file,'rb')
    outputn = clinvar_file.replace('.xml.gz','tagfiltered.xml.gz')
    outputf = gzip.open(pro_dir+outputn,'wb')

    for l in inputf:
        if re.search('xml|<ReleaseSet|</ReleaseSet|<ClinVarSet|</ClinVarSet|<ReferenceClinVarAssertion|</ReferenceClinVarAssertion|<ClinVarAccession|</ClinVarAccession|<MeasureSet|</MeasureSet|<Measure|</Measure|<SequenceLocation|</SequenceLocation|<Symbol|</Symbol|<ElementValue|</ElementValue|<TraitSet|</TraitSet|<Trait|</Trait|<Name|</Name|<ClinicalSignificance|</ClinicalSignificance|<Description|</Description',l.decode()):
            outputf.write(l)

    outputf.close()

    print('-->Creating ClinVar file...')
    #create tsv file
    global xml_section
    xml_section = pro_dir+'xml_temp.xml'
    outputf = open(pro_dir+'Clinvar.tsv','w')
    outputf.write('\t'.join(['Acc','Assembly','Chr','VCF_pos','VCF_ref','VCF_alt','Gene','Trait','Significance'])+'\n')
    outputf.close()

    #counter for ClinVarSet tag
    counter = 0
    xml = gzip.open(pro_dir+outputn, 'rb')

    for line in xml:
        line = line.decode()
        if re.search('^<\?xml version|^<ReleaseSet|^</ReleaseSet',line):
            continue
        elif re.search('^<ClinVarSet',line) and counter == 0:
            xml_section_file = open(xml_section,'w')
            xml_section_file.write('<ReleaseSet Dated="2021-10-30" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_1.64.xsd">'+'\n')
            xml_section_file.write(line)
        elif re.search('^</ClinVarSet',line):
            #limit of tags per processing
            if counter < 10000:
                xml_section_file.write(line)
                counter += 1
            else:
                xml_section_file.write(line)
                xml_section_file.write('</ReleaseSet>')
                xml_section_file.close()
                counter = 0
                process_xml_clinvar(pro_dir+'Clinvar.tsv')
        else:
            xml_section_file.write(line)

    xml.close()
    os.remove(xml_section)
    os.remove(pro_dir+outputn)

def process_pfam():

    print('Processing Pfam data...')

    pfam_file = 'Pfam-A.full.gz'

    print('-->Filtering Pfam...')
    #filter records
    inputf = gzip.open(dwn_dir+pfam_file,'rb')
    outputn = pfam_file.replace('.full.gz','.full.filtered.gz')
    outputf = gzip.open(pro_dir+outputn,'wb')

    for l in inputf:
        if re.search('^#=GF ID|^#=GF AC|^#=GF DE|_HUMAN|^//',l.decode('latin-1')):
            outputf.write(l)

    outputf.close()

    print('-->Creating Pfam file...')
    #create tsv file
    outputf = open(pro_dir+'Pfam-A.full.tsv','w')
    outputf.write('\t'.join(['DOMAIN_ID','PFAM_ACC','DOMAIN_DESCRIPT','PROTEIN_NAME','PROTEIN_ACC','START','END'])+'\n')

    filei = gzip.open(pro_dir+outputn,'rb')

    (domain_id, pfam_acc, domain_des, protein_name, protein_acc, start, end) = ('', '', '', '', '', '', '')

    for line in filei:
        line = line.decode()
        line = line.strip('\n')
        line_a = line.split(' ')
        if line_a[0] == '//':
            (domain_id, pfam_acc, domain_des, protein_name, protein_acc, start, end) = ('', '', '', '', '', '', '')
        else:
            if line_a[1] == 'ID': domain_id = ' '.join(line_a[4:])
            if line_a[1] == 'AC': pfam_acc = ' '.join(line_a[4:])
            if line_a[1] == 'DE': domain_des = ' '.join(line_a[4:])
            if line_a[0] == '#=GS' and len(line_a) > 10:
                if 'AC' in line_a:
                    protein_name = line_a[1].split('_')[0]
                    protein_acc = line_a[line_a.index('AC')+1].split('.')[0]
                    start = line_a[1].split('/')[1].split('-')[0]
                    end = line_a[1].split('/')[1].split('-')[1]

                    outputf.write('\t'.join([domain_id, pfam_acc, domain_des, protein_name, protein_acc, start, end])+'\n')

    filei.close()
    outputf.close()
    os.remove(pro_dir+outputn)

def process_interpro():

    print('Processing Interpro data...')

    interpro_file = 'match_complete.xml.gz'

    #filter records
    print('-->Filtering Interpro...')
    inputf = gzip.open(dwn_dir+interpro_file,'rb')
    outputn = interpro_file.replace('.xml.gz','.xml.filtered.gz')
    outputf = gzip.open(pro_dir+outputn,'wb')

    for line in inputf:
        line = line.decode().rstrip()
        if re.search('^<?xml', line) or re.search('interpromatch', line):
            outputf.write(bytes(line+'\n','utf-8'))
        elif re.search('<protein', line) and re.search('_HUMAN', line):
            while line != '</protein>':
                outputf.write(bytes(line+'\n','utf-8'))
                line = next(inputf).decode().rstrip()
            outputf.write(bytes(line+'\n','utf-8'))

    inputf.close()
    outputf.close()

    print('-->Creating Interpro file...')

    outputf = open(pro_dir+'Interpro.tsv','w')
    outputf.write('\t'.join(['DOMAIN_ID','DOMAIN_DESCRIPT','PROTEIN_NAME','PROTEIN_ACC','START','END'])+'\n')

    inputf = gzip.open(pro_dir+outputn,'rb')
    tree = ET.parse(inputf)
    root = tree.getroot()

    for protein in root.iter('protein'):
        (domain_id, pfam_acc, domain_des, protein_name, protein_acc, start, end) = ('', '', '', '', '', '', '')
        protein_acc = protein.attrib['id']
        protein_name = protein.attrib['name'].split('_')[0]

        for match in protein.iter('match'):
            for ipr in match.iter('ipr'):
                domain_id = ipr.attrib['id']
                domain_des = ipr.attrib['name']
                for lcn in match.iter('lcn'):
                    start = lcn.attrib['start']
                    end = lcn.attrib['end']
                    outputf.write('\t'.join([domain_id, domain_des, protein_name, protein_acc, start, end])+'\n')

    outputf.close()
    inputf.close()
    os.remove(pro_dir+outputn)

def process_uniprot():

    print('Processing Uniprot data...')

    uniprot_file = 'uniprot_sprot.xml.gz'

    print('-->Filtering Uniprot...')
    #filter records
    inputf = gzip.open(dwn_dir+uniprot_file,'rb')
    outputn = uniprot_file.replace('.xml.gz','.xml.filtered.gz')
    outputf = gzip.open(pro_dir+outputn,'wb')

    for l in inputf:
        l = l.decode()
        if re.search('\?xml|<uniprot|</uniprot|xmlns:|xsi:|<entry|</entry|<accession>|</accession>|<gene>|</gene>|<name .*type=\"primary\"|<organism>|<organism |</organism>|<name type=\"common\"',l):
            if re.search('<uniprot',l): l = '<uniprot>\n'
            if not re.search('^ xmlns:|^ xsi:',l):
                l = re.sub(r' xmlns\=\"http\:\/\/uniprot.org\/uniprot\"','',l)
                outputf.write(bytes(l,'utf-8'))

    outputf.close()

    print('-->Creating Uniprot file...')

    outputf = open(pro_dir+'Uniprot.tsv','w')
    outputf.write('\t'.join(['PROTEIN_ID','GENE_NAME','PUBMED_ID','dbSNP_ref'])+'\n')

    inputf = gzip.open(pro_dir+outputn,'rb')
    tree = ET.parse(inputf)
    root = tree.getroot()

#    for child in root:
#        print(child.tag, child.attrib)

    for entry in root.iter('entry'):
        (protein_id, gene_name, organism) = ([], '', '')
        for accession in entry.iter('accession'):
            protein_id.append(accession.text)
        for gene in entry.iter('gene'):
            for name in gene.iter('name'):
                gene_name = name.text
        for organism in entry.iter('organism'):
            for name in organism.iter('name'):
                organism = name.text
        if organism == 'Human':
            outputf.write('\t'.join([';'.join(protein_id),gene_name,'',''])+'\n')

    outputf.close()
    os.remove(pro_dir+outputn)

#process_DGIdb()
#process_sabdab()
#process_moalmanac()
#process_GDSC()
#process_KEGG_ATC()
#process_cmap()
#process_FDA()
#process_FDA_label()
#process_EMA()
#process_ct()
#process_cgc_oncovar()
#process_KEGG_ind()
#process_gene_names()
#corregir a mano NKX2.1 que se interpreta como NKX2 y 1 por separado
#process_civic()
#process_oncoKB()
#process_intogen()
#process_depmap()
#process_KEGG_pathways()
process_SL()

#exclusive for genomic annotation
#process_cosmic()
#process_clinvar()
#process_pfam()
#process_interpro()
#process_uniprot()
