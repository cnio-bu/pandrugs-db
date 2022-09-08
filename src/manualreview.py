#!/usr/bin/python

import re
import pdb
import pandas as pd
import progressbar

#Set to directory with processed files, pandrugs previous version and manual review files
pro_dir = '/home/epineiro/Analysis/PanDrugs/2.0/processed/'
mr_dir = '/home/epineiro/Analysis/PanDrugs/2.0/manual_review/'
pd_prev = '/home/epineiro/Analysis/PanDrugs/Update_Feb2020/Pandrugs_Feb2020.tsv'

def process_drug_names():

    print('Creating file with synonyms for review...')

    drugs = []

    files = ['DGIdb_interactions.tsv', 'oncokb.tsv', 'civic.tsv', 'moalmanac.tsv', 'GDSC.tsv']

    for f in files:
        inputf=open(pro_dir+f,'r')
        for line in inputf:
            line = line.rstrip('\n')
            line_a = line.split('\t')
            if 'drug_name' in line_a: drug_index = line_a.index('drug_name')
            drugs.append(line_a[drug_index])
        inputf.close()

    outputf = open(mr_dir+synonyms_mr, 'w')
    outputf.write('\t'.join(['drug_name', 'short_drug_name', 'standard_drug_name', 'short_standard_drug_name', 'show_drug_name', 'short_show_drug_name', 'synonyms', 'review'])+'\n')

    import pubchempy as pcp

    drugs = list(set(drugs))
    for i in progressbar.progressbar(range(len(drugs))):
        response = pcp.get_synonyms(drugs[i], 'name')
        INNname = ''
        if len(response) > 0:
            for e in response[0]['Synonym']:
                if '[' in e and 'INN' in e and 'INN-' not in e:
                    INNname = e.split(' [')[0]
                    break
            if INNname == '': INNname = response[0]['Synonym'][0]
            outputf.write('\t'.join([drugs[i], drugs[i].split(' ')[0], response[0]['Synonym'][0], response[0]['Synonym'][0].split(' ')[0], INNname, INNname.split(' ')[0], '::'.join(response[0]['Synonym']), ''])+'\n')
        else:
            outputf.write('\t'.join([drugs[i], drugs[i].split(' ')[0], drugs[i], drugs[i].split(' ')[0], drugs[i], drugs[i].split(' ')[0], '', ''])+'\n')

    outputf.close()

    drug_file = pd.read_csv(mr_dir+synonyms_mr, sep ='\t', low_memory=False)
    drug_file = drug_file.fillna('')

    counts = 0
    for index, row in drug_file.iterrows():
        gene = row['drug_name']
        sdn = row['standard_drug_name']
        sdns = []
        counts = counts + 1
        print(counts)
        for index2, row2 in drug_file.iterrows():
            review = ''
            sdn2 = row2['standard_drug_name']
            synonyms = row2['synonyms'].split('::')
            synonyms = [x.upper() for x in synonyms]
            if gene.upper() in synonyms and sdn != sdn2:
                sdns.append(sdn2)
            if review != '':
                pdb.set_trace()
                drug_file[index,'review'] = 'review'
                drug_file[index2,'review'] = 'review'
    
    drug_file.to_csv(mr_dir+synonyms_mr, index=False, sep='\t', header=True)

def process_FDA_EMA():

    print('Collecting FDA and EMA information about status...')

    fda_status_file = 'fda_status.tsv'
    approved_drug_data_mr = 'drug_approved_data_mr.tsv'
    ema_status_dwn = 'ema_status.tsv'
    fda_labels_dwn = 'fda_labels.tsv'

    fda_status = pd.read_csv(pro_dir+fda_status_file, sep ='\t', low_memory=False)
    fda_status.fillna('', inplace=True)

    ema_status = pd.read_csv(pro_dir+ema_status_dwn, sep ='\t', low_memory=False)
    ema_status.fillna('', inplace=True)
    ema_status['International non-proprietary name (INN) / common name'] = ema_status['International non-proprietary name (INN) / common name'].str.upper()

    fda_labels = pd.read_csv(pro_dir+fda_labels_dwn, sep ='\t', low_memory=False)
    fda_labels.fillna('', inplace=True)
    appl_list = fda_labels['ApplNo'].tolist()

    status_CT_file = pd.read_csv(pro_dir+'clinicaltrials.tsv', sep='\t', low_memory=False)
    status_CT_file['drug'] = status_CT_file['drug'].str.upper()
    status_CT_file[['drug']] = status_CT_file[['drug']].fillna('')

    fileo = open(mr_dir+approved_drug_data_mr, 'w')
    fileo.write('\t'.join(['standard_drug_name','ApplNo','ApplNo_labels','status_FDA','status_EMA','ct_condition','status_mc','PrevStat','indication_FDA','indication_EMA','indication_mc','indication_prev','pathology_mc','pathology_prev','cancer_type_mc','cancer_type_prev','drug_type_mc','drug_type_prev','short_name'])+'\n')

    drug_list = list(set(synonyms['standard_drug_name'].tolist()))
    for i in progressbar.progressbar(range(len(drug_list))):
        standard_drug_name = drug_list[i]
        #list of possible names in pandrugs database
        synonyms_list = synonyms.loc[synonyms['standard_drug_name'] == standard_drug_name][['drug_name', 'show_drug_name', 'standard_drug_name', 'short_drug_name', 'short_standard_drug_name', 'short_show_drug_name']].values.flatten().tolist()

        (fda_data,status_ema) = (['', ''],'')
        review = []

        #Data from FDA
        found = False
        for idx, syn in enumerate(synonyms_list):
            if syn in fda_status['drug'].tolist():
                fda_data = [str(x) for x in fda_status.loc[fda_status['drug'] == syn].values.flatten().tolist()[1:]]
                if idx > 2: review.append('fda')
                found = True
                break
        if not found:
            for idx, syn in enumerate(synonyms_list):
                if syn in fda_status['drug'].str.split(' ', n = 1, expand = True)[0].tolist():
                    index = fda_status['drug'].str.split(' ', n = 1, expand = True)[0].tolist().index(syn)
                    fda_data = [str(x) for x in fda_status.iloc[[index]].values.flatten().tolist()[1:]]
                    review.append('fda')
                    break

        #Data from EMA
        (status_ema, indication_ema) = ('', '')
        found = False
        for idx, syn in enumerate(synonyms_list):
            if syn in ema_status['International non-proprietary name (INN) / common name'].tolist():
                status_ema = list(set(ema_status.loc[ema_status['International non-proprietary name (INN) / common name'] == syn]['Authorisation status'].values.flatten().tolist()))
                indication_ema = [x.replace('\t',' ') for x in list(set(ema_status.loc[ema_status['International non-proprietary name (INN) / common name'] == syn]['Condition / indication'].values.flatten().tolist()))]
                if idx > 2: review.append('ema')
                found = True
                break
        if not found:
            for idx, syn in enumerate(synonyms_list):
                if syn in ema_status['International non-proprietary name (INN) / common name'].str.split(' ', n = 1, expand = True)[0].tolist():
                    index = ema_status['International non-proprietary name (INN) / common name'].str.split(' ', n = 1, expand = True)[0].tolist().index(syn)
                    status_ema = list(set(ema_status.iloc[[index]]['Authorisation status'].values.flatten().tolist()))
                    indication_ema = [x.replace('\t', ' ') for x in list(set(ema_status.iloc[[index]]['Condition / indication'].values.flatten().tolist()))]
                    review.append('ema')
                    break

        #Data from FDA labels
        indications = []
        applnos = []

        if fda_data[1] != '':
            applno = fda_data[1].split(';')

            for app in applno:
                app_label = ''
                if 'NDA'+app.zfill(6) in appl_list:
                    app_label = 'NDA'+app.zfill(6)
                elif 'ANDA'+app.zfill(6) in appl_list:
                    app_label = 'ANDA'+app.zfill(6)
                elif 'BLA'+app.zfill(6) in appl_list:
                    app_label = 'BLA'+app.zfill(6)
                if app_label != '':
                    for index2, row2 in fda_labels.loc[fda_labels['ApplNo'] == app_label].iterrows():
                        indications.append(row2.tolist()[0])
                        applnos.append(app)

        #Data from clinical trials
        ct_condition = ''
        found = False
        for idx, syn in enumerate(synonyms_list):
            if syn in status_CT_file['drug'].tolist():
                ct_condition = list(set(status_CT_file.loc[status_CT_file['drug'] == syn]['condition'].values.flatten().tolist()))
                if idx > 2: review.append('ct')
                found = True
                break
        if not found:
            for idx, syn in enumerate(synonyms_list):
                if syn in status_CT_file['drug'].str.split(' ', n = 1, expand = True)[0].tolist():
                    index = status_CT_file['drug'].str.split(' ', n = 1, expand = True)[0].tolist().index(syn)
                    ct_condition = list(set(status_CT_file.iloc[[index]]['drug'].values.flatten().tolist()))
                    review.append('ct')
                    break
        ct_condition = [str(x) for x in ct_condition]

        #Information from previous version
        prev_values = ['','','','','']
        if standard_drug_name in pandrugs_prev['standard_drug_name'].tolist():
            prev_values = pandrugs_prev.loc[pandrugs_prev['standard_drug_name'] == standard_drug_name][['status','pathology','cancer','extra','extra2']].values.tolist()[0]
        elif standard_drug_name in pandrugs_prev['show_drug_name'].tolist():
            prev_values = pandrugs_prev.loc[pandrugs_prev['show_drug_name'] == standard_drug_name][['status','pathology','cancer','extra','extra2']].values.tolist()[0]
        elif standard_drug_name in pandrugs_prev['source_drug_name'].tolist():
            prev_values = pandrugs_prev.loc[pandrugs_prev['source_drug_name'] == standard_drug_name][['status','pathology','cancer','extra','extra2']].values.tolist()[0]
        else:
            for syn in synonyms_list:
                if syn in pandrugs_prev['standard_drug_name'].tolist():
                    prev_values = pandrugs_prev.loc[pandrugs_prev['standard_drug_name'] == syn][['status','pathology','cancer','extra','extra2']].values.tolist()[0]
                    break
                elif syn in pandrugs_prev['show_drug_name'].tolist():
                    prev_values = pandrugs_prev.loc[pandrugs_prev['show_drug_name'] == syn][['status','pathology','cancer','extra','extra2']].values.tolist()[0]
                    break
                elif syn in pandrugs_prev['source_drug_name'].tolist():
                    prev_values = pandrugs_prev.loc[pandrugs_prev['source_drug_name'] == syn][['status','pathology','cancer','extra','extra2']].values.tolist()[0]
                    break

        #Type of therapy
        dtype = ''
        if re.search(('ib$|mab$|ib |mab '), standard_drug_name, re.IGNORECASE):
            dtype = 'TARGETED THERAPY'
        else:
            for syn in synonyms_list:
                if re.search(('ib$|mab$|ib |mab '), syn, re.IGNORECASE):
                    dtype = 'TARGETED THERAPY'
                    break
        fileo.write('\t'.join([standard_drug_name,fda_data[1],';'.join(applnos),fda_data[0],'::'.join(status_ema),'::'.join(ct_condition),'',prev_values[0],'::'.join(list(set(indications))),'::'.join(indication_ema),'',prev_values[3],'',prev_values[1],'',prev_values[2],dtype,prev_values[4],','.join(review)])+'\n')

    fileo.close()

def process_cgc_scores():

    print('Creating file for cancer type review for CGC...')

    cgc_dwn = 'cgc_scores.tsv'
    cgc_mr = 'cgc_scores_mr.tsv'

    cgc = pd.read_csv(pro_dir+cgc_dwn, sep ='\t', low_memory=False)
    cgc.fillna('', inplace=True)
    cgc['cancer_type'] = ''
 
    cgc.to_csv(pro_dir+cgc_mr, index=False, sep='\t', header=True)

synonyms_mr = 'drug_synonyms_mr.tsv'
#process_drug_names()

pandrugs_prev = pd.read_csv(pd_prev, sep ='\t', low_memory=False)
pandrugs_prev = pandrugs_prev.fillna('')
synonyms = pd.read_csv(mr_dir+synonyms_mr.replace('.tsv', 'd.tsv'), sep ='\t', low_memory=False)
synonyms['drug_name'] = synonyms['drug_name'].str.upper()
synonyms['standard_drug_name'] = synonyms['standard_drug_name'].str.upper()
synonyms['show_drug_name'] = synonyms['show_drug_name'].str.upper()

process_FDA_EMA() # run this after reviewing drug_names
#process_cgc_scores()
