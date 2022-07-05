#!/usr/bin/python
#conda activate pandrugs2

#import os
import re
#import subprocess
import pdb
#import numpy as np
#import xml.etree.ElementTree as ET
#import gzip
#import wget
import pandas as pd
#import itertools
#from datetime import datetime
#import multiprocessing
#from multiprocessing import Pool
import progressbar

##Revisar ficheros y pathways

#Directorio con ficheros procesados
pro_dir = '/home/epineiro/Analysis/PanDrugs/2.0/processed/'

#Diccionario con registros a duplicar en fda_status.tsv
#Revisar la lista de nombres de drogas del fichero fda_status para localizar aquellas que requieran ser disgregadas.
fda_mr = {'TRIPLE SULFA (SULFABENZAMIDE;SULFACETAMIDE;SULFATHIAZOLE)': ['SULFABENZAMIDE', 'SULFACETAMIDE', 'SULFATHIAZOLE'], 'TRISULFAPYRIMIDINES (SULFADIAZINE;SULFAMERAZINE;SULFAMETHAZINE)': ['SULFADIAZINE', 'SULFAMERAZINE', 'SULFAMETHAZINE'], 'ETHINYL ESTRADIOL;NORETHINDRONE': ['ETHINYL ESTRADIOL', 'NORETHINDRONE'], 'PANCRELIPASE (AMYLASE;LIPASE;PROTEASE)': ['AMYLASE', 'LIPASE', 'PROTEASE'], 'CONJUGATED ESTROGENS/MEDROXYPROGESTERONE ACETATE': ['CONJUGATED ESTROGENS', 'MEDROXYPROGESTERONE ACETATE'], 'LAMIVUDINE, NEVIRAPINE, AND STAVUDINE': ['LAMIVUDINE', 'NEVIRAPINE', 'STAVUDINE'], 'DEXTROMETHORPHAN HYDROBROMIDE;QUINIDINE SULFATE': ['DEXTROMETHORPHAN HYDROBROMIDE', 'QUINIDINE SULFATE'], 'NORETHINDRONE ACETATE;ETHINYL ESTRADIOL;FERROUS FUMARATE': ['NORETHINDRONE ACETATE', 'ETHINYL ESTRADIOL', 'FERROUS FUMARATE'], 'ABACAVIR;LAMIVUDINE': ['ABACAVIR', 'LAMIVUDINE'], 'EFAVIRENZ;LAMIVUDINE;TENOFOVIR DISOPROXIL FUMARATE': ['EFAVIRENZ', 'LAMIVUDINE', 'TENOFOVIR DISOPROXIL FUMARATE'], 'SODIUM SULFATE;POTASSIUM SULFATE;MAGNESIUM SULFATE': ['SODIUM SULFATE', 'POTASSIUM SULFATE', 'MAGNESIUM SULFATE'], 'SITAGLIPTIN PHOSPHATE;METFORMIN HYDROCHLORIDE': ['SITAGLIPTIN PHOSPHATE', 'METFORMIN HYDROCHLORIDE'], 'LOPINAVIR;RITONAVIR': ['LOPINAVIR', 'RITONAVIR'], 'ABACAVIR SULFATE;LAMIVUDINE': ['ABACAVIR SULFATE', 'LAMIVUDINE'], 'SAXAGLIPTIN;METFORMIN HYDROCHLORIDE': ['SAXAGLIPTIN', 'METFORMIN HYDROCHLORIDE'], 'LAMIVUDINE;ZIDOVUDINE;EFAVIRENZ': ['LAMIVUDINE', 'ZIDOVUDINE', 'EFAVIRENZ'], 'ATAZANAVIR SULFATE;RITONAVIR;LAMIVUDINE;ZIDOVUDINE': ['ATAZANAVIR SULFATE', 'RITONAVIR', 'LAMIVUDINE', 'ZIDOVUDINE'], 'EMTRICITABINE;TENOFOVIR DISOPROXIL FUMARATE;NEVIRAPINE': ['EMTRICITABINE', 'TENOFOVIR DISOPROXIL FUMARATE', 'NEVIRAPINE'], 'DEXTROAMPHETAMINE SACCHARATE;AMPHETAMINE ASPARTATE MONOHYDRATE;DEXTROAMPHETAMINE SULFATE;AMPHETAMINE SULFATE': ['DEXTROAMPHETAMINE SACCHARATE', 'AMPHETAMINE ASPARTATE MONOHYDRATE', 'DEXTROAMPHETAMINE SULFATE', 'AMPHETAMINE SULFATE'], 'OMBITASVIR, PARITAPREVIR, RITONAVIR': ['OMBITASVIR', 'PARITAPREVIR', 'RITONAVIR'], 'NORETHINDRONE ACETATE;ETHINYL ESTRADIOL;ETHINYL ESTRADIOL;FERROUS FUMARATE': ['NORETHINDRONE ACETATE', 'ETHINYL ESTRADIOL', 'ETHINYL ESTRADIOL', 'FERROUS FUMARATE'], 'AMPHETAMINE ASPARTATE/DEXTROAMPHETAMINE SULFATE': ['AMPHETAMINE ASPARTATE', 'DEXTROAMPHETAMINE SULFATE'], 'LINAGLIPTIN;METFORMIN HYDROCHLORIDE': ['LINAGLIPTIN', 'METFORMIN HYDROCHLORIDE'], 'EMTRICITABINE;RILPIVIRINE;TENOFOVIR DISOPROXIL FUMARATE': ['EMTRICITABINE', 'RILPIVIRINE', 'TENOFOVIR DISOPROXIL FUMARATE'], 'LINAGLIPTIN AND METFORMIN HYDROCHLORIDE': ['LINAGLIPTIN', 'METFORMIN HYDROCHLORIDE'], 'METFORMIN HYDROCHLORIDE;SITAGLIPTIN PHOSPHATE': ['METFORMIN HYDROCHLORIDE', 'SITAGLIPTIN PHOSPHATE'], 'BENZOYL PEROXIDE;CLINDAMYCIN PHOSPHATE': ['BENZOYL PEROXIDE', 'CLINDAMYCIN PHOSPHATE'], 'DOLUTEGRAVIR;LAMIVUDINE;TENOFOVIR DISOPROXIL FUMARATE': ['DOLUTEGRAVIR', 'LAMIVUDINE', 'TENOFOVIR DISOPROXIL FUMARATE'], 'IVACAFTOR, TEZACAFTOR': ['IVACAFTOR', 'TEZACAFTOR'], 'BUDESONIDE;FORMOTEROL FUMARATE DIHYDRATE': ['BUDESONIDE', 'FORMOTEROL FUMARATE DIHYDRATE'], 'EMPAGLIFLOZIN;METFORMIN HYDROCHLORIDE': ['EMPAGLIFLOZIN', 'METFORMIN HYDROCHLORIDE'], 'ELEXACAFTOR, IVACAFTOR, TEZACAFTOR': ['ELEXACAFTOR', 'IVACAFTOR', 'TEZACAFTOR'], 'EMPAGLIFLOZIN;LINAGLIPTIN': ['EMPAGLIFLOZIN', 'LINAGLIPTIN'], 'ELAGOLIX SODIUM,ESTRADIOL,NORETHINDRONE ACETATE': ['ELAGOLIX SODIUM', 'ESTRADIOL', 'NORETHINDRONE ACETATE'], 'DARATUMUMAB;HYALURONIDASE-FIHJ': ['DARATUMUMAB', 'HYALURONIDASE-FIHJ'], 'PERTUZUMAB;TRASTUZUMAB;HYALURONIDASE-ZZXF': ['PERTUZUMAB', 'TRASTUZUMAB', 'HYALURONIDASE-ZZXF']}

def process_FDA():

    print('Correcting FDA file...')

    fda_status = pd.read_csv(pro_dir+'fda_status.tsv', sep ='\t', low_memory=False)
    fda_status.fillna('', inplace=True)

    fileo = open(pro_dir+'fda_status_mr.tsv', 'w')
    fileo.write('\t'.join(fda_status.columns.tolist())+'\n')

    for index, row in fda_status.iterrows():
        drug = row['drug'].rstrip()
        if drug in fda_mr.keys():
            for d in fda_mr[drug]:
                fileo.write('\t'.join([d]+[str(x) for x in row[1:].tolist()])+'\n')
        else:
            fileo.write('\t'.join([str(x) for x in row.tolist()])+'\n')

    fileo.close()

def process_FDA_labels():

    print('Creating FDA label file for review...')

    fda_status = pd.read_csv(pro_dir+'fda_status_mr.tsv', sep ='\t', low_memory=False)
    fda_status.fillna('', inplace=True)
 
    fda_labels = pd.read_csv(pro_dir+'fda_labels.tsv', sep ='\t', low_memory=False)
    fda_labels.fillna('', inplace=True)
    appl_list = fda_labels['ApplNo'].tolist()

    fileo = open(pro_dir+'fda_labels_mr.tsv', 'w')
    fileo.write('\t'.join(fda_labels.columns.tolist()+['ApplNo_status','indication_mc','pathology','cancer_type','therapy_type'])+'\n')

    for index, row in fda_status.iterrows():
        if row['status'] == 'Approved':
            applno = row['ApplNo'].split(';')
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
                        fileo.write('\t'.join(row2.tolist()+[app,'','','',''])+'\n')

    fileo.close()

def process_EMA():

    print('Creating EMA file for review...')

    ema_status = pd.read_csv(pro_dir+'ema_status.tsv', sep ='\t', low_memory=False)
    ema_status.fillna('', inplace=True)
 
    fileo = open(pro_dir+'ema_status_mr.tsv', 'w')
    fileo.write('\t'.join(['drug', 'status', 'condition/indication', 'indication_mc', 'pathology', 'cancer_type'])+'\n')

    for index, row in ema_status.iterrows():
        if row['Authorisation status'] in ['Authorised', 'Withdrawn']:
            if row['Authorisation status'] == 'Authorised':
                status = 'Approved'
            if row['Authorisation status'] == 'Withdrawn':
                status = 'Withdrawn'
        for comp in re.split(', | / ', row['International non-proprietary name (INN) / common name']):
            fileo.write('\t'.join([comp, status, row['Condition / indication'], '', '', ''])+'\n')

    fileo.close()

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

    outputf = open(pro_dir+'drug_synonyms_mr.tsv', 'w')
    outputf.write('\t'.join(['drug_name', 'standard_drug_name', 'show_drug_name', 'synonyms', 'review'])+'\n')

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
            outputf.write('\t'.join([drugs[i], response[0]['Synonym'][0], INNname, '::'.join(response[0]['Synonym']), ''])+'\n')
        else:
            outputf.write('\t'.join([drugs[i], drugs[i], drugs[i], '', ''])+'\n')

    outputf.close()

#    drug_file = pd.read_csv(pro_dir+'drug_synonyms_mr.tsv', sep ='\t', low_memory=False)
#    drug_file = drug_file.fillna('')

#    counts = 0
#    for index, row in drug_file.iterrows():
#        gene = row['drug_name']
#        sdn = row['standard_drug_name']
#        sdns = []
#        counts = counts + 1
#        print(counts)
#        for index2, row2 in drug_file.iterrows():
#            review = ''
#            sdn2 = row2['standard_drug_name']
#            synonyms = row2['synonyms'].split('::')
#            synonyms = [x.upper() for x in synonyms]
#            if gene.upper() in synonyms and sdn != sdn2:
#                sdns.append(sdn2)
#            if review != '':
#                pdb.set_trace()
#                drug_file[index,'review'] = 'review'
#                drug_file[index2,'review'] = 'review'
    
#    drug_file.to_csv(pro_dir+'drug_synonyms_mr.tsv', index=False, sep='\t', header=True)


def process_drug_type():

    print('Creating file for review with drug types for FDA approved label drugs...')

    fda_status = pd.read_csv(pro_dir+'fda_status_mr.tsv', sep ='\t', low_memory=False)
    fda_status.fillna('', inplace=True)
    pdb.set_trace()
    ema_status = pd.read_csv(pro_dir+'ema_status_mr.tsv', sep ='\t', low_memory=False)
    ema_status.fillna('', inplace=True)

    fileo = open(pro_dir+'drug_type_mr.tsv', 'w')
    fileo.write('\t'.join(['drug','drug_type'])+'\n')

    for index, row in fda_status.iterrows():
        drug = row['drug']
        if row['status'] == 'Approved':
            dtype = ''
            if re.search(('ib$|mab$|ib |mab '), drug, re.IGNORECASE):
                dtype = 'TARGETED THERAPY'
            fileo.write('\t'.join([drug, dtype])+'\n')

    for index, row in ema_status.iterrows():
        drug = row['drug']
        if row['status'] == 'Approved':
            dtype = ''
            if re.search(('ib$|mab$|ib |mab '), drug, re.IGNORECASE):
                dtype = 'TARGETED THERAPY'
            fileo.write('\t'.join([drug, dtype])+'\n')

    fileo.close()

#recoger los farmacos aprobados
#si el farmaco termina en ib o mab poner targeted therapy
#escribir en el fichero de salida de manual review

#process_FDA()
#process_FDA_labels()
#process_EMA()
process_drug_names()
#process_drug_type() # run this after review of FDA and EMA files
