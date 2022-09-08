#!/usr/bin/python

import re
import pdb
import pandas as pd
from datetime import date
import progressbar
import math

#Directory with the downloaded files
dwn_dir = '/home/epineiro/Analysis/PanDrugs/2.0/downloads/'
#Directory with the processed files
pro_dir = '/home/epineiro/Analysis/PanDrugs/2.0/processed/'
#Directory with the manually reviewed files
mr_dir = '/home/epineiro/Analysis/PanDrugs/2.0/manual_review/'
#Directory with additional files
add_dir = '/home/epineiro/Analysis/PanDrugs/2.0/additional_files/'
#Output directory
out_dir = '/home/epineiro/Analysis/PanDrugs/2.0/'
#Output file
out_file = 'PanDrugs_'+date.today().strftime('%b_%d_%Y')+'.tsv'

#Dictionary with known targets - Mejor meter en un fichero de excepciones
knownT = { 'AKT1': ['MIRANSERTIB'], 'BCR': ['PONATINIB', 'PONATINIB HYDROCHLORIDE', 'AXITINIB', 'SARACATINIB', 'BOSUTINIB', 'DASATINIB', 'IMATINIB', 'IMATINIB MESYLATE', 'NILOTINIB', '923288-90-8'], 'ABL1': ['PONATINIB', 'PONATINIB HYDROCHLORIDE', 'AXITINIB', 'SARACATINIB', 'BOSUTINIB', 'DASATINIB', 'IMATINIB', 'IMATINIB MESYLATE', 'NILOTINIB', '923288-90-8'], 'KDR': ['SORAFENIB TOSYLATE'], 'PDGFRA': ['MASITINIB'], 'PIK3CA': ['OMIPALISIB', 'PILARALISIB', 'ZSTK474'], 'PIK3R1': ['OMIPALISIB', 'PILARALISIB'], 'PIK3R2': ['OMIPALISIB'], 'PML': ['ARSENIC TRIOXIDE', 'RETINOIC ACID'], 'RAF1': ['1096708-71-2', 'SORAFENIB TOSYLATE'], 'PIK3CB': ['OMIPALISIB'], 'EGFR': ['ERLOTINIB HYDROCHLORIDE', 'EGFR:AZD-9291 MESYLATE'], 'KRAS': ['2,6,10-TRIMETHYLDODECA-2,6,10-TRIENE'] }

def load_files():
    global checked_gene_symbol_file, drug_names_file, family_KEGG_file, family_cmap_file, drug_approved_data_file, ind_pathway_file, gene_dependency_file, moalmanac_sen_res, gdsc_sen_res, civic_sen_res, oncokb_sen_res, gene_names, gene_pathways

    print('Loading annotation files...')

    checked_gene_symbol_file = pd.read_csv(pro_dir+'genes_checked.tsv', sep ='\t', low_memory=False)
    checked_gene_symbol_file['gene_name'] = checked_gene_symbol_file['gene_name'].str.upper()

    drug_names_file = pd.read_csv(mr_dir+'drug_synonyms_mrd.tsv', sep ='\t', low_memory=False)
    drug_names_file['drug_name'] = drug_names_file['drug_name'].str.upper()
    drug_names_file['standard_drug_name'] = drug_names_file['standard_drug_name'].str.upper()
    drug_names_file['show_drug_name'] = drug_names_file['show_drug_name'].str.upper()
    drug_names_file = drug_names_file.fillna('')

    family_KEGG_file = pd.read_csv(pro_dir+'TargetBasedClassificationKEGG_formated.tsv', sep ='\t', low_memory=False)
    family_KEGG_file = family_KEGG_file.fillna('')
    family_KEGG_file['drug'] = family_KEGG_file['drug'].str.upper()
    family_cmap_file = pd.read_csv(pro_dir+'cmap_moa.tsv', sep ='\t', low_memory=False)
    family_cmap_file = family_cmap_file.fillna('')
    family_cmap_file['drug'] = family_cmap_file['drug'].str.upper()

    drug_approved_data_file = pd.read_csv(mr_dir+'drug_approved_data_mrd.tsv', sep='\t', low_memory=False)
    drug_approved_data_file = drug_approved_data_file.fillna('')

    ind_pathway_file = pd.read_csv(pro_dir+'KEGGmodeled/upstream_genes.tsv', sep='\t', low_memory=False)
    ind_pathway_file = ind_pathway_file.fillna('')
    ind_pathway_file['gene'] = ind_pathway_file['gene'].str.upper()

    gene_dependency_file = pd.read_csv(pro_dir+'genetic_dependencies.tsv', sep='\t', low_memory=False)
    gene_dependency_file = gene_dependency_file.fillna('')

    moalmanac_sen_res = pd.read_csv(pro_dir+'moalmanac.tsv', sep='\t', low_memory=False)
    moalmanac_sen_res = moalmanac_sen_res.fillna('')
    moalmanac_sen_res['drug_name'] = moalmanac_sen_res['drug_name'].str.upper()
    gdsc_sen_res = pd.read_csv(pro_dir+'GDSC.tsv', sep='\t', low_memory=False)
    gdsc_sen_res = gdsc_sen_res.fillna('')
    gdsc_sen_res['drug_name'] = gdsc_sen_res['drug_name'].str.upper()
    civic_sen_res = pd.read_csv(pro_dir+'civic.tsv', sep='\t', low_memory=False)
    civic_sen_res = civic_sen_res.fillna('')
    civic_sen_res['drug_name'] = civic_sen_res['drug_name'].str.upper()
    oncokb_sen_res = pd.read_csv(pro_dir+'oncokb.tsv', sep='\t', low_memory=False)
    oncokb_sen_res = oncokb_sen_res.fillna('')
    oncokb_sen_res['drug_name'] = oncokb_sen_res['drug_name'].str.upper()

    gene_names = pd.read_csv(dwn_dir+'custom', sep='\t', low_memory=False, dtype=str)
    gene_names = gene_names.fillna('')
    gene_pathways = pd.read_csv(pro_dir+'gene_pathway.tsv', sep='\t', low_memory=False, dtype=str)

def drug_gene_associations():

    print('Processing drug-gene associations...')
#CIVIC and OncoKB records are obtained from its web page in order to restrict to robut associations and avoid record of adverse efects
    files = ('DGIdb_interactions.tsv', 'oncokb.tsv', 'civic.tsv', 'moalmanac.tsv', 'GDSC.tsv')

    outputf = open(out_dir+out_file.replace('.tsv', '_prov.tsv'),'w')
    outputf.write('\t'.join(['gene_symbol', 'checked_gene_symbol', 'source', 'source_drug_name', 'standard_drug_name', 'show_drug_name', 'family', 'status', 'pathology', 'cancer', 'extra', 'extra2', 'pathways', 'target_marker', 'resistance', 'alteration', 'ind_pathway', 'gene_dependency', 'dscore', 'gscore', 'reviews'])+'\n')

    for fil in files:
        print('>Processing '+fil+'...')

        drug_gene_a = pd.read_csv(pro_dir+fil, sep='\t', low_memory=False)

        for i in progressbar.progressbar(range(len(drug_gene_a.index))):
            row = drug_gene_a.iloc[i]
            (gene_symbol,checked_gene_symbol, source, source_drug_name, standard_drug_name, show_drug_name, family, status, pathology, cancer, extra, extra2, pathways, target_marker, resistance, alteration, ind_pathway, gene_dependency, dscore, gscore, reviews) = ('', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', [])
            gene_symbol = re.split(';|,',row['gene_name'].upper())
            source = row['source']
            source_drug_name = row['drug_name'].upper()
            source_drug_name_reduced = source_drug_name.split(' ')[0]

            for gs in gene_symbol:
                #checked gene symbol
                checked_gene_symbol = checked_gene_symbol_file.loc[checked_gene_symbol_file['gene_name'] == gs]['checked_gene_symbol'].tolist()[0]

                standard_drug_name = drug_names_file.loc[drug_names_file['drug_name'] == source_drug_name]['standard_drug_name'].tolist()[0]
                standard_drug_name_reduced = standard_drug_name.split(' ')[0]

                #drug names file
                if drug_names_file.loc[drug_names_file['drug_name'] == source_drug_name]['show_drug_name'].tolist()[0] != '':
                    show_drug_name = drug_names_file.loc[drug_names_file['drug_name'] == source_drug_name]['show_drug_name'].tolist()[0]
                else:
                    show_drug_name = standard_drug_name

                show_drug_name_reduced = show_drug_name.split(' ')[0]
                drug_list = [standard_drug_name, show_drug_name, source_drug_name, standard_drug_name_reduced, show_drug_name_reduced, source_drug_name_reduced]

                #Drug_type
                extra2 = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['drug_type_mc'].tolist()[0]

                families = []
                #Drug family KEGG
                found = False
                for idx, d in enumerate(drug_list):
                    if d in family_KEGG_file['drug'].tolist():
                        families = families + [x+'(KEGG)' for x in family_KEGG_file.loc[family_KEGG_file['drug'] == d]['family2'].tolist()]
                        if idx > 2: reviews.append('KEGG')
                        found = True
                        break
                if not found:
                    for idx, d in enumerate(drug_list):
                        if d in family_KEGG_file['drug'].str.split(' ', n = 1, expand = True)[0].tolist():
                            index = family_KEGG_file['drug'].str.split(' ', n = 1, expand = True)[0].tolist().index(d)
                            families = families + [x+'(KEGG)' for x in family_KEGG_file.iloc[[index]]['family2'].tolist()]
                            reviews.append('KEGG')
                            break

                #Drug family cmap
                found = False
                for idx, d in enumerate(drug_list):
                    if d in family_cmap_file['drug'].tolist():
                        families = families + [x+'(Cmap)' for x in family_cmap_file.loc[family_cmap_file['drug'] == d]['moa'].tolist()]
                        if idx > 2: reviews.append('Cmap')
                        found = True
                        break
                if not found:
                    for idx, d in enumerate(drug_list):
                        if d in family_cmap_file['drug'].str.split(' ', n = 1, expand = True)[0].tolist():
                            index = family_cmap_file['drug'].str.split(' ', n = 1, expand = True)[0].tolist().index(d)
                            families = families + [x+'(Cmap)' for x in family_cmap_file.iloc[[index]]['moa'].tolist()]
                            reviews.append('Cmap')
                            break

                family = ', '.join(families)

                #Status
                status = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['status_mc'].tolist()[0]

                #Labels info
#                if source_drug_name == 'ERLOTINIB':
#                    pdb.set_trace() # probar con erlonitinb
                if status == 'Approved':
                    pathology = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['pathology_mc'].tolist()[0]
                    cancer = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['cancer_type_mc'].tolist()[0]
                    extra = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['indication_mc'].tolist()[0]

                #Target_marker
                if source in ['CancerCommons', 'ClearityFoundationClinicalTrial', 'DrugBank', 'MyCancerGenome', 'TALC', 'TEND', 'TTD']:
                    target_marker = 'target'
                else:
                    target_marker = 'marker'

                #Pathway members
                if len(ind_pathway_file.loc[ind_pathway_file['gene'] == gs]) > 0:
                    ind_pathway = ind_pathway_file.loc[ind_pathway_file['gene'] == gs]['upstream_genes'].tolist()[0]

                #Genetic dependencies
                if checked_gene_symbol in gene_dependency_file['gene'].tolist():
                    gene_dependency = gene_dependency_file[gene_dependency_file['gene'] == checked_gene_symbol]['genetic_dependency'].tolist()[0]

                #Resistance and alteration
                resistance = 'sensitivity'
                if source == 'CIViC':
                    if len(civic_sen_res[(civic_sen_res['drug_name'] == source_drug_name) & (civic_sen_res['gene_name'] == gs)]['response'].tolist()) > 0:
                        resistance = civic_sen_res[(civic_sen_res['drug_name'] == source_drug_name) & (civic_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = civic_sen_res[(civic_sen_res['drug_name'] == source_drug_name) & (civic_sen_res['gene_name'] == gs)]['variation'].tolist()[0]
                elif source == 'GDSC':
                    if len(gdsc_sen_res[(civic_sen_res['drug_name'] == source_drug_name) & (gdsc_sen_res['gene_name'] == gs)]['response'].tolist()) > 0:
                        resistance = gdsc_sen_res[(civic_sen_res['drug_name'] == source_drug_name) & (gdsc_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = gdsc_sen_res[(civic_sen_res['drug_name'] == source_drug_name) & (gdsc_sen_res['gene_name'] == gs)]['alteration'].tolist()[0]
                elif source == 'MOAlmanac':
                    if len(moalmanac_sen_res[(moalmanac_sen_res['drug_name'] == source_drug_name) & (moalmanac_sen_res['gene_name'] == gs)]['response'].tolist()[0]) > 0:
                        resistance = moalmanac_sen_res[(moalmanac_sen_res['drug_name'] == source_drug_name) & (moalmanac_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = moalmanac_sen_res[(moalmanac_sen_res['drug_name'] == source_drug_name) & (moalmanac_sen_res['gene_name'] == gs)]['alteration'].tolist()[0]
                elif source == 'OncoKB':
                    if len(oncokb_sen_res[(oncokb_sen_res['drug_name'] == source_drug_name) & (oncokb_sen_res['gene_name'] == gs)]['response'].tolist()) > 0:
                        resistance = oncokb_sen_res[(oncokb_sen_res['drug_name'] == source_drug_name) & (oncokb_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = oncokb_sen_res[(oncokb_sen_res['drug_name'] == source_drug_name) & (oncokb_sen_res['gene_name'] == gs)]['variation'].tolist()[0]

                #pathways
                if len(gene_names.loc[gene_names['Approved symbol']==checked_gene_symbol,]['NCBI Gene ID'].tolist()) > 0:
                    NCBI_name = gene_names.loc[gene_names['Approved symbol']==checked_gene_symbol,]['NCBI Gene ID'].tolist()[0]
                    if len(gene_pathways.loc[gene_pathways['KEGG Gene ID'] == NCBI_name]['KEGG Pathway ID'].tolist()) > 0:
                        pathways = gene_pathways.loc[gene_pathways['KEGG Gene ID'] == NCBI_name]['KEGG Pathway ID'].tolist()[0]

                outputf.write('\t'.join([gs, checked_gene_symbol, source, source_drug_name, str(standard_drug_name), str(show_drug_name), family, status, pathology, cancer, extra, str(extra2), pathways, target_marker, resistance, alteration, ind_pathway, gene_dependency, dscore, gscore, ';'.join(reviews)])+'\n')

    outputf.close()

    print('****IMPORTANT****'+'\n'+'Before creating the final PanDrugs file:')
    print('-> Check the records with review comments and modify them accordingly')
    print('-> Check the associations of KRAS, TP53, STK11 and APC (using checked_gene_symbol column) with drugs (using standard_drug_name column) with a target relationship (using target_marker column) in PanDrugs file. Update them in the control_records.tsv file along with the corresponding tag (keep/exclude/target/marker).')
    print('-> Recover EGFR (using checked_gene_symbol column) - drug (using standard_drug_name column) associations with a target relationship (using target_marker column) in PanDrugs file. For each record create a MET-drug association and update them in the control_records.tsv file along with the corresponding tag (amp_resistance).')
    print('Then re-run in update mode to create the final PanDrugs file.'+'\n'+'*****************')

def create_final_file():

        print('Creating final file...')

        pandrugs = pd.read_csv(out_dir+out_file.replace('.tsv', '_prov.tsv'), sep='\t', low_memory=False)
        pandrugs = pandrugs.fillna('')

        #update all target records
        target = pandrugs.loc[pandrugs['target_marker'] == 'target'].drop_duplicates(subset=['checked_gene_symbol', 'standard_drug_name'])

        for index, row in target.iterrows():
            pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name'])), 'target_marker'] = 'target'

        conres = pd.read_csv(add_dir+'controled_records.tsv', sep='\t', low_memory=False)

        for index, row in conres.loc[conres['reason'] == 'marker'].iterrows():
            pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name'])), 'target_marker'] = 'marker'

        for index, row in conres.loc[conres['reason'] == 'exclude'].iterrows():
            pandrugs = pandrugs.drop(pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name']))].index)

        for index, row in conres.loc[conres['reason'] == 'amp_resistance'].iterrows():
            dscore = -1
            pandrugs.loc[len(pandrugs.index)] = [row['checked_gene_symbol'], row['checked_gene_symbol'], 'Curated', row['standard_drug_name'], row['standard_drug_name']]+['','','','','','','','']+['marker','resistance', 'amplification']+['','','','']
#descomentar en cuanto los dscores están integrados
#            dscore = pandrugs.loc[pandrugs['standard_drug_name'] == row['standard_drug_name']][['dscore']].iloc[0].tolist() * -1
#            pandrugs.loc[len(pandrugs.index)] = [row['checked_gene_symbol'], row['checked_gene_symbol'], 'Curated', row['standard_drug_name'], row['standard_drug_name']]+pandrugs.loc[pandrugs['standard_drug_name'] == row['standard_drug_name']][['show_drug_name', 'family', 'status', 'pathology', 'cancer', 'extra', 'extra2', 'pathways']].iloc[0].tolist()+['marker','resistance', 'amplification']+pandrugs.loc[pandrugs['standard_drug_name'] == row['standard_drug_name']][['ind_pathway']].iloc[0].tolist()+dscore+pandrugs.loc[pandrugs['standard_drug_name'] == row['standard_drug_name']][['gscore']].iloc[0].tolist()+['']

        pandrugs.to_csv(out_dir+out_file, sep='\t', index=False, header=True)

def create_pubchem_ids_file():

    print('Creating PubChemID file...')

    pandrugs = pd.read_csv(out_dir+out_file, sep='\t', low_memory=False)
    pandrugs = pandrugs.fillna('')

    drugs = pandrugs['show_drug_name']

    outputf = open(out_dir+'show_drug_name_lk_pubchem.tsv', 'w')
    outputf.write('\t'.join(['show_drug_name', 'PubChemID'])+'\n')

    import pubchempy as pcp

    drugs = list(set(drugs))
    for i in progressbar.progressbar(range(len(drugs))):
        response = pcp.get_cids(drugs[i], 'name')
        CID = ''
        if len(response) > 0:
            CID = '|'.join([str(x) for x in response])
        outputf.write('\t'.join([drugs[i], CID])+'\n')

    outputf.close()

load_files()
drug_gene_associations()
#create_final_file()
#create_pubchem_ids_file()
