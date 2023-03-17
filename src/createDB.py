#!/usr/bin/python

import re
import pdb
import pandas as pd
from datetime import date
import progressbar
import math
import sys

if len(sys.argv) < 2:
    print("python createDB.py level [being 1(create intermediate file), and 2(create final file)]")
    sys.exit()
else:
    print(sys.argv)

#Directory with the downloaded files
dwn_dir = 'downloads/'
#Directory with the processed files
pro_dir = 'processed/'
#Directory with the manually reviewed files
mr_dir = 'manual_review/'
#Directory with additional files
add_dir = 'additional_files/'
#Output directory
out_dir = '2.0/'
#Output file
out_file = 'PanDrugs_'+date.today().strftime('%b_%d_%Y')+'.tsv'

def load_files():
    global checked_gene_symbol_file, drug_names_file, family_KEGG_file, family_cmap_file, drug_approved_data_file, ind_pathway_file, gene_dependency_file, moalmanac_sen_res, gdsc_sen_res, civic_sen_res, oncokb_sen_res, gene_names, gene_pathways, cgc, intogen_group, intogen_min, intogen_max, oncovar, oncovar_min, oncovar_max, chronos, chronos_min, hallmarks

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

    #Files for gscore calculation
    cgc = pd.read_csv(pro_dir+'cgc_scores.tsv', sep='\t', low_memory=False, dtype=str)
    intogen = pd.read_csv(pro_dir+'intogen.tsv', sep='\t', low_memory=False, dtype=str)
    intogen_group = intogen.groupby('gene_name')['qvalue_combination'].median().reset_index()
    intogen_min = sorted(intogen_group['qvalue_combination'][intogen_group['qvalue_combination'] > 0.05].tolist())[0]
    intogen_max = min(intogen_group['qvalue_combination'])
    oncovar = pd.read_csv(pro_dir+'oncovar_scores.tsv', sep='\t', low_memory=False, dtype=str)
    oncovar_max = float(max(oncovar[(oncovar['cancer_type'] == 'PanCancer')]['Consensus_Score'].tolist()))
    chronos = pd.read_csv(pro_dir+'chronos_skew.tsv', sep='\t', low_memory=False, dtype=str, header = None)
    chronos.columns = ['gene', 'value', 'min']
    chronos['value'] = chronos['value'].astype(float)
    chronos_min = sorted(chronos[chronos['value'] > -0.5]['value'].tolist())[0]
    hallmarks = pd.read_csv(pro_dir+'hallmarks.tsv', sep='\t', low_memory=False, dtype=str)

def drug_gene_associations():

    print('Processing drug-gene associations...')
    #CIVIC and OncoKB records are obtained from its web page in order to restrict to robut associations and avoid record of adverse efects
    files = ('DGIdb_interactions.tsv', 'oncokb.tsv', 'civic.tsv', 'DrugBank.tsv', 'moalmanac.tsv', 'GDSC.tsv', 'sabdab.tsv')

    outputf = open(out_dir+'PanDrugs_prov.tsv','w')
    outputf.write('\t'.join(['gene_symbol', 'checked_gene_symbol', 'source', 'source_drug_name', 'standard_drug_name', 'show_drug_name', 'family', 'status', 'pathology', 'cancer', 'extra', 'extra2', 'pathways', 'target_marker', 'resistance', 'alteration', 'ind_pathway', 'gene_dependency', 'dscore', 'gscore', 'reviews'])+'\n')

    for fil in files:
        print('>Processing '+fil+'...')

        drug_gene_a = pd.read_csv(pro_dir+fil, sep='\t', low_memory=False)
        drug_gene_a['drug_name'] = drug_gene_a['drug_name'].str.upper()
        drug_gene_a = drug_gene_a.fillna('')
        drug_gene_a.drop_duplicates(inplace=True)

        for i in progressbar.progressbar(range(len(drug_gene_a.index))):
            row = drug_gene_a.iloc[i]

            (gene_symbol, source, source_drug_name) = ('', '', '')

            gene_symbol = re.split(';|,',row['gene_name'].upper())
            source = row['source']
            source_drug_name = row['drug_name'].upper()
            source_drug_name_reduced = source_drug_name.split(' ')[0]

            for gs in gene_symbol:

                (checked_gene_symbol, standard_drug_name, show_drug_name, family, status, pathology, cancer, extra, extra2, pathways, target_marker, resistance, alteration, ind_pathway, gene_dependency, dscore, gscore, reviews) = ('', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', [])

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

                if len(families) > 0: family = ', '.join(families)
                else: family = 'Other'

                #Status
                status = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['status_mc'].tolist()[0]

                #Labels info
                if status == 'Approved':
                    pathology = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['pathology_mc'].tolist()[0]
                    extra = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['indication_mc'].tolist()[0]

                cancer = drug_approved_data_file.loc[drug_approved_data_file['standard_drug_name'] == standard_drug_name]['cancer_type_mc'].tolist()[0]

                #Target_marker
                if source in ['CancerCommons', 'ClearityFoundationClinicalTrial', 'DrugBank', 'MyCancerGenome', 'TALC', 'TEND', 'TTD', 'SAbDab']:
                    target_marker = 'target'
                else:
                    target_marker = 'marker'

                #Pathway members
                if len(ind_pathway_file.loc[ind_pathway_file['gene'] == checked_gene_symbol]) > 0:
                    ind_pathway = ind_pathway_file.loc[ind_pathway_file['gene'] == checked_gene_symbol]['upstream_genes'].tolist()[0]

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
                    if len(gdsc_sen_res[(gdsc_sen_res['drug_name'] == source_drug_name) & (gdsc_sen_res['gene_name'] == gs)]['response'].tolist()) > 0:
                        resistance = gdsc_sen_res[(gdsc_sen_res['drug_name'] == source_drug_name) & (gdsc_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = gdsc_sen_res[(gdsc_sen_res['drug_name'] == source_drug_name) & (gdsc_sen_res['gene_name'] == gs)]['alteration'].tolist()[0]
                elif source == 'MOAlmanac':
                    if len(moalmanac_sen_res[(moalmanac_sen_res['drug_name'] == source_drug_name) & (moalmanac_sen_res['gene_name'] == gs)]['response'].tolist()[0]) > 0:
                        resistance = moalmanac_sen_res[(moalmanac_sen_res['drug_name'] == source_drug_name) & (moalmanac_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = moalmanac_sen_res[(moalmanac_sen_res['drug_name'] == source_drug_name) & (moalmanac_sen_res['gene_name'] == gs)]['alteration'].tolist()[0]
                elif source == 'OncoKB':
                    if len(oncokb_sen_res[(oncokb_sen_res['drug_name'] == source_drug_name) & (oncokb_sen_res['gene_name'] == gs)]['response'].tolist()) > 0:
                        resistance = oncokb_sen_res[(oncokb_sen_res['drug_name'] == source_drug_name) & (oncokb_sen_res['gene_name'] == gs)]['response'].tolist()[0]
                        alteration = oncokb_sen_res[(oncokb_sen_res['drug_name'] == source_drug_name) & (oncokb_sen_res['gene_name'] == gs)]['variation'].tolist()[0]
                elif source == 'COSMIC':
                    resistance = 'resistance'

                #pathways
                if len(gene_names.loc[gene_names['Approved symbol']==checked_gene_symbol,]['NCBI Gene ID'].tolist()) > 0:
                    NCBI_name = gene_names.loc[gene_names['Approved symbol']==checked_gene_symbol,]['NCBI Gene ID'].tolist()[0]
                    if len(gene_pathways.loc[gene_pathways['KEGG Gene ID'] == NCBI_name]['KEGG Pathway ID'].tolist()) > 0:
                        pathways = gene_pathways.loc[gene_pathways['KEGG Gene ID'] == NCBI_name]['KEGG Pathway ID'].tolist()[0]

                gscore = str(compute_gscore(checked_gene_symbol))

                outputf.write('\t'.join([gs, checked_gene_symbol, source, source_drug_name, str(standard_drug_name), str(show_drug_name), family, status, pathology, cancer, extra, str(extra2), pathways, target_marker, resistance, alteration, ind_pathway, gene_dependency, dscore, gscore, ';'.join(reviews)])+'\n')

    outputf.close()

    print('****IMPORTANT****'+'\n'+'Before creating the final PanDrugs file:')
    print('-> Check the records with review comments and modify them accordingly')
    print('-> Check the associations of KRAS, TP53, STK11 and APC (using checked_gene_symbol column) with drugs (using standard_drug_name column) with a target relationship (using target_marker column) in PanDrugs file. Update them in the control_records.tsv file along with the corresponding tag (keep/exclude/target/marker).')
    print('-> Recover EGFR (using checked_gene_symbol column) - drug (using standard_drug_name column) associations with a target relationship (using target_marker column) in PanDrugs file. For each record without any described targeted MET association with the drug, create a MET-drug association and update them in the control_records.tsv file along with the corresponding tag (amp_resistance).')
    print('Then re-run in create final file mode to create the final PanDrugs file.'+'\n'+'*****************')

def create_final_file():

    print('Creating final file...')

    pandrugs = pd.read_csv(out_dir+'PanDrugs_prov.tsv', sep='\t', low_memory=False)
    pandrugs = pandrugs.fillna('')

    #update all target records
    print('Updating controled records...')
    target = pandrugs.loc[pandrugs['target_marker'] == 'target'].drop_duplicates(subset=['checked_gene_symbol', 'standard_drug_name'])

    for index, row in target.iterrows():
        pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name'])), 'target_marker'] = 'target'

    conres = pd.read_csv(add_dir+'controled_records.tsv', sep='\t', low_memory=False)

    for index, row in conres.loc[conres['reason'] == 'marker'].iterrows():
        pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name'])), 'target_marker'] = 'marker'

    for index, row in conres.loc[conres['reason'] == 'exclude'].iterrows():
        pandrugs = pandrugs.drop(pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name']))].index)

    for index, row in conres.loc[conres['reason'] == 'amp_resistance'].iterrows():
        pandrugs = pandrugs.reset_index(drop=True)
        pandrugs.loc[len(pandrugs.index)] = [row['checked_gene_symbol'], row['checked_gene_symbol'], 'Curated', row['standard_drug_name'], row['standard_drug_name']]+pandrugs.loc[pandrugs['standard_drug_name'] == row['standard_drug_name']][['show_drug_name', 'family', 'status', 'pathology', 'cancer', 'extra', 'extra2']].iloc[0].tolist()+pandrugs.loc[pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']][['pathways']].iloc[0].tolist()+['marker', 'resistance', 'amplification']+pandrugs.loc[pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']][['ind_pathway', 'gene_dependency']].iloc[0].tolist()+['']+pandrugs.loc[pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']][['gscore']].iloc[0].tolist()+['']

    for index, row in conres.loc[conres['reason'] == 'sensitivity'].iterrows():
        pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name'])), 'resistance'] = 'sensitivity'
    for index, row in conres.loc[conres['reason'] == 'resistance'].iterrows():
        pandrugs.loc[((pandrugs['checked_gene_symbol'] == row['checked_gene_symbol']) & (pandrugs['standard_drug_name'] == row['standard_drug_name'])), 'resistance'] = 'resistance'

    #homogenize target-marker for salt correction
    print('Homogenizing target-marker information...')
    genes = list(set(pandrugs['checked_gene_symbol'].tolist()))

    for g in genes:
        subset = pandrugs.loc[pandrugs['checked_gene_symbol'] == g]
        drugs = list(set(subset['show_drug_name'].tolist()))
        for d in drugs:
            drug = d.split(' ')
            if len(drug) > 1 and drug[0] in drugs:
                if subset.loc[(subset['checked_gene_symbol'] == g) & (subset['show_drug_name'] == d)]['target_marker'].tolist()[0] != subset.loc[(subset['checked_gene_symbol'] == g) & (subset['show_drug_name'] == drug[0])]['target_marker'].tolist()[0]:
                    idx = pandrugs.index[(pandrugs['checked_gene_symbol'] == g) & (pandrugs['show_drug_name'] == d)].tolist()
                    pandrugs.loc[idx, ['target_marker']] = 'target'

                    idx = pandrugs.index[(pandrugs['checked_gene_symbol'] == g) & (pandrugs['show_drug_name'] == drug[0])].tolist()
                    pandrugs.loc[idx, ['target_marker']] = 'target'

    #homogenize drug families
    print('Homogenizing drug families...')
    drugs = list(set(pandrugs['show_drug_name'].tolist()))
    for d in drugs:
        f_values = list(set(pandrugs.loc[pandrugs['show_drug_name'] == d]['family'].tolist()))
        if len(f_values) > 1:
            idx = pandrugs.index[pandrugs['show_drug_name'] == d].tolist()
            pandrugs.loc[idx, ['family']] = [x for x in f_values if x != 'Other'][0]

    #Drug-score calculation
    for index, row in pandrugs.iterrows():
        dscore = 0
        if pandrugs.iloc[index]['status'] == 'Approved':
            if pandrugs.iloc[index]['extra2'] != '':
                if pandrugs.iloc[index]['target_marker'] == 'target': dscore = 1
                else: dscore = 0.9
            else:
                if pandrugs.iloc[index]['cancer'] == 'clinical cancer':
                    if pandrugs.iloc[index]['target_marker'] == 'target': dscore = 0.8
                    else: dscore = 0.7
                else:
                    if pandrugs.iloc[index]['target_marker'] == 'target': dscore = 0.4
                    else: dscore = 0.3
        if pandrugs.iloc[index]['status'] == 'Clinical trials':
            if pandrugs.iloc[index]['cancer'] == 'cancer':
                if pandrugs.iloc[index]['target_marker'] == 'target': dscore = 0.6
                else: dscore = 0.5
            else:
                if pandrugs.iloc[index]['target_marker'] == 'target': dscore = 0.2
                else: dscore = 0.1
        if pandrugs.iloc[index]['status'] == 'Experimental':
            if pandrugs.iloc[index]['target_marker'] == 'target': dscore = 0.0008
            else: dscore = 0.0004
        if pandrugs.iloc[index]['resistance'] == 'resistance':
            dscore = dscore * (-1)

        pandrugs.loc[index, 'dscore'] = dscore

    #remove duplicates due to repeated genes in original databases
    pandrugs.drop_duplicates(inplace=True)

    #add GScore and gene role
    gscore = pd.read_csv(add_dir+ "gscore_Ene_2023.tsv", sep='\t', low_memory=False)
    driver = pd.read_csv(add_dir+ "drivers_Ene_2023.tsv", sep='\t', low_memory=False)

    pandrugs_gs = pd.merge(pandrugs, gscore, on="checked_gene_symbol", how='left')
    pandrugs_gs_dr = pd.merge(pandrugs_gs, driver, on="checked_gene_symbol", how='left')

    pandrugs_gs_dr.to_csv(out_dir+out_file, sep='\t', index=False, header=True)

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

def compute_gscore(gene):
    gscore = 0
    if not cgc.loc[cgc['gene'] == gene].empty:
        if '1' in list(set(cgc.loc[cgc['gene'] == gene]['Tier'])): gscore += 0.2
        else: gscore += 0.1

    if not intogen_group.loc[intogen_group['gene_name'] == gene].empty:
        qvalue = intogen_group.loc[intogen_group['gene_name'] == gene]['qvalue_combination'].tolist()[0]
        if qvalue <= 0.05: gscore += 0.2 * (qvalue - intogen_min) / (intogen_max - intogen_min)

    if not oncovar.loc[(oncovar['cancer_type'] == 'PanCancer') & (oncovar['gene'] == gene)].empty:
        score = float(oncovar.loc[(oncovar['cancer_type'] == 'PanCancer') & (oncovar['gene'] == gene)]['Consensus_Score'].tolist()[0])
        if score <= 3:
            gscore += 0.1 * (score - 0) / (3 - 0)
        else:
            gscore += 0.1 + (0.1 * (score - 3) / (oncovar_max - 3))

    if not chronos.loc[chronos['gene'] == gene].empty:
        value = chronos.loc[chronos['gene'] == gene]['value'].tolist()[0]
        if value < -2: gscore += 0.2
        elif value <= -0.5: gscore += 0.2 * (value - chronos_min) / (-2 - chronos_min)

    if not hallmarks.loc[hallmarks['Genes_list'] == gene].empty:
        n_hall = len(list((set(hallmarks.loc[hallmarks['Genes_list'] == gene]['Hallmark']))))
        if n_hall >= 5: gscore += 0.2
        else: gscore += 0.2 * (n_hall - 0) / (5 - 0)

    return round(gscore, 4)

load_files()

if sys.argv[1] == '1':
    drug_gene_associations()
else:
    create_final_file()
    create_pubchem_ids_file()
