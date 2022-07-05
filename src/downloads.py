#!/usr/bin/python
#conda activate pandrugs2

import os
import re
import subprocess
import pdb
import gzip
import wget
import pandas as pd
import progressbar

##Revisar ficheros y pathways

#Directorio de descargas
dwn_dir = '/home/epineiro/Analysis/PanDrugs/2.0/downloads/'

#URLs ficheros a descargar
gencode = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v39.annotation.gtf.gz'
gdsc = 'ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/ANOVA_results_GDSC2_20Feb20.xlsx'
features = 'https://www.cancerrxgene.org/downloads/download/genetic_feature'
kegg_ATC = 'https://www.genome.jp/kegg-bin/download_htext?htext=br08310&format=htext&filedir='
cmap = 'https://clue.io/fe0e0966-247c-4ad2-8258-498cf48acad5'
fda = 'https://www.fda.gov/media/89850/download'
fda_label = ['https://download.open.fda.gov/drug/label/drug-label-0001-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0002-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0003-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0004-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0005-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0006-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0007-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0008-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0009-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0010-of-0011.json.zip', 'https://download.open.fda.gov/drug/label/drug-label-0011-of-0011.json.zip']
ema = 'https://www.ema.europa.eu/sites/default/files/Medicines_output_european_public_assessment_reports.xlsx'
ct = 'https://clinicaltrials.gov/AllPublicXML.zip'
kegg_ind = 'http://rest.kegg.jp/get/$path/kgml' # $path sera sustituido por los nombres de los pathways
cgc = 'https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v95/cancer_gene_census.csv'
oncovar = 'https://oncovar.org/resource/download/All_genes_OncoVar_TCGA/TCGA.PanCancer.all.genes.OncoVar.tsv.gz'
therasabdab = 'http://opig.stats.ox.ac.uk/webapps/newsabdab/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv'
intogen = 'https://www.intogen.org/download?file=IntOGen-Drivers-20200201.zip'

#authentication string for COSMIC download: cambiar por el de una cuenta pandrugs
auth_str = ''

def download_DGIdb():
    print('Dowloading genecode annotations...')
    gencode_dwn = wget.download(gencode, out=dwn_dir)
    print(gencode_dwn)

    print('Obtaining list of human genes...')

    #Obtain the name of the downloaded file from gencode
    gencode_dwn = ''
    listdir = os.listdir(dwn_dir)
    for efile in listdir:
        if re.search('^gencode', efile) != None: gencode_dwn = efile

    #Retrieve the gene names
    genes_file = 'genes_fromGFT.tsv'
    genes = []
    inputf = gzip.open(dwn_dir+gencode_dwn,'rt')

    for line in inputf:
        line = line.rstrip("\n")
        line_a = line.split("\t")
        if re.search("##", line_a[0]) == None:
            names = line_a[8]
            pattern = re.compile('; gene_name "(.+?)";')
            gene = pattern.search(names)
            if gene.group(1) not in genes: genes.append(gene.group(1))
    inputf.close()

    print('Total number of genes: '+str(len(genes)))
    df = pd.DataFrame(genes)
    df.to_csv(dwn_dir+genes_file, index=False, sep='\t', header=False)
    print('Retrieving data from DGIdb...')
    #Por cada gen del listado recupero las asociaciones de la base de datos de DGIdb
    dgidb_dwn = 'DGIdb_interactions_dwn.tsv'
    open(dwn_dir+dgidb_dwn, 'w').close()
    filei = pd.read_csv(dwn_dir+genes_file, sep='\t', low_memory=False, header=None)
    genes = filei.iloc[:, 0].tolist()

#    import time
    import progressbar

    for i in progressbar.progressbar(range(len(genes))):
        command = "python ./python_example.py --genes='"+genes[i]+"' >> "+dwn_dir+dgidb_dwn
        subprocess.call(command, shell=True)

def download_therasabdab():

    print('Dowloading Thera-SAbDab annotations...')

    sabdab_dwn = wget.download(therasabdab, out=dwn_dir)
    print(sabdab_dwn)

def download_moalmanac():

    print('Dowloading moalmanac associations...')

    command = "curl -X 'GET' 'https://moalmanac.org/api/assertions' -H 'accept: */*' | jq -r '.[] | \"\\(.assertion_id)\t\\(.features[].attributes[].feature_type)\t\\(.features[].attributes[].gene)\t\\(.features[].attributes[].gene1)\t\\(.features[].attributes[].gene2)\t\\(.therapy_name)\t\\(.therapy_resistance)\t\\(.therapy_sensitivity)\t\\(.therapy_type)\t\\(.features[].attributes[].protein_change)\t\\(.predictive_implication)\t\\(.validated)\"' > "+dwn_dir+"moalmanac_dwn.tsv"
    subprocess.call(command,shell=True)

def download_GDSC():

    print('Dowloading GDSC associations and features...')

    gdsc_dwn = wget.download(gdsc, out=dwn_dir)
    print(gdsc_dwn)

    features_dwn = wget.download(features, out=dwn_dir+'GDSC_features.csv')
    print(features_dwn)

def download_KEGG_TB():

    print('Dowloading KEGG Target-based classification of drugs...')

    kegg_atc_dwn = wget.download(kegg_ATC, out=dwn_dir)
    print(kegg_atc_dwn)

def download_cmap():

    print('Downloading CLUE Repurposing data...')

    cmap_dwn = wget.download(cmap, out=dwn_dir)
    print(cmap_dwn)

def download_FDA():

    print('Downloading FDA data ...')

    fda_dwn = wget.download(fda, out=dwn_dir)
    print(fda_dwn)

def download_FDA_labels():

    print('Downloading FDA label data ...')

    for e in fda_label:
        fda_label_dwn = wget.download(e, out=dwn_dir)
        print(fda_label_dwn)

def download_EMA():

    print('Downloading EMA data ...')

    ema_dwn = wget.download(ema, out=dwn_dir)
    print(ema_dwn)

def download_ct():

    print('Downloading Clinical Trial records...')

    ct_dwn = wget.download(ct, out=dwn_dir)
    print(ct_dwn)

def download_KEGG_ind():

    print('Downloading KEGG pathways...')

    KEGG_dir = dwn_dir+'KEGGmodeled'

    if not os.path.exists(KEGG_dir): os.makedirs(KEGG_dir)

    pathways = ['hsa03320', 'hsa04010', 'hsa04012', 'hsa04014', 'hsa04015', 'hsa04020', 'hsa04022', 'hsa04024', 'hsa04062', 'hsa04064', 'hsa04066', 'hsa04068', 'hsa04071', 'hsa04110', 'hsa04114', 'hsa04115', 'hsa04150', 'hsa04151', 'hsa04152', 'hsa04210', 'hsa04261', 'hsa04270', 'hsa04310', 'hsa04330', 'hsa04340', 'hsa04350', 'hsa04370', 'hsa04390', 'hsa04510', 'hsa04520', 'hsa04530', 'hsa04540', 'hsa04611', 'hsa04620', 'hsa04621', 'hsa04622', 'hsa04630', 'hsa04650', 'hsa04660', 'hsa04662', 'hsa04664', 'hsa04666', 'hsa04668', 'hsa04670', 'hsa04722', 'hsa04910', 'hsa04912', 'hsa04914', 'hsa04915', 'hsa04916', 'hsa04919', 'hsa04920', 'hsa04921', 'hsa04922', 'hsa04971', 'hsa05010', 'hsa05012', 'hsa05160', 'hsa05200', 'hsa05205', 'hsa05212', 'hsa05214', 'hsa05218', 'hsa05231']

    for path in pathways:
        command = 'wget -O '+KEGG_dir+'/'+path+'.kgml.xml http://rest.kegg.jp/get/'+path+'/kgml'
        subprocess.call(command, shell=True)

def download_CGC():

    import ast

    print('Downloading CGC file...')

    command = 'curl -H "Authorization: Basic ' + auth_str + '" ' + cgc
    p = subprocess.Popen(['curl', '-H', 'Authorization: Basic ' + auth_str, cgc], stdout=subprocess.PIPE)
    out, err = p.communicate()

    cgc_dwn = wget.download(ast.literal_eval(out.decode('utf-8'))['url'], out=dwn_dir)
    print(cgc_dwn)

def download_oncovar():

    print('Downloading oncovar file...')

    oncovar_dwn = wget.download(oncovar, out=dwn_dir)
    print(oncovar_dwn)

def download_CIViC_evidence():

    print('Downloading CIViC evidence...')

    import requests
    import json

    url='https://civicdb.org/api/graphql'

    nodes = []

    query = """query evidenceItems($evidenceType: EvidenceType){
                 evidenceItems(evidenceType: $evidenceType) {
                 nodes {
                   id
                 }
                 pageCount
                 pageInfo {
                   endCursor
                   startCursor
                 }
                 totalCount
               }
            }"""

    variables = """{"evidenceType": "PREDICTIVE"}"""

    r = requests.post(url, json={'query': query, 'variables': variables})
    if r.status_code == 200:
        data = json.loads(r.text)
        for n in data['data']['evidenceItems']['nodes']:
            nodes.append(n['id'])

        query = """query evidenceItems($evidenceType: EvidenceType, $after: String){
                     evidenceItems(evidenceType: $evidenceType, after: $after) {
                     nodes {
                       id
                     }
                     pageCount
                     pageInfo {
                       endCursor
                       startCursor
                     }
                     totalCount
                   }
                }"""

        for i in range(1,data['data']['evidenceItems']['pageCount']):
            variables = """{"evidenceType": "PREDICTIVE", "after": \""""+data['data']['evidenceItems']['pageInfo']['endCursor']+"""\"}"""
            r = requests.post(url, json={'query': query, 'variables': variables})
            if r.status_code == 200:
                data = json.loads(r.text)
                for n in data['data']['evidenceItems']['nodes']:
                    nodes.append(n['id'])
            else:
                print('Error: ' + r.text)
    else:
        print('Error: ' + r.text)

    query = """query evidenceItem($id: Int!){
                 evidenceItem(id: $id) {
                 drugs { name }
                 gene { name }
                 variant { name }
                 status
                 clinicalSignificance
                 evidenceDirection
                 evidenceLevel
               }
            }"""

    print('length nodes: '+str(len(nodes)))
    count = 0
    output = open(dwn_dir+'civic_evidence.tsv','w')
    output.write('\t'.join(['drug', 'gene', 'variant', 'status', 'clinical_significance', 'evidence_direction', 'evidenceLevel'])+'\n')
    for n in progressbar.progressbar(range(len(nodes))):
        count = count + 1
        variables = """{"id": """+str(nodes[n])+"""}"""
        r = requests.post(url, json={'query': query, 'variables': variables})
        if r.status_code == 200:
            data = json.loads(r.text)
            gene = data['data']['evidenceItem']['gene']['name']
            variant = data['data']['evidenceItem']['variant']['name']
            status = data['data']['evidenceItem']['status']
            clinsig = data['data']['evidenceItem']['clinicalSignificance']
            evdirection = data['data']['evidenceItem']['evidenceDirection']
            evlevel = data['data']['evidenceItem']['evidenceLevel']
            drug = []
            for ele in data['data']['evidenceItem']['drugs']:
                drug.append(ele['name'])
            output.write('\t'.join([' + '.join(sorted(drug)), gene, variant, status, clinsig, evdirection, evlevel])+'\n')
        else:
            print('Error: ' + r.text)
    output.close()
    print('counts: '+str(count))

def download_oncoKB_evidence():

    print('Please, download oncoKB associations file from https://www.oncokb.org/actionableGenes and place it in downloads directory')

def download_intogen():

    print('Downloading IntOGen data...')

    print(intogen)
    intogen_dwn = wget.download(intogen, out=dwn_dir)
    print(intogen_dwn)

#download_DGIdb()
#download_therasabdab()
#download_moalmanac()
#download_GDSC()
#download_KEGG_TB()
#download_cmap()
#download_FDA()
#download_FDA_labels()
#download_EMA()
#download_ct()
#download_KEGG_ind()
#download_CGC()
#download_oncovar()
#download_CIViC_evidence()
#download_oncoKB_evidence()
download_intogen()
