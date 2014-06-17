import csv
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust
from numpy import linspace,ndarray,arange
from numpy.random import randn

remove_unmapped = False
just_ribosomes = False
use_LB = False
id_col_dict = { 'valgepea':'ko_num', 'heinmann':u'UP_AC' }
db_used = 'valgepea'

conf_fname_mod = '%s%s%s%s' % ('RibsOnly' if just_ribosomes else '', 'AnnotOnly' if remove_unmapped else '',"LB" if use_LB else '',db_used)

#Initialization of basic data containers, gene annotation data, growth rates and cell volumes and selection of conditions to analyze.
def ko_to_desc_dict():
    ko_annot_dict = {}

    with open('hierarchy_standardised.tms','rb') as annot_file:
         annots = csv.reader(annot_file,delimiter='\t') #,encoding='iso-8859-1')
         for row in annots:
             if len(row) == 2:
                cat = row [-1]
             elif len(row) == 3:
                subcat = row[-1]
             elif len(row) == 4:
                component = row[-1]
             elif len(row) == 5:
                ko_annot_dict[row[-1]]=(cat,subcat,component)
    return ko_annot_dict

def uni_ko_dict():
    uni_konum_dict = {}
    uni_to_konum = read_csv('eco_uniprot_mapping.csv',sep='[\t:]',encoding='iso-8859-1',header = None, names = ['ko','bla','uniprot'])
    for i,row in uni_to_konum.iterrows():
         uni_konum_dict[row['uniprot']]=row['ko']    
    return uni_konum_dict


def uniprot_to_desc_dict():
    uni_konum_dict = uni_ko_dict()
    #load the ko annotation tree:
    ko_annot_dict = ko_to_desc_dict()
    uni_to_annot = {}
    for uni in uni_konum_dict:
        if uni_konum_dict[uni] == 'NotMapped':
            uni_to_annot[uni]=['NotMapped']
        if uni_konum_dict[uni] == 'Not mapped':
            uni_to_annot[uni]=['NotMapped']
        elif uni_konum_dict[uni] not in ko_annot_dict:
	    uni_to_annot[uni]=['NotMapped']
        else:
            uni_to_annot[uni]=ko_annot_dict[uni_konum_dict[uni]]
    return uni_to_annot

def uniprot_to_offset():
    #load location information for genes:
    genome = SeqIO.read('U00096.2.gbk','genbank')

    locus_to_offset = {}
    for feat in genome.features:
        if feat.type == 'CDS':
           locus_to_offset[feat.qualifiers['locus_tag'][0]]=feat.location.start.real

    uniprot_to_locus = {}
    for row in open('all_ecoli_genes.txt','r'):
        uniprot_to_locus[row[48:54]]=row[0:5]

    uniprot_to_location = {}
    for uni in uniprot_to_locus.keys():
        if uniprot_to_locus[uni] in locus_to_offset.keys():
            uniprot_to_location[uni]= locus_to_offset[uniprot_to_locus[uni]]
        else:
            uniprot_to_location[uni]= 0
    return uniprot_to_location

# Define the list of conditions that will be relevant for the analysis, (and the description column), the growth rates and the cell volumes, according to the database used:
cond_list_dict = {'valgepea':[u'11', u'21', u'31', u'40', u'48'],
                  'heinmann':[
                      u'chemostat \u00b5=0.12', u'galactose',
                      u'chemostat \u00b5=0.20', u'acetate',
                      u'chemostat \u00b5=0.35', u'glucosamine',
                      u'pyruvate', u'glycerol', u'fumarate',
                      u'succinate', 
                      u'chemostat \u00b5=0.5',
                      u'anaerobic', u'glucose',],
                  'heinmann-chemo': [
                      u'chemostat \u00b5=0.12',
                      u'chemostat \u00b5=0.20',
                      u'chemostat \u00b5=0.35',
                      u'chemostat \u00b5=0.5']
                  }
if use_LB:
    cond_list_dict['heinmann'].append(u'LB')

gr_dict = {'valgepea':
    {u'11': 0.11, u'21':0.21, u'31':0.31, u'40':0.4, u'48':0.48},
           'heinmann':
    {u'chemostat \u00b5=0.12': 0.12, u'galactose':0.17, 
     u'chemostat \u00b5=0.20':0.2, u'acetate':0.29, 
     u'chemostat \u00b5=0.35':0.35, u'glucosamine':0.39, 
     u'pyruvate':0.4, u'glycerol':0.47, u'fumarate':0.47, 
     u'succinate':0.49, u'chemostat \u00b5=0.5':0.5, 
     u'anaerobic' : 0.55, u'glucose': 0.6, u'LB':1.61}
           }

def get_coli_data(db_used,use_weight):
    cond_list = cond_list_dict[db_used]
    if db_used == 'heinmann-chemo':
        db_used = 'heinmann'
    if db_used == 'heinmann':
        # As the file was exported from Excel, it uses Excel's encoding.
        ecoli_data = read_csv('coli_data.csv',header=1,encoding='iso-8859-1')

        #Split the data loaded into two sets - weight and count.
        idx = ecoli_data.columns
        ecoli_data_count = ecoli_data[idx[0:29]]
        ecoli_data_weight = ecoli_data[idx[0:10].append(idx[29:])]

        # Refine the DataFrames to include only these conditions (and the protein descriptions):
        desc_list = [id_col_dict[db_used]]
        count_cond = desc_list+cond_list
        ecoli_data_count = ecoli_data_count[count_cond]

        #duplicate headers are modified by read_csv and include a trailing '.1' string in their name.
        weight_cond = desc_list+[x+'.1' for x in cond_list] 
        ecoli_data_weight = ecoli_data_weight[weight_cond]

        #rename the columns to remove the trailing '.1'
        ecoli_data_weight.columns = count_cond

        #select the relevant data for analysis out of the two options:
        if use_weight:
            ecoli_data = ecoli_data_weight
        else:
            ecoli_data = ecoli_data_count

        #convert all data columns to floats, and description columns to strings.
        ecoli_data = ecoli_data[ecoli_data != 'below LOQ']
        ecoli_data[cond_list] = ecoli_data[cond_list].astype('float')

    if db_used == 'valgepea':
        ecoli_data = read_csv('valgepea.csv',header=0,encoding='iso-8859-1')

    #Normalize to get concentrations 
    ecoli_data[cond_list] = ecoli_data[cond_list] / ecoli_data[cond_list].sum()
    id_col = id_col_dict[db_used]
    ecoli_data[id_col] = ecoli_data[id_col].astype('string')
    return ecoli_data

def get_annotated_prots(db):
    coli_data = get_coli_data(db,use_weight=True)
    cond_list = cond_list_dict[db]
    if db == 'heinmann-chemo':
        db = 'heinmann'
    #annotate coli_data according to db.
    if db == 'heinmann':
        uni_to_konum = uni_ko_dict()
        coli_data['ko_num']=coli_data.apply(lambda x: 'NotMapped' if x[u'UP_AC'] not in uni_to_konum else uni_to_konum[x[u'UP_AC']],axis=1)
    id_to_annot = ko_to_desc_dict()
    id_col = 'ko_num'
    coli_data['group']=coli_data.apply(lambda x: 'NotMapped' if x[id_col] not in id_to_annot else (id_to_annot[x[id_col]])[0],axis=1)
    coli_data['func']=coli_data.apply(lambda x: 'NotMapped' if (x[id_col] not in id_to_annot) or (len(id_to_annot[x[id_col]]) < 3) else (id_to_annot[x[id_col]])[2],axis=1)

    if just_ribosomes:
        coli_data = coli_data[coli_data['func']=='Ribosome']
    if remove_unmapped:
        coli_data = coli_data[coli_data['group'] != 'NotMapped']

    coli_data = coli_data.dropna()
    gr = gr_dict[db]
    gr = pd.Series(gr)
    gr = gr[cond_list]
    return (cond_list,gr,coli_data)

def calc_gr_corr(df,cond_list,gr):
    df['gr_cov']=df[cond_list].apply(lambda x: x.corr(gr[cond_list]),axis=1)
    df['rsq']=df['gr_cov']**2
    return df

def add_loc_info(df):
    if db_used == 'heinmann':
        uni_to_loc = uniprot_to_offset()
        conc_data['loc']=conc_data.apply(lambda x: 0 if x[id_col_dict[db_used]] not in uni_to_loc else uni_to_loc[x[id_col_dict[db_used]]],axis=1)
