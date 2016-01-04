import csv
import pandas as pd
from pandas.io.parsers import read_csv
from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust,xscale
from numpy import linspace,ndarray,arange
from numpy.random import randn,shuffle,normal,seed

seed(123456)

remove_unmapped = False
just_ribosomes = False
use_LB = False
id_col_dict = { 'Valgepea':'B number identifier', 'Heinemann':u'UP_AC', 'Heinemann-chemo':u'UP_AC', 'HeinemannLB':u'UP_AC','Peebo':'B number identifier','Peebo-gluc':'B number identifier','Peebo-LB':'B number identifier','HuiAlim':'Gene','HuiClim':'Gene',  'HuiRlim':'Gene' }
db_used = 'Valgepea'
avg_conc_threshold = 0.00001

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

def name_to_desc_dict():
    n_annot_dict = {}

    with open('eco_hierarchy.tms','rb') as annot_file:
         annots = csv.reader(annot_file,delimiter='\t') #,encoding='iso-8859-1')
         for row in annots:
             if len(row) == 2:
                cat = (row [-1]).lower().capitalize()
             elif len(row) == 3:
                subcat = (row[-1]).lower().capitalize()
             elif len(row) == 4:
                component = (row[-1]).lower().capitalize()
             elif len(row) == 5:
                s=row[-1].split(':')[0]
                s = s[:1].lower()+s[1:] if s else ''
                n_annot_dict[s]=(cat,subcat,component)
    return n_annot_dict


def b_to_desc_dict():
    b_annot_dict = {}
    uni_to_b = uni_to_locus()[0]
    uni_to_cat = uniprot_to_category_dict()[0]
    for uni in uni_to_cat:
        if uni in uni_to_b:
            b_annot_dict[uni_to_b[uni]]=uni_to_cat[uni]
    return b_annot_dict

def uniprot_to_category_dict():
    uni_cat_dict = {}
    name_cat_dict = {}
    uni_name_dict = {}
    cats = set()
    df = read_csv('schmidt_prot_desc.csv')
    for i,row in df.iterrows():
        cat = row['Annotated functional COG class']
        if cat == '-' or cat == 'POORLY CHARACTERIZED':\
            cat = 'Unknown'
        else:
            cat = cat.lower().capitalize()
        gene = row['Gene']
        if "ribosomal protein" in row["Description"]:
            gene = "Ribosome"
        uni_cat_dict[row['Uniprot Accession']]=(cat,row['Annotated functional COG group (description)'],gene)
        name_cat_dict[row['Gene']]=(cat,row['Annotated functional COG group (description)'],gene)
        cats.add(cat)
        uni_name_dict[row['Uniprot Accession']]=row['Gene']
    print(cats)
    return uni_cat_dict,uni_name_dict,name_cat_dict

def uni_to_locus():
    uniprot_to_locus = {}
    name_to_uniprot = {}
    name_to_locus = {}
    locus_to_uniprot = {}
    for row in open('all_ecoli_genes.txt','r'):
        uniprot = row[48:54]
        locus = row[0:5]
        uniprot_to_locus[uniprot]=locus
        locus_to_uniprot[locus]=uniprot
        names = row[84:-1].split(';')
        for name in names:
            name_to_uniprot[name] = uniprot
            name_to_locus[name] = locus
    return (uniprot_to_locus,name_to_locus,name_to_uniprot,locus_to_uniprot)

# Define the list of conditions that will be relevant for the analysis, (and the description column), the growth rates and the cell volumes, according to the database used:
cond_list_dict = {'Valgepea':[u'0.11', u'0.21', u'0.31', u'0.4', u'0.49'],
                  'Heinemann':[
    u'Chemostat mu=0.12',
    u'Chemostat mu=0.20',
    u'Galactose',
    u'Acetate',
    u'Chemostat mu=0.35',
    u'Pyruvate',
    u'Fumarate',
    u'Succinate',
    u'Glucosamine',
    u'Glycerol',
    u'Mannose',
    u'Chemostat mu=0.5',
    u'Xylose',
    u'Osmotic-stress glucose',
    u'Glucose',
    u'pH6 glucose',
    u'Fructose',
    u'42C glucose'],
                  'HeinemannLB':[
    u'Chemostat mu=0.12',
    u'Chemostat mu=0.20',
    u'Galactose',
    u'Acetate',
    u'Chemostat mu=0.35',
    u'Pyruvate',
    u'Fumarate',
    u'Succinate',
    u'Glucosamine',
    u'Glycerol',
    u'Mannose',
    u'Chemostat mu=0.5',
    u'Xylose',
    u'Osmotic-stress glucose',
    u'Glucose',
    u'pH6 glucose',
    u'Fructose',
    u'42C glucose',
    u'Glycerol + AA',
    u'LB'
    ],
                  'Heinemann-chemo': [
    u'Chemostat mu=0.12',
    u'Chemostat mu=0.20',
    u'Chemostat mu=0.35',
    u'Chemostat mu=0.5'],
                  'Peebo-gluc': [
    u'0.21.1', u'0.31.1',
    u'0.41.1', u'0.51.2',
    u'0.22.3', u'0.26.1',
    u'0.36.1', u'0.46.1',
    u'0.51.3'],
                   'Peebo-LB': [
    u'0.22.4',
    u'0.25.1', u'0.35.1',
    u'0.45.1', u'0.55.1',
    u'0.65.1', u'0.74.1',
    u'0.82.1', u'0.22.5',
    u'0.42.1', u'0.53.1',
    u'0.63.1', u'0.73.1',
    u'0.78.1'],
                    'Peebo': [
    u'0.21.1', u'0.31.1',
    u'0.41.1', u'0.51.2',
    u'0.22.3', u'0.26.1',
    u'0.36.1', u'0.46.1',
    u'0.51.3', u'0.22.4',
    u'0.25.1', u'0.35.1',
    u'0.45.1', u'0.55.1',
    u'0.65.1', u'0.74.1',
    u'0.82.1', u'0.22.5',
    u'0.42.1', u'0.53.1',
    u'0.63.1', u'0.73.1',
    u'0.78.1'],
                  'HuiClim': [
    0.45205,0.57762,0.67079,0.86643,1.0397 ],
                  'HuiAlim': [
    0.45702, 0.60274, 0.71705, 0.88487, 0.96718],
                  'HuiRlim': [
    0.28292, 0.40773, 0.63983, 0.99021]
                      }
gr_dict = {'Valgepea': {
    u'0.11': 0.11, u'0.21':0.21, u'0.31':0.31, u'0.4':0.4, u'0.49':0.49},
           'Heinemann': {
    u'Chemostat mu=0.12': 0.12,
    u'Chemostat mu=0.20':0.2,
    u'Galactose':0.26,
    u'Acetate':0.3,
    u'Chemostat mu=0.35':0.35,
    u'Pyruvate':0.4,
    u'Fumarate':0.42,
    u'Succinate':0.44,
    u'Glucosamine':0.46,
    u'Glycerol':0.47,
    u'Mannose':0.47,
    u'Chemostat mu=0.5':0.5,
    u'Xylose':0.55,
    u'Osmotic-stress glucose':0.55,
    u'Glucose': 0.58,
    u'pH6 glucose':0.63,
    u'Fructose':0.65,
    u'42C glucose':0.66,
    },
           'HeinemannLB': {
    u'Chemostat mu=0.12': 0.12,
    u'Chemostat mu=0.20':0.2,
    u'Galactose':0.26,
    u'Acetate':0.3,
    u'Chemostat mu=0.35':0.35,
    u'Pyruvate':0.4,
    u'Fumarate':0.42,
    u'Succinate':0.44,
    u'Glucosamine':0.46,
    u'Glycerol':0.47,
    u'Mannose':0.47,
    u'Chemostat mu=0.5':0.5,
    u'Xylose':0.55,
    u'Osmotic-stress glucose':0.55,
    u'Glucose': 0.58,
    u'pH6 glucose':0.63,
    u'Fructose':0.65,
    u'42C glucose':0.66,
    u'Glycerol + AA':1.27,
    u'LB':1.9,
    },
             'Peebo-gluc': {
    u'0.21.1':0.21, u'0.31.1':0.31,
    u'0.41.1':0.41, u'0.51.2':0.51,
    u'0.22.3':0.22, u'0.26.1':0.26,
    u'0.36.1':0.36, u'0.46.1':0.46,
    u'0.51.3':0.51},
              'Peebo-LB': {
    u'0.22.4':0.22,
    u'0.25.1':0.25, u'0.35.1':0.35,
    u'0.45.1':0.45, u'0.55.1':0.55,
    u'0.65.1':0.65, u'0.74.1':0.74,
    u'0.82.1':0.82, u'0.22.5':0.22,
    u'0.42.1':0.42, u'0.53.1':0.53,
    u'0.63.1':0.63, u'0.73.1':0.73,
    u'0.78.1':0.78
    },
              'Peebo': {
    u'0.21.1':0.21, u'0.31.1':0.31,
    u'0.41.1':0.41, u'0.51.2':0.51,
    u'0.22.3':0.22, u'0.26.1':0.26,
    u'0.36.1':0.36, u'0.46.1':0.46,
    u'0.51.3':0.51, u'0.22.4':0.22,
    u'0.25.1':0.25, u'0.35.1':0.35,
    u'0.45.1':0.45, u'0.55.1':0.55,
    u'0.65.1':0.65, u'0.74.1':0.74,
    u'0.82.1':0.82, u'0.22.5':0.22,
    u'0.42.1':0.42, u'0.53.1':0.53,
    u'0.63.1':0.63, u'0.73.1':0.73,
    u'0.78.1':0.78
    },
              'HuiClim': {
    0.45205:0.45205,0.57762:0.57762,
    0.67079:0.67079,0.86643:0.86643,
    1.0397:1.0397 },
              'HuiAlim': {
    0.45702:0.45702, 0.60274:0.60274,
    0.71705:0.71705, 0.88487:0.88487,
    0.96718:0.96718 },
              'HuiRlim': {
    0.28292:0.28292, 0.40773:0.40773,
    0.63983:0.63983, 0.99021:0.99021}
   }

def get_coli_data(db_used,use_weight,rand):
    cond_list = cond_list_dict[db_used]
    if db_used == 'Heinemann-chemo':
        db_used = 'Heinemann'
    if db_used == 'Heinemann' or db_used == 'HeinemannLB':
        # As the file was exported from Excel, it uses Excel's encoding.
        ecoli_data = read_csv('matthias2.csv',header=1,encoding='iso-8859-1')
       # ecoli_data = read_csv('coli_data.csv',header=1,encoding='iso-8859-1')

        # Refine the DataFrames to include only these conditions (and the protein descriptions):
        desc_list = [id_col_dict[db_used]]
        count_cond = desc_list+cond_list
        ecoli_data_count = ecoli_data[count_cond]

        #duplicate headers are modified by read_csv and include a trailing '.1' string in their name.
        weight_cond = desc_list+[x+'.1' for x in cond_list]
        ecoli_data_weight = ecoli_data[weight_cond]

        cv_cond = desc_list+[x+'.2' for x in cond_list]
        ecoli_data_cv = ecoli_data[cv_cond]
        #rename the columns to remove the trailing '.1'
        ecoli_data_weight.columns = count_cond

        ecoli_data_cv.columns = count_cond

        #select the relevant data for analysis out of the two options:
        if use_weight:
            ecoli_data = ecoli_data_weight
        else:
            ecoli_data = ecoli_data_count

        #convert all data columns to floats, and description columns to strings.
        #ecoli_data = ecoli_data[ecoli_data != 'below LOQ']
        ecoli_data[cond_list] = ecoli_data[cond_list].astype('float')

    if db_used == 'HuiAlim':
        ecoli_data = pd.read_excel('hui.xlsx',0,skiprows=4)
    if db_used == 'HuiClim':
        ecoli_data = pd.read_excel('hui.xlsx',0,skiprows=4)
    if db_used == 'HuiRlim':
        ecoli_data = pd.read_excel('hui.xlsx',0,skiprows=4)
    if db_used == 'Valgepea':
        ecoli_data = read_csv('ecoli_Valgepea_et_al_2013.csv',header=0,encoding='iso-8859-1')
    if db_used == 'Peebo':
        ecoli_data = read_csv('valgepea2.csv',header=0,encoding='iso-8859-1')
    if db_used == 'Peebo-gluc':
        ecoli_data = read_csv('valgepea2.csv',header=0,encoding='iso-8859-1')
    if db_used == 'Peebo-LB':
        ecoli_data = read_csv('valgepea2.csv',header=0,encoding='iso-8859-1')

    #for rand_method in ["simulated"]:
    ## Randomize rows:
    if rand == "shuffle":
        for i in ecoli_data.index:
            y = ecoli_data[cond_list].loc[i]
            shuffle(y.values)
            ecoli_data.loc[i,cond_list] = y
    #Normalize to get concentrations
    if db_used not in ['HuiAlim','HuiClim','HuiRlim']:
        ecoli_data[cond_list] = ecoli_data[cond_list] / ecoli_data[cond_list].sum()
    #remove scarce proteins
    if db_used == 'Peebo' or db_used == 'Peebo-gluc' or db_used == 'Peebo-LB':
        ecoli_data = ecoli_data[ecoli_data[cond_list].mean(axis=1)>avg_conc_threshold]
    if db_used == 'Heinemann' or db_used == 'HeinemannLB':
        dropping = ecoli_data[ecoli_data_cv[cond_list].mean(axis=1)>20.0][cond_list].sum().mean()
        keeping = ecoli_data[ecoli_data_cv[cond_list].mean(axis=1)<20.0][cond_list].sum().mean()
        plottemp = ecoli_data[ecoli_data_cv[cond_list].mean(axis=1)>0]
        print("dropping:%f" % dropping)
        print("keeping:%f" % keeping)
        figure(figsize=(5,3))
        #hist(ecoli_data_cv[cond_list].mean(axis=1).dropna(),20)
        plot(plottemp[cond_list].mean(axis=1),ecoli_data_cv[cond_list].mean(axis=1).dropna(),'.')
        xscale('log')
        xlabel("mean fraction out of proteome")
        ylabel("mean CV")
        savefig('CVmean.pdf')
        ecoli_data = ecoli_data[ecoli_data_cv[cond_list].mean(axis=1)<20.0]
        for cond in cond_list:
            ecoli_data["%s.cv" % cond] = ecoli_data_cv[cond]

    #remove irrelevant proteins
    means = ecoli_data[cond_list].mean(axis=1)
    ecoli_data = ecoli_data[means>0]
    #renormalize
    if db_used not in ['HuiAlim','HuiClim','HuiRlim']:
        ecoli_data[cond_list] = ecoli_data[cond_list] / ecoli_data[cond_list].sum()
    ## create emulated data set based on actual average concentrations and noise.
    if rand == "simulated":
        ecoli_data['avg'] = ecoli_data[cond_list].mean(axis=1)
        ecoli_data = ecoli_data.sort('avg',ascending=False)
        # assume even entries scale and odd entries dont.
        gr = gr_dict[db_used]
        gr = pd.Series(gr)
        gr = gr[cond_list]
        response = gr/(2*gr.mean())+1.0/2
        print response
        even = True
        for i,row in ecoli_data.iterrows():
            avg = ecoli_data.ix[i,'avg']
            noise = normal(loc=0,scale = avg*0.25,size=len(cond_list)) # assume 20% noise
            if even:
                ecoli_data.ix[i,cond_list] = avg+noise
            else:
                ecoli_data.ix[i,cond_list] = avg*response+noise
            even = not even
# not renormalizing as that creates anti-correlation effect.
    id_col = id_col_dict[db_used]
    ecoli_data[id_col] = ecoli_data[id_col].astype('string')
    return ecoli_data

def get_annotated_prots(db,rand):
    coli_data = get_coli_data(db,use_weight=True,rand=rand)
    cond_list = cond_list_dict[db]

    uniprot_to_locus,name_to_locus,name_to_uniprot,locus_to_uni = uni_to_locus()
    uni_to_cat,uni_to_name,name_to_cat = uniprot_to_category_dict()
    if db == 'Heinemann-chemo':
        db = 'Heinemann'
    #annotate coli_data according to db.
    if db == 'Heinemann' or db == 'HeinemannLB':
        ##uni_to_konum = uni_ko_dict()
        x=0
        with open('unmappeduni.txt','w+') as f:
            for i,r in coli_data.iterrows():
                if r[u'UP_AC'] not in uniprot_to_locus or uniprot_to_locus[r[u'UP_AC']] == 'NotMapped':
                    f.write(r[u'UP_AC'])
                    f.write('\n')
                    x+=1
            print "unmapped uniprots:%d" % x
        ##coli_data['ko_num']=coli_data.apply(lambda x: 'NotMapped' if x[u'UP_AC'] not in uni_to_konum else uni_to_konum[x[u'UP_AC']],axis=1)
        coli_data['protName']=coli_data.apply(lambda x: 'NotMapped' if x[u'UP_AC'] not in uni_to_name else uni_to_name[x[u'UP_AC']],axis=1)
        id_to_annot = uni_to_cat
    if db == 'Valgepea':
        id_to_annot = b_to_desc_dict()
    if db == 'Peebo':
        id_to_annot = b_to_desc_dict()
    if db == 'Peebo-gluc':
        id_to_annot = b_to_desc_dict()
    if db == 'Peebo-LB':
        id_to_annot = b_to_desc_dict()
    if db == 'HuiAlim':
        id_to_annot = name_to_cat
    if db == 'HuiClim':
        id_to_annot = name_to_cat
    if db == 'HuiRlim':
        id_to_annot = name_to_cat
    id_col = id_col_dict[db]
    x=0
    y=0
    with open('unmappedko.txt','w+') as f:
        for i,r in coli_data.iterrows():
            if r[id_col] not in id_to_annot:
                if r[id_col] != "NotMapped":
                    f.write(r[id_col])
                    f.write("\n")
                    x+=1
                y+=1
    print "total unmapped ko's:%d" % x
    print "total unmapped :%d" % y
    coli_data['group']=coli_data.apply(lambda x: 'Unknown' if x[id_col] not in id_to_annot else (id_to_annot[x[id_col]])[0],axis=1)
    coli_data['func']=coli_data.apply(lambda x: 'NotMapped' if (x[id_col] not in id_to_annot) or (len(id_to_annot[x[id_col]]) < 3) else (id_to_annot[x[id_col]])[1],axis=1)
    coli_data['prot']=coli_data.apply(lambda x: 'NotMapped' if (x[id_col] not in id_to_annot) or (len(id_to_annot[x[id_col]]) < 3) else (id_to_annot[x[id_col]])[2],axis=1)
    coli_data['ID']=coli_data.apply(lambda x: 'NotMapped' if (x[id_col] not in id_to_annot) or (len(id_to_annot[x[id_col]]) < 3) else x[id_col],axis=1)
    coli_data['Temp']=coli_data[id_col_dict[db]]

    if just_ribosomes:
        coli_data = coli_data[coli_data['func']=='Ribosome']
    if remove_unmapped:
        coli_data = coli_data[coli_data['group'] != 'NotMapped']
    gr = gr_dict[db]
    gr = pd.Series(gr)
    gr = gr[cond_list]
    return (cond_list,gr,coli_data)

def calc_gr_corr(df,cond_list,gr):
    df['gr_cov']=df[cond_list].apply(lambda x: x.corr(gr[cond_list]),axis=1)
    df['rsq']=df['gr_cov']**2
    return df
