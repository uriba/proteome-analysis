import csv
import pandas as pd
import scipy.odr as od
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
from Bio import SeqIO

#Initialization of basic data containers, gene annotation data, growth rates and cell volumes and selection of conditions to analyze.
def uniprot_to_desc_dict():
    uni_konum_dict = {}
    uni_to_konum = read_csv('eco_uniprot_mapping.csv',sep='[\t:]',encoding='iso-8859-1',header = None, names = ['ko','bla','uniprot'])
    for i,row in uni_to_eco.iterrows():
         unidict[row['uniprot']]=row['ko']    

    #load the ko annotation tree:
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
                annot_dict[row[-1]]=(cat,subcat,component)
    uni_to_annot = {}
    for uni in uni_konum_dict:
        if uni_konum_dict[uni] == 'NotMapped':
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
    return uniprot_to_location

# define the growth rates:
gr = {
    u'chemostat \u00b5=0.12': 0.12, 
    u'galactose':0.17, 
    u'chemostat \u00b5=0.20':0.2, 
    u'acetate':0.29, 
    u'chemostat \u00b5=0.35':0.35,
    u'glucosamine':0.39, 
    u'pyruvate':0.4, 
    u'glycerol':0.47,
    u'fumarate':0.47, 
    u'succinate':0.49, 
    u'chemostat \u00b5=0.5':0.5, 
    u'anaerobic' : 0.55, 
    u'glucose': 0.6,
    u'LB':1.61
    }
gr = pd.Series(gr)
# define cell volumes:
volumes = {u'glucose': 3.2, u'anaerobic' : 2.9, u'chemostat \u00b5=0.12': 2.1, u'chemostat \u00b5=0.5':2.6, u'galactose':1.9, u'acetate':2.4, u'glycerol':2.3,
     u'LB':4.4,
     u'glucosamine':2.9, u'fumarate':2.4, u'succinate':2.4, u'pyruvate':2.1, u'chemostat \u00b5=0.20':2.2, u'chemostat \u00b5=0.35':2.4}
volumes = pd.Series(volumes)

# first load the data from the file.
# As the file was exported from Excel, it uses Excel's encoding.
ecoli_data = read_csv('coli_data.csv',header=1,encoding='iso-8859-1')

#Split the data loaded into two sets - weight and count.
idx = ecoli_data.columns
ecoli_data_weight = ecoli_data[idx[0:10].append(idx[29:])]
ecoli_data_count = ecoli_data[idx[0:29]]

# Define the list of conditions that will be relevant for the analysis, (and the description column):
cond_list = [
    u'chemostat \u00b5=0.12',
    u'galactose',
    u'chemostat \u00b5=0.20',
    u'acetate',
    u'chemostat \u00b5=0.35',
    u'glucosamine',
    u'pyruvate',
    u'glycerol',
    u'fumarate',
    u'succinate',
    u'chemostat \u00b5=0.5',
    u'anaerobic',
    u'glucose',
#    u'LB'
]
#Define the text fields that will be relevant for the analysis:
desc_list = ['Description','UP_AC']

# Refine the DataFrames to include only these conditions (and the protein descriptions):
count_cond = desc_list+cond_list
ecoli_data_count = ecoli_data_count[count_cond]

#duplicate headers are modified by read_csv and include a trailing '.1' string in their name.
weight_cond = desc_list+[x+'.1' for x in cond_list] 
ecoli_data_weight = ecoli_data_weight[weight_cond]

#rename the columns to remove the trailing '.1'
ecoli_data_weight.columns = count_cond

#select the relevant data for analysis out of the two options:
ecoli_data = ecoli_data_weight

#convert all columns, excluding 'Description', to floats, and 'Description' to string.
ecoli_data_weight = ecoli_data_weight[ecoli_data_weight != 'below LOQ']
ecoli_data_weight[cond_list] = ecoli_data_weight[cond_list].astype('float')
ecoli_data_weight[desc_list] = ecoli_data_weight[desc_list].astype('string')


ecoli_data = ecoli_data[ecoli_data != 'below LOQ']
ecoli_data[cond_list] = ecoli_data[cond_list].astype('float')
ecoli_data[desc_list] = ecoli_data[desc_list].astype('string')

#Normalize to get concentrations (Use the total protein weight as the normalizing factor to avoid errors in cell volume measurements)
ecoli_data[cond_list] = ecoli_data[cond_list] / ecoli_data_weight[cond_list].sum() #volumes[cond_list]

#get the ribosomal proteins.
ecoli_data_values = ecoli_data.drop(desc_list,1)
ecoli_data_values = ecoli_data_values.dropna()

ribs = ecoli_data['Description'].map(lambda x: 'ribosomal protein' in x)
ribs_data = ecoli_data_values[ribs]
ribs_avg = ribs_data.sum()
gr = gr[ribs_avg.index]
volumes = volumes[ribs_avg.index]
########################### Analysis of correlation of proteins to growth rate.
xs = linspace(-1,1,100)

covs = []
catcovs = {}
prots = {}
high_cor_prots = {}
threshold = 0.3
hc_indices = []
hc_prots = pd.DataFrame(columns = ecoli_data_values.columns)

loc_prots = {0:[],1:[],2:[]}
loc_indices = {0:[],1:[],2:[]}
loc_data = {}

for row_index, row in ecoli_data_values.iterrows():
    desc = []
    # calculate the correlation and add it to the statistics
    cov = row.cov(gr)/(row.std() * gr.std())
    # get the protein annotation information
    eco = ecoli_data['UP_AC'][row_index]
    uniprot = unidict[eco]
    if eco in uniprot_to_location:
	loc_prots[(uniprot_to_location[eco]/1600000)].append(cov)
	loc_indices[(uniprot_to_location[eco]/1600000)].append(row_index)
	loc_data[row_index]= uniprot_to_location[eco]
    if uniprot in annot_dict:
	desc = annot_dict[uniprot]
	covs.append(cov)
	# set the counters of all the annotation fields that this protein belong to.
	for i in range(0,3):
	    if desc[i] not in prots:
	       prots[desc[i]]=0
	       high_cor_prots[desc[i]]=0
	       catcovs[desc[i]]=[]
	    prots[desc[i]]+=1
	    catcovs[desc[i]].append(cov)
	    if cov > threshold: # threshold for gr correlated protein
		high_cor_prots[desc[i]]+=1

    if cov > 0.5:# and cov < 0.8:
       hc_indices.append(row_index)
       hc_prots = hc_prots.append(row)

loc_data = pd.Series(loc_data)
prots_num = len(covs)
density = gaussian_kde(covs)
plot(xs,density(xs),label='All proteins')
title('Estimated distributions of correlation with growth rate by functional group')
print "number of proteins compared: %d" % prots_num
print "correlation by categories:"
for i in loc_prots:
    print '%d, %d' % (i,len(loc_prots[i]))

for i in loc_prots:
    density = gaussian_kde(loc_prots[i])
    plot(xs,density(xs),label=str(i))

legend(loc='upper left')
savefig('some_subcategories_distributions.png')
