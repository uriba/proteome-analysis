import csv
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
#from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot
from numpy import linspace,ndarray,arange
from numpy.random import randn


#Initialization of basic data containers, gene annotation data, growth rates and cell volumes and selection of conditions to analyze.
def uniprot_to_desc_dict():
    uni_konum_dict = {}
    uni_to_konum = read_csv('eco_uniprot_mapping.csv',sep='[\t:]',encoding='iso-8859-1',header = None, names = ['ko','bla','uniprot'])
    for i,row in uni_to_konum.iterrows():
         uni_konum_dict[row['uniprot']]=row['ko']    

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
                ko_annot_dict[row[-1]]=(cat,subcat,component)
    uni_to_annot = {}
    for uni in uni_konum_dict:
        if uni_konum_dict[uni] == 'NotMapped':
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
    return uniprot_to_location

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
gr = gr[cond_list]

# define cell volumes:
volumes = {
    u'chemostat \u00b5=0.12': 2.1,
    u'galactose':1.9,
    u'chemostat \u00b5=0.20':2.2,
    u'acetate':2.4,
    u'chemostat \u00b5=0.35':2.4,
    u'glucosamine':2.9,
    u'pyruvate':2.1,
    u'glycerol':2.3,
    u'fumarate':2.4,
    u'succinate':2.4,
    u'chemostat \u00b5=0.5':2.6,
    u'anaerobic' : 2.9,
    u'glucose': 3.2,
    u'LB':4.4
}
volumes = pd.Series(volumes)
volumes = volumes[cond_list]

#Define the text fields that will be relevant for the analysis:
desc_list = ['Description','UP_AC']

#Convert dataframe types to standart types for analysis
def convert_types(df):
    df = df[df != 'below LOQ']
    df[cond_list] = df[cond_list].astype('float')
    df[desc_list] = df[desc_list].astype('string')
    return df

def get_coli_data(use_weight):
    # As the file was exported from Excel, it uses Excel's encoding.
    ecoli_data = read_csv('coli_data.csv',header=1,encoding='iso-8859-1')

    #Split the data loaded into two sets - weight and count.
    idx = ecoli_data.columns
    ecoli_data_count = ecoli_data[idx[0:29]]
    ecoli_data_weight = ecoli_data[idx[0:10].append(idx[29:])]

    # Refine the DataFrames to include only these conditions (and the protein descriptions):
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
    ecoli_data_weight = convert_types(ecoli_data_weight)
    ecoli_data = convert_types(ecoli_data)

    #Normalize to get concentrations (Use the total protein weight as the normalizing factor to avoid errors in cell volume measurements)
    ecoli_data[cond_list] = ecoli_data[cond_list] / ecoli_data_weight[cond_list].sum() #volumes[cond_list]
    return ecoli_data

### Results generation#####
### Figure 1 - Correlation to growth rate by functional group histogram.
uni_to_annot = uniprot_to_desc_dict()
ecoli_data = get_coli_data(use_weight=True)

conc_data = ecoli_data
conc_data = conc_data.dropna()
conc_data['gr_cov']=conc_data[cond_list].apply(lambda x: x.corr(gr[cond_list]),axis=1)
conc_data['rsq']=conc_data['gr_cov']**2
conc_data['group']=conc_data.apply(lambda x: (uni_to_annot[x['UP_AC']])[0],axis=1)

categories = set(conc_data['group'].values)
# Remove the unmapped proteins first and add them at the end so that they are stacked last in the histogram.
categories.remove("NotMapped")
bins = linspace(-1,1,20)
covs = ndarray(shape=(len(categories),len(bins)-1))
sets = [] 
figure(figsize=(5,3))
for x in categories:
    sets.append(conc_data[conc_data['group']==x].gr_cov)
    #sets.append(hist(conc_data[conc_data['group']==x].gr_cov,bins)[0])
sets.append(conc_data[conc_data['group']=="NotMapped"].gr_cov)
cats = list(categories)
cats.append("NotMapped")

hist(sets,bins = bins, stacked = True,label=cats)

tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
xlabel('Pearson correlation with growth rate',fontsize=10)
ylabel('Number of proteins',fontsize=10)

legend(loc=2,prop={'size':8})
tight_layout()
savefig('GrowthRateCorrelation.pdf')


### Global cluster analysis:
## The proteins that show a high correlation with growth rate have significant R^2 values.
## They change by xx fold across conditions measured.
## The correlation of each of the proteins with the global cluster is higher than with the GR (meaning it compensates for errors in GR measurements or degredation rates).
figure(figsize=(5,3))

high_corr_prots = conc_data[conc_data['gr_cov']>0.4]
high_corr_prots = high_corr_prots[high_corr_prots['gr_cov']<0.8]
high_corr_normed = high_corr_prots.copy()
high_corr_normed = high_corr_normed[cond_list].apply(lambda x: x/x.mean(),axis=1)

def cluster_corr(cluster):
    cluster_global = cluster.sum()
    cluster_global = cluster_global/cluster_global.mean()

    alpha,beta,r_val,p_val,std_err = linregress(gr,cluster_global)
    return (cluster_global,alpha,beta)

global_cluster = {}
global_weighted=cluster_corr(high_corr_prots[cond_list])
global_normed=cluster_corr(high_corr_normed[cond_list])

plot(gr.values,global_weighted[0].values,'o',label="Weighted")
plot(gr.values,global_weighted[1]*gr.values+global_weighted[2],color='blue',label=("Weighted Trend,$R^2$=%f" % (gr.corr(global_weighted[0])**2)))

plot(gr.values,global_normed[0].values,'o',label="Normalized")
plot(gr.values,global_normed[1]*gr.values+global_normed[2],color='green',label=("Normalized Trend,$R^2$=%f" % (gr.corr(global_normed[0])**2)))

xlim(0,0.7)
ylim(0,2)
xlabel('Growth rate',fontsize=10)
ylabel('Protein level (normalized)',fontsize=10)
legend(loc=2, prop={'size':8})
tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
tight_layout()
savefig('GlobalClusterGRFit.pdf')

## Figure 2, correlation inside global cluster
figure(figsize=(5,3))

high_corr_prots['weighted_cov']=high_corr_prots[cond_list].apply(lambda x: x.corr(global_weighted[0]),axis=1)
high_corr_prots['normed_cov']=high_corr_prots[cond_list].apply(lambda x: x.corr(global_normed[0]),axis=1)
sets = [high_corr_prots['weighted_cov'].values,high_corr_prots['normed_cov'].values]
hist(sets,bins = bins, stacked = False,label=['Weighted','Normalized'])
legend(loc=2, prop={'size':8})
xlabel('Pearson correlation with global cluster',fontsize=10)
ylabel('Number of proteins',fontsize=10)
tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
tight_layout()
savefig('GlobalClusterCorr.pdf')

## Figure 3, R^2 of proteins with global cluster
figure(figsize=(5,3))
sets = [(high_corr_prots['weighted_cov']**2).values,(high_corr_prots['normed_cov']**2).values]
hist(sets, stacked = False,label=['Weighted','Normalized'],bins=20)
legend(loc=2, prop={'size':8})
xlabel('R-square of protein with global cluster',fontsize=10)
ylabel('Number of proteins',fontsize=10)
tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
tight_layout()
savefig('GlobalClusterRSquare.pdf')

## Figure 4, coherent scaling of proteins in the global cluster - R^2 comparison between global cluster and specific fits.
def rsq(ys,xs,alpha,beta):
    n = len(xs)
    return 1.0-((ys-(alpha*xs + beta))**2).sum()/((n-1)*ys.var())

def rsq_self(ys,xs):
    alpha,beta,r_val,p_val,std_err = linregress(xs,ys)
    return rsq(ys,xs,alpha,beta)

rsq_global = high_corr_normed[cond_list].apply(lambda x: rsq(x,global_normed[0],1,0),axis=1)
rsq_selfs = high_corr_normed[cond_list].apply(lambda x: rsq_self(x,global_normed[0]),axis=1)
#alphas = high_corr_normed[cond_list].apply(lambda x: linregress(gr/gr.mean(),x/x.mean())[0],axis=1)

high_conc = conc_data[cond_list] #(conc_data[cond_list])[conc_data[cond_list]>conc_data[cond_list].median()].dropna()
#high_conc = (high_corr_prots[cond_list])[high_corr_prots[cond_list]>conc_data[cond_list].median()].dropna()
#high_conc = (conc_data[cond_list])[conc_data[cond_list]>conc_data[cond_list].median()].dropna()

high_conc['alpha'] = high_conc[cond_list].apply(lambda x: linregress(gr/gr.mean(),x/x.mean())[0],axis=1)
high_conc['inverse'] = high_conc[cond_list].apply(lambda x: linregress(x/x.mean(),gr/gr.mean())[0],axis=1)
high_conc['rsq'] = high_conc[cond_list].apply(lambda x: linregress(gr/gr.mean(),x/x.mean())[2]**2,axis=1)
linear_cluster = high_conc[high_conc['alpha']>0.9]
linear_cluster = linear_cluster[linear_cluster['alpha']<1.1]
#alphas = conc_data[cond_list].apply(lambda x: linregress((x/x.mean()),gr/gr.mean())[0],axis=1)

figure(figsize=(7,6))
p1 = subplot(221)
linsum = linear_cluster[cond_list].sum()
p1.plot(gr,linsum,'.')
#p1.set_xlim(0,0.7)
#p1.set_ylim(0,0.02)
p2=subplot(222)
p2.plot(high_conc[cond_list].mean(axis=1),high_conc['alpha'],'.',markersize=1)
p2.axhline(y=0,xmin=0,xmax=1)
p2.axhline(y=0.5,xmin=0,xmax=1)
p2.axhline(y=1,xmin=0,xmax=1)
#p2.axvline(x=0,ymin=0,ymax=100,color='r')
p2.set_xscale('log')
#p2.set_ylim(0,50)
#plot(rsq_global,rsq_selfs,'.',markersize=1)
#hist(rsq_selfs-rsq_global,30)
#xlim(-4,1)
m = high_conc['rsq'].median()
p3=subplot(223)
p3.plot(high_conc[cond_list].mean(axis=1),high_conc['rsq'],'.',markersize=1)
p3.set_xlabel('mean')
p3.set_ylabel('rsq')
p3.set_xscale('log')
p3.axhline(y=m,xmin=0,xmax=1)

p4=subplot(224)
p4.plot(high_conc['rsq'],high_conc['alpha'],'.',markersize=1)
p4.set_xlabel('rsq')
p4.set_ylabel('alpha')
#xlabel('$R^2$ of global cluster as predictor for normalized protein',fontsize=10)
#ylabel('$R^2$ of normalized protein with global cluster',fontsize=10)
#xlabel('$\Delta R^2$ between best fit and fit to global cluster',fontsize=10)
tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
tight_layout()
savefig('globalSlope.pdf')

## other figures: abundance vs, correlation, gene location vs. correlation, amount at two adjacent conditions with trendlines at equal amounts and global scaling amounts.

####################################### reference figure or R^2 distribution of random series
figure(figsize=(5,3))
xs = pd.DataFrame(randn(1000,10),columns=arange(0,10),index=arange(0,1000))
ys = pd.Series(arange(0,10))
cors = xs.apply(lambda x: x.corr(ys),axis=1)**2
hist(cors,bins=20)

tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
xlabel('R-square of random series with arbitrary series',fontsize=10)
ylabel('Number of proteins',fontsize=10)
tight_layout()
savefig('RandomRSquare.pdf')
#############################################################################################
