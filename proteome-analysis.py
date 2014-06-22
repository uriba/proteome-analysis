import csv
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust,text,subplots
from numpy import linspace,ndarray,arange
from numpy.random import randn
from analysis import *
import matplotlib

### Results generation#####
### Figure 1 - Correlation to growth rate by functional group histogram.
(cond_list_v,gr_v,ecoli_data_v) = get_annotated_prots('valgepea')
(cond_list_h,gr_h,ecoli_data_h) = get_annotated_prots('heinmann')
ecoli_data_h = calc_gr_corr(ecoli_data_h,cond_list_h,gr_h)
ecoli_data_v = calc_gr_corr(ecoli_data_v,cond_list_v,gr_v)

categories = set(ecoli_data_v['group'].values).union(set(ecoli_data_h['group'].values))

# Remove the unmapped proteins first and add them at the end so that they are stacked last in the histogram.
if not remove_unmapped and "NotMapped" in categories:
    categories.remove("NotMapped")
categories = list(categories)
categories.sort()
categories = ['Metabolism','Genetic Information Processing','Environmental Information Processing', 'Cellular Processes']

if not just_ribosomes:
    categories.append('NotMapped')

figure(figsize=(5,3))

p=subplot(111)
p1=subplot(121)
p2=subplot(122)

def plot_corr_hist(p,conc_data,categories):
    bins = linspace(-1,1,21)
    covs = ndarray(shape=(len(categories),len(bins)-1))
    sets = [] 

    for x in categories:
        sets.append(conc_data[conc_data['group']==x].gr_cov)

    p.hist(sets,bins = bins, stacked = True,label=categories)
    handles,labels=p.get_legend_handles_labels()
    p.tick_params(axis='both', which='major', labelsize=8)
    p.tick_params(axis='both', which='minor', labelsize=8)
    p.set_xlabel('Pearson correlation with growth rate',fontsize=8)
    p.set_ylabel('Number of proteins',fontsize=8)

    #legend(loc=2,prop={'size':8})
    tight_layout()
    return handles,labels

plot_corr_hist(p1,ecoli_data_v,categories)
plot_corr_hist(p2,ecoli_data_h,categories)

#assume both subplots have the same categories.
handles,labels=p1.get_legend_handles_labels()

figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.2,0.8,0.6,0.2),ncol=2)

text(0.03,0.8,"Valgepea",fontsize=8,transform=p.transAxes)
text(0.65,0.8,"Heinemann",fontsize=8,transform=p.transAxes)

subplots_adjust(top=0.83)
savefig('GrowthRateCorrelation.pdf')

### Global cluster analysis:
## The proteins that show a high correlation with growth rate have significant R^2 values.
## They change by xx fold across conditions measured.
## The correlation of each of the proteins with the global cluster is higher than with the GR (meaning it compensates for errors in GR measurements or degredation rates).
figure(figsize=(5,3))

def get_glob(db,df):
    if db == 'heinmann' and not use_LB:
        limits = (0.4,0.8)
    if db == 'heinmann' and use_LB:
        limits = (0.6,1.)
    if db == 'valgepea':
        limits = (0.8,1.)
    if db == 'heinmann-chemo':
        limits = (0.8,1.)
    glob = df[df['gr_cov']>limits[0]]
    glob = glob[glob['gr_cov']<limits[1]]
    print "for db %s global cluster is %d out of %d measured proteins" % (db, len(glob.index),len(df.index))
    if db == 'heinmann':
        print "for db %s annotated proteins in global are %d out of %d measured annotated proteins" % (db, len(glob[glob['group']!= "NotMapped"].index),len(df[df['group']!="NotMapped"].index))
    return glob
 
def get_high_corr(db,df,gr,conds):
    glob = get_glob(db,df)
    glob_tot = glob[conds].sum()
    alpha,beta,r_val,p_val,std_err = linregress(gr,glob_tot)
    return (glob_tot,alpha,beta)

(glob_h,alpha_h,beta_h) = get_high_corr('heinmann',ecoli_data_h,gr_h,cond_list_h)
(glob_v,alpha_v,beta_v) = get_high_corr('valgepea',ecoli_data_v,gr_v,cond_list_v)

plot(gr_h.values,glob_h.values,'o',label="Heinemann")
plot(gr_v.values,glob_v.values,'o',label="Valgepea")
plot(gr_h.values,alpha_h*gr_h.values+beta_h,color='blue',label=("Heinemann Trend,$R^2$=%.2f" % (gr_h.corr(glob_h)**2)))
plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("Valgepea Trend,$R^2$=%.2f" % (gr_v.corr(glob_v)**2)))

xlim(xmin=0.)
ylim(ymin=0.)
xlabel('Growth rate',fontsize=10)
ylabel('Protein fraction out of proteome',fontsize=10)
legend(loc=2, prop={'size':8},numpoints=1)
tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
tight_layout()
savefig('GlobalClusterGRFit%s.pdf' % conf_fname_mod)
savefig('GlobalClusterGRFit.pdf')

## Figure 3, global cluster slope vs. ribosomal slope
def set_alpha(db,df,gr):
    cond_list = cond_list_dict[db]
    df['alpha'] = df[cond_list].apply(lambda x: linregress(gr[cond_list]/gr[cond_list].mean(),x/x.mean())[0],axis=1)
    return df

def plot_response_hist(db,df,gr,p):
    bins = linspace(-1.7,1.7,35)
    ribs = df[df['func'] == 'Ribosome']
    glob_conc = get_glob(db,df)
    glob_conc = glob_conc[glob_conc['func'] != 'Ribosome']
    glob_conc = set_alpha(db,glob_conc,gr)
    ribs = set_alpha(db,ribs,gr)
    p.hist([glob_conc['alpha'].values,ribs['alpha'].values],bins=bins,stacked = True,label=['HC-proteins','Ribosomal proteins'])
    p.set_xlim(-1.7,1.7)
    p.set_xlabel('Normalized response',fontsize=8)
    p.set_ylabel('Number of proteins',fontsize=8)
    p.axvline(x=0,ymin=0,ymax=100,ls='--',color='black',lw=0.5)
    p.axvline(x=0.5,ymin=0,ymax=100,ls='--',color='black',lw=0.5)
    p.axvline(x=1,ymin=0,ymax=100,ls='--',color='black',lw=0.5)
    p.tick_params(axis='both', which='major', labelsize=8)
    p.tick_params(axis='both', which='minor', labelsize=8)

figure(figsize=(5,3))

p=subplot(111)
p1=subplot(121)
p2=subplot(122)
plot_response_hist('valgepea',ecoli_data_v,gr_v,p1)
plot_response_hist('heinmann',ecoli_data_h,gr_h,p2)

text(-0.01,0.9,"Valgepea",fontsize=8,transform=p.transAxes)
text(0.62,0.9,"Heinemann",fontsize=8,transform=p.transAxes)

tight_layout()
savefig('AllProtsVSRibosomalNormalizedSlopes.pdf')

#### plot figure of gr corr comparison by ko_num.
hgr = []
vgr = []
only_in_one = 0

v_ko_vals = set(ecoli_data_v['ko_num'].values)
h_ko_vals = set(ecoli_data_h['ko_num'].values)
ko_vals = v_ko_vals.union(h_ko_vals)

for ko in ko_vals:
    if ko == 'NotMapped':
        continue
    if len((ecoli_data_v[ecoli_data_v['ko_num']==ko])[['gr_cov']].values) >= 1 and len((ecoli_data_h[ecoli_data_h['ko_num']==ko])[['gr_cov']].values) >= 1:
        vgr.append((ecoli_data_v[ecoli_data_v['ko_num']==ko])[['gr_cov']].values[0][0])
        hgr.append((ecoli_data_h[ecoli_data_h['ko_num']==ko])[['gr_cov']].values[0][0])
    else:
        only_in_one +=1

figure(figsize=(5,3))

p=subplot(111)
p.plot(hgr,vgr,'.')
p.set_title('%d out of %d are only in one' % (only_in_one, len(ko_vals)))

savefig('vhcorrcomp.pdf')

# plot heinmann data only for chemostat conditions.
figure(figsize=(5,3))
(cond_list,gr_chemo,ecoli_data_chemo) = get_annotated_prots('heinmann-chemo')
ecoli_data_chemo = calc_gr_corr(ecoli_data_chemo,cond_list,gr_chemo)

categories = set(ecoli_data_chemo['group'].values)

# Remove the unmapped proteins first and add them at the end so that they are stacked last in the histogram.
if not remove_unmapped and "NotMapped" in categories:
    categories.remove("NotMapped")
categories = list(categories)
categories.sort()
categories = ['Metabolism','Genetic Information Processing','Environmental Information Processing', 'Cellular Processes']

if not just_ribosomes:
    categories.append('NotMapped')

p1=subplot(121)
p2=subplot(122)

def plot_corr_hist(p,conc_data,categories):
    bins = linspace(-1,1,21)
    covs = ndarray(shape=(len(categories),len(bins)-1))
    sets = [] 

    for x in categories:
        sets.append(conc_data[conc_data['group']==x].gr_cov)

    p.hist(sets,bins = bins, stacked = True,label=categories)
    handles,labels=p.get_legend_handles_labels()
    p.tick_params(axis='both', which='major', labelsize=8)
    p.tick_params(axis='both', which='minor', labelsize=8)
    p.set_xlabel('Pearson correlation with growth rate',fontsize=8)
    p.set_ylabel('Number of proteins',fontsize=8)

    #legend(loc=2,prop={'size':8})
    tight_layout()
    return handles,labels

plot_corr_hist(p1,ecoli_data_chemo,categories)

handles,labels=p1.get_legend_handles_labels()

figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.05,0.8,0.5,0.2),ncol=2)

(glob_chemo,alpha_chemo,beta_chemo) = get_high_corr('heinmann-chemo',ecoli_data_chemo,gr_chemo,cond_list)

p2.plot(gr_chemo.values,glob_chemo.values,'o',label="Heinemann Chemostat")
p2.plot(gr_chemo.values,alpha_chemo*gr_chemo.values+beta_chemo,color='blue',label=("Heinemann Chemostat Trend,$R^2$=%.2f" % (gr_chemo.corr(glob_chemo)**2)))
p2.plot(gr_v.values,glob_v.values,'o',label="Valgepea")
p2.plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("Valgepea Trend,$R^2$=%.2f" % (gr_v.corr(glob_v)**2)))

p2.set_xlim(xmin=0.)
p2.set_ylim(ymin=0.)
p2.set_xlabel('Growth rate',fontsize=8)
p2.set_ylabel('Protein fraction out of proteome',fontsize=8)
legend(loc=3, prop={'size':6},numpoints=1)
p2.tick_params(axis='both', which='major', labelsize=8)
p2.tick_params(axis='both', which='minor', labelsize=8)
tight_layout()

subplots_adjust(top=0.83)
savefig('HeinmannChemostatGr.pdf')
