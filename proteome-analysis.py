import csv
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
from scipy import stats
from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust,text,subplots
from numpy import linspace,ndarray,arange,sum,square,array,cumsum
from numpy.random import randn
from analysis import *
import matplotlib
from math import sqrt,isnan,log


### Results generation#####
def get_limits(db):
    if db == 'heinmann' and not use_LB:
        limits = (0.4,1.)
    if db == 'heinmann' and use_LB:
        limits = (0.6,1.)
    if db == 'valgepea':
        limits = (0.8,1.)
    if db == 'heinmann-chemo':
        limits = (0.5,1.)
    return limits

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

def plot_corr_hist(p,db,conc_data,categories):
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
    limits = get_limits(db)
    p.axvline(x=limits[0],ymin=0,ymax=250,ls='--',color='black',lw=0.5)
    p.axvline(x=limits[1],ymin=0,ymax=250,ls='--',color='black',lw=0.5)

    #legend(loc=2,prop={'size':8})
    tight_layout()
    return handles,labels

plot_corr_hist(p1,'valgepea',ecoli_data_v,categories)
plot_corr_hist(p2,'heinmann',ecoli_data_h,categories)

#assume both subplots have the same categories.
handles,labels=p1.get_legend_handles_labels()

figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.2,0.8,0.6,0.2),ncol=2)

text(0.03,0.8,"Valgepea et. al.",fontsize=8,transform=p.transAxes)
text(0.65,0.8,"Heinemann et. al.",fontsize=8,transform=p.transAxes)

subplots_adjust(top=0.83)
savefig('GrowthRateCorrelation.pdf')

### Global cluster analysis:
## The proteins that show a high correlation with growth rate have significant R^2 values.
## They change by xx fold across conditions measured.
## The correlation of each of the proteins with the global cluster is higher than with the GR (meaning it compensates for errors in GR measurements or degredation rates).
figure(figsize=(5,3))
def get_glob(db,df):
    limits = get_limits(db)
    glob = df[df['gr_cov']>limits[0]]
    glob = glob[glob['gr_cov']<limits[1]]
    print "for db %s global cluster is %d out of %d measured proteins" % (db, len(glob),len(df[df['gr_cov']>-1.]))
    if db == 'heinmann':
        print "for db %s annotated proteins in global are %d out of %d measured annotated proteins" % (db, len(glob[glob['group']!= "NotMapped"].index),len(df[df['group']!="NotMapped"].index))
    return glob
 
def get_high_corr(db,df,gr,conds):
    glob = get_glob(db,df)
    glob_tot = glob[conds].sum()
    alpha,beta,r_val,p_val,std_err = linregress(gr,glob_tot)
    print "global cluster sum follows alpha=%f, beta=%f" % (alpha,beta)
    print "horizontal intercept for %s is %f, corresponding to halflive %f" % (db,-beta/alpha, log(2)*alpha/beta)
    return (glob_tot,alpha,beta)

(glob_h,alpha_h,beta_h) = get_high_corr('heinmann',ecoli_data_h,gr_h,cond_list_h)
(glob_v,alpha_v,beta_v) = get_high_corr('valgepea',ecoli_data_v,gr_v,cond_list_v)

plot(gr_h.values,glob_h.values,'o',label="Heinemann et. al")
plot(gr_v.values,glob_v.values,'o',label="Valgepea et. al")
plot(gr_h.values,alpha_h*gr_h.values+beta_h,color='blue',label=("Heinemann Trend,$R^2$=%.2f" % (gr_h.corr(glob_h)**2)))
plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("Valgepea Trend,$R^2$=%.2f" % (gr_v.corr(glob_v)**2)))

xlim(xmin=0.)
ylim(ymin=0.)
xlabel('Growth rate',fontsize=10)
ylabel('Strongly correlated proteins\n fraction out of proteome',fontsize=10)
legend(loc=2, prop={'size':8},numpoints=1)
tick_params(axis='both', which='major', labelsize=8)
tick_params(axis='both', which='minor', labelsize=8)
tight_layout()
savefig('GlobalClusterGRFit%s.pdf' % conf_fname_mod)
savefig('GlobalClusterGRFit.pdf')

#gets values at cond_list and normalized in both axes
def std_err_fit(gr,s):
    alpha,beta,r,p,st = linregress(gr,s)
    n = sqrt(sum(square(s-(alpha * gr + beta)))/(len(s)-2))
    d = sqrt(sum(square(gr-gr.mean())))
    return n/d
    
def conf_int_min(degfr,s):
    res = stats.t.interval(0.95,degfr,loc=s['alpha'],scale=s['std_err'])
    return  res[0]

def conf_int_max(degfr,s):
    res = stats.t.interval(0.95,degfr,loc=s['alpha'],scale=s['std_err'])
    return  res[1]

def set_std_err(db,df,gr):
    cond_list = cond_list_dict[db]
    print "for db %s deg-free %d" %(db,len(cond_list)-2)
    df['std_err'] = df[cond_list].apply(lambda x: std_err_fit(gr[cond_list]/gr[cond_list].mean(),x/x.mean()),axis=1)
    
    df['conf_min'] = df.apply(lambda x: conf_int_min(len(cond_list)-2,x) ,axis=1)
    df['conf_max'] = df.apply(lambda x: conf_int_max(len(cond_list)-2,x) ,axis=1)
    return df

## Figure 3, global cluster slope vs. ribosomal slope
def set_alpha(db,df,gr):
    cond_list = cond_list_dict[db]
    df['alpha'] = df[cond_list].apply(lambda x: linregress(gr[cond_list]/gr[cond_list].mean(),x/x.mean())[0],axis=1)
    return df

def plot_response_hist(db,df,gr,p):
    bins = linspace(-1.7,1.7,35)
    xs = linspace(-1.75,1.75,100)
    glob_conc = get_glob(db,df)
    glob_conc = set_alpha(db,glob_conc,gr)
    glob_conc = set_std_err(db,glob_conc,gr)
    #glob_conc = glob_conc[glob_conc['std_err']<0.4]
    #print "for db %s, plotted slopes histogram includes %d proteins" % (db,len(glob_conc))
    avg = glob_conc['alpha'].mean()
    std_err = glob_conc['std_err'].mean()
    glob_conc_no_ribs = glob_conc[glob_conc['func'] != 'Ribosome']
    ribs = glob_conc[glob_conc['func'] == 'Ribosome']
    p.hist([glob_conc_no_ribs['alpha'].values,ribs['alpha'].values],bins=bins,stacked = True,label=['HC-proteins','Ribosomal proteins'])
    p.plot(xs,stats.t.pdf(xs,df=len(cond_list_dict[db])-2,loc=avg,scale=std_err)*len(glob_conc['alpha'])*0.1)
    p.set_xlim(-1.7,1.7)
    p.set_xlabel('Normalized slope',fontsize=8)
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

text(-0.01,0.9,"Valgepea et. al",fontsize=8,transform=p.transAxes)
text(0.62,0.9,"Heinemann et. al",fontsize=8,transform=p.transAxes)

tight_layout()
savefig('AllProtsVSRibosomalNormalizedSlopes.pdf')

#### plot figure of gr corr comparison by ko_num.
#hgr = []
#vgr = []
#only_in_one = 0

#v_ko_vals = set(ecoli_data_v['ko_num'].values)
#h_ko_vals = set(ecoli_data_h['ko_num'].values)
#ko_vals = v_ko_vals.union(h_ko_vals)

#for ko in ko_vals:
#    if ko == 'NotMapped':
#        continue
#    if len((ecoli_data_v[ecoli_data_v['ko_num']==ko])[['gr_cov']].values) >= 1 and len((ecoli_data_h[ecoli_data_h['ko_num']==ko])[['gr_cov']].values) >= 1:
#        vgr.append((ecoli_data_v[ecoli_data_v['ko_num']==ko])[['gr_cov']].values[0][0])
#        hgr.append((ecoli_data_h[ecoli_data_h['ko_num']==ko])[['gr_cov']].values[0][0])
#    else:
#        only_in_one +=1

#figure(figsize=(5,3))

#p=subplot(111)
#p.plot(hgr,vgr,'.')
#p.set_title('%d out of %d are only in one' % (only_in_one, len(ko_vals)))

#savefig('vhcorrcomp.pdf')

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

p2.plot(gr_chemo.values,glob_chemo.values,'o',label="Heinemann Chem.")
p2.plot(gr_chemo.values,alpha_chemo*gr_chemo.values+beta_chemo,color='blue',label=("Heinemann Chem. Trend,$R^2$=%.2f" % (gr_chemo.corr(glob_chemo)**2)))
p2.plot(gr_v.values,glob_v.values,'o',label="Valgepea")
p2.plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("Valgepea Trend,$R^2$=%.2f" % (gr_v.corr(glob_v)**2)))

p2.set_xlim(xmin=0.)
p2.set_ylim(ymin=0.)
p2.set_xlabel('Growth rate',fontsize=8)
p2.set_ylabel('Strongly correlated proteins\n fraction out of proteome',fontsize=8)
legend(loc=3, prop={'size':6},numpoints=1)
p2.tick_params(axis='both', which='major', labelsize=8)
p2.tick_params(axis='both', which='minor', labelsize=8)
tight_layout()

subplots_adjust(top=0.83)
savefig('HeinmannChemostatGr.pdf')

# plot slopes distribution for highly negatively correlated proteins from Valgepea dataset and sum of concentrations
#figure(figsize=(5,3))

#p1=subplot(121)
#p2=subplot(122)

#def get_low_corr(db,df,gr,conds):
#    if db == 'valgepea':
#        limits = (-1.0,-0.7)
#    glob = df[df['gr_cov']>limits[0]]
#    glob = glob[glob['gr_cov']<limits[1]]
#    print "for db %s anti-correlated cluster is %d out of %d measured proteins" % (db, len(glob.index),len(df.index))
#    glob_tot = glob[conds].sum()
#    alpha,beta,r_val,p_val,std_err = linregress(gr,-glob_tot)
#    return (glob_tot,alpha,beta)

#(neg_corr_v,alpha_neg,beta_neg) = get_low_corr('valgepea',ecoli_data_v,gr_v,cond_list_v)

#p2.plot(gr_v.values,neg_corr_v.values,'o',label="Valgepea anti correlated")
#p2.plot(gr_v.values,-alpha_neg*neg_corr_v.values+beta_neg,color='blue',label=("Valgepea anti correlated Trend,$R^2$=%.2f" % (gr_v.corr(neg_corr_v)**2)))
#p2.plot(gr_v.values,glob_v.values,'o',label="Valgepea")
#p2.plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("Valgepea Trend,$R^2$=%.2f" % (gr_v.corr(glob_v)**2)))

#p2.set_xlim(xmin=0.)
#p2.set_ylim(ymin=0.)
#p2.set_xlabel('Growth rate',fontsize=8)
#p2.set_ylabel('Protein fraction out of proteome',fontsize=8)
#legend(loc=3, prop={'size':6},numpoints=1)
#p2.tick_params(axis='both', which='major', labelsize=8)
#p2.tick_params(axis='both', which='minor', labelsize=8)
#tight_layout()

#subplots_adjust(top=0.83)
#savefig('Anticorrelated.pdf')

#plot rsq (or corr) vs slope.

figure(figsize=(5,3))
p1=subplot(121)
p2=subplot(122)
glob_conc = get_glob('heinmann',ecoli_data_h)
glob_conc = set_alpha('heinmann',glob_conc,gr_h)
p1.plot(glob_conc.gr_cov,glob_conc.alpha,'.')
glob_conc = get_glob('valgepea',ecoli_data_v)
glob_conc = set_alpha('valgepea',glob_conc,gr_v)
p2.plot(glob_conc.gr_cov,glob_conc.alpha,'.')
tight_layout()
savefig('slopecorr.pdf')

#plot slope vs std_err of estimate
def count(x,df):
    return float(len(df[(df['conf_min'] < x) & (df['conf_max'] > x)]))/len(df)

figure(figsize=(5,3))
p1=subplot(111)
#p2=subplot(122)
glob_conc_h = get_glob('heinmann',ecoli_data_h)
glob_conc_h = set_alpha('heinmann',glob_conc_h,gr_h)
glob_conc_h = set_std_err('heinmann',glob_conc_h,gr_h)
print "heinemann's data"
alphas = glob_conc_h['alpha'].values
mins = glob_conc_h['conf_min'].values
intervals = array(alphas)-array(mins)
alphas,intervals = (list (x) for x in zip(*sorted(zip(alphas,intervals))))
#p1.errorbar(range(0,len(alphas)),alphas,yerr=intervals,fmt='.',markersize='1',elinewidth=0.25)
p1.plot(glob_conc_h['std_err'],glob_conc_h[cond_list_dict['heinmann']].mean(axis=1),'.',markersize='1')
mins = sorted(glob_conc_h['conf_min'].values)
maxs = sorted(glob_conc_h['conf_max'].values)
allbound = sorted(mins + maxs)
fracs = [count(x,glob_conc_h) for x in allbound]
#p1.plot(allbound,fracs)
#p1.axhline(y=0.95,xmin=0,xmax=3,ls='--',color='black',lw=0.5)
p1.set_title('heinemann')
p1.set_yscale('log')
glob_conc_v = get_glob('valgepea',ecoli_data_v)
glob_conc_v = set_alpha('valgepea',glob_conc_v,gr_v)
glob_conc_v = set_std_err('valgepea',glob_conc_v,gr_v)
print "Valgepea's data"
mins = sorted(glob_conc_v['conf_min'].values)
maxs = sorted(glob_conc_v['conf_max'].values)
allbound = sorted(mins + maxs)
fracs = [count(x,glob_conc_v) for x in allbound]
#p2.plot(allbound,fracs)
#p2.set_title('valgepea')
#tight_layout()
savefig('slopestderr.pdf')

#check if for 95% of the slopes, the mean of all of the slopes lies in their 95% confidence interval

#Plot variability explained (R^2)/Var? in global cluster and in proteome as function of threshold for HC proteins.

corrs = linspace(-1,1,100)

def plotRsqRel():
    figure(figsize=(5,3))
    p1=subplot(121)
    p2=subplot(122)
    plots = {'valgepea':p1,'heinmann':p2}
    coli_data = {'valgepea':ecoli_data_v,'heinmann':ecoli_data_h}
    grs = {'valgepea':gr_v,'heinmann':gr_h}
    for db in ['valgepea','heinmann']:
        p=plots[db]
        r_sq = []
        gr = grs[db]
        conds = cond_list_dict[db]
        glob_conc = coli_data[db]
        tot_var = sqrt(glob_conc[conds].var(axis=1).sum())
        for threshold in corrs:
            glob_conc = glob_conc[glob_conc['gr_cov']>threshold]
            glob_tot = glob_conc[conds].sum()
            alpha,beta,r_val,p_val,std_err = linregress(gr,glob_tot)
            r_sq.append(gr.corr(glob_tot)**2)
        p.plot(corrs,r_sq)
        p.set_title(db)

    savefig('thresholdRsqrelation.pdf')

# calculate variability explained in proteome, take 1 (1 free parameter - selection of global cluster and scaling accordingly.
# calculate variability explained in global cluster, take 2 (1 free parameter - selection of global cluster and measurement of resulting variability reduction.
figure(figsize=(5,3))
p1=subplot(121)
p2=subplot(122)
plots = {'valgepea':p1,'heinmann':p2}
coli_data = {'valgepea':ecoli_data_v,'heinmann':ecoli_data_h}
grs = {'valgepea':gr_v,'heinmann':gr_h}
for db in ['valgepea','heinmann']:
    p=plots[db]
    gr = grs[db]
    conds = cond_list_dict[db]
    glob_conc = coli_data[db]
    glob_data = glob_conc[conds]
    tot_means = glob_data.mean(axis=1)
    tot_moved = glob_data.copy()
    for col in tot_moved.columns:
        tot_moved[col]=tot_moved[col]-tot_means
    tot_var = tot_moved **2
    tot_var = tot_var.sum()
    tot_var = tot_var.sum()
    print "tot_var is"
    print tot_var
    explained_glob = []
    explained_tot = []
    explained_normed = []
    explained_scaled = []
    for threshold in corrs:
        #print "threshold is %f" % threshold
        glob_cluster_idx = glob_conc['gr_cov']>threshold
        glob_compl_idx = glob_conc['gr_cov']<threshold
        glob_cluster_moved = tot_moved[glob_cluster_idx]
        glob_var = glob_cluster_moved ** 2
        glob_var = glob_var.sum().sum()
        glob_compl_moved = tot_moved[glob_compl_idx]
        glob_compl_var = glob_compl_moved ** 2
        glob_compl_var = glob_compl_var.sum().sum()
        #print "global cluster contains %d proteins and its variability is %f, compl variability is %f" % (len(glob_cluster_moved),glob_var,glob_compl_var)

        glob_cluster = glob_data[glob_cluster_idx]
        glob_tot = glob_cluster.sum()
        alpha,beta,r_val,p_val,std_err = linregress(gr,glob_tot)
        #print "alpha %f, beta %f" % (alpha,beta)
        vals = alpha*gr
        vals = vals+beta
        normed_vals = vals/vals.mean()
        compl = vals-1
        normed_compl = compl/compl.mean()
        scaled = glob_tot
        normed_scaled = scaled/scaled.mean()
        compl_scaled = scaled-1
        normed_compl_scaled = compl_scaled/compl_scaled.mean()

        glob_response = glob_cluster.copy()
        glob_scaled = glob_cluster.copy()
        means = tot_means[glob_cluster_idx]
        for col in glob_response.columns:
            glob_response[col]=means * normed_vals[col]
            glob_scaled[col]=means * normed_scaled[col]
        glob_cluster_response_diff = glob_cluster-glob_response
        glob_cluster_scaled_diff = glob_cluster-glob_scaled
        remaining_var = glob_cluster_response_diff **2
        remaining_var = remaining_var.sum().sum()
        remaining_scaled = glob_cluster_scaled_diff **2
        remaining_scaled = remaining_scaled.sum().sum()

        glob_compl = glob_data[glob_compl_idx]
        glob_compl_response = glob_compl.copy()
        glob_compl_scaled = glob_compl.copy()
        compl_means = tot_means[glob_compl_idx]
        for col in glob_compl_response.columns:
            glob_compl_response[col]=compl_means * normed_compl[col]
            glob_compl_scaled[col]=compl_means * normed_compl_scaled[col]
        glob_compl_response_diff = glob_compl-glob_compl_response
        glob_compl_scaled_diff = glob_compl-glob_compl_scaled
        remaining_compl_var = glob_compl_response_diff **2
        remaining_compl_var = remaining_compl_var.sum().sum()
        remaining_compl_scaled_var = glob_compl_scaled_diff **2
        remaining_compl_scaled_var = remaining_compl_scaled_var.sum().sum()

        explained_var = (glob_var-remaining_var)/glob_var
        explained_tot_frac = (glob_var-remaining_var)/tot_var
        explained_tot_compl = ((glob_var-remaining_var)+(glob_compl_var-remaining_compl_var))/tot_var
        explained_scaled_var = ((glob_var-remaining_scaled)+(glob_compl_var-remaining_compl_scaled_var))/tot_var

        explained_glob.append(explained_var)
        explained_tot.append(explained_tot_frac)
        explained_normed.append(explained_tot_compl)
        explained_scaled.append(explained_scaled_var)
        
    p.plot(corrs,explained_glob,markersize=1,label='Explained variability fraction of global cluster')
    p.plot(corrs,explained_tot,markersize=1,label='Explained variability fraction of total data')
    p.plot(corrs,explained_normed,markersize=1,label='Explained variability fraction when normalizing')
    p.plot(corrs,explained_scaled,markersize=1,label='Explained variability fraction when scaling')
    p.set_ylabel('Explained fraction of variability', fontsize=8)
    p.set_xlabel('global cluster correlation threshold', fontsize=8)
    p.tick_params(axis='both', which='major', labelsize=6)
    p.tick_params(axis='both', which='minor', labelsize=6)
    p.axhline(xmin=0,xmax=1,y=0.5,ls='--',color='black',lw=0.5)
    p.legend(loc=2,prop={'size':6})
    p.set_title(db)

tight_layout()
savefig('ExpVar2.pdf')

# compare to variability of modified proteome, where all proteins in global cluster are scaled(how) according to calculated fit, then calculating variability as above.


#Demonstrate predictive power of model by showing 3 normalized genes + trendline and how an unrelated gene "fits" the same line.

#For ron today - 6 panel graph - avg. exp. vs norm. slope, slope vs. r^2. non-global cluster avg. exp. vs. slope.
def plotMultiStats():
    figure(figsize=(5,3))
    p1=subplot(231)
    p2=subplot(232)
    p3=subplot(233)
    p4=subplot(234)
    p5=subplot(235)
    p6=subplot(236)

    glob_conc = ecoli_data_h
    gr=gr_h
    conds = cond_list_dict['heinmann']

    p1.plot(glob_conc[conds].mean(axis=1), glob_conc['rsq'],'.', markersize=1)
    p1.set_xlabel('Average concentraion', fontsize=6)
    p1.set_ylabel('$R^2$ with GR', fontsize=6)
    p1.set_xscale('log')
    p1.tick_params(axis='both', which='major', labelsize=6)
    p1.tick_params(axis='both', which='minor', labelsize=6)

    p2.plot(glob_conc[conds].mean(axis=1), glob_conc['gr_cov'],'.', markersize=1)
    p2.set_xlabel('Average concentraion', fontsize=6)
    p2.set_ylabel('Pearson corr. with GR', fontsize=6)
    p2.set_xscale('log')
    p2.tick_params(axis='both', which='major', labelsize=6)
    p2.tick_params(axis='both', which='minor', labelsize=6)

    glob_conc = glob_conc[glob_conc['gr_cov']>0.4]
    glob_conc = set_alpha('heinmann',glob_conc,gr)
    glob_conc = set_std_err('heinmann',glob_conc,gr)

    p3.plot(glob_conc[conds].mean(axis=1), glob_conc['alpha'],'.', markersize=1)
    p3.set_xlabel('Average concentraion (HC prots)', fontsize=6)
    p3.set_ylabel('Norm. Slope', fontsize=6)
    p3.set_xscale('log')
    p3.tick_params(axis='both', which='major', labelsize=6)
    p3.tick_params(axis='both', which='minor', labelsize=6)

    p4.plot(glob_conc[conds].mean(axis=1), glob_conc['std_err'],'.', markersize=1)
    p4.set_xlabel('Average concentraion (HC prots)', fontsize=6)
    p4.set_ylabel('std err of fit', fontsize=6)
    p4.set_xscale('log')
    p4.tick_params(axis='both', which='major', labelsize=6)
    p4.tick_params(axis='both', which='minor', labelsize=6)

    p5.plot(glob_conc['alpha'], glob_conc['std_err'],'.', markersize=1)
    p5.set_xlabel('Norm. slope (HC)', fontsize=6)
    p5.set_ylabel('std err of fit', fontsize=6)
    p5.tick_params(axis='both', which='major', labelsize=6)
    p5.tick_params(axis='both', which='minor', labelsize=6)

    p6.plot(glob_conc['alpha'], glob_conc['rsq'],'.', markersize=1)
    p6.set_xlabel('Norm. slope (HC)', fontsize=6)
    p6.set_ylabel('$R^2$ with GR', fontsize=6)
    p6.tick_params(axis='both', which='major', labelsize=6)
    p6.tick_params(axis='both', which='minor', labelsize=6)

    tight_layout()

    savefig('AvgConcStatsHein.pdf')

#comulative graph - x axis - avg. prot. conc. (or molecule count per cell), y axis, comulative % out of proteome.
def plotComulativeGraph():
    figure(figsize=(5,3))
    p1 = subplot(121)
    p2 = subplot(122)
    avgs = sorted(ecoli_data_h[conds].mean(axis=1).values)

    p1.plot(avgs,cumsum(avgs),'.',markersize=0.5)
    p1.set_xscale('log')
    p1.set_xlabel('Avg. prot. conc.',fontsize=6)
    p1.set_ylabel('accumulated fraction \n out of proteome',fontsize=6)
    p1.axhline(xmin=0,xmax=1,y=0.05,ls='--',color='black',lw=0.5)
    p1.axhline(xmin=0,xmax=1,y=0.01,ls='--',color='black',lw=0.5)
    p1.tick_params(axis='both', which='major', labelsize=6)
    p1.tick_params(axis='both', which='minor', labelsize=6)

    p2.plot(arange(0,len(avgs)),cumsum(avgs),'.',markersize=0.5)
    p2.set_xlabel('num. of prots',fontsize=6)
    p2.set_ylabel('accumulated fraction \n out of proteome',fontsize=6)
    p2.axhline(xmin=0,xmax=2000,y=0.05,ls='--',color='black',lw=0.5)
    p2.axhline(xmin=0,xmax=2000,y=0.01,ls='--',color='black',lw=0.5)
    p2.tick_params(axis='both', which='major', labelsize=6)
    p2.tick_params(axis='both', which='minor', labelsize=6)

    tight_layout()
    savefig('DistStatsHein.pdf')

#refactor all graphs to use for db in ['val','hein'].
