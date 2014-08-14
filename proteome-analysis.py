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
import random


### Results generation#####
def get_limits(db):
    if db == 'Heinemann' and not use_LB:
        limits = (0.4,1.)
    if db == 'Heinemann' and use_LB:
        limits = (0.6,1.)
    if db == 'Valgepea':
        limits = (0.8,1.)
    if db == 'Heinemann-chemo':
        limits = (0.5,1.)
    return limits

### Figure 1 - Correlation to growth rate by functional group histogram.
(cond_list_v,gr_v,ecoli_data_v) = get_annotated_prots('Valgepea')
(cond_list_h,gr_h,ecoli_data_h) = get_annotated_prots('Heinemann')
categories = set(ecoli_data_v['group'].values).union(set(ecoli_data_h['group'].values))

dbs = ['Heinemann','Valgepea']
cond_lists = {'Heinemann':cond_list_h,'Valgepea':cond_list_v}
grs = {'Heinemann':gr_h,'Valgepea':gr_v}
coli_datas = {'Heinemann':ecoli_data_h,'Valgepea':ecoli_data_v}

for db in dbs:
    coli_datas[db]= calc_gr_corr(coli_datas[db],cond_lists[db],grs[db])

def writeCorrsHist(db):
    conc_data = coli_datas[db]
    conds = cond_lists[db]
    limits = get_limits(db)
    threshold = limits[0]
    funcs = conc_data['func'].unique()
    func_stat = []
    for func in funcs:
        conc_func = conc_data[conc_data['func']==func]
        corred_idx = conc_func['gr_cov']>threshold
        tot = len(conc_func)
        tot_means = conc_func[conds].mean(axis=1)
        corr_means = tot_means[corred_idx]
        correlated = len(corr_means)
        func_stat.append(("{%s}" % func,tot,tot_means.sum()*100,correlated,corr_means.sum()*100))
    with open('funcs%s.csv' % db,'wb') as csvfile:
        csvwriter = csv.writer(csvfile,delimiter=';')
        csvwriter.writerow(['Function','Number of proteins','totPrctP','Correlated proteins','corPrctP'])
        for row in func_stat:
            csvwriter.writerow(row)

def writeTables():
    for db in dbs:
        writeCorrsHist(db)

categories = ['Metabolism','Genetic Information Processing','Environmental Information Processing', 'Cellular Processes','NotMapped']

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


def plotCorrelationHistograms():
    figure(figsize=(5,3))

    p=subplot(111)
    ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
    coords = {'Heinemann':0.03,'Valgepea':0.65}

    for db in dbs:
        plot_corr_hist(ps[db],db,coli_datas[db],categories)
        text(coords[db],0.8,"%s et. al." % db,fontsize=8,transform=p.transAxes)

    #assume both subplots have the same categories.
    handles,labels=ps['Heinemann'].get_legend_handles_labels()

    figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.2,0.8,0.6,0.2),ncol=2)

    subplots_adjust(top=0.83)
    savefig('GrowthRateCorrelation.pdf')
    savefig('GrowthRateCorrelation.png')

### Figure 3, Global cluster analysis:
def plotGlobalResponse():
    figure(figsize=(5,3))

    colors = {'Heinemann':'blue','Valgepea':'green'}

    for db in dbs:
        conds = cond_lists[db]
        coli_data = coli_datas[db]
        glob = get_glob(db,coli_data)
        gr = grs[db]
        #gr = gr[conds]

        print "%s global cluster is %d out of %d measured proteins" % (db, len(glob),len(coli_data[coli_data['gr_cov']>-1.]))

        glob_tot = glob[conds].sum()
        alpha,beta,r_val,p_val,std_err = linregress(gr,glob_tot)

        print "global cluster sum follows alpha=%f, beta=%f" % (alpha,beta)
        print "horizontal intercept for %s is %f, corresponding to halflive %f" % (db,-beta/alpha, log(2)*alpha/beta)
        plot(gr.values,glob_tot.values,'o',label="%s et. al" % db,color=colors[db])
        plot(gr.values,alpha*gr.values+beta,color=colors[db],label=("%s Trend,$R^2$=%.2f" % (db,gr.corr(glob_tot)**2)))

    xlim(xmin=0.)
    ylim(ymin=0.)
    xlabel('Growth rate',fontsize=10)
    ylabel('Strongly correlated proteins\n fraction out of proteome',fontsize=10)
    legend(loc=2, prop={'size':8},numpoints=1)
    tick_params(axis='both', which='major', labelsize=8)
    tick_params(axis='both', which='minor', labelsize=8)
    tight_layout()
    savefig('GlobalClusterGRFit.pdf')
    savefig('GlobalClusterGRFit.png')

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

## Figure 2, global cluster slope vs. ribosomal slope
def get_glob(db,df):
    limits = get_limits(db)
    return df[df['gr_cov']>limits[0]]
 
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
    glob_conc_no_ribs = glob_conc[glob_conc['prot'] != 'Ribosome']
    ribs = glob_conc[glob_conc['prot'] == 'Ribosome']
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
ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
coords = {'Heinemann':-0.01,'Valgepea':0.62}
for db in dbs:
    plot_response_hist(db,coli_datas[db],grs[db],ps[db])
    text(coords[db],0.9,"%s et. al" % db,fontsize=8,transform=p.transAxes)

tight_layout()
savefig('AllProtsVSRibosomalNormalizedSlopes.pdf')
savefig('AllProtsVSRibosomalNormalizedSlopes.png')

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

# plot Heinemann data only for chemostat conditions.
db = 'Heinemann-chemo'
figure(figsize=(5,3))
(cond_list,gr_chemo,ecoli_data_chemo) = get_annotated_prots(db)
ecoli_data_chemo = calc_gr_corr(ecoli_data_chemo,cond_list,gr_chemo)

categories = ['Metabolism','Genetic Information Processing','Environmental Information Processing', 'Cellular Processes','NotMapped']

p1=subplot(121)
p2=subplot(122)

plot_corr_hist(p1,db,ecoli_data_chemo,categories)

handles,labels=p1.get_legend_handles_labels()

figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.05,0.8,0.5,0.2),ncol=2)

glob_chemo = get_glob(db,ecoli_data_chemo)
print "%s global cluster is %d out of %d measured proteins" % (db, len(glob_chemo),len(ecoli_data_chemo[ecoli_data_chemo['gr_cov']>-1.]))

glob_tot_chemo = glob_chemo[cond_list].sum()
alpha,beta,r_val,p_val,std_err = linregress(gr_chemo,glob_tot_chemo)

print "global cluster sum follows alpha=%f, beta=%f" % (alpha,beta)
print "horizontal intercept for %s is %f, corresponding to halflive %f" % (db,-beta/alpha, log(2)*alpha/beta)
p2.plot(gr_chemo.values,glob_tot_chemo.values,'o',label="%s et. al" % db,color='blue')
p2.plot(gr_chemo.values,alpha*gr_chemo.values+beta,color='blue',label=("Heinemann Chem. Trend,$R^2$=%.2f" % (gr_chemo.corr(glob_tot_chemo)**2)))
#(glob_v,alpha_v,beta_v) = get_high_corr('Valgepea',coli_datas['Valgepea'],grs['Valgepea'],cond_lists['Valgepea'])
#p2.plot(gr_v.values,glob_v.values,'o',label="Valgepea")
#p2.plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("Valgepea Trend,$R^2$=%.2f" % (gr_v.corr(glob_v)**2)))

p2.set_xlim(xmin=0.)
p2.set_ylim(ymin=0.)
p2.set_xlabel('Growth rate',fontsize=8)
p2.set_ylabel('Strongly correlated proteins\n fraction out of proteome',fontsize=8)
legend(loc=3, prop={'size':6},numpoints=1)
p2.tick_params(axis='both', which='major', labelsize=8)
p2.tick_params(axis='both', which='minor', labelsize=8)
tight_layout()

subplots_adjust(top=0.83)
savefig('HeinemannChemostatGr.pdf')
savefig('HeinemannChemostatGr.png')

# plot slopes distribution for highly negatively correlated proteins from Valgepea dataset and sum of concentrations
#figure(figsize=(5,3))

#p1=subplot(121)
#p2=subplot(122)

#def get_low_corr(db,df,gr,conds):
#    if db == 'Valgepea':
#        limits = (-1.0,-0.7)
#    glob = df[df['gr_cov']>limits[0]]
#    glob = glob[glob['gr_cov']<limits[1]]
#    print "for db %s anti-correlated cluster is %d out of %d measured proteins" % (db, len(glob.index),len(df.index))
#    glob_tot = glob[conds].sum()
#    alpha,beta,r_val,p_val,std_err = linregress(gr,-glob_tot)
#    return (glob_tot,alpha,beta)

#(neg_corr_v,alpha_neg,beta_neg) = get_low_corr('Valgepea',ecoli_data_v,gr_v,cond_list_v)

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

#check if for 95% of the slopes, the mean of all of the slopes lies in their 95% confidence interval

#Plot variability explained (R^2)/Var? in global cluster and in proteome as function of threshold for HC proteins.

corrs = linspace(-1,1,100)

# calculate variability explained in proteome, take 1 (1 free parameter - selection of global cluster and scaling accordingly.
# calculate variability explained in global cluster, take 2 (1 free parameter - selection of global cluster and measurement of resulting variability reduction.
def calc_var(df):
    df = df.copy()
    means = df.mean(axis=1)
    for col in df.columns:
        df[col]=df[col]-means
    var = df**2
    var = var.sum().sum()
    return var

def calc_explained_var(df,gr):
    tot_var = calc_var(df)
    df = df.copy()
    means = df.mean(axis=1)
    response = df.sum()
    scaled_response = response/response.mean()
    alpha,beta,r_val,p_val,std_err = linregress(gr,response)
    normed_response = alpha*gr+beta
    normed_response = normed_response/normed_response.mean()
    pred = df.copy()
    scaled = df.copy()
    for col in pred.columns:
        pred[col]=means*normed_response[col]
        scaled[col]=means*scaled_response[col]
    remains = df-pred
    remained_var = calc_var(remains)
    scaled_remains = df-scaled
    remained_scaled_var = calc_var(scaled_remains)
    alpha,beta,r_val,p_val,std_err = linregress(gr,response/response.mean())
    return (alpha,tot_var,tot_var - remained_var,tot_var-remained_scaled_var)


def variabilityAndGlobClustSlopes():
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
    coli_data = {'Valgepea':ecoli_data_v,'Heinemann':ecoli_data_h}
    grs = {'Valgepea':gr_v,'Heinemann':gr_h}
    alphas = {'Valgepea':[],'Heinemann':[]}
    for db in ['Valgepea','Heinemann']:
        p=ps[db]
        conds = cond_list_dict[db]
        gr = grs[db]
        gr = gr[conds]
        glob_conc = coli_data[db]
        glob_data = glob_conc[conds]
        tot_var = calc_var(glob_data)
        print "tot_var is %f" % tot_var
        explained_glob = []
        explained_tot = []
        explained_compl_glob = []
        explained_compl_tot = []
        explained_normed = []
        explained_scaled = []
        for threshold in corrs:
            glob_cluster_idx = glob_conc['gr_cov']>threshold
            glob_compl_idx = glob_conc['gr_cov']<threshold
            alpha,glob_var,glob_explained,glob_scaled_explained = calc_explained_var(glob_data[glob_cluster_idx],gr)
            a,compl_var,compl_explained,compl_scaled_explained = calc_explained_var(glob_data[glob_compl_idx],gr)
            alphas[db].append(alpha)
            explained_var = glob_explained/glob_var
            explained_compl_var = compl_explained/compl_var
            explained_tot_frac = glob_explained/tot_var
            explained_compl_tot_frac = compl_explained/tot_var
            explained_tot_compl = (glob_explained+compl_explained)/tot_var
            explained_scaled_var = (glob_scaled_explained+compl_scaled_explained)/tot_var

            explained_glob.append(explained_var)
            explained_compl_glob.append(explained_compl_var)
            explained_tot.append(explained_tot_frac)
            explained_compl_tot.append(explained_compl_tot_frac)
            explained_normed.append(explained_tot_compl)
            explained_scaled.append(explained_scaled_var)

        p.plot(corrs,explained_glob,markersize=1,label='Explained variability fraction of global cluster')
        p.plot(corrs,explained_tot,markersize=1,label='Explained variability fraction of total data')
        p.plot(corrs,explained_compl_glob,markersize=1,label='Explained complementary variability fraction of global cluster')
        p.plot(corrs,explained_compl_tot,markersize=1,label='Explained complementary variability fraction of total data')
        p.plot(corrs,explained_normed,markersize=1,label='Explained variability fraction when normalizing')
        p.plot(corrs,explained_scaled,markersize=1,label='Explained variability fraction when scaling')
        p.set_ylabel('Explained fraction of variability', fontsize=8)
        p.set_xlabel('global cluster correlation threshold', fontsize=8)
        p.set_ylim(0,1)
        p.tick_params(axis='both', which='major', labelsize=6)
        p.tick_params(axis='both', which='minor', labelsize=6)
        p.axhline(xmin=0,xmax=1,y=0.5,ls='--',color='black',lw=0.5)
        p.legend(loc=2,prop={'size':6})
        p.set_title(db)

    tight_layout()
    savefig('ExpVar2.pdf')
    savefig('ExpVar2.png')

    figure(figsize=(5,3))
    p=subplot(111)
    ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
    for db in dbs:
        p = ps[db]
        p.plot(corrs,alphas[db])
        p.set_title(db)
        p.tick_params(axis='both', which='major', labelsize=6)
        p.tick_params(axis='both', which='minor', labelsize=6)
        p.set_ylim(0,2)
    tight_layout()
    savefig('ThresholdSlopes.pdf')
    savefig('ThresholdSlopes.png')

def variabilityAndGlobClustSlopesNormed():
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
    coli_data = {'Valgepea':ecoli_data_v,'Heinemann':ecoli_data_h}
    grs = {'Valgepea':gr_v,'Heinemann':gr_h}
    alphas = {'Valgepea':[],'Heinemann':[]}
    for db in ['Valgepea','Heinemann']:
        p=ps[db]
        conds = cond_list_dict[db]
        gr = grs[db]
        gr = gr[conds]
        glob_conc = coli_data[db]
        glob_data = glob_conc[conds]
        glob_data = glob_data.copy()
        tot_means = glob_data.mean(axis=1)
        for col in glob_data.columns:
            glob_data[col] = glob_data[col]/tot_means
        tot_var = calc_var(glob_data)
        print "tot_var is %f" % tot_var
        explained_glob = []
        explained_tot = []
        explained_compl_glob = []
        explained_compl_tot = []
        explained_normed = []
        explained_scaled = []
        for threshold in corrs:
            glob_cluster_idx = glob_conc['gr_cov']>threshold
            glob_compl_idx = glob_conc['gr_cov']<threshold
            alpha,glob_var,glob_explained,glob_scaled_explained = calc_explained_var(glob_data[glob_cluster_idx],gr)
            a,compl_var,compl_explained,compl_scaled_explained = calc_explained_var(glob_data[glob_compl_idx],gr)
            alphas[db].append(alpha)
            explained_var = glob_explained/glob_var
            explained_compl_var = compl_explained/compl_var
            explained_tot_frac = glob_explained/tot_var
            explained_compl_tot_frac = compl_explained/tot_var
            explained_tot_compl = (glob_explained+compl_explained)/tot_var
            explained_scaled_var = (glob_scaled_explained+compl_scaled_explained)/tot_var

            explained_glob.append(explained_var)
            explained_compl_glob.append(explained_compl_var)
            explained_tot.append(explained_tot_frac)
            explained_compl_tot.append(explained_compl_tot_frac)
            explained_normed.append(explained_tot_compl)
            explained_scaled.append(explained_scaled_var)

        p.plot(corrs,explained_glob,markersize=1,label='Explained variability fraction of global cluster')
        p.plot(corrs,explained_tot,markersize=1,label='Explained variability fraction of total data')
        p.plot(corrs,explained_compl_glob,markersize=1,label='Explained complementary variability fraction of global cluster')
        p.plot(corrs,explained_compl_tot,markersize=1,label='Explained complementary variability fraction of total data')
        p.plot(corrs,explained_normed,markersize=1,label='Explained variability fraction when normalizing')
        p.plot(corrs,explained_scaled,markersize=1,label='Explained variability fraction when scaling')
        p.set_ylabel('Explained fraction of variability', fontsize=8)
        p.set_xlabel('global cluster correlation threshold', fontsize=8)
        p.set_ylim(0,1)
        p.tick_params(axis='both', which='major', labelsize=6)
        p.tick_params(axis='both', which='minor', labelsize=6)
        p.axhline(xmin=0,xmax=1,y=0.5,ls='--',color='black',lw=0.5)
        p.legend(loc=2,prop={'size':6})
        p.set_title(db)

    tight_layout()
    savefig('ExpVar3.pdf')
    savefig('ExpVar3.png')

    figure(figsize=(5,3))
    p=subplot(111)
    ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
    for db in dbs:
        p = ps[db]
        p.plot(corrs,alphas[db])
        p.set_title(db)
        p.tick_params(axis='both', which='major', labelsize=6)
        p.tick_params(axis='both', which='minor', labelsize=6)
        p.set_ylim(0,2)
    tight_layout()
    savefig('ThresholdSlopes2.pdf')
    savefig('ThresholdSlopes2.png')


#6 panel graph - avg. exp. vs norm. slope, slope vs. r^2. non-global cluster avg. exp. vs. slope.
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
    conds = cond_list_dict['Heinemann']

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
    glob_conc = set_alpha('Heinemann',glob_conc,gr)
    glob_conc = set_std_err('Heinemann',glob_conc,gr)

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
    savefig('AvgConcStatsHein.png')

#comulative graph - x axis - avg. prot. conc. (or molecule count per cell), y axis, comulative % out of proteome.
def plotComulativeGraph():
    figure(figsize=(5,3))
    p1 = subplot(121)
    p2 = subplot(122)
    conds = cond_lists['Heinemann']
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
    savefig('DistStatsHein.png')

#plot the graphs for the 10 highest abundance proteins with their descriptions.
def plotHighAbundance():
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(121),'Valgepea':subplot(122)}
    for db in dbs:
        p = ps[db]
        conds = cond_lists[db]
        coli_data = coli_datas[db].copy()
        gr = grs[db]
        gr = gr[conds]
        coli_data['avg']=coli_data[conds].mean(axis=1)
        coli_data = coli_data.sort('avg',ascending=False)
        coli_data_conds = coli_data[conds]
        coli_data_conds = coli_data_conds.head(7)
        for i in coli_data_conds.index:
            desc = coli_data.ix[i]
            desc = "%s: %s: %s" % (desc['func'],desc['prot'],desc['ID'])
            p.plot(gr.values,coli_data_conds.ix[i].values,label=('%s' % desc))
        p.legend(loc=2, prop={'size':3},numpoints=1)
        p.set_ylim(0,0.1)
    tight_layout()
    savefig('highest.pdf')


#randomly select a few proteins and plot their prediction vs the actual concentration of a different protein in the HC prots.

def plotPrediction():
    for db in dbs:
        figure(figsize=(5,5))
        conds = cond_lists[db]
        coli_data = coli_datas[db]
        gr = grs[db]
        gr = gr[conds]
        glob = get_glob(db,coli_data)
        for i in range(1,10):
            p = subplot(330+i)
            samp = random.sample(glob.index,11)
            pred = samp[0:-1]
            est = samp[-1]
            pred = glob.ix[pred]
            est = glob.ix[est]
            pred = pred[conds]
            pred = pred.sum()
            pred = pred/pred.mean()
            est = est[conds]
            est = est/est.mean()
            alpha,beta,r_val,p_val,std_err = linregress(gr,pred)
            linpred = {}
            for c in gr.index:
                linpred[c]=alpha*gr[c]+beta
            linpred = pd.Series(linpred)
            linpred = linpred[conds]
            p.plot(gr.values,linpred.values,color='blue')
            p.plot(gr.values,pred,'o',color='blue',markersize=2)
            p.plot(gr.values,est,'o',color='green',markersize=2)
            p.set_ylim(0,3)
            p.set_xlim(0,0.7)
            p.tick_params(axis='both', which='major', labelsize=8)
            p.tick_params(axis='both', which='minor', labelsize=8)
            p.set_title("%.2f" % est.corr(linpred)**2,fontsize=8)
        tight_layout()    
        savefig('RandEstimate%s.pdf' % db)
        savefig('RandEstimate%s.png' % db)

# k-means

writeTables()
#plotCorrelationHistograms()
#plotGlobalResponse()
#plotMultiStats()
#plotComulativeGraph()
plotHighAbundance()
#plotPrediction()        
#variabilityAndGlobClustSlopes()
#variabilityAndGlobClustSlopesNormed()
