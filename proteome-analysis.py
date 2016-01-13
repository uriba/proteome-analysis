import csv
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
from scipy import stats
from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust,text,subplots,gcf,close,xscale,yscale
from matplotlib.markers import MarkerStyle
from numpy import linspace,ndarray,arange,sum,square,array,cumsum,ones,mean,std,cov
import numpy as np
from numpy.random import randn
from analysis import *
import matplotlib
from math import sqrt,isnan,log
import random
from numpy.random import randn,shuffle,normal
from matplotlib.ticker import FuncFormatter
from itertools import combinations

random.seed(123456)
#import plotly.plotly as py

#py.sign_in("uri.barenholz", "hvi3ma3m30")
### Results generation#####
def get_limits(db):
    if db == 'Heinemann':
        limits = (0.25,1.)
    if db == 'HeinemannLB':
        limits = (0.6,1.)
    if db == 'Valgepea':
        limits = (0.8,1.)
    if db == 'Heinemann-chemo':
        limits = (0.8,1.)
    #return limits
    return (0.5,1.)
    #return (-1.,-0.5)

#Initialize global data structures
#dbs = ['Heinemann','HeinemannLB','Peebo','Valgepea']
dbs = ['Heinemann','HeinemannLB','Heinemann-chemo','Peebo','Peebo-gluc','HuiAlim','HuiClim','HuiRlim','Valgepea']
datas = {}
rand_prefix = ""
db_name = { 'Heinemann':'Schmidt','HeinemannLB':'Schmidt','Heinemann-chemo':'Schmidt','Valgepea':'Valgepea','Peebo-gluc':'Peebo','Peebo':'Peebo','HuiAlim':'Hui','HuiClim':'Hui','HuiRlim':'Hui'}
db_suffix = { 'Heinemann':'','HeinemannLB':'Rich media','Heinemann-chemo':'chemo.','Valgepea':'','Peebo-gluc':'gluc.','Peebo':'','HuiAlim':'A-lim','HuiClim':'C-lim','HuiRlim':'R-lim'}
def init_datasets(rand_method):
    datas[rand_method] = {}
    for db in dbs:
        (conds,gr,coli_data) = get_annotated_prots(db,rand_method)
        gr = gr[conds]
        coli_data['avg']=coli_data[conds].mean(axis=1)
        coli_data['std']=coli_data[conds].std(axis=1)
        coli_data = calc_gr_corr(coli_data,conds,gr)
        CV = (coli_data['std']/coli_data['avg']).mean()
        datas[rand_method][db] = (conds,gr,coli_data)
        print "in %s, %s CV: %f" % (rand_method,db,CV)

#ecoli_data_h = ecoli_data_h[ecoli_data_h['prot']=='Ribosome']
#ecoli_data_v = ecoli_data_v[ecoli_data_v['prot']=='Ribosome']


#write tables data:
def writeCorrsHist(db):
    conds,gr,conc_data = datas[""][db]
    limits = get_limits(db)
    threshold = limits[0]
    funcs = conc_data['func'].unique()
    idx = ['Function','Number of proteins','totPrctP','Correlated proteins','corPrctP']
    func_stat = pd.DataFrame(columns = idx)
    for func in funcs:
        conc_func = conc_data[conc_data['func']==func]
        corred_idx = conc_func['gr_cov']>threshold
        tot = len(conc_func)
        tot_means = conc_func['avg']
        corr_means = tot_means[corred_idx]
        correlated = len(corr_means)
        stats = pd.Series([func,tot,tot_means.sum()*100,correlated,corr_means.sum()*100],index=idx)
        func_stat.loc['{%s}' % func]=stats
    func_stat.sort(columns = 'totPrctP',ascending=False,inplace = True)
    func_stat.to_csv('funcs%s.csv' % db,index=False,sep=';')

def writeTopProtsVar(db):
    conds,gr,conc_data = datas[""][db]
    conc_data = conc_data.copy()
    for cond in conds:
        conc_data[cond]=conc_data[cond]-conc_data['avg']
    conc_data_vars = (conc_data[conds]**2).sum(axis=1)
    conc_data['vars']=conc_data_vars
    tot_vars = conc_data['vars'].sum()
    conc_data = conc_data.sort('avg',ascending=False)
    if db == 'Heinemann' or db == 'HeinemannLB' or db == 'Heinemann-chemo':
        conc_data['Temp']=conc_data['protName']
    high_abdc = conc_data.head(20)
    with open('%svarsOfAbdcs%s.csv' % (rand_method,db),'wb') as csvfile:
        csvwriter = csv.writer(csvfile,delimiter=';')
        csvwriter.writerow(['Function','Sub Function','Name','totPrctP','prctOfVar','cov'])
        j = 0
        for i,row in high_abdc.iterrows():
            csvwriter.writerow((row['func'], row['prot'],row['Temp'],row['avg']*100,row['vars']*100/tot_vars,row['gr_cov']))
            #plot protein with second highest abundance in valgepea data set
            if (j == 1) and (db == 'Valgepea'):
                figure(figsize=(5,3))
                ax = subplot(111)
                ax.plot(gr,100*(row[conds]+row['avg']),'o',label="metE, Correlation: %.2f" % gr.corr(row[conds]))
                ax.set_ylim(0,5)
                ax.set_xlim(0,0.6)
                ax.set_xlabel("Growth rate [$h^{-1}$]")
                ax.set_ylabel("% of total proteome")
                legend(loc=2, prop={'size':8},numpoints=1)
                tight_layout()
                #fig = gcf()
                #py.plot_mpl(fig,filename="metE chemostat")
                savefig('SingleProt%s.pdf' % row['Temp'])
                close()
            j+=1
    conc_data = conc_data.sort('vars',ascending=False)
    high_vars = conc_data.head(20)
    with open('varsOfVars%s.csv' % db,'wb') as csvfile:
        csvwriter = csv.writer(csvfile,delimiter=';')
        csvwriter.writerow(['Function','Sub Function','Name','totPrctP','prctOfVar','cov'])
        for i,row in high_vars.iterrows():
            csvwriter.writerow((row['func'], row['prot'],row['Temp'],row['avg']*100,row['vars']*100/tot_vars,row['gr_cov']))

def writeTables():
    for db in dbs:
        writeCorrsHist(db)
        writeTopProtsVar(db)

### Figure 1 - Correlation to growth rate by functional group histogram.
categories = ['Metabolism','Genetic Information Processing','Environmental Information Processing', 'Cellular Processes','NotMapped']
categories = ['Metabolism','Information storage and processing', 'Cellular processes and signaling','Unknown']
def set_ticks(p,size):
    p.tick_params(axis='both', which='major', labelsize=size)
    p.tick_params(axis='both', which='minor', labelsize=size)

def plot_corr_hist(p,db,conc_data,categories):
    bins = linspace(-1,1,21)
    covs = ndarray(shape=(len(categories),len(bins)-1))
    sets = [] 

    for x in categories:
        sets.append(conc_data[conc_data['group']==x].gr_cov)
    if len(categories)>1:
        p.hist(sets,bins = bins, stacked = True,label=categories,zorder=0)
    else:
        p.hist(conc_data.gr_cov,bins = bins, color='0.75',zorder=0)

    set_ticks(p,6)
    p.set_xlabel('Pearson correlation with growth rate',fontsize=6)
    p.set_ylabel('Number of proteins',fontsize=6)
    for limit in get_limits(db):
        p.axvline(x=limit,ymin=0,ymax=250,ls='--',color='black',lw=0.5)

def plotCorrelationHistograms(dbs,suffix):
    figure(figsize=(6,2.6))

    coords = {'Heinemann':0.01,'Peebo':0.627,'Valgepea':0.627,'HuiAlim':0.01,'HuiClim':0.627,'HuiRlim':0.627}
    p=subplot(111)
    rands = [""]
    ps = {("",'Peebo'):subplot(132),("",'Valgepea'):subplot(132),("",'HuiAlim'):subplot(121),("",'HuiClim'):subplot(122),("",'HuiRlim'):subplot(122)}
    if(len(dbs)>1):
        ps['Heinemann'] = subplot(131)
        rands = ["","shuffle"]
        ps = {  ("",'Heinemann'):subplot(131),
                ("",'Peebo'):subplot(132),
                ("",'Valgepea'):subplot(132),
                ("",'HuiAlim'):subplot(131),
                ("",'HuiClim'):subplot(132),
                ("",'HuiRlim'):subplot(132),
                ("shuffle",'Heinemann'):subplot(133),
                ("shuffle",'Peebo'):subplot(133),
                ("shuffle",'Valgepea'):subplot(131)}
        horiz = {   ("",'Heinemann'):-0.027,
                    ("",'Peebo'):0.392,
                    ("shuffle",'Peebo'):0.81}

        for rand,db,panel in [("","Heinemann",'A'),('','Peebo','B'),('shuffle','Peebo','C')]:
            conds,gr,conc_data = datas[rand][db]
            if rand == 'shuffle':
                conc_data['group'] = ""
                plot_corr_hist(ps[(rand,db)],db,conc_data,[""])
                ps[(rand,db)].annotate("shuffled data",xy=(0.5,0.5),xytext=(-0.94,210),fontsize=6,zorder=10)
            else:
                plot_corr_hist(ps[(rand,db)],db,conc_data,categories)
                ps[(rand,db)].annotate("data from %s et. al. 2015" % db_name[db],xy=(0.5,0.5),xytext=(-0.94,210),fontsize=6,zorder=10)
            ps[(rand,db)].set_ylim(0,250)
            ps[(rand,db)].set_xlim(-1,1)
            ps[(rand,db)].annotate(panel,xy=(0.5,0.5),xytext=(-0.9,230),fontsize=10,zorder=10)

    #assume both subplots have the same categories.
    handles,labels=ps[("",dbs[0])].get_legend_handles_labels()

    tight_layout()
    figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.15,0.8,0.7,0.2),ncol=2)

    subplots_adjust(top=0.85)
    #fig = gcf()
    #py.plot_mpl(fig,filename="Growth rate Correlation histograms")
    savefig('%sGrowthRateCorrelation%s.pdf' % (rand_prefix,suffix))
    if rand_prefix == suffix:
        savefig('Fig2.eps')
    close()

##############temp 20 prots plot####################3
def tempprotsplot():
    db="Heinemann"
    size=6
    fignum=12
    rows=3
    columns=4
    conds,gr,coli_data = datas[db]
    slowconds= [ u'Chemostat mu=0.12', u'Chemostat mu=0.20', u'Galactose', u'Acetate', u'Chemostat mu=0.35', u'Pyruvate', u'Fumarate', 
     u'Succinate', u'Glucosamine', u'Glycerol', u'Mannose', u'Chemostat mu=0.5', u'Xylose', u'Osmotic-stress glucose', u'Glucose',
     u'pH6 glucose', u'Fructose', u'42C glucose', ]
    errorconds = ["%s.cv" % x for x in conds]
    glob = get_glob(db,coli_data)
    sampled = random.sample(glob.index,size)
    figure(figsize=(5,5))
    globprots = coli_data.loc[sampled,:]
    j=0
    for i,row in globprots.iterrows():
        p = subplot(rows,columns,(j % fignum)+1)
        y = row[conds]/row[slowconds].mean()
        ycv = row[errorconds]
        ycv.index = conds
        p.errorbar(gr[conds],y,fmt='b.',markersize=2,yerr=y*ycv/100,capthick=0.5,elinewidth=0.5,capsize=2)
        j+=1
    coli_data['fast_cov'] = coli_data['gr_cov']
    coli_data = calc_gr_corr(coli_data,slowconds,gr)
    glob = get_glob(db,coli_data)
    glob = glob[glob['fast_cov']<0.5]
    sampled = random.sample(glob.index,size)
    globprots = coli_data.loc[sampled,:]
    j=6
    for i,row in globprots.iterrows():
        p = subplot(rows,columns,(j % fignum)+1)
        y = row[conds]/row[slowconds].mean()
        ycv = row[errorconds]
        ycv.index = conds
        p.errorbar(gr[conds],y,fmt='r.',markersize=2,yerr=y*ycv/100,capthick=0.25,elinewidth=0.25,capsize=1)
        j+=1
    for j in range(fignum):
        p = subplot(rows,columns,j+1)
        p.set_ylim(0,3)
        set_ticks(p,8)

    tight_layout()
    savefig('slowvsfastcorrprotscomparison.pdf')
    close()
    coli_data = calc_gr_corr(coli_data,conds,gr)

### Figure 3, Global cluster analysis:
def plotGlobalResponse(dbs,rand_method):
    globalResponse[rand_method] = {}
    colors = {'Heinemann':'blue','Peebo':'green','Valgepea':'magenta','HuiAlim':'cyan','HuiClim':'gray','HuiRlim':'yellow'}

    for db in dbs:
        conds,gr,coli_data = datas[rand_method][db]
        glob = get_glob(db,coli_data)
        print "%s global cluster is %d out of %d measured proteins" % (db, len(glob),len(coli_data[coli_data['gr_cov']>-1.]))

        glob_tot = glob[conds].sum()
        alpha,beta,r_val,p_val,std_err = linregress(gr,glob_tot)
        print "global cluster sum follows alpha=%f, beta=%f" % (alpha,beta)
        print "horizontal intercept for %s is %f, corresponding to halflive %f" % (db,-beta/alpha, log(2)*alpha/beta)
        globalResponse[rand_method][db] = {}
        globalResponse[rand_method][db]['Rsq']=gr.corr(glob_tot)**2
        globalResponse[rand_method][db]['gr']=gr.values
        globalResponse[rand_method][db]['dots']=glob_tot.values
        globalResponse[rand_method][db]['line']=alpha*gr.values+beta

    if "shuffle" in globalResponse and "" in globalResponse:
        figure(figsize=(5,3))
        for db in dbs:
            print db_name
            plot(globalResponse[""][db]['gr'],globalResponse[""][db]['dots'],'o',label=("data from %s et. al, 2015 ($R^2$=%.2f)" % (db_name[db],globalResponse[""][db]['Rsq'])),color=colors[db])
            plot(globalResponse[""][db]['gr'],globalResponse[""][db]['line'],color=colors[db])
            plot(globalResponse["shuffle"][db]['gr'],globalResponse["shuffle"][db]['dots'],linestyle='None',marker="o",markerfacecolor='none',markeredgecolor=colors[db],label=("Shuffled protein amounts, based on %s" % db_name[db]))

        xlim(xmin=0.)
        ylim(ymin=-0.05,ymax=0.75)
        xlabel('Growth rate [$h^{-1}$]',fontsize=10)
        ylabel('Strongly correlated proteins\n fraction out of proteome',fontsize=10)
        legend(loc="upper left", prop={'size':6},numpoints=1)
        tick_params(axis='both', which='major', labelsize=8)
        tick_params(axis='both', which='minor', labelsize=8)
        tight_layout()
        #fig = gcf()
        #py.plot_mpl(fig,filename="Global cluster growth rate correlation")
        savefig('%sGlobalClusterGRFit.pdf' % rand_prefix)
        if rand_prefix == '':
            savefig('Fig4.eps')
        close()

#gets values at cond_list normalized in y axis
def std_err_fit(gr,s):
    alpha,beta,r,p,st = linregress(gr,s)
    return st
    
def conf_int_min(degfr,s):
    res = stats.t.interval(0.95,degfr,loc=s['alpha'],scale=s['std_err'])
    return  res[0]

def conf_int_max(degfr,s):
    res = stats.t.interval(0.95,degfr,loc=s['alpha'],scale=s['std_err'])
    return  res[1]

def set_std_err(df,gr,cond_list):
    if len(df)>0:
    #df['std_err'] = df[cond_list].apply(lambda x: std_err_fit(gr[cond_list]/gr[cond_list].mean(),x/x.mean()),axis=1)
        df['std_err'] = df[cond_list].apply(lambda x: std_err_fit(gr[cond_list],x/x.mean()),axis=1)
        
        df['conf_min'] = df.apply(lambda x: conf_int_min(len(cond_list)-2,x) ,axis=1)
        df['conf_max'] = df.apply(lambda x: conf_int_max(len(cond_list)-2,x) ,axis=1)
    return df

## Figure 2, global cluster slope vs. ribosomal slope
def get_glob(db,df):
    limits = get_limits(db)
    return df[(df['gr_cov']<limits[1]) & (df['gr_cov']>limits[0])]
 
def set_alpha(df,gr,cond_list):
    if len(df)>0:
    #df['alpha'] = df[cond_list].apply(lambda x: linregress(gr[cond_list]/gr[cond_list].mean(),x/x.mean())[0],axis=1)
        df['alpha'] = df[cond_list].apply(lambda x: linregress(gr[cond_list],x/x.mean())[0],axis=1)
    return df

def plot_response_hist(db,df,gr,conds,p,total,estimate):
    all_ribs = df[df['prot']=='Ribosome']
    all_ribs = set_alpha(all_ribs,gr,conds)
    all_ribs = set_std_err(all_ribs,gr,conds)
    if not total:
        print "total ribosomal proteins in db %s, %d, strongly positively correlated: %d" % (db,len(all_ribs),len(all_ribs[all_ribs['gr_cov']>get_limits(db)[0]]))

    bins = linspace(-5,5,41)
    xs = linspace(-5,5,200)
    glob_conc = get_glob(db,df)
    glob_conc = set_alpha(glob_conc,gr,conds)
    print "out of %d proteins, %d have slopes in the range 0.5 to 2" % (len(glob_conc),len(glob_conc[(glob_conc['alpha']>=0.5) & (glob_conc['alpha']<=2)]))
    glob_conc = set_std_err(glob_conc,gr,conds)
    avg = glob_conc['alpha'].mean()
    std_err = glob_conc['std_err'].mean()
    if not total:
        glob_conc_no_ribs = glob_conc[glob_conc['prot'] != 'Ribosome']
        ribs = glob_conc[glob_conc['prot'] == 'Ribosome']
        print "mean slope of strongly correlated prots: %.2f, std. dev. of slopes: %.2f, mean std.err of slope %.2f" % (glob_conc['alpha'].mean(), glob_conc['alpha'].std(), glob_conc['std_err'].mean())
        print "mean slope of ribosomal: %d proteins, %.2f, std. dev. of slopes: %.2f, mean std.err of slope %.2f" % (len(ribs),ribs['alpha'].mean(), ribs['alpha'].std(), ribs['std_err'].mean())
        idx = list(glob_conc.index)
        size = len(ribs)
        alphas = []
        for i in range(1000):
            shuffle(idx) 
            elems = idx[1:size]
            alphas.append(glob_conc.loc[elems,'alpha'].mean())
        print "for db %s, mean slope of randomly sampled data is %.2f, std_dev of slopes is %.2f. " % (db,mean(alphas),std(alphas))
        (alpha,b,r,pv,st)   = linregress(gr[conds], glob_conc[conds].sum()/glob_conc[conds].sum().mean())
        print "slope of sum of prots: %.2f, r-sq, %.2f" % (alpha,r**2)
        (ribs_alpha,b,ribs_r,pv,st)   = linregress(gr[conds], ribs[conds].sum()/ribs[conds].sum().mean())
        print "slope of sum of ribosomal prots: %.2f, r-sq, %.2f" % (ribs_alpha,ribs_r**2)
        p.hist([glob_conc_no_ribs['alpha'].values,ribs['alpha'].values],bins=bins,stacked = True,label=['High correlation proteins','Ribosomal proteins'],color=['blue','#20ff20'])
    else:
        p.hist(glob_conc['alpha'].values,bins=bins,label=['High correlation proteins'])
    if estimate:
        p.plot(xs,stats.t.pdf(xs,df=len(conds)-2,loc=avg,scale=std_err)*len(glob_conc['alpha'])*0.25,linestyle='-',color='0.5')
    p.set_xlim(-5,5)
    for x in (0.5,2):
        p.axvline(x=x,ymin=0,ymax=100,ls='--',color='black',lw=0.5)
    p.set_xlabel('Normalized slope',fontsize=8)
    p.set_ylabel('Number of proteins',fontsize=8)
    set_ticks(p,8)

def plot_response_hist_graphs(dbs):
    plots = {"AllProtsNormalizedSlopes":(True,False),"AllProtsVSRibosomalNoExpNormalizedSlopes":(False,False),"AllProtsVSRibosomalNormalizedSlopes":(False,True)}
    for (name,vals) in plots.iteritems():
        figure(figsize=(5,3))
        p = subplot(111)
        ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122),'HuiAlim':subplot(121),'HuiClim':subplot(122),'HuiRlim':subplot(122)}
        coords = {'Heinemann':0.0,'Peebo':0.62,'Valgepea':0.62,'HuiAlim':0.0,'HuiClim':0.62,'HuiRlim':0.62}
        for db in dbs:
            conds,gr,conc_data = datas[""][db]
            plot_response_hist(db,conc_data,gr,conds,ps[db],vals[0],vals[1])
            text(coords[db],0.93,"data from %s et. al" % db_name[db],fontsize=8,transform=p.transAxes)
            handles,labels=ps[db].get_legend_handles_labels()
            if db in ['Peebo','Valgepea']:
                ps[db].set_ylim(0,100)

        figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.25,0.8,0.5,0.2),ncol=2)
        tight_layout()
#fig = gcf()
#py.plot_mpl(fig,filename="Normalized slopes distribution")
        savefig('%s%s.pdf' % (rand_prefix,name))
        if name == 'AllProtsVSRibosomalNormalizedSlopes':
            savefig('Fig3.eps')
        close()


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
    
def corr_andGR_plot(db,ref):
    suffix = 'Chem'
    if db == "HeinemannLB":
        suffix = 'LB'
    rand = ''
    if db == 'Simulated':
        rand = 'simulated'
        db = 'Heinemann'
        suffix = 'simulated'
    figure(figsize=(5,3))
    (cond_list,gr_chemo,ecoli_data_chemo) = get_annotated_prots(db,rand)
    ecoli_data_chemo = calc_gr_corr(ecoli_data_chemo,cond_list,gr_chemo)

    p1=subplot(121)
    p2=subplot(122)
    plot_corr_hist(p1,db,ecoli_data_chemo,categories)

    handles,labels=p1.get_legend_handles_labels()
    figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.05,0.8,0.7,0.2),ncol=2)
    glob_chemo = get_glob(db,ecoli_data_chemo)
    print "%s global cluster is %d out of %d measured proteins" % (db, len(glob_chemo),len(ecoli_data_chemo[ecoli_data_chemo['gr_cov']>-1.]))

    glob_tot_chemo = glob_chemo[cond_list].sum()
    alpha,beta,r_val,p_val,std_err = linregress(gr_chemo,glob_tot_chemo)

    print "global cluster sum follows alpha=%f, beta=%f" % (alpha,beta)
    print "horizontal intercept for %s is %f, corresponding to halflive %f" % (db,-beta/alpha, log(2)*alpha/beta)
    p2.plot(gr_chemo.values,glob_tot_chemo.values,'o',label="%s et. al %s" % (db_name[db],suffix),color='blue')
    p2.plot(gr_chemo.values,alpha*gr_chemo.values+beta,color='blue',label=("%s %s. Trend,$R^2$=%.2f" % (db_name[db],suffix,gr_chemo.corr(glob_tot_chemo)**2)))

    cond_list,gr_v,conc_data = datas[""][ref]
    glob_v = get_glob(ref,conc_data)
    glob_tot_v = glob_v[cond_list].sum()
    alpha_v,beta_v,r_val,p_val,std_err = linregress(gr_v,glob_tot_v)
    p2.plot(gr_v.values,glob_tot_v.values,'o',label=db_name[ref],color='green')
    p2.plot(gr_v.values,alpha_v*gr_v.values+beta_v,color='green',label=("%s Trend,$R^2$=%.2f" % (db_name[ref],gr_v.corr(glob_tot_v)**2)))

    p2.set_xlim(xmin=0.)
    p2.set_ylim(ymin=0.)
    p2.set_xlabel('Growth rate',fontsize=8)
    p2.set_ylabel('Strongly correlated proteins\n fraction out of proteome',fontsize=8)
    legend(loc=3, prop={'size':6},numpoints=1)
    set_ticks(p2,8)
    tight_layout()

    subplots_adjust(top=0.83)
#fig = gcf()
#py.plot_mpl(fig,filename="Heinemann chemostat graphs")
    savefig('%s%ssummaryHistAndGr.pdf' % (db,rand))
    if db == 'HeinemannLB' and rand == '':
        savefig('FigS5.eps')
    if db == 'Heinemann' and rand == 'simulated':
        savefig('FigS9.eps')
    close()

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
#set_ticks(p2,8)
#tight_layout()

#subplots_adjust(top=0.83)
#savefig('Anticorrelated.pdf')

#check if for 95% of the slopes, the mean of all of the slopes lies in their 95% confidence interval

#Plot variability explained (R^2)/Var? in global cluster and in proteome as function of threshold for HC proteins.

corrs = linspace(-1,1,100)

# calculate variability explained in proteome, take 1 (1 free parameter - selection of global cluster and scaling accordingly.
# calculate variability explained in global cluster, take 2 (1 free parameter - selection of global cluster and measurement of resulting variability reduction.
def square_dist_func(df):
    return df**2

def abs_dist_func(df):
    return abs(df)

def calc_var(f,df,means):
    df = df.copy()
    for col in df.columns:
        df[col]=df[col]-means
    var = f(df)
    var = var.sum().sum()
    return var

def calc_explained_var(f,df,means,gr):
    tot_var = calc_var(f,df,means)
    df = df.copy()
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
    remained_var = calc_var(f,remains,remains.mean(axis=1))
    scaled_remains = df-scaled
    remained_scaled_var = calc_var(f,scaled_remains,scaled_remains.mean(axis=1))
    alpha,beta,r_val,p_val,std_err = linregress(gr,response/response.mean())
    return (alpha,tot_var,tot_var - remained_var,tot_var-remained_scaled_var)

def calc_var_stats(f,conds,gr,glob_conc):
    alphas = []
    glob_data = glob_conc[conds]
    tot_var = calc_var(f,glob_data,glob_conc['avg'])
    print "tot_var is %f" % tot_var
    explained_glob = []
    explained_tot = []
    explained_compl_glob = []
    explained_compl_tot = []
    explained_scaled = []
    glob_frac = []
    for threshold in corrs:
        glob_cluster_idx = glob_conc['gr_cov']>threshold
        glob_compl_idx = glob_conc['gr_cov']<threshold
        alpha,glob_var,glob_explained,glob_scaled_explained = calc_explained_var(f,glob_data[glob_cluster_idx],(glob_conc[glob_cluster_idx])['avg'],gr)
        a,compl_var,compl_explained,compl_scaled_explained = calc_explained_var(f,glob_data[glob_compl_idx],glob_conc[glob_compl_idx]['avg'],gr)
        alphas.append(alpha)
        explained_var = glob_explained/glob_var
        explained_compl_var = compl_explained/compl_var
        explained_tot_frac = glob_explained/tot_var
        explained_compl_tot_frac = compl_explained/tot_var
        explained_scaled_var = (glob_scaled_explained+compl_scaled_explained)/tot_var

        explained_glob.append(explained_var)
        explained_compl_glob.append(explained_compl_var)
        explained_tot.append(explained_tot_frac)
        explained_compl_tot.append(explained_compl_tot_frac)
        explained_scaled.append(explained_scaled_var)
        glob_frac.append(float(len(glob_data[glob_cluster_idx]))/len(glob_data))

    return (explained_glob,explained_tot,explained_compl_glob,explained_compl_tot,explained_scaled,alphas,glob_frac)

def variabilityAndGlobClustSlopes(dbs):
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122)}
    alphas = {}
    for db in dbs:
        alphas[db]=[]
        p=ps[db]
        conds,gr,glob_conc = datas[db]
        (explained_glob,explained_tot,explained_compl_glob,explained_compl_tot,explained_scaled,alphas[db],x) = calc_var_stats(square_dist_func,conds,gr,glob_conc)
        p.plot(corrs,explained_glob,markersize=1,label='Explained variability fraction of global cluster')
        p.plot(corrs,explained_tot,markersize=1,label='Explained variability fraction of total data')
        p.plot(corrs,explained_compl_glob,markersize=1,label='Explained complementary variability fraction of global cluster')
        p.plot(corrs,explained_compl_tot,markersize=1,label='Explained complementary variability fraction of total data')
        explained_normed = [x+y for x,y in zip(explained_tot,explained_compl_tot)]
        p.plot(corrs,explained_normed,markersize=1,label='Explained variability fraction when normalizing')
        p.plot(corrs,explained_scaled,markersize=1,label='Explained variability fraction when scaling')
        p.set_ylabel('Explained fraction of variability', fontsize=8)
        p.set_xlabel('global cluster correlation threshold', fontsize=8)
        p.set_ylim(0,1)
        set_ticks(p,6)
        p.axhline(xmin=0,xmax=1,y=0.5,ls='--',color='black',lw=0.5)
        p.legend(loc=2,prop={'size':6})
        p.set_title(db)

    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Non normalized variability statistics")
    savefig('%sExpVar2.pdf' % rand_prefix)
    close()

    figure(figsize=(5,3))
    p=subplot(111)
    ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122)}
    for db in dbs:
        p = ps[db]
        p.plot(corrs,alphas[db])
        p.set_title(db)
        set_ticks(p,6)
        p.set_ylim(0,2)
    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Global response slope dependence on threshold")
    savefig('%sThresholdSlopes.pdf' % rand_prefix)
    close()

def norm_glob_conc(glob_conc,conds):
    glob_conc = glob_conc.copy()
    tot_means = glob_conc['avg']
    for col in conds:
        glob_conc[col] = glob_conc[col]/tot_means
    glob_conc['avg'] = glob_conc[conds].mean(axis=1)
    return glob_conc

def keep_middle(glob_conc,conds):
    glob_conc = glob_conc.copy().sort('avg',ascending=False)
    num = len(glob_conc)
    #glob_conc = glob_conc[:-num/8]
    glob_conc = glob_conc[num/4:]
    return glob_conc

def drop_head(glob_conc,conds,num):
    glob_conc = glob_conc.copy().sort('avg',ascending=False)
    glob_conc = glob_conc[num:]
    return glob_conc

def variablityComparisonHein():
    figure(figsize=(8,5))
    db = 'Heinemann'
    conds,gr,glob_conc = datas[db]
    globs = []
    titles = []
    funcs = [square_dist_func,square_dist_func,square_dist_func,abs_dist_func,abs_dist_func,square]
    globs.append(glob_conc)
    titles.append('all prots')
    globs.append(norm_glob_conc(globs[0],conds))
    titles.append('all prots, normalized')
    globs.append(keep_middle(globs[0],conds))
    titles.append('prots excl. top $\\frac{1}{4}$')
    globs.append(globs[0])
    titles.append('all prots, abs')
    globs.append(globs[2])
    titles.append('prots excl. top $\\frac{1}{4}$, abs')
    globs.append(drop_head(globs[0],conds,10))
    titles.append('prots excl. top 10')
    for i in range(0,6):
        p = subplot(231+i)
        (explained_glob,explained_tot,explained_compl_glob,explained_compl_tot,explained_scaled,temp,x) = calc_var_stats(funcs[i],conds,gr,globs[i])

        p.plot(corrs,explained_glob,markersize=1,label='Explained variability fraction of global cluster')
        p.plot(corrs,explained_tot,markersize=1,label='Explained variability fraction of total data')
        p.plot(corrs,explained_compl_glob,markersize=1,label='Explained complementary variability fraction of global cluster')
        p.plot(corrs,explained_compl_tot,markersize=1,label='Explained complementary variability fraction of total data')
        explained_normed = [x+y for x,y in zip(explained_tot,explained_compl_tot)]
        p.plot(corrs,explained_normed,markersize=1,label='Explained variability fraction when normalizing')
        p.plot(corrs,explained_scaled,markersize=1,label='Explained variability fraction when scaling')
        p.set_ylabel('Explained fraction of variability', fontsize=8)
        p.set_xlabel('global cluster correlation threshold', fontsize=8)
        p.set_ylim(0,1)
        set_ticks(p,6)
        p.axhline(xmin=0,xmax=1,y=0.5,ls='--',color='black',lw=0.5)
        if i==0:
            p.legend(loc=2,prop={'size':6})
        p.set_title(titles[i],fontsize=8)

    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Various heuristics on explained variability for Heinemann data set")
    savefig('%sExpVarComp.pdf' % rand_prefix)
    close()

def variabilityAndGlobClustSlopesNormed(dbs,rand_method):
    figure(figsize=(5,3))
    p=subplot(111)
    ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122)}
    coords = {'Heinemann':0.03,'Peebo':0.03,'Valgepea':0.03}

    alphas = {}
    for db in dbs:
        alphas[db] = []
        p=ps[db]
        conds,gr,glob_conc = datas[rand_method][db]
        glob_conc = glob_conc.copy()
        tot_means = glob_conc['avg']
        for col in conds:
            glob_conc[col] = glob_conc[col]/tot_means
        glob_conc['avg'] = glob_conc[conds].mean(axis=1)
        (explained_glob,explained_tot,explained_compl_glob,explained_compl_tot,explained_scaled,alphas[db],glob_frac) = calc_var_stats(square_dist_func,conds,gr,glob_conc)

        p.plot(corrs,explained_tot,markersize=1,label='Explained variability fraction of total data')
        p.plot(corrs,explained_glob,markersize=1,label='Explained variability fraction of global cluster')
        p.plot(corrs,glob_frac,markersize=1,label='Correlated proteins fraction of proteome')
        #p.plot(corrs,explained_compl_glob,markersize=1,label='Explained complementary variability fraction of global cluster')
        #p.plot(corrs,explained_compl_tot,markersize=1,label='Explained complementary variability fraction of total data')
        explained_normed = [x+y for x,y in zip(explained_tot,explained_compl_tot)]
        #p.plot(corrs,explained_normed,markersize=1,label='Explained variability fraction when normalizing')
        #p.plot(corrs,explained_scaled,markersize=1,label='Explained variability fraction when scaling')
        p.set_ylabel('Explained fraction of variability', fontsize=8)
        p.set_xlabel('global cluster correlation threshold', fontsize=8)
        p.set_ylim(0,1)
        set_ticks(p,6)
        print "Maximum explained variability for db %s is %f" % (db,max(explained_tot))
        if db == 'Heinemann':
            exp_var = 0.09
        else:
            exp_var = 0.04
        p.axhline(xmin=0,xmax=1,y=exp_var,ls='--',color='black',lw=0.5)
        #p.axvline(ymin=0,ymax=1,x=get_limits(db)[0],ls='--',color='black',lw=0.5)
        text(coords[db],0.9,"data from %s et. al." % db_name[db],fontsize=8,transform=p.transAxes)

    handles,labels=ps['Heinemann'].get_legend_handles_labels()

    figlegend(handles,labels,fontsize=6,loc='upper left',bbox_to_anchor=(0.2,0.8,0.6,0.2))

    tight_layout()
    subplots_adjust(top=0.83)

    #fig = gcf()
    #py.plot_mpl(fig,filename="Explained variability statistics on normalized concentrations")
    savefig('%sExpVar3.pdf' % rand_method)
    if rand_method == '':
        savefig('FigS4.eps')
    if rand_method == 'shuffle':
        savefig('FigS8.eps')
    close()

    figure(figsize=(5,3))
    p=subplot(111)
    ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122)}
    for db in dbs:
        p = ps[db]
        p.plot(corrs,alphas[db])
        p.set_title(db)
        set_ticks(p,6)
        p.set_ylim(0,2)
    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Dependence on threshold of global response slopes for normalized concentrations")
    savefig('%sThresholdSlopes2.pdf' % rand_prefix)
    close()


#6 panel graph - avg. exp. vs norm. slope, slope vs. r^2. non-global cluster avg. exp. vs. slope.
def plotMultiStats(db):
    figure(figsize=(5,3))
    conds,gr,glob_conc = datas[db]
    sp = []
    for i in range(6):
        sp.append(subplot(231+i))

    sp[0].plot(glob_conc['avg'], glob_conc['rsq'],'.', markersize=1)
    sp[0].set_xlabel('Average concentraion', fontsize=6)
    sp[0].set_ylabel('$R^2$ with GR', fontsize=6)

    sp[1].plot(glob_conc['avg'], glob_conc['gr_cov'],'.', markersize=1)
    sp[1].set_xlabel('Average concentraion', fontsize=6)
    sp[1].set_ylabel('Pearson corr. with GR', fontsize=6)

    glob_conc = glob_conc[glob_conc['gr_cov']>get_limits(db)[0]]
    glob_conc = set_alpha(glob_conc,gr,conds)
    glob_conc = set_std_err(glob_conc,gr,conds)

    sp[2].plot(glob_conc['avg'], glob_conc['alpha'],'.', markersize=1)
    sp[2].set_xlabel('Average concentraion (HC prots)', fontsize=6)
    sp[2].set_ylabel('Norm. Slope', fontsize=6)

    sp[3].plot(glob_conc['avg'], glob_conc['std_err'],'.', markersize=1)
    sp[3].set_xlabel('Average concentraion (HC prots)', fontsize=6)
    sp[3].set_ylabel('std err of fit', fontsize=6)

    for i in range(4):
        sp[i].set_xscale('log')

    sp[4].plot(glob_conc['alpha'], glob_conc['std_err'],'.', markersize=1)
    sp[4].set_xlabel('Norm. slope (HC)', fontsize=6)
    sp[4].set_ylabel('std err of fit', fontsize=6)

    sp[5].plot(glob_conc['alpha'], glob_conc['rsq'],'.', markersize=1)
    sp[5].set_xlabel('Norm. slope (HC)', fontsize=6)
    sp[5].set_ylabel('$R^2$ with GR', fontsize=6)

    for i in range(6):
        set_ticks(sp[i],6)
    tight_layout()

    #fig = gcf()
    #py.plot_mpl(fig,filename="Proteins statistics for Heinemann dataset")
    glob_conc.to_csv('stats.csv')
    savefig('%sAvgConcStats%s.pdf' % (rand_prefix,db))
    close()

#comulative graph - x axis - avg. prot. conc. (or molecule count per cell), y axis, comulative % out of proteome.
def plotComulativeGraph():
    figure(figsize=(5,3))
    sp = [subplot(121),subplot(122)]

    conds,gr,coli_data = datas['Heinemann']
    avgs = sorted(coli_data['avg'].values)

    sp[0].plot(avgs,cumsum(avgs),'.',markersize=0.5)
    sp[0].set_xlabel('Avg. prot. conc.',fontsize=6)
    sp[0].set_xscale('log')

    sp[1].plot(arange(0,len(avgs)),cumsum(avgs),'.',markersize=0.5)
    sp[1].set_xlabel('num. of prots',fontsize=6)

    for i in range(2):
        sp[i].set_ylabel('accumulated fraction \n out of proteome',fontsize=6)
        sp[i].axhline(xmin=0,xmax=i*2000+1,y=0.05,ls='--',color='black',lw=0.5)
        sp[i].axhline(xmin=0,xmax=i*2000+1,y=0.01,ls='--',color='black',lw=0.5)
        set_ticks(sp[i],6)

    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Cumulative proteome concentration distribution for Heinemann")
    savefig('%sDistStatsHein.pdf' % rand_prefix)
    close()

#plot the graphs for the 10 highest abundance proteins with their descriptions.
def plotHighAbundance():
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(131),'Peebo':subplot(132),'Valgepea':subplot(133)}
    for db in dbs:
        p = ps[db]
        conds,gr,coli_data = datas[db]
        coli_data = coli_data.copy()
        if db == 'Heinemann':
            coli_data['ID']=coli_data['protName']
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
    #fig = gcf()
    #py.plot_mpl(fig,filename="Most abundant proteins concentration vs growth rate")
    savefig('%shighest.pdf' % rand_prefix)
    close()

def plotRibosomal(dbs):
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122),'HuiAlim':subplot(121),'HuiClim':subplot(122),'HuiRlim':subplot(122)}
    for db in dbs:
        p = ps[db]
        conds,gr,coli_data = datas[""][db]
        coli_data = coli_data.copy()
        if db == 'Heinemann':
            coli_data['ID']=coli_data['protName']
        coli_data = coli_data[coli_data['prot']=='Ribosome']
        coli_data_conds = coli_data[conds].copy()
        tot = coli_data_conds.sum()
        means = coli_data_conds.mean(axis=1)
        for col in conds:
            coli_data_conds[col] = coli_data_conds[col]/means
        for i in coli_data_conds.index:
            desc = coli_data.ix[i]
            desc = "%s" % (desc['ID'])
            p.plot(gr.values,coli_data_conds.ix[i].values,label=('%s' % desc))
        tot = tot/tot.mean()
        p.plot(gr.values,tot,'o')
        p.set_ylim(0,3)
    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Ribosomal proteins concentration vs growth")
    savefig('%sribosomal.pdf' % rand_prefix)
    close()


#randomly select a few proteins and plot their prediction vs the actual concentration of a different protein in the HC prots.
def plotPrediction():
    for db in dbs:
        figure(figsize=(5,5))
        conds,gr,coli_data = datas[""][db]
        glob = get_glob(db,coli_data)
        for i in range(1,10):
            p = subplot(330+i)
            samp = random.sample(glob.index,11)
            pred = samp[0:-1]
            est = samp[-1]
            pred = glob.ix[pred]
            est = glob.ix[est]
            pred = pred[conds].sum()
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
            set_ticks(p,8)
            p.set_title("$R^2$=%.2f" % est.corr(linpred)**2,fontsize=8)
        tight_layout()    
        #fig = gcf()
        #py.plot_mpl(fig,filename="Random proteins estimations, 10 proteins at a time, %s" % db)
        savefig('%sRandEstimate%s.pdf' % (rand_prefix,db))
        if rand_prefix == '' and db == 'Heinemann':
            savefig('FigS6.eps')
        close()

#plot ribosomal proteins vs. global cluster proteins with trendlines and R^2 estimates.
def plotRibosomalVsGlobTrend(dbs):
    figure(figsize=(5,3))
    ps = {'Heinemann':subplot(121),'Peebo':subplot(122),'Valgepea':subplot(122), 'HuiAlim':subplot(121),'HuiClim':subplot(122),'HuiRlim':subplot(122)}
    coords = {'Heinemann':0.03,'Peebo':0.03,'Valgepea':0.03,'HuiAlim':0.03,'HuiClim':0.03,'HuiRlim':0.03}
    for db in dbs:
        conds,gr,coli_data = datas[""][db]
        glob = get_glob(db,coli_data)
        no_ribs = glob[glob['prot'] != 'Ribosome']
        ribs = glob[glob['prot'] == 'Ribosome']
        colors = ['blue','green']
        p = ps[db]
        ser = ['Non ribosomal proteins','Ribosomal proteins']
        for j,d in enumerate([no_ribs,ribs]):
            c = colors[j]
            d = d[conds].sum()
            d = d/d.mean()
            p.plot(gr.values,d.values,'o',color=c,label=ser[j])
            alpha,beta,r_val,p_val,std_err = linregress(gr,d)
            p.plot(gr.values,alpha*gr.values+beta,color=c,label="Trend line $R^2$=%.2f" % (gr.corr(d)**2))
            p.set_xlim(xmin=0.)
            p.set_ylim(ymin=0.)
            p.set_xlabel('Growth rate',fontsize=10)
            p.set_ylabel('Normalized concentration',fontsize=10)
        p.legend(loc='lower left', prop={'size':8},numpoints=1)
        set_ticks(p,8)
        text(coords[db],0.93,"data from %s et. al" % db_name[db],fontsize=8,transform=p.transAxes)
    tight_layout()
    #fig = gcf()
    #py.plot_mpl(fig,filename="Ribosomal proteins vs global cluster")
    savefig('%sRibsVsGlob.pdf' % rand_prefix)
    if rand_prefix == '':
        savefig('FigS7.eps')
    close()

def model_effects_plot():
    grs = linspace(0.01,1,15)
    simple = grs
    neg = linspace(-0.4,0.01,6)
    simple = simple/simple.mean()
    degraded = grs+0.4
    degmean = degraded.mean()
    degraded = degraded/degmean
    neg_deg = neg+0.4
    neg_deg = neg_deg/degmean
    rate = 1/(1+0.2/grs)
    rate_effect = grs/rate
    rate_effect = rate_effect/rate_effect.mean()
    figure(figsize=(5,3))
    ax = subplot(111)
    ax.plot(grs,simple,'o',label="Unregulated protein - basic model")
    ax.plot(grs,degraded,'o',label="Unregulated protein - with degradation")
    ax.plot(neg,neg_deg,'--g')
    ax.plot(grs,rate_effect,'o',label="Unregulated protein - under decreasing biosynthesis rate")
    ax.plot(grs,rate,'--r',label="Biosynthesis rate")
    ax.annotate("degradation\nrate", xy=(-0.4,0),xytext=(-0.4,.6),arrowprops=dict(facecolor='black',shrink=0.05,width=1,headwidth=4),horizontalalignment='center',verticalalignment='center',fontsize=8)
    ax.set_xlim(xmin=-0.5)
    ax.set_ylim(ymin=0.)
    ax.set_xlabel('Growth rate [$h^{-1}$]',fontsize=8)
    set_ticks(ax,6)
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_ylabel('Normalized protein concentration',fontsize=8)
    tight_layout()
    subplots_adjust(top=0.83)
    handles,labels=ax.get_legend_handles_labels()
    figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.0,0.8,1,0.2),ncol=2,numpoints=1)
    savefig('TheoreticalModelEffects.pdf')
    savefig('FigS1.eps')
    close()

#plot heinemann/peebo protein correlation
def db_corr():
    print "db corr"
    pts = []
    not_in_peebo = 0
    not_in_hnm = 0
    figure(figsize=(5,5))
    uni_to_b,a,b,b_to_uni = uni_to_locus()
    pbo = datas[""]["Peebo"][2]
    hnm = datas[""]["Heinemann-chemo"][2]
    for i,r in hnm.iterrows():
        x = r["avg"]
        name = r["UP_AC"]
        y=0.1
        if name in uni_to_b:
            b = uni_to_b[name]
            r2 = pbo[pbo['B number identifier'] == b]
            if len(r2)>0:
                y = r2["avg"]
        if x>0:
            if name not in uni_to_b or len(pbo[pbo['B number identifier'] == uni_to_b[name]]) == 0:
                not_in_peebo+=1
        pts.append((x,y))
    for i,r in pbo.iterrows():
        name = r["B number identifier"]
        if r['avg']>0:
            if name not in b_to_uni or len(hnm[hnm['UP_AC'] == b_to_uni[name]]) == 0:
                not_in_hnm+=1
    xs,ys = zip(*pts)
    print ("Not in Peebo %d, not in Schmidt %d" % (not_in_peebo,not_in_hnm))
    plot(xs,ys,'.')
    diag = linspace(1e-5,0.1)
    plot(diag,diag,'r')
    xscale('log')
    yscale('log')
    savefig('dbscorr.pdf')
    close()
        
def plot_corr(p,db1,db2,id1_to_id2):
    print "plot corr %s,%s" % (db1,db2)
    pts = []
    conc1 = datas[""][db1][2]
    conc2 = datas[""][db2][2]
    id1 = id_col_dict[db1]
    id2 = id_col_dict[db2]

    for i,r in conc1.iterrows():
        x = r["gr_cov"]
        if isnan(x):
            continue
        name = r[id1]
        b=''
        if id1 == id2:
            b = name
        else:
            if name in id1_to_id2:
                b = id1_to_id2[name]
        r2 = conc2[conc2[id2] == b]
        if len(r2)>0:
            y = r2["gr_cov"].values[0]
            if not isnan(y):
                pts.append((x,y))
    xs,ys = zip(*pts)
    p.plot(xs,ys,'.',markersize=2)
    diag = linspace(-1,1)
    p.plot(diag,diag,'k--',markersize=0.5)
    p.set_xlabel(db_name[db1] + " " +db_suffix[db1],fontsize = 5)
    p.set_ylabel(db_name[db2] + " " +db_suffix[db2],fontsize = 5)
    p.set_xlim(-1,1)
    p.set_ylim(-1,1)
    p.set_title("cov = %.3f" %(cov(xs,ys)[0][1]),fontsize=6)
    print "covariance of %s and %s, %d proteins:" % (db1,db2,len(xs))
    print cov(xs,ys)
 
def plot_dbs_pairs_corr():
    fig = figure(figsize=(7,5))
    dbs = ['Heinemann-chemo','Peebo-gluc','Valgepea','HuiAlim','HuiClim','HuiRlim']
    dbext = {'Heinemann-chemo':'','Peebo-gluc':'','Valgepea':'','HuiAlim':', A-lim','HuiClim':', C-lim','HuiRlim':', R-lim'}
    (uniprot_to_locus,name_to_locus,name_to_uniprot,locus_to_uniprot) = uni_to_locus()
    for i,(db1,db2) in enumerate(combinations(dbs,2)):
        switch=False
        ps = fig.add_subplot(3,5,1+i)
        id1 = id_col_dict[db1]
        id2 = id_col_dict[db2]
        id1_to_id2 = {}
        if id1 == 'B number identifier' and id2 == 'UP_AC':
            id1_to_id2 = locus_to_uniprot
        if id2 == 'B number identifier' and id1 == 'UP_AC':
            id1_to_id2 = uniprot_to_locus
        if id1 == 'B number identifier' and id2 == 'Gene':
            id1_to_id2 = name_to_locus
            switch = True
        if id2 == 'B number identifier' and id1 == 'Gene':
            id1_to_id2 = name_to_locus
        if id1 == 'UP_AC' and id2 == 'Gene':
            id1_to_id2 = name_to_uniprot
            switch = True
        if id2 == 'UP_AC' and id1 == 'Gene':
            id1_to_id2 = name_to_uniprot
        if switch:
            print db2+db1
            plot_corr(ps,db2,db1,id1_to_id2)
        else:
            print db1+db2
            plot_corr(ps,db1,db2,id1_to_id2)
        set_ticks(ps,3)
    tight_layout()
    savefig('AllDbsCorrelation.pdf')
    savefig('FigS3.eps')
    close()

def plot_all_dbs_hist():
    figure(figsize=(6.5,5))
    dbs = ['Heinemann-chemo','Peebo-gluc','Valgepea','HuiAlim','HuiClim','HuiRlim']
    dbext = {'Heinemann-chemo':'','Peebo-gluc':'','Valgepea':'','HuiAlim':', A-lim','HuiClim':', C-lim','HuiRlim':', R-lim'}
    p=subplot(111)
    for i,db in enumerate(dbs):
        ps = subplot(231+i)
        conds,gr,conc_data = datas[""][db]
        plot_corr_hist(ps,db,conc_data,categories)
        if db == 'Valgepea':
            year = 2013
        else:
            year = 2015
        ps.annotate("data from %s et. al. %d%s" % (db_name[db],year,dbext[db]),xy=(0.5,0.5),xytext=(-0.87,303),fontsize=6,zorder=10)
        ps.set_ylim(0,350)
        ps.set_xlim(-1,1)
        ps.annotate(chr(65+i),xy=(0.5,0.5),xytext=(-0.87,320),fontsize=10,zorder=10)

    #assume both subplots have the same categories.
    handles,labels=ps.get_legend_handles_labels()

    tight_layout()
    figlegend(handles,labels,fontsize=6,mode='expand',loc='upper left',bbox_to_anchor=(0.15,0.8,0.7,0.2),ncol=2)

    subplots_adjust(top=0.9)
    #fig = gcf()
    #py.plot_mpl(fig,filename="Growth rate Correlation histograms")
    savefig('AllDbsGrowthRateCorrelation.pdf')
    savefig('FigS2.eps')
    close()

#init_datasets("")
model_effects_plot()
analyzed_dbs = ['Heinemann','Peebo']
#analyzed_dbs = ['HuiAlim','HuiClim']
special_dbs = ['Heinemann','Peebo','HeinemannLB']
globalResponse = {}
for rand_method in ["simulated","shuffle",""]:
#for rand_method in ["shuffle",""]:
#for rand_method in [""]:
    print "----------------------------------------------------------"
    print rand_method
#for rand_method in ["simulated"]:
    rand_prefix = rand_method
    init_datasets(rand_method)

plot_dbs_pairs_corr()
plot_all_dbs_hist()
db_corr()
print "plotting prediction"
#plotPrediction()        
print "plotting original data graphs"
#tempprotsplot()
corr_andGR_plot('Simulated','Heinemann')
writeTables()
#Single histogram for presentation
    #plotCorrelationHistograms(["Valgepea"],"Val")
 
plotRibosomal(analyzed_dbs)
plotRibosomalVsGlobTrend(analyzed_dbs)
plot_response_hist_graphs(analyzed_dbs)
plotCorrelationHistograms(analyzed_dbs,"")

for rand_method in ["simulated","shuffle",""]:
#for rand_method in ["shuffle",""]:
#for rand_method in [""]:
    plotGlobalResponse(analyzed_dbs,rand_method)
    #plotMultiStats('Valgepea')
    #plotComulativeGraph()
    #plotHighAbundance()
    #variabilityAndGlobClustSlopes(analyzed_dbs)
    variabilityAndGlobClustSlopesNormed(analyzed_dbs,rand_method) #This is the generating function for the variability analysis
    #variablityComparisonHein()
for db in ['Heinemann-chemo','HeinemannLB']:
    corr_andGR_plot(db,'Peebo')
