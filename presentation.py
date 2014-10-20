import csv
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.stats import gaussian_kde,linregress
from scipy import stats
from Bio import SeqIO
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust,text,subplots
from numpy import linspace,ndarray,arange,sum,square,array,cumsum,ones
from numpy.random import randn
from analysis import *
import matplotlib
from matplotlib.pyplot import gcf
from math import sqrt,isnan,log
import random
import plotly.plotly as py


py.sign_in("uri.barenholz", "hvi3ma3m30")

figure(figsize=(5,3))

p=subplot(111)
xs = linspace(0.1,0.8,12)
params = [(0.5,0.8),(1.5,0.5),(3,-0.3)]
xlim(0,1)
xlabel("Growth rate")
ylabel("Protein concentration")
for i,(a,b) in enumerate(params):
    plot(xs,a*xs+b,'o',label='Protein %d' % (i+1))

legend(loc=1, prop={'size':6},numpoints=1)
tight_layout()
fig = gcf()
#print py.plot_mpl(fig,filename="presentation1")
savefig('Noncorrelated.pdf')

figure(figsize=(6,3))
p = subplot(1,2,1)
labels = ['Education','Health','Social\nwelfare','Infrastructure']
fracs = [35,15,20,30]
patches, texts, autotexts = p.pie(fracs,labels=labels,autopct='%1.1f%%',radius=0.8)
for t in texts:
    t.set_size('smaller')
tight_layout()
savefig('SimpleBudget.pdf')

figure(figsize=(6,3))
p = subplot(1,2,1)
labels = ['Education','Health','Social\nwelfare','Infrastructure']
fracs = [35,15,20,30]
patches, texts, autotexts = p.pie(fracs,labels=labels,autopct='%1.1f%%',radius=0.8)
for t in texts:
    t.set_size('smaller')

p = subplot(1,2,2)
labels = ['Education','Health','Social\nwelfare','Infrastructure','Defence']
fracs = [35,15,20,30,25]
explode = [0,0,0,0,0.05]
patches, texts, autotexts = p.pie(fracs,labels=labels,autopct='%1.1f%%',radius=0.8,explode=explode)
for t in texts:
    t.set_size('smaller')
tight_layout()
savefig('SimpleBudget2.pdf')
