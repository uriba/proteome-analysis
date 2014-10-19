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
from math import sqrt,isnan,log
import random


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

savefig('Noncorrelated.pdf')

