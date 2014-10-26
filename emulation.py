import numpy
import pylab
import matplotlib
from scipy.stats import gaussian_kde,linregress,variation
import pandas as pd
from scipy import stats
from matplotlib.pyplot import hist, savefig, figure,figlegend,legend,plot,xlim,ylim,xlabel,ylabel,tight_layout,tick_params,subplot,subplots_adjust,text,subplots
from numpy import linspace,ndarray,arange,sum,square,array,cumsum,ones,mean
from numpy.random import randn
import random
import plotly.plotly as py


xs = arange(0.1,0.6,0.1)
slope = 1.6
b=1-1.6*0.3
ys = xs*slope+b
print linregress(xs,ys)
columns = ['1','2','3','4','5']
conc_data = pd.DataFrame(columns=columns)

figure(figsize=(16,12))
def put_alpha(x):
    r = linregress(xs,x.values)
    print x
    print r
    return r[0]

for i in range(0,1000):
    p = [random.gauss(y,y*0.2) for y in ys]
    conc_data.loc[i]=p
conc_data['avg'] = conc_data[columns].mean(axis=1)
conc_data['cv'] = conc_data[columns].apply(lambda x: variation(x), axis=1)
conc_data['alpha'] = conc_data[columns].apply(lambda x: linregress(xs,x/x.mean())[0],axis=1)
conc_data['corr'] = conc_data[columns].apply(lambda x: x.corr(pd.Series(xs,index=columns)),axis=1)
conc_data['r2']=conc_data['corr']**2

print conc_data.mean()
p = subplot(2,2,1)
p.hist(conc_data.alpha,20)
p.set_xlabel('normalized slope')

p1 = subplot(2,2,2)
p1.plot(conc_data['alpha'].values,conc_data['corr'].values,'.',markersize=2)
p1.set_xlabel('normalized slope')
p1.set_ylabel('correlation with GR')
p2 = subplot(223)
p2.hist(conc_data['corr'],bins=arange(-1.1,1.1,0.1))
p2.set_xlim(-1,1)
p2.set_xlabel('correlation with GR')
p3 = subplot(224)
p3.plot(conc_data.alpha.values,conc_data.r2.values,'.',markersize=2)
p3.set_xlabel('normalized slope')
p3.set_ylabel('$R^2$')
savefig('temp.pdf')
