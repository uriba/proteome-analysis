import csv
from pandas.io.parsers import read_csv

for rate in ['11','21','31','40','48']:
    read_csv('eco_Valgepea_2013_%s.csv' % rate,sep='[\t:]',encoding='iso-8859-1',header = None, names = ['loci','count','conc','weight_cnt','weight_conc','size','name','ko_num','function'])
