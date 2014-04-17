import csv
from pandas.io.parsers import read_csv

raw_data = {}

for rate in ['11','21','31','40','48']:
    raw_data[rate]=read_csv('eco_Valgepea_2013_%s.csv' % rate,sep='[\t:]',
                            encoding='iso-8859-1',
                            header = None, names = ['loci','count','conc','weight_cnt',
                                                    'weight_conc','size','name',
                                                    'ko_num','function'],skiprows = 1)

raw_base = raw_data['11'].dropna()

def get_weight_cnd(x,rate):
    df = raw_data[rate]
    row = df[df['loci']==x['loci']].dropna()
    return row['weight_cnt'].values[0]

for rate in ['21','31','40','48']:
    raw_base[rate]=raw_base.apply(lambda x: get_weight_cnd(x,rate),axis=1)
raw_base.rename(columns={'weight_cnt':'11'},inplace=True)
output_data = raw_base[['ko_num','11','21','31','40','48']]
output_data.to_csv('valgepea.csv',index=False)
