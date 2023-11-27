import numpy as np
import pandas as pd
import scipy.io;
import datetime
from datetime import *

def FLX2015_TO_STDNAMES(filename):
    df_fn_2015 = pd.read_csv(filename, header = 0);
    df_fn_2015['Time'] = pd.to_datetime(df_fn_2015['TIMESTAMP_END'], format='%Y%m%d%H%M');
    df_fn_2015 = df_fn_2015.set_index(['Time']);
    df_fn_2015 = df_fn_2015[df_fn_2015.index >= '2001-01-01']
    df_fn_2015['SWC_F_MDS_1']=df_fn_2015['SWC_F_MDS_1']/100.
    df_fn_2015['SWC_F_MDS_2']=df_fn_2015['SWC_F_MDS_2']/100.
    df_fn_2015 = df_fn_2015.rename(columns = {'TA_F':'TA', 'P_F_MDS':'P_mm','TS_F_MDS':'TS','NETRAD':'NETRAD','H_F_MDS':'H','LE_F_MDS':'LE','G_F_MDS':'G','VPD_F_MDS':'VPD', 'APAR_F_MDS':'PAR', 'SWC_F_MDS_1':'SM1','SWC_F_MDS_2':'SM2', 'WS_F':'Wind'})
    good_data = ((df_fn_2015.TA_F_MDS_QC==0)|(df_fn_2015.TA_F_MDS_QC==1))&((df_fn_2015.G_F_MDS_QC==0)|(df_fn_2015.G_F_MDS_QC==1))&((df_fn_2015.H_F_MDS_QC==0)|(df_fn_2015.H_F_MDS_QC==1))&((df_fn_2015.LE_F_MDS_QC==0)|(df_fn_2015.LE_F_MDS_QC==1))
    df_fn_2015[df_fn_2015==-9999]=np.nan;
    AA = np.ones(np.shape(df_fn_2015)); AA[:]=np.nan; AA = df_fn_2015[good_data];
    return df_fn_2015
