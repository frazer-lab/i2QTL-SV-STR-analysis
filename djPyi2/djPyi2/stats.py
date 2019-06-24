import pandas as pd
import scipy.stats as sps
import numpy as np

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

stats_r = importr('stats')

# from rpy2.robjects.packages import importr
utils = importr('utils')

def add_bh_fdr(top, col):
    top = top.copy()
    p_vals = top[col].tolist()
    p_adjust = stats_r.p_adjust(FloatVector(p_vals), method = 'fdr')
    top['fdr_corrected_p'] = list(p_adjust)
    return top

def per_variant_vc_unique(df, col1, col2, id_col, overlapping_sets = True):
    """ value counts that are mutually exclusive within first col T/F - second groupby groups-
    optionally, make sure the T/F sets are mutually exclusive entirely"""
    
    in_cat = df[(df[col1] == True)][id_col].unique()
    in_cat_sig = df[(df[col1] == True) & (df[col2] == True)][id_col].unique()
    in_cat_ns = df[(df[col1] == True) & (df[col2] == False)][id_col].unique()
    in_cat_ns = set(in_cat_ns).difference(in_cat_sig)
    num_in_cat_sig = len(in_cat_sig)
    num_in_cat_ns = len(in_cat_ns)
    
    if not overlapping_sets:
        ## make sure the two sets are totally mutually exclusive
        out_cat = df[(df[col1] == False)][id_col].unique()
        # remove things in the category from things out of the category 
        out_cat = set(out_cat).difference(in_cat)
        out_bin = df[(df[col1] == False) & (df[id_col].isin(out_cat))]
    else:
        out_cat = df[(df[col1] == False)][id_col].unique()
        out_bin = df[(df[col1] == False)]
        
    out_cat_sig =  out_bin[(out_bin[col2] == True)][id_col].unique()
    out_cat_ns =  out_bin[(out_bin[col2] == False)][id_col].unique()
    out_cat_ns = set(out_cat_ns).difference(out_cat_sig)
    num_out_cat_sig = len(out_cat_sig)
    num_out_cat_ns = len(out_cat_ns)

    v_in = [[num_in_cat_sig, num_in_cat_ns], [num_out_cat_sig, num_out_cat_ns]]
    return v_in



def vc_to_or(vc, v = False):
    def default_loc(df, a, b, default = 0):
        try:
            out = df.loc[a,b]
            return out
        except:
            return default

        
    if not v:
        
        t_g1 = [default_loc(vc, True, True), default_loc(vc, True, False)]
        f_g1 = [default_loc(vc, False, True), default_loc(vc, False, False)]
        v = [t_g1, f_g1]
    else:
        v = vc

    try:
        odds_ratio, p_fisher = sps.fisher_exact(v, )
    except:
        odds_ratio, p_fisher = (np.NaN, np.NaN)
    
    
    return v, odds_ratio, p_fisher

def gather_odds_ratio_data(df, gb1, gb2, bool_col, gb2_bool = True, unique_col = False,
                           overlapping_sets = False):
    
    gb1_cats = df[gb1].unique()
    if gb2_bool: # if this column is a bool- and not categorical
        data = []
        if not unique_col:
            vc = df.groupby((gb1, gb2))[bool_col].value_counts()
            for c1 in gb1_cats:
                tvc = vc.loc[c1]
                v, odds_ratio, p_fisher = vc_to_or(tvc)
                data.append([c1, gb2, v, odds_ratio, p_fisher])
        else:
            for c1 in gb1_cats:
                tdf = df[df[gb1] == c1]
                if tdf.shape[0] > 0:
                    vc = per_variant_vc_unique(tdf, gb2, 
                                               bool_col, unique_col, 
                                               overlapping_sets= overlapping_sets)
                    v, odds_ratio, p_fisher = vc_to_or(vc, v=True)
                    data.append([c1, gb2, v, odds_ratio, p_fisher]) 
            
        
        df_out = pd.DataFrame(data, columns = [gb1, gb2, 'contingency', 'odds_ratio', 
                                           'p_fisher']).pipe(add_bh_fdr, 'p_fisher')    
    else:
        data = []
        gb2_cats = df[gb2].unique()
        for c2 in gb2_cats:
            df['in_cat'] = (df[gb2] == c2)
            if not unique_col:
                vc = df.groupby((gb1, 'in_cat'))[bool_col].value_counts()     
            for c1 in gb1_cats:
                if unique_col:
                    tdf = df[df[gb1] == c1]
                    if tdf.shape[0] > 0:
                        v_in = per_variant_vc_unique(tdf, 'in_cat', 
                                               bool_col, unique_col, overlapping_sets= overlapping_sets)
#                         print v_in
                        v, odds_ratio, p_fisher = vc_to_or(v_in, v=True)
                        data.append([c1, c2, v, odds_ratio, p_fisher])   
                else:
                    tvc = vc.loc[c1]
                    v, odds_ratio, p_fisher = vc_to_or(tvc)
                    data.append([c1, c2, v, odds_ratio, p_fisher])
                    
        df_out = pd.DataFrame(data, columns = [gb1, gb2, 'contingency', 'odds_ratio', 
                                           'p_fisher']).pipe(add_bh_fdr, 'p_fisher') 
    
    return df_out



def annotate_tests_data(df, col = "significant"):
    df = df.copy()
    def safe_div(a, b):
        try:
            out = a/b
        except:
            out = np.NaN
        return out
            
#     df['frac_non_{}_ol_feat'.format(col)] = df['{}_False'.format(col)].apply(lambda x: safe_div(x[1], x[0]))
#     df['frac_{}_ol_feat'.format(col)] =  df['{}_True'.format(col)].apply(lambda x: safe_div(x[1], x[0]))
    try:
        df['-log10p_fisher'] = np.log10(df['p_fisher']) * -1
    except:
        pass
    
    try:
        df['log_odds_ratio'] = np.log10(df['odds_ratio'])
    except:
        pass
    
    try:
        df['log2_odds_ratio'] = np.log2(df['odds_ratio'])
    except:
        pass
    df = df.reset_index()
    
    df['log2_odds_ratio_raw'] = df['log2_odds_ratio']
    
    
    t_neg_inf = df.log2_odds_ratio == (-np.inf)
    t_pos_inf = (df.log2_odds_ratio == (np.inf))
    
    exclude = t_neg_inf[t_neg_inf].index.tolist() + t_pos_inf[t_pos_inf].index.tolist()
    if len(exclude) > 0:
        inds_non_inf = set(df.index.tolist()).difference(exclude)

        if t_neg_inf[t_neg_inf].shape[0] > 0:
            inds = t_neg_inf[t_neg_inf].index.tolist()
            try:
                m = df.loc[inds_non_inf].log2_odds_ratio.min()
            except:
                m = -1
            
            if m >= -0.5:
                m = -1
            df.loc[inds, 'log2_odds_ratio'] = m
            
        if t_pos_inf[t_pos_inf].shape[0] > 0:
            inds = t_pos_inf[t_pos_inf].index.tolist()
            try:
                m = df.loc[inds_non_inf].log2_odds_ratio.max()
            except:
                m = 2
            if m < 0:
                m = 2
            df.loc[inds, 'log2_odds_ratio'] = m
            
    return df



def compute_mannu_ttest_vs_group(df, col, comp_col = 'SVTYPE_NR_C', comp_group = 'STR'):
    """ Compute the Mann Whitney U and Tttest for difference against some group in long form df
    and Bonferonni Correct P values
    input: dataframe (df) with columnt to groupby (comp_col) that is categorical and a continous column (col)
    to compare
    return: dataframe with stats"""
    
    obs = df[df[comp_col] == comp_group][col].tolist()
    data = []
    for i,df in df.groupby(comp_col):
        if i != comp_group:
            comp = df[df[comp_col] == i][col].tolist()
            mu, p = sps.mannwhitneyu(comp, obs)
            mean_c = np.mean(comp)
            t, p_t = sps.ttest_ind(comp, obs)
            data.append([mu, p, i, mean_c, t, p_t, col])
    stats = pd.DataFrame(data, columns = ['mu', 'p_value', 'SVTYPE_NR_C', 'mean_c', 
                                          'T', 'p_value_t',
                                          "variable"])
    stats['p_bonf'] = stats['p_value'] * (stats.shape[0] -1)
    stats['p_bonf_t'] = stats['p_value_t'] * (stats.shape[0] -1)
    stats = stats.set_index('SVTYPE_NR_C', drop = False)
    return stats





