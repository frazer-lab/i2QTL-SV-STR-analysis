import pandas as pd



def concat_cols(df, cols, sep='_', formatting = False, name=False, inplace=False):
    if not name:
        name = sep.join(cols)
    if not inplace:
        df = df.copy()
        
    col_lists = zip(*[df[c].astype(str).tolist() for c in cols])
    
    if not formatting:
        concat = [sep.join(i) for i in col_lists]
    else:
        concat = [formatting.format(*i) for i in col_lists]
        
    df[name] = concat
    return df

def expand_col_split_join(df, col, sep = '_', col_names = False, subset = False, dtypes = False, overwrite=True):
    df = df.copy()
    tdf = df[col].str.split(sep, expand = True)
    if col_names:
        tdf.columns = col_names
        
    if subset:
        tdf = tdf[subset]
        col_names = subset
    if dtypes:
        for c, dt in zip(col_names, dtypes):
            tdf[c] = tdf[c].astype(dt)
        
    if overwrite:
        # if we don't care about duplicating columns and want to avoid throwing errors
        for c in col_names:
            df[c] = tdf[c]
    else:
        # if we wanna worry about duplicate columns
        df = df.join(tdf)
    return df


def limit_column_width(df, num_px = 100):   
    return df.style.set_table_styles([dict(selector="th",props=[('max-width', '{}px'.format(str(num_px)))])])

def vc_pivot_w_proportion(df, gb, col_unique, name_count_col = 'count', 
                         name_fraction_col = 'fraction', fraction_group = True):
    vc = df.groupby(gb)[col_unique].value_counts().to_frame(name_count_col).reset_index()
    vc[col_unique] =  vc[col_unique].apply(lambda x: '{}_{}'.format(col_unique, str(x)))
    vc = pd.pivot_table(vc, index = gb, columns = col_unique, values = name_count_col).fillna(0)
    vcfrac = (df.groupby(gb)[col_unique]
              .value_counts(normalize = True)
              .to_frame(name_fraction_col)
              .reset_index()
              .groupby(col_unique).get_group(fraction_group)
              .set_index(gb))

    vc = vc.join(vcfrac[name_fraction_col])
    vc[name_fraction_col] = vc[name_fraction_col].fillna(0)
    return vc



