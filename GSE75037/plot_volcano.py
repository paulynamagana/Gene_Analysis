from bioinfokit import analys, visuz


df = analys.get_data("FULL_RESULTS").data
df.head(2)

visuz.GeneExpression.volcano(df=df, lfc='logFC', pv='P.Value')
