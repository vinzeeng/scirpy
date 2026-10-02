import squidpy as sq

sdata = sq.datasets.visium_hne_sdata()

print(sdata)
table = sdata.tables["adata"]

print(table)
print()
print(table.obs.head())
print()
print(table.uns)

spots = sdata.shapes["spots"]

print(spots.head())
print(spots.columns)
