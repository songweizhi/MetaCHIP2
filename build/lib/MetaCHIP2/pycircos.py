
def circos_hgt():
    pass


'''

Rscript circos_HGT.R -m %s -p %s

'''
from pycirclize import Circos
from pycirclize.parser import Matrix
import pandas as pd

# Create from-to table dataframe & convert to matrix
fromto_table_df = pd.DataFrame(
    [
        ["A", "B", 10],
        ["A", "C", 5],
        ["A", "D", 15],
        ["A", "E", 20],
        ["A", "F", 3],
        ["B", "A", 3],
        ["B", "G", 15],
        ["F", "D", 13],
        ["F", "E", 2],
        ["E", "A", 20],
        ["E", "D", 6],
    ],
    columns=["from", "to", "value"],  # Column name is optional
)
matrix = Matrix.parse_fromto_table(fromto_table_df)

circos = Circos.initialize_from_matrix(
    matrix,
    space=3,
    cmap="viridis",
    ticks_interval=5,
    label_kws=dict(size=12, r=110),
    link_kws=dict(direction=1, ec="black", lw=0.5),
)

print(fromto_table_df.to_string(index=False))
fig = circos.plotfig()
