import pandas as pd
import numpy as np

vtx_file = "Meshes/yf17.xyz"
tri_file = "Meshes/yf17.tri"
tet_file = "Meshes/yf17.tet"

my_vtx = pd.read_csv(vtx_file, delim_whitespace=True, header=None,names=['x','y','z'])

lim = 20
my_vtx.loc[(abs(my_vtx['x']) > lim) | (abs(my_vtx['y']) > lim) | (abs(my_vtx['z']) > lim), 'in_range'] = 0
my_vtx['in_range'].fillna(1,inplace=True)

new_vtx=my_vtx[my_vtx['in_range'] ==1].reset_index()
# print(my_vtx.head(15))
my_vtx = None

new_vtx['index'] = new_vtx['index'].apply(lambda x: x + 1)
print(new_vtx.head(15))

pre_dict = new_vtx[['index']].to_dict()['index']
vtx_dict = {v : k+1 for k, v in pre_dict.items()}

def convert_db(df_file,ncols):
    print("\tImporting data frame from ",df_file)
    column_names=[]
    for i in range(0,ncols):
        column_names.append("n"+str(i+1))
    my_df = pd.read_csv(df_file, delim_whitespace=True, header=None,names=column_names)
#old
#    my_df.loc[(my_df['n1'].isin(new_vtx['index'])) & (my_df['n2'].isin(new_vtx['index'])) & (my_df['n3'].isin(new_vtx['index'])), 'is_valid'] = 1
#new
    cond = True
    for column in column_names:
        cond = cond & (my_df[column].isin(new_vtx['index']))    
    my_df.loc[cond, 'is_valid'] = 1    
    my_df['is_valid'].fillna(0,inplace=True)

    
    new_df = my_df[my_df['is_valid'] == 1].reset_index().drop(['index','is_valid'],axis=1)
    my_df = None
    print("\tRemapping vertices...")
    new_df = new_df.applymap(vtx_dict.get)
    new_df = new_df.astype(int)
    return new_df
print("Working on triangle files....")
my_tri=convert_db(tri_file,3)
print("Done, example:\n", my_tri.head(10))
print("Working on tetra files....")
my_tet=convert_db(tet_file,4)
print("Done, example:\n",my_tet.head(10))

new_vtx=new_vtx.drop(['index','in_range'],axis=1)
print("Exporting new csv files...")
new_vtx.to_csv("yf17.xyz",sep=" ",header=False,index=False)
my_tri.to_csv("yf17.tri",sep=" ",header=False,index=False)
my_tet.to_csv("yf17.tet",sep=" ",header=False,index=False)
