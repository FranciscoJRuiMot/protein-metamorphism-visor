import pandas as pd

'''
Código para crear el csv con la anotación de metamorfismo para importarlo en la app
'''

path = 'metamorphics_28_oct.csv'
data_df = pd.read_csv(path)
data_df['metamorphism'] = False
data_df.to_csv('data_with_annot.csv', index=False)