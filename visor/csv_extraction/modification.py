import pandas as pd

'''
Código para crear el csv con la anotación de metamorfismo para importarlo en la app
'''

path = 'mis_datos_07_12.csv'
data_df = pd.read_csv(path)
data_df['metamorphism'] = False
data_df['comments'] = 'Sin comentarios'
data_df.to_csv('mis_datos_07_12_mod2.csv', index=False)