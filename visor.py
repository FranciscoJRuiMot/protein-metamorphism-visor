import pymol
from pymol import cmd
import tkinter as tk
import pandas as pd

from function_csv import *

class VisorApp:

    def __init__(self, root):
        #Iniciar ventana princial
        self.root = root
        self.root.title('Visor de Pymol')

        #Análisis del csv
        self.path = 'example.csv' #cambiar aquí el fichero de interés a visualizar
        self.data_df = pd.read_csv(self.path)
        self.clusters_id = self.data_df['cluster_id'].unique()
        self.clusters_alignments = self.creat_clusters_alignments()

        #Initial values
        self.cluster_index = 0
        self.alignment_index = 0

        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.pdb_1, self.chain_1, self.pdb_2, self.chain_2 = get_estructures(self.alingments, self.alignment_index, self.data_df)
        self.sub1, self.sub2 = get_subclusters_id(self.alingments, self.alignment_index, self.data_df)
        self.annot = get_annotation(self.alingments, self.alignment_index, self.data_df)
        self.annot_metamor = tk.BooleanVar()
        self.annot_metamor.set(self.annot)

        #App Widgets
        self.btn_iniciar = tk.Button(root, text='Iniciar', command=self.iniciar_pymol)
        self.btn_iniciar.grid(row=0, column=0, columnspan=3, padx=10, pady=10)

        self.label_cluster = tk.Label(root, text='Clúster:')
        self.label_cluster.grid(row=1, column=1, padx=10, pady=10)

        self.label_subcluster1 = tk.Label(root, text='Sublcúster 1:')
        self.label_subcluster1.grid(row=2, column=1, padx=10, pady=10)

        self.label_subcluster2 = tk.Label(root, text='Sublcúster 2:')
        self.label_subcluster2.grid(row=3, column=1, padx=10, pady=10)

        self.btn_next_align = tk.Button(root, text='Siguiente alineamiento', command=self.next_sub)
        self.btn_next_align.grid(row=1, column=0, padx=10, pady=10)
        self.root.bind("<Right>", self.next_sub)

        self.btn_prev_align = tk.Button(root, text='Anterior alineamiento', command=self.prev_sub)
        self.btn_prev_align.grid(row=2, column=0, padx=10, pady=10)
        self.root.bind("<Left>", self.prev_sub)

        self.btn_next_cluster = tk.Button(root, text='Siguiente clúster', command=self.next_cluster)
        self.btn_next_cluster.grid(row=3, column=0, padx=10, pady=10)
        self.root.bind("<Down>", self.next_cluster)

        self.btn_prev_cluster = tk.Button(root, text='Anterior clúster', command=self.prev_cluster)
        self.btn_prev_cluster.grid(row=4, column=0, padx=10, pady=10)
        self.root.bind("<Up>", self.prev_cluster)

        self.label_annot = tk.Label(root, text='Anotación Polimorfismo')
        self.label_annot.grid(row=1, column=2, padx=10, pady=10)

        self.false_annot_btn = tk.Radiobutton(root, text="No hay metamorfismo", variable=self.annot_metamor, value=False, command=self.metamor_annotation)
        self.false_annot_btn.grid(row=2, column=2, padx=10, pady=10)

        self.true_annot_btn = tk.Radiobutton(root, text="Hay metamorfismo", variable=self.annot_metamor, value=True, command=self.metamor_annotation)
        self.true_annot_btn.grid(row=3, column=2, padx=10, pady=10)

        self.save_btn = tk.Button(root, text='Guardar fichero', command=self.save_file)
        self.save_btn.grid(row=4, column=2, padx=10, pady=10)

    def creat_clusters_alignments(self):
        clusters_alignments= {}
        # Iterar sobre cada valor único de 'cluster_id'
        for cluster_id in self.clusters_id:
            # Filtrar 'alignment_result_id' para el 'cluster_id' actual
            alignment_results = self.data_df.loc[self.data_df['cluster_id'] == cluster_id, 'alignment_result_id'].tolist()
            # Asignar la lista de resultados al diccionario
            clusters_alignments[cluster_id] = alignment_results

        return clusters_alignments

    def iniciar_pymol(self):
        pymol.finish_launching()

        self.act_label()
        self.cargar_pdb()

    def act_label(self):
        self.label_cluster.config(text=f'Clúster: {self.clusters_id[self.cluster_index]}')
        self.label_subcluster1.config(text=f'Sublcúster 1: {self.sub1}')
        self.label_subcluster2.config(text=f'Sublcúster 2: {self.sub2}')

    def cargar_pdb(self):
        cmd.fetch(self.pdb_1, async_=0)
        cmd.fetch(self.pdb_2, async_=0)
        # Enfocar en las cadenas de interés
        cmd.hide('everything')
        cmd.show('cartoon', f'{self.pdb_1} and chain {self.chain_1}')
        cmd.color('red', f'{self.pdb_1} and chain {self.chain_1}')
        cmd.show('cartoon', f'{self.pdb_2} and chain {self.chain_2}')
        cmd.color('blue', f'{self.pdb_2} and chain {self.chain_2}')

        # Alinear las estructuras y mostrar la distancia
        aln = cmd.cealign(f'{self.pdb_1} and chain {self.chain_1}', f'{self.pdb_2} and chain {self.chain_2}')

    def next_sub(self, event=None):

        self.del_data_pymol()
        
        self.alignment_index = next_alignment(self.alignment_index, self.alingments)
        print(f'El indice es {self.alignment_index}')
        self.load_data_pymol()

    def prev_sub(self, event=None):

        self.del_data_pymol()

        self.alignment_index = previous_alignment(self.alignment_index)
        print(f'El indice es {self.alignment_index}')
        self.load_data_pymol()

    def next_cluster(self, event=None):

        self.del_data_pymol()

        self.cluster_index = next_cluster(self.cluster_index, self.clusters_id)
        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.load_data_pymol()
    
    def prev_cluster(self, event=None):

        self.del_data_pymol()

        self.cluster_index = previous_cluster(self.cluster_index)
        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.load_data_pymol()

    def load_data_pymol(self):

        self.pdb_1, self.chain_1, self.pdb_2, self.chain_2 = get_estructures(self.alingments, self.alignment_index, self.data_df)
        self.sub1, self.sub2 = get_subclusters_id(self.alingments, self.alignment_index, self.data_df)
        print(self.sub1)
        print(self.sub2)

        self.act_label()
        self.cargar_pdb()

        self.annot = get_annotation(self.alingments, self.alignment_index, self.data_df)
        self.annot_metamor.set(self.annot)

    def del_data_pymol(self):
        try:
            cmd.delete(self.pdb_1) #importante poner esto antes de cargar otras estructuras
            cmd.delete(self.pdb_2)
        except:
            pass

    def metamor_annotation(self):       
        alignment = self.alingments[self.alignment_index]
        self.data_df.loc[self.data_df['alignment_result_id'] == alignment, 'metamorphism'] = self.annot_metamor.get()

    def save_file(self):
        self.data_df.to_csv(self.path, index=False)

root = tk.Tk()
app = VisorApp(root)
root.mainloop()