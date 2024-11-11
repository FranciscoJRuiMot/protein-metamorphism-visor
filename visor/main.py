import pymol
from pymol import cmd
import tkinter as tk
import pandas as pd
from tkinter import ttk

from csv_extraction.functions import *

class VisorApp:

    def __init__(self, root):
        """
        Inicializa la aplicación principal del visor.

        :param root: La ventana raíz de la aplicación tkinter.
        :type root: tkinter.Tk
        """
        self.root = root
        self.root.title('Visor de Pymol')

        self.initialize_data()
        self.initialize_ui_elements()
        self.initialize_data_trees()

    def initialize_data(self):
        """
        Inicializa y carga los datos iniciales de los clústeres y alineamientos
        desde un archivo CSV y define los índices y variables de estado.
        """
        self.path = 'data_with_annot.csv' #cambiar aquí el fichero de interés a visualizar
        self.data_df = pd.read_csv(self.path)
        self.clusters_id = self.data_df['cluster_id'].unique()
        self.clusters_alignments = self.create_clusters_alignments()
        self.cluster_index = 0
        self.alignment_index = 0

        self.load_initial_alignment_data()

    def load_initial_alignment_data(self):
        """
        Carga los datos de alineamiento específicos para el clúster y alineamiento
        inicial, incluyendo las estructuras y anotaciones.
        """
        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.pdb_1, self.chain_1, self.pdb_2, self.chain_2 = get_estructures(self.alingments, self.alignment_index, self.data_df)
        self.sub1, self.sub2 = get_subclusters_id(self.alingments, self.alignment_index, self.data_df)
        self.annot = get_annotation(self.alingments, self.alignment_index, self.data_df)
        self.annot_metamor = tk.BooleanVar()
        self.annot_metamor.set(self.annot)

    def initialize_ui_elements(self):
        """
        Inicializa y configura todos los elementos de la interfaz gráfica.
        """
        self.setup_navigation_buttons()
        self.setup_labels()
        self.setup_annotation_controls()
        self.setup_save_button()
        #self.create_df_for_trees()
        #self.create_my_tree()
        #self.create_data_trees()

    def setup_navigation_buttons(self):
        """
        Crea los botones de navegación y vincula las teclas para
        la navegación entre subclusters y clusters.
        """
        self.btn_iniciar = tk.Button(root, text='Iniciar', command=self.initiate_pymol)
        self.btn_iniciar.grid(row=0, column=0, columnspan=3, padx=10, pady=10)

        self.btn_next_align = tk.Button(root, text='Siguiente alineamiento', command=self.next_sub)
        self.btn_next_align.grid(row=1, column=0, padx=10, pady=10)
        
        self.btn_prev_align = tk.Button(root, text='Anterior alineamiento', command=self.prev_sub)
        self.btn_prev_align.grid(row=2, column=0, padx=10, pady=10)
        
        self.btn_next_cluster = tk.Button(root, text='Siguiente clúster', command=self.next_cluster)
        self.btn_next_cluster.grid(row=3, column=0, padx=10, pady=10)
        
        self.btn_prev_cluster = tk.Button(root, text='Anterior clúster', command=self.prev_cluster)
        self.btn_prev_cluster.grid(row=4, column=0, padx=10, pady=10)
               
        self.bind_navigation_keys()

    def bind_navigation_keys(self):
        """
        Asigna las teclas de flecha a las funciones de navegación.
        """
        self.root.bind("<Down>", self.next_sub)
        self.root.bind("<Up>", self.prev_sub)
        self.root.bind("<Right>", self.next_cluster)
        self.root.bind("<Left>", self.prev_cluster)

    def setup_labels(self):
        """
        Configura las etiquetas de la interfaz gráfica para mostrar información
        sobre el clúster y los subclusters.
        """
        self.label_cluster = tk.Label(root, text='Clúster:')
        self.label_cluster.grid(row=1, column=1, padx=10, pady=10)

        self.label_subcluster1 = tk.Label(root, text='Sublcúster 1:')
        self.label_subcluster1.grid(row=2, column=1, padx=10, pady=10)

        self.label_subcluster2 = tk.Label(root, text='Sublcúster 2:')
        self.label_subcluster2.grid(row=3, column=1, padx=10, pady=10)

        self.label_annot = tk.Label(root, text='Anotación Polimorfismo')
        self.label_annot.grid(row=1, column=2, padx=10, pady=10)

    def setup_annotation_controls(self):
        """
        Crea los controles de anotación para la interfaz gráfica, incluyendo
        las opciones de metamorfismo.
        """
        self.false_annot_btn = tk.Radiobutton(root, text="No hay metamorfismo", variable=self.annot_metamor, value=False, command=self.metamor_annotation)
        self.false_annot_btn.grid(row=2, column=2, padx=10, pady=10)

        self.true_annot_btn = tk.Radiobutton(root, text="Hay metamorfismo", variable=self.annot_metamor, value=True, command=self.metamor_annotation)
        self.true_annot_btn.grid(row=3, column=2, padx=10, pady=10)

    def metamor_annotation(self):
        """
        Actualiza la anotación de metamorfismo en el DataFrame para el alineamiento actual.

        :return: None
        """      
        alignment = self.alingments[self.alignment_index]
        self.data_df.loc[self.data_df['alignment_result_id'] == alignment, 'metamorphism'] = self.annot_metamor.get()

    def setup_save_button(self):
        """
        Crea el botón para guardar el archivo actualizado.
        """
        self.save_btn = tk.Button(root, text='Guardar fichero', command=self.save_file)
        self.save_btn.grid(row=4, column=2, padx=10, pady=10)

    def save_file(self):
        """
        Guarda el DataFrame actualizado en un archivo CSV en la ruta especificada.

        :return: None
        """
        self.data_df.to_csv(self.path, index=False)

    def initialize_data_trees(self):
        """
        Inicializa las vistas de árbol (TreeView) con los datos de clúster
        para visualizarlos en la interfaz gráfica.
        """
        self.create_df_for_trees()
        self.create_data_trees()

    def create_df_for_trees(self):
        """
        Crea dataframes filtrados por el clúster actual para mostrar
        en las vistas de árbol de la interfaz.
        """
        cluster = self.clusters_id[self.cluster_index]
        self.data_df_cluster = self.data_df.loc[self.data_df['cluster_id'] == cluster].copy() #poner cluster actual
        self.df1 = self.data_df_cluster.loc[:, ['alignment_result_id', 'pdb_id_1', 'chain_1' ,'pdb_id_2', 'chain_2', 'model_1', 'model_2', 'sequence_length_1', 'sequence_length_2', 'ce_rms', 'tm_rms', 'tm_seq_id']].copy()
        self.df1['alignment_result_id'] = self.df1['alignment_result_id'].astype(str)
        self.df2 = self.data_df_cluster.loc[:, ['alignment_result_id', 'tm_score_chain_1', 'tm_score_chain_2' ,'fc_rms', 'fc_identity', 'fc_similarity', 'fc_score', 'fc_align_len', 'cluster_identity_1', 'cluster_identity_2', 'subcluster_identity_1', 'subcluster_identity_2']].copy()
        self.df2['alignment_result_id'] = self.df2['alignment_result_id'].astype(str)

    def create_data_trees(self):
        """
        Crea y configura los widgets TreeView para visualizar los datos en la interfaz.
        """
        self.tree = self.create_treeview(self.root, self.df1, row=0, rowspan=2, column=3)
        self.tree2 = self.create_treeview(self.root, self.df2, row=2, rowspan=2, column=3)
        self.remarked_alignment = str(self.alingments[self.alignment_index])
        self.highlight_row(self.remarked_alignment)

    def create_treeview(self, parent, dataframe, **grid_options):
        """
        Crea y configura un TreeView para un dataframe específico.

        :param parent: Widget padre que contiene el TreeView.
        :param dataframe: Dataframe con los datos a mostrar en el TreeView.
        :param grid_options: Opciones adicionales para grid de tkinter.
        :return: TreeView configurado.
        :rtype: ttk.Treeview
        """
        tree = ttk.Treeview(parent, columns=list(dataframe.columns), show='headings', height=7)
        tree.bind('<Button-1>', self.prevent_selection)
        tree.grid(**grid_options)
        for col in dataframe.columns:
            tree.heading(col, text=col)
            tree.column(col, anchor=tk.CENTER, width=120)
        for _, row in dataframe.iterrows():
            tree.insert('', 'end', values=list(row))
        return tree
    
    def highlight_row(self, name):
        """
        Resalta una fila en las vistas de árbol de datos.

        :param name: Nombre del alineamiento a resaltar.
        :type name: str
        """
        self.highlight_tree_row(self.tree, name)
        self.highlight_tree_row(self.tree2, name)

    def highlight_tree_row(self, tree, name):
        """
        Aplica un color destacado a una fila en un TreeView dado.

        :param tree: TreeView en el que se aplicará el resaltado.
        :type tree: ttk.Treeview
        :param name: Nombre de la fila a resaltar.
        :type name: str
        """
        tree.tag_configure('highlighted', background='lightblue')
        for row in tree.get_children():
            item_values = tree.item(row, 'values')
            tree.item(row, tags=('highlighted',) if item_values[0] == name else ())

    def prevent_selection(self, event):
        """
        Previene la selección de filas en un TreeView.
        
        :param event: Evento de clic.
        :type event: tkinter.Event
        :return: "break" para evitar la selección.
        :rtype: str
        """
        return "break"

    def create_clusters_alignments(self):
        """
        Crea un diccionario con los alineamientos por cada clúster.

        :return: Diccionario de alineamientos por clúster.
        :rtype: dict
        """
        clusters_alignments= {}
        for cluster_id in self.clusters_id:
            alignment_results = self.data_df.loc[self.data_df['cluster_id'] == cluster_id, 'alignment_result_id'].tolist()
            clusters_alignments[cluster_id] = alignment_results

        return clusters_alignments

    def initiate_pymol(self):
        """
        Inicia la interfaz de PyMOL, carga las etiquetas y PDB iniciales.

        :return: None
        """
        pymol.finish_launching()
        self.update_labels()
        self.load_pdb()

    def update_labels(self):
        """
        Actualiza las etiquetas de la interfaz con los datos del clúster y subclúster actuales.

        :return: None
        """
        self.label_cluster.config(text=f'Clúster: {self.clusters_id[self.cluster_index]}')
        self.label_subcluster1.config(text=f'Sublcúster 1: {self.sub1}')
        self.label_subcluster2.config(text=f'Sublcúster 2: {self.sub2}')

    def load_pdb(self):
        """
        Descarga y oculta las estructuras PDB especificadas, y aplica coloración y alineación.

        :return: None
        """
        cmd.fetch(self.pdb_1, async_=0)
        cmd.fetch(self.pdb_2, async_=0)
        cmd.hide('everything')
        self.apply_pdb_colors_and_alignment()

    def apply_pdb_colors_and_alignment(self):
        """
        Aplica coloración y visualización en modo cartoon para las estructuras PDB cargadas. Luego alinea
        las cadenas especificadas.

        :return: None
        """
        cmd.show('cartoon', f'{self.pdb_1} and chain {self.chain_1}')
        cmd.color('red', f'{self.pdb_1} and chain {self.chain_1}')
        cmd.show('cartoon', f'{self.pdb_2} and chain {self.chain_2}')
        cmd.color('blue', f'{self.pdb_2} and chain {self.chain_2}')
        aln = cmd.cealign(f'{self.pdb_1} and chain {self.chain_1}', f'{self.pdb_2} and chain {self.chain_2}')

    def next_sub(self, event=None):
        """
        Cambia al siguiente alineamiento dentro del clúster actual.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        :return: None
        """
        self.change_alignment(next_alignment)

    def prev_sub(self, event=None):
        """
        Cambia al alineamiento anterior dentro del clúster actual.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        :return: None
        """
        self.change_alignment(previous_alignment)

    def next_cluster(self, event=None):
        """
        Cambia al siguiente clúster y actualiza las estructuras y alineamientos correspondientes.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        :return: None
        """
        self.change_cluster(next_cluster)
  
    def prev_cluster(self, event=None):
        """
        Cambia al clúster anterior y actualiza las estructuras y alineamientos correspondientes.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        :return: None
        """
        self.change_cluster(previous_cluster)

    def change_alignment(self, alignment_func):
        """
        Actualiza el índice de alineamiento actual y carga los datos correspondientes en PyMOL.

        :param alignment_func: Función para calcular el nuevo índice de alineamiento.
        :type alignment_func: Callable[[int, List[int]], int]
        :return: None
        """
        self.remove_pymol_data()
        self.alignment_index = alignment_func(self.alignment_index ,self.alingments) #cambiar las funciones para que sean compatibles con dos argumentos
        self.load_data_pymol()
        self.highlight_row(str(self.alingments[self.alignment_index]))

    def change_cluster(self, cluster_func):
        """
        Cambia el índice de clúster actual y carga los alineamientos correspondientes.

        :param cluster_func: Función para calcular el nuevo índice de clúster.
        :type cluster_func: Callable[[int, List[int]], int]
        :return: None
        """
        self.remove_pymol_data()
        self.cluster_index = cluster_func(self.cluster_index, self.clusters_id)
        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.load_data_pymol()
        self.initialize_data_trees()

    def load_data_pymol(self):
        """
        Carga los datos de estructuras PDB y subclústeres para el alineamiento actual en PyMOL y actualiza las etiquetas.

        :return: None
        """

        self.pdb_1, self.chain_1, self.pdb_2, self.chain_2 = get_estructures(self.alingments, self.alignment_index, self.data_df)
        self.sub1, self.sub2 = get_subclusters_id(self.alingments, self.alignment_index, self.data_df)
        print(self.sub1)
        print(self.sub2)

        self.update_labels()
        self.load_pdb()

        self.annot = get_annotation(self.alingments, self.alignment_index, self.data_df)
        self.annot_metamor.set(self.annot)

    def remove_pymol_data(self):
        """
        Elimina las estructuras PDB actuales de PyMOL antes de cargar otras nuevas.

        :return: None
        """
        try:
            cmd.delete(self.pdb_1) #importante poner esto antes de cargar otras estructuras
            cmd.delete(self.pdb_2)
        except:
            pass


root = tk.Tk()
app = VisorApp(root)
root.mainloop()