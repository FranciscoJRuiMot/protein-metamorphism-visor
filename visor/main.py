import pymol
from pymol import cmd
import tkinter as tk
import pandas as pd
from tkinter import ttk
from tkinter import StringVar

from csv_extraction.functions import *
from visor_alignment.manager import AlignmentManager

from logic import load_initial_alignment_data, create_clusters_alignments
from gui import initialize_ui_elements
from gui import create_df_for_trees, create_data_trees

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
        #self.create_align_manager()
        #self.initialize_ui_elements()
        initialize_ui_elements(self, root)
        self.initialize_data_trees()

    def initialize_data(self):
        """
        Inicializa y carga los datos iniciales de los clústeres y alineamientos
        desde un archivo CSV y define los índices y variables de estado.
        """
        self.path = '19feb.csv' #cambiar aquí el fichero de interés a visualizar
        self.data_df = pd.read_csv(self.path)
        self.data_df['metamorphism'] = False
        self.data_df['comments'] = 'Sin comentarios'
        print(self.data_df)
        self.clusters_id = self.data_df['cluster_id_1'].unique()
        self.clusters_alignments = create_clusters_alignments(self)
        self.cluster_index = 0
        self.alignment_index = 0

        load_initial_alignment_data(self)

    def save_comments(self):

        alignment = self.alingments[self.alignment_index]
        self.data_df.loc[self.data_df['alignment_result_id'] == alignment, 'comments'] = self.comment_text.get()

    def create_align_manager(self, alignment_type):
        self.remove_pymol_data_path()

        name_struct_1 = self.pdb_1 + '_' + self.chain_1
        name_struct_2 = self.pdb_2 + '_' + self.chain_2

        self.manager = AlignmentManager(prot_1=self.file_path_1, name_1=name_struct_1, prot_2=self.file_path_2, name_2=name_struct_2 ,job_name=str(self.alingments[self.alignment_index]))
        self.manager.call_alignment(alignment_type)

    def bind_navigation_keys(self, *args):
        """
        Asigna las teclas de flecha a las funciones de navegación.
        """
        self.root.bind("<Down>", self.next_sub)
        self.root.bind("<Up>", self.prev_sub)
        self.root.bind("<Right>", self.next_cluster)
        self.root.bind("<Left>", self.prev_cluster)

        self.root.bind("<c>", lambda event: self.create_align_manager('ce_alignment'))
        self.root.bind("<u>", lambda event: self.create_align_manager('US_alignment'))
        self.root.bind("<f>", lambda event: self.create_align_manager('fatcat_alignment'))

    def disable_hotkeys(self, *args):
        self.root.unbind("<Down>")
        self.root.unbind("<Up>")
        self.root.unbind("<Right>")
        self.root.unbind("<Left>")

        self.root.unbind("<c>")
        self.root.unbind("<u>")
        self.root.unbind("<f>")

    def metamor_annotation(self):
        """
        Actualiza la anotación de metamorfismo en el DataFrame para el alineamiento actual.

        :return: None
        """      
        alignment = self.alingments[self.alignment_index]
        self.data_df.loc[self.data_df['alignment_result_id'] == alignment, 'metamorphism'] = self.annot_metamor.get()

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
        create_df_for_trees(self)
        create_data_trees(self)

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

    def initiate_pymol(self):
        """
        Inicia la interfaz de PyMOL, carga las etiquetas y PDB iniciales.

        :return: None
        """
        pymol.finish_launching()
        self.update_labels()
        #Cambios para cargar desde path
        self.load_files_path()

    def update_labels(self):
        """
        Actualiza las etiquetas de la interfaz con los datos del clúster y subclúster actuales.

        :return: None
        """
        self.label_cluster.config(text=f'Clúster: {self.clusters_id[self.cluster_index]}')
        self.label_subcluster1.config(text=f'Sublcúster 1: {self.sub1}')
        self.label_subcluster2.config(text=f'Sublcúster 2: {self.sub2}')

    def load_files_path(self):
        print(self.file_path_1, self.pdb_1 + '_' + self.chain_1)
        cmd.load(self.file_path_1, self.pdb_1 + '_' + self.chain_1 + '_' + self.model_1) #cambiar a f string
        print(self.file_path_2, self.pdb_2 + '_' + self.chain_2)
        cmd.load(self.file_path_2, self.pdb_2 + '_' + self.chain_2 + '_' + self.model_2)

        cmd.hide('everything')
        self.apply_pdb_colors_and_alignment_path()

    def apply_pdb_colors_and_alignment_path(self): #aqui implementar el comando de object para que se vean las distancias (mirar como sería en biopyt)
        cmd.show('cartoon', self.pdb_1 + '_' + self.chain_1 + '_' + self.model_1)
        cmd.color('red', self.pdb_1 + '_' + self.chain_1 + '_' + self.model_1)
        cmd.show('cartoon', self.pdb_2 + '_' + self.chain_2 + '_' + self.model_2)
        cmd.color('blue', self.pdb_2 + '_' + self.chain_2 + '_' + self.model_2)

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
        self.remove_pymol_data_path()
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
        self.remove_pymol_data_path()
        self.cluster_index = cluster_func(self.cluster_index, self.clusters_id)
        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.load_data_pymol()
        self.initialize_data_trees()

    def load_data_pymol(self):
        """
        Carga los datos de estructuras PDB y subclústeres para el alineamiento actual en PyMOL y actualiza las etiquetas.

        :return: None
        """

        self.pdb_1, self.chain_1, self.pdb_2, self.chain_2, self.model_1, self.model_2 = get_estructures(self.alingments, self.alignment_index, self.data_df)
        self.sub1, self.sub2 = get_subclusters_id(self.alingments, self.alignment_index, self.data_df)
        ##
        self.file_path_1, self.file_path_2 = get_structures_path(self.alingments, self.alignment_index, self.data_df)
        ##
        print(self.sub1)
        print(self.sub2)

        self.update_labels()
        #cambiado para path
        self.load_files_path()

        self.annot = get_annotation(self.alingments, self.alignment_index, self.data_df)
        self.annot_metamor.set(self.annot)

        #Comentarios
        self.comment = get_comments(self.alingments, self.alignment_index, self.data_df)
        self.comment_text.set(self.comment)

    def remove_pymol_data_path(self):

        try:
            cmd.delete(self.pdb_1 + '_' + self.chain_1 + '_' + self.model_1) #importante poner esto antes de cargar otras estructuras
            cmd.delete(self.pdb_2 + '_' + self.chain_2 + '_' + self.model_2)
            cmd.delete("all")
        except:
            pass


if __name__ == '__main__':
    root = tk.Tk()
    root.bind_all("<Button-1>", lambda event: event.widget.focus_set())
    app = VisorApp(root)
    root.mainloop()