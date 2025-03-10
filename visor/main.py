import pymol
from pymol import cmd
import tkinter as tk
import pandas as pd

from csv_extraction.functions import *
from visor_alignment.manager import AlignmentManager

from logic import load_initial_alignment_data, create_clusters_alignments, apply_pdb_colors_and_alignment_path
from logic import load_config
from logic import create_go_graph
from logic import highlight_row
from gui import initialize_ui_elements
from gui import initialize_data_trees
from gui import show_go_graph

from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager
from sqlalchemy.sql import text

from sqlalchemy import update, select
from protein_metamorphisms_is.sql.model.operational.structural_alignment.group import AlignmentGroup
from protein_metamorphisms_is.sql.model.operational.structural_alignment.result import AlignmentResult 

from protein_metamorphisms_is.sql.model.model import (
    AccessionManager,
    PDBExtractor,
    UniProtExtractor,
    Structure3DiManager,
    SequenceClustering,
    StructuralSubClustering,
    SequenceGOAnnotation,
    StructuralAlignmentManager,
    GoMultifunctionalityMetrics,
    SequenceEmbeddingManager
)

class VisorApp:
    """
    Clase principal para la aplicación de visualización en PyMOL.
    
    :param root: Ventana principal de la interfaz gráfica.
    :type root: tkinter.Tk
    :param config_path: Ruta del archivo de configuración.
    :type config_path: str
    """

    def __init__(self, root, config_path):
        """
        Inicializa la aplicación principal del visor.
        """
        self.root = root
        self.root.title('Visor de Pymol')
        self.config = load_config(config_path)
        self.db_manager = DatabaseManager(self.config)
        self.session = self.db_manager.get_session()

        self.initialize_data()
        initialize_ui_elements(self, root)
        initialize_data_trees(self)

    def initialize_data(self):
        """
        Carga los datos iniciales de los alineamientos desde la base de datos.
        """
        sql_file = self.config['consult_path']
        with open(sql_file, 'r') as file:
          sql_query = file.read()
        result = self.session.execute(text(sql_query))
        data = result.fetchall()
        columns = result.keys()
        self.data_df = pd.DataFrame(data, columns=columns)
        print(self.data_df)

        self.clusters_id = self.data_df['cluster_id_1'].unique()
        self.clusters_alignments = create_clusters_alignments(self)
        self.cluster_index = 0
        self.alignment_index = 0

        load_initial_alignment_data(self)
    
    def initiate_pymol(self):
        """
        Inicia la interfaz de PyMOL y carga los archivos cif.
        """
        pymol.finish_launching()
        self.update_labels()
        #Cambios para cargar desde path
        self.load_files_path()

    def create_align_manager(self, alignment_type):
        """
        Crea un gestor de alineamiento para las estructuras.

        :param alignment_type: Tipo de alineamiento a utilizar.
        :type alignment_type: str
        """
        self.remove_pymol_data_path()
        self.manager = AlignmentManager(prot_1=self.file_path_1, name_1=self.name_1, prot_2=self.file_path_2, name_2=self.name_2 ,job_name=str(self.alingments[self.alignment_index]))
        self.manager.call_alignment(alignment_type)

    def bind_navigation_keys(self, *args):
        """
        Asigna las teclas de flecha a las funciones de navegación, alineamiento de estructuras y generación de gráfico GO.
        """
        self.root.bind("<Down>", self.next_sub)
        self.root.bind("<Up>", self.prev_sub)
        self.root.bind("<Right>", self.next_cluster)
        self.root.bind("<Left>", self.prev_cluster)

        self.root.bind("<c>", lambda event: self.create_align_manager('ce_alignment'))
        self.root.bind("<u>", lambda event: self.create_align_manager('US_alignment'))
        self.root.bind("<f>", lambda event: self.create_align_manager('fatcat_alignment'))

        self.root.bind("<g>", self.generate_go_graph)

    def disable_hotkeys(self, *args):
        """
        Deshace la asignación de las teclas para poder escribir comentarios sin generar problemas.
        """
        self.root.unbind("<Down>")
        self.root.unbind("<Up>")
        self.root.unbind("<Right>")
        self.root.unbind("<Left>")

        self.root.unbind("<c>")
        self.root.unbind("<u>")
        self.root.unbind("<f>")
        self.root.unbind("<g>")

    def save_db(self):
        """
        Actualiza los valores de 'metamorphism' y 'comments' en la tabla alignment_group
        usando los datos modificados en el DataFrame.
        """
        try:
            for _, row in self.data_df.iterrows():
                alignment_result_id = row['alignment_result_id']

                # Obtener alignment_group_id desde alignment_result
                alignment_group_id = self.session.scalar(
                    select(AlignmentResult.alignment_group_id).where(AlignmentResult.id == alignment_result_id)
                )

                if alignment_group_id:
                    stmt = (
                        update(AlignmentGroup)
                        .where(AlignmentGroup.id == alignment_group_id)
                        .values(
                            is_metamorphic = row['metamorphism'],
                            comments = row['comments']
                        )
                    )
            
                    self.session.execute(stmt)
            
            self.session.commit()
            print("Base de Datos actualizada")

        except Exception as e:
            print(f"Error al actualizar la BD: {e}")

    def update_labels(self):
        """
        Actualiza las etiquetas de la interfaz con los datos del clúster y subclúster actuales.
        """
        self.label_cluster.config(text=f'Clúster: {self.clusters_id[self.cluster_index]}')
        self.label_subcluster1.config(text=f'Sublcúster 1: {self.sub1}')
        self.label_subcluster2.config(text=f'Sublcúster 2: {self.sub2}')

    def load_files_path(self):
        """
        Carga los ficheros cif en PyMol y nombra las estructuras.
        """
        print(self.file_path_1, self.name_1)
        cmd.load(self.file_path_1, self.name_1)
        print(self.file_path_2, self.name_2)
        cmd.load(self.file_path_2, self.name_2)

        cmd.hide('everything')
        apply_pdb_colors_and_alignment_path(self.name_1, self.name_2)

    def next_sub(self, event=None):
        """
        Cambia al siguiente alineamiento dentro del clúster actual.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        """
        self.change_alignment(next_alignment)

    def prev_sub(self, event=None):
        """
        Cambia al alineamiento anterior dentro del clúster actual.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        """
        self.change_alignment(previous_alignment)

    def next_cluster(self, event=None):
        """
        Cambia al siguiente clúster y actualiza las estructuras y alineamientos correspondientes.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        """
        self.change_cluster(next_cluster)
  
    def prev_cluster(self, event=None):
        """
        Cambia al clúster anterior y actualiza las estructuras y alineamientos correspondientes.

        :param event: Evento opcional para la navegación (como teclas de flecha).
        :type event: tkinter.Event, optional
        """
        self.change_cluster(previous_cluster)

    def change_alignment(self, alignment_func):
        """
        Actualiza el índice de alineamiento actual y carga los datos correspondientes en PyMOL.

        :param alignment_func: Función para calcular el nuevo índice de alineamiento.
        :type alignment_func: Callable[[int, List[int]], int]
        """
        self.remove_pymol_data_path()
        self.alignment_index = alignment_func(self.alignment_index ,self.alingments)
        self.load_data_pymol()
        highlight_row(self, str(self.alingments[self.alignment_index]))

    def change_cluster(self, cluster_func):
        """
        Cambia el índice de clúster actual y carga los alineamientos correspondientes.

        :param cluster_func: Función para calcular el nuevo índice de clúster.
        :type cluster_func: Callable[[int, List[int]], int]
        """
        self.remove_pymol_data_path()
        self.cluster_index = cluster_func(self.cluster_index, self.clusters_id)
        self.alingments, self.alignment_index = get_alingments(self.clusters_id, self.cluster_index, self.clusters_alignments, self.alignment_index)
        self.load_data_pymol()
        initialize_data_trees(self)

    def load_data_pymol(self):
        """
        Carga los datos de estructuras y subclústeres para el alineamiento actual en PyMOL y actualiza las etiquetas.
        """

        self.pdb_1, self.chain_1, self.name_1, self.pdb_2, self.chain_2, self.name_2 = get_structures(self.alingments, self.alignment_index, self.data_df)
        self.sub1, self.sub2 = get_subclusters_id(self.alingments, self.alignment_index, self.data_df)
        self.c_go_term_1_p1, self.c_go_term_2_p1, self.p_go_term_1_p1, self.p_go_term_2_p1, self.f_go_term_1_p1, self.f_go_term_2_p1, self.c_go_term_1_p2, self.c_go_term_2_p2, self.p_go_term_1_p2, self.p_go_term_2_p2, self.f_go_term_1_p2, self.f_go_term_2_p2 = get_go_terms(self.alingments, self.alignment_index, self.data_df)
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
        """
        Elimina las estructuras actuales del PyMol para poder cargar las siguientes.
        """
        try:
            cmd.delete("all")
        except:
            pass

    #Funcionalidad
    def generate_go_graph(self, event=None):
        """
        Genera un gráfico en forma de red que representa y situa la jerarquía de los términos GO de las proteinas que se están visualizando
        """
        go_categories = {'BP': 'p', 'CC' :'c', 'MF':'f'}
        print(f"Elegiste: {self.combo.get()}")
        go_category = go_categories.get(self.combo.get())
        print(f'Que se corresponde con: {go_category}')
        term = f'{go_category}_go_term'

        go_terms_p1 = [getattr(self, term + '_1_p1'), getattr(self, term + '_2_p1')] #poner excepciones de None
        create_go_graph(self.config['obo'], go_terms_p1, "go_path/go_plot_1.png", self.config['color_1'], self.config['color_2'])
        show_go_graph("go_path/go_plot_1.png", "Protein 1 GO Terms")

        go_terms_p2 = [getattr(self, term + '_1_p2'), getattr(self, term + '_2_p2')] #poner excepciones de None
        create_go_graph(self.config['obo'], go_terms_p2, "go_path/go_plot_2.png", self.config['color_1'], self.config['color_2'])
        show_go_graph("go_path/go_plot_2.png", "Protein 2 GO Terms")

def ignore_focus_error(event):
    """
    Previene errores de foco en la interfaz de tkinter.
    """
    try:
        event.widget.focus_set()
    except AttributeError:
        pass

if __name__ == '__main__':
    root = tk.Tk()
    root.bind_all("<Button-1>", lambda event: ignore_focus_error(event))
    app = VisorApp(root, "config/config.yaml")
    root.mainloop()