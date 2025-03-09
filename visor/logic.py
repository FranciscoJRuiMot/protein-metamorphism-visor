#Lógica
import tkinter as tk
from csv_extraction.functions import *
from pymol import cmd
import yaml

from goatools.obo_parser import GODag
from goatools.cli.gosubdag_plot import PlotCli
import os

def load_config(file_path):
    """Carga la configuración desde un archivo YAML."""
    with open(file_path, "r") as file:
        return yaml.safe_load(file)

def create_clusters_alignments(app_visor):
    """
    Crea un diccionario con los alineamientos por cada clúster.

    :return: Diccionario de alineamientos por clúster.
    :rtype: dict
    """
    clusters_alignments= {}
    for cluster_id in app_visor.clusters_id:
        alignment_results = app_visor.data_df.loc[app_visor.data_df['cluster_id_1'] == cluster_id, 'alignment_result_id'].tolist()
        clusters_alignments[cluster_id] = alignment_results

    return clusters_alignments

def load_initial_alignment_data(app_visor):
    """
    Carga los datos de alineamiento específicos para el clúster y alineamiento
    inicial, incluyendo las estructuras y anotaciones.
    """
    app_visor.alingments, app_visor.alignment_index = get_alingments(app_visor.clusters_id, app_visor.cluster_index, app_visor.clusters_alignments, app_visor.alignment_index)
    app_visor.pdb_1, app_visor.chain_1, app_visor.name_1, app_visor.pdb_2, app_visor.chain_2, app_visor.name_2 = get_structures(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    ##Cambios para cargar desde path
    app_visor.file_path_1, app_visor.file_path_2 = get_structures_path(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    ##
    app_visor.sub1, app_visor.sub2 = get_subclusters_id(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    app_visor.annot = get_annotation(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    app_visor.annot_metamor = tk.BooleanVar()
    app_visor.annot_metamor.set(app_visor.annot)

def apply_pdb_colors_and_alignment_path(name_1, name_2): #aqui implementar el comando de object para que se vean las distancias (mirar como sería en biopyt)
    cmd.show('cartoon', name_1)
    cmd.color('blue', name_1)
    cmd.show('cartoon', name_2)
    cmd.color('red', name_2)

#funcionalidad
def create_go_graph(obo_file, go_terms, out_file, color_1, color_2):
    godag = GODag(obo_file)

    with open("go_path/go_colors.txt", "w") as color_file:
        color_term_1 = f"{go_terms[0]} {color_1}"
        color_term_2 = f"{go_terms[1]} {color_2}"
        color_file.write(color_term_1 + "\n" + color_term_2)

    kws_plt = {
    "obo": obo_file,
    "GO": go_terms,
    "go_color_file": "go_path/go_colors.txt",
    "outfile": out_file}
        
    plotter = PlotCli()
    fout_img, objplt = plotter.plot(godag, kws_plt)
    os.rename(fout_img, out_file)

#comentarios
def save_comments(app_visor):
    alignment = app_visor.alingments[app_visor.alignment_index]
    app_visor.data_df.loc[app_visor.data_df['alignment_result_id'] == alignment, 'comments'] = app_visor.comment_text.get()

def metamor_annotation(app_visor):
    """
    Actualiza la anotación de metamorfismo en el DataFrame para el alineamiento actual.

    :return: None
    """      
    alignment = app_visor.alingments[app_visor.alignment_index]
    app_visor.data_df.loc[app_visor.data_df['alignment_result_id'] == alignment, 'metamorphism'] = app_visor.annot_metamor.get()

#trees
def create_df_for_trees(app_visor):
    """
    Crea dataframes filtrados por el clúster actual para mostrar
    en las vistas de árbol de la interfaz.
    """
    cluster = app_visor.clusters_id[app_visor.cluster_index]
    app_visor.data_df_cluster = app_visor.data_df.loc[app_visor.data_df['cluster_id_1'] == cluster].copy() #poner cluster actual

    config = app_visor.config
    app_visor.trees_df = {}

    for df_name, df_config in config['trees_dataframes'].items():
        columns = df_config['columns']

        df_filtered = app_visor.data_df_cluster[columns].copy()
        df_filtered['alignment_result_id'] = df_filtered['alignment_result_id'].astype(str)

        app_visor.trees_df[df_name] = df_filtered

def prevent_selection(app_visor, event=None):
    """
    Previene la selección de filas en un TreeView.
        
    :param event: Evento de clic.
    :type event: tkinter.Event
    :return: "break" para evitar la selección.
    :rtype: str
    """
    return "break"

def highlight_row(app_visor, name):
    """
    Resalta una fila en las vistas de árbol de datos.

    :param name: Nombre del alineamiento a resaltar.
    :type name: str
    """
    for tree in app_visor.trees:
        highlight_tree_row(tree, name)

def highlight_tree_row(tree, name):
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