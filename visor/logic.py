#Lógica
import tkinter as tk
from csv_extraction.functions import *

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
    app_visor.pdb_1, app_visor.chain_1, app_visor.pdb_2, app_visor.chain_2 = get_estructures(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    ##Cambios para cargar desde path
    app_visor.file_path_1, app_visor.file_path_2 = get_structures_path(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    ##
    app_visor.sub1, app_visor.sub2 = get_subclusters_id(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    app_visor.annot = get_annotation(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    app_visor.annot_metamor = tk.BooleanVar()
    app_visor.annot_metamor.set(app_visor.annot)