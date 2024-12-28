#User interface
from csv_extraction.functions import *
import tkinter as tk
from tkinter import StringVar

def setup_navigation_buttons(app_visor, root):
    """
    Crea los botones de navegación y vincula las teclas para
    la navegación entre subclusters y clusters.
    """
    app_visor.btn_iniciar = tk.Button(root, text='Iniciar', command=app_visor.initiate_pymol)
    app_visor.btn_iniciar.grid(row=0, column=0, columnspan=3, padx=10, pady=10)

    app_visor.btn_next_align = tk.Button(root, text='Siguiente alineamiento', command=app_visor.next_sub)
    app_visor.btn_next_align.grid(row=1, column=0, padx=10, pady=10)
        
    app_visor.btn_prev_align = tk.Button(root, text='Anterior alineamiento', command=app_visor.prev_sub)
    app_visor.btn_prev_align.grid(row=2, column=0, padx=10, pady=10)
        
    app_visor.btn_next_cluster = tk.Button(root, text='Siguiente clúster', command=app_visor.next_cluster)
    app_visor.btn_next_cluster.grid(row=3, column=0, padx=10, pady=10)
        
    app_visor.btn_prev_cluster = tk.Button(root, text='Anterior clúster', command=app_visor.prev_cluster)
    app_visor.btn_prev_cluster.grid(row=4, column=0, padx=10, pady=10)

    ##Nuevos botones alineamientos
    app_visor.btn_ce_align = tk.Button(root, text='CE-Align', command=lambda: app_visor.create_align_manager('ce_alignment'))#faltan que se borren las estructuras
    app_visor.btn_ce_align.grid(row=5, column=0, padx=10, pady=10)

    app_visor.btn_us_align = tk.Button(root, text='US-Align', command=lambda: app_visor.create_align_manager('US_alignment'))
    app_visor.btn_us_align.grid(row=5, column=1, padx=10, pady=10)

    app_visor.btn_fatcat_align = tk.Button(root, text='FATCAT-Align', command=lambda: app_visor.create_align_manager('fatcat_alignment'))
    app_visor.btn_fatcat_align.grid(row=5, column=2, padx=10, pady=10)
               
    app_visor.bind_navigation_keys()

def setup_labels(app_visor, root):
    """
    Configura las etiquetas de la interfaz gráfica para mostrar información
    sobre el clúster y los subclusters.
    """
    app_visor.label_cluster = tk.Label(root, text='Clúster:')
    app_visor.label_cluster.grid(row=1, column=1, padx=10, pady=10)

    app_visor.label_subcluster1 = tk.Label(root, text='Sublcúster 1:')
    app_visor.label_subcluster1.grid(row=2, column=1, padx=10, pady=10)

    app_visor.label_subcluster2 = tk.Label(root, text='Sublcúster 2:')
    app_visor.label_subcluster2.grid(row=3, column=1, padx=10, pady=10)

    app_visor.label_annot = tk.Label(root, text='Anotación Polimorfismo')
    app_visor.label_annot.grid(row=1, column=2, padx=10, pady=10)

def setup_annotation_controls(app_visor, root):
    """
    Crea los controles de anotación para la interfaz gráfica, incluyendo
    las opciones de metamorfismo.
    """
    app_visor.false_annot_btn = tk.Radiobutton(root, text="No hay metamorfismo", variable=app_visor.annot_metamor, value=False, command=app_visor.metamor_annotation)
    app_visor.false_annot_btn.grid(row=2, column=2, padx=10, pady=10)

    app_visor.true_annot_btn = tk.Radiobutton(root, text="Hay metamorfismo", variable=app_visor.annot_metamor, value=True, command=app_visor.metamor_annotation)
    app_visor.true_annot_btn.grid(row=3, column=2, padx=10, pady=10)

def setup_save_button(app_visor, root):
    """
    Crea el botón para guardar el archivo actualizado.
    """
    app_visor.save_btn = tk.Button(root, text='Guardar fichero', command=app_visor.save_file)
    app_visor.save_btn.grid(row=4, column=2, padx=10, pady=10)

def setup_entry_comments(app_visor, root):
    app_visor.comment = get_comments(app_visor.alingments, app_visor.alignment_index, app_visor.data_df)
    app_visor.comment_text = StringVar()
    app_visor.comment_text.set(app_visor.comment) #añadir que se vaya cambiando en cada movimiento
        
    app_visor.entry_comment = tk.Entry(root, textvariable=app_visor.comment_text)
    app_visor.entry_comment.grid(row=6, column=3, padx=10, pady=10, sticky="ew")
    app_visor.entry_comment.bind("<FocusIn>", app_visor.disable_hotkeys)
    app_visor.entry_comment.bind("<FocusOut>", app_visor.bind_navigation_keys)

    app_visor.save_comment_btn = tk.Button(root, text='Guardar comentarios actuales', command=app_visor.save_comments)
    app_visor.save_comment_btn.grid(row=6, column=0, padx=10, pady=10, columnspan=2)

    app_visor.comment_label = tk.Label(root ,text='Comentarios:')
    app_visor.comment_label.grid(row=6, column=2)

def initialize_ui_elements(app_visor, root):
    """
    Inicializa y configura todos los elementos de la interfaz gráfica.
    """
    setup_navigation_buttons(app_visor, root)
    setup_labels(app_visor, root)
    setup_annotation_controls(app_visor, root)
    setup_save_button(app_visor, root)
    setup_entry_comments(app_visor, root)

##Trees

def create_df_for_trees(app_visor):
    """
    Crea dataframes filtrados por el clúster actual para mostrar
    en las vistas de árbol de la interfaz.
    """
    cluster = app_visor.clusters_id[app_visor.cluster_index]
    app_visor.data_df_cluster = app_visor.data_df.loc[app_visor.data_df['cluster_id_1'] == cluster].copy() #poner cluster actual
    app_visor.df1 = app_visor.data_df_cluster.loc[:, ['alignment_result_id', 'pdb_id_1', 'chain_1' ,'pdb_id_2', 'chain_2', 'model_1', 'model_2', 'sequence_length_1', 'sequence_length_2', 'ce_rms', 'tm_rms', 'tm_seq_id']].copy()
    app_visor.df1['alignment_result_id'] = app_visor.df1['alignment_result_id'].astype(str)
    app_visor.df2 = app_visor.data_df_cluster.loc[:, ['alignment_result_id', 'tm_score_chain_1', 'tm_score_chain_2' ,'fc_rms', 'fc_identity', 'fc_similarity', 'fc_score', 'fc_align_len']].copy()
    app_visor.df2['alignment_result_id'] = app_visor.df2['alignment_result_id'].astype(str)

def create_data_trees(app_visor):
        """
        Crea y configura los widgets TreeView para visualizar los datos en la interfaz.
        """
        app_visor.tree = app_visor.create_treeview(app_visor.root, app_visor.df1, row=0, rowspan=2, column=3)
        app_visor.tree2 = app_visor.create_treeview(app_visor.root, app_visor.df2, row=2, rowspan=2, column=3)
        app_visor.remarked_alignment = str(app_visor.alingments[app_visor.alignment_index])
        app_visor.highlight_row(app_visor.remarked_alignment)