'''
Código para crear el csv con la anotación de metamorfismo para importarlo en la app

path = 'metamorphics_28_oct.csv'
data_df = pd.read_csv(path)
data_df['metamorphism'] = False
print(data_df)
data_df.to_csv('data_with_annot.csv', index=False)
'''

##Funciones para la extracción de datos del csv
def get_alingments(clusters_id, cluster_index, clusters_alignments, alignment_index):
    cluster_id = int(clusters_id[cluster_index])
    print(cluster_id)

    alingments = clusters_alignments.get(cluster_id, [])
    #print(alingments)
    alignment_index = 0

    return alingments, alignment_index

def get_estructures(alingments, alignment_index, data_df):
    alignment = alingments[alignment_index]

    pdb_1 = data_df.loc[data_df['alignment_result_id'] == alignment, 'pdb_id_1'].values[0]
    chain_1 = data_df.loc[data_df['alignment_result_id'] == alignment, 'chain_1'].values[0]
    pdb_2 = data_df.loc[data_df['alignment_result_id'] == alignment, 'pdb_id_2'].values[0]
    chain_2 = data_df.loc[data_df['alignment_result_id'] == alignment, 'chain_2'].values[0]
    print(pdb_1)
    print(chain_1)
    print(pdb_2)
    print(chain_2)

    return pdb_1, chain_1, pdb_2, chain_2

def get_subclusters_id(alingments, alignment_index, data_df):
    alignment = alingments[alignment_index]

    subcluster_id_1 = data_df.loc[data_df['alignment_result_id'] == alignment, 'subcluster_id_1'].values[0]
    subcluster_id_2 = data_df.loc[data_df['alignment_result_id'] == alignment, 'subcluster_id_2'].values[0]
    print(subcluster_id_1)
    print(subcluster_id_2)

    return subcluster_id_1, subcluster_id_2

def get_annotation(alingments, alignment_index, data_df):
    alignment = alingments[alignment_index]

    annot = bool(data_df.loc[data_df['alignment_result_id'] == alignment, 'metamorphism'].values[0])

    return annot

def next_cluster(cluster_index, clusters_id):
    if cluster_index < len(clusters_id)-1:
        cluster_index += 1
    else:
        print("Es el último cluster")

    return cluster_index

def previous_cluster(cluster_index):
    if cluster_index > 0:
        cluster_index -= 1
    else:
        print("Es el primer cluster")

    return cluster_index

def next_alignment(alignment_index, alingments):
    if alignment_index < len(alingments)-1:
        alignment_index += 1
    else:
        print("Es el último alineamiento")

    return alignment_index

def previous_alignment(alignment_index):
    if alignment_index > 0:
        alignment_index -= 1
    else:
        print("Es el primer alineamiento")

    return alignment_index

#previous_cluster()
#next_cluster(clusters_id)
#align = get_alingments(clusters_id)
#previous_alignment()
#next_alignment(alingments=align)
#get_estructures(alingments=align)
#get_subclusters_id(alingments=align)


#sort_df = data_df.sort_values('tm_rms', ascending=False)

#print(sort_df['tm_rms'].values)



