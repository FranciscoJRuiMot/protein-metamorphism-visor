SELECT 
    -- Información de las cadenas
    c1.name AS chain_1,
    s1.id AS pdb_id_1,
    se1.sequence_length AS sequence_length_1,
    st1.model_id AS model_1,
    seq1.sequence AS sequence_1,       -- Secuencia del primer pdb
    
    c2.name AS chain_2,
    s2.id AS pdb_id_2,
    se2.sequence_length AS sequence_length_2,
    st2.model_id AS model_2,
    seq2.sequence AS sequence_2,       -- Secuencia del segundo pdb

    -- Información de las proteínas
    p1.id AS protein_1,
    p1.gene_name AS gene_name_1,
    p1.description AS description_1,
    
    p2.id AS protein_2,
    p2.gene_name AS gene_name_2,
    p2.description AS description_2,
    
    -- Información del resultado de alineación
    ar.id AS alignment_result_id,
    ar.ce_rms,
    ar.tm_rms,
    ar.tm_seq_id,
    ar.tm_score_chain_1,
    ar.tm_score_chain_2,
    ar.fc_rms,
    ar.fc_identity,
    ar.fc_similarity,
    ar.fc_score,
    ar.fc_align_len,

    -- Información de los clusters
    sub1.id AS subcluster_id_1,         -- ID del primer subcluster
    sub1.cluster_id AS cluster_id_1,    -- Cluster ID del primer subcluster
    sub2.id AS subcluster_id_2,         -- ID del segundo subcluster
    sub2.cluster_id AS cluster_id_2,    -- Cluster ID del segundo subcluster

    -- Embeddings
    s31.embedding AS embedding_1,       -- Embedding para el 3di del primer pdb
    s32.embedding AS embedding_2,       -- Embedding para el 3di del segundo pdb

    -- Rutas de archivo
    st1.file_path AS file_path_1,
    st2.file_path AS file_path_2,

    -- Nuevas columnas de alignment_group
    ag.is_metamorphic AS metamorphism,
    ag.comments

FROM 
    alignment_result ar
JOIN 
    alignment_group ag ON ar.alignment_group_id = ag.id
JOIN 
    alignment_group_entry age1 ON ag.id = age1.alignment_group_id
JOIN 
    alignment_group_entry age2 ON ag.id = age2.alignment_group_id AND age1.id < age2.id
JOIN 
    subcluster_entry se1 ON age1.subcluster_entry_id = se1.id
JOIN 
    subcluster_entry se2 ON age2.subcluster_entry_id = se2.id
JOIN 
    subcluster sub1 ON se1.subcluster_id = sub1.id
JOIN 
    subcluster sub2 ON se2.subcluster_id = sub2.id
JOIN 
    structure_3di s31 ON se1.structure_3di_id = s31.id
JOIN 
    structure_3di s32 ON se2.structure_3di_id = s32.id
JOIN 
    state st1 ON s31.state_id = st1.id
JOIN 
    state st2 ON s32.state_id = st2.id
JOIN 
    structure s1 ON st1.structure_id = s1.id
JOIN 
    structure s2 ON st2.structure_id = s2.id
JOIN 
    chain c1 ON st1.chain_id = c1.id
JOIN 
    chain c2 ON st2.chain_id = c2.id
JOIN 
    sequence seq1 ON c1.sequence_id = seq1.id      -- Secuencia de la primera cadena
JOIN 
    sequence seq2 ON c2.sequence_id = seq2.id      -- Secuencia de la segunda cadena

-- Uniones adicionales para obtener información de las proteínas
JOIN 
    protein p1 ON s1.protein_id = p1.id            -- Información de la primera proteína
JOIN 
    protein p2 ON s2.protein_id = p2.id            -- Información de la segunda proteína

ORDER BY 
    cluster_id_1 ASC,
    subcluster_id_1 ASC,
    subcluster_id_2 ASC,
    fc_rms DESC;
