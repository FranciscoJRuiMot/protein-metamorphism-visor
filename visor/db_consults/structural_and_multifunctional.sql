-- Definir las CTEs
WITH max_mbl_per_protein_category AS (
    SELECT 
        p.id AS protein_id,
        gt.category,
        gr.mbl,
        gte.go_term_pair_id,
        gte.go_term_id,
        gte_other.go_term_id AS other_go_term_id,
        ROW_NUMBER() OVER (
            PARTITION BY p.id, gt.category 
            ORDER BY gr.mbl DESC
        ) AS rn
    FROM 
        public.protein p
    JOIN 
        public.go_term_pair_protein gtp 
        ON p.id = gtp.protein_id
    JOIN 
        public.go_term_pair_result gr 
        ON gtp.go_term_pair_id = gr.go_term_pair_id
    JOIN 
        public.go_term_pair_entry gte 
        ON gr.go_term_pair_id = gte.go_term_pair_id
    JOIN 
        public.go_terms gt 
        ON gte.go_term_id = gt.go_id
    JOIN 
        public.go_term_pair_entry gte_other 
        ON gr.go_term_pair_id = gte_other.go_term_pair_id 
        AND gte_other.go_term_id <> gte.go_term_id
    WHERE 
        gr.mbl IS NOT NULL  -- Asegura que mbl no sea nulo
),
pivoted_categories AS (
    SELECT 
        mm.protein_id,
        
        -- Categoría C
        MAX(mm.go_term_pair_id) FILTER (WHERE mm.category = 'C') AS C_go_term_pair_id,
        MAX(mm.mbl) FILTER (WHERE mm.category = 'C') AS C_max_mbl,
        MAX(mm.go_term_id) FILTER (WHERE mm.category = 'C') AS C_go_term_1,
        MAX(gt1.description) FILTER (WHERE mm.category = 'C') AS C_go_term_1_description,
        MAX(mm.other_go_term_id) FILTER (WHERE mm.category = 'C') AS C_go_term_2,
        MAX(gt2.description) FILTER (WHERE mm.category = 'C') AS C_go_term_2_description,
        
        -- Categoría P
        MAX(mm.go_term_pair_id) FILTER (WHERE mm.category = 'P') AS P_go_term_pair_id,
        MAX(mm.mbl) FILTER (WHERE mm.category = 'P') AS P_max_mbl,
        MAX(mm.go_term_id) FILTER (WHERE mm.category = 'P') AS P_go_term_1,
        MAX(gt1.description) FILTER (WHERE mm.category = 'P') AS P_go_term_1_description,
        MAX(mm.other_go_term_id) FILTER (WHERE mm.category = 'P') AS P_go_term_2,
        MAX(gt2.description) FILTER (WHERE mm.category = 'P') AS P_go_term_2_description,
        
        -- Categoría F
        MAX(mm.go_term_pair_id) FILTER (WHERE mm.category = 'F') AS F_go_term_pair_id,
        MAX(mm.mbl) FILTER (WHERE mm.category = 'F') AS F_max_mbl,
        MAX(mm.go_term_id) FILTER (WHERE mm.category = 'F') AS F_go_term_1,
        MAX(gt1.description) FILTER (WHERE mm.category = 'F') AS F_go_term_1_description,
        MAX(mm.other_go_term_id) FILTER (WHERE mm.category = 'F') AS F_go_term_2,
        MAX(gt2.description) FILTER (WHERE mm.category = 'F') AS F_go_term_2_description

    FROM 
        max_mbl_per_protein_category mm
    JOIN 
        public.go_terms gt1 
        ON mm.go_term_id = gt1.go_id
    JOIN 
        public.go_terms gt2 
        ON mm.other_go_term_id = gt2.go_id
    WHERE 
        mm.rn = 1  -- Mantener solo el máximo mbl por proteína y categoría
    GROUP BY 
        mm.protein_id
)

-- Consulta principal combinada
SELECT 
    /* 1. Información de las Proteínas */
    p1.id AS protein_1_id,
    p1.gene_name AS protein_1_gene_name,
    p1.description AS protein_1_description,
    
    p2.id AS protein_2_id,
    p2.gene_name AS protein_2_gene_name,
    p2.description AS protein_2_description,
    
    /* 2. Información de los Clusters */
    sub1.id AS subcluster_id_1,         -- ID del primer subcluster
    sub1.cluster_id AS cluster_id_1,    -- Cluster ID del primer subcluster
    sub2.id AS subcluster_id_2,         -- ID del segundo subcluster
    sub2.cluster_id AS cluster_id_2,    -- Cluster ID del segundo subcluster,
    
    /* 3. Información de la Secuencia */
    seq1.sequence AS sequence_1,         -- Secuencia del primer pdb
    seq2.sequence AS sequence_2,         -- Secuencia del segundo pdb
    se1.sequence_length AS sequence_length_1,
    se2.sequence_length AS sequence_length_2,
    
    /* 4. Información Estructural y Métricas */
    -- Información de las cadenas
    c1.name AS chain_1,
    s1.id AS pdb_id_1,
    st1.model_id AS model_1,
    st1.file_path AS file_path_1,
    
    c2.name AS chain_2,
    s2.id AS pdb_id_2,
    st2.model_id AS model_2,
    st2.file_path AS file_path_2,
    
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
    
    /* 5. Información Funcional y Métricas */
    -- Categorías para la proteína 1
    p1_cat.C_go_term_pair_id AS c_go_term_pair_id_p1,
    p1_cat.C_max_mbl AS c_max_mbl_p1,
    p1_cat.C_go_term_1 AS c_go_term_1_p1,
    p1_cat.C_go_term_1_description AS c_go_term_1_description_p1,
    p1_cat.C_go_term_2 AS c_go_term_2_p1,
    p1_cat.C_go_term_2_description AS c_go_term_2_description_p1,

    p1_cat.P_go_term_pair_id AS p_go_term_pair_id_p1,
    p1_cat.P_max_mbl AS p_max_mbl_p1,
    p1_cat.P_go_term_1 AS p_go_term_1_p1,
    p1_cat.P_go_term_1_description AS p_go_term_1_description_p1,
    p1_cat.P_go_term_2 AS p_go_term_2_p1,
    p1_cat.P_go_term_2_description AS p_go_term_2_description_p1,

    p1_cat.F_go_term_pair_id AS f_go_term_pair_id_p1,
    p1_cat.F_max_mbl AS f_max_mbl_p1,
    p1_cat.F_go_term_1 AS f_go_term_1_p1,
    p1_cat.F_go_term_1_description AS f_go_term_1_description_p1,
    p1_cat.F_go_term_2 AS f_go_term_2_p1,
    p1_cat.F_go_term_2_description AS f_go_term_2_description_p1,

    -- Categorías para la proteína 2
    p2_cat.C_go_term_pair_id AS c_go_term_pair_id_p2,
    p2_cat.C_max_mbl AS c_max_mbl_p2,
    p2_cat.C_go_term_1 AS c_go_term_1_p2,
    p2_cat.C_go_term_1_description AS c_go_term_1_description_p2,
    p2_cat.C_go_term_2 AS c_go_term_2_p2,
    p2_cat.C_go_term_2_description AS c_go_term_2_description_p2,

    p2_cat.P_go_term_pair_id AS p_go_term_pair_id_p2,
    p2_cat.P_max_mbl AS p_max_mbl_p2,
    p2_cat.P_go_term_1 AS p_go_term_1_p2,
    p2_cat.P_go_term_1_description AS p_go_term_1_description_p2,
    p2_cat.P_go_term_2 AS p_go_term_2_p2,
    p2_cat.P_go_term_2_description AS p_go_term_2_description_p2,

    p2_cat.F_go_term_pair_id AS f_go_term_pair_id_p2,
    p2_cat.F_max_mbl AS f_max_mbl_p2,
    p2_cat.F_go_term_1 AS f_go_term_1_p2,
    p2_cat.F_go_term_1_description AS f_go_term_1_description_p2,
    p2_cat.F_go_term_2 AS f_go_term_2_p2,
    p2_cat.F_go_term_2_description AS f_go_term_2_description_p2,
    
    /* 6. Embeddings y 3DI */
    -- Embeddings
    s31.embedding AS embedding_1,       -- Embedding para el 3di del primer pdb
    s32.embedding AS embedding_2,       -- Embedding para el 3di del segundo pdb

    /* 7. Información del visor */
    ag.is_metamorphic AS metamorphism,
    ag.comments

FROM 
    -- Consulta de alineación principal
    alignment_result ar
JOIN 
    alignment_group ag 
    ON ar.alignment_group_id = ag.id
JOIN 
    alignment_group_entry age1 
    ON ag.id = age1.alignment_group_id
JOIN 
    alignment_group_entry age2 
    ON ag.id = age2.alignment_group_id AND age1.id < age2.id
JOIN 
    subcluster_entry se1 
    ON age1.subcluster_entry_id = se1.id
JOIN 
    subcluster_entry se2 
    ON age2.subcluster_entry_id = se2.id
JOIN 
    subcluster sub1 
    ON se1.subcluster_id = sub1.id
JOIN 
    subcluster sub2 
    ON se2.subcluster_id = sub2.id
JOIN 
    structure_3di s31 
    ON se1.structure_3di_id = s31.id
JOIN 
    structure_3di s32 
    ON se2.structure_3di_id = s32.id
JOIN 
    state st1 
    ON s31.state_id = st1.id
JOIN 
    state st2 
    ON s32.state_id = st2.id
JOIN 
    structure s1 
    ON st1.structure_id = s1.id
JOIN 
    structure s2 
    ON st2.structure_id = s2.id
JOIN 
    chain c1 
    ON st1.chain_id = c1.id
JOIN 
    chain c2 
    ON st2.chain_id = c2.id
JOIN 
    sequence seq1 
    ON c1.sequence_id = seq1.id      -- Secuencia de la primera cadena
JOIN 
    sequence seq2 
    ON c2.sequence_id = seq2.id      -- Secuencia de la segunda cadena

-- Uniones adicionales para obtener información de las proteínas
JOIN 
    protein p1 
    ON s1.protein_id = p1.id            -- Información de la primera proteína
JOIN 
    protein p2 
    ON s2.protein_id = p2.id            -- Información de la segunda proteína

-- Unir las categorías para las proteínas 1 y 2
LEFT JOIN 
    pivoted_categories p1_cat 
    ON p1.id = p1_cat.protein_id
LEFT JOIN 
    pivoted_categories p2_cat 
    ON p2.id = p2_cat.protein_id

ORDER BY 
    /* 1. Ordenar por Proteína */
    p1.id ASC,
    p2.id ASC,
    
    /* 2. Ordenar por Clusters */
    sub1.cluster_id ASC,
    sub1.id ASC,
    sub2.id ASC,
    
    /* 3. Ordenar por Métricas */
    ar.fc_rms DESC;