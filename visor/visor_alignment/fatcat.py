#FATCAT

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBIO
import os
import subprocess

def cif_to_pdb(cif_path, pdb_path): #mudar al archivo de fatcat
    """Convierte un archivo CIF a PDB."""

    parser = MMCIFParser()
    structure = parser.get_structure('ID', cif_path)
    io = PDBIO()
    io.set_structure(structure)

    for model in structure:
        for chain in model:
            if len(chain.id) > 1:
                chain.id = "X"
    io.save(pdb_path)

def fatcat_align_task(prot_1, name_1, prot_2, name_2, job_name):
    resul_dir = 'visor_alignment/results/' +  job_name
    cif_to_pdb(prot_1, f'visor_alignment/results/{name_1}.pdb')
    cif_to_pdb(prot_2, f'visor_alignment/results/{name_2}.pdb')

    command = [os.path.abspath('binaries/FATCAT'),'-p1', f'visor_alignment/results/{name_1}.pdb', '-p2', f'visor_alignment/results/{name_2}.pdb', '-o', resul_dir ,'-m', '-ac', '-t']
    print(command)
    subprocess.run(command)