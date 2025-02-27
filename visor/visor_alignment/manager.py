#Manager

from visor_alignment.fatcat import fatcat_align_task
from visor_alignment.universal import US_align_task
from logic import apply_pdb_colors_and_alignment_path

import os
from pymol import cmd

class AlignmentManager(): #cambiar el nombre al propuesto
    def __init__(self, prot_1, name_1, prot_2, name_2 ,job_name):
        self.prot_1, self.name_1, self.prot_2, self.name_2, self.job_name = prot_1, name_1, prot_2, name_2 ,job_name

    def US_alignment(self):
        self.del_results_files() #esto seguramente se podría meter en la llamada a los alineamientos
        US_align_task(self.prot_1, self.prot_2, self.job_name)
        #lógica para cargar el fichero
        file_path = os.path.abspath(f'visor_alignment/results/{self.job_name}_all_atm_lig.pml')
        cmd.load(file_path)
        #Añadir como cambiar el nombre de las estructuras leyendo la salida del principio
        ##Parece que la estructura 1 es la 1 y la 2 es la 2
        cmd.bg_color("black")
        cmd.set_name("structure1", f"{self.name_1}")
        cmd.set_name("structure2", f"{self.name_2}")
        cmd.zoom()

    def fatcat_alignment(self):
        self.del_results_files()
        fatcat_align_task(self.prot_1, self.name_1, self.prot_2, self.name_2, self.job_name)
        #lógica para cargar el fichero
        file_path = os.path.abspath(f'visor_alignment/results/{self.job_name}.opt.twist.pdb')
        cmd.load(file_path, f"fatcat_align_{self.job_name}") #

        ##Cambios para manipular mejor las cadenas
        ##La cadena A es la prot_1 y la B es la 2
        cmd.hide('everything')
        cmd.extract(f"{self.name_1}", f"fatcat_align_{self.job_name} and chain A")
        cmd.show('cartoon', f"{self.name_1}")
        cmd.color("red" ,f"{self.name_1}")

        cmd.extract(f"{self.name_2}", f"fatcat_align_{self.job_name} and chain B")
        cmd.show('cartoon', f"{self.name_2}")
        cmd.color("blue" ,f"{self.name_2}")
        cmd.zoom()
        

    def ce_alignment(self):
        self.del_results_files()
        cmd.load(self.prot_1, self.name_1)
        cmd.load(self.prot_2, self.name_2)

        cmd.cealign(self.name_1, self.name_2)
        cmd.zoom()

        apply_pdb_colors_and_alignment_path(self.name_1, self.name_2)

    def call_alignment(self, alignment_type):
        align = getattr(self, alignment_type, None)
        if callable(align):
            align()
        else:
            print('Esta prueba no ha funcionado')

    def del_results_files(self):
        for file in os.listdir('visor_alignment/results'):
            os.remove(os.path.abspath(f'visor_alignment/results/{file}'))
        