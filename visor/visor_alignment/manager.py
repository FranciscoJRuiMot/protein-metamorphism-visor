#Manager

from visor_alignment.fatcat import fatcat_align_task
from visor_alignment.universal import US_align_task

import os
from pymol import cmd

class AlignmentManager(): #cambiar el nombre al propuesto
    def __init__(self, prot_1, prot_2, job_name):
        self.prot_1, self.prot_2, self.job_name = prot_1, prot_2, job_name

    def US_alignment(self):
        self.del_results_files() #esto seguramente se podría meter en la llamada a los alineamientos
        US_align_task(self.prot_1, self.prot_2, self.job_name)
        #lógica para cargar el fichero
        file_path = os.path.abspath(f'visor_alignment/results/{self.job_name}_all_atm_lig.pml')
        cmd.load(file_path)
        #cmd.zoom()
        cmd.bg_color("black")


    def fatcat_alignment(self):
        self.del_results_files()
        fatcat_align_task(self.prot_1, self.prot_2, self.job_name)
        #lógica para cargar el fichero
        file_path = os.path.abspath(f'visor_alignment/results/{self.job_name}.opt.twist.pdb')
        cmd.load(file_path)

        cmd.color("red", "chain A")
        cmd.color("blue", "chain B")

    def ce_alignment(self):
        self.del_results_files()
        cmd.load(self.prot_1, 'struc_1')
        cmd.load(self.prot_2, 'struc_2')

        cmd.cealign('struc_1', 'struc_2')
        #cmd.zoom()

    def call_alignment(self, alignment_type):
        align = getattr(self, alignment_type, None)
        if callable(align):
            align()
        else:
            print('Esta prueba no ha funcionado')

    def del_results_files(self):
        for file in os.listdir('visor_alignment/results'):
            os.remove(os.path.abspath(f'visor_alignment/results/{file}'))
        