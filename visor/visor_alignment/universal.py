#US-align
import os
import subprocess

def US_align_task(prot_1, prot_2, job_name):
    resul_dir = 'visor_alignment/results/' +  job_name
    command = [os.path.abspath('binaries/USalign'), prot_1, prot_2, '-o', resul_dir]
    #print(command)

    subprocess.run(command)