# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:38:16 2024

@author: sarah
"""
import csv
import numpy as np

from pathlib import Path


#------------------------------------------------------------------------------
def read_stokes(filename, nm):
    u = np.zeros(nm)
    v = np.zeros(nm)
    psi = np.zeros(nm)

    with open(filename, newline='') as file:
        reader = csv.reader(file)

        for i in range(nm):
            line = next(reader)
            ui, vi, psii = line[0].split(' ')
            u[i] = float(ui)
            v[i] = float(vi)
            psi[i] = float(psii)
        cnvg_err = int(next(reader)[0])
        file.close()
    return u, v, psi, cnvg_err

def write_stokes(dmn, u, v, psi, cnvg_err):
    nm = dmn.Nx * dmn.Ny
    filename = dmn.filestr + ".csv"
    Path(dmn.dirstr).mkdir(parents=True, exist_ok=True)
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(nm):
            writer.writerow([u[i], v[i], psi[i]])
        
        writer.writerow([cnvg_err])
        file.close()    


    
    
    
    
    