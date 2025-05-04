import numpy as np
indir = "/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/CPC_geometry_1024/"
file = f"{indir}/phi_1024_19660_0.0001__CPC_40_cohesin_16_eps_0.007504684956431058.txt"
def LastNlines(file_name, N):
    # opening file using with() method
    # so that file get closed
    # after completing work
    lines = np.empty((0,1024))
    
    with open(file_name) as file:
        # loop to read iterate 
        # last n lines and print it
        for line in (file.readlines() [-N:]):
            line = np.array([float(i) for i in line[0:-2].split(" ")]).reshape(1,-1)
            lines = np.append(lines, line, axis = 0)
    return lines

lines = (LastNlines(file, 1024))

np.savetxt(f"{indir}/tmp.txt",lines)