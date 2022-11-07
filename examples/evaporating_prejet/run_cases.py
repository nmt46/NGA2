import numpy as np
import os
import re


# Process the input file to obtain simulation parameters
# with open('input') as fobj:
#     for line in fobj:
#         line_data = re.split(':',line)
#         if line_data[0].rstrip()=='Gravity':
#             gx=float(line_data[1].split()[0])
#             gy=float(line_data[1].split()[1])
#             gz=float(line_data[1].split()[2])
#         if line_data[0].rstrip()=='Liquid density':
#             rho_l=float(line_data[1])
#         if line_data[0].rstrip()=='Gas density':
#             rho_g=float(line_data[1])
#         if line_data[0].rstrip()=='Liquid dynamic viscosity':
#             mu_l=float(line_data[1])
#         if line_data[0].rstrip()=='Gas dynamic viscosity':
#             mu_g=float(line_data[1])
#         if line_data[0].rstrip()=='Surface tension coefficient':
#             sigma=float(line_data[1])

# Data will be stored in 'result.txt' file
# os.system('cp result.txt oldresult.txt')
# os.system('rm result.txt')

# Loop over cases, generate sigma, and run
n = 41
temps = np.linspace(300,700,n)

testSeries = 'F' # identified for series of cases

for i in range(n):
    print(f'Running simulation #{i}')
    command=f'mpiexec -n 1 nga.dp.gnu.dbg.mpi.exe -i input --"Gas inlet temperature={temps[i]}"'

    os.system(command)
    storeSim = f'cp monitor/simulation testData/simulation{testSeries}{i}'
    os.system(storeSim)
    i+=1
