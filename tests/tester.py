#######################################################################################################
# This python script does automatic test runs of MCML.
# For some sets of input parameters it runs all combinations.
# Parameters with physical meaning are taken from a Medical & Biological Engineering & Computing paper.
#######################################################################################################

import sys
import os
import subprocess
import tempfile








# Modelling the sampling volume for skin blood oxygenation measurements (Meglinsky, 2001)
# Table 2 Optical properties of skin model (lambda=600nm, Caucasian skin)
# k Skin layer             mus mua   g    n
# 1 stratum corneum        100 0.02  0.9  1.53
# 2 living epidermis       40  0.015 0.85 1.34
# 3 papillary dermis       30  0.07  0.8  1.4
# 4 upper blood net dermis 35  0.1   0.9  1.39
# 5 dermis                 20  0.07  0.76 1.4
# 6 deep blood net dermis  35  0.1   0.95 1.39
# 7 subcutaneous fat       15  0.03  0.8  1.44
# The refractive index of air is n=1


# Input values to test:

mua = [0.1, 1, 10] # absorption events per unit
mus = [1, 10, 100] # scattering events per unit
g = [0, 0.5, 0.9]  # scattering anisotropy https://omlc.org/classroom/ece532/class3/hg.html
n = [1.34, 1.4, 1.44] # refractive indices of Caucasian skin (from table 2 in )

PHOTON_MILLIONS = 10

boundary_presets = { # user can either choose one from here or use implicit flat boundaries
	'flat50': {
		'description': 'A trivial flat boundary with 50 samples and uniform spacing 0.1.',
		'data': [0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1]
	} #TODO smaller spacing
}

dz = 0.01 #TODO more resolution
dr = 0.01
nz = 200
nr = 500
na = 10

dl = 0.5 # layer thickness







useMPI = True
mpiDebugLevel = 0
mpiMachinefile = 'mpi.txt'
mpiCmd = ['mpiexec', '-debug', str(mpiDebugLevel), '-lines', '-machinefile', mpiMachinefile, '-pwd', 'S>v2zv]08ct/lg4+9{lvi-[bcv-Qc4[|]$ne1NE~']







# Interpret options:

if len(sys.argv) < 2:
	print('Usage: python tester.py <mcml executable> <optional executable args> -b <optional list of explicit boundaries(*)> -o <optional output folder>')
	
	print('\nThis script will run multiple MCML simulations with a user-specified implementation, producing mco files.')
	print('Plots will be created at the end from the mco files, for which it depends on another script plotter.py.')
	
	print('\nYou can specify an executable in any folder, which the exe will recognize as its current working directory, e.g. to find kernel files to compile.')
	
	print('\nHardcoded at the top of the script you can find lists of values for MCML input parameters.')
	print('The script will iterate these lists to start a simulation for all possible combinations.')
	
	print('\nOutput files have the following naming convention:')
	print('<number of photons>_<absorption coefficient mua>_<scattering coefficient mus>_<scattering anisotropy g>_<refractive index n>.mco')
	print('For each mco file, a folder with the same name is created to store the plots.')
	
	print('\n(*) Explicit boundary presets:')
	for name in boundary_presets: print(name + ' - ' + boundary_presets[name]['description'])
	
	print('\nIf no explicit boundaries are specified, implicit flat boundaries are used.')
	print('If one explicit boundary is specified, it is used above and below one layer.')
	print('If two explicit boundaries are specified, the first is used above and the second below one layer.')
	print('If more than two explicit boundaries are specified, the same layer is used multiple times so that all boundaries are used.')
	exit()


mcml_exe = []
while len(sys.argv)>1 and sys.argv[1]!='-b' and sys.argv[1]!='-o':
	mcml_exe.append(sys.argv.pop(1))

explicit_boundaries = []
if len(sys.argv)>2 and sys.argv[1]=='-b':
	sys.argv.pop(1)
	while len(sys.argv)>1 and sys.argv[1]!='-o':
		b = sys.argv.pop(1)
		try:
			explicit_boundaries.append(boundary_presets[b]['data'])
		except KeyError:
			print('Preset boundary ' + b + ' does not exist')
			exit()

out_folder = os.getcwd() + '/'
if len(sys.argv)>2 and sys.argv.pop(1)=='-o':
	out_folder += sys.argv.pop(1) + '/'
if not os.path.exists(out_folder): os.makedirs(out_folder)







# Run simulations:

mcos = []

for mua_i in mua:
	for mus_i in mus:
		for g_i in g:
			for n_i in n:

				mco = os.path.abspath(out_folder + str(PHOTON_MILLIONS)+'M_'+str(mua_i)+'mua_'+str(mus_i)+'mus_'+str(g_i)+'g_'+str(n_i)+'n.mco')

				if (os.path.isfile(mco)):
					print('Skip existing ' + str(mco))
					mcos.append(mco)
					continue

				mci = os.path.join(tempfile.mkdtemp(), 'test.mci')
				with open(mci, 'w') as f:
					f.write('1.0\n') # file version
					f.write('1\n') # number of sims

					# sim 1
					f.write(mco+' A\n') # output filename, ASCII
					f.write(str(PHOTON_MILLIONS*1000000) + '\n') # number of photons
					f.write(str(dz)+' '+str(dr)+'\n') # detection deltas in cm (dz, dr)
					f.write(str(nz)+' '+str(nr)+' '+str(na)+'\n') # detection array dims (z, r, a)

					# number of layers
					if len(explicit_boundaries) > 1:
						f.write(str(len(explicit_boundaries)-1) + '\n')
					else:
						f.write('1\n')

					f.write('1.0\n') # n above

					if len(explicit_boundaries) > 1:
						for i in range(0, len(explicit_boundaries)-1):
							f.write(' '.join(str(f) for f in explicit_boundaries[i]) + '\n')
							f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(dl)+'\n')
						f.write(' '.join(str(f) for f in explicit_boundaries[len(explicit_boundaries)-1]) + '\n')
					elif len(explicit_boundaries) > 0:
						f.write(' '.join(str(f) for f in explicit_boundaries[0]) + '\n')
						f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(dl)+'\n')
						f.write(' '.join(str(f) for f in explicit_boundaries[0]) + '\n')
					else:
						f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(dl)+'\n') # layer 1

					f.write('1.0') # n below


				print('===========================================================')
				print('[STARTING] Simulation ' + str(len(mcos)+1) + ' of ' + str(len(mua)*len(mus)*len(g)*len(n)))
				print('===========================================================')
				print('[INPUT] ' + mci + ':')
				with open(mci, 'r') as f:
					print(f.read())
				print('===========================================================')

				success = False

				while not success:

					pathPair = os.path.split(mcml_exe[0]) #(head, tail) where tail is the last pathname component and head is everything leading up to that

					cmd = [*mcml_exe, mci]

					if useMPI:
						cmd[0] = pathPair[1]
						cmd = mpiCmd + cmd

					print('[EXEC] ' + str(cmd) + ':')
					result = subprocess.run(cmd, cwd=pathPair[0])

					print('===========================================================')
					try:
						result.check_returncode()
						print('[OUTPUT] ' + mco)

						print('[PLOTTING]')
						subprocess.run(['python', 'plotter.py', mco], stdin=subprocess.DEVNULL, universal_newlines=True)

						mcos.append(mco)

						success = True

					except subprocess.CalledProcessError as err:
						print('[FAILED] ' + str(err))
						print('Retrying...')
					print('\n\n')


print('[FINISHED] ' + str(len(mcos)) + ' mco files produced')