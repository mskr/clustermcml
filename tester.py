#######################################################################################################
# This python script does automatic test runs of MCML.
# For some sets of input parameters it runs all combinations.
# Parameters with physical meaning are taken from a Medical & Biological Engineering & Computing paper.
#######################################################################################################

import sys
import os
import subprocess
import tempfile
from os import listdir

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

mua = [0.1, 1, 10] # absorption events per unit
mus = [1, 10, 100] # scattering events per unit
g = [0, 0.5, 0.9]  # scattering anisotropy https://omlc.org/classroom/ece532/class3/hg.html
n = [1.34, 1.4, 1.44] # refractive indices of Caucasian skin (from table 2 in )

PHOTON_MILLIONS = 10

boundaries = {
	'flat50': [0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1]
}




if len(sys.argv) < 2:
	print('Usage: python tester.py <mcml executable> -eb <optional list of explicit boundaries*>')
	print('\n* Preset explicit boundaries are:')
	for b in boundaries: print(b)
	print('\nIf no explicit boundaries are specified, implicit flat boundaries are used.')
	print('If one explicit boundary is specified, it is used above and below one layer.')
	print('If two explicit boundaries are specified, the first is used above and the second below one layer.')
	print('If more than two explicit boundaries are specified, the same layer is used multiple times so that all boundaries are used.')
	exit()


mcml_exe = []
while len(sys.argv)>1 and sys.argv[1]!='-eb':
	mcml_exe.append(sys.argv.pop(1))

explicit_boundaries = []
if len(sys.argv)>2 and sys.argv[1]=='-eb':
	for i in range(2, len(sys.argv)):
		b = sys.argv[i]
		try:
			explicit_boundaries.append(boundaries[b])
		except KeyError:
			print('Preset boundary ' + b + ' does not exist')
			exit()


mcos = []

for mua_i in mua:
	for mus_i in mus:
		for g_i in g:
			for n_i in n:

				mco = str(PHOTON_MILLIONS)+'M_'+str(mua_i)+'mua_'+str(mus_i)+'mus_'+str(g_i)+'g_'+str(n_i)+'n.mco'

				mci = os.path.join(tempfile.mkdtemp(), 'test.mci')
				with open(mci, 'w') as f:
					f.write('1.0\n') # file version
					f.write('1\n') # number of sims

					# sim 1
					f.write(mco+' A\n') # output filename, ASCII
					f.write(str(PHOTON_MILLIONS*1000000) + '\n') # number of photons
					f.write('0.01 0.01\n') # detection deltas in cm (dz, dr)
					f.write('200 500 10\n') # detection array dims (z, r, a)

					# number of layers
					if len(explicit_boundaries) > 1:
						f.write(str(len(explicit_boundaries)-1) + '\n')
					else:
						f.write('1\n')

					f.write('1.0\n') # n above

					if len(explicit_boundaries) > 1:
						for i in range(0, len(explicit_boundaries)-1):
							f.write(' '.join(str(f) for f in explicit_boundaries[i]) + '\n')
							f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(10)+'\n')
						f.write(' '.join(str(f) for f in explicit_boundaries[len(explicit_boundaries)-1]) + '\n')
					elif len(explicit_boundaries) > 0:
						f.write(' '.join(str(f) for f in explicit_boundaries[0]) + '\n')
						f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(10)+'\n')
						f.write(' '.join(str(f) for f in explicit_boundaries[0]) + '\n')
					else:
						f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(10)+'\n') # layer 1

					f.write('1.0') # n below


				print('===========================================================')
				print('[INPUT] ' + mci + ':')
				with open(mci, 'r') as f:
					print(f.read())
				print('===========================================================')
				print('[EXEC] ' + str([*mcml_exe, mci]) + ':')
				result = subprocess.run([*mcml_exe, mci],
										stdin=subprocess.DEVNULL,
										universal_newlines=True)  # Will be "text" in Python 3.7
				print('===========================================================')
				try:
					result.check_returncode()
					print('[OUTPUT] ' + mco)
				except CalledProcessError as err:
					print('[FAILED] ' + str(err))
				print('\n\n')

				mcos.append(mco)


print('[FINISHED] Now plotting...')
result = subprocess.run(['python', 'plotter.py', *mcos], stdin=subprocess.DEVNULL, universal_newlines=True)