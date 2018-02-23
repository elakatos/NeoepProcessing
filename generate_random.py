from numpy import random

aaOneLetter = 'GPAVLIMCFYWHKRQNEDST'

with open('Samples/random_epitopes.fasta', 'w') as outFile:
	for i in range(500):
		fastaLine = ''.join(random.choice(list(aaOneLetter), 19, replace=True))
		outFile.write('>line'+str(i)+'\n'+fastaLine+'\n')

