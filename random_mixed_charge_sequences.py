from __future__ import division
import random

r2rk = 0 #proportion of R
d2de = 0 # proportion of D
length = 100
half = int(length/2)

numberof_r=int(r2rk*half)
numberof_d=int(d2de*half)
positive = ['R' for i in range(numberof_r)]+['K' for i in range(half- numberof_r)]
negative = ['D' for i in range(numberof_d)]+['E' for i in range(half- numberof_d)]

positive_idx = [i for i in range(half)]
negative_idx = [i for i in range(half)]

random.shuffle(positive_idx)
random.shuffle(negative_idx)

with open('seqs.txt','w') as o:
	o.write('length = 100, FC = 0.7, number of R = '+str(numberof_r)+', number of D = '+str(numberof_d)+'\n')
	for i in range(8):
		random.shuffle(positive_idx)
		random.shuffle(negative_idx)
		start = random.randint(0,1)
		o.write('>'+str(i+1)+'\n')
		sequence = ''
		if start%2:
			for j in range(len(positive_idx)):
				sequence+=positive[positive_idx[j]]
				sequence+=negative[negative_idx[j]]
		else:
			for j in range(len(positive_idx)):
				sequence+=negative[negative_idx[j]]
				sequence+=positive[positive_idx[j]]
		print sequence
		o.write(sequence+'\n')