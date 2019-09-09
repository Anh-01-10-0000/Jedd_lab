from __future__ import division
from Bio import SeqIO
import csv,itertools,re, numpy
from localcider.sequenceParameters import SequenceParameters

kappa = 1
charge_bal = 0.7
window_sizes = [60]
cutoffs = [0.6]

def score_polyamp(s, window, cutoff):
    #s: amino acid seq, window: window size, cutoff: the minimum score to be considered as significant
    #returns a string combining all the mixed charge regions
    polyamph_regions = []
    polyamph_start = 0
    polyamph_end = 0
    last_start = 0
    last_end = 0
    for j in range(window//2,1+len(s)-window//2):
        region = ''
        start = j-(window//2)
        end = j+(window//2)
        region = s[start:end]

        pos = region.count('R')+region.count('K')
        neg = region.count('E')+region.count('D')
        score = pos+neg

        if score>=cutoff:
            if polyamph_start == 0:
                extended_reg = s[polyamph_start:end]
                if extended_reg.count('R')+extended_reg.count('K')+extended_reg.count('E')+extended_reg.count('D')>=cutoff:
                    polyamph_start = start
                    polyamph_end = end
                    last_start = polyamph_start
                    last_end = polyamph_end
            elif polyamph_end+1 > start: #if next start is within 10 from last end, join the 2
                if extended_reg.count('R')+extended_reg.count('K')+extended_reg.count('E')+extended_reg.count('D')>=cutoff:
                    polyamph_end = end
                    last_start = polyamph_start
                    last_end = polyamph_end
            else:
                if extended_reg.count('R')+extended_reg.count('K')+extended_reg.count('E')+extended_reg.count('D')>=cutoff:
                    polyamph_regions.append(s[polyamph_start:polyamph_end+1])
                    polyamph_start = 0
                    polyamph_end = 0
    if last_end!=0: 
        temp1 = s[last_start:last_end+1]
        if temp1 not in polyamph_regions:
            polyamph_regions.append(temp1)
    return polyamph_regions

def get_polyamp_regions(inputFasta,outputFile):
    # look for mixed charge regions in each protein. Each entry is preceeded by its accession
    seqs,out,added = [],[],{}
    with open(inputFasta,'rU') as infile, open(outputFile,'wb') as outf:
        reader = SeqIO.parse(infile,'fasta')
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(['accession','name','highly charged region','r/(r+k)','d/(d+e)','length','charge_density','kappa'])
        for s in reader:
            tmp = s.id.split('|')
            acc = tmp[1]
            tmp2 = s.description.split()[-3]
            name = tmp2[tmp2.find('=')+1:]
            t2=score_polyamp(str(s.seq),window_size,cut_off*window_size)
            if len(t2)>0:
                for i in range(len(t2)):
                    to_write = [acc,name, '', 0,0, 0,0,1]
                    r,k = t2[i].count('R'), t2[i].count('K')
                    d,e = t2[i].count('D'), t2[i].count('E')
                    neg = d+e #+t2[i].count('X')
                    pos = r+k
                    if pos>0 and neg>0 and pos/neg>=charge_bal and neg/pos>=charge_bal:
                        l = len(t2[i])
                        r_ratio = r/pos
                        d_ratio = d/neg
                        region_noX = ''
                        for char in t2[i]:
                            if char=='X':
                                region_noX = region_noX+'S'
                            elif char=='S' or char == 'U':
                                region_noX = region_noX+'T'
                            else:
                                region_noX = region_noX+char
                        SeqOb = SequenceParameters(region_noX)
                        kap = SeqOb.get_kappa_X(grp1 = ['E','D'], grp2 = ['K','R'])
                        if kap<=kappa:
                            to_write = [acc,name, t2[i], str(r_ratio),str(d_ratio), l,(pos+neg)/l ,str(kap)]
                        writer.writerow(to_write)


path = "C:\\Users\\NTAnh\\Dropbox\\TLL\\polyamph\\no_domain_removed\\"
path2 = "C:\\Users\\NTAnh\\Dropbox\\TLL\\polyamph\\"
human_proteome = path2 + 'uniprot_human_phosphorylated.fasta'#'uniprot_human_phosphorylated_no_dna_rna.fasta'



for window_size in window_sizes:
    for cut_off in cutoffs:
        outfile = path + str(window_size)+'_'+str(cut_off)+'_k'+str(kappa)+'bal'+str(charge_bal)+"polyamph_no_phos_separated.txt"
        get_polyamp_regions(human_proteome,outfile)

### get a list of those annotated to be in speckles
speck = {}
with open(path2+'panther_speckle.txt','r') as inf:
    reader = csv.reader(inf,delimiter = '\t')
    for line in reader:
        speck[line[0]] = line[1]




import seaborn as sns
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
sns.set()
plt.style.use('ggplot')

###plot speckle
count = 0
for window_size in window_sizes:
    for cut_off in cutoffs:
        df = pd.read_csv(path + str(window_size)+'_'+str(cut_off)+'_k'+str(kappa)+'bal'+str(charge_bal)+"polyamph_no_phos_separated.txt",sep='\t')
        cmap = sns.cubehelix_palette(start=0,as_cmap=True, reverse=True, dark=0.05, light = 0.7)
        f, ax = plt.subplots()
        plt.axis([-0.1,1.1,-0.1,1.1])
        # ax.set_autoscale_on(False)

        x = []
        y = []
        colors = []
        sizes = []

        for j in range(len(df)):
            if df['accession'][j] in speck and df['charge_density'][j]>0:
                count+=1
                x.append(df['r/(r+k)'][j])
                region = df['highly charged region'][j]
                d = region.count('D')
                e = region.count('E')
                if d+e==0:
                    y.append(0)
                else:
                    y.append(d/(d+e))
                colors.append(df['kappa'][j])
                sizes.append(df['length'][j]/5)

        points = ax.scatter(x,y, c=colors, s=sizes, alpha = 0.5, cmap=cmap)
        plt.xlabel('R/(R+K)')
        plt.ylabel('D/(D+E)')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))

        plt.title('Window size: '+str(window_size)+', FC: '+str(cut_off))
        f.colorbar(points,boundaries=numpy.linspace(0,1,11),label = 'kappa') 
        # plt.show()
        plt.savefig(path+str(window_size)+'_'+str(cut_off)+'bal_'+str(charge_bal)+'k'+str(kappa)+'_speck_no_phos.pdf')
        plt.clf()
        plt.close()
