import sys
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True

import fitDF.models as models
import fitDF.analyse as analyse

from methods import switch_samples

import flares
fl = flares.flares('')

model = models.DoubleSchechter()

# ====== Custom Variables ========
sim = 'flares'

if len(sys.argv) > 1:
    tag_idx = int(sys.argv[1])
else:
    tag_idx = 3

if sim == 'flares':
    tag = fl.tags[tag_idx]
else:
    tag = fl.ref_tags[tag_idx]
# =================================

print(tag)
z = float(tag[5:].replace('p','.'))

#sample_ID = '%s_gsmf_%s'%(sim,tag)
sample_ID = '%s_sfrf_%s'%(sim,tag)

a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False)
_samples = switch_samples(a.samples) 
a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False, samples=_samples)


# delete fixed parameter
a.parameters = list(a.parameters)
a.parameters.remove("alpha_2")
for ip, p in enumerate(a.parameters):
    a.median_fit[p] = np.percentile(a.samples[p], 50)


fig = a.triangle(hist2d=True, ccolor='0.5')

# plt.show()
fname = 'images/posteriors_%s_%i.png'%(sim,int(z))
print(fname)
plt.savefig(fname,dpi=150,bbox_inches='tight')
