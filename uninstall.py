import os, shutil, glob

for path in ['build', 'dist', 'HAFFET.egg-info/']:
    if os.path.isdir(path):  shutil.rmtree(path)

for filename in glob.glob('*~'):
    os.remove(filename)
    
for filename in glob.glob('*/*~'):
    os.remove(filename)

for filename in glob.glob('*/*/*~'):
    os.remove(filename)
