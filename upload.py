from __future__ import print_function
from shutil import rmtree
import os,sys
from src.__version__ import version

if __name__ == '__main__':
    here = os.path.abspath(os.path.dirname(__file__))
    try:  rmtree(os.path.join(here, 'dist'))
    except OSError:  pass

    #print ('Building Source and Wheel (universal) distribution...')
    #os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

    #print ('Uploading the package to PyPI via Twine...')
    #os.system('twine upload --skip-existing --verbose dist/*')

    print ('Pushing git tags...')
    os.system('git add .')
    os.system('git commit -m {0}'.format(version))
    os.system('git push')
