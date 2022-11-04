from setuptools import find_packages, setup
from sdapy. __version__ import name, version, description, python_requires, author, author_email, url

setup(
    name=name,
    version=version,
    description=description,
    python_requires=python_requires,
    url=url,
    author=author,
    author_email=author_email,    
    packages=find_packages(),
    package_data={
        '': [
            'data/*.txt',
            'data/**/*',
        ],
    },        
    entry_points = {
        'console_scripts' :
        ['sdapy_run = sdapy.snerun:main',
         'sdapy_gui = sdapy.sneGui:main',],
    },
    install_requires=[
        'astropy',
        'numpy',
        'pandas',
        'scipy',
        'joblib',
        'matplotlib',
    ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    zip_safe = False,
    include_package_data=True
)
