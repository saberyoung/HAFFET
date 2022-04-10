from setuptools import setup
from sdapy. __version__ import version, description

setup(
    name='sn_data_analysis',
    version=version,
    description=description,
    python_requires='>=2.7',
    packages=[
        'sdapy',
    ],   
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
    zip_safe = False
)
