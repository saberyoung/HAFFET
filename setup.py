from setuptools import setup
from src. __version__ import version

setup(
    name='sn_data_analysis',
    version=version,
    python_requires='>=2.7',
    packages=[
        'src',        
    ],   
    zip_safe = False
)
