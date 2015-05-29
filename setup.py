__author__ = 'svfilhol'

from distutils.core import setup

setup(
    name='dempy',
    version='0.1dev',
    packages=['dempy',],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    description='Useful tool for DEM and point cloud processing',
    long_description=open('README.txt').read(),
    author='Simon Filhol',
    author_email='svfilhol@alaska.edu',
    scripts=['demtool/fftdem.py','demtool/linearDetrend.py','demtool/MatrixGenerator.py','demtool/morphometry.py',
             'demtool/plotting.py','demtool/rasterTool.py','demtool/smooth2d.py',
             'pctool/pointcloud.py',
             'misceallenous/Qgis_tool.py'],
    url='https://github.com/ArcticSnow/dempy',
    install_requires=[
        "pandas >= 0.16.1",
        "numpy >= 1.9.2",
        "gdal >= 1.9.0",
        "matplotlib >= 1.4.2",
        "pyfftw >= 0.9.2"
    ],
)