import setuptools
from numpy.distutils.core import setup, Extension

sources = ['source/lapjv.hpp', 'source/lapjv.cpp', 'source/lapjv.f90', 'source/lap.f90', 'source/utils.f90', 'source/tiling.f90', 'source/potential.f90', 'source/transform.f90']
functions = ["only:","munkres", "free_trans", "rot_mat", "center", "eye", "norm", "split", "det", "sort", "sphere", "circle", "distance", "derivative", "intoptimization", "fastoptimization", "lapjv", ":"]

mod = Extension(name='p2ptrans.fmodules', sources=sources, extra_f90_compile_args=["-fopenmp","-O3"], extra_link_args=['-lgomp', '-lgfortran'], f2py_options=functions)

setup(name='p2ptrans',
      version='2.0.2',
      description='An algorithm to match crystal structures',
      url='https://github.com/ftherrien/p2ptrans',
      author='Felix Therrien',
      author_email='felix.therrien@gmail.com',
      license='MIT',
      packages=['p2ptrans'],
      ext_modules=[mod],
      python_requires='>=3',
      install_requires=[
          'numpy',
          'matplotlib',
          'spglib',
          'pylada @ git+https://github.com/pylada/pylada-light#egg=pylada',
        ],
      scripts=['bin/p2ptrans', 'bin/p2pint', 'bin/mview'],)
