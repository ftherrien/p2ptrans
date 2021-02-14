import setuptools
from numpy.distutils.core import setup, Extension

mod = Extension(name='p2ptrans.fmodules', sources=['source/lap.f90', 'source/utils.f90', 'source/tiling.f90', 'source/potential.f90', 'source/transform.f90'], extra_f90_compile_args=["-fopenmp","-O3"], extra_link_args=['-lgomp'], f2py_options=["only:","munkres", "free_trans", "rot_mat", "center", "eye", "norm", "split", "det", "sort", "sphere", "circle", "distance", "derivative", "closest", "fixed_tmat", "fixed_tmat_int", "intoptimization", "fastoptimization", ":"])

setup(name='p2ptrans',
      version='2.0 beta 6',
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
      scripts=['bin/p2ptrans', 'bin/p2pint'],)
