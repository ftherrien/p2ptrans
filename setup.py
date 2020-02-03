import setuptools

setuptools.setup(
    setup_requires=[
        'numpy'
    ],)

from numpy.distutils.core import setup, Extension

mod = Extension(name='p2ptrans.fmodules', sources=['source/lap.f90', 'source/tiling.f90', 'source/transform.f90'], extra_f90_compile_args=["-fopenmp","-O3"], extra_link_args=['-lgomp'], f2py_options=["only:", "sphere", "parallelepiped", "circle", "center", "fastoptimization",":"])

setup(name='p2ptrans',
      version='0.0',
      description='An algorithm to match crystal structures',
      url='https://github.com/ftherrien/p2ptrans',
      author='Felix Therrien',
      author_email='felix.therrien@gmail.com',
      license='MIT',
      packages=['p2ptrans'],
      ext_modules=[mod],
      install_requires=[
          'matplotlib',
          'spglib',
        ],
      dependency_links=['git+https://github.com/pylada/pylada-light#egg=pylada'],
      scripts=['bin/p2ptrans'],)
