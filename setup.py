from skbuild import setup

setup(name='p2ptrans',
      version='2.0.2',
      description='An algorithm to match crystal structures',
      url='https://github.com/ftherrien/p2ptrans',
      author='Felix Therrien',
      author_email='felix.therrien@gmail.com',
      license='MIT',
      packages=['p2ptrans'],
      python_requires='>=3',
      install_requires=[
          'numpy',
          'matplotlib',
          'spglib',
          'pylada @ git+https://github.com/pylada/pylada-light#egg=pylada',
        ],
      scripts=['bin/p2ptrans', 'bin/p2pint', 'bin/mview'],)
