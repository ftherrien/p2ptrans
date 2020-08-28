A Simple Guide to p2ptrans
==========================

Welcome to the documentation for p2ptrans and thank you for your interest! I will try my best to add as much information as possible to this guide; if you encounter problems or if you have any suggestion to improve the package or the documentation, please feel free to raise an issue or--even better--to create a pull request. 

.. _install:

Installation
^^^^^^^^^^^^

To install p2ptrans simply run:

.. code-block:: console

   pip install git+https://github.com/ftherrien/p2ptrans

If you do not have `pylada
<https://github.com/pylada/pylada-light>`_, you will need to install the `py` module first:

.. code-block:: console

   pip install py

On certain systems, the pylada installation fails with ``error: ‘v’ does not name a type``. If you encounter this error retry the installation with:

.. code-block:: console

   CXXFLAGS="-std=c++11" pip install git+https://github.com/ftherrien/p2ptrans

.. toctree::
   :maxdepth: 1
   :caption: Documentation

   p2ptrans
   p2pint

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   martensite
   interfaces
   pymodule
