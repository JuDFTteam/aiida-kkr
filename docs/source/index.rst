#########################################
`AiiDA`_ plugin for the `Jülich KKRcode`_
#########################################

.. image:: images/juKKR_logo_square_new.jpg
    :width: 20%
.. image:: images/AiiDA_transparent_logo.png
    :width: 40%

.. _AiiDA: http://www.aiida.net
.. _Jülich KKRcode: http://jukkr.fz-juelich.de
.. _juDFT: http://www.judft.de/pm/index.php



Welcome to documentation of the AiiDA plugin for the `Jülich KKRcode`_!
===========================================================================

The plugin is available at https://github.com/JuDFTteam/aiida-kkr

If you use this plugin for your research, please cite the following work:

.. highlights:: Philipp Rüßmann, Fabian Bertoldo, and Stefan Blügel, *The AiiDA-KKR plugin and its application to high-throughput impurity embedding into a topological insulator*, in preparation (2020).

Also please cite the original `AiiDA`_ paper:

.. highlights:: Giovanni Pizzi, Andrea Cepellotti, Riccardo Sabatini, Nicola Marzari,
  and Boris Kozinsky, *AiiDA: automated interactive infrastructure and database
  for computational science*, Comp. Mat. Sci 111, 218-230 (2016);
  http://dx.doi.org/10.1016/j.commatsci.2015.09.013; http://www.aiida.net.

Requirements
------------

- Installation of `aiida-core`_
- Installation of KKR codes (*kkrhost*, *kkrimp*, *voronoi*) of the `JuKKR package`_
- Installation of `aiida-kkr`_

Once all requirements are installed you need to `set up the computers and codes`_ before you can submit KKR calcutions using the *aiida-kkr* plugin.


.. _`aiida-core`: https://aiida-core.readthedocs.io/en/stable/installation/index.html
.. _`aiida-kkr`: https://github.com/JuDFTteam/aiida-kkr/blob/master/README.md
.. _`JuKKR package`: https://iffgit.fz-juelich.de/kkr/jukkr
.. _`set up the computers and codes`: https://aiida-core.readthedocs.io/en/stable/get_started/index.html#setup-of-computers-and-codes


User's guide
++++++++++++

.. toctree::
   :maxdepth: 4

   user_guide/index

Modules provided with aiida-kkr (API reference)
+++++++++++++++++++++++++++++++++++++++++++++++

.. toctree::
   :maxdepth: 4

   module_guide/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

