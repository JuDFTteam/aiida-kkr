===============
Calculations
===============

Some intro...

Prerequisite: aiida knowledge, KKR knowledge, installation of KKRcode, Voronoi, aiida-kkr, code, computer setup
output of ``verdi calculation plugins`` should contain::
    $ verdi calculation plugins
    * kkr.kkr
    * kkr.kkrimp
    * kkr.kkrimporter
    * kkr.voro
    
In the following first set up calculation creating shapefunctions and starting potential in voronoi calculation.
Example: bulk Cu



Voronoi starting potential generator
++++++++++++++++++++++++++++++++++++

``kkr.voro``

Create an aiida structure
-------------------------

Create KKR parameter set
------------------------

Get Voronoi code
----------------

Create calculation
------------------

Submit calculation and retrieve results
---------------------------------------



KKR calculation for bulk and interfaces
+++++++++++++++++++++++++++++++++++++++

``kkr.kkr``

Update KKR parameter set
------------------------

Get code and reate calculation
------------------------------

Submit calculation and retrieve results
---------------------------------------



KKR calculation importer
++++++++++++++++++++++++

``kkr.kkrimporter``

Get remote code and set remote working directory
------------------------------------------------

Set file names
--------------

Submit KKR importer calculation and retrieve results
----------------------------------------------------



KKR impurity calculation
++++++++++++++++++++++++

``kkr.kkrimp``

Create impurity cluster
-----------------------

Writeout host Green function with KKRcode
-----------------------------------------

Create impurity starting potential
----------------------------------

Set input parameter for KKRimp
------------------------------

Submit impurity calculation and retrieve results
------------------------------------------------





