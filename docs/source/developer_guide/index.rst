===============
Developer guide
===============

Running the tests
+++++++++++++++++

AiiDA-KKR comes with a set of tests for its functionality. The tests are run through `pytest` and they are defined in `tests/` and the sub directories therein. 

The following will discover and run the unit test::

    # install aiida-kkr with testing extra
    pip install -e .[testing]
    # go to path where tests are defined
    cd tests
    # create fake executables
    mkdir -p jukkr; cd jukkr && export PATH="$PWD:$PATH"; touch kkr.x; touch voronoi.exe; touch kkrflex.exe; chmod +x kkr.x voronoi.exe kkrflex.exe
    # run tests (-h shows help)
    ./run_all.sh -h
    
The coverage of the tests is controlled via environment variables (see `-h` option of `run_all.sh`), e.g.::

    RUN_VORONOI=1 RUN_KKRHOST=1 ./run_all.sh
    
If you use aiida-core >= v2.0 you should first migrate the input data::

    python migrate_exports.py
    
In order to recreate test export files you need real executables instead of the fakes we create above::

    cd tests
    # this will download the code and compile voronoi, kkrhost and kkrimp
    ./jukkr_installation.sh -f
    # make sure the executables are found in the PATH
    cd jukkr && export PATH="$PWD:$PATH" && cd ..
    
If your changes require updates to reference data (checked via the `pytest-regressions` package) you should add the `--force-regen` option to the pytest run::

    pytest --force-regen workflows/test_bs_wc.py

Automatic coding style checks
+++++++++++++++++++++++++++++

Enable enable automatic checks of code sanity and coding style::

    pip install -e .[pre-commit]
    pre-commit install

After this, the `yapf <https://github.com/google/yapf>`_ formatter,
the `pylint <https://www.pylint.org/>`_ linter
and the `pylint <https://www.pylint.org/>`_ code analyzer will
run at every commit.

To run the pre-commit hooks without making a commit use::

    pre-commit run --all-files

If you ever need to skip these pre-commit hooks, just use::

    git commit -n


Continuous integration
++++++++++++++++++++++

``aiida-kkr`` comes with a ``.github`` folder that contains continuous integration tests on every commit using `GitHub Actions <https://github.com/features/actions>`_. It will:

#. run all tests for the ``django`` ORM
#. build the documentation
#. check coding style and version number (not required to pass by default)

Building the documentation
++++++++++++++++++++++++++

 #. Install the ``docs`` extra::

        pip install -e .[docs]

 #. Edit the individual documentation pages::

        docs/source/index.rst
        docs/source/developer_guide/index.rst
        docs/source/user_guide/index.rst
        ...

 #. Use `Sphinx`_ to generate the html documentation::

        cd docs
        make html

Check the result by opening ``build/html/index.html`` in your browser.

PyPI release
++++++++++++

With every tag that is pushed a continuous deployment github action runs that uploads the code to pypi.
Note that this will only be done if the tests (see Continuous integration) pass.

The latest release is therefore able to be installed via::

    pip install aiida-kkr


.. note::

   When updating the plugin package to a new version, remember to update the version number both in ``setup.json`` and ``aiida_kkr/__init__.py``.


.. _ReadTheDocs: https://readthedocs.org/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
