===============
Developer guide
===============

Running the tests
+++++++++++++++++

The following will discover and run the unit test::

    pip install -e .[testing]
    cd aiida_kkr/tests
    ./run_all.sh

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
