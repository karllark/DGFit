#################
Dust Grain Fitter
#################

``DGFit`` is a python package to derive dust grain size and composition
distributions based on fitting observations of interstellar dust.

User Documentation
==================

.. toctree::
   :maxdepth: 2

   Running DGFit <dgfit/running_dgfit.rst>
   Size Distributions <dgfit/size_distributions.rst>

Installation
============

.. toctree::
  :maxdepth: 2

  How to install <dgfit/install.rst>

Repository
==========

GitHub: `DGFit <https://github.com/karllark/DGFit>`_

Quick Start
===========

Material needed.

Reporting Issues
================

If you have found a bug in ``DGFit`` please report it by creating a
new issue on the ``DGFit`` `GitHub issue tracker
<https://github.com/karllark/DGFit/issues>`_.

Please include an example that demonstrates the issue sufficiently so that the
developers can reproduce and fix the problem. You may also be asked to provide
information about your operating system and a full Python stack trace.  The
developers will walk you through obtaining a stack trace if it is necessary.

Contributing
============

Like the `Astropy`_ project, ``DGFit`` is made both by and for its
users.  We accept contributions at all levels, spanning the gamut from fixing a
typo in the documentation to developing a major new feature. We welcome
contributors who will abide by the `Python Software Foundation Code of Conduct
<https://www.python.org/psf/conduct/>`_.

``DGFit`` follows the same workflow and coding guidelines as
`Astropy`_.  The following pages will help you get started with contributing
fixes, code, or documentation (no git or GitHub experience necessary):

* `How to make a code contribution <https://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_

* `Coding Guidelines <https://docs.astropy.org/en/stable/development/codeguide.html>`_

* `Developer Documentation <https://docs.astropy.org/en/latest/index_dev.html#developer-docs>`_

For the complete list of contributors please see the `DGFit
contributors page on Github
<https://github.com/karllark/DGFit/graphs/contributors>`_.



Reference API
=============

.. automodapi:: dgfit.obsdata

.. automodapi:: dgfit.dustgrains

.. automodapi:: dgfit.dustmodel
    :inherited-members:
