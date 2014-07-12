.. _getting_started:


***************
Getting started
***************

.. _installing-docdir:

Installing your doc directory
=============================

You may already have sphinx `sphinx <http://sphinx.pocoo.org/>`_
installed -- you can check by doing::

  python -c 'import sphinx'

If that fails grab the latest version of and install it with::

  > sudo easy_install -U Sphinx

Now you are ready to build a template for your docs, using
sphinx-quickstart::

  > sphinx-quickstart

accepting most of the defaults.  I choose "sampledoc" as the name of my
project.  cd into your new directory and check the contents::
