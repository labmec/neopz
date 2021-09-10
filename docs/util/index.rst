Utilities
=========

.. toctree::
   :maxdepth: 1
   
.. contents:: Table of Contents
   :local:
  
      


      
TPZSimpleTimer
--------------

.. doxygenclass:: TPZSimpleTimer
   :members:

TPZLogger
---------

Logging in NeoPZ is currently done through the
`log4cxx <https://logging.apache.org/log4cxx/latest_stable/>`_ library, and its
linkage with NeoPZ can be configured through the ``USING_LOG4CXX`` CMake option.

A user familiarised with the log4cxx library should expect no trouble creating their
own custom logger, but NeoPZ has several loggers available with sane defaults. However,
one is able to easily change the verbosity levels (and adding their own loggers) through
the configuration file.

Configuring the logger for an application
+++++++++++++++++++++++++++++++++++++++++



The configuration of the loggers can be done through a call to

.. doxygenfunction:: TPZLogger::InitializePZLOG()

The default configuration file is ``$INSTALL_DIR/include/Util/log4cxx.cfg``
(or ``$BUILD_DIR/Util/log4cxx.cfg`` when inside the NeoPZ buildtree). This file is
generated based on the ``$SOURCE_DIR/Util/log4cxx.cfg``, which should **not** be modified
by the user.

.. note::

   If the ``log4cxx.cfg`` in the source directory is modified
   (by a user that does not care about this documentation,
   or by means of a ``git pull``),
   then it will be copied to both the build and install directories.
   Any existing configuration file will be renamed with a ``-backup`` suffix.

Alternatively, the user can provide a custom configuration file and call

.. doxygenfunction:: TPZLogger::InitializePZLOG(const std::string &configfile)

which is the recommended method for external projects, since each might have quite
distinct verbosity levels for their loggers.