Python API
==========

To use the plugin from your Python script, you can do:

.. code-block:: python

    import openmm as mm
    import openmmlab as mmlab
    system = mm.System()
    force = mmlab.SlicedNonbondedForce(2)
    system.addForce(force)

This is the implemented subclass of :OpenMM:`Force`:

.. toctree::
    :titlesonly:

    SlicedNonbondedForce


.. testsetup::

    from openmmlab import *
