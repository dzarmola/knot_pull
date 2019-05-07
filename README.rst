| |version| |versions||wheel|


.. |version| image:: http://img.shields.io/pypi/v/intspan.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/knot_pull

.. |versions| image:: https://img.shields.io/pypi/pyversions/intspan.svg
    :alt: Supported versions
    :target: https://pypi.org/project/knot_pull

.. |wheel| image:: https://img.shields.io/pypi/wheel/intspan.svg
    :alt: Wheel packaging support
    :target: https://pypi.org/project/knot_pull


KnotPull reduces a user provided 3D structure, to simplify it,
while preserving the topology of the chain. It has been successfully
used for knot detection in proteins and chromatine chains.

It can return either a final, simplified structure, or a .pdb formatted
trajectory of consecutive simplifications performed.

Finally, it can give a knot type assignment (including connected knots and links)
of the simplified structure (using the Dowker notation for topology detection).
