
.. image:: http://img.shields.io/pypi/l/knot_pull.svg?style=flat
    :alt: License
    :target: https://pypi.org/project/knot_pull

.. image:: http://img.shields.io/pypi/v/knot_pull.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/knot_pull

..  image:: https://img.shields.io/pypi/pyversions/knot_pull.svg
    :alt: Supported versions
    :target: https://pypi.org/project/knot_pull

..  image:: https://travis-ci.org/dzarmola/knot_pull.svg?branch=master
    :alt: Travis CI build passing
    :target: https://pypi.org/project/knot_pull

..  image:: https://codecov.io/gh/dzarmola/knot_pull/branch/master/graph/badge.svg
    :alt: Code coverage
    :target: https://pypi.org/project/knot_pull

..  image:: https://img.shields.io/pypi/wheel/knot_pull.svg
    :alt: Wheel packaging support
    :target: https://pypi.org/project/knot_pull


KnotPull reduces a user provided 3D structure, to simplify it,
while preserving the topology of the chain. It has been successfully
used for knot detection in proteins and chromatine chains.

It can return either a final, simplified structure, or a .pdb formatted
trajectory of consecutive simplifications performed.

Finally, it can give a knot type assignment (including connected knots and links)
of the simplified structure (using the Dowker notation for topology detection).
