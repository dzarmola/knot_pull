<table>
    <tr>
        <td>License</td>
        <td><img src='https://img.shields.io/pypi/l/knot_pull.svg'></td>
        <td>Version</td>
        <td><img src='https://img.shields.io/pypi/v/knot_pull.svg'></td>
    </tr>
    <tr>
        <td>Travis CI</td>
        <td><img src='https://travis-ci.org/dzarmola/knot_pull.svg?branch=master'></td>
        <td>Coverage</td>
        <td><img src='https://codecov.io/gh/dzarmola/knot_pull/branch/master/graph/badge.svg'></td>
    </tr>
    <tr>
        <td>Status</td>
        <td><img src='https://img.shields.io/pypi/status/knot_pull.svg'></td>
        <td>Downloads</td>
        <td><img src='https://img.shields.io/pypi/dm/knot_pull.svg'></td>
    </tr>
    <tr>
        <td>Wheel</td>
        <td><img src='https://img.shields.io/pypi/wheel/knot_pull.svg'></td>
        <td>Supported versions</td>
        <td><img src='https://img.shields.io/pypi/pyversions/knot_pull.svg'></td>
    </tr>
</table>

# KnotPull - a simplifier for 3D structures

KnotPull reduces a user provided 3D structure, to simplify it,
while preserving the topology of the chain. It has been successfully
used for knot detection in proteins and chromatin chains.

It can return either a final, simplified structure, or a .pdb formatted
trajectory of consecutive simplifications performed.

## Install
`pip install knot_pull `

(older versions) `pip install --index-url https://test.pypi.org/simple/ knot_pull`

If installing with `pip --user` knot_pull_check and knot_pull_show
will be installed in your user base
(`python -c "import site; print(site.USER_BASE)"`), e.g.
`$HOME/.local/bin/knot_pull_check`

## Running

To run on a specific chain of a local pdb file:
```
knot_pull_check -o -k -c A 3ris.pdb
```
To run on a local xyz file (with format id<whitespace>x<whitespace>y<whitespace>z; chains separated by a line with "END"):
```
knot_pull_check -o from_xyz.pdb -k 05.xyz
```
To run on by structure id from RCSB PDB (without creating an output file):
```
knot_pull_check -k 4mcb
```

To show coloured trajectory in PyMOL (if installed):
```
knot_pull_show my_output_file.pdb
```
If PyMOL is not available in system path it can be specified as the second argument

### Help:
```
usage: knot_pull_check [-h] [-p] [-t] [-f {pdb,xyz,guess}] [-k] [-c CHAIN]
                       [-o [OUTPUT]] [-s [SAVE]]
                       [infile]

Simplify the 3D structure while preserving the topology

positional arguments:
  infile                Structure file which will be used

optional arguments:
  -h, --help            show this help message and exit
  -p, --preserve_resnums
                        Keep original numbering
  -t, --trajectory      Write out all steps of the simplification [default is
                        1]
  -f {pdb,xyz,guess}, --format {pdb,xyz,guess}
                        Input file format: xyz or pdb [default is guess]
  -k, --detect_knot     Calculate the knot Dowker notation and assign the knot
                        type
  -c CHAIN, --chain CHAIN
                        Specify which chain(s) you want to simplify. Multiple
                        chains should be separated by commas.
  -b BEGIN, --begin BEGIN
                        Smooth and analyze only part of the structure: specify
                        the first bead to be considered. '-b 4' will start
                        from the 4th bead, '-b i4' will start from the
                        residue/coordinate with id 4. When specified, only the
                        first chain from selected/present will be used.
  -e END, --end END     Smooth and analyze only part of the structure: specify
                        the first bead to be considered. '-e 104' will
                        disregard beads from 105th forward, '-e i104' will
                        disregard residue/coordinates with id >=105. When
                        specified, only the first chain from selected/present
                        will be used.
  -o [OUTPUT], --output [OUTPUT]
                        Specify the output file. Format is guessed based on
                        extension [default pdb]. For stdout use "-". Empty
                        argument makes an '_out.pdb' file. If this option is
                        missing just the last frame will be shown (sets
                        trajectory to 0)
  -q, --quiet           Nothing except the knot type (if calculated) will be
                        written to the screen


```

### Copyright (c) Aleksandra Jarmolińska
