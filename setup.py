#Copyright (c) Aleksandra Jarmolinska 2018

import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

VERSION = "0.2.5"

setuptools.setup(
    name="knot_pull",
    version=VERSION,
    author="Aleksandra Jarmolinska",
    author_email="dzarmola@gmail.com",
    description="Simplifier of 3D structures of molecules",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/dzarmola/knot_pull",
    packages=setuptools.find_packages(),
    python_requires='>=2.7.*',
    classifiers=[
        'Development Status :: 4 - Beta',
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    keywords='protein chromatine topology knot dowker',
    install_requires=[
          'numpy',
      ],
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=['bin/knot_pull_check','bin/knot_pull_show'],
    include_package_data=True
)
