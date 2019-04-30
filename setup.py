#Copyright (c) Aleksandra Jarmolinska 2018

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="knot_pull",
    version="0.1.1",
    author="Aleksandra Jarmolinska",
    author_email="dzarmola@gmail.com",
    description="Simplifier of 3D structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dzarmola/knot_pull",
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    keywords='protein chromatine topology knot dowker',
    install_requires=[
          'numpy',
      ],
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=['bin/knot_pull_check'],
    include_package_data=True
)
