# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:33:58 2020

@author: Dylan Agius
"""

import setuptools


setuptools.setup(
    name="lengmorph-DylanAgius", 
    version="0.0.1",
    author="Dylan Agius, Solid Mechanics Research Group, Univeristy of Bristol",
    author_email="dylan.j.agius@gmail.com",
    description="Create binary files containing arrays required in the implementation of a length scale modification for crystal plasticity simulations",
    long_description_content_type="text/markdown",
    url="https://github.com/DylanAgius/LengMorph",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)