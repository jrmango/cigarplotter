#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 19:20:27 2022

@author: cmr1
"""

import re

matches = 0
cigartally = 0
mismatches = 1
identity = 0

match = (re.findall(r'\d+[A-Z]', "100S84M16S"))
if match != None:
    for element in match:
        print(element)
        if element[-1] == "M":
            print("matches")
            cigartally += int(element[:-1])
        else:
            mismatches += int(element[:-1])
#code lifted straight from linked resource "on calculating percent identity from SAM files"
identity = 100 * (matches / (matches + mismatches))

print(identity)
print(cigartally)