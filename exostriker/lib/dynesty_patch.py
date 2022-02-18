#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities for handling results.

"""
import dynesty
   
class Results(dynesty.results.Results):


    def __setattr__(self,name, value):
        if name[0] != '_' and self._initialized:
            pass
       #     #raise RuntimeError("Cannot set attributes directly")
        super().__setattr__(name, value)
  
