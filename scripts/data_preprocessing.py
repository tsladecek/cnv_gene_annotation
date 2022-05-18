#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data Preprocessing
"""

import numpy as np
from constants import IMPUTE_VALUE, SCALING
from constants import BINARIZE, COUNT, DUMMY, DBNSFP_NUMERICAL
from copy import copy

class cat:
    def binarize(col):
        """Is missing (0) or not (1)

        :param col: pd.dataframe column
        return np.array
        """
        return np.array([i is not np.nan for i in col]) * 1

    def count(col, sep=';'):
        """Count values in a cell. NaN = 0

        :param col: pd.dataframe column
        :return np.array
        """
        col = col.replace(np.nan, '')
        return num.scale(np.array([len(i.strip(sep)) for i in col]), how='minmax')

    def dummy(col):
        """Create dummy variables for categorical variables. Rows with 0 only are NaN

        :param col: pd.dataframe column
        :return tuple: (dummy np.array, column names)
        """

        classes = np.unique(col.dropna())

        # we dummy it for all classes, since we want to preserve info about missing values
        dummied_array = np.empty((len(col), len(classes)))

        for i, c in enumerate(classes):
            dummied_array[:, i] = np.array((col == c) * 1)
        return dummied_array, [f'{col.name}_{c}' for c in classes]

class num:
    def impute(col, value="mean"):
        """
        Impute values with a predefined function or a value
        
        :param col: pd.dataframe column
        :param value: if float, then NaN are replaced with this value.
                      Other options include: "mean", "median", "min", "max"
        """
        
        if type(value) is float or type(value) is int:
            return col.replace(np.nan, value)
        
        elif value == "mean":
            return col.replace(np.nan, np.mean(col))
        elif value == "median":
            return col.replace(np.nan, np.median(col).dropna())
        elif value == "min":
            return col.replace(np.nan, np.min(col))
        elif value == "max":
            return col.replace(np.nan, np.max(col))
        else:
            print("Unknown value/function")
            return
    
    def scale(col, how='standardize', _mean=None, _std=None, _min=None, _max=None):
        """What type of scaling is desired for numerical variables
        
        :param how: ['standardize', 'mean_norm', 'minmax']
        
        """
        if how == 'standardize':
            if _mean is not None and _std is not None:
                return (col - _mean) / _std
            else:
                return (col - np.mean(col)) / np.std(col)
        
        elif how == 'mean_norm':
            if _mean is not None and _min is not None and _max is not None:
                return (col - _mean) / (_max - _min)
            else:
                return (col - np.mean(col)) / (np.max(col) - np.min(col))
        
        elif how == 'minmax':
            if _min is not None and _max is not None:
                return (col - _min) / (_max - _min)
            else:
                return (col - np.min(col)) / (np.max(col) - np.min(col))
        
        else:
            print("Unknown scaling method.")
            print("Please choose one of: ['standardize', 'mean_norm', 'minmax']")
            return

def preprocess(df0,
               categorical=True,
               scale_numerical=True,
               impute_numerical=True):
    """
    Preprocess dataframe
    """
    
    df = copy(df0)
    
    # CATEGORICAL PREPROCESSING  
    if categorical:

        # BINARIZE
        for b in BINARIZE:
            df.loc[:, b] = cat.binarize(df.loc[:, b])
        
        # COUNT
        for c in COUNT:
            df.loc[:, c] = cat.count(df.loc[:, c])
        
        # DUMMY
        for d in DUMMY:
            col = df.loc[:, d]
            df.drop(d, axis=1, inplace=True)
            
            dummied, colnames = cat.dummy(col)
            for i in range(len(colnames)):
                df[colnames[i]] = dummied[:, i]
                
    # NUMERICAL PREPROCESSING    
    if scale_numerical and impute_numerical:
        for n in DBNSFP_NUMERICAL:
            df.loc[:, n] = num.scale(num.impute(df.loc[:, n], IMPUTE_VALUE), how=SCALING)
    
    elif scale_numerical and not impute_numerical:
        for n in DBNSFP_NUMERICAL:
            df.loc[:, n] = num.scale(df.loc[:, n], how=SCALING)

    elif not scale_numerical and impute_numerical:
        for n in DBNSFP_NUMERICAL:
            df.loc[:, n] = num.impute(df.loc[:, n], IMPUTE_VALUE)
    
    return df

