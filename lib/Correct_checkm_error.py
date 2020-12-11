import os

class PathNotExist(Exception):
    """error of path not exist"""
    def __init__(self,errorfile:str):
        super().__init__(self)
        self.ErrorFile=errorfile

    def __str__(self):
        return "%s is not exist!"%self.ErrorFile

class BinidInconsistent(Exception):
    '''error of binid inconsistent between bins and marker set files'''
    def __init__(self,errorfile:str):
        super().__init__(self)
        self.ErrorFile=errorfile

    def __str__(self):
        return "%s exist in bins but not in stat"%self.ErrorFile


class MarkerGeneInconsistent(Exception):
    def __init__(self,errorgene:str):
        super().__init__(self)
        self.ErrorGene=errorgene

    def __str__(self):
        return "gene %s not exist in the marker gene file! Is the file is right?"%self.ErrorGene

