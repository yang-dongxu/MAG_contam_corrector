import os
import sys
from Correct_checkm_error import PathNotExist


class Values(object):

    __PATH=['work_directory','marker_gene_exclude','select_map_file','storage_path','bins_path','marker_file','mark_gene_stats']
    """values to transport"""
    def __init__(self,options):
        self.work_directory=options.work_directory
        self.marker_gene_exclude_file=options.marker_gene_exclude_file
        self.select_map_file=options.select_map_file

        self.storage_path=os.path.join(options.work_directory,options.storage_path)
        self.bins_path=os.path.join(options.work_directory,options.bins_path)
        self.marker_file=os.path.join(options.work_directory,options.marker_file)
        self.mark_gene_stats=os.path.join(self.storage_path,options.marker_gene_stat)
        self.binids=self.__get_binids()
        keys=self.binids
        genes_faa=[os.path.join(self.bins_path,i,options.genes_faa) for i in keys]
        self.genes_faa={}
        for i in range(len(keys)):
            self.genes_faa[keys[i]]=genes_faa[i]
        self.__check_exist()

        self.evalue=options.evalue
        self.min_len=options.min_len
        self.identity_low=options.identity_low
        self.identity_high=options.identity_high
        self.web=options.web
        self.hitlist=options.hitlist
        self.db=options.db
        self.blast_threads=options.blast_threads
        self.select_by_direction = options.select_by_direction
        self.semi_partial=options.semi_partial
        self.overlap_len=options.overlap_len
        self.exact_same_len=options.exact_same_len
        self.hit_len_miss=options.hit_len_miss/100
        self.identity_diff=options.identity_diff

        self.out_path=options.out_path
        self.out_file=os.path.join(self.out_path,options.out_file)
        self.blast_out=os.path.join(self.out_path,options.blast_out)
        self.debug=options.debug
        self.log=os.path.join(self.out_path,options.log)
        self.binids_file=options.binids_file

        self.diamond=options.diamond


    def __get_binids(self):
        binids=os.listdir(self.bins_path)
        return binids

    def __check_exist(self):
        for i in dir(self):
            if '__' in i:
                continue
            elif i=='genes_faa':
                for t,j in getattr(self,i,{}).items():
                    if not os.path.exists(j):
                        raise (PathNotExist(j))
            elif i  in self.__PATH:
                if  not os.path.exists(getattr(self,i,i)):
                    raise(PathNotExist(getattr(self,i,i)))
        return  True


        


