from Values import Values
import Correct_checkm_error 

class Marker_sets_parser(object):
    """a parser for marker gene sets"""
    def __init__(self,values:Values):
        self.file=values.marker_file
        self.binids=values.binids
        self.values=values
        self.type=self.classify_type()
        self.binid_to_marker_sets={}
        self.parse_marker_sets()

    def __getitem__(self,binid):
        return self.binid_to_marker_sets[binid]

    def __len__(self):
        return len(self.binid_to_marker_sets)

    def classify_type(self):
        with open(self.file) as f:
            head=f.readline()
        if head.strip()=='# [Taxon Marker File]':
            self.type=2
        elif head.strip()=='# [Lineage Marker File]':
            self.type=1
        else:
            raise(TypeError("type label of %s is not supported!"%self.file))
        return self.type

    def parse_marker_sets(self):
        if self.type==1:
            self.lineage_specific_parse_marker_sets()
        elif self.type==2:
            self.taxonomy_parse_marker_sets()
        return self.binid_to_marker_sets

    def taxonomy_parse_marker_sets(self):
        binid_to_marker_sets={}
        Marker_set.select_map_parser(self.values.select_map_file)
        Marker_set.exclude_genes(self.values.marker_gene_exclude_file)
        with open(self.file) as f:
            f.readline()
            marker_sets=Marker_set(self.type,f.readline())
        for binid in self.binids:
            binid_to_marker_sets[binid]=marker_sets
        self.binid_to_marker_sets=binid_to_marker_sets
        return self.binid_to_marker_sets

    
    def lineage_specific_parse_marker_sets(self):
        binid_to_marker_sets={}
        Marker_set.select_map_parser(self.values.select_map_file)
        Marker_set.exclude_genes(self.values.marker_gene_exclude_file)
        with open(self.file) as f:
            f.readline()
            for line in f:
                marker_sets=Marker_set(self.type,line)
                binid_to_marker_sets[marker_sets.binid]=marker_sets
        for binid in self.binids:
            if binid not in binid_to_marker_sets:
                raise(Correct_checkm_error.BinidInconsistent(binid))
        self.binid_to_marker_sets=binid_to_marker_sets
        return self.binid_to_marker_sets


class Marker_set(object):

    @classmethod
    def select_map_parser(cls,map_name):
        selectedMarkerSetMap = {}
        for line in open(map_name):
            lineSplit = line.split('\t')
            internalID = lineSplit[0]
            selectedID = lineSplit[1].rstrip()
            selectedMarkerSetMap[internalID] = selectedID
        cls.select_map=selectedMarkerSetMap
        return selectedMarkerSetMap

    @classmethod
    def exclude_genes(cls,exclude_file_name):
        with open(exclude_file_name) as f:
            cls.exclude=eval(f.read())
        return cls.exclude

    def __init__(self,type:int,line:str):
        self.type=type
        self.raw_data=line
        self.get_marker_set()

    def __iter__(self):
        return self.marker_sets

    def __len__(self):
        return len(self.marker_sets)

    def get_marker_set(self):
        if self.type==1:
            return self.lineage_marker_sets()
        elif self.type==2:
            return self.taxon_marker_sets()
        else:
            raise(TypeError("your type code has to be 1(lineage_specific) or 2(taxonomy)\n"))

    def taxon_marker_sets(self):
        """Construct bin marker set data from line."""
        lineSplit = self.raw_data.split('\t')
        self.uid=lineSplit[2]
        self.lineageStr=lineSplit[3]
        self.numGenomes=int(lineSplit[4])
        self.marker_sets=eval(lineSplit[5])
        return self.marker_sets

    def lineage_marker_sets(self):
        """Construct bin marker set data from line."""
        total_info=[]
        lineSplit = self.raw_data.split('\t')
        self.binid=lineSplit[0]
        numMarkerSets = int(lineSplit[1])
        for i in range(0, numMarkerSets):
            uid = lineSplit[i * 4 + 2]
            lineageStr = lineSplit[i * 4 + 3]
            numGenomes = int(lineSplit[i * 4 + 4])
            markerSet = eval(lineSplit[i * 4 + 5])
            total_info.append((uid, lineageStr, numGenomes, markerSet))
        raw_uid=total_info[0][0]
        self.uid=self.select_map[raw_uid]

        #筛选出对应的uid和marker sets
        for info in total_info:
            if self.uid==info[0]:
                self.lineageStr=info[1]
                self.numGenomes=info[2]
                self.marker_sets=info[3]
                break
        self.exclude_useless()
        return self.marker_sets

    def exclude_useless(self):
        clean_sets=[]
        for aset in self.marker_sets:
            p=aset-self.exclude
            if len(p) !=0:
                clean_sets.append(p)
        self.marker_sets=clean_sets
        return self.marker_sets

    def get_set_len(self,gene):
        for aset in self.marker_sets:
            if gene in aset:
                return len(aset)
        raise(Correct_checkm_error.MarkerGeneInconsistent(gene))



