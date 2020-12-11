from Values import Values
import Correct_checkm_error 
from Markersets_parser import Marker_sets_parser
from Parser_of_blast import Parser_of_blast
from copy import deepcopy

class Sample(object):
    """info of sample"""
    def __init__(self,values:Values,marker_sets:Marker_sets_parser):

        self.values=values
        self.binids_need_to_correct=[]

        self.marker_gene_stats=values.mark_gene_stats
        self.binids=values.binids
        self.marker_sets=marker_sets
        self.binid_to_rawInfo={}
        self.get_rawInfo()
        self.binid_to_gene_to_seg={}
        self.get_binid_to_gene_to_seg()
        self.completeness={}
        self.get_completeness()
        self.contamination={}
        self.get_contamination()
        self.binid_to_redundantGene={}
        self.get_binid_to_redundantGene()



    def get_rawInfo(self):
        if len(self.values.binids_file):
            f=open(self.values.binids_file)
            for line in f:
                self.binids_need_to_correct.append(line.strip())
            f.close()

        for line in open(self.marker_gene_stats):
            binid=line.split('\t')[0]
            info=eval(line.split('\t')[1])
            if len(self.values.binids_file):
                if binid in self.binids_need_to_correct:
                    self.binid_to_rawInfo[binid]=info
            else:
                self.binid_to_rawInfo[binid]=info
        for i in self.binid_to_rawInfo:
            if i not in self.binids:
                raise(Correct_checkm_error.BinidInconsistent(i))

        return self.binid_to_rawInfo

    def get_binid_to_gene_to_seg(self):
        for binid,rawInfo in self.binid_to_rawInfo.items():
            seg_to_gene={}
            for seg,genes_info in rawInfo.items():
                for gene in genes_info:
                    if gene in seg_to_gene:
                        seg_to_gene[gene].append(seg)
                    else:
                        seg_to_gene[gene]=[seg]
            self.binid_to_gene_to_seg[binid]=seg_to_gene
        return self.binid_to_gene_to_seg

    def get_completeness(self):
        for binid,genes_info in self.binid_to_gene_to_seg.items():
            comp=0
            for gene,segs in genes_info.items():
                set_length=self.marker_sets[binid].get_set_len(gene)
                sets_num=len(self.marker_sets[binid])
                comp+=(1.0/set_length)/sets_num
            self.completeness[binid]=comp
        return deepcopy(self.completeness)

    def get_contamination(self):
        for binid,genes_info in self.binid_to_gene_to_seg.items():
            cont=0
            for gene,segs in genes_info.items():
                set_length=self.marker_sets[binid].get_set_len(gene)
                sets_num=len(self.marker_sets[binid])
                cont+=max(0,len(segs)-1.0)/set_length/sets_num
            self.contamination[binid]=cont
        return deepcopy(self.contamination)

    def get_binid_to_redundantGene(self):
        self.binid_to_redundantGene={}
        for binid,genes_info in self.binid_to_gene_to_seg.items():
            redundant={}
            for gene , segs in genes_info.items():
                if len(segs)>1:
                    redundant[gene]=segs
            self.binid_to_redundantGene[binid]=redundant
        return deepcopy(self.binid_to_redundantGene)

    def __edit(self,change_info:Parser_of_blast):
        binid_to_gene_to_selectResult=change_info.binid_to_gene_to_selectResult
        for binid,gene_to_selectResult in binid_to_gene_to_selectResult.items():
            for gene, selectResult in gene_to_selectResult.items():
                for i in range(abs(len(selectResult))):
                    self.binid_to_gene_to_seg[binid][gene].pop()

    def output(self,change_info:Parser_of_blast):
        old_binid_to_gene_to_seg=deepcopy(self.binid_to_gene_to_seg)
        old_completeness=deepcopy(self.get_completeness())
        old_contamination=deepcopy(self.get_contamination())
        old_binid_to_redundantGene=deepcopy(self.binid_to_redundantGene)
        self.__edit(change_info)
        new_completeness=self.get_completeness()
        new_contamination=self.get_contamination()
        new_binid_to_redundantGene=deepcopy(self.binid_to_redundantGene)
        self.binid_to_gene_to_seg=old_binid_to_gene_to_seg
        with open(self.values.out_file,'w',encoding='utf8') as f:
            head=['binid','raw_completeness','raw_contamination','correct_completeness','correct_contamination','']
            f.write(','.join(head))
            f.write('\n')
            for binid in old_completeness:
                data=[old_completeness[binid],old_contamination[binid],new_completeness[binid],new_contamination[binid]]
                f.write('%s,'%binid)
                for word in data:
                    f.write('%.2f,'%(word*100))
                f.write('\n')

        if self.values.debug:
            f=open(self.values.log+'_repair','w',encoding='utf8')
            header=['binid','gene_name','r_seg_id','f_seg_id','same_accessions','highest_identity_accession','r_accession_index_1','f_accession_index_1','r_identity_1','f_identity_1','r_start_1','r_end_1','f_start_1','f_end_1','nearest_identity_accession','r_accession_index_2','f_accession_index_2','r_identity_2','f_identity_2','r_start_2','r_end_2','f_start_2','f_end_2']
            f.write(','.join(header))
            f.write('\n')
            binid_to_gene_to_selectResult=change_info.binid_to_gene_to_selectResult
            for binid,gene_to_selectResult in binid_to_gene_to_selectResult.items():
                for gene, selectResult in gene_to_selectResult.items():
                    for eachone in selectResult:
                        p=eachone.output()
                        f.write('%s,%s,%s\n'%(binid,gene,','.join(p)))
            f.close()

            f=open(self.values.log+'_raw','w',encoding='utf8')
            for binid , redundant in old_binid_to_redundantGene.items():
                for gene,segs in redundant.items():
                    f.write('>%s|%s\n'%(binid,gene))
                    f.write(str(segs))
                    f.write('\n')
            f.close()

            f=open(self.values.log+'_repairResult','w',encoding='utf8')
            for binid , redundant in new_binid_to_redundantGene.items():
                for gene,segs in redundant.items():
                    f.write('>%s|%s\n'%(binid,gene))
                    f.write(str(segs))
                    f.write('\n')
            f.close()


        return True




