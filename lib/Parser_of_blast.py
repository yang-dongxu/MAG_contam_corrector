from Redundant_gene_parser import *
from Values import Values
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Dot_aligner import Dot_aligner

import os
import sys
import Correct_checkm_error 

class Parser_of_blast(object):
    """description of class"""
    def __init__(self,values:Values,redundant_gene_parser:Redundant_gene_parser):
        self.blast_out = values.blast_out
        self.binids = values.binids
        self.values = values
        self.binid_to_redundantGene_info = redundant_gene_parser.binid_to_redundantGene_info
        self.make_path_exist()
        self.binid_to_gene_to_selectResult={}

    def make_path_exist(self):
        for binid,genes in self.binid_to_redundantGene_info.items():
            for gene in genes:
                if os.path.exists(os.path.join(self.blast_out,binid,gene)):
                    continue
                else:
                    os.makedirs(os.path.join(self.blast_out,binid,gene))
        return True

    def blast(self):
        for binid, genes in self.binid_to_redundantGene_info.items():
            for gene,segs in genes.items():
                if len(segs.get_pairs()) :
                    segs.blast(os.path.join(self.blast_out,binid,gene))
                else:
                    pass
        return True

    def pairs_and_select(self):
        for binid, genes in self.binid_to_redundantGene_info.items():
            gene_to_result={}
            for gene, segs in genes.items():
                seg_pairs=Seg_pair_container()
                for p in segs.get_pairs():
                    f=Seg_for_cmp(p[0],p[0].blast_out_file)
                    r=Seg_for_cmp(p[1],p[1].blast_out_file)
                    seg_pair=Seg_pair(f,r,self.values)
                    if binid=="1278073.fasta_200_1" and gene=="PF01281.14" :
                        #print("nothing")
                        pass
                    seg_pairs.append(seg_pair)

                    #here is something used for debug
                    if binid=="1278073.fasta_200_1" and gene=="PF01281.14" :
                        pass
                gene_to_result[gene]=seg_pairs.get_edit_info()
            self.binid_to_gene_to_selectResult[binid]=gene_to_result
                    

class Seg_for_cmp(Seg_info):
    def __init__(self,seg:Seg_info,blast_out_file):
        super().__init__(id=seg.id,partial=seg.partial,direction=seg.direction,seq=seg.seq)
        self.blast_out_file = blast_out_file
        self.accession_hsp = {}
        self.parse_blast_xml()
        self.ranks=[]
        self.init_accession_rank()



    def parse_blast_xml(self):
        with open(self.blast_out_file) as f:
            handle = NCBIXML.read(f)
        for alignment in handle.alignments:
            self.accession_hsp[alignment.accession] = alignment.hsps[0]
        return self.accession_hsp

    def get_identity(self,accesion):
        a=self.accession_hsp[accesion].identities
        b=self.accession_hsp[accesion].align_length
        return (a/b)*100

    def init_accession_rank(self):
        for accession,hit in self.accession_hsp.items():
            if hit.bits:
                if hit.bits not in self.ranks:
                    self.ranks.append(hit.bits)
        self.ranks=sorted(self.ranks)[::-1]

    def get_accession_rank(self,accession):
        bit_score=self.accession_hsp[accession].bits
        return str(self.ranks.index(bit_score)+1)

    def __getitem__(self, key):
        return self.accession_hsp[key]


class Seg_pair(object):
    def __init__(self,f_seg:Seg_for_cmp,r_seg:Seg_for_cmp,values:Values):
        self.f_seg = f_seg
        self.r_seg = r_seg
        self.Flag=True
        self.accessions=[]
        self.values=values
        self.select()

    def select(self):
        self.select_by_accesion()
        if self.values.select_by_direction:
            self.select_by_direction()
        self.select_by_identity()
        self.select_by_identity_diff()
        self.select_by_lenCoverage()
        self.select_by_site()
             

    def select_by_direction(self):
        if  not self.f_seg.is_direction_consistent(self.r_seg):
            self.Flag=False
        pass

    def select_by_accesion(self):
        f_accesion=self.f_seg.accession_hsp.keys()
        r_accesion=self.r_seg.accession_hsp.keys()
        self.accessions=list(set(f_accesion) & set(r_accesion))

    def select_by_site(self):
        new_accesions=[]
        for accession in self.accessions:
            if self.r_seg.accession_hsp[accession].sbjct_end < self.f_seg.accession_hsp[accession].sbjct_start:
                new_accesions.append(accession)
            else:
                overlap=self.r_seg.accession_hsp[accession].sbjct_end - self.f_seg.accession_hsp[accession].sbjct_start 

                r_shift=len(self.r_seg.seq)+1-self.f_seg.accession_hsp[accession].query_end
                r_overlap_part=str(self.r_seg.seq)[self.f_seg.accession_hsp[accession].sbjct_start-self.r_seg.accession_hsp[accession].sbjct_end-1-r_shift:]

                f_shift=self.f_seg.accession_hsp[accession].query_start
                f_overlap_part=str(self.f_seg.seq)[:self.r_seg.accession_hsp[accession].sbjct_end-self.f_seg.accession_hsp[accession].sbjct_start+1+f_shift]
                exact_same_len=Dot_aligner(r_overlap_part,f_overlap_part).longest_hit_len
                if overlap < self.values.overlap_len and overlap > self.values.exact_same_len:
                    if exact_same_len >self.values.exact_same_len:
                        new_accesions.append(accession)
                elif overlap <= self.values.exact_same_len and overlap >self.values.exact_same_len/2:
                    if exact_same_len >self.values.exact_same_len/2:
                        new_accesions.append(accession)
                elif overlap <=self.values.exact_same_len/2:
                    new_accesions.append(accession)

                #below is the part for debug
            if  self.values.debug:
                print('!,%s,%s,%d,%d,%d,%d,%f,%f,%s,%s,%s\n'%(self.r_seg.id,self.f_seg.id,self.r_seg.accession_hsp[accession].sbjct_start,self.r_seg.accession_hsp[accession].sbjct_end,self.f_seg.accession_hsp[accession].sbjct_start,self.f_seg.accession_hsp[accession].sbjct_end,self.r_seg.get_identity(accession),self.f_seg.get_identity(accession),str(self.r_seg.seq),str(self.f_seg.seq),self.r_seg.accession_hsp[accession].match[self.f_seg.accession_hsp[accession].sbjct_start-self.r_seg.accession_hsp[accession].sbjct_end:]))
        self.accessions=new_accesions

    def select_by_identity(self):
        new_accession=[]
        for accession in self.accessions:
            if (self.f_seg.get_identity(accession)>=self.values.identity_low) and (self.f_seg.get_identity(accession)<=self.values.identity_high):
                if self.r_seg.get_identity(accession)>= self.values.identity_low and self.r_seg.get_identity(accession) <= self.values.identity_high:
                    new_accession.append(accession)
            else:
                pass
            
        self.accessions=new_accession

    def select_by_identity_diff(self):
        new_accession=[]
        for accession in self.accessions:
            if abs(self.f_seg.get_identity(accession) -self.r_seg.get_identity(accession))<self.values.identity_diff:
                new_accession.append(accession)
        self.accessions=new_accession

    def select_by_lenCoverage(self):
        new_accession=[]
        for accession in self.accessions:
            f_query_len=self.f_seg.accession_hsp[accession].query_end-self.f_seg.accession_hsp[accession].query_start+1
            r_query_len=self.r_seg.accession_hsp[accession].query_end-self.r_seg.accession_hsp[accession].query_start+1
            if (self.f_seg.length-f_query_len)/self.f_seg.length <self.values.hit_len_miss:
                if (self.r_seg.length-r_query_len)/self.r_seg.length <self.values.hit_len_miss:
                    new_accession.append(accession)
            if self.values.debug:
                print('!!,%s,%s,%s,%.4f,%.4f\n'%(self.r_seg.id,self.f_seg.id,accession,(self.r_seg.length-r_query_len)/self.r_seg.length,(self.r_seg.length-r_query_len)/self.r_seg.length))
        self.accessions=new_accession


    def is_redundant(self,seg_pair):
        assert isinstance(seg_pair,Seg_pair)
        if self.f_seg.id== seg_pair.f_seg.id or self.r_seg.id==seg_pair.r_seg.id:
            return True
        else:
            return False

    def __nonzero__(self):
        if len(self.accessions) and self.Flag:
            return True
        else :
            return False
    def is_not_empty(self):
        if len(self.accessions) and self.Flag:
            return True
        else :
            return False

    def output(self):
        highest_identity_accession=0
        highest_identity_sum=0
        for accession in self.accessions:
            if self.f_seg.get_identity(accession)+self.r_seg.get_identity(accession)>highest_identity_sum:
                highest_identity_sum=self.f_seg.get_identity(accession)+self.r_seg.get_identity(accession)
                highest_identity_accession=accession
        r_identity_1=str(self.r_seg.get_identity(highest_identity_accession))
        f_identity_1=str(self.f_seg.get_identity(highest_identity_accession))
        r_start_1=str(self.r_seg.accession_hsp[highest_identity_accession].sbjct_start)
        r_end_1=str(self.r_seg.accession_hsp[highest_identity_accession].sbjct_end)
        f_start_1=str(self.f_seg.accession_hsp[highest_identity_accession].sbjct_start)
        f_end_1=str(self.f_seg.accession_hsp[highest_identity_accession].sbjct_end)
        f_accession_index_1=self.f_seg.get_accession_rank(highest_identity_accession)
        r_accession_index_1=self.r_seg.get_accession_rank(highest_identity_accession)


        nearest_identity_accession=0
        nearest_identity_minus=100
        for accession in self.accessions:
            if abs(self.f_seg.get_identity(accession)-self.r_seg.get_identity(accession))<nearest_identity_minus:
                nearest_identity_accession=accession
                nearest_identity_minus=abs(self.f_seg.get_identity(accession)-self.r_seg.get_identity(accession))
        r_identity_2=str(self.r_seg.get_identity(nearest_identity_accession))
        f_identity_2=str(self.f_seg.get_identity(nearest_identity_accession))
        r_start_2=str(self.r_seg.accession_hsp[nearest_identity_accession].sbjct_start)
        r_end_2=str(self.r_seg.accession_hsp[nearest_identity_accession].sbjct_end)
        f_start_2=str(self.f_seg.accession_hsp[nearest_identity_accession].sbjct_start)
        f_end_2=str(self.f_seg.accession_hsp[nearest_identity_accession].sbjct_end)
        f_accession_index_2=self.f_seg.get_accession_rank(nearest_identity_accession)
        r_accession_index_2=self.r_seg.get_accession_rank(nearest_identity_accession)




        return (self.r_seg.get_id(),self.f_seg.get_id(),str(len(self.accessions)),highest_identity_accession,r_accession_index_1,f_accession_index_1,r_identity_1,f_identity_1,r_start_1,r_end_1,f_start_1,f_end_1,nearest_identity_accession,r_accession_index_2,f_accession_index_2,r_identity_2,f_identity_2,r_start_2,r_end_2,f_start_2,f_end_2)

class Seg_pair_container(object):
    def __init__(self):
        self.edit_info=[]
        self.seg_pairs=[]
        pass

    def append(self,seg_pair:Seg_pair):
        if  not seg_pair.is_not_empty():
            return False
        for eachone in self.seg_pairs:
            if eachone.is_redundant(seg_pair):
                return False
        self.seg_pairs.append(seg_pair)
        self.get_edit_info()
        return True


    def get_edit_info(self):
        self.edit_info=self.seg_pairs
        return self.seg_pairs

    def output(self):
        return [seg_pair.output() for seg_pair in self.seg_pairs]

        










