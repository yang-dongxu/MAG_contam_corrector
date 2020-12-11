from Values import Values
import  Correct_checkm_error  
from Bio import SeqIO,Seq
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import time

class Redundant_gene_parser(object):
    """this is the parser of redundant genes,  further process to prepare for blast and correct"""
    def __init__(self,values:Values,binid_to_redundantGene:dict):
        self.values=values
        self.genes_faa=values.genes_faa
        self.binid_to_redundantGene=binid_to_redundantGene
        self.binid_to_redundantGene_info={}
        self.get_genes_info()

    def get_genes_info(self):
        self.binid_to_redundantGene_info={}
        binid_to_genes_faa={}
        for binid,file_name in self.genes_faa.items():
            segid_to_info={}
            seqs=SeqIO.parse(file_name,'fasta')
            for seq in seqs:
                seq_id=seq.id.strip()
                direction=seq.description.split('#')[3]
                partial=seq.description.split('#')[4].split(';')[1].split('=')[1]
                segid_to_info[seq_id]=Seg_info(id=seq_id,direction=direction,partial=partial,seq=seq.seq)
            binid_to_genes_faa[binid]=segid_to_info
        for binid,genes in self.binid_to_redundantGene.items():
            id_info={}
            for gene,segs in genes.items():
                segs_info=Segs(self.values)
                for i in segs:
                    try:
                        if '&&' in i:
                            segs_info.append(binid_to_genes_faa[binid][i.split('&&')[0]])
                            segs_info.append(binid_to_genes_faa[binid][i.split('&&')[1]])
                        else:
                            segs_info.append(binid_to_genes_faa[binid][i])
                    except Exception:
                        print(binid)
                        print('\n')
                        print(i)
                        raise KeyError(Exception)
                segs_info.classify()
                if segs_info:
                    id_info[gene]=segs_info
            self.binid_to_redundantGene_info[binid]=id_info
        return self.binid_to_redundantGene_info



class Seg_info(object):
    def __init__(self,id:str,direction,partial:str,seq:Seq):
        '''this contains information of a gene in genes_faa,id is the first part,while direction is True when it is plus strand'''
        self.id=id
        self.partial=partial
        self.seq=seq
        self.length=len(self.seq)
        if type(direction) ==bool:
            self.direction=direction
        else:
            self.direction=bool(int(direction)+1)
        self.blast_out_file=''

    def __len__(self):
        return self.length

    def is_semipartial(self):
        if self.partial=='00':
            return False
        elif self.partial=='11':
            return False
        elif self.partial=='10' or self.partial=='01':
            return True
        else:
            raise(TypeError("partial of %s is %s, which is not supported type"%(self.id,self.partial)))

    def is_partial(self):
        if self.partial=='00':
            return False
        elif self.partial=='11':
            return True
        elif self.partial=='10' or self.partial=='01':
            return True
        else:
            raise(TypeError("partial of %s is %s, which is not supported type"%(self.id,self.partial)))


    def is_l_incomplete_seg(self):
        '''is True when the gene in the seg is incomplete in the left side of seg'''
        if self.partial=='10' or self.partial=='11':
            return True
        else:
            return False

    def is_r_incomplete_seg(self):
        '''is True when the gene in the seg is incomplete in the right side of seg'''
        if self.partial=='01' or self.partial=='11':
            return True
        else:
            return False
    
    def is_direction_consistent(self,seg):
        return not( self.direction ^ seg.direction)

    def is_strand_plus(self):
        return self.direction

    def is_f_incomplete_gene_seg(self):
        '''is True only when the gene on this seg is incomplete in the front,for example its partial is 01, which means it is  a l_incomplete_seg, while it's a minus-strand, so the gene is incomplete in the front, so the seg has a f_incomplete_gene'''
        return (self.is_r_incomplete_seg() ^self.is_strand_plus())

    def set_blast_out_file(self,file_name):
        self.blast_out_file=file_name

    def get_id(self):
        return self.id


class Segs(object):
    def __init__(self,values:Values):
        self.segs=[]
        self.min_len=values.min_len
        self.values=values

    def append(self,seg:Seg_info):
        if self.values.semi_partial:
            if seg.is_semipartial() and seg.length>self.min_len:
                self.segs.append(seg)
                flag = True
            else :
                flag = False
            return flag
        else:
            if seg.is_partial() and seg.length>self.min_len:
                self.segs.append(seg)
                flag = True
            else :
                flag = False
            return flag

    def classify(self):
        self.f_incomplete_gene_segs=[]
        self.r_incomplete_gene_segs=[]
        for seg in self.segs:
            if seg.is_partial() and not seg.is_semipartial():
                self.f_incomplete_gene_segs.append(seg)
                self.r_incomplete_gene_segs.append(seg)
            if seg.is_f_incomplete_gene_seg():
                self.f_incomplete_gene_segs.append(seg)
            else:
                self.r_incomplete_gene_segs.append(seg)

        return self.f_incomplete_gene_segs,self.r_incomplete_gene_segs

    def __nonzero__(self):
        self.classify()
        if len(self.f_incomplete_gene_segs) and len(self.r_incomplete_gene_segs):
            return True
        else :
            return False

    def is_empty(self):
        self.classify()
        if len(self.f_incomplete_gene_segs) and len(self.r_incomplete_gene_segs):
            return False
        else :
            return True

    def blast(self,path):
        for seg_info in self.r_incomplete_gene_segs:
            file_name=self.__blast(path=path,seg=seg_info)
            seg_info.set_blast_out_file(file_name)
        for seg_info in self.f_incomplete_gene_segs:
            file_name=self.__blast(path=path,seg=seg_info)
            seg_info.set_blast_out_file(file_name)


    def __blast(self,path,seg:Seg_info):
        if os.path.exists(os.path.join(path,seg.id)) and os.path.getsize(os.path.join(path,seg.id)):
            print('%s hasn\'t been blasted because it was exsiting! '%(os.path.join(path,seg.id))) 
            return os.path.join(path,seg.id)
        if self.values.web:
            blast_handler=qblast('blastp',database=self.values.db,sequence=seg.seq,expect=self.values.evalue,hitlist_size=self.values.hitlist)
            with open(os.path.join(path,seg.id),'w') as f:
                f.write(blast_handler.read())
        else:
            file_name=os.path.join(path,seg.id)
            fasta_file_name=file_name+'.fasta'
            f=open(fasta_file_name,'w')
            f.write('>%s\n'%seg.id)
            f.write(str(seg.seq))
            f.write('\n')
            f.close()
            time.sleep(0.1)
            if os.path.exists(file_name):
                os.remove(file_name)
            if not self.values.diamond:
                blastp_cline=NcbiblastpCommandline(query=fasta_file_name,db=self.values.db,evalue=self.values.evalue,outfmt=5,out=file_name,num_alignments=self.values.hitlist,num_threads=self.values.blast_threads)
                blastp_cline()
            else:
                cmd="diamond blastp -d %s -q %s -f 5 -o %s -e %f -k %d -p %d --sensitive"%(self.values.db,fasta_file_name,file_name,self.values.evalue,self.values.hitlist,self.values.blast_threads)
                os.system(cmd)
            os.remove(fasta_file_name)
        print('%s is got\n'%os.path.join(path,seg.id))
        return os.path.join(path,seg.id)

    def get_pairs(self):
        pairs=[]
        for f in self.f_incomplete_gene_segs:
            for r in self.r_incomplete_gene_segs:
                if f.id!=r.id:
                    pairs.append((f,r))
        return pairs

