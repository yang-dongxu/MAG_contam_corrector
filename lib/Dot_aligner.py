import numpy as np

class Dot_aligner(object):
    """this class is used for dot align """

    def __init__(self,f_seg,r_seg):
        self.f_seg=str(f_seg)
        self.r_seg=str(r_seg)

        self.score_matrix=np.zeros(shape=(len(self.f_seg)+1,len(self.r_seg)+1),dtype=int)
        self.dot_score_matrix=np.zeros(shape=(len(self.f_seg)+1,len(self.r_seg)+1),dtype=int)
        self.longest_hit_len=0

        self.dot_align()
        self.calcul()
        self.get_longest_len()

    def dot_align(self):
        for i in range(1,len(self.f_seg)+1):
            for j in range(1,len(self.r_seg)+1):
                if self.f_seg[i-1]==self.r_seg[j-1]:
                    self.score_matrix[i,j]=1
                else:
                    self.score_matrix[i,j]=0
    
    def calcul(self):
        
        for i in range(1,len(self.f_seg)+1):
            for j in range(1,len(self.r_seg)+1):
                if self.score_matrix[i,j]==0:
                    self.dot_score_matrix[i,j]=0
                else:
                    self.dot_score_matrix[i,j]=self.dot_score_matrix[i-1,j-1]+1
    
    def get_longest_len(self):
        self.longest_hit_len=self.dot_score_matrix.max()
        return self.longest_hit_len



