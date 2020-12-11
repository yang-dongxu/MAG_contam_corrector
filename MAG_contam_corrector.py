import sys
import os

sys.path.append(os.path.join(os.getcwd(),'lib'))

from optparse import OptionParser
import Correct_checkm_error
from Values import Values
from Markersets_parser import Marker_sets_parser
from Sample import Sample
from Redundant_gene_parser import Redundant_gene_parser
from Parser_of_blast import Parser_of_blast

if __name__ == '__main__':

    VERSIOM='0.0.1'

    parser = OptionParser("usage:-D checkm_result_directory -m marker_gene_file_name -o out_path -f out_file_name\n")

    #输入选项
    parser.add_option('-D','--work_directory',action='store',type="string",dest='work_directory',help="this is the directory of checkm result\n")
    #parser.add_option('-s','--summary_file',action='store',type='string',dest='summary_file',help="the
    #summary result of checkm");
    parser.add_option('-B','--bins_path',action='store',type='string',dest='bins_path',default='bins',help="this is the name of checkm result bins\n")
    parser.add_option('-S','--storage_path',action='store',type='string',dest='storage_path',default='storage',help="this is the path of checkm result storage\n")
    parser.add_option('-m','--marker_file',action='store',type='string',dest='marker_file',default='lineage.ms',help="the marker file of checkm\n")
    parser.add_option('-s','--marker_gene_stat',action="store",type='string',dest='marker_gene_stat',default='marker_gene_stats.tsv',help="this is name of checkm storage marker_gene_stat file, which contains marker gene detail info in each bin\n")
    parser.add_option('-g','--genes_faa',action='store',type='string',dest='genes_faa',default='genes.faa',help="this is the fasta contains protein seq, default is the bins/BinID/genes.faa\n")
    parser.add_option('--selected_marker_sets',type='string',dest="select_map_file",default=os.path.join('lib','selected_marker_sets.tsv'),help="this is the map file for lineage_specific mode ,usually don't need to change.\n")
    parser.add_option('--marker_gene_exclude',type='string',dest="marker_gene_exclude_file",default=os.path.join('lib','marker_gene_to_exclude.dat'),help="this file contains the marker gene don't need to calculate.\n")

    #筛选条件相关
    parser.add_option('-e','--evalue',action='store',type='float',dest='evalue',default=0.00001,help="evalue of blast\n")
    parser.add_option('-l','--min_len',action='store',type='int',dest='min_len',default=20,help="the min length of seq to blast \n")
    parser.add_option('-i','--identity',action='store',type='float',dest='identity_low',default=50,help="the identity lower bound\n")
    parser.add_option('--identity_high',action='store',type='float',dest='identity_high',default=100,help="the identity higher bound\n")
    parser.add_option('--hitlist',action='store',type='int',dest='hitlist',default=100,help='the length of blast result list\n')
    parser.add_option('-w','--web_blast',action='store_true',dest='web',default=False,help='if use web blast instead of blast on your computer\n')
    parser.add_option('--db',action='store',type='string',dest='db',default='nr',help='the default blast database\n')
    parser.add_option('--blast_threads',action='store',type='int',dest='blast_threads',default=1,help='the threads used in blast when you choose to do it on your local computer\n')
    parser.add_option('--select_by_direction',action='store_true',dest='select_by_direction',default=False,help='if you want select by direction of gene strand, choose this, but it\'s recommended that run without it.\n')
    parser.add_option('--semi',action='store_true',dest='semi_partial',help='when select, the gene has to be semi-partial so that it can be a part of mis-count contaminant gene, which is predicted by prodigal of 10 or 01')
    parser.add_option('--overlap_len',action='store',type='int',dest='overlap_len',default=50,help='this is the length which is  the longest situation allowed to tolerent overlap between twp incomplete protein seqs, so that they could been seen as a broken one. the default is 50aa, from the length of sequencing reads. \n')
    parser.add_option('--exact_same_len',action='store',type='int',dest='exact_same_len',default=20,help='this is set for determine the overlap seqs whether similiar enough. \n')
    parser.add_option('--hit_len_miss',action='store',type='int',dest='hit_len_miss',default=25,help='the value is set to avoid the situation that query-hit sequence is too smaller than your raw data\n')
    parser.add_option('--identity_diff',action='store',type='float',dest='identity_diff',default=18,help='the max difference you allow between the identities of your seg pair\n')


    #输出相关选项
    parser.add_option('-o','--out_path',action='store',type='string',dest='out_path',default='out',help='the out path of result\n')
    parser.add_option('-f','--out_file',action='store',type="string",dest='out_file',default='result.csv',help="the name of output file\n")
    parser.add_option('-b','--blast_out',action='store',type='string',dest='blast_out',default='blast_out',help="the out path of blast in the out dir\n")


    #其他选项
    parser.add_option('--version',action='store_true',dest='version',default=False,help="show version\n")
    parser.add_option('--debug',action='store_true',dest='debug',default=False,help='debug mode\n')
    parser.add_option('--log_file',action='store',type='string',dest='log',default='log',help='the file name of log\n')
    parser.add_option('--binids_file',action='store',type='string',dest='binids_file',default='',help='this is the file contains the binid which you want to correct\n')

    parser.add_option('--diamond',action='store_true',dest='diamond',default=False,help='use diamond instead of blast+ if you choose this,but hasn\'t been supported now! \n')

    #解析参数
    options,args = parser.parse_args(sys.argv)

    #检测帮助、版本
    if options.version:
        print("version: %s\n"%VERSIOM)
        sys.exit(0)
    


values=Values(options)
ms=Marker_sets_parser(values)
sample=Sample(values,ms)
p=Redundant_gene_parser(values,sample.binid_to_redundantGene)
bp=Parser_of_blast(values,p)
bp.blast()
bp.pairs_and_select()
z=sample.output(bp)


print('success!')

