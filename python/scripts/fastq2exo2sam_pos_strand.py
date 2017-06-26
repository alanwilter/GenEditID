#!/usr/bin/python3
import os
import re
import argparse


########################## START OF FUNCTIONS ################################

# print sam header

def print_sam_header( sam_header_file, samp_name, out_file ):
   sam_header = open( sam_header_file, 'r')
   for sam_h_line in sam_header:
      #sam_h_line=sam_h_line.rstrip()
      out_file.write( sam_h_line )
   #print('@RG', '\t', 'ID:',samp_name, '\t',  'LB:',samp_name, '\t', 'SM:',samp_name, '\t', 'PL:ILLUMINA')
   rg_line ="{}\t{}\t{}\t{}\t{}\n".format( '@RG', 'ID:'+ samp_name,  'LB:'+ samp_name, 'SM:' + samp_name,  'PL:ILLUMINA')
   out_file.write( rg_line)
   #print( '@PG', '\t', 'ID:exonorate2sam', '\t', 'VN:1', '\t', 'CL:python3 fastq2exo2sam.py', '-', '\t', 'PN:exonorate2sam')
   pg_line = "{}\t{}\t{}\t{}\t{}\n".format( '@PG', 'ID:exonorate2sam', 'VN:1', 'CL:python3 fastq2exo2sam.py', 'PN:exonorate2sam')
   out_file.write( pg_line)
   sam_header.close()
   
def make_fasta_and_get_info( file_name, lines ):
   temp_fa= open( file_name, 'w')
   fa_header = '>' + lines[0].split( sep=' ' )[0]
   fa_seq = lines[1]
   fa_qal = lines[3]
   read_length = len(fa_seq)
   temp_fa.write( fa_header + '\n' + fa_seq )
   temp_fa.close()
   return[ read_length, fa_seq, fa_qal ]


def run_exonerate ( query_file, out_file, target_file ):
   cmd = '/home/chilam01/Apps/exonerate/exonerate-2.2.0-x86_64/bin/exonerate --model affine:local --showcigar --bestn 1 \
            --exhaustive TRUE  --forwardcoordinates FALSE --alignmentwidth 200 ' + '--query ' + query_file + ' --target ' + target_file + ' 1>' + out_file + ' 2> /dev/null'
   os.system( cmd )
             
def revcomp(seq):
   return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]

def reversed_string(a_string):
   return a_string[::-1]

def get_soft_clipped_bases( cigar_str):
   if re.search( 'S$',cigar_str ):
      spl_cigar_str = re.split( '[A-Z]', cigar_str )
      soft_clipped_bases = int( spl_cigar_str[-2])
      return soft_clipped_bases
   else:
      soft_clipped_bases = 0
      return soft_clipped_bases
   

def process_exonerate( exo_out_file, read_length  ):
   # read exonerate out file
   exo_out = open( exo_out_file, 'r' )
   for out_line in exo_out:
      if 'Query:' in out_line:
         qname= out_line.split( sep='Query: ' )[1].rstrip().replace( '@', '' ).replace( ' [revcomp]', '')
         #print( '################ ', qname )

      if 'cigar:' in out_line:
         # Query alignment start and end
         q_aln_start = int(out_line.split( sep=' ')[2] ) + 1
         q_aln_end = int(out_line.split( sep=' ')[3] )

         # Target alignment start and end
         t_aln_start = int(out_line.split( sep=' ')[6] ) + 1
         t_aln_end = int(out_line.split( sep=' ')[7] )

         # cigar - string
         #cigar_str = out_line.split( sep='+')[-1].split( sep='  ').pop().replace( ' ', '').rstrip()
         cigar_str = out_line.split( sep='+')[-1].split( sep='  ').pop().rstrip()


         '''
         from cigar string need to reverse pairwise elemnts to fit SAM format
         old string = ['M', '109', 'I', '3', 'D', '10']
         new string looks like 109M3I10D 
         '''
         spl_cigar_str = cigar_str.split( sep=' ')

         half_list_size = len( spl_cigar_str ) / 2
         half_list_size = int(half_list_size)
         cigar_text=0
         cigar_num=1
         new_cigar_str=''
         for i in range( half_list_size):
            pair_str =  spl_cigar_str[ cigar_num] +  spl_cigar_str[ cigar_text]
            new_cigar_str = new_cigar_str + pair_str
            cigar_text +=2
            cigar_num += 2

         cigar_str = new_cigar_str
 
         # add soft clipping info if aln not end to end
         if q_aln_start > 1:
            #cigar_str = 'S' + str(q_aln_start - 1)  + cigar_str
            cigar_str =  str(q_aln_start - 1) + 'S' + cigar_str

         if q_aln_end < read_length:
            end_soft_clip = read_length - q_aln_end
            #cigar_str = cigar_str + 'S' + str( end_soft_clip )
            cigar_str = cigar_str + str( end_soft_clip ) + 'S'

         # alignment start position with respect to chromosome
         # amplicion or target is small region with in chromosome
         # sam requires chr start of the alignment
         aln_pos = amp_chr_start + (t_aln_start -1 )
   return [ qname, aln_pos, cigar_str ]
      


def alignment_test( exo_out_file ):
   exo_out = open( exo_out_file, 'r' )
   aligned = False
   for out_line in exo_out:
      if 'C4 Alignment:' in out_line:
         aligned = True 
   return aligned



########################## END OF FUNCTIONS   ################################


#############################################################################
# amplicon start and end
#path='/mnt/scratcha/bioinformatics/chilam01/projects/CRISPR/NGS/20150522_ChanD_DF_CRISPR/pairwiseToSam/pten/SLX-13774/FLD0224/test_on_real_data/GEP00002/FLD0277'
#sample_name='FLD0277'
#amp_chr_start = 89653730
#amp_chr_end = 89653854
#strand = '+'
#chromosome='chr10'
#r1_sam_flag = 99   # read paired, read mapped in proper pair, mate reverse strand,  first in pair
#r2_sam_flag = 147  # read paired, read mapped in proper pair, read reverse strand,  second in pair
#mqal_score = 60
#rnext = '='

parser = argparse.ArgumentParser()
parser.add_argument('-path', action='store', dest='path', type=str, help='data path')
parser.add_argument( '-sample_name', action='store', dest='sample_name', type=str, help='sample name' )
parser.add_argument( '-amp_start', action='store', dest='amp_chr_start', type=int, help='amplicon start position' )
parser.add_argument( '-amp_end', action='store', dest='amp_chr_end', type=int, help='amplicon end position')
parser.add_argument( '-strand', action='store', dest='strand', type=str, help='amplicon strand')
parser.add_argument( '-chr', action='store', dest='chromosome', type=str, help='amplicon chromosome name')
parser.add_argument( '-r1_flag', action='store', dest='r1_sam_flag', type=int, help='read 1 sam flag')
parser.add_argument( '-r2_flag', action='store', dest='r2_sam_flag', type=int, help= 'read 2 sam flag')
parser.add_argument( '-mapq', action='store', dest='mqal_score', type=int, help='sam maping quality score')
parser.add_argument( '-rnext', action='store', dest='rnext', type=str, help='sam column 7 value')
parser.add_argument( '-target', action='store', dest='target_file', type=str, help='target sequence file name')
parser.add_argument( '-r1_fq', action='store', dest='r1_fq', type=str, help='read 1 fastq file name')
parser.add_argument( '-r2_fq', action='store', dest='r2_fq', type=str, help='read 2 fastq file name')
parser.add_argument( '-sam_out', action='store', dest='sam_out', type=str, help='sam output file name')
parser.add_argument( '-sam_header', action='store', dest='sam_header_file', type=str, help='sam header file name')

args = parser.parse_args()

path = args.path
sample_name = args.sample_name
amp_chr_start = args.amp_chr_start
amp_chr_end = args.amp_chr_end
strand = args.strand
chromosome = args.chromosome
r1_sam_flag = args.r1_sam_flag
r2_sam_flag = args.r2_sam_flag
mqal_score = args.mqal_score
rnext = args.rnext
target_file = args.target_file
r1_fq =  args.r1_fq
r2_fq = args.r2_fq
sam_out_file = args.sam_out
sam_header_file = args.sam_header_file

#target_file = path + '/' + 'pten.txt'
#sam_out_file= path + '/' + sample_name + '.exonerate.sam'
sam_out = open(sam_out_file, 'w' )

############################################################################
#print_sam_header( 'bam_header_v2.0.txt', sample_name, out_file=sam_out )
print_sam_header( sam_header_file, sample_name, out_file=sam_out )

# read in fastq files
#fq1=open(  'SLX-13774.FLD0277.000000000-B85M8.s_1.r_1.fq', 'r')
fq1=open(  r1_fq, 'r')
fq2=open(  r2_fq, 'r')

# lines list
fq1_lines = fq1.read().splitlines()
fq2_lines = fq2.read().splitlines()

noOfLines = len(fq1_lines) 

for i in range( noOfLines ):
   if (i % 4) == 0:
      print( i)
      ### r1 process
      # temp R1 fasta 
      f_temp_r1_fa = path + '/' + sample_name  + '_temp_r1_fa'
      r1_fq_info = make_fasta_and_get_info( file_name = f_temp_r1_fa, lines = fq1_lines[i : i+4] ) 
      seq_read_length = r1_fq_info[0]
      r1_seq = r1_fq_info[1]
      r1_b_qal = r1_fq_info[2]


      # run R1 exonerate 
      r1_exo_out = path + '/' + sample_name + '_temp_r1_out'
      run_exonerate( query_file = f_temp_r1_fa, out_file = r1_exo_out, target_file = target_file ) 


      # read one alignment test
      r1_aln_test = alignment_test( r1_exo_out ) 



      ### r2 process
      # temp R2 fasta
      f_temp_r2_fa = path + '/' + sample_name + '_temp_r2_fa'
      r2_fq_info = make_fasta_and_get_info( file_name = f_temp_r2_fa, lines = fq2_lines[i : i+4] )

      # if reads from the gene on +ve strand
      # R1 aligns to forward starnd
      # R2 reverse complement align to forward starnd
      # ref seq is always forward strand
      r2_seq = r2_fq_info[1]
      r2_seq_rev_comp = revcomp( r2_seq )

      # if read is reverse complementaed then need to reverse the quality score string
      r2_b_qal = r2_fq_info[2]
      r2_b_qal_rev = reversed_string(r2_b_qal)

      # run R2 exonerate
      r2_exo_out = path + '/' + sample_name + '_temp_r2_out'
      run_exonerate( query_file = f_temp_r2_fa, out_file = r2_exo_out, target_file = target_file )

      # read2 alignment test
      r2_aln_test = alignment_test( r2_exo_out ) 


      #print( '######r1_aln_test == ', r1_aln_test, '\t', 'r2_aln_test == ', r2_aln_test, '\n' )


      # both r1 and r2 shold align to be in sam 
      # if any one or none aligned read pair is ignored

      if r1_aln_test == True and r2_aln_test == True:
         # process R1 exonorate output
         r1_exo_pro_out  = process_exonerate( exo_out_file=r1_exo_out, read_length = seq_read_length )
         qname = r1_exo_pro_out[0]
         r1_aln_start = r1_exo_pro_out[1]
         r1_cigar = r1_exo_pro_out[2]
   
         
   
         # process R2 exonorate output
         r2_exo_pro_out  = process_exonerate( exo_out_file=r2_exo_out, read_length = seq_read_length )
         r2_aln_start = r2_exo_pro_out[1]
         r2_cigar = r2_exo_pro_out[2]
   
   
         # get TLEN - template length
   
         '''           R1                               R2
              Start ---------> End                 S<-----------End
         ref<--------------------------------------------------------->
   
         TLEN = number of bases between leftmost mapped base to rightmost mapped base
         TLEN    = End of R2 - Start of R1
   
         Be aware of soft clippings at the ends
   
         '''  
         no_soft_clip_bases = get_soft_clipped_bases( r2_cigar )
         l_start = r1_aln_start
         r_end = (r2_aln_start + seq_read_length ) - no_soft_clip_bases
         t_len = r_end - l_start
   
         r1_tlen = str( t_len )
         r2_tlen = '-' +  str( t_len )       
   
   
         # r1 sam
         r1_sam_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(qname, r1_sam_flag, chromosome, r1_aln_start, mqal_score, r1_cigar, rnext, r2_aln_start, r1_tlen, r1_seq, r1_b_qal)
              
         # r2 sam
         r2_sam_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( qname, r2_sam_flag, chromosome, r2_aln_start, mqal_score, r2_cigar, rnext, r1_aln_start, r2_tlen, r2_seq_rev_comp, r2_b_qal_rev )
         sam_out.write( r1_sam_line )
         sam_out.write( r2_sam_line)
      else:
         next
