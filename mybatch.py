from argparse import ArgumentParser
import math

parser = ArgumentParser()


parser.add_argument('--f_1', default='' ,help='txt file with paths of NanoAOD data')
parser.add_argument('--f_2', default='',help='txt file with paths of NanoAOD mc mu')
parser.add_argument('--f_3', default='',help='txt file with paths of NanoAOD mc tau')
parser.add_argument('--f_4', default='',help='txt file with paths of NanoAOD mc onia')

#,nargs='+',
#nargs='?',
#parser.add_argument('f_out', help='file path')
#parser.add_argument('nchunks', type=int)
#parser.add_argument('ichunk', type=int)
#parser.add_argument('--mc', action = 'store_true')
args = parser.parse_args()



#all_files = [i.strip() for i in open(args.f_in) if i.strip() and not i.startswith('#')]
#chunk_size = math.ceil(len(all_files) / args.nchunks)
#infiles = all_files[chunk_size*args.ichunk:chunk_size*(args.ichunk+1)]
#infiles=[args.f_in]
#if infiles:
#    print(len(infiles), "to be processed.")
#else:
#    print('Nothing to be done, no files left')
#    exit(0)

