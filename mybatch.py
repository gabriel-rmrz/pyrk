from argparse import ArgumentParser
import math

parser = ArgumentParser()

#averlo in input come file (prendi spunto da samples di crab)
#ho come input 3 file diversi, uno con i path dei dati, uno con i path mc tau e l'altro con path mc mu. Lui li prende e li divide in 3 gruppi. io faccio un for in quei tre gruppi e poi il for all'interno di ognuno.
#Un dizionario data: array con tutti i file data; mc mu: array file mu etc e nel codice principale giro prima su data mc e mc e per ognuno sui file

#fare sul codice principale 	f=open("guru99.txt", "r") e poi f1=f.readlines() e un ciclo for x in f1 e legge riga per riga. così posso leggere riga per riga i path dei file. quindi farei tutto nel programma principale, senza questo mybatch, che tanto non fa niente praticamente


parser.add_argument('f_data', help='txt file with paths of NanoAOD files of data')
parser.add_argument('f_mc_tau', help='txt file with paths of NanoAOD files of data of mc 1 (tau)')
parser.add_argument('f_mc_mu', help='txt file with paths of NanoAOD files of data of mc 2 (mu)')

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

