from argparse import ArgumentParser
import math

parser = ArgumentParser()


parser.add_argument('--data', default='' ,help='txt file with paths of NanoAOD data')
parser.add_argument('--mc_mu', default='',help='txt file with paths of NanoAOD mc mu')
parser.add_argument('--mc_tau', default='',help='txt file with paths of NanoAOD mc tau')
parser.add_argument('--mc_x', default='',help='txt file with paths of NanoAOD mc that needs to delete signal')
parser.add_argument('--mc_onia', default='',help='txt file with paths of NanoAOD mc that needs to delete signal')
parser.add_argument('--mc_gen', default='',help='txt file with paths of NanoAOD generic mc')

args = parser.parse_args()

