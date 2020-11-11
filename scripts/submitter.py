import os
import datetime
import sys

dataset_dict = {'BcToJpsiMuNu':['--mc_mu /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToJpsiMuNu_files_path.txt'],
                'BcToJpsiTauNu':['--mc_tau /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToJpsiTauNu_files_path.txt'],
                'OniaX':['--mc_x /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018/OniaX_files_path.txt'],
                'data':['--data /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/data_files_path_pichannel.txt']}

dataset = 'data'
dataset_opt = dataset_dict[dataset][0]

#splitting of files
file_name = dataset_opt.split(' ')[1]
count_files = len(open(file_name).readlines(  ))

files_per_job = 200
njobs = count_files//files_per_job + 1  

print("Submitting %s jobs" %(njobs))
print("Submitting in quick.")
production_tag = datetime.date.today().strftime('%Y%b%d')

date_dir = 'dataframes_'+ production_tag


out_dir = date_dir + '/%s' %(dataset)

##########################################################################################
##########################################################################################

# make output dir
if not os.path.exists(date_dir):
# if there si already one it does not delete it
    os.makedirs(date_dir)

#do not write in the same folder if I forget to change the name of the folder!
if os.path.exists(out_dir):
    sys.exit("WARNING: the folder "+ out_dir + " already exists!")
    
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    os.makedirs(out_dir + '/logs')
    os.makedirs(out_dir + '/errs')
os.system('cp nanoframe.py '+ out_dir+ '/.')
os.system('cp mybatch.py '+ out_dir+ '/.')

if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' + date_dir):
    os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' + date_dir)

if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir):
    os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir)
else:
    sys.exit("WARNING: the folder "+ out_dir + " already exists in the SE!")

for ijob in range(njobs):

    #input file
    fin = open("Resonant_Rjpsi_crab.py", "rt")
    #output file to write the result to
    fout = open("%s/Resonant_Rjpsi_chunk%d.py" %(out_dir, ijob), "wt")
    #for each line in the input file
    for line in fin:
        #read replace the string and write to output file
        if   'HOOK_MAX_FILES' in line: fout.write(line.replace('HOOK_MAX_FILES' , '%s' %files_per_job))
        elif 'HOOK_FILE_OUT'   in line: fout.write(line.replace('HOOK_FILE_OUT'   , '/scratch/friti/%s/%s_UL_%d' %(dataset, dataset,ijob)))
        elif 'HOOK_SKIP_FILES'in line: fout.write(line.replace('HOOK_SKIP_FILES', '%d' %(files_per_job*ijob)))
        else: fout.write(line)
    #close input and output files
    fout.close()
    fin.close()

    flauncher = open("%s/submitter_chunk%d.sh" %(out_dir, ijob), "wt")
    flauncher.write(
'''#!/bin/bash
cd {dir}
#scramv1 runtime -sh
mkdir -p /scratch/friti/{scratch_dir}
ls /scratch/friti/
python {cfg} {option}
ls /scratch/friti/{scratch_dir}
xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/.
rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}.root

'''.format(dir='/'.join([os.getcwd(), out_dir]), scratch_dir= dataset, cfg='Resonant_Rjpsi_chunk%d.py' %(ijob), option= dataset_opt, dat = dataset,ijob=ijob, se_dir=out_dir)
    )
    flauncher.close()
    
    command_sh_batch = 'sbatch -p quick --account=t3 -o %s/logs/chunk%d.log -e %s/errs/chunk%d.err --job-name=%s --time=60 --mem=6GB %s/submitter_chunk%d.sh' %(out_dir, ijob, out_dir, ijob, out_dir, out_dir, ijob)

    os.system(command_sh_batch)
    
    
    
