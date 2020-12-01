import os
import datetime
import sys

dataset_dict = {'BcToJpsiMuNu':['--mc_mu /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToJpsiMuNu_files_path.txt'],
                'BcToJpsiTauNu':['--mc_tau /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToJpsiTauNu_files_path.txt'],
                'OniaX':['--mc_onia /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018/OniaX_files_path.txt'],
                'data':['--data /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/data_files_path_17Nov.txt'],
                'BcToXToJpsi':['--mc_x /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToXToJpsi_files_path.txt']
}

dataset = 'BcToXToJpsi'
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

if ("BcToXToJpsi") in dataset:
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_hc_mu"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_hc_mu")

    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_chic0_mu"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_chic0_mu")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_chic1_mu"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_chic1_mu")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_chic2_mu"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_chic2_mu")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_jpsi_3pi"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_jpsi_3pi")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_jpsi_hc"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_jpsi_hc")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_jpsi_mu"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_jpsi_mu")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_jpsi_pi"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_jpsi_pi")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_jpsi_tau"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_jpsi_tau")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_psi2s_mu"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_psi2s_mu")
    if not os.path.exists('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir+ "/is_psi2s_tau"):
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/friti/' +out_dir + "/is_psi2s_tau")

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
    if not ("BcToXToJpsi") in dataset:
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
    else:
        flauncher.write(
            '''#!/bin/bash
            cd {dir}
            #scramv1 runtime -sh
            mkdir -p /scratch/friti/{scratch_dir}
            ls /scratch/friti/
            python {cfg} {option}
            ls /scratch/friti/{scratch_dir}
            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_chic0_mu.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_chic0_mu/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_chic0_mu.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_chic1_mu.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_chic1_mu/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_chic1_mu.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_chic2_mu.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_chic2_mu/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_chic2_mu.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_hc_mu.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_hc_mu/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_hc_mu.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_3pi.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_jpsi_3pi/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_3pi.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_hc.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_jpsi_hc/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_hc.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_mu.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_jpsi_mu/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_mu.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_pi.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_jpsi_pi/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_pi.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_tau.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_jpsi_tau/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_jpsi_tau.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_psi2s_mu.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_psi2s_mu/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_psi2s_mu.root

            xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_is_psi2s_tau.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/friti/{se_dir}/is_psi2s_tau/.
            rm /scratch/friti/{scratch_dir}/{dat}_UL_{ijob}_psi2s_tau.root


            
            '''.format(dir='/'.join([os.getcwd(), out_dir]), scratch_dir= dataset, cfg='Resonant_Rjpsi_chunk%d.py' %(ijob), option= dataset_opt, dat = dataset,ijob=ijob, se_dir=out_dir)
        )
        
        flauncher.close()
        
    command_sh_batch = 'sbatch -p quick --account=t3 -o %s/logs/chunk%d.log -e %s/errs/chunk%d.err --job-name=%s --time=60 --mem=6GB %s/submitter_chunk%d.sh' %(out_dir, ijob, out_dir, ijob, out_dir, out_dir, ijob)

    os.system(command_sh_batch)
    
