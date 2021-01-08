import os
import datetime
import sys

dataset_dict = {'BcToJpsiMuNu':['--mc_mu /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToJpsiMuNu_files_path.txt'],
                'BcToJpsiTauNu':['--mc_tau /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToJpsiTauNu_files_path.txt'],
                'OniaX':['--mc_onia /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018/OniaX_files_path.txt'],
                'data':['--data /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/data_files_path_7Dic.txt'],
                'BcToXToJpsi':['--mc_x /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/pathFiles/nanoAOD/Run2018UL/BcToXToJpsi_files_path.txt']
}

dataset = 'BcToJpsiTauNu'
dateFolder = '2020Dec09'
total_files = 2
new_files_per_job = 1
dataset_opt = dataset_dict[dataset][0]

files_ok = os.popen('ls /pnfs/psi.ch/cms/trivcat/store/user/friti/dataframes_' + dateFolder + '/'+ dataset).readlines()

numbers_ok = []

for fil in files_ok:
    print(fil)
    numbers_ok.append(int(fil.split("_")[2]))
print(numbers_ok)
numbers_ok.sort()
print(numbers_ok)

for i in range(total_files):
#for i in range(1):
    if(i not in numbers_ok):
        fin = open("/work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/dataframes_"+ dateFolder+ "/"+dataset+"/Resonant_Rjpsi_chunk"+str(i)+".py", "rt")
        lines = fin.readlines()
        files_per_job = 0
        skip_files = 0
        #        print(i)
        for line in lines[24:26]:
            #            print(line)
            if 'nMaxFiles' in line: files_per_job = int(line.split(" ")[2])
            if 'skipFiles' in line: skip_files = int(line.split(" ")[2])
            print("files per job",files_per_job,skip_files)
        #input file
        fin.close()
        #we have to send new jobs!
        fin2 = open("Resonant_Rjpsi_crab.py", "rt")
        lines = fin2.readlines()
        #output file to write the result to
        print(files_per_job,new_files_per_job)
        if(new_files_per_job <= files_per_job and files_per_job%new_files_per_job==0 ):
            nNewJobs= (files_per_job//new_files_per_job) 
            print("new jobs",nNewJobs)
        for j in range(nNewJobs):
            fout = open("/work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/dataframes_"+ dateFolder+ "/"+dataset+"/Resonant_Rjpsi_chunk%d_%d.py" %(i, j), "wt")
            #for each line in the input file
            print("Creating /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/dataframes_"+ dateFolder+ "/"+dataset+"/Resonant_Rjpsi_chunk%d_%d.py" %(i, j))
            for line in lines:
                #read replace the string and write to output file
                if 'HOOK_MAX_FILES' in line: fout.write(line.replace('HOOK_MAX_FILES' , '%s' %(new_files_per_job)))
                elif 'HOOK_FILE_OUT'   in line: fout.write(line.replace('HOOK_FILE_OUT'   , '/scratch/friti/%s/%s_UL_%d_%d' %(dataset, dataset,i,j)))
                elif 'HOOK_SKIP_FILES'in line: fout.write(line.replace('HOOK_SKIP_FILES', '%d' %(skip_files + new_files_per_job*j)))
                else: fout.write(line)
                #close input and  output files

            fout.close()
            fin2.close()
            
            flauncher = open("/work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/dataframes_"+ dateFolder+ "/"+dataset+"/submitter_chunk%d_%d.sh" %(i,j), "wt")
            if not ("BcToXToJpsi") in dataset:
                flauncher.write(
                    '''#!/bin/bash
                    cd {dir}
                    #scramv1 runtime -sh
                    mkdir -p /scratch/friti/{scratch_dir}
                    ls /scratch/friti/
                    python {cfg} {option}
                    ls /scratch/friti/{scratch_dir}
                    xrdcp /scratch/friti/{scratch_dir}/{dat}_UL_{i}_{j}_ptmax.root root://t3dcachedb.psi.ch:1094//{se_dir}/.
                    rm /scratch/friti/{scratch_dir}/{dat}_UL_{i}_{j}_ptmax.root'''.format(dir="/work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/dataframes_"+ dateFolder+ "/"+dataset, scratch_dir= dataset, cfg='Resonant_Rjpsi_chunk%d_%d.py' %(i,j), option= dataset_opt, dat = dataset,i=i,j=j, se_dir='/pnfs/psi.ch/cms/trivcat/store/user/friti/dataframes_' + dateFolder + '/'+ dataset)
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
        
            #command_sh_batch = 'sbatch -p wn --account=t3 -o %s/logs/chunk%d.log -e %s/errs/chunk%d.err --job-name=%s --time=60 --mem=6GB %s/submitter_chunk%d.sh' %(out_dir, ijob, out_dir, ijob, out_dir, out_dir, ijob)
            out_dir = "/work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/dataframes_"+ dateFolder+ "/"+dataset
            command_sh_batch = 'sbatch -p wn --account=t3 -o %s/logs/chunk%d_%d.log -e %s/errs/chunk%d_%d.err --job-name=%s  --mem=6GB %s/submitter_chunk%d_%d.sh' %(out_dir, i,j, out_dir, i,j, dataset, out_dir, i,j)
            print(command_sh_batch)
            os.system(command_sh_batch)
