
#Script that saves a .txt file with the paths of nanoAOD files created with CRAB. 
#In input datasets names separated by a comma, without blanck space

import subprocess
import os
import sys

nanotools = True

tier = 'CSCS'
crab_path = '/work/friti/CMSSW_10_6_5/src/PhysicsTools/NanoAODTools/crab/RJPsiNanoTools_2020Nov26'
#crab_path = '/work/friti/UL/CMSSW_10_6_14/src/PhysicsTools/BParkingNano/production/RJPsiNANO_2020Nov25'
#crab_path = '/work/friti/CMSSW_10_6_5/src/PhysicsTools/NanoAODTools/crab/'

if(nanotools):
    dataset_dict = {'BcToJpsiMuNu':['BcToJpsiMuNu_nanotools'],
                    'BcToJpsiTauNu':['BcToJpsiTauNu_nanotools'],
                    'OniaX':['OniaX_nanotools'],
                    'BcToXToJpsi':['BcToXToJpsi_2_nanotools'],
                }
else:
    dataset_dict = {
        #'data':['data_Run2018D_UL'],
        'BcToXToJpsi':['BcToXToJpsi2'],
        'data':['data_Run2018A_UL','data_Run2018B_UL','data_Run2018C_UL','data_Run2018D_UL'],
        'BcToJpsiMuNu':['BcToJpsiMuNu'],
        'BcToJpsiTauNu':['BcToJpsiTauNu'],
        'OniaX':['OniaX'],
        'BToJpsi_ToMuMu':['BToJpsi_ToMuMu']
                }
#in input i nomi dei dataset separati da una virgola, senza spazio
datasets = map(str,sys.argv[1].split(','))
 
for dataset in datasets:
    print(" ")
    print("Dataset: %s"%dataset)
    if not os.path.exists("txt_files/"):
        os.makedirs("txt_files/")
    f = open("txt_files/"+dataset+"_files_path.txt", "w")
    folders = dataset_dict[dataset]
    for folder in folders:
        url = os.popen('crab getoutput --xrootd --jobids=1 -d ' + crab_path + '/crab_' + folder + '/').readlines()[0]
        print(crab_path + '/crab_' + folder)
        print(url)
        s1 = url.split('friti')
        s2 = s1[1].split('0000')
        if(tier == 'PSI'):
            newurl = 'root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/friti' + s2[0]
        elif(tier == 'CSCS'):
            newurl = 'srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/friti' + s2[0]
            #newurl = 'root://cms-xrd-global.cern.ch//store/user/friti' + s2[0]


        i = 0
        print('Checking files in the folder '+newurl.strip('\n'))
        while True:
            files_name = os.popen('eval `scram unsetenv -sh`; gfal-ls '+ newurl.strip('\n')+'000'+str(i)).readlines()
            path_name = newurl.strip('\n')+'000'+str(i)
            if(len(files_name)==0):
                print("The folder does not exist: "+ str(path_name))
                break
            print('subfolder: '+'000'+str(i))
            for file in range(len(files_name)):
                if(files_name[file].strip('\n') == 'log'):
                        continue
                f.write(path_name+'/'+files_name[file]) 
            i+=1

    f.close()
    print("The file txt_files/"+dataset+"_files_path.txt has been created.")
    
print('\nGoodbye\n')
