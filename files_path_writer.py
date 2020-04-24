import subprocess
import os

crab_path = '/work/friti/new/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/production/BParkingNANO_2020Apr21'
datasets = ['data','BcToJpsiMuNu','OniaX']



for dataset in datasets:
    print(" ")
    print("Dataset: %s"%dataset)
    if not os.path.exists("txt_files/"):
        os.makedirs("txt_files/")
    f=open("txt_files/"+dataset+"_files_path.txt", "w")
    
    if(dataset=='data'):
        folders=[]
        #folders=['data_Run2018A_sep_v1', 'data_Run2018B_sep_v1', 'data_Run2018C_sep_v1', 'data_Run2018D_pr_v2']                                                                                         
    if(dataset=='BcToJpsiMuNu'):
        folders=['BcToJpsiMuNu_JpsiToMuMu']
    if(dataset=='OniaX'):
        folders=['OniaX_JpsiToMuMu']
    
    for folder in folders:
        url = os.popen('crab getoutput --xrootd --jobids=1 -d ' + crab_path + '/crab_' + folder + '/').readlines()[0]
        s1=url.split('friti')
        s2=s1[1].split('0000')
        newurl='root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/friti' + s2[0]

        i=0
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
    
    #print('Setting cmsenv again')
    #os.popen("eval `scram runtime -csh`")
print('\nGoodbye\n')
