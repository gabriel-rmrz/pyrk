import sys
import os
import datetime

var = 'mu'
xmin = 0
xmax = 1
flag = 'UL_newctau'
cut = '04'
FRValue = 0.85373
production_tag = datetime.date.today().strftime('%Y%b%d')

path_dir = var+ "_" + production_tag + "_cut"+cut
if not os.path.exists(path_dir):
# if there si already one it does not delete it
    os.makedirs(path_dir)
    print("Made directory "+ path_dir)
else:
    print("Directory already exists! ")
    sys.exit("WARNING: the folder "+ path_dir + " already exists!")

passrootfileName =  var + "Pass_cut" + cut + ".root"
os.system("cp /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/rootFiles/" + passrootfileName +" "+ path_dir + "/")
print("Copied root File " + passrootfileName)

failrootfileName = var + "Fail_cut" + cut + ".root"
os.system("cp /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/rootFiles/" + failrootfileName +" " + path_dir+ "/")
print("Copied root File " + failrootfileName)

datacardfailName = "datacard_" + var + "_Fail.txt"
os.system("cp /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/datacards/datacard_" + var + "_Fail.txt " + path_dir+ "/")
print("Copied datacard   datacard_" + var + "_Fail.txt" )

datacardpassName = "datacard_" + var + "_Pass.txt"
os.system("cp /work/friti/new/CMSSW_10_2_15/src/pyrk/scripts/datacards/datacard_" + var + "_Pass.txt " + path_dir+ "/")
print("Copied datacard   datacard_" + var + "_Fail.txt" )

#workspace
fin = open("workspace.C", "rt")
fout = open("%s/workspace.C"%(path_dir),"wt")
for line in fin:
    if 'REPLACE' in line:
        newline = line.replace('REPLACE_FILE_FAIL', '%s'%(var + "Fail_cut" + cut + ".root"))
        newline = newline.replace('REPLACE_FILE_PASS', '%s'%(var + "Pass_cut" + cut + ".root"))
        newline = newline.replace('REPLACE_FR', '%s'%(FRValue))
        newline = newline.replace('REPLACE_VAR', '%s'%(var))
        newline = newline.replace('REPLACE_XMIN', '%s'%(xmin))
        newline = newline.replace('REPLACE_XMAX', '%s'%(xmax))
        fout.write(newline)
       
    else:
        fout.write(line)
fout.close()
fin.close()
print("Created workspace workspace.C")

#entering the right directory
os.chdir(path_dir)

#combine datacards and fitDiagnostics
os.system("root -l -q workspace.C")
os.system("combineCards.py " + datacardpassName + " " + datacardfailName + ">& datacard.txt")
os.system("combine -M FitDiagnostics --plots --robustFit 1 --saveShapes --cminDefaultMinimizerStrategy 0 --saveNormalizations --maxFailedSteps 20 --saveWithUncertainties --robustHesse 1  --ignoreCovWarning datacard.txt")

#yields
os.system("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py -u fitDiagnostics.root")
output = os.popen("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py -u fitDiagnostics.root").readlines()


for line in output:
    if 'ch1' in line:
        if 'mc_comb' in line:
            splittedLine = line.split(" ")
            smartLine = []
            for piece in splittedLine:
                if piece != '':
                    smartLine.append(piece)
            preComb = smartLine[2]
            preCombUnc = smartLine[4]
            postComb  = smartLine[5]
            postCombUnc = smartLine[7]

        if 'mc_mu' in line:
            splittedLine = line.split(" ")
            smartLine = []
            for piece in splittedLine:
                if piece != '':
                    smartLine.append(piece)
            preMu = smartLine[2]
            preMuUnc = smartLine[4]
            postMu  = smartLine[5]
            postMuUnc = smartLine[7]

        if 'mc_tau' in line:
            splittedLine = line.split(" ")
            smartLine = []
            for piece in splittedLine:
                if piece != '':
                    smartLine.append(piece)
            preTau = smartLine[2]
            preTauUnc = smartLine[4]
            postTau  = smartLine[5]
            postTauUnc = smartLine[7]

ratioComb = float(postComb)/ float(preComb)
ratioMu = float(postMu) / float(preMu)
ratioTau = float(postTau) / float(preTau)

print("The ratios are: ")
print("- Factor for mu sample = %s"%(ratioMu))
print("- Factor for tau sample = %s"%(ratioTau))
print("- Factor for comb sample = %s"%(ratioComb))

