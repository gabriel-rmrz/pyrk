def create_datacard(root_file_name,data_int,mu_int,tau_int,x_int):


    f= open("datacard.txt","w")
    f.write("imax 1 number of channels \n")
    f.write("jmax 2 number of backgrounds \n")
    f.write("kmax * number of independent systematical uncertainties \n") 
    f.write("--------------------------------------------------------------------------- \n") 
    f.write("shapes * * /work/friti/new/CMSSW_10_2_15/src/pyrk/root_files/%s $PROCESS \n"%root_file_name) 
    f.write("--------------------------------------------------------------------------- \n") 
    f.write("bin rjpsi\n") 
    f.write("observation %s\n"%data_int) 
    f.write("----------------------------------------------------------------------------\n") 
    f.write("bin \t rjpsi \t rjpsi \t rjpsi\n") 
    f.write("process \t mc_mu \t mc_tau \t mis_id\n") 
    f.write("process \t 1 \t 0 \t 2\n") 
    f.write("rate \t %s \t %s \t %s \n"%(mu_int,tau_int,x_int)) 
    f.write("----------------------------------------------------------------------------\n") 
    f.write("mu_norm rateParam rjpsi mc_mu 1\n") 
    f.write("mis_id_unc lnN - - 1.5\n") 
    f.write("#rjpsi autoMCStats 1\n") 
    f.close()

#create_datacard("Bm_miss_sq2.root",16626,14863.04568741178,1040.9416810714274,721.9011999999963)
