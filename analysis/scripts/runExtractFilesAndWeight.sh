#! /bin/sh
#run extractJSONS.py for all samples in file

# FLASHgg version 1_0_0 for 25ns
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/Data/DoubleEG.json		-o DoubleEG		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py	-i lists_25ns_v1_1_0/MC/ZpZH.json		-o ZpZH			-d lists_25ns_v1_1_0
python extractFilesAndWeight.py	-i lists_25ns_v1_1_0/MC/ttHJetToGG_M125.json		-o ttHJetToGG		-d lists_25ns_v1_1_0
python extractFilesAndWeight.py	-i lists_25ns_v1_1_0/MC/VBFHToGG_M-125.json		-o VBFHToGG		-d lists_25ns_v1_1_0

#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/GJet_Pt-20to40.json	-o GJet_Pt-20to40	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/VHToGG_M125.json	-o VH			-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/QCD_Pt-30to40.json	-o QCD_Pt-30to40	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/QCD_Pt-40toInf.json	-o QCD_Pt-40toInf	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/QCD_Pt-30toInf.json	-o QCD_Pt-30toInf	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/DiPhoton.json		-o DiPhoton		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/GluGluHToGG_M-125.json	-o GluGluHToGG		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/DYJetsToLL.json		-o DYJetsToLL		-d lists_25ns_v1_1_0
#
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-600.json	-o 2HDM_MZP600		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-1000.json	-o 2HDM_MZP1000		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-1200.json	-o 2HDM_MZP1200		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-1400.json	-o 2HDM_MZP1400		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-1700.json	-o 2HDM_MZP1700		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-2000.json	-o 2HDM_MZP2000		-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-2500.json	-o 2HDM_MZP2500		-d lists_25ns_v1_1_0
##python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZp-800.json	-o 2HDM_MZP800		-d lists_25ns_v1_1_0
#
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP600.json	-o privMC_2HDM_MZP600	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP800.json	-o privMC_2HDM_MZP800	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP1000.json	-o privMC_2HDM_MZP1000	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP1200.json	-o privMC_2HDM_MZP1200	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP1400.json	-o privMC_2HDM_MZP1400	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP1700.json	-o privMC_2HDM_MZP1700	-d lists_25ns_v1_1_0
#python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP2500.json	-o privMC_2HDM_MZP2500	-d lists_25ns_v1_1_0
##python extractFilesAndWeight.py -i lists_25ns_v1_1_0/MC/2HDM_MZP2000.json	-o privMC_2HDM_MZP2000	-d lists_25ns_v1_1_0


# all 25ns, Spring15BetaV7 
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP600.json		-o 2HDM_MZP600		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP800.json		-o 2HDM_MZP800		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP1000.json		-o 2HDM_MZP1000		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP1200.json		-o 2HDM_MZP1200		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP1400.json		-o 2HDM_MZP1400		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP1700.json		-o 2HDM_MZP1700		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP2000.json		-o 2HDM_MZP2000		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/2HDM_MZP2500.json		-o 2HDM_MZP2500		-d lists_25ns_v7

#python extractFilesAndWeight.py -i lists_25ns_v7/MC/GJet_Pt-20to40.json		-o GJet_Pt-20to40	-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf	-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/VHToGG_M125.json		-o VH			-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/QCD_Pt-30to40.json		-o QCD_Pt-30to40	-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/QCD_Pt-40toInf.json		-o QCD_Pt-40toInf	-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/QCD_Pt-30toInf.json		-o QCD_Pt-30toInf	-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/DiPhoton.json		-o DiPhoton		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/GluGluHToGG_M-125.json	-o GluGluHToGG		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/MC/DYJetsToLL.json		-o DYJetsToLL		-d lists_25ns_v7
#
#python extractFilesAndWeight.py -i lists_25ns_v7/Data/DoubleEG.json		-o DoubleEG		-d lists_25ns_v7
#python extractFilesAndWeight.py -i lists_25ns_v7/Data/DoubleEG_RunD.json	-o DoubleEG_RunD	-d lists_25ns_v7





#all FLASHgg Spring15BetaV4 
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/GJet_Pt-20to40.json	-o GJet_Pt-20to40	-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf	-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/VHToGG_M125.json	-o VH			-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/QCD_Pt-30to40.json	-o QCD_Pt-30to40	-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/QCD_Pt-40toInf.json	-o QCD_Pt-40toInf	-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/QCD_Pt-30toInf.json	-o QCD_Pt-30toInf	-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/DiPhoton.json		-o DiPhoton		-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/GluGluHToGG_M-125.json	-o GluGluHToGG		-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/MC/DYJetsToLL.json		-o DYJetsToLL		-d lists_50ns_betaV4
#python extractFilesAndWeight.py -i lists_50ns_betaV4/Data/DoubleEG.json		-o DoubleEG		-d lists_50ns_betaV4

#mix of Flashgg versions (spring15betaV2 and betaV4)
#python extractFilesAndWeight.py -i lists_50ns/MC/DiPhoton.json  -o DiPhoton     -d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ZH.json	-o ZH 		-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/WplusH.json	-o WplusH 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/WminusH.json	-o WminusH 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/Higgs_scalar_nohdecay_gg_1000GeV.json	-o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/Higgs_scalar_nohdecay_gg_100GeV.json	-o Higgs_scalar_nohdecay_gg_100GeV	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/Higgs_scalar_nohdecay_gg_10GeV.json	-o Higgs_scalar_nohdecay_gg_10GeV	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/Higgs_scalar_nohdecay_gg_1GeV.json	-o Higgs_scalar_nohdecay_gg_1GeV	-d lists_50ns

#signal MC
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_1000GeV.json	-o Higgs_scalar_nohdecay_gg_1000GeV	-d lists_sig1
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_100GeV.json	-o Higgs_scalar_nohdecay_gg_100GeV	-d lists_sig1
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_10GeV.json	-o Higgs_scalar_nohdecay_gg_10GeV	-d lists_sig1
#python extractFilesAndWeight.py -i lists_sig1/Higgs_scalar_nohdecay_gg_1GeV.json	-o Higgs_scalar_nohdecay_gg_1GeV	-d lists_sig1

#MC
#python extractFilesAndWeight.py -i lists_50ns/MC/GJet_Pt-20to40.json	-o GJet_Pt-20to40	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/GJet_Pt-40toInf.json	-o GJet_Pt-40toInf 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/QCD_Pt-30to40.json	-o QCD_Pt-30to40  	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/QCD_Pt-30toInf.json	-o QCD_Pt-30toInf 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/QCD_Pt-40toInf.json	-o QCD_Pt-40toInf 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/GluGluHToGG.json	-o GluGluHToGG  	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ttHJetToGG_M120.json	-o ttHJetToGG_M120 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ttHJetToGG_M125.json	-o ttHJetToGG_M125 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/MC/ttHJetToGG_M130.json	-o ttHJetToGG_M130 	-d lists_50ns
##Data
#python extractFilesAndWeight.py -i lists_50ns/Data/SinglePhoton.json	-o SinglePhoton 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/Data/DoubleEG.json	-o DoubleEG	 	-d lists_50ns
#python extractFilesAndWeight.py -i lists_50ns/Data/EGamma.json 		-o EGamma	 	-d lists_50ns
