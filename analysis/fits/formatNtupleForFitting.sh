#! /bin/sh

#############################################
#
# fitterFormatting input arguements:
#   1st: file to convert
#   2nd: type (sig, bkg, data)
#   3rd: sample name
#   4th: output file name
#
#############################################




root -l -b << EOF
.L fitterFormatting.cc++

fitterFormatting("../macro/data/25ns_v7_LooseSel/2HDM_mZP600.root","sig","2HDM_mZP600_mA0300","2HDM_mZP600_new.root")
fitterFormatting("../macro/data/25ns_v7_LooseSel/2HDM_mZP800.root","sig","2HDM_mZP800_mA0300","2HDM_mZP800_new.root")
fitterFormatting("../macro/data/25ns_v7_LooseSel/2HDM_mZP1000.root","sig","2HDM_mZP1000_mA0300","2HDM_mZP1000_new.root")
fitterFormatting("../macro/data/25ns_v7_LooseSel/2HDM_mZP1200.root","sig","2HDM_mZP1200_mA0300","2HDM_mZP1200_new.root")
fitterFormatting("../macro/data/25ns_v7_LooseSel/2HDM_mZP1400.root","sig","2HDM_mZP1400_mA0300","2HDM_mZP1400_new.root")

.q

EOF
echo "Done"


#fitterFormatting("test.root","sig","2HDM_mZP600_mA0300","outputtest.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/VH.root","bkg","VH","VH_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/QCD.root","bkg","QCD","QCD_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/GJets.root","bkg","GJets","GJets_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/DYJetsToLL.root","bkg","DYJetsToLL","DYJetsToLL_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/DiPhoton.root","bkg","DiPhoton","DiPhoton_new.root")

