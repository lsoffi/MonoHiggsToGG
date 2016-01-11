#! /bin/sh

#############################################
#
# fitterFormatting input arguements:
#   1st: input directory
#   2nd: output directory
#   3rd: type (sig, bkg, data)
#   4th: input filename 
#   5th: sample name
#   6th: outfile name
#
# final files used for the fit are:
#   outdir/Output_MC.root
#   outdir/Output_Data.root
#
#############################################


indir="../macro/data/25ns_v7_LooseSel/"
outdir="ntuples4fit"

# Sidenote: Safely ignore warning message "tab completion not implemented for this context" 
# This comes from the tabs included for making the inputs to fitterFormatting easier to read below 

root -l -b << EOF
.L fitterFormatting.cc++

fitterFormatting("$indir","$outdir","sig",	"2HDM_mZP600.root",	"2HDM_mZP600_mA0300",	"2HDM_mZP600_new.root")
fitterFormatting("$indir","$outdir","sig",	"2HDM_mZP800.root",	"2HDM_mZP800_mA0300",	"2HDM_mZP800_new.root")
fitterFormatting("$indir","$outdir","sig",	"2HDM_mZP1000.root",	"2HDM_mZP1000_mA0300",	"2HDM_mZP1000_new.root")
fitterFormatting("$indir","$outdir","sig",	"2HDM_mZP1200.root",	"2HDM_mZP1200_mA0300",	"2HDM_mZP1200_new.root")
fitterFormatting("$indir","$outdir","sig",	"2HDM_mZP1400.root",	"2HDM_mZP1400_mA0300",	"2HDM_mZP1400_new.root")
fitterFormatting("$indir","$outdir","bkg",	"GluGluHToGG.root",	"GluGluHToGG",		"GluGluHToGG_new.root")
fitterFormatting("$indir","$outdir","bkg",	"VH.root",		"VH",			"VH_new.root")
fitterFormatting("$indir","$outdir","data",	"DoubleEG.root",	"DoubleEG",		"Output_Data.root")

.q

EOF
echo "Done"

echo "Adding MC Files Together"
hadd $outdir/Output_MC.root $outdir/2HDM_mZP600_new.root $outdir/2HDM_mZP800_new.root $outdir/2HDM_mZP1000_new.root $outdir/2HDM_mZP1200_new.root $outdir/2HDM_mZP1400_new.root $outdir/GluGluHToGG_new.root $outdir/VH_new.root 

#fitterFormatting("../macro/data/25ns_v7_LooseSel/QCD.root","bkg","QCD","QCD_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/GJets.root","bkg","GJets","GJets_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/DYJetsToLL.root","bkg","DYJetsToLL","DYJetsToLL_new.root")
#fitterFormatting("../macro/data/25ns_v7_LooseSel/DiPhoton.root","bkg","DiPhoton","DiPhoton_new.root")
