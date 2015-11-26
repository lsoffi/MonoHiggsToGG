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

fitterFormatting("test.root","sig","2HDM_mZP600_mA0300","outputtest.root")

.q

EOF
echo "Done"
