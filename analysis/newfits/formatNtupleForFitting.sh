#! /bin/sh

#############################################
#
# fitterFormatting input arguements:
#   1st: file to convert
#   2nd: type (sig, bkg, data)
#   3rd: sample name
#
#############################################



root -l -b << EOF
.L fitterFormatting.cc++

fitterFormatting("test.root","sig","2HDM_mZP600_mA0300")

.q

EOF
echo "Done"
