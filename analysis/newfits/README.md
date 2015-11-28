----------------------------------------------------------
# To make fits:
----------------------------------------------------------

## Step 1) 
Convert ntuples from FLASHgg format to format for fits by using fitterFormatting.cc 
This can be called by:  

```
  ./formatNtupleForFitting.sh
```
  
specifying fitterFormatting input arguements:
  1st: file to convert
  2nd: type (sig, bkg, dat)
  3rd: sample name
  4th: output file name

## Step 2) 
Run the fitter (combine_maker.py & templates_maker.py) called by 
`./combine_maker.sh <analysis_version> <list-of-options>`

- runs template_maker.py to generate the input workspace from the analysis trees
- runs background model, signal model and datacard creation according to the options
- all the outuput goes to a dedicate folder named `<analysis_version>_<fitname>_lumi_<luminosity>_<background_model>[_bias][_use_templates][_<extra_label>]`
    the extra label can be specified through the `--label` option


Example:
```
./combine_maker.sh  full_analysis_spring15_7412v2_sync_idv2_v2 --luminosity 0.82 --data-file output.root --default-model pow --label mytest  --fit-name cic2
```
   


