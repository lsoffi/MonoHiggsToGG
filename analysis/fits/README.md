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
 fitterFormatting input arguements:
   1st: input directory
   2nd: output directory
   3rd: type (sig, bkg, data)
   4th: input filename 
   5th: sample name
   6th: outfile name

 final files used for the fit are:
   `outdir/Output_MC.root` and 
   `outdir/Output_Data.root`

Nb: MET categorization no longer implemented in fitterFormatting, but if it is implemented in the future... inside fitterFormatting you need to specify the MET bin criteria.

## Step 2)
Check that the sample shapes don't change in the different MET bins (this only works when MET cat is implemented -- i.e. not working at the moment):

```root -l ProduceWorkspaces.C```

which calls runfits() from `FitTools.cc`. In FitTools, need to specify the number of MET cat (nMetCat) and the number of PHO cat (nPhoCat).
And in the runfits function specify files to run over with format:
`AddSigData(workspace,sample,type)`

where if (type != 0) name = sample.
Or if (type == 1) then script assumes using 2HDM sample, so just put in sample=mZP mass.

## Step 3) 
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
   


