# Folder with script and data for paper "A partial prehistory of the Southwest Silk Road: Archaeometallurgical data from the sub-Himalayan corridor."

List of files and folders:

- `data/`
  + `data/DATA isotopie_Pryce_V13092022.csv` Frozen verison used in `scriptMeetingWednesday16.Rmd`
  + `data/DATA isotopie_Pryce_V4.2.1.csv` Last version of the dataset
  + `data/*.RDS` RDS of the object as used for final publication
  + `data/RG_table_s1.csv` table with elemental compition from [Radivojević & Grujić (2017)](https://academic.oup.com/comnet/article/6/1/106/4030792)
  + `data/RG_table_s2.csv` table with communities found in [Radivojević & Grujić (2017)](https://academic.oup.com/comnet/article/6/1/106/4030792)
- `SSRgeoloc/` `kmz` file with single coordinates collected by OP
- `maps/` multiple useful shapefiles used for visualisation and analysis
  + `maps/allsea.shp` big shapefile.
  + `maps/allseaborder_lowq.shp` lower quality general shapefile.
  + `maps/riverbasin_sitelimited.shp` river basin where sites falls in.
- `vignettes/` R vignettes that describe series of analysis and tool
  + `vignettes/reproduceRG2017.Rmd` a vignette that reproduce some of the analysis done in [Radivojević & Grujić (2017)](https://academic.oup.com/comnet/article/6/1/106/4030792)
  + `vignettes/exploreSEALIPelementalNetwork.Rmd` a vignette reproduce the analysis of [Radivojević & Grujić (2017)](https://academic.oup.com/comnet/article/6/1/106/4030792) using SEALIP data
  + `vignettes/recreateNetworkSSR.Rmd` this is the main vignette used throughout the paper 
- `scripts/` small script that concatenate various R functions to quickly generate useful tools/outputs
  + `scripts/shapefileope.Rmd` a script that was used to create the shapefiles avalable in `maps/`
  + `scripts/paperFigure.R` a script  to generate high quality figures used for the paper
  + `scripts/cleanDataset.R` a script to clean the csv and generate labels automatically 
- `R/` set of R-methods and functions used to run the analysis


To render a vignette with R in your command line ; from the root folder of this reporitory

```
Rscript  -e  "rmarkdown::render('doc/reproduceRG2017.Rmd',knit_root_dir='..')"
```

if you wanter to render them all:

```
for v in vignettes/*.Rmd; 
do
    Rscript  -e  "rmarkdown::render('$v',knit_root_dir='..')"
done
```
