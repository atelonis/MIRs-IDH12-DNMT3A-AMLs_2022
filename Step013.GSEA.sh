#!/bin/bash 

[ ! -d output/GSEA ] && mkdir output/GSEA
[ ! -d output/GSEA/Quant3.SelectKEGG ] && mkdir output/GSEA/Quant3.SelectKEGG

Scripts/GSEA.sh selKEGG_Entrez output/GSEA/Quant3/Genelist.rnk output/GSEA/Quant3.SelectKEGG/output > output/GSEA/Quant3.SelectKEGG/Log.txt


