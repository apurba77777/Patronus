# Project for Analysis of Transients in Radio Observations with Not-so Unique Software

Patronus is a collection of scripts conjured for interferometric data reduction and image-based transient search in the LIGER project using some simple muggle-friendly spells. Patronus also supports some advanced and dangerous magic. 

Going through the following references is highly recommended before using this software: 
* [https://www.bloomsbury.com/uk/discover/harry-potter/fun-facts/harry-potter-glossary/](url)
* [https://harrypotter.fandom.com/wiki/List_of_spells](url).

Patronus is neither automated nor an end-to-end pipeline. It is _NOT_ designed for general interferometric data analysis. If you are not sure if it is suitable for your purpose, it is probably _NOT_ so.

Patronus uses CASA ([https://casa.nrao.edu/](url)) for calibration and imaging of GMRT interferometric data. aNKflag ([https://github.com/devojyoti96/ankflag](url)) is used for identification and flagging of RFI-affected visibilities from data.

Additionally, Patronus uses a number of Python packages listed in `templates/basilisk.txt`. 

Basic instructions for using Patronus are given in the github Wiki. 
