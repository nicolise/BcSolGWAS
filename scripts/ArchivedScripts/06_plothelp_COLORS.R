06_plottingCOLORS

#VENN DIAGRAMS
#http://www.eulerdiagrams.org/eulerAPE/

#colorblind-friendly palettes can be found at
#http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
colourpicker::colourPicker()

#many-color set
myColors <- c("#999999", "#292929","#684800" ,"#CBA22A", "#63B2D3", "#1FA69D", "#57B761", "#DAD94C","#2B869D", "#EE82EE", "#D2652D", "#CC79A7")
#greyscale set
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80")
#AND fill = grey20 dark for domesticated, grey60 pale for wild
#domestication 2 set
#domestication is first color, pale green. Wild is lilac.
myColors <- c("#9EFA6C", "#AB82FF")
#colorblind-friendly 3 set
myColors <- c("#2F4F4F", "#9EFA6C", "#AB82FF")

#another 3 set (blue, orange, black)
#domesticated is blue, wild is orange, sensitivity is black
myColors <- c("#1C86EE", "#EE7600", "#050505")

names(myColors) <- levels(ModDat$Species)
colScale <- scale_colour_manual(name = "Species",values = myColors)


#paper dimensions: 7.5 by 5.4
#mini paper dimensions: 3.5 by 5.4
#presentation dimensions: 6 by 4

#after loading dataframe have to reorder domestication levels to plot
ModDat$SpLabs <- factor(ModDat$SpLabs, levels = c("Wild", "Domesticated"))
