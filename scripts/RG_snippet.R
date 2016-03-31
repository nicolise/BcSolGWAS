#Nicole E Soltis
#031116
#single-isolate linear models

#-------------------------------------------------------------

#split dataset by isolate
out <- split( MYDATA , f = MYDATA$ISOLATE )
head(out[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file='FILENAME.txt')
#skip 13: blank, 59: 94.1, 68: blank, 77: Gallo3, 99: blank
#if you get error messages in the next step, you must figure out which particular isolates the script does not work for. Then skip those in the c(1:96) list.
for (i in c(1:96)) {
  print(unique(out[[i]]$ISOLATE))
  Mod <- lmer(LESIONSIZE ~ SPECIES/PLANTGENO + (1|EXPERIMENT), data=out[[i]])
  result <- anova(Mod)
  random <- rand(Mod)
  print(result)
  print(random)
}
sink()

#YOU CAN ALSO TRY ADDING A COUPLE MORE TERMS TO THE MODEL SUCH AS BLOCK IF YOU WANT

#-----------------------------------------------------------
#FDR cutoff

#NEXT, make a table with the p-values from each of your outputs from part 1. You can then read that table in to R and change the p-values to correct for the false discovery rate (FDR)

#p is a vector of p values

MyPvals <- read.csv("PVALUEDATA.csv")
names(MyPvals)
p <- MyPvals$SPECIESPVALUES
MyPvals$SPECIESFDR <- p.adjust(p, method = "fdr", n = length(p))

write.csv(MyPvals, "MYFILE.csv")
