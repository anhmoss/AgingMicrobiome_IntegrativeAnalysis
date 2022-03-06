#function to log10 normalize according to the Fodor lab's method
#there is a step to also remove samlples with no counts at all (rowSum ==0)
#input: table with counts to be normalized, usually a df
#output: table of lognormalized data

lognorm_function = function(table) {

tab1 = table
tab1 =tab1[which(rowSums(tab1)!=0),]
tab2=tab1/rowSums(tab1)*mean(rowSums(tab1))
lognorm_table =log10(tab2+1)

return(lognorm_table)

}

