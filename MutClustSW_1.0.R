#-------------------------------------------
# MutClustSW V1.0
#-------------------------------------------

# the example mutation file includes all the COSMIC missense mutations on KRAS (except for those with 'TCGA-##-####')
# the header should be "symbol	position (and additional columns separated by tab)"
# the actual file may have numerous lines such as "KRAS(tab)12" indicating one incidence of KRAS harbor one somatic mutation on 12th codon 
mutation.missense = read.table("mutation_positions_KRAS.txt", header=TRUE, sep="\t")

# the example gene length file includes the gene symbol of KRAS and the length (the number of amino acids length)
# the head should be "symbol(tab)length (and additional columns separated by tab)"
gene.lengths = read.table("gene_length_KRAS.txt", header=TRUE, sep="\t")

output.file1 = c("mutClustSW_results_v1.0.txt")
write.table(cbind("Gene", "Start", "End", "Scoresum", "Recursion", "Mutation", "NES", "Pvalue"), file = output.file1, row.names=F, col.names=F, quote=F, sep="\t")

# repeat the procedure for all the genes with missense mutations
genes.to.test = unique(mutation.missense$symbol)

for(k in 1:length(genes.to.test)){	

	print(paste("Processing ", k, "/", length(genes.to.test)))

	target.symbol = genes.to.test[[k]] 

	### obtain mutation.positions (missense) and mutation.positions.silent (silent) - ignore the genes without enough mutations 
	mutation.gene = mutation.missense[which(mutation.missense$symbol== as.character(target.symbol)),]
	gene.length = gene.lengths$length[which(gene.lengths$symbol ==  as.character(target.symbol))]
	mutation.positions = mutation.gene$position
	if (length(mutation.positions) < 10) next # (recommended) minmutation = 10 for missense and 5 for nonsense/silent

	### swarray.random to estimate the mean and standard deviation of scores in permuted data
	swscore.random = c()
	for(perm in 1:100){
		mutation.positions.perm = sample(1:gene.length, length(mutation.positions), replace=T)
		mut.pos.table = table(mutation.positions.perm)

		# prepare one-dimensional score - delta value of '1/gene.length' is introduced to ensure the negative sum of the score
		swarray = rep(-1/(gene.length - length(mut.pos.table))-1/gene.length, gene.length)
		swarray[as.numeric(names(mut.pos.table))] = 1/length(mutation.positions.perm) * mut.pos.table - 1/gene.length

		# sw-array for permuted mutations
		sws = 0; swb = 0; swscore = 0; swstart = 0; swend = 0;
		for (i in 1:gene.length) {
			if (sws + swarray[i] > 0) {sws = sws + swarray[i]} else {sws = 0; swb = i + 1;}
			if (sws > swscore) {swscore = sws; swstart = swb; swend = i;}
		}
		swscore.random = c(swscore.random, swscore)
	}

	### sw-score/segment discovery - lasts until no segment with positive sw-score
	mut.pos.table = table(mutation.positions)
	swarray = rep(-1/(gene.length - length(mut.pos.table)) - 1/gene.length, gene.length)
	swarray[as.numeric(names(mut.pos.table))] = 1/length(mutation.positions) * mut.pos.table - 1/gene.length

	while (1) {

		sws = 0; swb = 1; swscore = 0; swstart = 0; swend = 0;
		for (i in 1:gene.length) {
			if (sws + swarray[i] > 0) {sws = sws + swarray[i]} else {sws = 0; swb = i + 1;}
			if (sws > swscore) {swscore = sws; swstart = swb; swend = i;}
		}

		if (swscore == 0) break # exit when no more positive-score segments

		# Recursive (single-round) - takes the recursive segments when swscore.recur > swscore (size-adjusted)
		gene.length.recur = swend - swstart + 1
		mutation.positions.recur = subset(mutation.positions, mutation.positions >= swstart & mutation.positions <= swend) - swstart + 1
		mut.pos.table.recur = table(mutation.positions.recur)

		# prepare one-dimensional score
		swarray.recur = rep(-1/(gene.length.recur - length(mut.pos.table.recur)), gene.length.recur)
		swarray.recur[as.numeric(names(mut.pos.table.recur))] = 1/length(mutation.positions.recur) * mut.pos.table.recur

		sws.recur = 0; swb.recur = 1; swscore.recur = 0; swstart.recur = 0; swend.recur = 0;
		for (i in 1:gene.length.recur) {
			if (sws.recur + swarray.recur[i] > 0) {sws.recur = sws.recur + swarray.recur[i]} else {sws.recur = 0; swb.recur = i + 1;}
			if (sws.recur > swscore.recur) {swscore.recur = sws.recur; swstart.recur = swb.recur; swend.recur = i;}
		}

		# recursive-segment is selected
		if (swscore/sqrt(swend-swstart+1) < swscore.recur/sqrt(swend.recur - swstart.recur + 1) & gene.length.recur > 3) {
			swscore = sum(swarray[(swstart + swstart.recur - 1):(swstart + swend.recur - 1)])
			nes = (swscore-mean(swscore.random))/sd(swscore.random)
			nomP = sum(swscore.random >= swscore)/(length(swscore.random) + 1)
			mutcount = sum(mutation.positions >= (swstart + swstart.recur - 1) & mutation.positions <= (swstart + swend.recur - 1))
			write.table(cbind(as.character(target.symbol), swstart + swstart.recur - 1, swstart + swend.recur - 1, swscore, "Y", mutcount, nes, nomP), file = output.file1, row.names=F, col.names=F, quote=F, sep="\t", append = TRUE)
			swarray[(swstart + swstart.recur - 1):(swstart + swend.recur - 1)] = 0
			mutation.positions = mutation.positions[!mutation.positions %in% (swstart + swstart.recur - 1):(swend.recur + swstart - 1)]
				
		} else {
			nes = (swscore-mean(swscore.random))/sd(swscore.random)
			nomP = sum(swscore.random >= swscore)/(length(swscore.random) + 1)
			mutcount = sum(mutation.positions >= swstart & mutation.positions <= swend)
			write.table(cbind(as.character(target.symbol), swstart , swend, swscore, "N", mutcount, nes, nomP), file = output.file1, row.names=F, col.names=F, quote=F, sep="\t", append = TRUE)
			swarray[swstart:swend] = 0
			mutation.positions = mutation.positions[!mutation.positions %in% swstart:swend]

		}
	}
}
