options(stringsAsFactors = F)

q_score <- read.table("q_score_file.txt", header=T)
ranges <- read.table("ranges.txt")

files <- rev(ranges$V1)

numbers <- c()
for (i in files)
{
	name <- paste("snp_score_",i,".txt", sep="")
	number <- sub("S0*","",sub('*_(.*)',"", i))

	infos  <- tryCatch(read.table(name), error = function(e) NULL)
	if (is.null(infos))
	{
		print(paste("Error: ", i ,"Failed Reading"))

	}
	else
	{
		numbers <- c(numbers, number)
		q_score$P[q_score$SNP %in% infos$V2] <- number
	}

}

write.table(q_score, "q_score_to_new_data.txt", col.names=T, row.names=F)

new_ranges <- data.frame(files, rep(0,length(files)),c(numbers,1))

write.table(new_ranges, "ranges_to_new_data.txt", row.names=F, col.names=F)