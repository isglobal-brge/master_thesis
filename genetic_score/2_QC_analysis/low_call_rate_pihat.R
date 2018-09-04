pihat <- read.table("pihat_min0.2_in_founders.genome", header =T)
miss <- read.table("plink.imiss", header =T)
# Read data into R

ind_rel <- data.frame()
for(i in 1:nrow(pihat))
{
  # For each pair of related individuals we save the pihat and missingness values
	i1 <- pihat[i,2]
	i2 <- pihat[i,4]

	i1_miss <- miss[grep(paste("^",i1,"$",sep=""), miss[,2]), 6]
	i2_miss <- miss[grep(paste("^",i2,"$",sep=""), miss[,2]), 6]
  
	# Then we compare them to see which individual have lower call rate and we save those individuals to remove them
	
	if(i1_miss > i2_miss)
		{
			ind_rel <- rbind(ind_rel, miss[grep(paste("^",i2,sep=""), miss[,2]), 1:2])
		}
	if(i1_miss < i2_miss)
		{
			ind_rel <- rbind(ind_rel, miss[grep(paste("^",i1,sep=""), miss[,2]), 1:2])
		}
}
write.table(ind_rel, file="0.2_low_call_rate_pihat.txt", col.names =T, row.names=F)
