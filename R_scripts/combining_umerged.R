#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/r/3.5.0/bin/Rscript

args = commandArgs(TRUE)
start = args[1] # Specify which sequences in "list_ind" file you want to align, directly from the shell. Alternatively, you can do this from the "alignments" function itself.

setwd(start)

###files to work with
system(paste("ls -1 *unassembled.forward.fastq >umerged.forward.files",sep = ""))
system(paste("ls -1 *unassembled.reverse.fastq >umerged.reverse.files",sep = ""))
system(paste("ls -1 *.assembled.fastq >assembled.files",sep = ""))

files_f = read.table("umerged.forward.files",stringsAsFactors = F)
files_r = read.table("umerged.reverse.files",stringsAsFactors = F)
files_a = read.table("assembled.files",stringsAsFactors = F)


for(i in 1:nrow(files_f))
  {
#grep 1st and 2rd lines
  	unassembled.fastq = paste("awk 'NR % 2 == 1' ",files_f[i,1]," >unassembled.names.seq.fastq",sep="")
  	system(unassembled.fastq)

#grep sequences forward...
  	unassembled.forward.fastq = paste("awk 'NR % 2 == 0' ",files_f[i,1]," >unassembled.seq.forward.fastq",sep="")
  	system(unassembled.forward.fastq)

#grep sequences reverse...
  	unassembled.reverse.fastq = paste("awk 'NR % 2 == 0' ",files_r[i,1]," >unassembled.seq.reverse.fastq",sep="")
  	system(unassembled.reverse.fastq)

#N and E quality file
 	system("wc -l unassembled.seq.reverse.fastq >wc")
 	wc = read.table("wc")

 	write.table(c(rbind(rep('NNNNNNNNNNNNNNNNNNNN',wc[1,1]/2),rep('EEEEEEEEEEEEEEEEEEEE',wc[1,1]/2))),"NE_file",row.names = F, col.names = F, quote = F)

#paste the unassembled sequences into a single file
	system("paste -d '\\0' unassembled.seq.forward.fastq NE_file unassembled.seq.reverse.fastq >unassembledN")

#put everything back into a single fastq
	back = paste("paste -d '\\n' unassembled.names.seq.fastq unassembledN >",gsub("merged.unassembled.forward","cat",files_f[i,1]),sep= "")
	system(back)

#add the merged sequences
	all = paste("cat ",files_a[i,1]," ",gsub("merged.unassembled.forward","cat",files_f[i,1])," >",gsub("merged.assembled","merged.assembled2",files_a[i,1]),sep = "")
	system(all)

#remove the clutter (file specific)
	remove = c("rm NE_file unassembledN wc unassembled.names.seq.fastq unassembled.seq.forward.fastq unassembled.seq.reverse.fastq")
	system(remove)
}

#remove the clutter (listing of all files)
system("rm assembled.files umerged.forward.files umerged.reverse.files")
