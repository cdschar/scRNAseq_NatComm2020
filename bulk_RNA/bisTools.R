library("Biostrings");
library("BSgenome");
library("rtracklayer");
library("Rsamtools");
library("ShortRead");
library("chipseq")
library("data.table")
library("GenomeInfoDb")


# global constants
gSigDigits = 5;
gNormReads = 1e6

## bamToBigWig constants
gFragSize = NA;
gBinSize = 300;
gBwFile = "gsub('.bam', paste('.frag', paste(fragSizes, collapse = '_'), '.norm', format(normReads, scientific = FALSE), '.bw', sep = ''), bamFiles)";

gGenome = "BSgenome.Mmusculus.UCSC.mm9"

#### plotting constants
gMainSize = 1.5;
gLabSize = 1;
gAxisSize = 0.8;
gPch = 19;
gCex = 0.25
gHeight = 3;
gWidth = 3;

gSdLwd = 0.5;
gSdCol = rgb(0.5, 0.5, 0.5);

gRegLwd = 3;
gRegCol = rgb(0, 0, 0.5);
gRegCex = 2;

gMai = c(0.5, 0.5, 0.5, 0.25);
gMgp = c(1.5, 0.5, 0);

#annotates a bisulfite coverage file

annotBisCov = function(covFile, txDb, org, prefix, removeChrs = c(NA), covExt = c(".cov$"), annotExt = c(".annot.cov")) {

	#Debug	covFile = paste0(covDir, "Sample.common.coverage.", minCov, ".cov"); txDb = tx; org = org; prefix = "mm9"; removeChrs = NA; covExt = ".cov$"; annotExt = ".annot.cov";

	cov = fread(covFile, header = T);
	setnames(cov, "rname", "chr")
	setkey(cov, chr, pos)

	cov = cov[!(chr %in% removeChrs)]

	covR = GRanges(seqnames = cov$chr, ranges = IRanges(start = cov$pos, end = cov$pos))
	annot = annotBedtoTxDbGene(bed = covR, tx = txDb, prefix = prefix)

	#convert annotation to data table
	annot = as.data.table(as.data.frame(annot))
	setnames(annot, c("seqnames", "start"), c("chr", "pos"))
	setkey(annot, chr, pos)

	symbol = rep(NA, dim(annot)[1])
	symbol[!is.na(as.character(annot[[paste0(prefix, ".tsKg")]]))] = select(org, keytype = "ENTREZID", keys = as.character(annot[[paste0(prefix, ".tsKg")]])[!is.na(as.character(annot[[paste0(prefix, ".tsKg")]]))], columns = c("SYMBOL"))$SYMBOL
	annot[, eval(parse(text = paste0(prefix, ".symbol:=symbol")))]

	cov = merge(cov, annot)

	write.table(cov, gsub(covExt, annotExt, covFile), sep = "\t", quote = F, col.names = T, row.names = F)

}

#Compiles sample coverage summarizing by group and enforcing a minimium coverage (minCov) at the group level
compileBisGrpCov = function(samples, grps, covFiles, covFile = c(NA), minCov = c(10), ...) {
	#Debug samples = files$sample; grps = files$group; covFiles = paste0(files$dir, files$covFile); covFile = paste0(homeDir, "Group.common.coverage.", minCov, ".cov"); minCov = minCov;

	bedExt = ".bed";

	if (length(samples) == length(covFiles) & length(samples) == length(grps)) {
		for (i in 1:length(samples)) {

			print(paste(samples[i], Sys.time()))
			cov = readBistoolsCov(covFiles[i])

			#set appropriate names
			cov = setnames(cov, names(cov)[!names(cov) %in% key(cov)], paste(samples[i], names(cov)[!names(cov) %in% key(cov)], sep = "."))

			#aggregate
			if (i == 1) covs = cov else covs = merge(covs, cov, all = T)
			rm(cov); gc()
		}

		#Summarize each group
		for (grp in unique(grps)) {

			grpSamples = samples[grps == grp]
			covs[, eval(parse(text = paste0(grp, ".u:=rowSums(covs[, list(", paste0(grpSamples, ".u", collapse = ","), ")], na.rm = T)"))), ]
			covs[, eval(parse(text = paste0(grp, ".m:=rowSums(covs[, list(", paste0(grpSamples, ".m", collapse = ","), ")], na.rm = T)"))), ]
			covs[, eval(parse(text = paste0(grp, ".meth:=", grp, ".m/(", grp, ".u+", grp, ".m)"))), ]
		}

		#Filter summarized groups
		grpCrit	= list(); #character logical expressions to be evaluated
		for (grp in unique(grps)) grpCrit[[grp]] = paste0(paste(paste0(grp, c(".u", ".m")), collapse = "+"), ">=", minCov)
		covs = covs[eval(parse(text = paste(grpCrit, collapse = "&")))]

		#make grp bedGraphs
		for (grp in unique(grps)) writeBistoolsCovToBw(covs, file = paste0(dirname(covFile), "/", grp, ".minCov.", minCov, ".bw"), chrCol = "rname", scoreCol = paste0(grp, ".meth"), ...)

		if (!is.na(covFile)) {
			writeBistoolsCov(covs, covFile);
			writeBistoolsCovToBed(covs, paste0(covFile, bedExt))
		} else return(covs)

	} else warning("sample length is not the same as coverage files length");
}


compileBisSampleCov = function(samples, covFiles, covFile = c(NA), minCov = c(10), ...) {
	#Debug samples = files$sample; covFiles = paste0(files$dir, files$covFile); covFile = paste0(homeDir, "Sample.common.coverage.", minCov, ".cov"); minCov = minCov;

	bedExt = ".bed";
	bwExt = ".bw"

	if (length(samples) == length(covFiles)) {
		for (i in 1:length(samples)) {

			print(paste(samples[i], Sys.time()))
			cov = readBistoolsCov(covFiles[i])

			#set appropriate names & limit to minCov
			cov = setnames(cov[m + u >= minCov], names(cov)[!names(cov) %in% key(cov)], paste(samples[i], names(cov)[!names(cov) %in% key(cov)], sep = "."))

			#aggregate
			if (i == 1) covs = cov else covs = merge(covs, cov)
			rm(cov); gc()
		}

		if (!is.na(covFile)) {
			writeBistoolsCov(covs, covFile);
			writeBistoolsCovToBed(covs, paste0(covFile, bedExt))
			for (i in 1:length(samples)) writeBistoolsCovToBw(covs, file = paste0(covFile, ".", samples[i], bwExt), chrCol = "rname", scoreCol = paste0(samples[i], ".meth"), ...)


		} else return(covs)

	} else warning("sample length is not the same as coverage files length");

}

#Extract sra file
sraExtract = function(sraFile) {

	print(paste("Extracting sra files", Sys.time()))

	sraCmd = "sudo fastq-dump --split-3 "
	fq.ext = ".fastq"

	sraCmdFile = paste0(sraCmd, files$dir[i], files$sraFile[i], " -O ", files$dir[i]);
	system(sraCmdFile)

	#bashFile = "sraExt.bash"
	#write(sraCmdFile, file = bashFile, append = if (i == 1) F else T)
	# ssh in and run 'nohup ./sraExt.bash & exit'
}

#Get fastq file names from extracted sra files
getSraFastqNames = function(sraFile) {

	mate1Tag = "_1";
	mate2Tag = "_1";

	dir = dirname(sraFile);
	sraFile = basename(sraFile);

	print(paste("Getting fastq file names", sraFile, Sys.time()))

	fqExt = ".fastq";
	sraExt = ".sra";

	fqFiles = list.files(dir)
	fqMate1 = paste(fqFiles[grepl(paste0(mate1Tag, fqExt, "$"), fqFiles) & grepl(paste0("^", gsub(sraExt, "", sraFile)), fqFiles)], collapse = "|")
	fqMate2 = paste(fqFiles[grepl(paste0(mate2Tag, fqExt, "$"), fqFiles) & grepl(paste0("^", gsub(sraExt, "", sraFile)), fqFiles)], collapse = "|")
	fqFile = paste(fqFiles[grepl(paste0(fqExt, "$"), fqFiles) & grepl(paste0("^", gsub(sraExt, "", sraFile)), fqFiles) & !grepl(paste0(mate1Tag, fqExt, "$"), fqFiles) & !grepl(paste0(mate2Tag, fqExt, "$"), fqFiles)], collapse = "|")

	return(c(fqFile, fqMate1, fqMate2))

}

#mapping STAR index selection function
select_ref_STAR = function(genome="hg38") {

	if ( genome == "mm9"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/mm9.ERCC"
	} else if ( genome == "hg19"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/hg19.nohap.ERCC"
	} else if ( genome == "mm10"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/mm10.ERCC"
	} else if ( genome == "hg38"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/hg38.nohap.ERCC"
	} else if ( genome == "rn6"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/rn6.ERCC"
	} else if ( genome == "MacaM"){
		starGenome = "/Volumes/GRAID2/seqTools/STAR/MacaM.ERCC"
	} else if ( genome == "mm10.ERCC.GFP"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/mm10.ERCC.GFP"
	} else if ( genome == "mm10.Nr5a2"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/mm10.Nr5a2"
	} else if ( genome == "BALBc"){
		starGenome = "/Volumes/GRAID/seqTools/STAR/BALBc"
	}

	return(starGenome)
}

#gene counts reference txdb selection function
select_ref_txdb = function(genome="hg38") {

	if ( genome == "mm9"){
		library("TxDb.Mmusculus.UCSC.mm9.knownGene")
		tx = TxDb.Mmusculus.UCSC.mm9.knownGene
	} else if ( genome == "hg19"){
		library("TxDb.Hsapiens.UCSC.hg19.knownGene")
		tx = TxDb.Hsapiens.UCSC.hg19.knownGene
	} else if ( genome == "mm10"){
		library("TxDb.Mmusculus.UCSC.mm10.knownGene")
		tx = TxDb.Mmusculus.UCSC.mm10.knownGene
	} else if ( genome == "hg38"){
		library("TxDb.Hsapiens.UCSC.hg38.knownGene")
		tx = TxDb.Hsapiens.UCSC.hg38.knownGene
	} else if ( genome == "rn6"){
		library("TxDb.Rnorvegicus.UCSC.rn6.refGene")
		tx = TxDb.Rnorvegicus.UCSC.rn6.refGene
	} else if ( genome == "MacaM"){
		tx = makeTxDbFromGFF("/Volumes/GRAID2/seqTools/genomes/MacaM/MacaM_Rhesus_Genome_Annotation_v7.8.2.gtf", dataSource="UNMC", organism="Macaca mulatta",taxonomyId=9544)
	}

	return(tx)
}

#gene counts reference org selection function
select_ref_org = function(genome="hg38") {

	if ( genome == "mm9"){
		library("org.Mm.eg.db")
		org = org.Mm.eg.db
	} else if ( genome == "hg19"){
		library("org.Hs.eg.db")
		org = org.Hs.eg.db
	} else if ( genome == "mm10"){
		library("org.Mm.eg.db")
		org = org.Mm.eg.db
	} else if ( genome == "hg38"){
		library("org.Hs.eg.db")
		org = org.Hs.eg.db
	} else if ( genome == "rn6"){
		library("org.Rn.eg.db")
		org = org.Rn.eg.db
	} else if ( genome == "MacaM"){
		library("org.Mmu.eg.db")
		org = org.Mmu.eg.db
	}

	return(org)
}

#mapping bowtie index selection function
select_ref_bw = function(genome="hg38") {

	if ( genome == "mm9"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie-1.1.1/indexes/mm9"
	} else if ( genome == "hg19"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie-1.1.1/indexes/hg19"
	} else if ( genome == "mm10"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie-1.1.1/indexes/mm10"
	} else if ( genome == "hg38"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie-1.1.1/indexes/hg38"
	}

	return(bowtieGenome)
}

#mapping bowtie2 index selection function
select_ref_bw2 = function(genome="hg38") {

	if ( genome == "mm9"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie2.indices/mm9/mm9"
	} else if ( genome == "hg19"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie2.indices/hg19/hg19"
	} else if ( genome == "mm10"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie2.indices/mm10/mm10"
	} else if ( genome == "hg38"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie2.indices/hg38.nohap/hg38.nohap"
	} else if ( genome == "mm10_S288C"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie2.indices/mm10.S288C/mm10.S288C"
	} else if ( genome == "hg38_S288C"){
		bowtieGenome = "/Volumes/GRAID/seqTools/bowtie/bowtie2.indices/hg38.S288C/hg38.S288C"
	}

	return(bowtieGenome)
}

#mapping macs genome selection function
select_ref_macs = function(genome="hg38") {

	if ( genome == "mm9"){
		macsGenome = "mm"
	} else if ( genome == "hg19"){
		macsGenome = "hs"
	} else if ( genome == "mm10"){
		macsGenome = "mm"
	} else if ( genome == "hg38"){
		macsGenome = "hs"
	}

	return(macsGenome)
}

#create dir and move files
mvFiles = function(files, dir, mvDir, split = c("\\|")) {

	#Debug	files = paste0(files$dir[i], strsplit(files$fqFile[i], "\\|")[[1]]); dir = paste0(files$dir[i], files$sample[i], "/");

	if (!is.null(files) && !is.na(files) && files != "") {

		files = paste0(dir, strsplit(files, split)[[1]])

		if (!file.exists(mvDir)) dir.create(mvDir); #create new directory if necessary
		file.copy(files, paste0(mvDir, basename(files)))

		return(mvDir)
	}
}

formatFastq = function(fqFiles, dir, sampleName = c(NA), split = c("\\|"), threads=c(1)) {

	#Debug	fqFiles = files$fqMate1[i]; dir = files$dir[i]; sampleName = paste0(files$sample[i], "_1"); split = "\\|";
	#Debug  fqFiles = files$fqFile[i]; dir = files$dir[i]; sampleName = files$sample[i]; split = "\\|";

	if (!is.null(fqFiles) && !is.na(fqFiles) && fqFiles != "") {

		fqFiles = paste0(dir, strsplit(fqFiles, split)[[1]])

		gzExt = ".gz";
		fqExt = ".fastq";

		print(paste("Formatting fastq files", Sys.time()))

		#if multiple fastq files are provided decompress and combine for ...
		if (length(fqFiles) > 1) {

			if (is.na(sampleName)) sampleName = gsub(fqExt, "", gsub(gzExt, "", basename(fqFiles[1])))

			fqFile = paste0(dir, sampleName, fqExt); #combined fastq File name

			for (i in 1:length(fqFiles)) if (grepl(paste0(gzExt, "$"), fqFiles[i])) system(paste("pigz -p", threads, "-cd", fqFiles[i], ">", gsub(gzExt, "", fqFiles[i]))); #decompress files

			system(paste("cat", paste(gsub(gzExt, "", fqFiles), collapse = " "), ">", fqFile))

			for (i in 1:length(fqFiles)) if (grepl(paste0(gzExt, "$"), fqFiles[i])) file.remove(gsub(gzExt, "", fqFiles[i])); #delete files after concatenated

		} else {

			if (grepl(paste0(gzExt, "$"), fqFiles)) {
				system(paste("pigz -p", threads, "-cd", fqFiles, ">", gsub(gzExt, "", fqFiles))); #decompress files
				fqFile = gsub(gzExt, "", fqFiles)
			} else fqFile = fqFiles;
		}
		return(as.character(basename(fqFile)))

	} else return(NULL)

}

#check to see if fastq file is gzipped and if not compress and return name
gzipFastq = function(fqFile,  overwrite = c(F), rmIfExists = c(F), threads=8) {

	#Debug	fqFile = paste0(files$dir[i], files$fqMate1[i]);
	gz.ext = ".gz";
	fastq.ext = ".fastq"

	if (grepl(paste0(fastq.ext, "$"), fqFile) & (!file.exists(paste0(fqFile, gz.ext)) | overwrite)) {
		system(paste("pigz -p", threads, fqFile))
		return(paste0(basename(fqFile), gz.ext))
	} else if (grepl(paste0(fastq.ext, "$"), fqFile) & file.exists(paste0(fqFile, gz.ext)) & rmIfExists) {
		file.remove(fqFile)
		return(paste0(basename(fqFile), gz.ext))
	} else return(basename(fqFile))
}

#QC file
qcFastq = function(fqFile, dir, sampleName = c("")) {

	#Debug	fqFile = files$fqFile[i]; dir = files$dir[i]; sampleName = files$sample[i];

	fqExt = ".fastq"
	qual.ext = ".qual.txt";

	fxQual = "fastx_quality_stats";
	fxQualPlot = "fastq_quality_boxplot_graph.sh";
	fxBasePlot = "fastx_nucleotide_distribution_graph.sh";

	qualPlot.ext = ".qual.png";
	basePlot.ext = ".nucleotideDist.png"

	print(paste("Fastq QC", sampleName, Sys.time()))

	if (!is.null(fqFile) && !is.na(fqFile) && fqFile != "") {

		#Make Quality Plots
		system(paste0(fxQual, " -i ", dir, fqFile, " -o ", dir, gsub(fqExt, qual.ext, fqFile, fixed = T)))
		system(paste0(fxQualPlot, " -i ", dir, gsub(fqExt, qual.ext, fqFile, fixed = T), " -t ", sampleName, " -o ", dir, gsub(fqExt, qualPlot.ext, fqFile, fixed = T)))
		system(paste0(fxBasePlot, " -i ", dir, gsub(fqExt, qual.ext, fqFile, fixed = T), " -t ", sampleName, " -o ", dir, gsub(fqExt, basePlot.ext, fqFile, fixed = T)))
	}
}


#Function to quality trim fastq
qtrimFastq = function(fqFile, dir, minQ = c(30)) {

	#Debug	fqFile = files$fqFile[i]; dir = files$dir[i]; sampleName = files$sample[i];

	if (!is.null(fqFile) && !is.na(fqFile) && fqFile != "") {

		fqExt = ".fastq"
		fxTrim = "fastq_quality_trimmer";
		trimExt = paste0(".qt.", minQ, ".fastq");

		print(paste("Fastq Quality Trimmer", Sys.time()))

		#Make Quality Plots
		system(paste0(fxTrim, " -t ", minQ, " -i ", dir, fqFile, " -o ", dir, gsub(fqExt, trimExt, fqFile, fixed = T)))

		gzipFastq(paste0(dir, fqFile), rmIfExists = T)

		return(gsub(fqExt, trimExt, fqFile, fixed = T))
	}
}

#Function to quality trim fastq
trimFastq = function(fqFile, dir, f = c(NA), l = c(NA)) {

	#Debug	fqFile = files$fqFile[i]; dir = files$dir[i]; sampleName = files$sample[i];

	if (!is.null(fqFile) && !is.na(fqFile) && fqFile != "") {

		fqExt = ".fastq"
		fxTrim = "fastx_trimmer";
		trimExt = paste0(".qt.f.", f, ".l.", l, ".fastq");

		print(paste("Fastq Trimmer", Sys.time()))

		#Make Quality Plots
		if (!is.na(f)) fxTrim = paste0(fxTrim, " -f ", f)
		if (!is.na(l)) fxTrim = paste0(fxTrim, " -l ", l)

		system(paste0(fxTrim, " -i ", dir, fqFile, " -o ", dir, gsub(fqExt, trimExt, fqFile, fixed = T)))

		gzipFastq(paste0(dir, fqFile), rmIfExists = T)

		return(gsub(fqExt, trimExt, fqFile, fixed = T))
	}
}

#Function to cut adapter seqs using cutadapt
fqCutadapt = function(fqFile, dir, platform="nextera", threads=c(1)) {

	if( platform == "nextera"){
		adapter_seq = "CTGTCTCTTATA"
	} else if( platform == "illumina"){
		adapter_seq = "AGATCGGAAGAGC"
	} else if( platform == "small_rna"){
		adapter_seq = "TGGAATTCTCGG"
	}

	if( grepl(".gz", fqFile) ){
		fqname = gsub(".fastq.gz", "", fqFile)
		fout = paste0(fqname, "-trimmed.fastq.gz")
	} else {
		fqname = gsub(".fastq", "", fqFile)
		fout = paste0(fqname, "-trimmed.fastq")
	}

	system(paste0("skewer -t ",threads, " -x ", adapter_seq, " -q 3 -o ", dir,fqname, " ", dir,fqFile))
	system(paste0("rm ", dir, fqFile))
	return(fout)

}

#Function to cut adapter seqs for paired end reads
fqPECutadapt = function(fqMate1, fqMate2, dir, platform="nextera", threads=c(1)){

	if( platform == "nextera"){
		adapter_seq = "CTGTCTCTTATA"
	} else if( platform == "illumina"){
		adapter_seq = "AGATCGGAAGAGC"
	} else if( platform == "small_rna"){
		adapter_seq = "TGGAATTCTCGG"
	}

	if( grepl(".gz", fqMate1) ){
		fqname = gsub("_1.fastq.gz", "", fqMate1)
		fout1 = paste0(fqname, "-trimmed-pair1.fastq.gz")
		fout2 = paste0(fqname, "-trimmed-pair2.fastq.gz")
	} else {
		fqname = gsub("_1.fastq", "", fqMate1)
		fout1 = paste0(fqname, "-trimmed-pair1.fastq")
		fout2 = paste0(fqname, "-trimmed-pair2.fastq")
	}

	system(paste0("skewer -t ",threads, " -x ", adapter_seq, " -q 3 -o ", dir,fqname, " ", dir,fqMate1, " ", dir,fqMate2))
	system(paste0("rm ", dir, fqMate1, " ", dir, fqMate2))
	return(list(fout1, fout2))

}


#Sort bam file
sortBam = function(bamFile, bamSortFile = c(NA), delBam = c(F), threads = c(1)) {

	print(paste("sorting Bam", basename(bamFile), Sys.time()))

	idxExt = ".bai";
	bamExt = ".bam";
	bamSortExt = ".sort";

	sortCmd = "samtools sort "

	if (is.na(bamSortFile)) bamSortFile = gsub(bamExt, bamSortExt, bamFile)

	#sort BAM file
    system(paste0(sortCmd, "-o ", bamSortFile, ".bam ", bamFile, " -@ ", threads, "-m", mem))

	#Remove unsorted bam file
	if (delBam) file.remove(bamFile)
	if (delBam & file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))

	#Create BAM index for fast access
	indexBam(paste0(bamSortFile, bamExt))

	return(basename(paste0(bamSortFile, bamExt)))

}

#Use picard to mark duplicates in bam file
markDups = function(bamFile, picardCmd, delBam = c(F)) {

	print(paste("Marking Duplicates", bamFile, Sys.time()))

	idxExt = ".bai";
	bamExt = ".bam";
	bamDupMarkExt = ".dupMark.bam";

	bamDupFile = gsub(paste0(bamExt, "$"), bamDupMarkExt, bamFile)

	#Mark duplicates using picard
	system(paste0(picardCmd, "INPUT=", bamFile, " OUTPUT=", bamDupFile, " METRICS_FILE=", bamDupFile, ".metrics"))

	#Create BAM index for fast access
	indexBam(bamDupFile)

	#Remove BAM without duplicates marked
	if (delBam) file.remove(bamFile)
	if (delBam & file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))

	return(basename(bamDupFile))
}

#Get read counts for bam file
getBamCts = function(bamFile, unMappedBamFile = c(NA)) {

	print(paste("Getting Bam File Read Counts", basename(bamFile), Sys.time()))

	baiExt = ".bai";

	#create index for bamFile if it doesn't exist
	if (!file.exists(paste0(bamFile, baiExt))) indexBam(bamFile);
	if (!is.na(unMappedBamFile)) {
		if (!file.exists(paste0(unMappedBamFile, baiExt))) indexBam(unMappedBamFile);
		unmappedCt = countBam(unMappedBamFile)$records;
	} else {
		sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = T))
		unmappedCt = countBam(bamFile, param = sbp)$records;
	}

	sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F))
	mappedCt = countBam(bamFile, param = sbp)$records;

	sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F))
	uniqueCt = countBam(bamFile, param = sbp)$records;

	sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F, isPaired = T))
	pairedCt = countBam(bamFile, param = sbp)$records;

	sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F, isPaired = T, isProperPair = T))
	properPairedCt = countBam(bamFile, param = sbp)$records;

	return(c(unmappedCt, mappedCt, uniqueCt, pairedCt, properPairedCt))

}

#Function that reads Bistools table and returns % methylation data for spike in genomes
getSpikeInMeth = function(callFile = c(NA), spikeIn = c(NA)) {

	#callFile = "/Volumes/GRAID/SLE/RRBS/data/SLE969.rN/SLE969.rN.CpG.call.gz"
	#spikeIn = c("Lamda", "PhiX");

	print(paste("Calculating Spike In Methylation Values", Sys.time()))

	#spikeIn genomes chromosome IDs
	Lamda = "gi|215104|gb|J02459.1|LAMCG";
	PhiX = "gi|9626372|ref|NC_001422.1|";
	Yeast = "S288C"

	#read in call file
	call = readBistoolsCall(callFile);

	#cycle through spike in genomes and compile coverage/% meth
	for(i in 1:length(spikeIn)){

		if(spikeIn[i] == "Lamda"){
			l = call[call$rname == Lamda];
			l.reads = dim(l)[1]
			l.convert = round((1-mean(l$meth))*100, 2)
			}
		else if(spikeIn[i] == "PhiX"){
			p = call[call$rname == PhiX];
			p.reads = dim(p)[1]
			p.convert = round((1-mean(p$meth))*100, 2)
			}

	}
	return(c(l.reads, l.convert, p.reads, p.convert))
}

#Function that calculates read coverage for spike in genomes
getSpikeInCts = function(bamFile = c(NA), spikeIn = c(NA)) {

	#bamFile = "/Volumes/GRAID2/CutNRun/mouse/data/Cre1_62K_K27me3/Cre1_62K_K27me3.sort.dupMark.bam"
	#spikeIn = c("Yeast");

	print(paste("Calculating Spike In Coverage", Sys.time()))

	#spikeIn genomes chromosome IDs
	Lamda = "gi|215104|gb|J02459.1|LAMCG";
	PhiX = "gi|9626372|ref|NC_001422.1|";
	Yeast = "S288C"

	#set spikeIn
		if(spikeIn == "Lamda"){
			spike = Lamda
				}
		else if(spikeIn == "PhiX"){
			spike = PhiX
			}
		else if(spikeIn == "Yeast"){
			spike = Yeast
			}

	#get chrs list from bamfile
	si = seqinfo(BamFile(bamFile))
	cts = data.frame(chr = si@seqnames, start = 1, end = si@seqlengths)

	#filter for Spike In genome
	cts = cts[grepl(spike, cts$chr), ]

	#make GRanges object
	bed = makeGRangesFromDataFrame(cts)

	#set up parameters
	sbp = ScanBamParam(flag = scanBamFlag(isDuplicate = F), which = bed)
	reads = readGAlignmentPairs(bamFile, param = sbp)

	#Get counts
	so = summarizeOverlaps(bed, reads, mode = "Union", singleEnd = F, ignore.strand = T)

	#sumarize coverage
	spikeCov = sum(assays(so)$counts)

	return(c(spikeCov))
}


#Function that reads Bismark Call table and returns data
readBismarkCall = function(file, gz = c(NA)) {

	#read in ot and ob strand sort strand CpG calls
	covNames = c("read", "fauxStrand", "chr", "start", "meth")
	covClasses = c("NULL", "NULL", "character", "numeric", "character")

	if ((is.na(gz) & grepl("\\.gz$", file)) | gz) call = gzfile(file, "r") else call = file(file, "r")
	data = suppressWarnings(read.table(call, skip = 1, header = F, sep = "\t", comment.char = "", colClasses = covClasses, col.names = covNames))
	close(call)

	return(data)
}

#Function that reads Bistools table and returns data
readBistoolsCall = function(file, gz = c(NA), covClasses = c(NA)) {

	#file = "/home/bbarwick/Documents/Boss/Dnmt3ab/RRBS/RRBS_fasta/Dnmt3_702/t.call.gz"

	#read in ot and ob strand sort strand CpG calls
	covNames = c("readPos", "qname", "flag", "rname", "strand", "pos", "qwidth", "meth")
	if (all(is.na(covClasses))) covClasses = c("NULL", "NULL", "NULL", "character", "NULL", "numeric", "NULL", "integer")

	if ((is.na(gz) & grepl("\\.gz$", file)) | gz) file = paste("zcat <", file)
	data = fread(file, colClasses = covClasses)

	return(data)
}


#Function that reads Bistools Coverage data table
readBistoolsCov = function(file, gz = c(NA), covClasses = c(NA)) {

	#Debug
	#file = "/home/bbarwick/Documents/Boss/Dnmt3ab/RRBS/RRBS_fasta/Dnmt3_702/t.cov.gz"; covClasses = NA; gz = NA;

	#read in ot and ob strand sort strand CpG calls
	if (any(is.na(covClasses))) covClasses = c("chr", "numeric", "numeric", "integer", "integer")

	if ((is.na(gz) & grepl("\\.gz$", file)) | gz) file = paste("zcat <", file)
	cov = fread(file, colClasses = covClasses)
	setkey(cov, rname, pos)

	return(cov)
}


#Function that writes a bed formatted data.table to a bed file, potentially compressed
writeBed = function(bed, file, gz = c(NA), sep = "\t", col.names = c(F), quote = c(F), row.names = c(F)) {

	#Debug
	#cov = sum; file = paste0(covDir, grp, ".cov.", minCov, ".cov"); gz = NA; sep = "\t"

	if ((is.na(gz) & grepl("\\.gz$", file)) | (!is.na(gz) & gz)) {
		con = gzfile(file, open = "w")
		write.table(bed, con, sep = sep, row.names = row.names, col.names = col.names, quote = quote)
		close(con)
	} else write.table(bed, file, sep = sep, row.names = row.names, col.names = col.names, quote = quote)
}



#Function that writes Bistools Coverage to tab delimited txt
writeBistoolsCov = function(cov, file, gz = c(NA), sep = c("\t")) {

	#Debug
	#cov = sum; file = paste0(covDir, grp, ".cov.", minCov, ".cov"); gz = NA; sep = "\t"

	if ((is.na(gz) & grepl("\\.gz$", file)) | (!is.na(gz) & gz)) {
		con = gzfile(file, open = "w")
		write.table(cov, con, sep = sep, row.names = F)
		close(con)
	} else write.table(cov, file, sep = sep, row.names = F)
}

#Function that writes Bistools Coverage to tab delimited bed file
writeBistoolsCovToBed = function(cov, file, gz = c(NA), sep = c("\t"), chrCol = c("rname"), posCol = c("pos")) {

	#Debug
	#cov = sum; file = paste0(covDir, grp, ".cov.", minCov, ".cov"); gz = NA; sep = "\t"

	#handle gz compression
	if ((is.na(gz) & grepl("\\.gz$", file)) | (!is.na(gz) & gz)) con = gzfile(file, open = "w")

	options(scipen = 999)
	covBed = data.frame(chr = cov[[chrCol]], start = cov[[posCol]] - 1, end = cov[[posCol]], strand = rep("*", dim(cov)[1]))
	write.table(covBed, file = file, sep = sep, col.names = F, row.names = F, quote = F)

	#handle gz compression
	if ((is.na(gz) & grepl("\\.gz$", file)) | (!is.na(gz) & gz)) close(con)


}

writeBistoolsCovToBw = function(cov, file, chrCol = c("chr"), posCol = c("pos"), scoreCol = c("meth"), removeChrs = c(NA), removeBg = c(T), ...) {

	#Debug
	#cov = covs; file = "test.bw"; 	chrCol = "rname"; posCol = "pos"; scoreCol = "D30_pos.meth"; genome = genome

	bgExt = ".bedGraph";

	covBg = data.frame(chr = cov[[chrCol]], start = cov[[posCol]] - 1, end = cov[[posCol]], score = cov[[scoreCol]])
	if (any(!is.na(removeChrs))) covBg = covBg[!(covBg$chr %in% removeChrs), ]

	options(scipen = 999)
	write.table(covBg, file = paste0(file, bgExt), sep = "\t", col.names = F, row.names = F, quote = F)

	bedgraph2bw(paste0(file, bgExt), bwFile = file, ...);

	if (removeBg) file.remove(paste0(file, bgExt))
}



#Converts a BisTools Cov data.table to bedGraph
bistoolsCov2Bedgraph = function(bg, bgFile = c(NA), minCov = c(1), removeChrs = c(NA)) {

	#bg = cov; bgFile = "test.bedGraph.gz"; minCov = 10; removeChrs = c("gi|215104|gb|J02459.1|LAMCG")

	#Calculate total coverage and filter for minCov
	bg[, start:=pos-1]
	bg[, end:=pos]
	bg = bg[m+u >= minCov]

	#remove extraneous columns
	bg[, m:=NULL]
	bg[, u:=NULL]
	bg[, pos:=NULL]

	#reorder columns
	setcolorder(bg, c("rname", "start", "end", "meth"))

	#remove chromsomes contained in the removeChrs variable
	if (any(!is.na(removeChrs))) bg = bg[!(rname %in% removeChrs)]

	#sort by chromosome and position
	setkeyv(bg, c("rname", "start"))

	#if bgFile is given write to bgFile else return bg
	if (!is.na(bgFile)) {
		options(scipen=999)
		if (grepl("\\.gz$", bgFile)) bgCon = gzfile(bgFile, open = "w") else bgCon = file(bgFile, open = "w")
		write.table(bg, file = bgCon, sep = "\t", row.names = F, quote = F, col.names = F)
		close(bgCon)
	} else return(bg)
}

##Wrapper function to convert bedgraph to bigwig
bedgraph2bw = function(bgFile, genome = c(NA), seqLengths = c(NA), chromFile = c(NA), bwFile = c(gsub(".bedGraph.*$", ".bw", bgFile))) {

	#bgFile = paste0(files$dir[i], files$bgFile[i]); genome = genome; seqLengths = NA; chromFile = NA; bwFile = gsub(".bedGraph.*$", ".bw", bgFile)

	#make bigwig file. decompress bedgraph if necessary
	if (grepl(".gz$", bgFile)) system(paste("gzip -cd", bgFile, ">", gsub(".gz$", "", bgFile)))

	#write chromosome sizes to txt file for bedGraphToBigWig
	if (suppressWarnings(!is.na(genome))) {
		chromFile = paste0(genome@pkgname, ".chrom.txt")
		write.table(seqlengths(genome), file = chromFile, sep = "\t", quote = F, col.names = F)
	} else if (!is.na(seqLengths)) {
		if (is.na(chromFile)) chromFile = "chrom.txt"
		write.table(seqLengths, file = chromFile, sep = "\t", quote = F, col.names = F);
	}

	if (!is.na(chromFile)) {
		bwCmd = paste("bedGraphToBigWig", gsub(".gz$", "", bgFile), chromFile, bwFile)
		system(bwCmd)
	} else warning("No compatible chromosome sequence length data");

	if (grepl(".gz$", bgFile)) file.remove(gsub(".gz$", "", bgFile))
}

#Wrapper function to convert bedgraph to bigBed
bed2bb = function(bedFile, genome, bbFile = c(gsub(".bed.*$", ".bb", bedFile))) {

	#make bigwig file. decompress bedgraph if necessary
	if (grepl(".gz$", bedFile)) system(paste("gzip -cd", bedFile, ">", gsub(".gz$", "", bedFile)))

	#write chromosome sizes to txt file for bedGraphToBigWig
	write.table(seqlengths(genome), file = paste0(genome@pkgname, ".chrom.txt"), sep = "\t", quote = F, col.names = F)

	bbCmd = paste("bedToBigBed", gsub(".gz$", "", bedFile), paste0(genome@pkgname, ".chrom.txt"), bbFile)
	system(bbCmd)

	if (grepl(".gz$", bedFile)) file.remove(gsub(".gz$", "", bedFile))
}



#Summarizes Bismark Call data to coverage
bistoolsCall2Cov = function(call, covFile = c(NA)) {

	#use data.table to summarize average methylation 'meth', and total calls of methylated 'm' and 'u' states
	mean = data.table(call)[, list(meth = mean(as.numeric(meth))), by = "rname,pos"]
	u = data.table(call)[meth == 0, list(u = length(meth)), by = "rname,pos"]
	m = data.table(call)[meth == 1, list(m = length(meth)), by = "rname,pos"]

	#merge data
	cov = merge(mean, u, by = c("rname", "pos"), all = T)
	cov = merge(cov, m, by = c("rname", "pos"), all = T)

	#replace NA values with 0 for the m and u counts
	cov$u[is.na(cov$u)] = 0
	cov$m[is.na(cov$m)] = 0

	#if bgFile is given write to bgFile else return bg
	if (!is.na(covFile)) writeBistoolsCov(cov, covFile)

	return(cov)
}

#Summarizes Bismark Call data to coverage
bismarkCall2Cov = function(call) {

	#convert methylation calls into numeric
	call$meth[call$meth == "Z"] = 1
	call$meth[call$meth == "z"] = 0

	#use data.table too summarize average methylation 'meth', and total calls of methylated 'm' and 'u' states
	mean = data.table(call)[, list(meth = mean(as.numeric(meth))), by = "chr,start"]
	u = data.table(call)[meth == "0", list(u = length(meth)), by = "chr,start"]
	m = data.table(call)[meth == "1", list(m = length(meth)), by = "chr,start"]

	#merge data
	sum = merge(mean, u, by = c("chr", "start"), all = T)
	sum = merge(sum, m, by = c("chr", "start"), all = T)

	#replace NA values with 0 for the m and u counts
	sum$u[is.na(sum$u)] = 0
	sum$m[is.na(sum$m)] = 0

	return(sum)
}


##### Wrapper function for bismarkBam2call specifically for CpGs
bismarkBam2callCpg = function(...) bismarkBam2call(callChar = c("z"), ...)
bismarkBam2callChh = function(...) bismarkBam2call(callChar = c("x"), ...)
bismarkBam2callChg = function(...) bismarkBam2call(callChar = c("h"), ...)

##### Function to read in bamFile and make MethCalls
bismarkBam2call = function(bamFile, callFile = c(NA), trimFile = c(NA), collapseCpgFile = c(NA), trim1 = c(0, 0), trim2 = c(0, 0), callChar = c("z"), ... ) {

	#bamFile = "/home/bbarwick/Documents/Boss/Dnmt3ab/RRBS/RRBS_fasta/Dnmt3_702/ts.bam"; callFile = NA; collapseCpgFile = NA; trim1 = c(0, 0); trim2 = c(3, 0); callChar = "z";
	#bamFile = "/home/bbarwick/Documents/Boss/Dnmt3ab/RRBS/RRBS_fasta/Dnmt3_702/ts.bam"; callFile = callCollapseFile; trimFile = NA; collapseCpgFile = NA; trim1 = c(0, 0); trim2 = c(3, 0); callChar = "z";

	#bamFile = bamSortFile; callFile = callFile; trimFile = NA; collapseCpgFile = callCollapseFile; trim1 = c(0, 0); trim2 = c(0, 0); callChar = "z";

	if (!exists("bWhat") || is.na(bWhat)) bWhat = c("qname", "rname", "strand", "pos", "flag", "qwidth");

	#if the bam file is from paired end reads
	if (bamPaired(bamFile)) {

		m1t = Sys.time()
		print("Bam file is paired-end ... Processing mate 1 ... ")
		sbp1 = ScanBamParam(flag = scanBamFlag(isFirstMateRead = T), what = bWhat, tag = "XM")
		m1call = bismarkBam2callChr(bamFile, callChar = callChar, bParam = sbp1)

		print(paste("Done with mate 1 - ", format(Sys.time() - m1t)))

		#trim left and right sides of mate pair 1
		if (exists("trim1") && length(trim1) == 2 && !all(trim1 == 0)) {
			t = Sys.time()
			cat(paste("\nTrimming mate 1 by ", trim1[1], "5' bases and", trim1[2], "3' bases ..."))
			m1trim = getCallDtTrimmings(m1call, trim = trim1)
			m1call = trimCallDt(m1call, trim = trim1)
			cat(paste("Done - ", format(Sys.time() - t), "\n"))
		} else m1trim = NULL

		#Read in second mate from bam file
		m2t = Sys.time()
		print("Processing mate 2 ... ")
		sbp2 = ScanBamParam(flag = scanBamFlag(isSecondMateRead = T), what = bWhat, tag = "XM")
		m2call = bismarkBam2callChr(bamFile, callChar = callChar, bParam = sbp2)
		m2call[strand=="-",pos:=pos+1]; #adjust second mate by 1
		m2call[strand=="+",pos:=pos-1]; #adjust second mate by 1

		print(paste("Done with mate2 - ", format(Sys.time() - m2t)))

		#trim left and right sides of mate pair 1
		if (exists("trim2") && length(trim2) == 2 && !all(trim2 == 0)) {
			t = Sys.time()
			cat(paste("\nTrimming mate 2 by ", trim2[1], "5' bases and", trim2[2], "3' bases ..."))
			m2trim = getCallDtTrimmings(m2call, trim = trim2)
			m2call = trimCallDt(m2call, trim = trim2)
			cat(paste("Done - ", format(Sys.time() - t), "\n"))
		} else m2trim = NULL

		trim = rbind(m1trim, m2trim)
		methCall = rbind(m1call, m2call);

		suppressWarnings(rm(m1call, m2call, m1trim, m2trim)); gc()

	#If reads are single end
	} else {

		t = Sys.time()
		print("Bam file is single-end ... ")
		sbp = ScanBamParam(what = bWhat, tag = "XM")
		methCall = bismarkBam2callChr(bamFile, callChar = callChar, bParam = sbp)
		print(paste("Done - ", format(Sys.time() - t)))

		#trim left and right sides of mate pair 1
		if (exists("trim1") && length(trim1) == 2 && !all(trim1 == 0)) {
			t = Sys.time()
			cat(paste("\nTrimming by ", trim1[1], "5' bases and", trim1[2], "3' bases ..."))
			trim = getCallDtTrimmings(methCall, trim = trim1)
			methCall = trimCallDt(methCall, trim = trim1)
			cat(paste("Done - ", format(Sys.time() - t), "\n"))
		} else trim = NULL

		gc();
	}

	#order by chrom & pos
	methCall = methCall[order(rname, pos)]
	if (!is.null(trim)) trim = trim[order(rname, pos)]

	if (!is.na(trimFile) & !is.null(trim)) {
		if  (grepl("\\.gz$", trimFile)) {
			trimCon = gzfile(trimFile, open = "w")
			write.table(trim, file = trimCon, sep = "\t", row.names = F, quote = F)
			close(trimCon)
		} else 	write.table(trim, file = trimFile, sep = "\t", row.names = F, quote = F)
	}

	if (!is.na(callFile)) {
		if  (grepl("\\.gz$", callFile)) {
			callCon = gzfile(callFile, open = "w")
			write.table(methCall, file = callCon, sep = "\t", row.names = F, quote = F)
			close(callCon)
		} else write.table(methCall, file = callFile, sep = "\t", row.names = F, quote = F)
	}


	if (!is.na(collapseCpgFile)) {
		methCall$pos[methCall$strand == "-"] = methCall$pos[methCall$strand == "-"] - 1
		if  (grepl("\\.gz$", collapseCpgFile)) {
			collapseCpgCon = gzfile(collapseCpgFile, open = "w")
			write.table(methCall, file = collapseCpgCon, sep = "\t", row.names = F, quote = F)
			close(collapseCpgCon)
		} else write.table(methCall, file = collapseCpgFile, sep = "\t", row.names = F, quote = F)
	}

	rm(methCall); gc();
}


#Function to check if bam file is paired
bamPaired = function(bamFile, yieldSize = c(1000)) {
	#Debug;
	#bamFile = bamFile; yieldSize = 1000;

	bam = scanBam(BamFile(bamFile, yieldSize = yieldSize))[[1]]
	sum = colSums(bamFlagAsBitMatrix(bam$flag))
	if (sum[["isPaired"]] > 0) return(T) else return(F)
}


##### Helper function which checks to see if bamFile is indexed and if so reads in bam file 1 chromosome at a time converting it to a methylation calls data table
bismarkBam2callChr = function(bamFile, callChar = c("z"), bParam = c(ScanBamParam(what = c("qname", "rname", "strand", "pos", "qwidth"), tag = "XM")), ... ) {

	#Debug
	#bamFile = bamFile; callChar = "z"; bParam = sbp1

	bf = BamFile(bamFile)

	t = Sys.time()
	if (length(index(bf)) == 0 | is.na(index(bf))) {
		print(paste(bamFile, "is not indexed reading in entire bam file ..."));
		call = bismarkBam2callDt(bamFile, callChar = callChar, bParam = bParam)
	} else {
		call = NULL
		print(paste(bamFile, "is indexed reading in one chromosome at a time ..."));
		chrLen = scanBamHeader(bf)$targets
		chrs = names(chrLen)[order(chrLen)]
		for (chr in chrs) {
			cat(paste(chr, "...", "\n"))
			bamWhich(bParam) = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = chrLen[which(names(chrLen) == chr)]))
			call = rbind(call, bismarkBam2callDt(bamFile, callChar = callChar, bParam = bParam))
		}
	}
	return(call);
}

##### Helper function for converting a bamFile to methylation calls data table
bismarkBam2callDt = function(bamFile, callChar = c("z"), bParam = c(ScanBamParam(what = c("qname", "rname", "strand", "pos", "qwidth"), tag = "XM")), ... ) {

	#bamFile = bamFile; callChar = "z"; bParam = bParam

	#Read in from bam file
	bam = scanBam(bamFile, param = bParam)[[1]]

	#If there are no records return NULL else format methyl calls data table
	if (length(bam[[1]]) == 0) {
		return(NULL)
	} else {
		### Get CpG methylation Calls
		m = strsplit(gsub(paste0("[^", callChar, "]"), "", bam$tag$XM, ignore.case = T), ""); #Parse CpG Methylation Calls "Z" and "z" from bam file

		#if there are no methylation calls return NULL
		if (all(lapply(m, length) == 0)) return(NULL) else {
			mLoc = gregexpr(callChar, bam$tag$XM, ignore.case = T); #Parse CpG Methylation Calls and save locations of Z and z to list of lists
			mLoc = mLoc[unlist(lapply(mLoc, function(x) return(!all(x == -1))))]; #Remove records with no call (i.e. neither z or Z) for which gregexpr returns -1
			call = data.table(mCall = unlist(m), readPos = unlist(mLoc))

			#add all other field in bamWhat
			fields = names(bam)[names(bam) != "tag"]
			for (field in fields) call[, eval(parse(text = paste(field, ":=rep(bam[[field]], lapply(m, length))")))]

			rm(m, mLoc); gc()

			#adjust position by read position
			call[, pos:=pos+readPos-1]

			#Convert mCall (e.g. "Z" and "z") into boolean (i.e. 1 and 0)
			call[, meth:=0]
			call[mCall==toupper(callChar), meth:=1]
			call[, mCall:=NULL];		#eliminate unnecessary columns

			return(call);
		}
	}
}


##### Helper function for directionally trimming a methylation calls data table
trimCallDt = function(call, trim = c(0, 0)) {

	#Debug
	#call = m1call; trim = trim1;

	call = call[((readPos > trim[1] & strand == "+") | (qwidth - readPos >= trim[1] & strand == "-")) & ((qwidth - readPos >= trim[2] & strand == "+") | (readPos > trim[2] & strand == "-"))]
	return(call)
}



##### Helper function for getting directionally trimmed a methylation calls data table
getCallDtTrimmings = function(call, trim = c(0, 0)) {

	#Debug
	#call = m1call; trim = trim1;
	t1 = call[(readPos <= trim[1] & strand == "+") | (qwidth - readPos < trim[1] & strand == "-")]
	t2 = call[(qwidth - readPos < trim[2] & strand == "+") | (readPos <= trim[2] & strand == "-")]

	return(rbind(t1, t2))
}



#Function gets estimated fragment sizes by chromosome for a bam file and weights it by the # of reads / chromosome
getBamFragSizes = function(bamFiles) {

	#Debug
	#bamFiles = bamFile;

	fragSizes = NULL;

	for (i in 1:length(bamFiles)) {

		fragSize = NULL;
		chrReads = NULL;
		chrs = unique(getBamChrs(bamFiles[i]))

		si = seqinfo(BamFile(bamFiles[i])); #Stupid Encode Bam files have duplicated sequence names
		#si = Seqinfo(seqnames = names(scanBamHeader(bamFiles[i])[[1]][1]$targets[!duplicated(scanBamHeader(bamFiles[i])[[1]][1]$targets)]), seqlengths = scanBamHeader(bamFiles[i])[[1]][1]$targets[!duplicated(scanBamHeader(bamFiles[i])[[1]][1]$targets)])

		for (chr in chrs) {
			sbf = scanBamFlag(isDuplicate = FALSE);
			sbp = ScanBamParam(flag = sbf, which = GRanges(seqnames = chr, ranges = IRanges(start = 1, end = si@seqlengths[si@seqnames == chr])))
			bam = granges(readGAlignments(bamFiles[i], param = sbp));

			t = readGAlignments(bamFiles[i], param = sbp);

			chrReads = c(chrReads, length(bam));

			fragSize = c(fragSize, estimate.mean.fraglen(bam, method = "SISSR"));

			print(paste(bamFiles[i], chr, "estimated fragment size: ", fragSize[length(fragSize)] , Sys.time()));
		}

		#calc average fragsize weighted by chromosomal reads and add that to the list fragSizes
		fragSizes = c(fragSizes, sum(fragSize * chrReads / sum(chrReads), na.rm = TRUE));
	}

	return(fragSizes)
}

coupleFactor = function(bamFile, pattern, fragSize = c(gFragSize), binSize = c(gBinSize), genome = c(gGenome), normReads = c(gNormReads), normPat = c(NA), normWeight = c(NA)) {

	#Debug: bamFile = "/home/bbarwick/Documents/Boss/Bcell/Plasmablast.input/input.unique.bam"; pattern = "CG"; fragSize = gFragSize; binSize = gFragSize; normReads = gNormReads; genome = c(gGenome); normReads = gNormReads;
	#Debug norm: normPat = "CG"; normWeight = sampleNormWeight

	library(genome, character.only = TRUE)
	org = gsub("^[[:alnum:]]+\\.([[:alnum:]]+)\\.[[:alnum:]]+\\.[[:alnum:]]+$", "\\1", genome);	#Parse organism common name
	si = SeqinfoForBSGenome(genome); #Get genome version and seqinfo

	readCount = 0

	cfs = NULL; #couple factor of bin
	reads = NULL; # reads in bin

	for (chr in si@seqnames) {
		print(paste("Chr", chr, " start:", Sys.time()));

		sbp = ScanBamParam(which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))
		bam = suppressWarnings(resize(granges(readGAlignments(as.character(bamFile), param = sbp)), width = fragSize, fix = "start"));

		readCount = readCount + length(bam);

		#Bin chromosome
		starts = seq(from = 0, to = si@seqlengths[si@seqnames == chr] - binSize, by = binSize);
		chrBins = IRanges(start = starts, end = starts + binSize - 1)

		#Find chr patterns and calculate coupling factor of bins
		chrPat = matchPattern(pattern, eval(parse(text = org))[[chr]], fixed = FALSE);
		cfs = c(cfs, countOverlaps(chrBins, ranges(chrPat), type= "any"));

		rm(chrBins); rm(chrPat); gc();

		#if sample normalization is required
		if (!is.na(normPat) & any(!is.na(normWeight))) {
			chrNormPat = matchPattern(normPat, eval(parse(text = org))[[chr]], fixed = FALSE);
			cpl = countOverlaps(ranges(bam), ranges(chrNormPat)); #memory hog
			#Calc representative normalization
			cplNorm = normWeight[as.character(cpl), ];
			cplNorm[is.na(cplNorm)] = 1;
			#Calc weighted coverage add to existing coverage
			cv = coverage(bam, weight = cplNorm)[[chr]];
			rm(cpl); rm(cplNorm); rm(chrNormPat); gc();
		} else {
			cv = coverage(bam)[[chr]]
		}

		#Calculate running means to calculate mean read density
		means = runmean(cv, k = binSize);
		reads = c(reads, as.numeric(means)[(starts + 1)]);

		rm(means); rm(bam); gc();
	}

	#Data frame of coupling factors (cf) and reads per bin
	cfReads = data.frame(cfs = cfs, reads = reads);

	rm(cfs);
	rm(reads);
	gc();

	#Mean # of reads in bins w/ a specific coupling factor
	mean = aggregate(cfReads$reads, by = list(cfReads$cfs), mean); #Capture the sum of read counts per cf

	#Variance of the number of reads in bins w/ a specific coupling factor
	var = aggregate(cfReads$reads, by = list(cfReads$cfs), var); #Capture the variance of the read counts per cf

	#Calculate the mean and standard deviation
	return(data.frame(cf = mean$Group.1, mean = mean$x * (normReads / readCount), sd = sqrt(var$x) * (normReads / readCount)));
}

plotCoupleFactor = function(data, binSize, normReads, pattern, output, outFile, yMax = c(NULL), xMax = c(NULL), title = c(NULL), plotReg = c(TRUE), mainSize = c(gMainSize), labSize = c(gLabSize), axisSize = c(gAxisSize), pch = c(gPch), cex = c(gCex), plotSd = c(TRUE), sdLwd = c(gSdLwd), sdCol = c(gSdCol), regLwd = c(gRegLwd), regCol = c(gRegCol), regCex = c(gRegCex), height = c(gHeight), width = c(gWidth), mai = gMai, mgp = gMgp, ...) {

	#Debug
	#mainSize = 2.5; labSize = 2; axisSize = 1.5; pch = 19; cex = 1; sdLwd = 1; sdCol = rgb(0.5, 0.5, 0.5); regLwd = 3; reg.col = rgb(0, 0, 0.5); regCex = 2; sigDigits = 3; yMax = NULL; xMax = NULL; output = "screen"; title = "Debug"
	#data = cfs; binSize = binSize; normReads = normReads; pattern = pattern; output = output; outFile = outFile; yMax = yMax; xMax = NULL; title = samples$sample[i]; plotReg = TRUE; mainSize = gMainSize; labSize = gLabSize; axisSize = gAxisSize; pch = gPch; cex = gCex; sdLwd = gSdLwd; sdCol = gSdCol; regLwd = gRegLwd; regCol = gRegCol; regCex = gRegCex;

	if(is.null(list(...)$sigDigits)) sigDigits = gSigDigits else sigDigits = list(...)$sigDigits;

	#Calculate first local max
	lm = localMax(data$mean, trend = 3);
	reg = lm(data$mean[1:lm] ~ data$cf[1:lm]);

	if (is.null(yMax)) yMax = max(c(max(data$mean) * reg$coef[[2]] + reg$coef[[1]], max(data$mean, na.rm = TRUE) + max(data$sd, na.rm = TRUE)));

	if (is.null(xMax)) xMax = max(data$cf);

	if (output == "pdf") cairo_pdf(outFile, family = "Arial", height = height, width = width);

	par(mai = mai, mgp = mgp, bty = "n")
	plot(data$cf, data$mean, main = title, xlab = paste(pattern, "coupling Factor"), ylab = paste("Reads per ", normReads, " per ", binSize, "bp bin", sep = ""), pch = pch, cex = cex, cex.main = mainSize, cex.lab = labSize, cex.axis = axisSize, ylim = c(0, yMax), xlim = c(0, xMax))
	if (plotSd) for (cf in data$cf) lines(c(cf, cf), c(data$mean[which(data$cf == cf)] - data$sd[which(data$cf == cf)], data$mean[which(data$cf == cf)] + data$sd[which(data$cf == cf)]), lwd = sdLwd, col = sdCol);

	if (plotReg) {
		abline(reg, lwd = regLwd, col = regCol);
		legendText = paste("y = ", round(reg$coef[[2]], sigDigits), "x", ifelse(reg$coef[[1]] > 0, " + ", " - "), round(reg$coef[[1]], sigDigits), sep = "");
		legend("topleft", legendText, lwd = regLwd, col = regCol, bty = "n", cex = regCex)
	}

	if (output != "screen" | output != "X11") dev.off();
}

readComp = function(bamFile, pat, normReads = c(gNormReads), fragSize = c(NA), genome = c(gGenome), normPat = c(NA), normWeight = c(NA), isDup = c(NA)) {

	#Debug: bamFile = "/home/bbarwick/Documents/Boss/Bcell/MeDIPseq/Bcell1/Bcell1.MeDIP.sort.dupMark.bam"; pat = "CG"; fragSize = NA; normReads = gNormReads; yieldSize = NA; genome = c(gGenome); normPat = NA; normWeight = NA; isDup = FALSE;
	#Debug norm: normPat = "CG"; normWeight = sampleNormWeight

	library(genome, character.only = TRUE)
	org = gsub("^[[:alnum:]]+\\.([[:alnum:]]+)\\.[[:alnum:]]+\\.[[:alnum:]]+$", "\\1", genome);	#Parse organism common name
	si = SeqinfoForBSGenome(genome); #Get genome version and seqinfo

	readCount = 0
	comp = NA;

	for (chr in getBamChrs(bamFile)) {
		print(paste("Chr", chr, " start:", Sys.time()));

		sbf = scanBamFlag(isDuplicate = isDup);
		sbp = ScanBamParam(flag = sbf, which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))

		if (is.na(fragSize)) {
			bam = granges(readGAlignments(bamFile, param = sbp));
		} else {
			bam = suppressWarnings(resize(granges(readGAlignments(bamFile, param = sbp)), width = fragSize, fix = "start"));
		}

		readCount = readCount + length(bam);

		#Find chr patterns and calculate coupling factor of bins
		chrPat = matchPattern(pat, eval(parse(text = org))[[chr]], fixed = FALSE);
		chrComp = as.data.frame(table(countOverlaps(ranges(bam), ranges(chrPat), type = "any")), responseName = chr);

		if (all(is.na(comp))) comp = chrComp else comp = merge(comp, chrComp, by = "Var1", all = TRUE);

		rm(chrPat); rm(bam); gc();


		#if sample normalization is required
		#if (!is.na(normPat) & any(!is.na(normWeight))) {
		#	chrNormPat = matchPattern(normPat, eval(parse(text = org))[[chr]], fixed = FALSE);
		#	cpl = countOverlaps(ranges(bam), ranges(chrNormPat)); #memory hog
			#Calc representative normalization
		#	cplNorm = normWeight[as.character(cpl), ];
		#	cplNorm[is.na(cplNorm)] = 1;
			#Calc weighted coverage add to existing coverage
		#	cv = coverage(bam, weight = cplNorm)[[chr]];
		#	rm(cpl); rm(cplNorm); rm(chrNormPat); gc();
		#} else {
		#	cv = coverage(bam)[[chr]]
		#}

	}

	rownames(comp) = comp$Var1
	colnames(comp)[1] = paste("comp", pat, sep = ".");

	compSum = data.frame(comp = comp[, 1], reads = apply(comp[, -1], 1, sum, na.rm = TRUE) * normReads / readCount);
	compSum = compSum[order(compSum$comp), ];

	#Calculate the mean and standard deviation
	return(compSum);
}


writeCf = function(cfs, file) write.table(cfs, file, sep = ",", row.names = FALSE, quote = FALSE);

#Read coupling factor file
readCf = function(file) return(read.table(file, sep = ",", header = TRUE));

localMax = function(x, trend) {

	inc.trend = 0; #Times the trend has increased
	inc = FALSE; #Is the trend increasing

	for (i in 2:length(x)) {
		if (x[i] > x[i-1]) {
			inc.trend = inc.trend + 1
			dec.trend = 0
		} else {
			inc.trend = 0;
			dec.trend = dec.trend + 1;
		}

		if (inc.trend >= trend) inc = TRUE;
		if (inc & dec.trend >= trend) return(i - trend);
	}
}

#Function to return chromosomes included in a Bam file
getBamChrs = function(bamFile) return(names(scanBamHeader(bamFile)[[1]][1]$targets));

getBamChrLengths = function(bamFile) {

	#Debug: bamFile = "/home/bbarwick/Documents/Boss/Bcell/MeDIPseq/Bcell1/Bcell1.MeDIP.unique.bam";

	h = scanBamHeader(bamFile)
	cl = h[[1]][2]$text[names(h[[1]][2]$text) == "@SQ"]

	chrs = gsub("SN:", "", unlist(cl)[grepl("SN:", unlist(cl))]);
	len = as.numeric(gsub("LN:", "", unlist(cl)[grepl("LN:", unlist(cl))]));
	names(len) = chrs;

	return(len);
}



#bamToBigWig function: converts a bam file into a bigWig coverage file.  Removes chromosomes that match the variable removeChrs by grepl.  Normalizes for read count .Only reads not in removeChrs are counted. The normReads variable sets the sequencing coverage. The default normalization is reads per million (e.g. normRead = 1e6).
bamToBigWig = function(bamFiles, bwFile = c(gBwFile), fragSizes = c(gFragSize), sigDigits = c(gSigDigits), normReads = c(gNormReads), extCall = c(T), flag = scanBamFlag(), removeChrs = c(NA), ...) {

	#Debug
	#bamFiles = bamFile; bwFile = "test.bw"; fragSizes = gFragSize; sigDigits = gSigDigits; normReads = gNormReads; extCall = T; flag = scanBamFlag(); removeChrs = c("random");

	#Get genome version and seqinfo
	si = list(); for (i in 1:length(bamFiles)) si[[i]] = seqinfo(BamFile(bamFiles[i]))

	#if there are more bamFiles than fragSizes then recycle fragSizes
	if (!all(is.na(fragSizes)) & length(bamFiles) > length(fragSizes)) fragSizes = rep(fragSizes, ceiling(length(bamFiles) / length(fragSizes)));

	#Get readcounts for each bamFile
	#readCounts = list(); for (i in 1:length(bamFiles)) readCounts[[i]] = countBam(bamFiles[i], param = ScanBamParam(flag = flag))$records

	#Create bigwig and bedfile output file names
	if (bwFile == gBwFile) bwFile = eval(parse(text = gBwFile))
	bgFile = paste0(bwFile, ".bedGraph")

	#get all chromosomes in bamFiles
	chrs = unique(unlist(lapply(si, seqnames)))

	#remove chromosomes that match removeChrs
	if (!is.na(removeChrs)) chrs = chrs[!grepl(removeChrs, chrs)]

	#coverage and read count variables
	cv = RleList(compress = F);
	ct = list(); for (bamFile in bamFiles) ct[[bamFile]] = 0

	#Recurse through each  chromosome
	for (chr in chrs) {
		print(paste(chr, Sys.time()))

		#Recursively build coverage for each bamFile
		for (i in 1:length(bamFiles)) {

			#Check to make sure the bam file has reads on chromosome chr
			if (chr %in% si[[i]]@seqnames) {

				#Set up chr specific bam files
				sbp = ScanBamParam(flag = flag, which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si[[1]]@seqlengths[si[[1]]@seqnames == chr])))

				#Read in chromosome from each bam file
				bam = readGAlignments(as.character(bamFiles[i]), param = sbp);

				#Extend reads to fragment size
				if (!is.na(fragSizes[i])) bam = suppressWarnings(resize(granges(bam), width = fragSizes[i], fix = "start"));

				#cumulatively add up coverages and reads
				if (i == 1) cv[[chr]] = coverage(bam)[[chr]] else cv[[chr]] = cv[[chr]] + coverage(bam)[[chr]];

				#add up read counts
				ct[[bamFiles[i]]] = ct[[bamFiles[i]]] + length(bam)

				#clear memory
				rm(bam); gc();
			}
		}
	}

	#normalize coverage to reads
	for (chr in chrs) cv[[chr]] = round(cv[[chr]] * normReads / sum(unlist(ct)), sigDigits)

	#if external calll write to bedgraph then call Kent tools
	if(extCall) {
		export.bedGraph(cv, bgFile)
		writeSiChromSizes(si[[1]], "temp.chrom.size")
		system(paste("bedGraphToBigWig", bgFile, "temp.chrom.size", bwFile))
		file.remove(bgFile); file.remove("temp.chrom.size");
	} else {
		export.bw(rleCv, bwFile)
	}
}

writeSiChromSizes = function(si, file) write.table(cbind(si@seqnames, si@seqlengths), file = file, sep = "\t", quote = F, row.names = F, col.names = F)

#### function to get bam pileup values

#### Constants for Hist functions
gHistBin = 10;
gHistOutFile = "paste(bamFiles, '.', gsub('^.*/', '', bedFile), '.range', range, '.bin', bin, '.csv', sep = '')"

makeBamHist = function(bamFiles, bed, range, outFile = c(gHistOutFile), bin = c(gHistBin), fragSizes = c(gFragSize), normReads = c(gNormReads), sigDigits = c(gSigDigits), sampleNames = c(NA), flag = c(NA), readCounts = c(NA), ...) {

	#Debug
	#bamFiles = paste0(files$dir[i], files$bamFile[i]); bed = proms; range = 100; outFile = outFile; bin = 10; fragSizes = c(gFragSize); normReads = gNormReads; sampleNames = NA; sigDigits = gSigDigits; flag = NA; readCounts = NA;

	bins = seq(from = -range, to = range, by = bin)

	if (!is.null(bed@elementMetadata$name) & length(bed@elementMetadata$name) == length(unique(bed@elementMetadata$name))) {
		histRows = bed@elementMetadata$name
	} else {
		histRows = paste(as.character(bed@seqnames), start(bed), sep = "_");
	}

	if (outFile == gHistOutFile)  outFile = eval(parse(text = outFile))

	bedCenters = round(apply(cbind(start(bed), end(bed)), MARGIN = 1, mean))

	#Make hist matrix index
	histIdx = matrix(bins, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);
	histIdx = histIdx + bedCenters
	histIdx[histIdx <= 0] = NA

	#make hist matrix
	hist = matrix(0, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);

	hists = list();
	for (i in 1:length(bamFiles)) hists[[i]] = hist;

	if (all(is.na(readCounts))) {
		readCounts = list();
		for (i in 1:length(bamFiles)) readCounts[[i]] = countBam(bamFiles[i], param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F)))$records
	}

	for (chr in levels(bed@seqnames)) {

		print(paste("chr", chr, Sys.time()))

		#Recursively build coverage for each bamFile
		for (i in 1:length(bamFiles)) {

			si = seqinfo(BamFile(bamFiles[i]))

			#Check to make sure the bam file has reads on chromosome chr
			if (chr %in% unique(getBamChrs(as.character(bamFiles[i])))) {

				bedChr = bed[bed@seqnames == chr];

				#Read in chromosome from each bam file, extend reads to fragment size
				sbp = ScanBamParam(which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))

				if (is.na(fragSizes[i])) {
					bam = granges(readGAlignments(as.character(bamFiles[i]), param = sbp));
				} else 	bam = suppressWarnings(resize(granges(readGAlignments(as.character(bamFiles[i]), param = sbp)), width = fragSizes[i], fix = "start"));

				#Calculate base level coverage
				cv = coverage(bam)[[chr]];

				#Calculate the running mean to put into bins
				means = runmean(cv, k = bin)

				#Get the names of rows for chromosome chr
				chrRows = rownames(hists[[i]])[as.character(seqnames(bed)) == chr];

				#Replace regions that are larger than than the genome with NA
				hists[[i]][chrRows, ][hists[[i]][chrRows, ] > si@seqlengths[si@seqnames == chr]] = NA

				#Replace hist rows with mean
				hists[[i]][chrRows, ] = as.numeric(means)[histIdx[chrRows, ]]

				rm(bam); rm(cv); rm(means); gc();
			}
		}
	}

	#add up hists for each bam file
	histSum = NA;
	for (i in 1:length(bamFiles)) if (all(is.na(histSum))) histSum = hists[[i]] else histSum = histSum + hists[[i]];

	#Adjust for strandness
	histSum[as.character(bed@strand) == "-", ] = histSum[as.character(bed@strand) == "-", rev(seq_len(ncol(histSum)))];

	#Normalize
	histAvg = histSum * normReads / sum(unlist(readCounts));

	histAvg = cbind(row = row.names(histAvg), round(histAvg, sigDigits))

	if (!is.na(outFile)) write.table(histAvg, file = outFile, sep = ",", row.names = F, quote = F) else return(histSum)

	rm(histSum); rm(histAvg); gc();
}


makeBamHistAvg = function(bamFiles, bed, range, outFile = c(gHistOutFile), bin = c(gHistBin), fragSizes = c(gFragSize), normReads = c(gNormReads), sigDigits = c(gSigDigits), sampleNames = c(NA), flag = c(NA), readCounts = c(NA), ...) {

	#Debug
	#bamFiles = paste0(files$dir[i], files$bamFile[i]); bed = proms; range = 100; outFile = outFile; bin = 10; fragSizes = c(gFragSize); normReads = norm; sampleNames = NA; sigDigits = gSigDigits; flag = NA; readCounts = NA;

	bins = seq(from = -range, to = range, by = bin)

	if (!is.null(bed@elementMetadata$name) & length(bed@elementMetadata$name) == length(unique(bed@elementMetadata$name))) {
		histRows = bed@elementMetadata$name
	} else {
		histRows = paste(as.character(bed@seqnames), start(bed), sep = "_");
	}

	if (outFile == gHistOutFile)  outFile = eval(parse(text = outFile))

	bedCenters = round(apply(cbind(start(bed), end(bed)), MARGIN = 1, mean))

	#Make hist matrix index
	histIdx = matrix(bins, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);
	histIdx = histIdx + bedCenters
	histIdx[histIdx <= 0] = NA

	#make hist matrix
	hist = matrix(0, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);

	hists = list();
	for (i in 1:length(bamFiles)) hists[[i]] = hist;

	if (all(is.na(readCounts))) {
		readCounts = list();
		for (i in 1:length(bamFiles)) readCounts[[i]] = countBam(bamFiles[i], param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F)))$records
	}

	for (chr in levels(bed@seqnames)) {

		print(paste("chr", chr, Sys.time()))

		#Recursively build coverage for each bamFile
		for (i in 1:length(bamFiles)) {

			si = seqinfo(BamFile(bamFiles[i]))

			#Check to make sure the bam file has reads on chromosome chr
			if (chr %in% unique(getBamChrs(as.character(bamFiles[i])))) {

				bedChr = bed[bed@seqnames == chr];

				#Read in chromosome from each bam file, extend reads to fragment size
				sbp = ScanBamParam(which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))

				if (is.na(fragSizes[i])) {
					bam = granges(readGAlignments(as.character(bamFiles[i]), param = sbp));
				} else 	bam = suppressWarnings(resize(granges(readGAlignments(as.character(bamFiles[i]), param = sbp)), width = fragSizes[i], fix = "start"));

				#Calculate base level coverage
				cv = coverage(bam)[[chr]];

				#Calculate the running mean to put into bins
				means = runmean(cv, k = bin)

				#Get the names of rows for chromosome chr
				chrRows = rownames(hists[[i]])[as.character(seqnames(bed)) == chr];

				#Replace regions that are larger than than the genome with NA
				hists[[i]][chrRows, ][hists[[i]][chrRows, ] > si@seqlengths[si@seqnames == chr]] = NA

				#Replace hist rows with mean
				hists[[i]][chrRows, ] = as.numeric(means)[histIdx[chrRows, ]]

				rm(bam); rm(cv); rm(means); gc();
			}
		}
	}

	#add up hists for each bam file
	histSum = NA;
	for (i in 1:length(bamFiles)) if (all(is.na(histSum))) histSum = hists[[i]] else histSum = histSum + hists[[i]];

	#Adjust for strandness
	histSum[as.character(bed@strand) == "-", ] = histSum[as.character(bed@strand) == "-", rev(seq_len(ncol(histSum)))];

	#Normalize
	histSum = histSum * normReads / sum(unlist(readCounts));

	histAvg = colMeans(histSum, na.rm = T)

	write.table(round(histAvg, sigDigits), file = outFile, sep = ",", quote = F, col.names = F)

	rm(histSum); gc();
}

makeBamHistScaled = function(bamFiles, bed, range, outFile = c(gHistOutFile), bin = c(gHistBin), fragSizes = c(gFragSize), normReads = c(gNormReads), sigDigits = c(gSigDigits), sampleNames = c(NA), flag = c(NA), readCounts = c(NA), ...) {

	#Debug
	#bamFiles = bamFile; bed = motifr; range = 100; outFile = NA; bin = 1; fragSizes = c(gFragSize); normReads = gNormReads; sampleNames = NA; sigDigits = gSigDigits; flag = NA; readCounts = NA;
	#normReads = 1e6 / files$frip[i]; fragSizes = 1; outFile = statFile

	#Determine mean width of bed range
	meanWidth = round(mean(width(bed)))

	#binned region to calculate enrichment
	bins = seq(from = -range, to = range + meanWidth, by = bin)

	#Create bin labels, normalize distance within bed region
	binlabs = bins
	binlabs[binlabs > 0 & binlabs <= meanWidth] = binlabs[binlabs > 0 & binlabs <= meanWidth] / meanWidth
	binlabs[binlabs > meanWidth] = binlabs[binlabs > meanWidth] - meanWidth + 1

	if (!is.null(bed@elementMetadata$name) & length(bed@elementMetadata$name) == length(unique(bed@elementMetadata$name))) {
		histRows = bed@elementMetadata$name
	} else histRows = paste(as.character(bed@seqnames), start(bed), sep = "_");

	if (outFile == gHistOutFile)  outFile = eval(parse(text = outFile))

	#bedCenters = round(apply(cbind(start(bed), end(bed)), MARGIN = 1, mean))

	#Make hist matrix index
	histIdx = matrix(bins, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);
	histIdx = histIdx + start(bed)
	histIdx[histIdx <= 0] = NA

	#make hist matrix
	hist = matrix(0, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);

	hists = list();
	for (i in 1:length(bamFiles)) hists[[i]] = hist;

	if (all(is.na(readCounts))) {
		readCounts = list();
		for (i in 1:length(bamFiles)) readCounts[[i]] = countBam(bamFiles[i], param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F)))$records
	}

	for (chr in levels(bed@seqnames)) {

		print(paste("chr", chr, Sys.time()))

		#Recursively build coverage for each bamFile
		for (i in 1:length(bamFiles)) {

			si = seqinfo(BamFile(bamFiles[i]))

			#Check to make sure the bam file has reads on chromosome chr
			if (chr %in% unique(getBamChrs(as.character(bamFiles[i])))) {

				bedChr = bed[bed@seqnames == chr];

				#Read in chromosome from each bam file, extend reads to fragment size
				sbp = ScanBamParam(which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))

				if (is.na(fragSizes[i])) {
					bam = granges(readGAlignments(as.character(bamFiles[i]), param = sbp));
				} else 	bam = suppressWarnings(resize(granges(readGAlignments(as.character(bamFiles[i]), param = sbp)), width = fragSizes[i], fix = "start"));

				#Calculate base level coverage
				cv = coverage(bam)[[chr]];

				#Calculate the running mean to put into bins
				means = runmean(cv, k = bin)

				#Get the names of rows for chromosome chr
				chrRows = rownames(hists[[i]])[as.character(seqnames(bed)) == chr];

				#Replace regions that are larger than than the genome with NA
				hists[[i]][chrRows, ][hists[[i]][chrRows, ] > si@seqlengths[si@seqnames == chr]] = NA

				#Replace hist rows with mean
				hists[[i]][chrRows, ] = as.numeric(means)[histIdx[chrRows, ]]

				rm(bam); rm(cv); rm(means); gc();
			}
		}
	}

	#add up hists for each bam file
	histSum = NA;
	for (i in 1:length(bamFiles)) if (all(is.na(histSum))) histSum = hists[[i]] else histSum = histSum + hists[[i]];

	#Adjust for strandness
	histSum[as.character(bed@strand) == "-", ] = histSum[as.character(bed@strand) == "-", rev(seq_len(ncol(histSum)))];

	colnames(histSum) = binlabs; #label bins correctly
	#Normalize and round
	histAvg = histSum * normReads / sum(unlist(readCounts));
	histAvg = cbind(row = row.names(histAvg), round(histAvg, sigDigits))

	if (!is.na(outFile)) write.table(histAvg, file = outFile, sep = ",", row.names = F, quote = F)
	invisible(histSum)

	#rm(histSum); rm(histAvg); gc();
}

gListOutFile = "paste(bwFile, '.', gsub('^.*/', '', bedFile), '.range', range, '.csv', sep = '')"

makeBwScoreList = function(bwFile, bed, range = c(NA), outFile = c(gListOutFile), sigDigits = c(gSigDigits), flag = c(NA), ...) {

	#Debug
	#bwFile = bwFile; bed = ts; range = NA; outFile = paste0(homeDir, files$sample[i], ".", bedPrefix, ".txt");
	#bwFile = paste0(conDir, conFile); bed = cov; range = 100; outFile = "test.txt"; sigDigits = gSigDigits; flag = NA
	#bwFile = bwFile; bed = deg; outFile = paste0(homeDir, files$sample[i], ".", "mm9.txt"); sigDigits = gSigDigits; flag = NA

	if (is.character(bed)) bed = unique(import.bed(bed)) else if (!inherits(bed, what = "GRanges")) print("ERROR: bed not bed file or GRanges object");

	if (outFile == gListOutFile) outFile = eval(parse(text = outFile))

	scores = data.frame(chr = NULL, start = NULL, score = NULL, relLoc = NULL)

	si = seqinfo(BigWigFile(bwFile))

	chrs = intersect(si@seqnames, unique(seqnames(bed)))

	for (chr in chrs) {

		print(paste("chr", chr, Sys.time()))

		bedChr = bed[bed@seqnames == chr];

		if (!is.na(range)) {
			iRanges = IRanges(start(bedChr) - range, end(bedChr) + range)
			start(iRanges)[start(iRanges) < 1] = 1
			end(iRanges)[end(iRanges) > si@seqlengths[si@seqnames == chr]] = si@seqlengths[si@seqnames == chr]
		} else iRanges = IRanges(1, si@seqlengths[si@seqnames == chr]);

		bw = import.bw(as.character(bwFile), which = GRanges(seqnames = as.character(chr), ranges = iRanges, seqinfo = si))

		if (length(bw) > 0) {
			scoresChr = data.frame(chr = chr, start = start(bw), score = round(score(bw), sigDigits), relStart = NA, relEnd = NA, relLoc = NA)

			seqlevels(bw) = chr
			seqlevels(bedChr) = chr;
			dist = distanceToNearest(bw, bedChr)

			bedChrPos = as.character(strand(bedChr[dist@subjectHits])) == "+"

			scoresChr$relStart[bedChrPos] = start(bw)[bedChrPos] - start(bedChr[dist@subjectHits])[bedChrPos];
			scoresChr$relStart[!bedChrPos] = end(bedChr[dist@subjectHits])[!bedChrPos] - start(bw)[!bedChrPos];

			scoresChr$relEnd[bedChrPos] = start(bw)[bedChrPos] - end(bedChr[dist@subjectHits])[bedChrPos];
			scoresChr$relEnd[!bedChrPos] = start(bedChr[dist@subjectHits])[!bedChrPos] - start(bw)[!bedChrPos];

			#set relLoc (relative location) to distance. Scale bw elements with bed files from 0 to 1
			beforeBed = scoresChr$relStart <= 0
			scoresChr$relLoc[beforeBed] = scoresChr$relStart[beforeBed];

			#if bw location is within the bed file then scale the relLoc from 0 (start) to 1 (end) directionally
			inBed = scoresChr$relStart > 0 & scoresChr$relEnd < 0
			scoresChr$relLoc[inBed] = scoresChr$relStart[inBed] / (scoresChr$relStart[inBed] - scoresChr$relEnd[inBed]);

			afterBed = scoresChr$relEnd >= 0
			scoresChr$relLoc[afterBed] = scoresChr$relEnd[afterBed];

			scores = rbind(scores, scoresChr)

			rm(bw); rm(scoresChr); gc();
		}
	}

	write.table(scores, outFile, sep = "\t", row.names = F)
	return(scores)
}


annotBedtoTxDbGene = function(bed, tx, org = c(NULL), prefix = c(NULL), promUp = 2500, promDown = 0, ...) {

	#Debug
	#bed = testR; tx = tx; prefix = "mm9"; promUp = 2500; promDown = 0; org = org
	#bed = motifr; tx = tx; org = org; prefix = "mm9"; promUp = 2500; promDown = 0;

	#Columns to add to bed Granges
	annotCols = paste(prefix, c("promoter", "exon", "intron", "5utr", "3utr", "ts", "tsKg", "tsDist", "tsStart", "tsEnd", "tsStrand", "tsRelPos", "tss", "tssKg", "tssDist", "tssRelPos", "tssStart", "tssEnd", "tssStrand"), sep = ".");
	bed@elementMetadata[annotCols] = rep(NA, dim(bed@elementMetadata)[1])

	#get Ts database
	ts = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")); #, ...); #get transcripts

	#Determine the the closest transcript
	closeTs = distanceToNearest(bed, ts, ignore.strand = T)

	#Assign ts name to annotation column (AnnotCol) "tx"
	bed@elementMetadata[[paste0(prefix, ".ts")]][closeTs@queryHits] = ts@elementMetadata$TXNAME[closeTs@subjectHits]
 	bed@elementMetadata[paste0(prefix, ".tsKg")] = select(tx, keytype = "TXNAME", keys = as.character(bed@elementMetadata[[paste0(prefix, ".ts")]]), columns = "GENEID")$GENEID

	bed@elementMetadata[[paste0(prefix, ".tsStart")]][closeTs@queryHits] = start(ts)[closeTs@subjectHits]
	bed@elementMetadata[[paste0(prefix, ".tsEnd")]][closeTs@queryHits] = end(ts)[closeTs@subjectHits]
	bed@elementMetadata[[paste0(prefix, ".tsStrand")]][closeTs@queryHits] = as.character(strand(ts)[closeTs@subjectHits])

	 #Determine the orientation
	sign = sign(start(bed[closeTs@queryHits]) - start(ts[closeTs@subjectHits]))
	sign[as.character(strand(ts[closeTs@subjectHits])) == "-"] = -sign[as.character(strand(ts[closeTs@subjectHits])) == "-"]; #adjust for strand

	#Assign ts dist
	bed@elementMetadata[[paste(prefix, "tsDist", sep = ".")]][closeTs@queryHits] = closeTs@elementMetadata$distance * sign

	#Determine ts rel pos
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][closeTs@queryHits] = closeTs@elementMetadata$distance * sign

	#Determine if the bed is within or after a transcript, if so scale it to the transcript
	inTx = bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][closeTs@queryHits] == 0
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][closeTs@queryHits][inTx] = (start(bed[closeTs@queryHits][inTx]) - start(ts[closeTs@subjectHits])[inTx]) / (end(ts[closeTs@subjectHits])[inTx] - start(ts[closeTs@subjectHits])[inTx])
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][closeTs@queryHits][inTx & as.character(strand(ts[closeTs@subjectHits])) == "-"] = (1 - bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][closeTs@queryHits][inTx & as.character(strand(ts[closeTs@subjectHits])) == "-"])

	tss = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")); #, ...);
	#make strand-specific ts range object
	end(tss)[as.character(strand(tss)) == "+"] = start(tss)[as.character(strand(tss)) == "+"]
	start(tss)[as.character(strand(tss)) == "-"] = end(tss)[as.character(strand(tss)) == "-"]

	#Determine the relative location
	closeTss = distanceToNearest(bed, tss, ignore.strand = T)
	#closeTss = closeTss[!is.na(closeTss@subjectHits), ]

	#Assign ts name to annotation column (AnnotCol) "ts"
	bed@elementMetadata[[paste0(prefix, ".tss")]][closeTs@queryHits] = tss@elementMetadata$TXNAME[closeTss@subjectHits]
	bed@elementMetadata[[paste0(prefix, ".tssKg")]] = select(tx, keytype = "TXNAME", keys = bed@elementMetadata[[paste0(prefix, ".tss")]], columns = "GENEID")$GENEID

	bed@elementMetadata[[paste0(prefix, ".tssStart")]][closeTss@queryHits] = start(ts)[closeTss@subjectHits]
	bed@elementMetadata[[paste0(prefix, ".tssEnd")]][closeTss@queryHits] = end(ts)[closeTss@subjectHits]
	bed@elementMetadata[[paste0(prefix, ".tssStrand")]][closeTss@queryHits] = as.character(strand(ts)[closeTss@subjectHits])

	#Determine the orientation
	sign = sign(start(bed[closeTss@queryHits]) - start(tss[closeTss@subjectHits]))
	sign[as.character(strand(tss[closeTss@subjectHits])) == "-"] = -sign[as.character(strand(tss[closeTss@subjectHits])) == "-"]; #adjust for strand

	#Assign tss dist
	bed@elementMetadata[[paste0(prefix, ".tssDist")]][closeTss@queryHits] = closeTss@elementMetadata$distance * sign

	#Determine the position relative to the closest tss
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits] = closeTss@elementMetadata$distance * sign

	#Determine if the bed is within or after a transcript, if so scale it to the transcript
	inTx = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits] >= 0 & bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits] < (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][closeTss@queryHits] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][closeTss@queryHits])
	afterTx = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits] >= (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][closeTss@queryHits] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][closeTss@queryHits])

	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits][inTx] = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits][inTx] / (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][closeTss@queryHits][inTx] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][closeTss@queryHits][inTx])
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits][afterTx] = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][closeTss@queryHits][afterTx] - bed@elementMetadata[[paste0(prefix, ".tssEnd")]][closeTss@queryHits][afterTx] + bed@elementMetadata[[paste0(prefix, ".tssStart")]][closeTss@queryHits][afterTx]

	#Assign promoter overlap
	proms = suppressWarnings(promoters(tx, upstream = promUp, downstream = promDown, columns = c("GENEID")))
	overlaps = findOverlaps(bed, proms)
	overlaps = cbind(as.data.frame(overlaps), id = as.character(proms@elementMetadata$GENEID[overlaps@subjectHits]))
	overlaps = data.table(unique(overlaps[, c(1, 3)]))
	overlaps = overlaps[!is.na(overlaps$id), ]
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "promoter", sep = ".")]][ag$queryHits] = ag$id

	#Assign exon overlap
	exons = exonsBy(tx, by = "tx")
	if (exists("vals") && all(!is.na(vals$tx_name))) exons = exons[names(exons) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset Exon")
	overlaps = findOverlaps(bed, exons)
	overlaps = cbind(as.data.frame(overlaps), id = names(exons)[overlaps@subjectHits]);
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXID)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[, c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "exon", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign intron overlap
	introns = intronsByTranscript(tx, use.name = T)
	if (exists("vals") && all(!is.na(vals$tx_name))) introns = introns[names(introns) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset Intron")
	overlaps = findOverlaps(bed, introns)
	overlaps = cbind(as.data.frame(overlaps), id = names(introns)[overlaps@subjectHits]);
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "intron", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign 5'utr classification
	fiveUtr = fiveUTRsByTranscript(tx, use.name = T)
	if (exists("vals") && all(!is.na(vals$tx_name))) fiveUtr = fiveUtr[names(fiveUtr) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset 5' UTR")
	overlaps = findOverlaps(bed, fiveUtr)
	overlaps = cbind(as.data.frame(overlaps), id = names(fiveUtr)[overlaps@subjectHits]);
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "5utr", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign 3'utr classification
	threeUtr = threeUTRsByTranscript(tx, use.name = T)
	if (exists("vals") && all(!is.na(vals$tx_name))) threeUtr = threeUtr[names(threeUtr) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset 3' UTR")
	overlaps = findOverlaps(bed, threeUtr)
	overlaps = cbind(as.data.frame(overlaps), id = names(threeUtr)[overlaps@subjectHits]);
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "3utr", sep = ".")]][ag$queryHits] = as.character(ag$id)

	if (!is.null(org)) {
		bed@elementMetadata[paste0(prefix, ".tsSymbol")] = rep(NA, length(bed))
		#bed@elementMetadata[[paste0(prefix, ".tsSymbol")]][!is.na(as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]]))] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]])[!is.na(as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]]))], columns = c("SYMBOL"))$SYMBOL
		bed@elementMetadata[[paste0(prefix, ".tsSymbol")]][closeTs@queryHits] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]][closeTs@queryHits]), columns = c("SYMBOL"))$SYMBOL

		bed@elementMetadata[paste0(prefix, ".tssSymbol")] = rep(NA, length(bed))
		#bed@elementMetadata[[paste0(prefix, ".tssSymbol")]][!is.na(as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]]))] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[[paste0(prefix, ".tssKg")]])[!is.na(as.character(bed@elementMetadata[[paste0(prefix, ".tssKg")]]))], columns = c("SYMBOL"))$SYMBOL
		bed@elementMetadata[[paste0(prefix, ".tssSymbol")]][closeTss@queryHits] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[[paste0(prefix, ".tssKg")]][closeTss@queryHits]), columns = c("SYMBOL"))$SYMBOL

	}

	return(bed);
}



annotBedtoGffStart = function(bed, gff, prefix = c(NULL)) {

	#Debug
	#bed = bed; gff = gff; prefix = "test"; fileName = NA;

	#Columns to add to bed Granges
	gffCols = paste(prefix, c("start", "end", "strand", "Name", "distStart", "distEnd"), sep = ".");

	bed@elementMetadata[gffCols] = rep(NA, dim(bed@elementMetadata)[1])

	# get the chromosomes of the bed bed file
	chrs = as.character(seqnames(bed)@values);

	for (chr in chrs) {

		print(paste(chr, Sys.time()))

		closeGffs = NULL;

		#limit data sets to relevant chromsome chr
		bedChr = bed[bed@seqnames == chr];
		gffChr = gff[gff@seqnames == chr];

		#set up strand specific start / ends
		gffStarts = start(gffChr);
		gffEnds = end(gffChr);
		gffStarts[as.character(strand(gffChr)) == "-"] = end(gffChr)[as.character(strand(gffChr)) == "-"];
		gffEnds[as.character(strand(gffChr)) == "-"] = start(gffChr)[as.character(strand(gffChr)) == "-"];

		#set up bed starts / ends
		bedStarts = start(bedChr);
		bedStarts[as.character(strand(bedChr)) == "-"] = end(bedChr)[as.character(strand(bedChr)) == "-"];
		bedEnds = end(bedChr);
		bedEnds[as.character(strand(bedChr)) == "-"] = start(bedChr)[as.character(strand(bedChr)) == "-"];

		centers = (start(bedChr) + end(bedChr)) / 2;

		for (i in 1:length(bedChr)) closeGffs = c(closeGffs, which(abs(centers[i] - gffStarts) == min(abs(centers[i] - gffStarts)))[1]);

		bed@elementMetadata[bed@seqnames == chr, gffCols[1]] = gffStarts[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[2]] = gffEnds[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[3]] = as.character(strand(gffChr[closeGffs]))
		bed@elementMetadata[bed@seqnames == chr, gffCols[4]] = gffChr@elementMetadata$Name[closeGffs];
		bed@elementMetadata[bed@seqnames == chr, gffCols[5]] = (bedStarts - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))
		bed@elementMetadata[bed@seqnames == chr, gffCols[6]] = (bedEnds - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))

	}

	return(bed);
}

annotBedtoGffEnd = function(bed, gff, prefix = c(NULL)) {

	#Debug
	#bed = bed; gff = gff; prefix = "test"; fileName = NA;

	#Columns to add to bed Granges
	gffCols = paste(prefix, c("start", "end", "strand", "Name", "distStart", "distEnd"), sep = ".");

	bed@elementMetadata[gffCols] = rep(NA, dim(bed@elementMetadata)[1])

	# get the chromosomes of the bed bed file
	chrs = as.character(seqnames(bed)@values);

	for (chr in chrs) {

		print(paste(chr, Sys.time()))

		closeGffs = NULL;

		#limit data sets to relevant chromsome chr
		bedChr = bed[bed@seqnames == chr];
		gffChr = gff[gff@seqnames == chr];

		#set up strand specific start / ends
		gffStarts = start(gffChr);
		gffEnds = end(gffChr);
		gffStarts[as.character(strand(gffChr)) == "-"] = end(gffChr)[as.character(strand(gffChr)) == "-"];
		gffEnds[as.character(strand(gffChr)) == "-"] = start(gffChr)[as.character(strand(gffChr)) == "-"];

		#set up bed starts / ends
		bedStarts = start(bedChr);
		bedStarts[as.character(strand(bedChr)) == "-"] = end(bedChr)[as.character(strand(bedChr)) == "-"];
		bedEnds = end(bedChr);
		bedEnds[as.character(strand(bedChr)) == "-"] = start(bedChr)[as.character(strand(bedChr)) == "-"];

		centers = (start(bedChr) + end(bedChr)) / 2;

		for (i in 1:length(bedChr)) closeGffs = c(closeGffs, which(abs(centers[i] - gffEnds) == min(abs(centers[i] - gffEnds)))[1]);

		bed@elementMetadata[bed@seqnames == chr, gffCols[1]] = gffStarts[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[2]] = gffEnds[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[3]] = as.character(strand(gffChr[closeGffs]))
		bed@elementMetadata[bed@seqnames == chr, gffCols[4]] = gffChr@elementMetadata$Name[closeGffs];
		bed@elementMetadata[bed@seqnames == chr, gffCols[5]] = (bedStarts - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))
		bed@elementMetadata[bed@seqnames == chr, gffCols[6]] = (bedEnds - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))

	}

	return(bed);
}

#Function to annotate a bed object to another bed object
annotBedtoBed = function(bedA, bed, prefix = c(NULL)) {

	#Debug
	#bedA = covR; bed = bedAnnot; prefix = "cpgi";

	chrs = as.character(intersect(seqnames(bed), seqnames(bedA)));
	bedACom = bedA[seqnames(bedA) %in% chrs]
	bedCom = bed[seqnames(bed) %in% chrs]

	seqlevels(bedACom) = chrs;
	seqlevels(bedCom) = chrs;

	#Columns to add to bed Granges
	bedCols = paste(prefix, c("start", "end", "strand", "name", "dist", "relDist"), sep = ".");
	bedACom@elementMetadata[bedCols] = rep(NA, dim(bedACom@elementMetadata)[1])

	closeBed = distanceToNearest(bedACom, bedCom)

	#set the closest
	bedACom@elementMetadata[, bedCols[1]] = start(bedCom)[subjectHits(closeBed)]
	bedACom@elementMetadata[, bedCols[2]] = end(bedCom)[subjectHits(closeBed)]
	bedACom@elementMetadata[, bedCols[3]] = as.character(strand(bedCom)[subjectHits(closeBed)])

	if (!is.null(bedCom@elementMetadata$name)) bedACom@elementMetadata[, bedCols[4]] = bedCom@elementMetadata$name[subjectHits(closeBed)];

	bedACom@elementMetadata[, bedCols[5]] = closeBed@elementMetadata$distance

	#orient distance relative to location and strand (e.g. upstream or downstream)
	upstream = (bedACom@elementMetadata[, bedCols[1]] > end(bedACom) & (bedACom@elementMetadata[bedCols[3]][[1]] == "*" | bedACom@elementMetadata[bedCols[3]][[1]] == "+")) | (bedACom@elementMetadata[, bedCols[2]] < start(bedACom) & bedACom@elementMetadata[bedCols[3]][[1]] == "-")
	bedACom@elementMetadata[upstream, bedCols[5]] = -bedACom@elementMetadata[upstream, bedCols[5]]

	#Set the relative distance
	bedACom@elementMetadata[, bedCols[6]] = bedACom@elementMetadata[, bedCols[5]]
	overlapPos = bedACom@elementMetadata[, bedCols[6]] == 0 & (bedACom@elementMetadata[, bedCols[3]] == "+" | bedACom@elementMetadata[bedCols[3]][[1]] == "*")
	overlapNeg = bedACom@elementMetadata[, bedCols[6]] == 0 & bedACom@elementMetadata[, bedCols[3]] == "-"
	bedACom@elementMetadata[overlapPos, bedCols[6]] = (start(bedACom)[overlapPos] - bedACom@elementMetadata[overlapPos, bedCols[1]]) / (bedACom@elementMetadata[overlapPos, bedCols[2]] - bedACom@elementMetadata[overlapPos, bedCols[1]])
	bedACom@elementMetadata[overlapNeg, bedCols[6]] = (bedACom@elementMetadata[overlapNeg, bedCols[2]] - start(bedACom)[overlapNeg]) / (bedACom@elementMetadata[overlapNeg, bedCols[2]] - bedACom@elementMetadata[overlapNeg, bedCols[1]])

	bedA@elementMetadata[bedCols] = rep(NA, dim(bedA@elementMetadata)[1])
	for (col in bedCols) bedA@elementMetadata[seqnames(bedA) %in% chrs, col] = bedACom@elementMetadata[col][[1]]

	return(bedA);
}

gMetadata = c("type", "Name");

bedOverlapGff = function(bed, gff, prefix = c(NULL), metadata = gMetadata) {

	#Debug
	#bed = bed; gff = gff; prefix = "test"; metadata = gMetadata;

	#Columns to add to bed Granges
	bedCols = paste(prefix, metadata, sep = ".");
	bed@elementMetadata[bedCols] = rep(NA, dim(bed@elementMetadata)[1])

	laps = findOverlaps(bed, gff)

	for (i in unique(laps@queryHits)) for (j in 1:length(bedCols)) bed@elementMetadata[i, bedCols[j]] = paste(gff@elementMetadata[metadata[j]][[1]][laps@subjectHits[laps@queryHits == i]], collapse = ";")

	return(bed);
}


#Function to consolidate bed intersection
bedInt = function(bed1, bed2, outFile = c(NA)) {

	if (class(bed1) == "character") bed1 = import.bed(bed1)
	if (class(bed2) == "character") bed2 = import.bed(bed2)

	int = intersect(bed1, bed2)

	score1 = aggregate(bed1@elementMetadata$score[findOverlaps(bed1, int)@queryHits], by = list(findOverlaps(bed1, int)@subjectHits), FUN = mean, na.rm = T)$x
	score2 = aggregate(bed2@elementMetadata$score[findOverlaps(bed2, int)@queryHits], by = list(findOverlaps(bed2, int)@subjectHits), FUN = mean, na.rm = T)$x

	score = apply(cbind(score1, score2), MARGIN = 1, mean, na.rm = TRUE)
	values(int) = data.frame(score = score);

	if (!is.na(outFile)) write.table(int, outFile = outFile, sep = "\t", row.names = F);

	return(int);
}



#Function to consolidate bed union
bedUnion = function(bed1File, bed2File) {

	#Debug

	bed1 = import.bed(bed1File)
	bed2 = import.bed(bed2File)

	bed = union(bed1, bed2)

	score1 = rep(NA, length(bed));
	score2 = rep(NA, length(bed));

	score1[findOverlaps(bed1, bed)@subjectHits] = bed1@elementMetadata$score
	score2[findOverlaps(bed2, bed)@subjectHits] = bed2@elementMetadata$score

	score = apply(cbind(score1, score2), MARGIN = 1, mean, na.rm = TRUE)

	name1 = rep(NA, length(bed));
	name2 = rep(NA, length(bed));

	name1[findOverlaps(bed1, bed)@subjectHits] = bed1@elementMetadata$name
	name2[findOverlaps(bed2, bed)@subjectHits] = bed2@elementMetadata$name

	name = paste(name1, name2, sep = "_")

	values(bed) = data.frame(score = score, name = name);

	return(bed);
}


#Function to consolidate bed union
exportBed = function(bed, con) {

	bed = amp; con = "test.bed"

	#Debug
	out = as.data.frame(bed);
	if (any(grepl("strand", names(out)))) out$strand = gsub("\\*", "\\.", out$strand)
	write.table(out, con, sep = "\t", quote = F, col.names = F, row.names = F)
}

#function to generate sam file for HOMER that is missing chromosome data
removeChrSamFile = function(bamFile, dir, removeChrs = c(NA)) {
	#bhanu edits

	bamFile = paste0(dir,bamFile)
	#Get genome version and seqinfo
	si =  seqinfo(BamFile(bamFile))
	#get all chromosomes in bamFiles
	chrs = unique(seqnames(si))
	#remove chromosomes that match removeChrs
	if (!is.na(removeChrs)) chrs = chrs[!grepl(removeChrs, chrs)]
	chrs = noquote(chrs)
	chrList = paste(chrs, collapse = " ")

	tempbamSortFile = bamFile
	samFile = paste0(dir,"rmChr.sam")
	system(paste("samtools view -h ", tempbamSortFile, chrList, " >", samFile))
	return(samFile)

}

#Funtion to filter for detected genes based on groups
filter_detected = function(x, group, cutoff){

	one_gene = data.frame(group=group, rpm=as.vector(t(x)))
	one_gene = aggregate(rpm~group, data=one_gene, FUN=function(x) sum(x>cutoff))

	if( any( one_gene$rpm >= table(group) ) )
		{ return(TRUE) }
	else
		{ return(FALSE) }
}

#Function to summarize reads for each chromosome
readsByChr = function(bamFile){

	samtoolsCall = paste("samtools idxstats", bamFile, "> temp.txt")
	system(samtoolsCall)
	counts = read.table("temp.txt", sep="\t", row.names=1)
	counts = counts[grepl("chr", row.names(counts)), ]
	colnames(counts) = c("length", "mapped", "unmapped")
	counts = t(counts)
	if (file.exists("temp.txt")) file.remove("temp.txt")

	return(counts)

}

#####annot function for peak data and will limit annotation to a custom RNA file
annotPeak = function(peak, tx, org = c(NULL), rna = c(NULL), prefix = c(NULL), promUp = 2000, promDown = 500, both = FALSE, ...) {

	#Debug
	#peak = data; tx = tx; prefix = "DEG_"; promUp = 2500; promDown = 0; org = org; rna = rna; both = FALSE;

	#establish GRanges object
	bed = GRanges(seqnames = peak$chr, ranges = IRanges(start = peak$start, end = peak$end), strand = Rle(peak$strand), peak = peak$peak)

	if(is.null(rna) | both == TRUE) {
		#Columns to add to bed Granges / kg is the ENTREZID and needs updating
		annotCols = c("ENTREZID", "SYMBOL", "GeneName", "UCSCID", "ENTREZ_TSS_dist");
		bed@elementMetadata[annotCols] = rep(NA, dim(bed@elementMetadata)[1])
		#changed column names
		#tsKg = ENTREZID
		#ts = UCSCID
		#tsRelPos = ENTREZ_TSS_dist
		#tsSymbol = SYMBOL

		#get Ts database
		ts = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")); #, ...); #get transcripts

		#Determine the the closest transcript
		closeTs = distanceToNearest(bed, ts, ignore.strand = T)

		#Assign ts name to annotation column (AnnotCol) "tx"
		bed@elementMetadata[["UCSCID"]][queryHits(closeTs)] = ts@elementMetadata$TXNAME[subjectHits(closeTs)]
		bed@elementMetadata["ENTREZID"] = select(tx, keytype = "TXNAME", keys = as.character(bed@elementMetadata[["UCSCID"]]), columns = "GENEID")$GENEID

		#Determine the orientation
		sign = sign(start(bed[queryHits(closeTs)]) - start(ts[subjectHits(closeTs)]))
		sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"] = -sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"]; #adjust for strand

		#Determine the relative location
		tss = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")); #, ...);
		#make strand-specific ts range object
		end(tss)[as.character(strand(tss)) == "+"] = start(tss)[as.character(strand(tss)) == "+"]
		start(tss)[as.character(strand(tss)) == "-"] = end(tss)[as.character(strand(tss)) == "-"]

		closeTss = distanceToNearest(bed, tss, ignore.strand = T)

		#Determine the position relative to the closest tss
		bed@elementMetadata[["ENTREZ_TSS_dist"]][queryHits(closeTss)] = closeTss@elementMetadata$distance * sign

		#annnotate SYMBOL
		if (!is.null(org)) {
			bed@elementMetadata["SYMBOL"] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[["ENTREZID"]]), columns = c("SYMBOL"))$SYMBOL
			bed@elementMetadata["GeneName"] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[["ENTREZID"]]), columns = "GENENAME")$GENENAME
		}

		#if no other annotation return bed
		if(is.null(rna) & both == FALSE) {
		return(bed); 
		}
	} 
	
	#if RNA limit annotation to RNA file
	if(!is.null(rna)) {
	
		annotCols = c(paste0(prefix, "ENTREZID"), paste0(prefix, "SYMBOL"), paste0(prefix, "GeneName"), paste0(prefix, "UCSCID"), paste0(prefix, "ENTREZ_TSS_dist"));
		bed@elementMetadata[annotCols] = rep(NA, dim(bed@elementMetadata)[1])
		#changed column names
		#tsKg = ENTREZID
		#ts = UCSCID
		#tsRelPos = ENTREZ_TSS_dist
		#tsSymbol = SYMBOL

		#get Ts database and limit to ENTREZIDS in RNA table
		ts = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID"), filter = list(gene_id = rna$ENTREZID)); #get transcripts

		#Determine the the closest transcript
		closeTs = distanceToNearest(bed, ts, ignore.strand = T)

		#Assign ts name to annotation column (AnnotCol) "tx"
		bed@elementMetadata[[paste0(prefix, "UCSCID")]][queryHits(closeTs)] = ts@elementMetadata$TXNAME[subjectHits(closeTs)]
		bed@elementMetadata[paste0(prefix, "ENTREZID")] = select(tx, keytype = "TXNAME", keys = as.character(bed@elementMetadata[[paste0(prefix, "UCSCID")]]), columns = "GENEID")$GENEID

		#Determine the orientation
		sign = sign(start(bed[queryHits(closeTs)]) - start(ts[subjectHits(closeTs)]))
		sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"] = -sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"]; #adjust for strand

		#Determine the relative location and limit to ENTREZIDs in RNA table
		tss = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID"), filter = list(gene_id = rna$ENTREZID));
		#make strand-specific ts range object
		end(tss)[as.character(strand(tss)) == "+"] = start(tss)[as.character(strand(tss)) == "+"]
		start(tss)[as.character(strand(tss)) == "-"] = end(tss)[as.character(strand(tss)) == "-"]

		closeTss = distanceToNearest(bed, tss, ignore.strand = T)

		#Determine the position relative to the closest tss
		bed@elementMetadata[[paste0(prefix, "ENTREZ_TSS_dist")]][queryHits(closeTss)] = closeTss@elementMetadata$distance * sign

		#annnotate SYMBOL
		if (!is.null(org)) {
			bed@elementMetadata[paste0(prefix, "SYMBOL")] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[[paste0(prefix, "ENTREZID")]]), columns = c("SYMBOL"))$SYMBOL
			bed@elementMetadata[paste0(prefix, "GeneName")] = select(org, keytype = "ENTREZID", keys = as.character(bed@elementMetadata[[paste0(prefix, "ENTREZID")]]), columns = "GENENAME")$GENENAME
		}

		return(bed); 
	}
}
