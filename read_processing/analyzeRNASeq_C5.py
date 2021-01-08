## same RNA-seq pipeline specific to using C5 reference genome fasta and gff
import glob, sys, os, string, datetime
now = datetime.datetime.now()
timeprint = now.strftime("%Y-%m-%d %H:%M")

# Create genome index
def genomeIndex(genomeDir,starDir,genomeFasta):
    indexCmd = '%s --runMode genomeGenerate --runThreadN 4 --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s  --genomeSAindexNbases 7 --sjdbOverhang 74 --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript locus_tag' %(starDir, genomeDir, genomeFasta, genomeGFF)
    print indexCmd

    print "\033[34m %s Indexing genome... \033[0m" %(timeprint)
    if os.path.exists('%s/SAindex' %(genomeDir)):
        print 'Genome indexes exist. Not creating!'
    else:
        print 'Creating genome indexes'
    #os.system(indexCmd)


# Function to run and control pipeline
def runPipeline():
    fastqFolders = glob.glob('%s/*/' %(dataDir)) # Modify for location of your read folders
    print
    print 'FASTQ Folders: %s' %fastqFolders
    print
    # Start processing each folder
    for folder in fastqFolders:
        sampleFolder = folder.split("/")[-2] # Folder containing reads
        print(folder)
        print(sampleFolder)
        fastqFilesFirst = glob.glob('%s/*_R1*.fastq.gz' %(folder)) # Get list of 1st pair files
        #print 'FASTQ Files: %s' %fastqFilesFirst
        # Start processing each file
        for file in fastqFilesFirst:
            print
            print "\033[33m Processing %s \033[0m" %(file)
            fileName = file.split("_R1")[0] # Read file name

            print "FileName: %s" %fileName
            sampleResultsDir = resultsDir+ '/'+sampleFolder # folder name in the results folder

            print "sampleResultsDir: %s" %sampleResultsDir
            sampleTitle = fileName.split("/")[-2]

            libraryName = sampleTitle.split("_")[0]

            lane = fileName.split("/")[-1].split("_")[-1]

            sampleName = fileName.split("/")[-1]

            #sampleId = fileName.split("/")[2].split("_")[0].split("-")[3] # uncomment for actual run
            sampleId = fileName.split("/")[-1].split("_")[0] # uncomment for test run
            print "FileName: %s Lane: %s sample: %s ID: %s Library name: %s" %(fileName, lane, sampleName, sampleId, libraryName)
            print

            # Create read names
            firstPair = fileName + "_R1_001.fastq.gz" # or "_R1_001.fastq.gz"
            secondPair = fileName + "_R2_001.fastq.gz"

            print "First Pair: %s, Second Pair: %s" %(firstPair, secondPair)
            print
            #sys.exit()
            ############ Run Functions ##############
            # Run Fastqc
            #runQC(firstPair, secondPair)
            # Run trimmomatic
            #irstPairTrimmedPaired, secondPairTrimmedPaired = runTrim(firstPair, secondPair)
            # Run trimgalore
            firstPairTrimmedPaired, secondPairTrimmedPaired = trimgalore(firstPair, secondPair,sampleName,fileName)

        # Run STAR
        runStar(firstPairTrimmedPaired, secondPairTrimmedPaired, libraryName, sampleResultsDir, folder)

        # Run HTSeq count
        runHTseq(libraryName, sampleResultsDir)



    return firstPair, secondPair, libraryName, sampleResultsDir, folder

# Quality control
def trimgalore(firstPair, secondPair,sampleName,fileName):
    print
    print "\033[34m %s Running TrimGalore \033[0m" %(timeprint)

    # define result files
    print firstPair
    filesFolder = fileName.split('/%s' %sampleName)[0]
    print filesFolder
    firstPairTrimmedPaired = filesFolder+"/trimmed/"+ sampleName + "/" + firstPair.split('/')[-1].split(".fastq.gz")[0] + "_val_1.fq.gz"

    secondPairTrimmedPaired = filesFolder+"/trimmed/" + sampleName + "/" + secondPair.split('/')[-1].split(".fastq")[0] + "_val_2.fq.gz"
    trimDir = filesFolder + "/trimmed/"+ sampleName + "/"
    fqcDir = filesFolder + "/fastqc/"+ sampleName + "/"
    # create fastqc folder
    if not os.path.exists('%s' %(fqcDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(fqcDir)
        os.makedirs('%s' %(fqcDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(fqcDir)

    # create results folder
    if not os.path.exists('%s' %(trimDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(trimDir)
        os.makedirs('%s' %(trimDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(trimDir)
    # run Command and write output into both screen and logfile with 2>&1 | tee -a %s
    cmd = 'trim_galore --fastqc_args "--outdir %s" --paired --output_dir %s %s %s' %(fqcDir, trimDir,firstPair, secondPair)
    print
    print 'FastQC Command:', cmd
    #os.system(cmd)

    print firstPairTrimmedPaired
    print secondPairTrimmedPaired
    print
    return firstPairTrimmedPaired, secondPairTrimmedPaired



# Run STAR alignment
def runStar(firstPairTrimmedPaired, secondPairTrimmedPaired, libraryName, sampleResultsDir, folder):
    print
    print "\033[34m %s Running STAR alignment... \033[0m" %(timeprint)
    # collect list of recalibrated bam files
    print folder
    # We need to comma seperate all firt read pairs and then second read pairs for ReadFilesIn
    readFilesInFirst = glob.glob('%strimmed/*/*_R1_*_val_1.fq.gz' %(folder))
    print(readFilesInFirst)

    readFilesInSecond = []
    # Get the first read pair filename and replace R1 with R2 to create second read file name
    for file in readFilesInFirst:
        secondPair = file.replace("_R1_", "_R2_")
        secondPair = secondPair.replace("_val_1", "_val_2")
        readFilesInSecond.append(secondPair)
    readFilesIn1stJoined = ",".join(readFilesInFirst)
    readFilesIn2ndJoined = ",".join(readFilesInSecond)

    # create results folder
    if not os.path.exists('%s' %(sampleResultsDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(sampleResultsDir)
        os.makedirs('%s' %(sampleResultsDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(sampleResultsDir)

    starOptions ="--runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --alignSJoverhangMin 5 --alignSJDBoverhangMin 3 --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --quantMode TranscriptomeSAM GeneCounts"
    outFileNamePrefix = sampleResultsDir+"/"+libraryName+"_star_"

    cmd = '%s --genomeDir %s %s --readFilesIn %s %s --outFileNamePrefix %s' %(starPath, genomeDir, starOptions, readFilesIn1stJoined, readFilesIn2ndJoined, outFileNamePrefix )
    print "Star Command: ", cmd
    #os.system(cmd)



# Run htseq counts
def runHTseq(libraryName, sampleResultsDir):
    print "\033[34m Running HTSeq counts... \033[0m"

    # create results folder
    htseqFolder = '%s/htseq-counts' %(resultsDir)
    if not os.path.exists(htseqFolder):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(htseqFolder)
        os.makedirs(htseqFolder)
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(htseqFolder)

    htseqInputFile = sampleResultsDir+"/"+libraryName+"_star_Aligned.sortedByCoord.out.bam"
    cmd = '/users/sturkars/anaconda2/bin/python -m HTSeq.scripts.count -s "reverse" -t "CDS" -i "locus_tag" -r pos -f bam %s %s > %s/%s_htseqcounts.txt' %(htseqInputFile,genomeGFF,htseqFolder,libraryName)
    print cmd
    os.system(cmd)

################# Main #################

# Input files
dataDir = "/proj/omics4tb/sturkarslan/icalvum_assembly/data"
genomeDir = "/proj/omics4tb/sturkarslan/icalvum_assembly/reference_c5"
resultsDir = "/proj/omics4tb/sturkarslan/icalvum_assembly/results_c5"
fastqcDir = '%s/fastqc' %(resultsDir)
#trimDir = '%s/trimmed' %(resultsDir)
genomeFasta = glob.glob('%s/*_genome.fna' %(genomeDir))
genomeFasta = ' '.join(genomeFasta)
genomeGFF = glob.glob('%s/*_genome.gff' %(genomeDir))
genomeGFF = ' '.join(genomeGFF)

print genomeGFF, genomeFasta,fastqcDir
#sys.exit()


#Create Results directory
if not os.path.exists('%s' %(resultsDir)):
    print '\033[34m %s directory does NOT exists.Creating one! \033[0m' %(resultsDir)
    os.makedirs('%s' %(resultsDir))
else:
    print '\033[34m %s directory exists. Not creating. \033[0m' %(resultsDir)



# Programs
starPath = "/users/sturkars/STAR-2.6.0a/bin/Linux_x86_64/STAR" # path to STAR executable
trimmomaticPath = "/users/sturkars/Trimmomatic-0.35/trimmomatic-0.35.jar"
fastqc = "/users/sturkars/FastQC/fastqc" # path to fastqc executable


# Genome index
genomeIndex(genomeDir,starPath,genomeFasta)
#sys.exit()
firstPair, secondPair, libraryName, sampleResultsDir, folder = runPipeline()
