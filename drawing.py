#!/usr/bin/env python3
#SBATCH --partition=long       ### Partition (like a queue in PBS)
#SBATCH --job-name=RRPS7          ### Job Name
#SBATCH --time=1-20:01:00       ### Wall clock 0ime limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=28     ### Number of tasks to be launched per Node
#SBATCH --mail-user=rarichardson92@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Don't forget to load modules (python/3.6.0) before running code!

####
# Required imports
####

import argparse
import re
import cairo

####
# User input reading
####

def getarguments():
    parser=argparse.ArgumentParser(description = "Reads in an RNA/DNA FASTA file or UCSC coordinates along with desired RNA/DNA motifs and outputs one of the following: visually projects where motifs are in the sequence, provides a tabular list of indexes for motif and exon locations in the sequence on the command line, tsv with output similar to commandline printout. First base index position = 1, reverse positions based on reverse compliment with the same orientation as the initial sequence (From the left of the reverse compliment). Motif matches are insensitive to RNA/DNA file types (i.e., an RNA motif will match to a DNA FASTA and vice-versa). Defaults to FASTA input and an image output. Exons are read as capitalized sequence from FASTA. Motifs should be seperated by newline characters. Use of UCSC options requires download of UCSC 2bit file and twoBitToFa function (if not already present), thus requires internet access. Any pulled UCSC sequence can be found in UCSC.fa file after program completion. Imports argparse, re, and cairo, please install packages appropriately before use.")
    parser.add_argument("-m", "--motiffile", help = "Motif file name and path. Motifs should be seperated by newline characters.", required = True, type = str)
    parser.add_argument("-inF", "--infile", help = "String that determines if a FASTA file or UCSC coordinates are used. Defaults to FASTA. Set to \"UCSC\" for UCSC coordinates.", required = False, type = str, default = "FASTA")
    parser.add_argument("-FASTA", help = "Defines name and path of FASTA file to use in program. Must be a string.", required = False, type = str)
    parser.add_argument("-UCSC", help = "Defines UCSC coordinates use in program. Read as a string in format of database.chr##:#-#. Example: -UCSC hg19.chr21:1000-2000", required = False, type = str)
    parser.add_argument("-o", "--output", help = "String that determines if output is a tsv, on the command line, or an image. Defaults to an image. Options are \"tsv\", \"cmd\", and \"image\".", required = False, type = str, default = "image")
    return parser.parse_args()

args=getarguments()
motiffile=str(args.motiffile)
inF=str(args.infile) #FASTA/UCSC option
FASTA=str(args.FASTA)
UCSC=str(args.UCSC)
outtype=str(args.output)

#Raises Exceptions for incorrect entries
if inF=="FASTA" and FASTA=="None":
    raise Exception('Must input a FASTA value. Set -inF flag to UCSC for UCSC coordinates.')
if inF=="UCSC":
    if UCSC=="None":
        raise Exception('Must input a UCSC value. Set -inF flag to FASTA (default) for FASTA files.')
    elif "-" not in UCSC or ":" not in UCSC or "." not in UCSC:
        raise Exception('Incorrect UCSC input format. Example: -UCSC hg19.chr21:1000-2000')
if inF!="UCSC" and inF!="FASTA":
    raise Exception('inF type must be UCSC or FASTA. See help for setting -inF flag.')
if outtype!="cmd" and outtype!="image" and outtype!="tsv":
    raise Exception('Output type must be \"tsv\", \"cmd\", or \"image\". See help for setting -o flag.')

####
# Handling UCSC input
####

#Calls sequence from UCSC to use as a fasta
if inF=="UCSC":
    import subprocess           #Allows for UCSC calling via bash commands
    from pathlib import Path    #Checks for file existance, no need to redownload if file exists

    #Splitting correctly inputted UCSC coordinates
    UCSC = UCSC.split(".")
    database = UCSC[0]
    pos = UCSC[1].split(":")
    seq = pos[0]
    stend = pos[1].split("-")
    start = stend[0]
    end = stend[1]

    #Creating bash command from user input
    inp = "./twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/"+database+"/"+database+".2bit UCSC.fa -seq="+seq+" -start="+start+" -end="+end

    #Downloads twoBitToFa if its not already existing and makes it executible,  calls the UCSC query and makes a FASTA, then sets as FASTA for rest of program
    if not Path("./twoBitToFa").is_file():
        subprocess.call("wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa", shell=True)
        subprocess.call("chmod a+x twoBittoFa", shell=True)
    subprocess.call(inp, shell = True)
    FASTA = "UCSC.fa"

####
#Functions
####

def motifs(motiffile):
    """Generates a list of motifs, list variable denoted as motiflist"""
    motiflist = []                          #holder for motifs
    with open(motiffile, "rt") as motifs:
        for line in motifs:                 #Each motif is on a seperate line of code
            line=line.strip('\n')               #Strip newline characters
            motiflist.append(line)              #Add to list
    return(motiflist)

def reverse(seq):
    """Generates the reverse complimentary sequence of input sequence. Case sensitive."""
    comp=""                             #holds reverse sequence
    for i in reversed(range(len(seq))): #reverses seq string; then for each letter, adds it compliment to comp variable
        if seq[i] == "A":
            comp=comp+"T"
        elif seq[i] == "T":
            comp=comp+"A"
        elif seq[i] == "C":
            comp=comp+"G"
        elif seq[i] == "G":
            comp=comp+"C"
        elif seq[i] == "N":
            comp=comp+"N"
        elif seq[i] == "a":
            comp=comp+"t"
        elif seq[i] == "t":
            comp=comp+"a"
        elif seq[i] == "c":
            comp=comp+"g"
        elif seq[i] == "g":
            comp=comp+"c"
        elif seq[i] == "n":
            comp=comp+"n"
    return(comp)

def createdict(motiflist):
    """Generates a dictionary for use in finding motifs and exons, dictionary variable denoted as features."""
    features = dict()
    features["header"] = ""     #holds header value for each sequence, initialized as empty
    features["exon"] = ""       #holds pairs of exon start/end indexes, initialized as empty
    features["length"] = ""     #holds length of sequence, initialized as empty
    for motif in motiflist:
        features[motif] = ""        #holds pairs of motif start/end indexes for each motif (forward and reverse), initialized as empty
    return(features)

def exonfind(currentseq):
    """Finds exons in a string, marked in captital letters. Output updates features dictionary. Note: N is not considered in exon search. Exons will be considered to start before and after any encountered N."""
    indexlist = []                              #holder for tuple of exon start/end indexes
    exons = re.finditer(r'[ATCG]+', currentseq) #Picks out indexes of capital bases in file.
    for match in exons:                             #Each match is the location of an exon
        startend = (match.start()+1, match.end())   #Makes a tuple of match start and end index. Add one to start to reflect true start site in sequence (first base = 1)
        indexlist.append(startend)                  #Adds tuple to list
    return(indexlist)

def motifadjust(motiflist):
    """Handles ambiguous motifs, substitutes for string that can be used in search functions. Case insensitive."""
    adjustlist = []             #Holder for adjusted motifs

    #Substitutes any ambiguous notatation for a string thats readable in regex; see "ambiguous bases" on wikipedia for list reference
    for motif in motiflist:     #adjusts each motif, re.sub makes sure to do a global substitute for each ambiguous notation
        adjust = re.sub('N', '[ATCG]', motif, flags=re.IGNORECASE)
        adjust = re.sub('R', '[AG]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('Y', '[TC]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('K', '[GT]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('M', '[AC]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('B', '[CTG]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('D', '[ATG]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('H', '[ACT]', adjust, flags=re.IGNORECASE)
        adjust = re.sub('V', '[ACG]', adjust, flags=re.IGNORECASE)
        adjustlist.append(adjust)
    return(adjustlist)

def motiffind(currentseq, motif):
    """Finds the adjusted motif in the forward and reverse complement of a given sequence. Case insensitive. Uses the reverse function."""
    fwdindexlist = []   #Holder for as-is matches in the input Sequence
    revindexlist = []   #Holder for reverse matches in the input sequence. NOTE: that these are indexed based on reversed input, and therefore the same directionality as the initial input (either 5' - 3' in both or 3' to 5' in both)
    revseq = reverse(currentseq)                                    #Generates reverse compliment of input sequence

    if "U" in currentseq or "u" in currentseq:
        if "T" in motif or "t" in motif:         #Converts for DNA/RNA file mismatches
            motif = re.sub('[Tt]', 'U', motif, flags=re.IGNORECASE)
    elif "U" in motif or "u" in motif:
        if "T" in currentseq or "t" in currentseq:       #Converts for DNA/RNA file mismatches
            motif = re.sub('[Uu]', 'T', motif, flags=re.IGNORECASE)

    fwdmatch = re.finditer(motif, currentseq, flags=re.IGNORECASE)  #finds motifs in original input
    revmatch = re.finditer(motif, revseq, flags=re.IGNORECASE)      #finds motifs in reverse compliment input. See above note.
    for match in fwdmatch:                          #each match contains motif indexes
        startend = (match.start()+1, match.end())   #Add one to start to reflect true start site in sequence
        fwdindexlist.append(startend)               #Adds to list
    for match in revmatch:                          #each match contains motif indexes
        startend = (match.start()+1, match.end())   #Add one to start to reflect true start site in sequence
        revindexlist.append(startend)               #Adds to list
    return(fwdindexlist, revindexlist)

def featureupdate(features, currentseq, motifdict):
    """Updates features dictionary. Uses the following functions: exonfind, motiffind."""
    features["length"] = str(len(currentseq))                    #Holds length of input sequence
    features["exon"] = exonfind(currentseq)                 #Initiates exon search on sequence, holds exon start/end indexes from search output
    for motif in motifdict:                                    #For each motif,
        fwdmotif, revmotif = motiffind(currentseq, motifdict[motif])       #Initiates motif search on sequence, variable holds motif start/end indexes from search output for original input and its reverse complement
        features[motif] = (fwdmotif, revmotif)                  #stores variable in features dictionary
    return(features)

def draw(features, width, height, margin, fontsize):
    """Creates svg image output."""
    length = features.pop("length")                                             #Assigns the length from features dictionary input as a variable. NOTE: Pop is used throughout to reduce dictionary to motifs only by the end of the function call. (Easy for drawing based on each motif in dictionary)
    adjustedwidth = (width - 2*margin)/int(length)                              #For scaling the image based on input sequence length; this varible is multiplied by index coordinates
    header = features.pop("header")                                             #Assigns header variable for printing in image
    title = header                                                              #Seperate varible from header; manipulated if header cannot be used as a file name

    #Makes sure the svg name is able to be saved.
    forbid = ["\\", "/", ":", "*", "?", "\"", "<", ">", "|"]                    #All forbidden characters for file naming
    for item in forbid:
        if item in title:
            title = re.sub(item, '.', title, flags=re.IGNORECASE)               #Replaces all forbidden characters with "."; ensures no overwite occurs from shortening since headers should be unique
    surface = cairo.SVGSurface(title+".motifoutput.svg", width, height)         #opens up svg based on sequence header

    #Print text on image for header, sequence length, and notes about the displayed image
    context = cairo.Context(surface)    #generate surface
    context.set_source_rgb(1.0, 1.0, 1.0)
    context.rectangle(0, 0, width, height) #x initial, y inital (upper), width of box, height of box
    context.fill() ## fill rectangle
    context.set_source_rgb(0.0, 0.0, 0.0)
    context.move_to(margin, fontsize+margin)
    context.set_font_size(fontsize)
    context.show_text("Sequence: "+header+", length: "+length+" bases")
    context.move_to(margin, 3*fontsize+margin)
    context.show_text("Motifs:")
    context.move_to(margin, height-2*fontsize-margin)
    context.show_text("Exons depicted as black boxes. Introns depicted as thin black lines. Motifs are colored lines above and below the sequence.")
    context.move_to(margin, height-1*fontsize-margin)
    context.show_text("Motifs above the sequence are aligned to the input sequence. Motifs below align to the reverse compliment of the input sequence.")

    #Draw line for full sequence size in the middle of the image
    context.move_to(margin, height/2) #initial x, y coordinates
    context.line_to((margin+int(length)*adjustedwidth), height/2) # Line destination equal to full sequence length
    context.stroke() #draws lines

    #For each exon, makes a box over sequence line at scaled coordinates; exon[0] is the exon start coordinate, exon[1] is the end coordinate
    exons = features.pop("exon") #See above note about popping, stores exons as a variable.
    for exon in exons:
        context.rectangle((margin+int(exon[0])*adjustedwidth),240,(int(exon[1])-int(exon[0]))*adjustedwidth,20) #x initial, y inital (upper), width of box, height of box
        context.fill() ## fill rectangle

    #Creates increments based on the number of motifs left in the features dictionary
    colorminus = 2.5/(len(features)+1)      #Subtracts from a light coloration; since color minus takes is divided by 2.5 and not 2.8 (0.9*3, initial set of colorlist), colors should not be too light or dark.
    loop = 0
    colorlist = [0.90, 0.90, 0.90]          #inital color setting

    for motif in features:
        colorlist[loop%3] -= colorminus                                         #adjust red, green, or blue color based on loop and previous setting
        loop += 1
        context.set_source_rgb(colorlist[0], colorlist[1], colorlist[2])        #set color based on colorlist

        #Populates text for motif, same color as motif marks
        context.move_to(margin, (loop+3)*fontsize+margin)
        context.show_text(motif)
        #Assigns variables
        motif = features[motif]
        fwdmotif = motif[0]
        revmotif = motif[1]

        if fwdmotif != ['']:        #Where there are alignments
            for span in fwdmotif:
                context.move_to(margin+(int(span[0])*adjustedwidth), height/2-20-(10*loop)) #initial x, y coordinates (from left)
                context.line_to(margin+(int(span[1])*adjustedwidth), height/2-20-(10*loop)) # Line destination equal to full sequence length
                context.stroke() #draws lines

        if revmotif != ['']:
            for span in revmotif:
                context.move_to(width-(margin+(int(span[0])*adjustedwidth)), height/2+20+(10*loop)) #initial x, y coordinates (from right, revcomp coordinates)
                context.line_to(width-(margin+(int(span[1])*adjustedwidth)), height/2+20+(10*loop)) # Line destination equal to full sequence length
                context.stroke() #draws lines
    surface.finish() # end, closes image

def output(features, out, output):
    """Creates and saves desired output as a svg, tsv, or on the command line."""
    if output == "cmd":                                                         #prints on the cmd line, does length, header, motif into from the features dictionary
        print(features["header"])
        print("Sequence length: ", features["length"])
        print("Exon positions", features["exon"], sep = "\t")
        print('\n')
        print("Motif", "Occurances", "Forward positions", "Reverse positions", sep = "\t")
        for item in features:
            if item != "exon" and item != "header" and item != "length":        #motifs only
                print(item, len(features[item][0])+len(features[item][1]), features[item][0], features[item][1], sep = "\t")
        print("\n")
        print("\n")
        out.write('Selected output: cmd. Output information can be found on the command line.') #tsv message informing output location

    elif output == "tsv":                                                       #writes to tsv, does length, header, motif into from the features dictionary
        out.write(features["header"])
        out.write('\n')
        out.write("Sequence length: %s" % str(features["length"]))
        out.write('\n')
        out.write("Exon positions\t%s" % str(features["exon"]))
        out.write('\n')
        out.write('\n')
        out.write("Motif    Occurances  Forward positions   Reverse positions")
        out.write('\n')
        for item in features:
            if item != "exon" and item != "header" and item != "length":        #motifs only
                out.write("%s\t%s\t%s\t%s" % (item, str(len(features[item][0])+len(features[item][1])), features[item][0], features[item][1]))
                out.write('\n')
        out.write('\n')
        out.write('\n')

    elif output == "image":
        out.write('Selected output: image. Output information can be found in motif.svg file.')  #tsv message informing output location
        width, height, margin, fontsize = 1000, 500, 20, 13   #image parameters
        draw(features, width, height, margin, fontsize)       #draws image for sequence

def readFASTA(FASTA, features, motifdict, outtype):
    """Finds motifs and exons in fasta files. Outputs a dictionary with features corresponding to each header.  Uses the featureupdate function."""
    with open(FASTA, "rt") as FASTA, open("Motifoutput.tsv", "w") as out:
        currentseq = ""                 #Holder for sequence for the current header, String initialized as empty
        for line in FASTA:
            if ">" in line:                 #Finds header lines
                if currentseq != "":            #When its not the first header, process prior header dictionary info
                    features = featureupdate(features, currentseq, motifdict)   #Updates features dictionary
                    output(features, out, outtype)                              #Creates output
                    currentseq = ""                                             #resets current sequence holder
                line=line.strip('\n')
                features["header"]=line.strip('>') #Sets dictionary header; stripped of newline and >
            else:
                currentseq+=line.strip('\n')       #Adds sequence line to current seq; stripped of newline
        features = featureupdate(features, currentseq, motifdict)               #Updates features dictionary
        output(features, out, outtype)                                          #Creates output

####
# Main code
####

motiflist=motifs(motiffile)                         #Makes motif list from motifs file
adjustlist=motifadjust(motiflist)                    #Adjusts motifs for ambiguous motifs

motifdict={}                                          #Creates dict with original motif as keys and adjusted motif as values
i = 0
for motif in motiflist:
    motifdict[motif]= adjustlist[i]
    i+=1
features=createdict(motiflist)                      #Creates features dictionary
readFASTA(FASTA, features, motifdict, outtype)      #Reads fasta, creating output for each header based on updating the features dictionary over the lines of fasta
