#! usr/bin/python

from sys import argv
import math

script, spikealignment, spikestats = argv
# when running the script, you have to enter the filenames
# for the bowtie-output (.sam file: spikealignment), and
# for the outputfile to which you want to write the stats

def spikeAnalysis(alignmentfile, output):
    ''' Function that analyses the spike-ins in your sequencing
    reads by reading through the .sam file that bowtie(2) created.
    
    Input: .sam alignment file with sequences aligned to
    spike-ins AND outputfile, the file to which you want to
    write the output of the analysis (statistics of spike-ins
    in your sample.
    
    Output: writes a file (txt) with all required
    statistics of the spike-ins in your sample.
    Statistics include: spike lengths, number of correct/wrong
    sequences, number of mismatches and Phred score.
    
    Return: list with all different spike-in sequences
    in your sample, that can be used to filter your samples.
    '''
    
    with open(alignmentfile, 'r') as samfile:
            
        spikecounter = {}
        difseqs = {}
        bases = 'ACGT'
    
        for line in samfile:
            if line[0] == '@':
                if line[1:3] == 'SQ':
                    line = line.rstrip('\n').split('\t')
                    spikename = line[1][3:]
                    if len(spikename) == 10:
                        spikename = spikename[:-1] + '0' + spikename[-1:]
                    spikelength = int(line[2][3:])
                    spikecounter[spikename] = [spikelength, 0, 0, 999, 0, -999, 0, 0, 0, 0]
                    # this makes a dictionary with spikenames as
                    # keys and a list of [length, count correct,
                    # count wrong, minimal mismatches, average
                    # mismatches, maximal mismatches, seq too short
                    # , seq right length, seq too long, avg seqlen] as value
            else:
                # if line does not start with '@' -> non-header
                line = line.rstrip('\n').split('\t')
                # remove the '\n's at end of lines and
                # split by tabs -> separate different elements
                flag = int(line[1])
                spikename = line[2]
                if len(spikename) == 10:
                    spikename = spikename[:-1] + '0' + spikename[-1:]
                #convert spikenumbers in 2-digit numbers (i.e. 1 -> 01)
                #this helps with sorting
                sequence = line[9]
                cigar = line[5]
                mismatches = ''
                
                # flags are 0, 4 or 16 -> mapped on +, unmapped, mapped on -
                if flag in [0, 16]:
                #check if sequences have been mapped
                    if line[17][:5] == "MD:Z:":
                        mdz = line[17]
                    elif line[18][:5] == "MD:Z:":
                        mdz = line[18]
                    
                    if sequence in difseqs:
                    #check if sequences have been recorded yet
                        difseqs[sequence] += 1
                    else:
                        difseqs[sequence] = 1
                    
                    spikecounter[spikename][9] += len(sequence)
                    #keep track of the sequence lengths, so that
                    #you can calculate an average length
                    
                    if len(sequence) < spikecounter[spikename][0]:
                    #if the seq is shorter than it should be, count:
                        spikecounter[spikename][6] += 1
                    elif len(sequence) > spikecounter[spikename][0]:
                    #if it is too long, count:
                        spikecounter[spikename][8] += 1
                    else:
                    #if the length is correct:
                        spikecounter[spikename][7] += 1
                        
                    if cigar[:-1] == mdz[5:]:
                    #check if sequence is correct: no mismatches
                        spikecounter[spikename][1] += 1
                        # count correct spikes
                    else:
                    # if they do have mismatches:
                        for character in mdz[5:]:
                            if character in bases:
                                mismatches += character
                    # count the number of mismatches
                        spikecounter[spikename][3] = min(spikecounter[spikename][3], len(mismatches))
                        spikecounter[spikename][4] += len(mismatches)
                        spikecounter[spikename][5] = max(spikecounter[spikename][5], len(mismatches))
                        spikecounter[spikename][2] += 1
                    # and write some statistics about those mismatches


        spikecountersorter = list(spikecounter.keys())
        spikecountersorter.sort()
        
        #make a sorted list of spikenames, to print them 
        
        for spikename in spikecountersorter:
            spikecounter[spikename][4] /= float(spikecounter[spikename][1] + spikecounter[spikename][2])
            spikecounter[spikename][4] = (spikecounter[spikename][4] * 100)
            #calculate the average number of mismatches, rather than keeping the sum
            # -> as a percentage!
        
        header = "Spike\tLength\tAverageLengthSequenced\t#TooShort\t#RightLength\t#TooLong"
        header2 = "\t#Correct\t#Incorrect\tPercentageCorrect"
        header3 = "\tMinimum#Mismatches\tAverageMismatches%\tMaximum#Mismatches\tMismatchesPerBase\tPhredScore"
        
        # separate headers for different statistics

### HERE COMES WRITING THE OUTPUT FILE, COMMENT THIS OUT IF
### YOU DO NOT NEED IT
        with open(output, 'w') as spikestats:
            spikestats.write("%s%s%s\n" % (header, header2, header3))
            
            for key in spikecountersorter:
                Length = spikecounter[key][0]
                TotalLength = spikecounter[key][9]
                TooShort = spikecounter[key][6]
                RightLength = spikecounter[key][7]
                TooLong = spikecounter[key][8]
                CorrectNumber = spikecounter[key][1]
                IncorrectNumber = spikecounter[key][2]
                percentageCorrect = float(CorrectNumber) / (float(CorrectNumber) + float(IncorrectNumber)) * 100
                MinimumMismatches = spikecounter[key][3]
                AverageMismatches = spikecounter[key][4]
                MaximumMismatches = spikecounter[key][5]
                mismatchPerLength = AverageMismatches / (CorrectNumber + IncorrectNumber)
                mismatchPerBase = int(Length) / mismatchPerLength
                AverageLength = float(TotalLength) / (CorrectNumber + IncorrectNumber)
                Phred = -10 * math.log((1 / mismatchPerBase), 10)
                
                spikestats.write("%s\t%s\t%2.2f\t%i\t%i\t%i\t%i\t%i\t%2.2f\t%i\t%2.2f\t%i\t1 in %i\t%i\n" % (key, Length, AverageLength, TooShort, RightLength, TooLong, CorrectNumber, IncorrectNumber, percentageCorrect, MinimumMismatches, AverageMismatches, MaximumMismatches, mismatchPerBase, Phred))

### END OF FILE WRITING
### NOW COMES THE PRINTING, AGAIN COMMENT OUT IF YOU DO NOT
### NEED THIS
        
        #print (header, header2, header3)
#        for key in spikecountersorter:
#        
#            Length = spikecounter[key][0]
#            TotalLength = spikecounter[key][9]
#            TooShort = spikecounter[key][6]
#            RightLength = spikecounter[key][7]
#            TooLong = spikecounter[key][8]
#            CorrectNumber = spikecounter[key][1]
#            IncorrectNumber = spikecounter[key][2]
#            percentageCorrect = float(CorrectNumber) / (float(CorrectNumber) + float(IncorrectNumber)) * 100
#            MinimumMismatches = spikecounter[key][3]
#            AverageMismatches = spikecounter[key][4]
#            MaximumMismatches = spikecounter[key][5]
#            mismatchPerLength = AverageMismatches / (CorrectNumber + IncorrectNumber)
#            mismatchPerBase = int(Length) / mismatchPerLength
#            AverageLength = float(TotalLength) / (CorrectNumber + IncorrectNumber)
#            Phred = -10 * math.log((1 / mismatchPerBase), 10)
            
            #print (key, "\t%s\t%2.2f\t%i\t%i\t%i\t%i\t%i\t%2.2f\t%i\t%2.2f\t%i\t1 in %i\t%i" % (Length, AverageLength, TooShort, RightLength, TooLong, CorrectNumber, IncorrectNumber, percentageCorrect, MinimumMismatches, AverageMismatches, MaximumMismatches, mismatchPerBase, Phred))

        print ("\nThe statistics have now been written to the file you named!")

### END OF PRINT

        seqslist = difseqs.keys()
        return seqslist
                
if __name__ == "__main__":
    spikeAnalysis(spikealignment, spikestats)
