#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Take in quality score string. Return average quality score of the entire input string."""
    sum = 0
    for i in phred_score:
        sum+=convert_phred(i)
    return sum/len(phred_score)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    #pass
    assert validate_base_seq(DNA), "String contains invalid characters"  #this is the error message that will be shown
    seq = DNA.upper()                        #ensure string is all uppercase
    GC = seq.count('G') + seq.count('C')     #count number of Gs + Cs
    return GC/len(seq)                       #calculate GC content of entire sequence

def oneline_fasta(input_file:str, output_file:str) -> str:
    '''Take in string parameters of the input file name where sequences have '\n' in the middle, and the output file name. 
        Function removes '\n' from the middle of sequences in the input file. Output file will have headers on one line, 
        entire sequence on one line. Will return the string 'Successful' when the function is complete'''
    with open(input_file) as input:
        with open(output_file, 'w') as output:
            counter = 0 
            for line in input:
                if ">" in line and counter!=0:
                    output.write('\n')
                if ">" not in line:
                    line=line.strip('\n')

                output.write(line)
                counter+=1
    return('Successful')         

def init_list(lst: list, length: int, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it as a list of length 101, 
    each position storing the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    
    counter = 0
    while counter < length:
        if(value==[]):
            lst.append([])
        else:
            lst.append(value)

        counter+=1
    return lst

def populate_list(file: str) -> tuple[list, int]:
    """Take in file name as a string. Return a list of converted Phred quality score its corresponding number and 
    added to an ongoing sum inside your list of the quality scores for each base position. Also return an int of the number
    of lines in your file."""
    
    #call init_list to create an empty list
    qual_score: list = []
    
    qual_score = init_list(qual_score)
    
    #open FASTQ file
    with open(file, "r") as opened: 
    
    #loop through every record
    #Keep a counter called to keep track of the total number of lines in the file.
        counter = 0
        while True:
            line = opened.readline()
            counter+=1
            if(line == ""):
                counter-=1
                break
            elif("@HWI" in line):
                line2 = opened.readline()
                line3 = opened.readline()
                line4 = opened.readline()
                counter+=3
                line4 = line4.strip() #get rid of \n
            
    #convert the Phred quality score from a letter to its corresponding number 
    #and add it to an ongoing sum inside your list of the quality scores for each base pair.
                counter_phred = 0
                for i in line4:
                    phred_pos = bioinfo.convert_phred(i)
                    qual_score[counter_phred] += phred_pos
                    counter_phred+=1
    return(qual_score,counter)

def rev_comp(seq: str) -> str: 
    """Take in output file and """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} #the complements for each DNA pair, N has no complement so it stays as N
    forward = list(seq)
    reverse = ""
    for i in reversed(forward):
        reverse += complement[i]
    return reverse  

def write_output(file, array, ind1: str, ind2: str) -> str:
    """Take in output file and """
    #add the indexes to the end of the header line
    #ind2_rev = rev_comp(ind2)
    new_header = array[0] + ' ' + ind1 + '_' + ind2
    #write the entire record (updated header) to desired output file
    file.write(new_header+ '\n' + array[1]+ '\n' + array[2]+ '\n' + array[3]+ '\n')
    return("Success")

if __name__ == "__main__":
    # write tests for functions above
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Convert_phred function is working!")
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calculated the correct average phred score")
    assert validate_base_seq("AATAGAT") == True, "DNA string not recognized"
    assert validate_base_seq("AGCTACTGNNNNCTACTG") == True, "DNA string not recognized"
    print("Correctly identified DNA")
    assert validate_base_seq("Coding is fun") == False, "Non-DNA identified as DNA"
    assert validate_base_seq("This week was exhausting") == False, "Non-DNA identified as DNA"
    print("Correctly determined non-DNA")
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")
    #I will include this bioinfo_test.txt file with my submission
    #testing = oneline_fasta('bioinfo_test.txt', 'bioinfo_result.txt')
    #assert(testing == 'Successful')
    print('Converted FASTA correctly')
    assert rev_comp('AATGCTA') == "TAGCATT"
    assert rev_comp('ATCATGCG') == "CGCATGAT"
    print('Correctly found reverse complement')
        #my bioinfo_test.txt file:
        # >seq1
        # A
        # C
        # T
        # G
        # >seq2
        # AAAAAAA
        # >seq3
        # ATCTCA
        # GTACAG
        # TGACST
        # >seq4
        # shfgosk;l
