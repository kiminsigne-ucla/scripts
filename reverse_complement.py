def reverse_complement(sequence):
    """Return the reverse complement of the supplied sequence string """ 
    
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', '.':'.', '*':'*'} 
    #Reverse
    sequence = sequence[::-1]
    #Complement
    letters = list(sequence) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters) 