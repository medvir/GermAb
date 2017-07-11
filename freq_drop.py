#!/usr/local/bin/python3.4
'''
freq_drop.py
Input is a list of sorted counts (descending).
Determines where the biggest relative change between two consecutive elements occurs.
Change has to be at least four-fold.
Returns the number of elements before this drop.
Returns 1 if only one element is in the list.
'''

import sys

def freq_drop(reads):
    n = len(reads)
    
    if n == 1:
        d = 1
    else:
        for i in range(1, n):
            if float(reads[i-1]) >= 3 * float(reads[i]):
                d = i
                break
            else:
                d = n
                
    return(d)
    
    
    #else:
    #    changes = []
    #    for i in range(1, n):
    #        change = float(reads[i]) / float(reads[i-1])
    #        changes.append(change)
    #        #print >> sys.stderr, changes
    #        d = changes.index(min(changes)) + 1
    #       if reads[d-1] < 3 * reads[d]:
    #           d = n
    
    

def main():
    reads = sys.argv[1:]
    reads = [int(r) for r in reads] 
    d = freq_drop(reads)
    print(d)

if __name__ == '__main__':
    main()