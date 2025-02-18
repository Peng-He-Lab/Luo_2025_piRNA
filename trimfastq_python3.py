#!/usr/bin/env python3
import sys
import os

def main(argv):
    if len(argv) < 3:
        print("usage: python {} <inputfilename> <bpToKeep | max> [-trim5 bp] [-flowcellID flowcell] [-addEnd 1 | 2] [-replace string newstring | blank] [-renameIDs prefix] [-stdout]".format(argv[0]))
        print("\tthe -trim5 option will trim additional bp from the 5 end, i.e. if you want the middle 36bp of 38bp reads, use 36 as bp to keep and 1 as the trim5 argument")
        print("\tUse - to specify standard input, and the -stdout option to tell the script to print to standard output")
        print("\tThe script can read compressed files as long as they have the correct suffix - .bz2 or .gz")
        print("\tReplace inputfilename with - if you want to read from standard input")
        sys.exit(1)

    inputfilename = argv[1]
    doMax = False
    if argv[2] == 'max':
        doMax = True
        trim = 'max'
    else:
        try:
            trim = int(argv[2])
        except ValueError:
            print("Error: bpToKeep must be an integer or 'max'")
            sys.exit(1)

    outputfilename = inputfilename.split('/')[-1].split('.fastq')[0] + '.' + str(trim) + 'mers.fastq'
    doFlowcellID = False

    doStdOut = False
    if '-stdout' in argv:
        doStdOut = True

    if '-flowcellID' in argv:
        doFlowcellID = True
        flowcellID = argv[argv.index('-flowcellID')+1]
        if not doStdOut:
            print("will include flowcell ID", flowcellID, "in reads headers")

    doRenameIDs = False
    if '-renameIDs' in argv:
        doRenameIDs = True
        RID = '@' + argv[argv.index('-renameIDs') + 1]

    dotrim5 = False
    if '-trim5' in argv:
        dotrim5 = True
        try:
            trim5 = int(argv[argv.index('-trim5')+1])
        except ValueError:
            print("Error: -trim5 argument must be an integer")
            sys.exit(1)
        if not doStdOut:
            print("will trim", trim5, "bp from the 5'-end")
        outputfilename = inputfilename.split('.fastq')[0] + '.' + str(trim) + 'bp-5prim-trim.fastq'

    doAddEnd = False
    if '-addEnd' in argv:
        doAddEnd = True
        END = argv[argv.index('-addEnd')+1]
        if not doStdOut:
            print("will add", '/' + END, "to read IDs")

    doReplace = False
    if '-replace' in argv:
        doReplace = True
        oldstring = argv[argv.index('-replace')+1]
        newstring = argv[argv.index('-replace')+2]
        if newstring == 'blank':
            newstring = ''
        if not doStdOut:
            print("will replace", oldstring, "with", newstring, "in read IDs")

    i = 0 
    shorter = 0

    if not doStdOut:
        outfile = open(outputfilename, 'w')

    doStdIn = False
    if inputfilename != '-':
        if inputfilename.endswith('.bz2'):
            cmd = 'bzip2 -cd ' + inputfilename
        elif inputfilename.endswith('.gz'):
            cmd = 'gunzip -c ' + inputfilename
        else:
            cmd = 'cat ' + inputfilename
        p = os.popen(cmd, "r")
    else:
        doStdIn = True

    line = 'line'

    if dotrim5:
        i = 1
        j = 0
        while True:
            if doStdIn:
                line = sys.stdin.readline()
            else:
                line = p.readline()
            if line == '':
                break
            
            if i == 1 and line[0] == '@':
                if doFlowcellID and flowcellID not in line:
                    ID = '@' + flowcellID + '_' + line.replace(' ', '_')[1:-1] + '\n'
                else:
                    ID = line.replace(' ', '_')
                if doReplace:
                    ID = ID.replace(oldstring, newstring)
                if doRenameIDs:
                    ID = RID + str(j)
                if doAddEnd:
                    ID = ID.strip() + '/' + END + '\n'
                i = 2
                continue
            
            if i == 2:
                i = 3
                sequence = line[trim5:].strip()
                continue
            
            if i == 3 and line[0] == '+':
                plus = '+\n'
                i = 4
                continue
            
            if i == 4:
                i = 1
                scores = line[trim5:].strip()
                scores = scores[0:trim]
                j += 1
                if j % 5000000 == 0:
                    if not doStdOut:
                        print(str(j // 1000000) + 'M reads processed')
                if doMax:
                    sequence = sequence.replace('.', 'N')
                else:
                    sequence = sequence[0:trim].replace('.', 'N') + '\n'
                if doStdOut:
                    print(ID.strip())
                    print(sequence.strip())
                    print(plus.strip())
                    print(scores)
                else:
                    outfile.write(ID.strip() + '\n')
                    outfile.write(sequence.strip() + '\n')
                    outfile.write(plus.strip() + '\n')
                    outfile.write(scores + '\n')
                continue
    else:
        i = 1
        j = 0
        while True:
            if doStdIn:
                line = sys.stdin.readline()
            else:
                line = p.readline()
            if line == '':
                break
            if i == 1 and line[0] == '@':
                if doFlowcellID and flowcellID not in line:
                    ID = '@' + flowcellID + '_' + line.replace(' ', '_')[1:-1] + '\n'
                else:
                    ID = line.replace(' ', '_')
                if doReplace:
                    ID = ID.replace(oldstring, newstring)
                if doRenameIDs:
                    ID = RID + str(j)
                if doAddEnd:
                    ID = ID.strip() + '/' + END + '\n'
                i = 2
                continue
            if i == 2:
                i = 3
                j += 1
                if j % 5000000 == 0:
                    if not doStdOut:
                        print(str(j // 1000000) + 'M reads processed')
                if doMax:
                    sequence = line
                else:
                    if len(line.strip()) < trim:
                        shorter += 1
                        sequence = line.strip().replace('.', 'N') + '\n'
                    else:
                        sequence = line[0:trim].replace('.', 'N') + '\n'
                continue
            if i == 3 and line[0] == '+':
                plus = '+\n'
                i = 4
                continue
            if i == 4:
                i = 1
                if doMax:
                    scores = line
                    if doStdOut:
                        print(ID.strip())
                        print(sequence.strip())
                        print(plus.strip())
                        print(line.strip())
                    else:
                        outfile.write(ID)
                        outfile.write(sequence)
                        outfile.write(plus)
                        outfile.write(line)
                else:
                    if len(line.strip()) < trim:
                        continue
                    scores = line[0:trim] + '\n'
                    if doStdOut:
                        print(ID.strip())
                        print(sequence.strip())
                        print(plus.strip())
                        print(scores.strip())
                    else:
                        outfile.write(ID)
                        outfile.write(sequence)
                        outfile.write(plus)
                        outfile.write(scores)
                continue

    if not doStdOut:
        outfile.close()

    if shorter > 0:
        print(str(shorter) + ' sequences shorter than desired length')

if __name__ == '__main__':
    main(sys.argv)

