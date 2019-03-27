with open('samplesheet.csv') as input:
    for line in input:
        sline = line.strip().split(',')
        if sline[0] == '"SampleName"':
            print('"SampleName","FSFSampleName"')
        else:
            print('{},{}'.format(sline[0], sline[0]))
