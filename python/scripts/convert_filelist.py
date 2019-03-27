with open('filelist.csv') as input:
    for line in input:
        sline = line.strip().split(',')
        if sline[0] == '"SampleName"':
            print('"SampleName","SLXId","Barcode","Filename","Flowcell","Lane","SequencingPlatform","PlatformModel","SequencingCentre","SequencingDate","EndType","ReadLength"')
        else:
            ssource = sline[1].split('.')
            slxid = ssource[0][1:]
            barcode = ssource[1]
            print('{},"{}","{}",{},{},{},{},{},{},{},{},{}'.format(sline[0], slxid, barcode, sline[2], sline[3], sline[4], sline[5], sline[6], sline[7], sline[8], sline[9], sline[10]))

# "SampleName","SourceID",       "Filename","Flowcell","Lane","SequencingPlatform","PlatformModel","SequencingCentre","SequencingDate","EndType","ReadLength"
# "SampleName","SLXId","Barcode","Filename","Flowcell","Lane","SequencingPlatform","PlatformModel","SequencingCentre","SequencingDate","EndType","ReadLength"
