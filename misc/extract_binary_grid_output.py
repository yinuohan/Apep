'''
Extracts outputs from binary_grid. 
'''

filepath = '/import/silo4/snert/TINTAGEL_DATA/NACO/2019_03/Analysis/Apep/binary_grid_output.txt'

files = []
params = []
errors = []

lines = []
clean = []

with open(filepath) as f:  
    line = f.readline()
    i = 1
    while line:
        #print("Line {}: {}".format(i, line.strip()))
        line = f.readline()
        
        '''
        if 'oifits file used' in line:
            files.append(line)
        elif 'Final best solution' in line:
            params.append(line)
        elif 'Errors: ' in line:
            errors.append(line)
        elif ' VS ' in line or ' HOLES' in line:
            params.append('\n' + line + '\n')
        '''
        
        if 'Final best solution' in line:
            clean.append(line)
            lines.append(line)
        elif 'Errors: ' in line or 'Significance:' in line or 'Contrast (mags):' in line or 'low contrast ratio' in line:
            lines.append(line)
        elif ' HOLES' in line or ' POS' in line:
            lines.append('\n' + line + '\n')
            clean.append('\n' + line)
    
        i += 1
       
savepath = '/import/silo4/snert/TINTAGEL_DATA/NACO/2019_03/Analysis/Apep/output_extract.txt'
savepath_clean = '/import/silo4/snert/TINTAGEL_DATA/NACO/2019_03/Analysis/Apep/output_extract_clean.txt'

savef = open(savepath, 'w')
savef.writelines('                        Sep(mas) PA(degs) Contrast\n\n')
savef.writelines(lines)
savef.close()

savef = open(savepath_clean, 'w')
savef.writelines('                        Sep(mas) PA(degs) Contrast\n\n')
savef.writelines(clean)
savef.close()


