import os

def ReadTfToM(input_file):
'''
Reads a file where each line contains a transcription factor and its associated motifs. Returns a dictionary mapping each transcription factor to a list of its motifs.
'''
    with open(input_file, 'r') as f:

        dic = {}

        for line in f.readlines():

            line = line.split()
            tf = line[0]
            m = line[1].split()
            dic[tf] = m

    return dic

def CheckQvalue(pwms_dir):
'''
Concatenates all .meme files in a directory and runs Tomtom to compare motifs. If any comparison has a q-value above 0.001, it returns False. Otherwise, it returns the name of the longest PWM file.
'''
    os.system('cat ' + pwms_dir + '*.meme > multipwm.meme')
    os.system('tomtom multipwm.meme multipwm.meme -thresh 1.0 -verbosity 1')
    dic = {}

    with open('./tomtom_out/tomtom.tsv', 'r') as tomtom_output:
        header = tomtom_output.readline()
        line = tomtom_output.readline()

        while line != '\n':

            line = line.split()
            qvalue = float(line[5])

            if qvalue > 0.001: return False

            line = tomtom_output.readline()

    longest = ''
    l = 0
    for pwm in os.listdir('pwms'):
        f = open('pwms/'+pwm,'r')
        n = len(f.readlines())
        if l < n: longest = pwm
    return longest[:-5]

def CheckMotifs(input_file, output_file, pwms_dir, dummy_dir):
'''
For each transcription factor and its associated motifs, checks motif similarity using Tomtom. If motifs are similar (based on q-value), writes the transcription factor and selected motif to an output file. Uses a temporary directory to run the checks.
'''
    tfs = ReadTfToM(input_file)

    if not os.exists(dummy_dir): os.mkdir(dummy_dir)

    with open(output_file, 'w') as out_file:

        for tf in tfs:
            
            os.system('rm ' + dummy_dir + '*')

            motifs = [pwms_dir + m + '.meme' for m in tfs[tf]]
            os.system('cp ' + ' '.join(motifs) + ' ' + dummy_dir)

            if len(motifs) > 1:
                check = check_qvalue()
                if check:
                    out_file.write(f"{tf} {check}\n")

            else:
                out_file.write(f"{tf} {' '.join(tfs[tf])}\n")
