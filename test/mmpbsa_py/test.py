import subprocess

OUTPUT_NAME = "TEST_OUTPUT.dat"

def change_inp_igb(igb):
    new_lines=[]
    f = open('mmpbsa.in')
    for l in f:
        if l.__contains__('igb='):
            new_lines.append('   igb=' + igb + ',\n')
        else:
            new_lines.append(l)
    f.close()
    f = open('mmpbsa.in', 'w')
    for l in new_lines:
        f.write(l)
    f.close()

for i, igb in enumerate(['2', '5', '7', '8', '66']):
    if i > 0:
        subprocess.run(['rm', OUTPUT_NAME])
    change_inp_igb(igb)
    subprocess.run(['./Run.sh', igb, OUTPUT_NAME])