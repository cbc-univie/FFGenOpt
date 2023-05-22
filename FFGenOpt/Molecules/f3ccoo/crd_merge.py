import sys

try:
    cgen = sys.argv[1]
    dgen = sys.argv[2]
except IndexError:
    print("Usage: python crd_merge.py nonpolarizable.crd polarizable.crd")
    sys.exit()

f1 = open(cgen, 'r')
f2 = open(dgen, 'r')
f3 = open(cgen[:-4]+'_conv.crd', 'w')

x1 = f1.readline()
x2 = f2.readline()
f3.write(x1)

while not x1.split()[0][0].isdigit():
    x1 = f1.readline()
    x2 = f2.readline()
    f3.write(x1)

num_atoms = int(x1.split()[0])
for i in range(num_atoms):
    x1 = f1.readline()
    x2 = f2.readline()
    x3 = x1[0:46]+x2[46:101]+x1[101:]
    f3.write(x3)

f1.close()
f2.close()
f3.close()
