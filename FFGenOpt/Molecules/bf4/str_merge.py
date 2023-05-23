import sys

try:
    cgen = sys.argv[1]
    dgen = sys.argv[2]
except IndexError:
    print("Usage: python str_merge.py nonpolarizable.str polarizable.str")
    sys.exit()

f1 = open(cgen, 'r')
f2 = open(dgen, 'r')
f3 = open(cgen[:-4]+'_conv.str', 'w')

atom_types = {}

x1 = f1.readline()
x2 = f2.readline()
f3.write(x1)

while not x1.startswith("GROUP"):
    x1 = f1.readline()
    f3.write(x1)

while not x2.startswith("GROUP"):
    x2 = f2.readline()

x1 = f1.readline()
x2 = f2.readline()
while len(x1.split()) != 0:
    if x1.split()[0] == "ATOM" and x2.split()[0] == "ATOM" and \
    x1.split()[1] == x2.split()[1]:
        old_atom = x1.split()[2]
        new_atom = x2.split()[2]
        atom_types[old_atom] = new_atom
    x3 = x1.replace(old_atom, new_atom)
    f3.write(x3)
    x1 = f1.readline()
    x2 = f2.readline()

f3.write(x1)

while len(x1) != 0:
    x1 = f1.readline()
    x3 = x1
    for i in x1.split():
        if i in atom_types:
            x3 = x3.replace(i, atom_types[i])
    f3.write(x3)

f1.close()
f2.close()
f3.close()
