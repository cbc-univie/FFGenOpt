import sys

try:
    cgen = sys.argv[1]
    dgen = sys.argv[2]
except IndexError:
    print("Usage: python psf_merge.py nonpolarizable.psf polarizable.psf")
    sys.exit()

f1 = open(cgen, 'r')
f2 = open(dgen, 'r')
f3 = open(cgen[:-4]+'_pol.psf', 'w')

num_count = 0

while num_count < 2:
    x1 = f1.readline()
    x2 = f2.readline()
    f3.write(x1)
    if len(x1) > 1 and x1.split()[0].isdigit():
        num_count += 1
num_atom = int(x1.split()[0])
for i in range(num_atom):
    x1 = f1.readline()
    x2 = f2.readline()
    while x1.split()[4] != x2.split()[4]:
        x2 = f2.readline()
    old_atom = x1.split()[5]
    new_atom = x2.split()[5]
    len_diff = len(old_atom)-len(new_atom)
    if len_diff < 0:
        for i in range(-len_diff):
            old_atom += " "
    elif len_diff > 0:
        for i in range(len_diff):
            new_atom += " "
    x3 = x1.replace(old_atom, new_atom)
    f3.write(x3)
while len(x1) != 0:
    x1 = f1.readline()
    f3.write(x1)


f1.close()
f2.close()
f3.close()