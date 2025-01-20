import sys

_atom=[]
_posit=[]
for rl in open(sys.argv[1],'r'):
    _atom.append(int(rl[6:11]))
    _x=float(rl[30:38])
    _y=float(rl[38:46])
    _z=float(rl[46:54])
    _posit.append([_x,_y,_z])

print('[ moleculetype ]')
print('; Name         Exclusions')
print('Protein_DON   1')
print('\n')


print("[ atoms ]")
#    1    Qd     1   LYS    BB     1  1.0000 ; C
atomline="%4d    DD %5d   DON    DO %5d  0.0000 ; O"
for i in range(len(_atom)):
    print(atomline%(i+1,i+1,i+1))
print('\n')


_cutoff_min=float(sys.argv[2])
_cutoff_max=float(sys.argv[3])
def EN(p1,p2):
    dist=(((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)**0.5)*0.1
    return dist
#enline="%5d %5d 1 %9.6f RUBBER_FC*1.000000"
enline="%5d %5d 6 %9.6f 1000"
print("[ bonds ]")
for i in range(len(_posit)):
    p1=_posit[i]
    for j in range(i+1,len(_posit)):
        p2=_posit[j]
        dist=EN(p1,p2)
       # print(dist)
        if dist >= _cutoff_min and dist <= _cutoff_max:
            print(enline%(i+1,j+1,dist) )
print('\n')

print("[ position_restraints ]")
#   1    1  POSRES_FC    POSRES_FC    POSRES_FC
posreline="%5s    1   1000  1000  1000"
for i in range(len(_atom)):
    print(posreline%(i+1))
print('\n')
