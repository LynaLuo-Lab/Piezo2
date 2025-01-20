import sys
import math
import numpy as np
from scipy.optimize import leastsq

def groGxyz(line):
    _x=float(rl[20:28])*10
    _y=float(rl[28:36])*10
    _z=float(rl[36:44])*10
    _xyz=[_x,_y,_z]
    return _xyz
def pdbGxyz(line):
    _x=float(rl[30:38])
    _y=float(rl[38:46])
    _z=float(rl[46:54])
    _xyz=[_x,_y,_z]
    return _xyz

# Input Parameters
pdb=sys.argv[1]   ### input PDB from OPM database
_length=float(sys.argv[2])   ### lipid length, width
_width=float(sys.argv[3])



_height       = 50
_mp_posit_ot  = "outer"
_mp_posit_ul  = "upper"
_pdb_data     = []
_xyz_data     = []
_ddo_xyz_data = []
_ddi_xyz_data = []

for rl in open(pdb,'r'):
    if len(rl)>=50:
        if rl[:4]=="ATOM" or rl[:6]=="HETATM":
            _pdb_data.append(rl)
            _xyz=pdbGxyz(rl)
            _xyz_data.append(_xyz)
            if rl[:6]=="HETATM" and rl[17:20]=="DUM":
                if rl[13:14]=="O":
                    _ddo_xyz_data.append(_xyz)
                elif rl[13:14]=="N":
                    _ddi_xyz_data.append(_xyz)

### according to the membrane protein position, to get the SSMP postion
_dd_xyz_data=[]
if _mp_posit_ot=="outer":
    _dd_xyz_data=_ddo_xyz_data
elif _mp_posit_ot=="inner":
    _dd_xyz_data=_ddi_xyz_data
#print(_dd_xyz_data)
_np_dd_xyz_data=np.array(_dd_xyz_data)
_dd_min=_np_dd_xyz_data.min(axis=0)
_dd_max=_np_dd_xyz_data.max(axis=0)
_box_z=40+_height*2+_height*2+(_dd_max[2]-_dd_min[2])*2
if _mp_posit_ul=="upper":
    _delta_z=40+_height+(_height*2+(_dd_max[2]-_dd_min[2]))+(0-_dd_min[2])
elif _mp_posit_ul=="lower":
    _delta_z=(_dd_max[2]-_dd_min[2])+_height+(0-_dd_max[2])
#print(_delta_z,_box_z)
_dd_x_col=int(max(_dd_max[0],-(_dd_min[0]))/4)
print("simulation box:",_length,_width,int(_box_z))
# print(_delta_z*0.1)

'''
### write a new pdb file
fwp=open(pdb[:4]+"_box.pdb",'w')
pdbline="%s%8.3f%8.3f%8.3f%s"
for rl in _pdb_data:
    _x=float(rl[30:38])+(_length/2)
    _y=float(rl[38:46])+(_width/2)
    _z=float(rl[46:54])+_delta_z
    rl=pdbline%(rl[:30],_x,_y,_z,rl[54:])
    fwp.write(rl)
fwp.close()
'''

### generate SSMB pdb file 
def yFit(y):
    _y_max=max(y)
    _y_min=min(y)
    _y_li=int(max(_y_max,-(_y_min))/4)
    _y_fit=[]
    for i in range(-(_y_li),_y_li+1):
        _y_fit.append(float(4*i))
    return _y_fit
def residuals(p):
    a,b,r=p
    return r**2-(_fz-b)**2-(_fy-a)**2
fwd=open("dome.pdb",'w')
boxline="CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1 \n"
fwd.write(boxline%(_length,_width,int(_box_z)))
pdbline_1="ATOM  %5d  ON  DUM A%4d    %8.3f%8.3f%8.3f  1.00  1.00 \n"
num=0
_dd_fit_xyz_data=[]
for i in range(-(_dd_x_col),_dd_x_col+1):
    _x_posit=float(4*i)
    _fy=[]
    _y_fit=[]
    _fz=[]
    for j in range(len(_dd_xyz_data)):
        if _x_posit==_dd_xyz_data[j][0]:
            _fy.append(_dd_xyz_data[j][1])      ## y z values to fit function
            _fz.append(_dd_xyz_data[j][2])
    _y_fit=yFit(_fy)
    result=leastsq(residuals,[1,1,1])
    a,b,r=result[0]
    if _delta_z >= 0:
        _z_fit=-((r**2-(_y_fit-a)**2)**0.5+b)   
    else:
        _z_fit=((r**2-(_y_fit-a)**2)**0.5+b)   
    for k in range(len(_y_fit)):
        num=num+1
        _line=pdbline_1%(num,num,_x_posit+(_length/2),_y_fit[k]+(_width/2),_z_fit[k]+_delta_z)
        _dd_fit_xyz_data.append([_x_posit+(_length/2),_y_fit[k]+(_width/2),_z_fit[k]+_delta_z])
        fwd.write(_line)
        #fwd.write(pdbline_1%(num,num,_x_posit,_y_fit[k],_z_fit[k]))
pdbline_2="ATOM  %5d  DT  DD1 A%4d    %8.3f%8.3f%8.3f  1.00  1.00 \n"
pdbline_4="ATOM  %5d  DT  DD2 A%4d    %8.3f%8.3f%8.3f  1.00  1.00 \n"
if _delta_z >= 0:
    #fwd.write(pdbline_2%(num+1,num+1,_length/2,_width/2,_height*2+(_dd_max[2]-_dd_min[2])*0.5))
    fwd.write(pdbline_2%(num+1,num+1,_length/2,_width/2,_box_z-_height*2.5 - 1.5*(_dd_max[2]-_dd_min[2])))
    fwd.write(pdbline_4%(num+1,num+1,_length/2,_width/2,_box_z-_height*2-(_dd_max[2]-_dd_min[2])*0.5))
else:
    fwd.write(pdbline_2%(num+1,num+1,_length/2,_width/2,_box_z-_height))

    
pdbline_3="ATOM  %5d  DO  DON A%4d    %8.3f%8.3f%8.3f  1.00  1.00 \n"
center=160
sp=4
cn=5
_r_z=sp*cn
radius=[]
for i in range(-(cn),cn+1):
    radius.append(center+i*sp)
#print(radius)
_dno_num=0
for r in radius:    
    _r_x = float(center-abs(r))
    z1 =  (_r_z**2 - _r_x**2)**0.5
    z2 = -(_r_z**2 - _r_x**2)**0.5
    for angle in range(0,360,4):
        #print(angle)
        _x_d = r * math.cos(angle * 3.14 / 180)
        _y_d = r * math.sin(angle *3.14 /180)
        if z1==z2:
            num=num+1
            _dno_num=_dno_num+1
            fwd.write(pdbline_3%(num+1,num+1,_x_d + _length*0.5, _y_d + _width*0.5,z1+(_box_z*0.5-20-_height-_r_z)))
        else:
            fwd.write(pdbline_3%(num+1,num+1,_x_d + _length*0.5, _y_d + _width*0.5,z1+(_box_z*0.5-20-_height-_r_z)))
            num=num+1
            _dno_num=_dno_num+1
            fwd.write(pdbline_3%(num+1,num+1,_x_d + _length*0.5, _y_d + _width*0.5,z2+(_box_z*0.5-20-_height-_r_z)))
            num=num+1
            _dno_num=_dno_num+1
fwd.close()


### generate topology file
fwt=open("ssmb.itp",'w')
fwt.write("[ moleculetype ] \n")
fwt.write("; Name         Exclusions \n")
fwt.write("Protein_DUM   1 \n")
fwt.write("\n")

### atoms
fwt.write("[ atoms ] \n")
atomline="%4d    DD %5d   DUM    ON %5d  0.0000 ; O \n"
for i in range(len(_dd_fit_xyz_data)):
    fwt.write(atomline%(i+1,i+1,i+1))
fwt.write("\n")

### bonds
_cutoff_min=0.5   #float(sys.argv[2])
_cutoff_max=3.0   #float(sys.argv[3])
def Bond(p1,p2):
    dist=(((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)**0.5)*0.1
    return dist
bondline="%5d %5d 6 %9.6f 1000 \n"
fwt.write("[ bonds ] \n")
for i in range(len(_dd_fit_xyz_data)):
    _p1=_dd_fit_xyz_data[i]
    for j in range(i+1,len(_dd_fit_xyz_data)):
        _p2=_dd_fit_xyz_data[j]
        dist=Bond(_p1,_p2)
        if dist >= _cutoff_min and dist <= _cutoff_max:
            fwt.write(bondline%(i+1,j+1,dist))
fwt.write("\n")

### position restraints
fwt.write("[ position_restraints ] \n")
#   1    1  POSRES_FC    POSRES_FC    POSRES_FC
posreline="%5s    1   1000  1000  1000 \n"
for i in range(len(_dd_fit_xyz_data)):
    fwt.write(posreline%(i+1))
fwt.close



### write system.top
fwst=open("system.top",'w')
fwst.write("#include \"ff/martini_v2.2.itp\" \n")
fwst.write("#include \"ff/martini_v2.0_lipids_all_201506.itp\" \n")
fwst.write("#include \"ff/martini_v2.0_ions.itp\" \n")
fwst.write("#include \"ssmb.itp\" \n \n")
fwst.write("#include \"ff/martini_v2.2_ddd.itp\" \n")
fwst.write("#include \"don.itp\" \n \n")

fwst.write("[ system ] \n")
fwst.write("; name \n")
fwst.write("sphere shaped membranes \n \n")

fwst.write("[ molecules ] \n")
fwst.write("; name  number \n")
fwst.write("Protein_DUM        1  \n")
fwst.write("DD1        1  \n")
fwst.write("DD2        1  \n")
fwst.write("Protein_DON        1  \n")


