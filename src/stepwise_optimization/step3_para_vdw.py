##エネルギー等計算の対象を確認　Rt:t-shaped Rp:slipped-parallel
import numpy as np
import pandas as pd
from utils import Rod

def get_monomer_xyzR1(monomer_name,Ta,Tb,Tc,A2,A3):## rotation and translation operations
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('~/path/to/monomer/{}.csv'.format(monomer_name))
    df_mono=df_mono[df_mono['Z']>0]
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)## A2 is a torsion angle in step2
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)## A3 is a half of dihedral angle
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    return np.concatenate([xyz_array,R_array],axis=1)
    
def get_monomer_xyzR2(monomer_name,Ta,Tb,Tc,A2,A3):## rotation and translation operations
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('~/path/to/monomer/{}.csv'.format(monomer_name))
    df_mono=df_mono[df_mono['Z']<0]
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)## A2 is a torsion angle in step2
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)## A3 is a half of dihedral angle
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    return np.concatenate([xyz_array,R_array],axis=1)
    
def get_c_vec_vdw(monomer_name,Rt,Rp,a_,b_,theta):#,name_csv
    
    i=np.zeros(3); a=np.array([a_,0,2*Rt-Rp]); b=np.array([0,b_,Rp]); t1=(a+b)/2;t2=(a-b)/2 ##ずらし方の定義
    monomer_array_i0 = get_monomer_xyzR1(monomer_name,0.,0.,0.,0.,theta);monomer_array_t0 = get_monomer_xyzR1(monomer_name,0.,0.,0.,0.,-theta)
    monomer_array_i = get_monomer_xyzR2(monomer_name,0.,0.,0.,0.,theta);monomer_array_t = get_monomer_xyzR2(monomer_name,0.,0.,0.,0.,-theta)
    
    arr_list1=[[i,'p'],[b,'p'],[-b,'p'],[a,'p'],[-a,'p'],[t1,'t'],[-t1,'t'],[t2,'t'],[-t2,'t']]##molecular layer with 9 atoms
    arr_list2=[[i,'t'],[b,'t'],[-b,'t'],[a,'t'],[-a,'t'],[t1,'p'],[-t1,'p'],[t2,'p'],[-t2,'p']]##molecular layer with 9 atoms type2
    cy_list=[np.round(cy,1) for cy in np.linspace(-np.round(b_/2,1),np.round(b_/2,1),int(np.round(2*np.round(b_/2,1)/0.1))+1)]
    cx_list=[np.round(cx,1) for cx in np.linspace(-np.round(a_/2,1),np.round(a_/2,1),int(np.round(2*np.round(a_/2,1)/0.1))+1)]
    para_list=[]
    for cx in cx_list:
        for cy in cy_list:
            shift = (2*Rt-Rp)*cx/a_ + Rp*cy/b_

            z_max1=0
            for R,arr in arr_list1:
                vx,vy,vz=R
                if arr=='t':
                    monomer_array1=monomer_array_t0
                elif arr=='p':
                    monomer_array1=monomer_array_i0
                x1 = monomer_array1[:,0];y1 = monomer_array1[:,1];z1 = monomer_array1[:,2];R1 = monomer_array1[:,3]
                x1=x1[:,None]+vx;y1=y1[:,None]+vy;z1=z1[:,None]+vz
                x2 = monomer_array_i[:,0];y2 = monomer_array_i[:,1];z2 = monomer_array_i[:,2];R2 = monomer_array_i[:,3]
                dx = x1[:,None] - (x2[None,:] + cx)
                dy = y1[:,None] - (y2[None,:] + cy)

                z_sq = (R1[:, None] + R2[None, :])**2 - dx**2 - dy**2

                z_clps = np.where(z_sq > 0,np.sqrt(z_sq) + z1[:, None] - (z2[None, :] + shift),0.0)
                z_max1=max(z_max1,np.max(z_clps))
            z_max2=0
            for R,arr in arr_list2:
                vx,vy,vz=R
                if arr=='t':
                    monomer_array1=monomer_array_t0
                elif arr=='p':
                    monomer_array1=monomer_array_i0
                x1 = monomer_array1[:,0];y1 = monomer_array1[:,1];z1 = monomer_array1[:,2];R1 = monomer_array1[:,3]
                x1=x1[:,None]+vx;y1=y1[:,None]+vy;z1=z1[:,None]+vz
                x2 = monomer_array_t[:,0];y2 = monomer_array_t[:,1];z2 = monomer_array_t[:,2];R2 = monomer_array_t[:,3]
                dx = x1[:,None] - (x2[None,:] + cx)
                dy = y1[:,None] - (y2[None,:] + cy)

                z_sq = (R1[:, None] + R2[None, :])**2 - dx**2 - dy**2

                z_clps = np.where(z_sq > 0,np.sqrt(z_sq) + z1[:, None] - (z2[None, :] + shift),0.0)
                z_max1=max(z_max1,np.max(z_clps))
            z_max=max(z_max1,z_max2)
            cz=round(z_max+(2*Rt-Rp)*cx/a_+Rp*cy/b_,1)
            para_list.append([cx,cy,cz,z_max*a_*b_])
    df=pd.DataFrame(para_list,columns=['cx','cy','cz','V'])
    return df

