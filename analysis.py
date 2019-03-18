#!/usr/bin/env python
import numpy as np
import os

lx = 30.0
ly = 30.0
lz = 30.0

def loadConfiguration(filename):
    f = open(filename,'r')
    c = 0
    for line in f:
        c += 1
    conf = np.zeros((c,18))
    f = open(filename,'r')
    for i, line in enumerate(f):
        if i > c-2: break
        conf[i] = np.array(line.split())
    return conf

def loadTrajectory(filename,save=1,load=0):
    l = np.array([lx,ly,lz])
    rg = np.array([0.0,0.0,0.0])
    if not load:
        conf = loadConfiguration(filename)
        poss = conf[:,0:3]
        traj = np.zeros((len(poss)-1,3))
        for i, r in enumerate(poss[1:]):
            d = poss[i]-poss[i-1]
            s = np.sign(d)
            for j in range(3):
                if np.abs(d[j]) > l[j]/2.0:
                    rg[j] -= l[j]*s[j];
            traj[i] = poss[i]+rg
    else:
        traj = np.load('traj_big.npy')

    if save:
        np.savetxt('traj_big',traj)
        np.save('traj_big',traj)

    return traj

def msd(dti=1,save=1,load=1):
    '''
    Mean squared displacement
    '''

    print('msd')

    l = np.array([lx,ly,lz])
    rg = np.array([0.0,0.0,0.0])

    traj = loadTrajectory('config_big.txt',save=1,load=0)

    dt_i = 0
    dt_f = 1000000
    ddt = 100

    pn = len(traj)

    ti_i = 0

    sc = np.zeros(int((dt_f-dt_i)/ddt))
    for i, dt in enumerate(range(dt_i,dt_f,ddt)):
        print(dt)
        dr = traj[ti_i:pn-dt]-traj[ti_i+dt:pn]
        sc[i] = np.mean(dr*dr)
    
    if save:
        np.savetxt('msd_{}'.format(dti),sc)
        np.save('msd_{}'.format(dti),sc)
    
    print(sc)    

    return sc


def vsc(angular=0,dti=1,save=1):
    '''
    Velocity self correlation
    '''

    print('vsc')

    conf = loadConfiguration('config_big.txt')
    vels = conf[:,3:6]

    dt_i = 0
    dt_f = 30000
    ddt = 1

    ti_i = 0
    ti_f = len(vels)

    sc = np.zeros(dt_f)
    c = 0
    for dt in range(dt_i,dt_f,ddt):
        print(dt)
        sc[dt] += np.dot(vels[ti_i:ti_f-dt,0],vels[ti_i+dt:ti_f,0])
        sc[dt] += np.dot(vels[ti_i:ti_f-dt,1],vels[ti_i+dt:ti_f,1])
        sc[dt] += np.dot(vels[ti_i:ti_f-dt,2],vels[ti_i+dt:ti_f,2])
    sc /= (ti_f-ti_i)
    
    if save:
        if angular:
            np.savetxt('wsc_{}'.format(dti),sc)
        else:
            np.savetxt('vsc_{}'.format(dti),sc)

    return sc

#msd(save=1,load=1)
vsc()