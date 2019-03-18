def msd(dti=1,save=1,traj=1):
    '''
    Mean squared displacement
    '''

    print('msd')

    if not traj:
        conf = sim.loadConfiguration('config_big.npy')
        poss = np.array(conf[:,0,2:5])
    else:
        conf = np.load(sim.path / 'pt_0.npy')
        poss = conf
    
    poss = poss[:-1]

    dt_i = 0
    dt_f = 50000

    pn = len(poss)

    ti_i = len(poss)-100000
    ti_f = len(poss)-dt_f

    '''
    sc = np.zeros(dt_f)
    c = 0
    for ti in range(ti_i,ti_f,dti):
        c += 1
        for dt in range(dt_i,dt_f):
            dx = poss[ti,0]-poss[ti+dt,0]
            dy = poss[ti,1]-poss[ti+dt,1]
            dz = poss[ti,2]-poss[ti+dt,2]
            sc[dt] += dx*dx+dy*dy+dz*dz
    sc = sc/(c)
    '''
    sc = np.zeros(dt_f)
    for dt in range(dt_i,dt_f):
        dr = poss[0:pn-dt]-poss[dt:pn]
        sc[dt] = np.mean(dr*dr)
    
    if save:
        np.savetxt(sim.path / 'msd_{}'.format(dti),sc)
        np.save(sim.path / 'msd_{}'.format(dti),sc)
        
    return sc