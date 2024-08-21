import happi
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["mathtext.fontset"]='cm'
fig,ax = plt.subplots(4,4,sharey='col',sharex='col',squeeze=1,figsize=[13,5.5])
for simName, axCol,time_idx,label in zip(['2Dpinch13','2Dpinch17','2Dpinch10','2Dpinch11'],
                                         [ax[:,3],ax[:,2],ax[:,1],ax[:,0]],
                                         [[0,4,10,40],[0,8,30,90],[0,4,80,150],[0,50,100,199]],
                                         ['(d)','(c)','(b)','(a)']):
    S = happi.Open(simName)
    Jz = S.Field(0,'Jz_ions+Jz_electrons')
    tsteps = Jz.getTimesteps()
    times  = Jz.getTimes()
    xaxis = Jz.getAxis('x')
    yaxis = Jz.getAxis('y')
    n0    = S.namelist.mass_ratio
    for axes,t in zip(axCol,time_idx):
        c = axes.imshow(Jz.getData(tsteps[t])[0]/n0,aspect='equal',
                        extent = [yaxis[0],yaxis[-1],xaxis[0],xaxis[-1]],
                        origin = 'lower',
                        cmap = 'Reds',
                        vmin = 0,
                        )
        ###### Calculate flux functions and plot their contours
        By = S.Field(0,'By')
        Bx = S.Field(0,'Bx')
        dx = S.namelist.dx
        dy = S.namelist.dy
        Nx = S.namelist.xmax
        Ny = S.namelist.ymax
        Az = -np.cumsum(By.getData(tsteps[t])[0],axis=0)*dx
        f  = np.cumsum(Bx.getData(tsteps[t])[0][0],axis=0)*dy
        f  = np.tile(f,(int(Nx+1),1))
        Az += f 
        axes.contour(yaxis,xaxis,Az,levels=10,
                     colors = 'k',
                     linewidths=0.5,
                     linestyles='solid')
        #######
        # axes.set_ylim([xaxis[-1]/4,xaxis[-1]*3/4])
        plt.colorbar(c,ax=axes)
        axes.text(0.01,1.05,'$t\omega_{pi}=%.2f$'%times[t],transform=axes.transAxes,c='k')
    axCol[0].text(0.01,1.3,label,transform=axCol[0].transAxes)
[ax.set_ylabel(r'$x/d_i$') for ax in ax[:,0]]
[ax.set_xlabel(r'$y/d_i$') for ax in ax[-1,:]]
plt.show()
# plt.savefig('/Users/young/Dropbox/APCTP/Research/2023/2023_ReconnectionOnset/manuscript/Jz_evolution.pdf',dpi=300,bbox_inches='tight')
