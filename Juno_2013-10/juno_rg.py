"""
http://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/pck/pck00010.tpc
http://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/spk_merge_110805_171017_130515.bsp
http://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/lsk/naif0011.tls
\begindata

KERNELS_TO_LOAD = (
  '$K/pck00010.tpc'
  '$K/spk_merge_110805_171017_130515.bsp'
  '$K/naif0011.tls'
)

PATH_VALUES = ( '.' )
PATH_SYMBOLS = ( 'K' )

ETNEARPRG = @2013-10-09T00:00:00

SPACECRAFT = 'JUNO'
TARGET = 'EARTH'

\begintext
"""

import SpiceyPy as spice

dpr = spice.dpr()

spice.furnsh('.'.join(__file__.split('.')[:-1]+['py']))

target = spice.gcpool('TARGET',0,1,99)[0]
spacecraft = spice.gcpool('SPACECRAFT',0,1,99)[0]
etnearca = spice.gdpool('ETNEARPRG',0,1)[0]

targetID = spice.bods2c(target)

et = etnearca
state,lt = spice.spkezr(spacecraft,et,'J2000','NONE',target)

dt = spice.vdot(state[:3],state[3:]) / spice.vdot(state[3:],state[3:])
oldet = et
et = et - dt
while abs(et-oldet)>0.:
  state,lt = spice.spkezr(spacecraft,et,'J2000','NONE',target)
  dt = spice.vdot(state[:3],state[3:]) / spice.vdot(state[3:],state[3:])
  oldet = et
  et = et - dt

utcPRG = spice.et2utc(et,'ISOC',3,99)

mtxJ2kToBF = spice.tipbod('J2000', targetID, et)

velocityPRG_BF = spice.mxv(mtxJ2kToBF,state[3:])
V_prg,ra,Theta = spice.recrad(velocityPRG_BF)

Theta_deg = Theta * dpr

print((dt,dt/et,utcPRG,V_prg,Theta_deg))

import matplotlib.pyplot as plt

days = 7
hours = xrange(-24*days,24*days + 1)

speeds = [ spice.vnorm(spice.spkezr(spacecraft,et+(hour*3600.),'J2000','NONE',target)[0][3:]) for hour in hours]

Vmin = min(speeds)
plt.axhline(y=V_prg,label='V_prg = %.3fkm/s' % (V_prg,))
plt.axhline(y=Vmin,label='Vmin = %.3fkm/s' % (Vmin,))
plt.plot(hours,speeds,'r',label='%s-relative speed' % (target,))

plt.title('Perigee at %s (UTC); Theta=%.3fdeg' % (utcPRG,Theta_deg,))
plt.xlabel('t - T_prg, h')
plt.ylabel('%s-relative speed, km/s' % (target,))

plt.legend()

plt.show()
