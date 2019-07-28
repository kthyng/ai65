'''
Plot results from calcs
'''

import matplotlib.pyplot as plt
import pandas as pd

basesrc = 'savedoutput/'
bc = pd.read_csv(basesrc + 'power_baroclinic_anomaly.csv', index_col=0, parse_dates=True)
bt = pd.read_csv(basesrc + 'power_barotropic_anomaly.csv', index_col=0, parse_dates=True)
o = pd.read_csv(basesrc + 'power_otherterms_correct.csv', index_col=0, parse_dates=True)
z = pd.read_csv(basesrc + 'tidal.csv', index_col=0, parse_dates=True)

# combine other terms from o into bc and bt
bc = pd.concat([bc,o[['Pmix','Pmom','Pfriction']],z['zeta']], axis=1)
bt = pd.concat([bt,o[['Pmix','Pmom','Pfriction']],z['zeta']], axis=1)

# # take absolute value of P
# bc['P'] = abs(bc['P'])
# bt['P'] = abs(bt['P'])

# add total diss terms
bc['terms'] = bc['Pmix'] + bc['Pmom'] + bc['Pfriction']

# plots in time
bc[['Pin','Pout','Pmix','Pmom','Pfriction']].plot(figsize=(12,4))
bc[['P','Pmix','Pmom','Pfriction']].plot(figsize=(12,4))
bc[['P','terms']].plot(figsize=(12,4))

sdate = '2006-09-01 18:00'
edate = '2006-09-30 17:00'
fig = plt.figure()
ax = fig.add_subplot(211)
bc[sdate:edate]['zeta'].plot(ax=ax)
# ax = fig.add_subplot(412)
# (bc[sdate:edate][['P','Pin','Pout']]/1e6).plot(ax=ax)
# ax = fig.add_subplot(413)
# # (bc[:date][['PEwest','PEsouth', 'KEsouth']]/1e6).plot(ax=ax)
# (bc[sdate:edate]['P']/1e6-bc[sdate:edate]['terms']/1e6).plot(ax=ax)
ax = fig.add_subplot(212)
(bc[sdate:edate][['Pmix','Pmom','Pfriction']]/1e6).plot(ax=ax)
# (bc[:'2006-09-02']['P']/1e6-bc[:'2006-09-02']['terms']/1e6).plot(ax=ax)
plt.show()


fig = plt.figure()
ax = fig.add_subplot(311)
bt[:'2006-09-02']['zeta'].plot(ax=ax)
ax = fig.add_subplot(312)
(bt[:'2006-09-02'][['KEwest']]/1e6).plot(ax=ax)
ax = fig.add_subplot(313)
(bt[:'2006-09-02'][['PEwest']]/1e6).plot(ax=ax)
plt.show()
