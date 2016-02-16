import numpy as np
import matplotlib.pyplot as plt

# import data
obs_mags = np.loadtxt(fname='SDSSPS1_de.txt', dtype=float)
syn_mags = np.loadtxt(fname='fkcor_SDSS_5_15_2_ab.csv', dtype=float, delimiter=',')

# plot boundaries
xmin = -1.0
xmax = 3
ymin = -0.4
ymax = 0.3

# grab data for observed magnitudes
ps1g = np.copy(obs_mags[:, 0])
ps1r = np.copy(obs_mags[:, 1])
ps1i = np.copy(obs_mags[:, 2])
ps1z = np.copy(obs_mags[:, 3])
sdssg = np.copy(obs_mags[:, 4])
sdssr = np.copy(obs_mags[:, 5])
sdssi = np.copy(obs_mags[:, 6])
sdssz = np.copy(obs_mags[:, 7])

# grab data for synthetic magnitudes
syn_ps1g = np.copy(syn_mags[:, 6])
syn_ps1r = np.copy(syn_mags[:, 7])
syn_ps1i = np.copy(syn_mags[:, 8])
syn_ps1z = np.copy(syn_mags[:, 9])
syn_sdssg = np.copy(syn_mags[:, 2])
syn_sdssr = np.copy(syn_mags[:, 3])
syn_sdssi = np.copy(syn_mags[:, 4])
syn_sdssz = np.copy(syn_mags[:, 5])


def plot(obslist1, obslist2, synlist1, synlist2, name1, name2):
    # set up plots
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    obslist1_obslist2 = []
    synlist1_synlist2 = []
    ps1g_ps1i = []
    syn_ps1g_ps1i = []

    for i in range(len(obslist1)):
        if abs(obslist1[i]-obslist2[i]) > 10.0:
            pass
        else:
            obslist1_obslist2.append(obslist1[i]-obslist2[i])
            ps1g_ps1i.append(ps1g[i]-ps1i[i])

    for i in range(len(synlist1)):
        if abs(synlist1[i]-synlist2[i]) > 10.0:
            pass
        else:
            synlist1_synlist2.append(synlist1[i]-synlist2[i])
            syn_ps1g_ps1i.append(syn_ps1g[i]-syn_ps1i[i])

    # calculate best fit line for both observed and synthetic magnitude
    obs_linfit = np.polyfit(x=ps1g_ps1i, y=obslist1_obslist2, deg=1)
    syn_linfit = np.polyfit(x=syn_ps1g_ps1i, y=synlist1_synlist2, deg=1)
    x_list = np.linspace(xmin, xmax, 10)
    obsy_list = []
    syny_list = []
    for f in x_list:
        obsy_list.append(obs_linfit[0] * f + obs_linfit[1])
        syny_list.append(syn_linfit[0] * f + syn_linfit[1])

    # plot ps1g-sdssg
    plt.plot(x_list, obsy_list, 'blue', linestyle='-', linewidth=2, label='Observed fit')
    plt.plot(x_list, syny_list, 'm', linestyle='-', linewidth=2, label='Synthetic fit')
    plt.scatter(x=ps1g_ps1i, y=obslist1_obslist2, c='black', marker='o', label='Observed')
    plt.scatter(x=syn_ps1g_ps1i, y=synlist1_synlist2, c='red', marker='x', label='Synthetic')
    # make sure label is not covering a bunch of data points
    if obs_linfit[0] <= 0:
        plt.legend(loc='upper right')
    else:
        plt.legend(loc='lower right')
    plt.ylabel('%s-%s (mag)' % (name1, name2))
    plt.xlabel('PS1g-PS1i (mag)')
    plt.show()

    # histogram and bins
    count = xmin
    interval = 0.5
    obsdiff_list = []
    syndiff_list = []
    obsbinned_list = []
    synbinned_list = []

    # create bins of y-offset data per interval
    while count < xmax:
        obsbin = []
        synbin = []
        for k in range(len(ps1g_ps1i)):
            if count <= ps1g_ps1i[k] < count + interval:
                y = obslist1_obslist2[k]
                y0 = obs_linfit[0] * ps1g_ps1i[k] + obs_linfit[1]
                obsdiff_list.append(y-y0)
                obsbin.append(obslist1_obslist2[k])
        for l in range(len(syn_ps1g_ps1i)):
            if count <= syn_ps1g_ps1i[l] < count + interval:
                z = synlist1_synlist2[l]
                z0 = syn_linfit[0] * syn_ps1g_ps1i[l] + syn_linfit[1]
                syndiff_list.append(z-z0)
                synbin.append(synlist1_synlist2[l])

        obsbinned_list.append(obsbin)
        synbinned_list.append(synbin)
        obsbin = []
        synbin = []
        count += interval
    # bins = len(np.arange(-np.max(obslist1_obslist2), np.max(obslist1_obslist2), 0.01))

    # create histogram of all differences
    # plt.hist(obsdiff_list, 40, histtype='step')
    plt.hist(syndiff_list, 40, color='green', histtype='step')
    plt.xlabel('%s-%s (mag)' % (name1, name2))
    plt.show()
    return [obsbinned_list, synbinned_list]

plot(ps1g, sdssg, syn_ps1g, syn_sdssg, 'ps1g', 'sdssg')
plot(ps1i, sdssi, syn_ps1i, syn_sdssi, 'ps1i', 'sdssi')
plot(ps1r, sdssr, syn_ps1r, syn_sdssr, 'ps1r', 'sdssr')
plot(ps1z, sdssz, syn_ps1z, syn_sdssz, 'ps1z', 'sdssz')


