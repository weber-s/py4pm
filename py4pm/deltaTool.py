import os
import pandas as pd
import numpy as np
import sqlite3
import glob
from zipfile import ZipFile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import scipy.stats as st
import seaborn as sns
import collections
from matplotlib.colors import TABLEAU_COLORS
import matplotlib.lines as mlines
import matplotlib.text as mtext
import itertools

from pyutilities.chemutilities import *

class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, orig_handle, **self.text_props)
        handlebox.add_artist(title)
        return title

# def sourcesColor(s):
#     colors = {
#         'Aged seasalt': "#00b0ff",
#         'Aged seasalt/HFO': "#8c564b",
#         'Biomass burning': "#92d050",
#         'Dust': "#dac6a2",
#         'Fresh seasalt': "#00b0f0",
#         'HFO': "#70564b",
#         'HFO (stainless)': "#8c564b",
#         'Industries': "#7030a0",
#         'Mineral dust': "#dac6a2",
#         'Nitrate rich': "#ff7f2a",
#         'Primary biogenic': "#ffc000",
#         'Resuspended dust': "#dac6a2",
#         'Road salt': "#00b0f0",
#         'Road traffic': "#000000",
#         'Seasalt': "#00b0f0",
#         'Secondary biogenic': "#8c564b",
#         'Secondary biogenic/HFO': "#8c564b",
#         'Sulfate rich': "#ff2a2a",
#         'Sulfate rich/HFO': "#ff2a2a",
#         'Anthropogenic SOA': "#8c564b",
#     }
#     if s not in colors.keys():
#         print("WARNING: no {} found in colors".format(s))
#         return "#666666"
#     return colors[s]

def to_relativeMass(df, totalVar="PM10"):
    """
    Normalize the profile df to the relative mass with regard to the totalVar
    (=PM10 or PM2.5).
    """
    if totalVar not in df.index:
        # print("WARNING: totalVar {} not in index.".format(totalVar))
        totalVar = df[df.index.str.contains("PM")].index[0]
        # print("Taking {}".format(totalVar))

    dftmp = df.copy()
    if isinstance(dftmp, pd.DataFrame):
        dftmp.dropna(axis=1, how="all", inplace=True)
    dftmp /= dftmp.loc[totalVar]
    dftmp.drop(totalVar, inplace=True)
    return dftmp

def compute_SID(df1, df2, source1, source2=None, isRelativeMass=True):
    """ 
    Compute the SID of the sources `source1` and `source2` in the profile
    `df1`and `df2`. 
    """
    if not(source2):
        source2 = source1
    if source1 not in df1.dropna(axis=1, how="all").columns:
        return np.nan
    if source2 not in df2.dropna(axis=1, how="all").columns:
        return np.nan
    if not isRelativeMass:
        p1 = to_relativeMass(df1.loc[:, source1])
        p2 = to_relativeMass(df2.loc[:, source2])
    else:
        p1 = df1.loc[:, source1]
        p2 = df2.loc[:, source2]
    sp = p1.index.intersection(p2.index)
    if len(sp)>3:
        ID = pd.np.abs(p1[sp]-p2[sp])
        MAD = p1[sp] + p2[sp]
        SID = pd.np.sqrt(2)/len(sp) * (ID/MAD).sum()
    else:
        SID = pd.np.nan
    return SID

def compute_PD(df1, df2, source1, source2=None, isRelativeMass=True):
    """
    Compute the PD of the sources `source1` and `source2` in the profile
    `df1`and `df2`. 
    """
    if not(source2):
        source2 = source1
    if source1 not in df1.dropna(axis=1, how="all").columns:
        return np.nan
    if source2 not in df2.dropna(axis=1, how="all").columns:
        return np.nan
    if not isRelativeMass:
        p1 = to_relativeMass(df1.loc[:, source1])
        p2 = to_relativeMass(df2.loc[:, source2])
    else:
        p1 = df1.loc[:, source1]
        p2 = df2.loc[:, source2]
    p1.dropna(inplace=True)
    p2.dropna(inplace=True)
    sp = p1.index.intersection(p2.index)
    if len(sp)>3:
        PD = 1 - pd.np.corrcoef(p1[sp], p2[sp])[1,0]**2
    else:
        PD = pd.np.nan
    return PD

def plot_deltatool_pretty(ax):
    """
    Format the given ax to conform with the "deltatool-like" visualization.
    """
    rect = mpatch.Rectangle((0, 0),
                            width=1, height=0.4,
                            facecolor="green",
                            alpha=0.1, zorder=-1)
    ax.add_patch(rect)
    ax.set_xlim(0, 1.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("SID")
    ax.set_ylabel("PD")

def plot_similarityplot(PMF_profile, station1, station2, source1, source2=None,
                        SID=None, PD=None, isRelativeMass=False, ax=None, plot_kw={}):
    """
    Plot the distance in the SID/PD space of 2 profiles for 2 stations.
    """
    if not source2:
        source2 = source1

    for station, source in itertools.product([station1, station2], [source1, source2]):
        if PMF_profile.loc[station, source].isnull().all():
            print("There is no {} in {}".format(source, station))
            return

    newFig = False
    if not ax:
        f = plt.figure()
        newFig = True
        ax = plt.gca()

    if (SID is None) and (PD is None):
        SID = compute_SID(
            PMF_profile.loc[station1, :], PMF_profile.loc[station2, :],
            source1, source2,
            isRelativeMass=isRelativeMass
        )
        PD = compute_PD(
            PMF_profile.loc[station1, :], PMF_profile.loc[station2, :],
            source1, source2,
            isRelativeMass=isRelativeMass
        )

    ax.plot(SID, PD, "o", **plot_kw)
    
    if newFig:
        plot_deltatool_pretty(ax=ax)

# def plot_similarityplot_all(PMF_profile, stationRef, stations, factors):
#
#     similarity = pd.DataFrame(columns=["x", "y", "xerr", "yerr", "n"],
#                               index=factors)
#     
#     for p in factors:
#         SIDtmp = pd.Series(index=stations)
#         PDtmp = pd.Series(index=stations)
#         for station in stations:
#             if station == stationRef:
#                 continue
#             SIDtmp.loc[station] = compute_SID(
#                 PMF_profile.loc[stationRef],
#                 PMF_profile.loc[station],
#                 source1=p
#             )
#             PDtmp.loc[station] = compute_PD(
#                 PMF_profile.loc[stationRef],
#                 PMF_profile.loc[station],
#                 source1=p
#             )
#         SIDtmp.dropna(inplace=True)
#         PDtmp.dropna(inplace=True)
#         similarity.loc[p, "x"] = SIDtmp.mean()
#         similarity.loc[p, "xerr"] = SIDtmp.mean() \
#                 - st.t.interval(0.95, len(SIDtmp)-1,
#                                 loc=np.mean(SIDtmp),
#                                 scale=st.sem(SIDtmp)
#                               )[0]
#         similarity.loc[p, "y"] = PDtmp.mean()
#         similarity.loc[p, "yerr"] = PDtmp.mean() \
#                         - st.t.interval(0.95, len(PDtmp)-1,
#                                 loc=np.mean(PDtmp),
#                                 scale=st.sem(PDtmp)
#                                )[0]
#         similarity.loc[p, "n"] = SIDtmp.notnull().sum()
#     # ---- plot part
#     f = plt.figure(figsize=(7,5))
#     ax = plt.gca()
#     for p in factors:
#         if similarity.loc[p, :].isnull().any():
#             continue
#         ax.errorbar(similarity.loc[p, "x"], similarity.loc[p, "y"],
#                     fmt="o", markersize=similarity.loc[p, "n"], 
#                     color=sourcesColor(p), alpha=0.5,
#                     xerr=similarity.loc[p, "xerr"],
#                     yerr=similarity.loc[p, "yerr"],
#                     label=p
#                    )
#     rect = mpatch.Rectangle((0, 0), width=1, height=0.4, facecolor="green", alpha=0.1, zorder=-1)
#     ax.add_patch(rect)
#     ax.set_xlim(0, 1.5)
#     ax.set_ylim(0, 1)
#     ax.set_xlabel("SID")
#     ax.set_ylabel("PD")
#     ax.set_title('Similarity between {station} profiles\nand all profiles'.format(station=stationRef))
#     plt.subplots_adjust(top=0.88, bottom=0.11, left=0.100, right=0.720)
#     ax.legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False)
#     return similarity

def get_all_SID_PD(PMF_profile, stations, source2=None, isRelativeMass=False):
    """
    Compute the SID and PD for all profiles in PMF_profile for the stations
    `stations`.
    """
    profiles = PMF_profile.dropna(axis=1, how='all').columns
    SID = pd.DataFrame(index=pd.MultiIndex.from_product((profiles, stations)),
                       columns=stations)
    PD = pd.DataFrame(index=pd.MultiIndex.from_product((profiles, stations)),
                      columns=stations)
    MAD = pd.DataFrame()

    list_stations1 = []
    if source2 is None:
        sameSource = True
    else:
        sameSource = False

    for p in profiles:
        if sameSource:
            source2 = p
        else:
            source2 = source2
        print(p)
        for station1 in stations:
            list_stations1.append(station1)
            if all(PMF_profile.loc[station1, p].isnull()):
                continue
            for station2 in stations:
                # if station2 in list_stations1:
                #     continue
                if all(PMF_profile.loc[station2, p].isnull()):
                    continue
                SID.loc[(p, station1), station2] = compute_SID(
                    PMF_profile.loc[station1, :],
                    PMF_profile.loc[station2, :],
                    p,
                    source2=source2,
                    isRelativeMass=isRelativeMass
                )
                PD.loc[(p, station1), station2] = compute_PD(
                    PMF_profile.loc[station1, :],
                    PMF_profile.loc[station2, :],
                    p,
                    source2=source2,
                    isRelativeMass=isRelativeMass
                )
        list_stations1 = []
    return (SID, PD)

def plot_similarity_profile(PMF_profile, SID, PD, err="ci", plotAll=False):
    """
    Plot a point in the SID/PD space (+/-err) for all profile in PMF_profile.

    PMF_profile: DataFrame in long-form. Column are the PMF factors.
    SID : DataFrame with index (factor, station) and column (station) : the SID matrix
    PD :  DataFrame with index (factor, station) and column (station) : the PD matrix
    err : "ci" or "sd", the type of error for xerr and yerr.
    plotAll: boolean. Either or not plot each pair of profile.

    Returns
    -------
    similarity : pd.DataFrame
        columns: x, y, xerr, yerr, n
        index: profiles
    handles_labels: tuple of handles and labels
        legend of the plot
    """
    from pyutilities.chemutilities import get_sourceColor

    similarity = pd.DataFrame(columns=["x", "y", "xerr", "yerr", "n"],
                              index=PMF_profile.columns)
    for p in PMF_profile.columns:
        if p not in SID.index:
            continue
        SIDtmp = SID.loc[p].dropna(axis=1, how='all').dropna(how='all')
        keep = np.triu(np.ones(SIDtmp.shape))
        SIDtmp = SIDtmp * keep
        SIDtmp = SIDtmp[SIDtmp > 0]
        SIDtmp = pd.Series(SIDtmp.values.reshape((SIDtmp.shape[0]*SIDtmp.shape[1])))
        SIDtmp = SIDtmp.dropna()
        PDtmp = PD.loc[p].dropna(axis=1, how='all').dropna(how='all')
        keep = np.triu(np.ones(PDtmp.shape))
        PDtmp = PDtmp * keep
        PDtmp = PDtmp[PDtmp > 0]
        PDtmp = pd.Series(PDtmp.values.reshape((PDtmp.shape[0]*PDtmp.shape[1])))
        PDtmp = PDtmp.dropna()

        similarity.loc[p, "x"] = SIDtmp.mean()
        similarity.loc[p, "y"] = PDtmp.mean()
        if err == "ci":
            similarity.loc[p, "xerr"] = SIDtmp.mean() \
                    - st.t.interval(0.95, len(SIDtmp)-1,
                                    loc=np.mean(SIDtmp),
                                    scale=st.sem(SIDtmp)
                                   )[0]
            similarity.loc[p, "yerr"] = PDtmp.mean() \
                    - st.t.interval(0.95, len(PDtmp)-1,
                                    loc=np.mean(PDtmp),
                                    scale=st.sem(PDtmp)
                                   )[0]
        elif err == "sd":
            similarity.loc[p, "xerr"] = SIDtmp.std()
            similarity.loc[p, "yerr"] = PDtmp.std()

        similarity.loc[p, "n"] = SIDtmp.notnull().sum()

        if plotAll:
            f = plt.figure(figsize=(7, 5))
            ax = plt.gca()
            ax.plot(SID.loc[p], PD.loc[p], "o", color=get_sourceColor(p))
            ax.set_title(p)                                       
            plot_deltatool_pretty(ax)
            # ax.set_xlabel("SID")                                          
            # ax.set_ylabel("PD") 
            # rect = mpatch.Rectangle((0, 0), width=1, height=0.4, facecolor="green", alpha=0.1, zorder=-1)
            # ax.add_patch(rect)                                                                           
            # ax.set_xlim(0, 1.5)                                                                          
            # ax.set_ylim(0, 1)
            plt.savefig("distance_all_profile_{p}.png".format(
                p=p.replace(" ", "-").replace("/", "-")
            ))
    # ---- plot part
    f = plt.figure(figsize=(8,5))
    ax = plt.gca()
    maxNumber = similarity.loc[:, "n"].max()
    for p in PMF_profile.columns:
        print(p)
        if similarity.loc[p, :].isnull().any():
            continue
        ax.errorbar(similarity.loc[p, "x"], similarity.loc[p, "y"],
                    fmt="o", markersize=14*similarity.loc[p, "n"]/maxNumber, 
                    color=get_sourceColor(p), alpha=0.5,
                    xerr=similarity.loc[p, "xerr"],
                    yerr=similarity.loc[p, "yerr"],
                    label=p
                   )
    plot_deltatool_pretty(ax)
    ax.set_title('Similarity between all pairs of profile')
    plt.subplots_adjust(top=0.88, bottom=0.11, left=0.100, right=0.700)
    handles, labels = ax.get_legend_handles_labels()
    newLabels = []
    for l in labels:
        newLabels.append("{l} ({n})".format(l=l.replace("_", " "), n=similarity.loc[l, "n"]))
    ax.legend(handles, newLabels, bbox_to_anchor=(1, 1), loc="upper left",
              frameon=False, fontsize=12)
    return (similarity, (handles, newLabels))

def plot_all_stations_similarity_by_source(PMF_profile):
    """
    Plot all individual pair of profile for each common source.
    """
    stations = PMF_profile.index.get_level_values("station").unique()
    sources = PMF_profile.columns

    for source in sources:
        if PMF_profile.loc[:, source].isnull().all():
            continue
        f, ax = plt.subplots()
        color = sourcesColor(source)
        for station1, station2 in itertools.product(stations, stations):
            if station1 == station2:
                continue
            plot_similarityplot(
                PMF_profile,
                station1, station2, source,
                ax=ax, plot_kw={"color":color}
            )
        plot_deltatool_pretty(ax=ax)
        f.suptitle(source)

def plot_relativeMass(PMF_profile, source="Biomass burning",
                      isRelativeMass=True, totalVar="PM10", naxe=1,
                      site_typologie=None):
    if not site_typologie:
        site_typologie = get_site_typology()
        # site_typologie = collections.OrderedDict()
        # site_typologie["Urban"] = ["Talence", "Lyon", "Poitiers", "Nice", "MRS-5av",
        #                            "PdB", "Aix-en-provence", "Nogent",
        #                            "Lens-2011-2012", "Lens-2013-2014"]
        # site_typologie["Valley"] = ["Chamonix", "GRE-fr"]
        # site_typologie["Traffic"] = ["Roubaix", "STG-cle"]
        # site_typologie["Rural"] = ["Revin"]


    carboneous = ["OC*", "EC"]
    ions = ["Cl-", "NO3-", "SO42-", "Na+", "NH4+", "K+", "Mg2+", "Ca2+"]
    organics = [
        "MSA", "Polyols", "Levoglucosan", "Mannosan",
    ]
    metals = [
        "Al", "As", "Ba", "Cd", "Co", "Cr", "Cs", "Cu", "Fe", "La", "Mn",
        "Mo", "Ni", "Pb", "Rb", "Sb", "Se", "Sn", "Sr", "Ti", "V", "Zn"
    ]
    keep_index = ["PM10"] + carboneous + ions + organics + metals
    if naxe==2:
        keep_index1 = carboneous + [
            "NO3-", "SO42-", "NH4+", "Na+", "K+", "Ca2+", "Mg2+", "Cl-",
            "Polyols", "Fe", "Al", "Levoglucosan", "Mannosan", "MSA"
        ]
        keep_index2 = []#"MSA"] 
        keep_index2tmp = list(set(keep_index) - set(keep_index1) - set(["PM10",
            "Mg2+"]))
        keep_index2tmp.sort()
        keep_index2 += keep_index2tmp
    dfperµg = pd.DataFrame(columns=keep_index)
    for station, df in PMF_profile.reset_index().groupby("station"):
        df = df.set_index("specie").drop("station", axis=1)
        if isRelativeMass:
            dfperµg.loc[station, :] = df.loc[:, source]
        else:
            dfperµg.loc[station, :] = to_relativeMass(df.loc[:, source],
                                                       totalVar=totalVar).T
    dfperµg = dfperµg.convert_objects(convert_numeric=True)

    # FIGURE
    f, axes = plt.subplots(nrows=naxe, ncols=1, figsize=(12.67, 6.54))
    if naxe == 1:
        axes = [axes, None]

    d = dfperµg.T.copy()
    d["specie"] = d.index
    # for i, keep_index in enumerate([keep_index1, keep_index2]):
    if naxe == 1:
        xtick_list = [keep_index]
    elif naxe == 2:
        xtick_list = [keep_index1, keep_index2]
    for i, keep_index in enumerate(xtick_list):
        dd = d.reindex(keep_index)
        # dd.rename(rename_station, inplace=True, axis="columns")
        dd = dd.melt(id_vars=["specie"])
        # if not percent:
        dd.replace({0: pd.np.nan}, inplace=True)
        axes[i].set_yscale("log")
        sns.boxplot(data=dd, x="specie", y="value", ax=axes[i], color="white",
                    showcaps=False,
                    showmeans=False, meanprops={"marker": "d"})
        ntypo = len(site_typologie.keys())
        colors = list(TABLEAU_COLORS.values())
        # for sp, specie in enumerate(keep_index):
        for t, typo in enumerate(site_typologie.keys()):
            if typo == "Urban":
                marker = "*"
            elif (typo == "Valley") or typo == ("Urban+Alps"):
                marker = "o"
            elif typo == "Traffic":
                marker = "p"
            else:
                marker = "d"
            step = 0.1
            j = 0
            for site in site_typologie[typo]:
                if site not in PMF_profile.index.get_level_values("station").unique():
                    continue
                axes[i].scatter(
                    pd.np.arange(0,len(keep_index))-ntypo*step/2+step/2+t*step,
                    d.loc[keep_index, site],
                    marker=marker, color=colors[j], alpha=0.8
                )
                j += 1

        # sns.swarmplot(data=dd, x="specie", y="value", color=".2", alpha=0.5, size=4,
        #              ax=axes[i]) 
    if naxe == 1:
        axes[0].set_ylim([1e-5,2])
    else:
        axes[0].set_ylim([1e-3,1])
        axes[1].set_ylim([1e-5,1e-2])
    for ax in axes:
        if ax:
            ax.set_ylabel("µg/µg of PM$_{10}$")
            ax.set_xlabel("")
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
            # ax.legend(loc="center", 
            #           ncol=(len(dfperµg.columns))//2,
            #           bbox_to_anchor=(0.5, -.55),
            #           frameon=False)
    # create custom legend
    labels = []
    artists = []
    for typo in site_typologie.keys():
        if typo == "Urban":
            marker = "*"
        elif (typo == "Valley") or (typo == "Urban+Alps"):
            marker = "o"
        elif typo == "Traffic":
            marker = "p"
        else:
            marker = "d"
        noSiteYet = True
        j = 0
        for site in site_typologie[typo]:
            if site not in PMF_profile.index.get_level_values("station").unique():
                continue
            if noSiteYet:
                artist = typo
                label = ""
                artists.append(artist)
                labels.append(label)
                noSiteYet = False
            artist = mlines.Line2D([], [], ls='', marker=marker, color=colors[j],
                                   label=site)
            artists.append(artist)
            labels.append(site)
            j += 1

    axes[0].legend(artists, labels, bbox_to_anchor=(1.1, 1),
                   handler_map={str: LegendTitle()},
                   frameon=False
              )
    # ax.legend('', frameon=False)

    f.suptitle(source)
    # plt.subplots_adjust(left=0.07, right=0.83, bottom=0.15, top=0.900, hspace=0.5)
    plt.subplots_adjust(
        top=0.925,
        bottom=0.07,
        left=0.053,
        right=0.899,
        hspace=0.489,
        wspace=0.2
    )

def save4deltaTool(contrib, profile):
    stations = profile.index.get_level_values('station').unique()
    for station in stations:
        F = profile.loc[station,:].dropna(how='all', axis=1)\
                .rename({"OC*":"OC", "SO42-": "SO4="})
        G = contrib.loc[station, :].dropna(how='all', axis=1)
        dfcontrib = G * F.loc['PM10']
        dfcontrib2specie = ((F.T/F.sum(axis=1)).T * 100).copy()

        fname = "{station}_4deltaTool.csv".format(station=station)
        F.index.name = "CON (µg/m3)"
        dfcontrib.index.name = "TREND (µg/m3)"
        dfcontrib2specie.index.name = "CONTR"

        F.to_csv("CONC.csv")
        dfcontrib.to_csv("TREND.csv", date_format='%m/%d/%Y')
        dfcontrib2specie.to_csv("CONTR.csv")
        with open(fname, "w") as f:
            f.write("MATRIX 1\n")
            f.write(F.to_csv())
            f.write("\nMATRIX 3\n")
            f.write(dfcontrib.to_csv(date_format='%m/%d/%Y'))
            f.write("\nMATRIX 4\n")
            f.write(dfcontrib2specie.to_csv())
        list_files = ["CONC.csv", "TREND.csv", "CONTR.csv"]
        with ZipFile(station+".zip", "w") as zip:
            for file in list_files:
                zip.write(file)
        for file in list_files:
            os.remove(file)

## TO ADAPT! ==================================================================
# conn = sqlite3.connect("/home/webersa/Documents/BdD/BdD_PM/db.sqlite")
# contrib = pd.read_sql(
#     "SELECT * FROM PMF_contributions WHERE programme IN ('SOURCES');",
#     con=conn,
#     parse_dates=["date"]
# )
# contrib.set_index(['station', 'date'], drop=True, inplace=True)
# contrib.drop(["programme", "index"], inplace=True, axis=1)
# profile = pd.read_sql(
#     "SELECT * FROM PMF_profiles WHERE programme IN ('SOURCES');",
#     con=conn
# )
# profile.set_index(['station', 'specie'], drop=True, inplace=True)
# profile.drop(["programme", "index"], inplace=True, axis=1)
### ===========================================================================

# stations = profile.index.get_level_values('station').unique()

