import os
import pandas as pd
import numpy as np
from zipfile import ZipFile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import scipy.stats as st
import seaborn as sns
from matplotlib.colors import TABLEAU_COLORS
import matplotlib.lines as mlines
import matplotlib.text as mtext
import itertools

from py4pm.chemutilities import *

class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, orig_handle, **self.text_props)
        handlebox.add_artist(title)
        return title

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

def compute_SID(df1, df2, factor1=None, factor2=None, isRelativeMass=True):
    """
    Compute the SID of the factors `factor1` and `factor2` in the profile
    `df1`and `df2`.
    """

    if isinstance(df1, pd.Series) and isinstance(df2, pd.Series):
        p1 = df1
        p2 = df2
    else:
        if not(factor2):
            factor2 = factor1

        if factor1 not in df1.dropna(axis=1, how="all").columns:
            return np.nan

        if factor2 not in df2.dropna(axis=1, how="all").columns:
            return np.nan

        p1 = df1.loc[:, factor1]
        p2 = df2.loc[:, factor2]

    if not isRelativeMass:
        p1 = to_relativeMass(p1)
        p2 = to_relativeMass(p2)

    sp = p1.index.intersection(p2.index)
    if len(sp) > 3:
        ID = pd.np.abs(p1[sp]-p2[sp])
        MAD = p1[sp] + p2[sp]
        SID = pd.np.sqrt(2)/len(sp) * (ID/MAD).sum()
    else:
        SID = pd.np.nan

    return SID

def compute_PD(df1, df2, factor1=None, factor2=None, isRelativeMass=True):
    """
    Compute the PD of the factors `factor1` and `factor2` in the profile
    `df1`and `df2`.
    """

    if isinstance(df1, pd.Series) and isinstance(df2, pd.Series):
        p1 = df1
        p2 = df2
    else:
        if not(factor2):
            factor2 = factor1

        if factor1 not in df1.dropna(axis=1, how="all").columns:
            return np.nan

        if factor2 not in df2.dropna(axis=1, how="all").columns:
            return np.nan

        p1 = df1.loc[:, factor1]
        p2 = df2.loc[:, factor2]

    if not isRelativeMass:
        p1 = to_relativeMass(p1)
        p2 = to_relativeMass(p2)

    if p1.index.str.contains("PM").any():
        p1 = p1.loc[~p1.index.str.contains("PM")]
    if p2.index.str.contains("PM").any():
        p2 = p2.loc[~p2.index.str.contains("PM")]

    # Remove this part (impact performance)
    # p1.dropna(inplace=True)
    # p2.dropna(inplace=True)
    sp = p1.index.intersection(p2.index)
    if len(sp) > 3:
        PD = 1 - pd.np.corrcoef(p1[sp], p2[sp])[1, 0]**2
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

def get_all_SID_PD(PMF_profile, stations, factor2=None, isRelativeMass=False):
    """
    Compute the SID and PD for all profiles in PMF_profile for the stations
    `stations`.
    """
    factors = PMF_profile.dropna(axis=1, how='all').columns
    SID = pd.DataFrame(index=pd.MultiIndex.from_product((factors, stations)),
                       columns=list(stations))
    PD = pd.DataFrame(index=pd.MultiIndex.from_product((factors, stations)),
                      columns=list(stations))
    MAD = pd.DataFrame()

    list_stations1 = []
    if factor2 is None:
        sameFactor = True
    else:
        sameFactor = False

    for factor1 in factors:
        if sameFactor:
            factor2 = factor1
        else:
            factor2 = factor2
        print(factor1)
        for station1 in stations:
            list_stations1.append(station1)
            if all(PMF_profile.loc[station1, factor1].isnull()):
                continue
            for station2 in stations:
                # print(station1, station2)
                # if station2 in list_stations1:
                #     continue
                if all(PMF_profile.loc[station2, factor1].isnull()):
                    continue
                profile1 = PMF_profile.loc[station1, factor1]
                profile2 = PMF_profile.loc[station2, factor2]
                SID.loc[(factor1, station1), station2] = compute_SID(
                    profile1,
                    profile2,
                    factor1=factor1,
                    factor2=factor2,
                    isRelativeMass=isRelativeMass
                )
                PD.loc[(factor1, station1), station2] = compute_PD(
                    profile1,
                    profile2,
                    factor1=factor1,
                    factor2=factor2,
                    isRelativeMass=isRelativeMass
                )
        list_stations1 = []
    return (SID, PD)

def plot_similarity_profile(SID, PD, err="ci", plotAll=False):
    """
    Plot a point in the SID/PD space (+/-err) for all profile in SID and PD.

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
    from py4pm.chemutilities import get_sourceColor

    
    stations = list(SID.index.get_level_values("station").unique())
    profiles = list(SID.index.get_level_values("profile").unique())
    similarity = pd.DataFrame(columns=["x", "y", "xerr", "yerr", "n"],
                              index=profiles)
    for p in profiles:
        if p not in SID.index:
            continue
        x = SID.loc[(p, stations), stations]  # .sort_index(axis=1)
        x = x.where(np.triu(x, k=1).astype(bool))
        x = x.reset_index().melt(id_vars=["profile", "station"]).dropna()
        x = x.round(3)
        y = PD.loc[(p, stations), stations]  # .sort_index(axis=1)
        y = y.where(np.triu(y, k=1).astype(bool))
        y = y.reset_index().melt(id_vars=["profile", "station"]).dropna()
        y = y.round(3)

        similarity.loc[p, "x"] = x["value"].mean()
        similarity.loc[p, "y"] = y["value"].mean()
        if err == "ci":
            similarity.loc[p, "xerr"] = x["value"].mean() \
                    - st.t.interval(
                        0.95, len(x["value"])-1,
                        loc=np.mean(x["value"]),
                        scale=st.sem(x["value"])
                    )[0]
            similarity.loc[p, "yerr"] = y["value"].mean() \
                    - st.t.interval(
                        0.95, len(y["value"])-1,
                        loc=np.mean(y["value"]),
                        scale=st.sem(y["value"])
                    )[0]
        elif err == "sd":
            similarity.loc[p, "xerr"] = x["value"].std()
            similarity.loc[p, "yerr"] = y["value"].std()
        similarity.loc[:, "xerr"] = similarity.loc[:, "xerr"].fillna(0)
        similarity.loc[:, "yerr"] = similarity.loc[:, "yerr"].fillna(0)

        similarity.loc[p, "n"] = x["value"].notnull().sum()

        if plotAll:
            f = plt.figure(figsize=(7, 5))
            ax = plt.gca()
            ax.plot(x["value"], y["value"], "o", color=get_sourceColor(p))
            ax.set_title(p)
            plot_deltatool_pretty(ax)
            # ax.set_xlabel("SID")
            # ax.set_ylabel("PD")
            # rect = mpatch.Rectangle((0, 0), width=1, height=0.4, facecolor="green", alpha=0.1, zorder=-1)
            # ax.add_patch(rect)
            # ax.set_xlim(0, 1.5)
            # ax.set_ylim(0, 1)
            plt.savefig("distance_all_profile_{p}.pdf".format(
                p=p.replace(" ", "-").replace("/", "-")
            ))
    # ---- plot part
    f = plt.figure(figsize=(8,5))
    ax = plt.gca()
    maxNumber = similarity.loc[:, "n"].max()
    for p in profiles:
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

def get_profile_from_PMF(pmfs):
    """Get a profile matrix from a list of PMF object

    :pmfs: TODO
    :returns: TODO

    """
    profiles = pd.DataFrame()
    for pmf in pmfs:
        dftmp = pmf.to_relative_mass()
        possible_sources = {
            p: get_sourcesCategories([p])[0]
            for p in dftmp.columns
        }
        dftmp.rename(columns=possible_sources, inplace=True)
        dftmp.sort_index(axis=1, inplace=True)
        dftmp["station"] = pmf._site+"_"+pmf._program
        profiles = pd.concat([profiles, dftmp])
    return profiles


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

