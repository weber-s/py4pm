import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd
import seaborn as sns
import sqlite3
from py4pm.dateutilities import add_season

MAPPER_METALS_NAME_TO_SYMBOLE = {
    "Potassium": "K",
    "Calcium": "Ca",
    "Titanium": "Ti",
    "Vanadium": "V",
    "Chromium": "Cr",
    "Manganese": "Mn",
    "Iron": "Fe",
    "Cobalt": "Co",
    "Nickel": "Ni",
    "Copper": "Cu",
    "Zinc": "Zn",
    "Gallium": "Ga",
    "Germanium": "Ge",
    "Arsenic": "As",
    "Selenium": "Se",
    "Bromine": "Br",
    "Yttrium": "Yt",
    "Molybdenum": "Mo",
    "Cadmium": "Cd",
    "Tin": "Sn",
    "Antimony": "Sb",
    "Mercury": "Hg",
    "Thallium": "Tl",
    "Lead": "Pb",
    "Bismuth": "Bi",
}

def replace_QL(dftmp, species=None, conn=None):
    """Replace the -1 and -2 in the dataframe by the appropriate DL and QL
    values

    The change are done inplace.

    :dftmp: pandas DataFrame
    """
    stations = dftmp.station.unique()

    if species is None:
        species = dftmp.columns

    if conn is None:
        conn = sqlite3.connect("/home/webersa/Documents/BdD/BdD_PM/aerosols.db")
    
    sqlquery = """
        SELECT {sp} FROM QL 
        WHERE station IN ("{stations}")
        AND "sample ID" LIKE "%QL%";
        """.format(
            sp='", "'.join(species),
            stations='", "'.join(stations)
        )
    print(sqlquery)
    QLtmp = pd.read_sql(sqlquery, con=conn)
    print(QLtmp)
    conn.close()
    QLtmp = QLtmp.apply(pd.to_numeric, errors='ignore').dropna(how="all", axis=1)
    for station in stations:
        QLtmpmean = QLtmp[QLtmp.station==station].mean()
        to_replace = {
            c: {-2: QLtmpmean[c]/2, -1: QLtmpmean[c]/2} for c in QLtmpmean.index
        }
        for c in dftmp.columns:
            if c not in species: continue
            if (c in to_replace.keys()) and (pd.notna(to_replace[c][-1])):
                idx = dftmp.station == station
                dftmp.loc[idx, c] = dftmp.loc[idx, c].clip_lower(to_replace[c][-1])

def get_sourceColor(source=None):
    """Return the hexadecimal color of the source(s) 

    If no option, then return the whole dictionary
    
    Optional Parameters
    ===================

    source : str
        The name of the source
    """
    color = {
        "Traffic": "#000000",
        "Traffic 1": "#000000",
        "Traffic 2": "#102262",
        "Road traffic": "#000000",
        "Primary traffic": "#000000",
        "Traffic_ind": "#000000",
        "Traffic_exhaust": "#000000",
        "Traffic_dir": "#444444",
        "Traffic_non-exhaust": "#444444",
        "Resuspended_dust": "#444444",
        "Oil/Vehicular": "#000000",
        "Road traffic/oil combustion": "#000000",
        "Biomass_burning": "#92d050",
        "Biomass burning": "#92d050",
        "Biomass_burning1": "#92d050",
        "Biomass_burning2": "#92d050",
        "Sulfate-rich": "#ff2a2a",
        "Sulfate_rich": "#ff2a2a",
        "Sulfate rich": "#ff2a2a",
        "Nitrate-rich": "#217ecb", # "#ff7f2a",
        "Nitrate_rich": "#217ecb", # "#ff7f2a",
        "Nitrate rich": "#217ecb", # "#ff7f2a",
        "Secondary_inorganics": "#0000cc",
        "MSA_rich": "#ff7f2a", # 8c564b",
        "Secondary_oxidation": "#ff87dc",
        "Marine SOA": "#ff7f2a", # 8c564b",
        "Biogenic SOA": "#8c564b",
        "Anthropogenic SOA": "#8c564b",
        "Marine/HFO": "#a37f15", #8c564b",
        "Aged seasalt/HFO": "#8c564b",
        "Marine_biogenic": "#fc564b",
        "HFO": "#70564b",
        "HFO (stainless)": "#70564b",
        "Oil": "#70564b",
        "Vanadium rich": "#70564b",
        "Cadmium rich": "#70564b",
        "Marine": "#33b0f6",
        "Marin": "#33b0f6",
        "Salt": "#00b0f0",
        "Seasalt": "#00b0f0",
        "Sea-road salt": "#209ecc",
        "Sea/road salt": "#209ecc",
        "Fresh seasalt": "#00b0f0",
        "Aged_salt": "#97bdff", #00b0f0",
        "Aged seasalt": "#97bdff", #00b0f0",
        "Fungal spores": "#ffc000",
        "Primary_biogenic": "#ffc000",
        "Primary biogenic": "#ffc000",
        "Biogenique": "#ffc000",
        "Biogenic": "#ffc000",
        "Dust": "#dac6a2",
        "Mineral dust": "#dac6a2",
        "Crustal_dust": "#dac6a2",
        "Industrial": "#7030a0",
        "Industries": "#7030a0",
        "Indus/veh": "#5c304b",
        "Industry/traffic": "#5c304b", #7030a0",
        "Arcellor": "#7030a0",
        "Siderurgie": "#7030a0",
        "Plant debris": "#2aff80",
        "Plant_debris": "#2aff80",
        "Débris végétaux": "#2aff80",
        "Choride": "#80e5ff",
        "PM other": "#cccccc",
        "Traffic/dust (Mix)": "#333333",
        "SOA/sulfate (Mix)": "#6c362b",
        "Sulfate rich/HFO": "#8c56b4",
        "nan": "#ffffff"
    }
    color = pd.DataFrame(index=["color"], data=color)
    if source:
        if source not in color.keys():
            print("WARNING: no {} found in colors".format(source))
            return "#666666"
        return color.loc["color", source]
    else:
        return color

def get_sourcesCategories(profiles):
    """Get the sources category according to the sources name.

    Ex. Aged sea salt → Aged_sea_salt

    :profiles: list
    :returns: list

    """
    possible_sources = {
        "Vehicular": "Traffic",
        "VEH": "Traffic",
        "VEH ind": "Traffic_ind",
        "Traffic_exhaust": "Traffic_exhaust",
        "Traffic_non-exhaust": "Traffic_non-exhaust",
        "VEH dir": "Traffic_dir",
        "Oil/Vehicular": "Traffic",
        "Oil": "Oil",
        "Vanadium rich": "Vanadium rich",
        "Road traffic/oil combustion": "Traffic",
        "Traffic": "Road traffic",
        "Traffic 1": "Traffic 1",
        "Traffic 2": "Traffic 2",
        "Primary traffic": "Road traffic",
        "Road traffic": "Road traffic",
        "Road trafic": "Road traffic",
        "Road traffic/dust": "Traffic/dust (Mix)",
        "Bio. burning": "Biomass_burning",
        "Bio burning": "Biomass_burning",
        "Comb fossile/biomasse": "Biomass_burning",
        "BB": "Biomass_burning",
        "Biomass_burning": "Biomass_burning",
        "Biomass Burning": "Biomass_burning",
        "Biomass burning": "Biomass_burning",
        "BB1": "Biomass_burning1",
        "BB2": "Biomass_burning2",
        "Sulfate-rich": "Sulfate_rich",
        "Sulphate-rich": "Sulfate_rich",
        "Nitrate-rich": "Nitrate_rich",
        "Sulfate rich": "Sulfate_rich",
        "Sulfate_rich": "Sulfate_rich",
        "Nitrate rich": "Nitrate_rich",
        "Nitrate_rich": "Nitrate_rich",
        "Secondary inorganics": "Secondary_inorganics",
        "Secondaire": "MSA_rich",
        "Secondary bio": "MSA_rich",
        "Secondary biogenic": "MSA_rich",
        "Secondary organic": "MSA_rich",
        "Secondary oxidation": "Secondary_oxidation",
        "Secondaire organique": "MSA_rich",
        # "Marine SOA": "Marine SOA",
        "Marine SOA": "MSA_rich",
        "MSA_rich": "MSA_rich",
        "MSA rich": "MSA_rich",
        "Secondary biogenic/sulfate": "SOA/sulfate (Mix)",
        "Marine SOA/SO4": "SOA/sulfate (Mix)",
        "Marine/HFO": "Marine/HFO",
        "Marine biogenic/HFO": "Marine/HFO",
        "Secondary biogenic/HFO": "Marine/HFO",
        "Marine bio/HFO": "Marine/HFO",
        "Marin bio/HFO": "Marine/HFO",
        "Sulfate rich/HFO": "Marine/HFO",
        "Marine secondary": "MSA_rich",
        "Marin secondaire": "MSA_rich",
        "HFO": "HFO",
        "HFO (stainless)": "HFO",
        "Marin": "MSA_rich",
        "Sea/road salt": "Sea-road salt",
        "Sea-road salt": "Sea-road salt",
        "sea-road salt": "Sea-road salt",
        "Road salt": "Salt",
        "Sea salt": "Salt",
        "Seasalt": "Salt",
        "Salt": "Salt",
        "Fresh seasalt": "Salt",
        "Sels de mer": "Salt",
        "Aged_salt": "Aged_salt",
        "Aged sea salt": "Aged_salt",
        "Aged seasalt": "Aged_salt",
        "Aged seasalt": "Aged_salt",
        "Aged salt": "Aged_salt",
        "Primary_biogenic": "Primary_biogenic",
        "Primary bio": "Primary_biogenic",
        "Primary biogenic": "Primary_biogenic",
        "Biogénique primaire": "Primary_biogenic",
        "Biogenique": "Primary_biogenic",
        "Biogenic": "Primary_biogenic",
        "Mineral dust": "Dust",
        "Mineral dust ": "Dust",
        "Resuspended_dust": "Resuspended_dust",
        "Resuspended dust": "Resuspended_dust",
        "Dust": "Dust",
        "Crustal dust": "Dust",
        "Dust (mineral)": "Dust",
        "Dust/biogénique marin": "Dust",
        "AOS/dust": "Dust",
        "Industrial": "Industrial",
        "Industry": "Industrial",
        "Industrie": "Industrial",
        "Industries": "Industrial",
        "Industry/vehicular": "Industry/traffic",
        "Industry/traffic": "Industry/traffic",
        "Industries/trafic": "Industry/traffic",
        "Cadmium rich": "Cadmium rich",
        "Fioul lourd": "HFO",
        "Arcellor": "Industrial",
        "Siderurgie": "Industrial",
        "Débris végétaux": "Plant_debris",
        "Chlorure": "Chloride",
        "PM other": "Other"
        }
    s = [possible_sources[k] for k in profiles]
    return s

def get_site_typology():
    import collections
    
    site_typologie = collections.OrderedDict()
    site_typologie["Urban"] = ["Talence", "Lyon", "Poitiers", "Nice", "MRS-5av",
                               "PdB", "Aix-en-provence", "Nogent", "Poitiers",
                               "Lens-2011-2012", "Lens-2013-2014", "Lens", "Rouen"]
    site_typologie["Valley"] = ["Chamonix", "Passy", "Marnaz", "GRE-cb", "VIF",
                                "GRE-fr", "Passy_decombio"]
    site_typologie["Traffic"] = ["Roubaix", "STG-cle"]
    site_typologie["Rural"] = ["Revin", "Peyrusse", "ANDRA-PM10", "ANDRA-PM2.5"]

    site_typologie_SOURCES = collections.OrderedDict()
    site_typologie_SOURCES["Urban"] = [
        "LEN", "LY", "MRS", "NGT", "NIC", "POI", "PdB", "PROV", "TAL", "ROU"
    ]
    site_typologie_SOURCES["Valley"] = ["CHAM", "GRE"]
    site_typologie_SOURCES["Traffic"] = ["RBX", "STRAS"]
    site_typologie_SOURCES["Rural"] = ["REV"]

    for typo in site_typologie.keys():
        site_typologie[typo] += site_typologie_SOURCES[typo]

    return site_typologie

def get_OC_from_OC_star_and_organic(df):
    """
    Re-compute OC taking into account the organic species

    OC = OC* + sum(eqC_sp)
    """
    OC = df.loc['OC*'].copy()
    equivC = {
        'Oxalate': 0.27,
        'Arabitol': 0.40,
        'Mannitol': 0.40,
        'Sorbitol': 0.40,
        'Polyols': 0.40,
        'Levoglucosan': 0.44,
        'Mannosan': 0.44,
        'Galactosan': 0.44,
        'MSA': 0.12,
        'Glucose': 0.44,
        'Cellulose': 0.44,
        'Maleic': 0.41,
        'Succinic': 0.41,
        'Citraconic': 0.46,
        'Glutaric': 0.45,
        'Oxoheptanedioic': 0.48,
        'MethylSuccinic': 0.53,
        'Adipic': 0.49,
        'Methylglutaric': 0.49,
        '3-MBTCA': 0.47,
        'Phtalic': 0.58,
        'Pinic': 0.58,
        'Suberic': 0.55,
        'Azelaic': 0.57,
        'Sebacic': 0.59,
    }
    for sp in equivC.keys():
        if sp in df.index:
            OC += df.loc[sp] * equivC[sp]
    return OC

def get_sample_where(sites=None, date_min=None, date_max=None, species=None,
                     min_sample=None, particle_size=None, con=None):
    """Get dataframe that meet conditions

    :sites: TODO
    :date_min: TODO
    :date_max: TODO
    :min_sample: int, minimum samples size
    :particle_size:
    :con: sqlite3 connection
    :returns: TODO

    """
    df = pd.read_sql("SELECT * FROM values_all;", con=con)

    df["Date"] = pd.to_datetime(df["Date"])
    if date_min:
        df = df.loc[date_min < df["Date"]]
    if date_max:
        df = df.loc[df["Date"] < date_max]
    if species:
        df = df.loc[df[species].notnull().all(axis=1)]
    if particle_size:
        df["Station"] = df["Station"]+"—"+df["Particle_size"]
    if min_sample:
        keep_stations = df.groupby("Station").size()
        keep_stations = list(keep_stations.loc[keep_stations > min_sample].index)
        df = df.loc[df["Station"].isin(keep_stations)]
    return df


def _format_ions(text):
    map_ions = {
        "Cl-": "Cl$^-$",
        "Na+": "Na$^+$",
        "K+": "K$^+$",
        "NO3-": "NO$_3^-$",
        "NH4+": "NH$_4^+$",
        "SO42-": "SO$_4^{2-}$",
        "Mg2+": "Mg$^{2+}$",
        "Ca2+": "Ca$^{2+}$",
        "nss-SO42-": "nss-SO$_4^{2-}$"
    }
    if text in map_ions.keys():
        return map_ions[text]
    else:
        return text

def format_ions(text):
    if isinstance(text, list):
        mapped = [_format_ions(x) for x in text]
    elif isinstance(text, str):
        mapped = _format_ions(text)
    else:
        raise KeyError(
            "`text` must be a {x,y}ticklabels, a list of string or string"
        )
    return mapped

class plot():
    
    def _mainComponentOfPM(dff, station):
        COLORS = {
            "OM": "#008000",
            "EC": "#000000",
            "Cl-": "#59B2B2",
            "NO3-": "#0000FF",
            "SO42-": "#FF0000",
            "NH4+": "#FF8000",
            "Ca2+": "#CED770",
            "Other ions": "#710077",
            "Metals": "#804000",
            "Anhydrous monosaccharides": "#004000",
            "Organic acids": "#CE9E8E",
            "Polyols": "#A0A015",
            "Oxalate": "#7D0000",
            "MSA": "#2D00BB",
            "Glucose": "#4B8A08",
            "Cellulose": "#0B3B0B",
            "HULIS": "#58ACFA"
        }
        TEXTCOLORS = {
            "OM": "#000000",
            "EC": "#FFFFFF",
            "Cl-": "#000000",
            "NO3-": "#FFFFFF",
            "SO42-": "#000000",
            "NH4+": "#000000",
            "Ca2+": "#000000",
            "Other ions": "#FFFFFF",
            "Metals": "#FFFFFF",
            "Anhydrous monosaccharides": "#FFFFFF",
            "Organic acids": "#000000",
            "Polyols": "#000000",
            "Oxalate": "#FFFFFF",
            "MSA": "#FFFFFF",
            "Glucose": "#FFFFFF",
            "Cellulose": "#FFFFFF",
            "HULIS": "#000000"
        }

        ORGANICS = ["HULIS", "Anhydrous monosaccharides", "Polyols", "Organic acids", "Oxalate",
                    "MSA", "Glucose", "Cellulose"]

        # 2 dataframes: one for the 'main' components, one for the organics
        df_proportion_perday = pd.DataFrame()
        nonorganics = list(set(dff.columns)-set(ORGANICS))
        for c in nonorganics:
            df_proportion_perday[c] = dff[c]/dff[nonorganics].sum(axis=1)

        df_proportion_OM_perday = pd.DataFrame()
        for c in dff.columns:
            if c in ORGANICS:
                df_proportion_OM_perday[c] = dff[c]/dff["OM"]

        d = pd.DataFrame(index=list(df_proportion_OM_perday.columns) +
                         list(df_proportion_perday.columns))

        d["other"] = df_proportion_perday.median()
        d["organics"] = df_proportion_OM_perday.median()
        d.loc[ORGANICS, "other"] = pd.np.nan
        df_mg_per_gOM = df_proportion_OM_perday.median() * 1000

        # Plot part
        order1 = ["OM", "EC", "Cl-", "NO3-", "SO42-", "NH4+", "Ca2+", "Other ions",
                  "Metals"]
        order2 = ORGANICS.copy()

        d = d.reindex(order1+order2)
        d.dropna(axis=0, how="all", inplace=True)

        OMidentified = df_proportion_OM_perday.median().sum() * 100
        dnormalize = d/d.sum() * 100

        # d1 = d["other"].reindex(order1, axis=0)
        # d2 = d["organics"].reindex(order2, axis=0)
        # d1 = d1/d1.sum()
        # d2 = d2/d2.sum()


        f, ax = plt.subplots(figsize=(9.5,7.5))
        dnormalize.T.plot.bar(
            stacked=True,
            color=dnormalize.index.map(COLORS).dropna(),
            rot=0,
            ax=ax,
        )

        xpos = {"other": 0, "organics": 1}
        texts = {"other": [], "organics": []}
        for xvar in ["other", "organics"]:
            val = dnormalize[xvar].reset_index().melt(id_vars=["index"])
            cumsum = 0
            for i, v in zip(val["index"], val["value"]):
                if pd.np.isnan(v):
                    continue
                cumsum += v
                if xvar == "other":
                    annot = ax.annotate("{}".format(format_ions(i)), 
                                       (xpos[xvar]-0.28, (cumsum -v/2) ),
                                       ha="right",
                                       va="center"
                                      )
                else:
                    text = "{}\n({:.2f} mg.g$_{{OM}}^{{-1}}$)".format(i, df_mg_per_gOM.loc[i])
                    if len(text)<40:
                        text = text.replace("\n", " ")

                    annot = ax.annotate(text, 
                                        (xpos[xvar]+0.28, (cumsum -v/2) ),
                                        ha="left",
                                        va="center"
                                       )
                texts[xvar].append(annot)

                
                ax.annotate("{:.0f}%".format(v),
                            (xpos[xvar], (cumsum - v/2) ),
                            ha="center",
                            va="center",
                            color=TEXTCOLORS[i],
                            fontweight="bold"
                           )
        # texts = pd.Series(plt.gcf().get_children()[1].get_children())
        # idx = [type(t)==matplotlib.text.Annotation for t in texts]
        # texts = texts[idx].tolist()

        # adjust_text(
        #     texts["organics"],
        #     # arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
        #     autoalign='', only_move={'points': 'y', 'text': 'y'}
        # )

        yOMidentified = OMidentified * dnormalize.loc["OM", "other"]/100
        ax.annotate("{:.0f}% identified".format(OMidentified),
                    (xpos["other"], yOMidentified/2),
                    ha="center",
                    va="center",
                    color="#FFFFFF",
                    fontweight="bold"
                   )
        ax.plot([0.25, 0.75], [yOMidentified, 100], "-k")
        ax.plot([-0.25, 0.25], [yOMidentified, yOMidentified], '-w')

        ax.set_title(station, fontsize=16)
        ax.set_xticklabels([])
        f.subplots_adjust(top=0.88,
                         bottom=0.11,
                         left=0.125,
                         right=0.85,
                         hspace=0.2,
                         wspace=0.2)
        ax.legend('', frameon=False)
        ax.yaxis.set_major_formatter(FuncFormatter('{0:.0f}%'.format))
        sns.despine()


    def mainCompentOfPM(station, dateStart, dateEnd, seasonal=False,
                        savefig=False, savedir=None):
        """
        Plot a stacked bar plot of the different constitutant of the PM

        Parameters
        ----------

        station : str
            name of the station
        dateStart, dateEnd : str
            starting and ending date
        seasonal : boolean, default False
            Either to make separate graph per season
        savefig : boolean, default False
            Save the fig in png and pdf
        savedir : str path, default None
            Where to save the figures
        """
        TO_GROUP = {
            "Metals": [
                "Al", "As", "Cd", "Cr", "Cu", "Fe", "Mn", "Mo", "Ni", "Pb", "Rb", "Sb",
                "Se", "Sn", "Ti", "V", "Zn"
            ],
            "Anhydrous monosaccharides":  ["Levoglucosan", "Mannosan", "Galactosan"],
            "Polyols":  ["Arabitol", "Sorbitol", "Mannitol"],
            "Organic acids": [
                "Maleic", "Succinic", "Citraconic", "Glutaric", "Oxoheptanedioic",
                "MethylSuccinic", "Adipic", "Methylglutaric", "3-MBTCA", "Phtalic",
                "Pinic", "Suberic", "Azelaic", "Sebacic"
            ],
            "Other ions": [
                "Na+", "K+", "Mg2+",
            ]
        }

        TO_MICROGRAMME = ["OM", "EC", "HULIS"]

        conn = sqlite3.connect("/home/webersa/Documents/BdD/BdD_PM/aerosols.db")
        df = pd.read_sql(
            "SELECT * FROM values_all WHERE station IN ('{}');".format(station),
            con=conn
        )

        df.date = pd.to_datetime(df.date)
        df.set_index("date", inplace=True, drop=True)

        df = df[(dateStart < df.index) & (df.index < dateEnd)]

        if seasonal:
            df = add_season(df)

        # Metals = [
        #     "Al", "As", "Ba", "Cd", "Co", "Cr", "Cs", "Cu", "Fe", "La", "Mn",
        #     "Mo", "Ni", "Pb", "Rb", "Sb", "Se", "Sn", "Sr", "Ti", "V", "Zn"
        # ]


        dff = pd.DataFrame()
        for k in TO_GROUP.keys():
            df[k] = df[TO_GROUP[k]].sum(axis=1, min_count=1)
        
        # Get only the columns we have
        dff = df.reindex(TO_GROUP.keys(), axis=1)
        dff["OM"] = df["OC"]*1.8
        to_keep = ["EC", "NO3-", "NH4+", "Cl-", "SO42-", "Ca2+", "Oxalate", "MSA",
                   "Glucose", "Cellulose", "HULIS"]
        for k in to_keep:
            if k in df.columns:
                dff[k] = df[k]
        dff.apply(pd.to_numeric)

        if seasonal:
            dff["season"] = df["season"]
        
        # Convert ng to µg
        for i in TO_MICROGRAMME:
            dff[i] *= 1000
        
        
        DF = []
        seasonName = []
        if seasonal:
            for season in df["season"].unique():
                DF.append(dff[dff["season"] == season].drop("season", axis=1))
                seasonName.append(season)
        else:
            DF = [dff]
            seasonName = ["annual"]

        for dfff, season in zip(DF, seasonName):
            plot._mainComponentOfPM(dfff, station)
            ax = plt.gca()
            if season:
                title = ax.get_title()
                plt.title(title+" "+season)
            if savefig:
                plt.savefig(
                    "{BDIR}/{station}_{temp}.png".format(
                        BDIR=savedir, station=station, temp=season
                    )
                )
                plt.savefig(
                    "{BDIR}/{station}_{temp}.pdf".format(
                        BDIR=savedir, station=station, temp=season
                    )
                )


    def what_do_we_have(sites=None, date_min=None, date_max=None, species=None,
                        min_sample=None, particle_size=None, con=None):
        """TODO: Docstring for what_do_we_have.

        :sites: TODO
        :date_min: TODO
        :date_max: TODO
        :species: TODO
        :min_sample: TODO
        :con: TODO
        :returns: TODO

        """
        df = get_sample_where(
            sites=sites,
            date_min=date_min,
            date_max=date_max,
            species=species,
            min_sample=min_sample,
            particle_size=particle_size,
            con=con
        )

        df.set_index(["Station", "Date"], inplace=True)
        stations = df.index.get_level_values("Station").unique()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
        for i, station in enumerate(stations):
            date = df.loc[station].index
            ax.plot(date, [i]*len(date), "-o", label=station)

        ax.set_yticks(range(len(stations)))
        ax.set_ylim(-0.5, len(stations)-0.5)
        ax.set_yticklabels(stations)

