import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from py4pm.chemutilities import get_sourceColor, get_sourcesCategories, format_ions

class CachedAccessor:
     """
     Custom property-like object (descriptor) for caching accessors.
     Parameters
     ----------
     name : str
         The namespace this will be accessed under, e.g. ``df.foo``
     accessor : cls
         The class with the extension methods. The class' __init__ method
         should expect one of a ``Series``, ``DataFrame`` or ``Index`` as
         the single argument ``data``
     """

     def __init__(self, name, accessor):
         self._name = name
         self._accessor = accessor

     def __get__(self, obj, cls):
         if obj is None:
             # we're accessing the attribute of the class, i.e., Dataset.geo
             return self._accessor
         accessor_obj = self._accessor(obj)
         # Replace the property with the accessor object. Inspired by:
         # http://www.pydanny.com/cached-property.html
         # We need to use object.__setattr__ because we overwrite __setattr__ on
         # NDFrame
         object.__setattr__(obj, self._name, accessor_obj)
         return accessor_obj


class ReaderAccessor():
    def __init__(self, data):
        self._parent = data
    
    def read_metadata(self):
        """Get profiles, species and co

        It add a totalVariable (by default one of "PM10", "PM2.5", "PMrecons" or
        "PM10recons", "PM10rec"]. Otherwise, try to guess.

        """
        pmf = self._parent
        if pmf.dfprofiles_b is None:
            pmf.read.read_base_profiles()

        pmf.profiles = pmf.dfprofiles_b.columns.tolist()
        pmf.nprofiles = len(pmf.profiles)
        pmf.species = pmf.dfprofiles_b.index.tolist()
        pmf.nspecies = len(pmf.species)

        TOTALVAR = ["PM10", "PM2.5", "PMrecons", "PM10rec", "PM10recons"]
        for x in TOTALVAR:
            if x in pmf.species:
                pmf.totalVar = x
        if pmf.totalVar is None:
            print("Warning: trying to guess total variable.")
            pmf.totalVar = [x for x in pmf.species if "PM" in x]
            if len(pmf.totalVar) >= 1:
                print("Warning: several possible total variable: {}".format(pmf.totalVar))
                print("Watning: taking the first one {}".format(pmf.totalVar[0]))
            pmf.totalVar = pmf.totalVar[0]
        print("Total variable set to: {}".format(pmf.totalVar))

    def read_base_profiles(self):
        """TODO: Docstring for read_base_profiles.

        Returns
        -------
        TODO

        """
        pmf = self._parent

        dfbase = pd.read_excel(
                pmf._basename+"_base.xlsx",
                sheet_name=['Profiles'],
                header=None,
                )["Profiles"]

        idx = dfbase.iloc[:, 0].str.contains("Factor Profiles").fillna(False)
        idx = idx[idx].index.tolist()

        dfbase = dfbase.iloc[idx[0]:idx[1], 1:]
        dfbase.dropna(axis=0, how="all", inplace=True)
        factor_names = list(dfbase.iloc[0, 1:])
        dfbase.columns = ["specie"] + factor_names
        dfbase = dfbase\
                .drop(dfbase.index[0])\
                .set_index("specie")

        # check correct number of column
        idx = dfbase.columns.isna().argmax()
        if idx > 0:
            dfbase = dfbase.iloc[:, :idx]
            dfbase.dropna(how="all", inplace=True)
        # avoid 10**-12 possible concentration...
        dfbase = dfbase.infer_objects()
        dfbase[dfbase < 10e-6] = 0

        pmf.dfprofiles_b = dfbase

        pmf.read.read_metadata()

    def read_constrained_profiles(self):
        """TODO: Docstring for read_constrained_profiles.

        :returns: TODO

        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        dfcons = pd.read_excel(
                    pmf._basename+"_Constrained.xlsx",
                    sheet_name=['Profiles'],
                    header=None,
                )["Profiles"]

        idx = dfcons.iloc[:, 0].str.contains("Factor Profiles").fillna(False)
        idx = idx[idx].index.tolist()

        dfcons = dfcons.iloc[idx[0]:idx[1], 1:]
        dfcons.dropna(axis=0, how="all", inplace=True)

        # check correct number of column
        idx = dfcons.columns.isna().argmax()
        if idx > 0:
            dfcons = dfcons.iloc[:, :idx]
            dfcons.dropna(how="all", inplace=True)
        nancolumns = dfcons.isna().all()
        if nancolumns.any():
            dfcons = dfcons.loc[:, :nancolumns.idxmax()]
            dfcons.dropna(axis=1, how="all", inplace=True)

        dfcons.columns = ["specie"] + pmf.profiles
        dfcons = dfcons.set_index("specie")
        dfcons = dfcons[dfcons.index.notnull()]
        # avoid 10**-12 possible concentration...
        dfcons = dfcons.infer_objects()
        dfcons[dfcons < 10e-6] = 0

        pmf.dfprofiles_c = dfcons

    def read_base_contributions(self):
        """TODO: Docstring for read_base_contributions.

        :returns: TODO

        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        dfcontrib = pd.read_excel(
            pmf._basename+"_base.xlsx",
            sheet_name=['Contributions'],
            parse_date=["date"],
            header=None,
        )["Contributions"]

        try:
            idx = dfcontrib.iloc[:, 0].str.contains("Factor Contributions").fillna(False)
            idx = idx[idx].index.tolist()
            if len(idx) > 1:
                dfcontrib = dfcontrib.iloc[idx[0]:idx[1], :]
            else:
                dfcontrib = dfcontrib.iloc[idx[0]+1:, :]
        except AttributeError:
            print("WARNING: no total PM reconstructed in the file")

        dfcontrib.dropna(axis=1, how="all", inplace=True)
        dfcontrib.dropna(how="all", inplace=True)
        dfcontrib.drop(columns=dfcontrib.columns[0], inplace=True)
        dfcontrib.columns = ["Date"] + pmf.profiles
        dfcontrib.set_index("Date", inplace=True)
        dfcontrib = dfcontrib[dfcontrib.index.notnull()]

        dfcontrib.replace({-999: pd.np.nan}, inplace=True)

        pmf.dfcontrib_b = dfcontrib

    def read_constrained_contributions(self):
        """TODO: Docstring for read_constrained_contributions.

        :returns: TODO

        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        dfcontrib = pd.read_excel(
            pmf._basename+"_Constrained.xlsx",
            sheet_name=['Contributions'],
            parse_date=["date"],
            header=None,
        )["Contributions"]

        idx = dfcontrib.iloc[:, 0].str.contains("Factor Contributions").fillna(False)
        idx = idx[idx].index.tolist()

        if len(idx) > 1:
            dfcontrib = dfcontrib.iloc[idx[0]+1:idx[1], 1:]
        else:
            dfcontrib = dfcontrib.iloc[idx[0]+1:, 1:]

        nancolumns = dfcontrib.isna().all()
        if nancolumns.any():
            dfcontrib = dfcontrib.loc[:, :nancolumns.idxmax()]
        dfcontrib.dropna(axis=0, how="all", inplace=True)
        dfcontrib.dropna(axis=1, how="all", inplace=True)
        dfcontrib.columns = ["Date"] + pmf.profiles
        dfcontrib.replace({-999:pd.np.nan}, inplace=True)
        dfcontrib.set_index("Date", inplace=True)
        dfcontrib = dfcontrib[dfcontrib.index.notnull()]

        pmf.dfcontrib_c = dfcontrib.infer_objects()

    def read_base_bootstrap(self):
        """TODO: Docstring for read_base_bootstrap.

        :returns: TODO

        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        dfprofile_boot = pd.read_excel(
            pmf._basename+"_boot.xlsx",
            sheet_name=['Profiles'],
            header=None,
        )["Profiles"]

        dfbootstrap_mapping_b = dfprofile_boot.iloc[2:2+pmf.nprofiles, 0:pmf.nprofiles+2]
        dfbootstrap_mapping_b.columns = ["mapped"] + pmf.profiles + ["unmapped"]
        dfbootstrap_mapping_b.set_index("mapped", inplace=True)
        dfbootstrap_mapping_b.index = ["BF-"+f for f in pmf.profiles]

        idx = dfprofile_boot.iloc[:, 0].str.contains("Columns are:").fillna(False)
        idx = idx[idx].index.tolist()

        # 13 is the first column for BS result
        dfprofile_boot = dfprofile_boot.iloc[idx[0]+1:, 13:]
        df = pmf._split_df_by_nan(dfprofile_boot)

        df = pd.concat(df)
        df.index.names = ["specie", "profile"]
        # handle nonconvergente BS
        nBSconverged = dfbootstrap_mapping_b.sum(axis=1)[0]
        nBSnotconverged = len(df.columns)-1-nBSconverged
        if nBSnotconverged > 0:
            print("Warging: trying to exclude non-convergente BS")
            idxStrange = (df.loc[pmf.totalVar]>100)
            colStrange = df[idxStrange]\
                    .dropna(axis=1, how="all")\
                    .dropna(how="all")\
                    .columns
            print("BS eliminated:")
            print(df[colStrange])
            df = df.drop(colStrange, axis=1)

        # handle BS without totalVariable
        if pmf.totalVar:
            lowmass = (df.loc[pmf.totalVar, :] < 10**-3)
            if lowmass.any().any():
                print("Warning: BS with totalVar < 10**-3 encountered ({})".format(lowmass.any().sum()))
                df = df.loc[:, ~lowmass.any()]

        pmf.dfBS_profile_b = df
        pmf.dfbootstrap_mapping_b = dfbootstrap_mapping_b

    def read_constrained_bootstrap(self):
        """Read the '_Gcon_profile_boot.xlsx' and format it in dataframe

        Add:
        - self.dfBS_profile_c: all mapped profile
        - self.dfbootstrap_mapping_c: table of mapped profiles
        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        dfprofile_boot = pd.read_excel(
            pmf._basename+"_Gcon_profile_boot.xlsx",
            sheet_name=['Profiles'],
            header=None,
        )["Profiles"]

        dfbootstrap_mapping_c = dfprofile_boot.iloc[2:2+pmf.nprofiles, 0:pmf.nprofiles+2]
        dfbootstrap_mapping_c.columns = ["mapped"] + pmf.profiles + ["unmapped"]
        dfbootstrap_mapping_c.set_index("mapped", inplace=True)
        dfbootstrap_mapping_c.index = ["BF-"+f for f in pmf.profiles]

        idx = dfprofile_boot.iloc[:, 0].str.contains("Columns are:").fillna(False)
        idx = idx[idx].index.tolist()
        # 13 is the first column for BS result
        dfprofile_boot = dfprofile_boot.iloc[idx[0]+1:, 13:]
        df = pmf._split_df_by_nan(dfprofile_boot)

        df = pd.concat(df)
        df.index.names = ["specie", "profile"]
        # handle nonconvergente BS
        nBSconverged = dfbootstrap_mapping_c.sum(axis=1)[0]
        nBSnotconverged = len(df.columns)-1-nBSconverged
        if nBSnotconverged > 0:
            print("Warging: trying to exclude non-convergente BS")
            idxStrange = (df.loc[pmf.totalVar]>100)
            colStrange = df[idxStrange]\
                    .dropna(axis=1, how="all")\
                    .dropna(how="all")\
                    .columns
            print("BS eliminated: ", colStrange)
            df = df.drop(colStrange, axis=1)

        # handle BS without totalVariable
        if pmf.totalVar:
            lowmass = (df.loc[pmf.totalVar, :] < 10**-3)
            if lowmass.any().any():
                print("Warning: BS with totalVar < 10**-3 encountered ({})".format(lowmass.any().sum()))
                df = df.loc[:, ~lowmass.any()]

        pmf.dfBS_profile_c = df
        pmf.dfbootstrap_mapping_c = dfbootstrap_mapping_c

    def read_base_uncertainties_summary(self):
        """Read the _BaseErrorEstimationSummary.xlsx file

        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        if pmf.species is None:
            pmf.read.read_base_profiles()

        rawdf = pd.read_excel(
            pmf._basename+"_BaseErrorEstimationSummary.xlsx",
            sheet_name=["Error Estimation Summary"],
            header=None,
            )["Error Estimation Summary"]
        rawdf = rawdf.dropna(axis=0, how="all").reset_index().drop("index", axis=1)


        # ==== DISP swap
        idx = rawdf.iloc[:, 1].str.contains("Swaps").fillna(False)
        if idx.sum() > 0:
            df = pd.DataFrame()
            df = rawdf.loc[idx, :]\
                    .dropna(axis=1)\
                    .iloc[:, 1:]\
                    .reset_index(drop=True)
            df.columns = pmf.profiles
            df.index = ["swap count"]

            pmf.df_disp_swap_b = df

        # ==== uncertainties summary
        # get only the correct rows
        idx = rawdf.iloc[:, 0].str.contains("Concentrations for").astype(bool)
        idx = rawdf.loc[idx]\
                .iloc[:, 0]\
                .dropna()\
                .index
        df = pd.DataFrame()
        df = rawdf.loc[idx[0]+1:idx[-1]+1+pmf.nspecies, :]
        idx = df.iloc[:, 0].str.contains("Specie|Concentration").astype(bool)
        df = df.drop(idx[idx].index)
        df = df.dropna(axis=0, how='all')
        df["profile"] = pd.np.repeat(pmf.profiles, len(pmf.species)).tolist()

        df.columns = ["specie", "Base run", 
                "BS 5th", "BS 25th", "BS median", "BS 75th", "BS 95th", "tmp1",
                "BS-DISP 5th", "BS-DISP average", "BS-DISP 95th", "tmp2",
                "DISP Min", "DISP average", "DISP Max",
                "profile"
                ]
        df["specie"] = pmf.species * len(pmf.profiles)
        df.set_index(["profile", "specie"], inplace=True)
        df.drop(["tmp1", "tmp2"], axis=1, inplace=True)

        pmf.df_uncertainties_summary_b = df.infer_objects()

    def read_constrained_uncertainties_summary(self):
        """Read the _ConstrainedErrorEstimationSummary.xlsx file

        """
        pmf = self._parent

        if pmf.profiles is None:
            pmf.read.read_base_profiles()

        if pmf.species is None:
            pmf.read.read_base_profiles()

        rawdf = pd.read_excel(
            pmf._basename+"_ConstrainedErrorEstimationSummary.xlsx",
            sheet_name=["Constrained Error Est. Summary"],
            header=None,
            )["Constrained Error Est. Summary"]
        rawdf = rawdf.dropna(axis=0, how="all").reset_index().drop("index", axis=1)

        # ==== DISP swap
        idx = rawdf.iloc[:, 1].str.contains("Swaps").fillna(False)
        if idx.sum() > 0:
            df = pd.DataFrame()
            df = rawdf.loc[idx, :]\
                    .dropna(axis=1)\
                    .iloc[:, 1:]\
                    .reset_index(drop=True)
            df.columns = pmf.profiles
            df.index = ["swap count"]

            pmf.df_disp_swap_c = df

        # ==== uncertainties summary
        # get only the correct rows
        idx = rawdf.iloc[:, 0].str.contains("Concentrations for").astype(bool)
        idx = rawdf.loc[idx]\
                .iloc[:, 0]\
                .dropna()\
                .index
        df = pd.DataFrame()
        df = rawdf.loc[idx[0]+1:idx[-1]+1+pmf.nspecies, :]
        idx = df.iloc[:, 0].str.contains("Specie|Concentration").astype(bool)
        df = df.drop(idx[idx].index)
        df = df.dropna(axis=0, how='all')
        df["profile"] = pd.np.repeat(pmf.profiles, len(pmf.species)).tolist()

        df.columns = ["specie", "Constrained base run",
                "BS 5th", "BS median", "BS 95th", "tmp1",
                "BS-DISP 5th", "BS-DISP average", "BS-DISP 95th", "tmp2",
                "DISP Min", "DISP average", "DISP Max",
                "profile"
                ]
        df["specie"] = pmf.species * len(pmf.profiles)
        df.set_index(["profile", "specie"], inplace=True)
        df.drop(["tmp1", "tmp2"], axis=1, inplace=True)

        pmf.df_uncertainties_summary_c = df.infer_objects()


class PlotterAccessor():
    def __init__(self, data):
        self._parent = data

    def __call__(self):
        print("Called!")


    def _save_plot(self, formats=["png"], name="plot", DIR=""):
        """Save plot in a given format.
        
        :formats: list, format of the figure (see plt.savefig)
        :name: string, default "plot". File name.
        :DIR: string, default "". Directory.
        """
        for fmt in formats:
            plt.savefig("{DIR}{name}.{fmt}".format(DIR=DIR,
                                                   name=name.replace("/", "-"), fmt=fmt))

    def _plot_per_microgramm(self, constrained=True, df=None, profile=None, species=None,
                             new_figure=False, **kwargs):
        """
        """
        import seaborn as sns
        pmf = self._parent

        if new_figure:
            plt.figure(figsize=(12, 4))
            ax = plt.gca()
        elif "ax" in kwargs:
            ax = kwargs["ax"]

        if constrained:
            dfprofiles = pmf.dfprofiles_c
        else:
            dfprofiles = pmf.dfprofiles_b

        d = df.xs(profile, level="profile") \
                / (df.xs(profile, level="profile").loc[pmf.totalVar])
        d = d.reindex(species).unstack().reset_index()
        dref = dfprofiles[profile] / dfprofiles.loc[pmf.totalVar, profile]
        dref = dref.reset_index()
        sns.boxplot(data=d.replace({0: pd.np.nan}), x="specie", y=0,
                    color="grey", ax=ax)
        sns.stripplot(data=dref.replace({0: pd.np.nan}), x="specie", y=profile,
                      ax=ax, jitter=False, color="red")
        ax.set_yscale('log')
        ax.set_xticklabels(
            format_ions([t.get_text() for t in ax.get_xticklabels()]),
            rotation=90
        )
        ax.set_ylim((10e-6, 3))
        ax.set_ylabel("µg/µg")
        ax.set_xlabel("")
        ax.set_title(profile)

        #Create custom artists
        refArtist = plt.Line2D((0, 1),(0, 0), color='red', marker='o', linestyle='')
        BSArtist = plt.Rectangle((0, 0), 0, 0, color="grey")
        handles = [refArtist, BSArtist]
        labels = ["Ref. run", "BS"]
        ax.legend(handles=handles, labels=labels, loc="upper left", bbox_to_anchor=(1., 1.), frameon=False)

    def _plot_totalspeciesum(self, constrained=True, df=None, profile=None,
                             species=None, sumsp=None, new_figure=False,
                             **kwargs):
        import seaborn as sns
        """TODO: Docstring for _plot_totalspeciesum.

        :df: TODO
        :constrained: Boolean, either to use the constrained run or the base run
        :profile: TODO
        :species: TODO
        :sumsp: dataframe with the sum of each species
        :new_figure: TODO
        :returns: TODO

        """
        pmf = self._parent

        if new_figure:
            plt.figure(figsize=(12, 4))
            ax = plt.gca()
        elif "ax" in kwargs:
            ax = kwargs["ax"]

        if constrained:
            dfprofiles = pmf.dfprofiles_c
        else:
            dfprofiles = pmf.dfprofiles_b

        if sumsp is None:
            sumsp = pd.DataFrame(columns=species, index=['sum'])
            for sp in species:
                sumsp[sp] = df.loc[(sp, slice(None)), :].mean(axis=1).sum()

        d = df.xs(profile, level="profile").divide(sumsp.iloc[0], axis=0) * 100
        d = d.reindex(species).unstack().reset_index()
        dref = dfprofiles[profile].divide(dfprofiles.sum(axis=1)) * 100
        dref = dref.reset_index()
        sns.barplot(data=d, x="level_1", y=0, color="grey", ci="sd", ax=ax,
                    label="BS (sd)")
        sns.stripplot(data=dref, x="specie", y=0, color="red", jitter=False,
                      ax=ax, label="Ref. run")
        ax.set_xticklabels(
            format_ions([t.get_text() for t in ax.get_xticklabels()]),
            rotation=90
        )
        ax.set_ylim((0, 100))
        ax.set_ylabel("% of total specie sum")
        ax.set_xlabel("")
        ax.set_title(profile)

        h, l = ax.get_legend_handles_labels()
        h = h[-2:]
        l = l[-2:]
        ax.legend(handles=h, labels=l, loc="upper left", bbox_to_anchor=(1., 1.), frameon=False)

    def _plot_contrib(self, constrained=True, dfBS=None, dfDISP=None, dfcontrib=None,
                      profile=None, specie=None, BS=True, DISP=True,
                      BSDISP=False, new_figure=False, **kwargs):
        """TODO: Docstring for _plot_contrib.

        Parameters
        ----------

        dfBS : pd.DataFrame
        dfDISP: TODO
        :dfcontrib: TODO
        :profile: TODO
        :specie: TODO
        :BS: TODO
        :DISP: TODO
        :BSDISP: TODO
        :new_figure: TODO
        :returns: TODO

        """
        import seaborn as sns
        pmf = self._parent
        
        if new_figure:
            plt.figure(figsize=(12, 4))
            ax = plt.gca()
        elif "ax" in kwargs:
            ax = kwargs["ax"]

        if constrained:
            dfprofiles = pmf.dfprofiles_c
        else:
            dfprofiles = pmf.dfprofiles_b

        fill_kwarg = dict(
            alpha=0.5,
            edgecolor="black",
            linewidth=0,
        )
        with sns.axes_style("ticks"):
            if BS:
                d = pd.DataFrame(
                    columns=dfBS.columns,
                    index=dfcontrib.index
                )
                for BS in dfBS.columns:
                    d[BS] = dfcontrib[profile] * dfBS.xs(profile, level="profile").loc[specie][BS]
                mstd = d.std(axis=1)
                ma = d.mean(axis=1)
                plt.fill_between(
                    ma.index, ma-mstd, ma+mstd,
                    label="BS (sd)", **fill_kwarg
                )
                # d.mean(axis=1).plot(marker="*")
            if DISP:
                d = pd.DataFrame(
                    columns=dfDISP.columns,
                    index=dfcontrib.index
                )
                for DISP in ["DISP Min", "DISP Max"]:
                    d[DISP] = dfcontrib[profile] * dfDISP.xs(profile, level="profile").loc[specie][DISP]
                plt.fill_between(
                    d.index, d["DISP Min"], d["DISP Max"],
                    label="DISP (min-max)", **fill_kwarg
                )
            plt.plot(
                dfcontrib.index, dfcontrib[profile] * dfprofiles.loc[specie, profile],
                color="#888a85", marker="*", label="Ref. run"
            )
            ax.set_ylabel("Contribution to {} ($µg.m^{{-3}}$)".format(specie))
            ax.set_xlabel("")
            ax.set_title(profile)
            ax.legend(loc="upper left", bbox_to_anchor=(1., 1.), frameon=False)

    def _plot_profile(self, constrained=True, dfcontrib=None, dfBS=None, dfDISP=None, profile=None,
                      specie=None, BS=False, DISP=False, BSDISP=False):
        """TODO: Docstring for _plot_profile.

        :constrained: Boolean, either to use the constrained run or the base one
        :dfcontrib: TODO
        :dfprofile: TODO
        :profile: TODO
        :specie: TODO
        :BS: TODO
        :DISP: TODO
        :BSDISP: TODO
        :returns: TODO

        """
        pmf = self._parent

        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 12))
        axes[0].get_shared_x_axes().join(axes[0], axes[1])
        self._plot_per_microgramm(
            df=dfBS, constrained=constrained, profile=profile, species=pmf.species,
            new_figure=False, ax=axes[0]
        )

        self._plot_totalspeciesum(
            df=dfBS, constrained=constrained, profile=profile, species=pmf.species,
            new_figure=False, ax=axes[1]
        )

        self._plot_contrib(
            constrained=constrained,
            dfcontrib=dfcontrib,
            dfBS=dfBS, dfDISP=dfDISP,
            BS=BS, DISP=DISP,
            profile=profile, specie=specie,
            new_figure=False, ax=axes[2]
        )

        # axes[0].xaxis.tick_top()

        for ax in axes:
            ax.set_title("")
        fig.suptitle(profile)

        fig.subplots_adjust(
            top=0.95,
            bottom=0.05,
            left=0.125,
            right=0.865,
            hspace=0.40,
            wspace=0.015
        )

    def plot_per_microgramm(self, df=None, constrained=True, profiles=None, species=None,
            plot_save=False, BDIR=None):
        """Plot profiles in concentration unique (µg/m3).

        :df: DataFrame with multiindex [species, profile] and an arbitrary
            number of column.  Default to dfBS_profile_c.
        :constrained: Boolean, either to use the constrained run or the base run
        :profiles: list, profile to plot (one figure per profile)
        :species: list, specie to plot (x-axis)
        :plot_save: boolean, default False. Save the graph in BDIR.
        :BDIR: string, directory to save the plot.
        """
        pmf = self._parent

        if df is None:
            if constrained:
                if pmf.dfBS_profile_c is None:
                    pmf.read.read_constrained_bootstrap()
                df = pmf.dfBS_profile_c
                if pmf.dfprofiles_c is None:
                    pmf.read.read_constrained_profiles()
            else:
                if pmf.dfBS_profile_b is None:
                    pmf.read.read_base_bootstrap()
                df = pmf.dfBS_profile_b
                if pmf.dfprofiles_b is None:
                    pmf.read.read_base_profiles()
        elif not(isinstance(df, pd.DataFrame)):
            raise TypeError("df should be a pandas DataFrame.")
        

        if profiles is None:
            if pmf.profiles is None:
                pmf.read.read_metadata()
            profiles = pmf.profiles
        elif not(isinstance(profiles, list)):
            raise TypeError("profiles should be a list.")

        if species is None:
            if pmf.species is None:
                pmf.read.read_metadata()
            species = pmf.species
        elif not(isinstance(species, list)):
            raise TypeError("species should be a list.")

        if BDIR is None:
            BDIR = pmf._BDIR

        for p in profiles:
            self._plot_per_microgramm(df, constrained=constrained, profile=p, species=species,
                                      new_figure=True)
            plt.subplots_adjust(left=0.1, right=0.9, bottom=0.3, top=0.9)
            if plot_save: self._save_plot(DIR=BDIR, name=p+"_profile_perµg")

    def plot_totalspeciesum(self, df=None, profiles=None, species=None,
                            constrained=True,
                            plot_save=False, BDIR=None, **kwargs):
        """Plot profiles in percentage of total specie sum (%).

        :df: DataFrame with multiindex [species, profile] and an arbitrary
            number of column.  Default to dfBS_profile_c.
        :profiles: list, profile to plot (one figure per profile)
        :species: list, specie to plot (x-axis)
        :plot_save: boolean, default False. Save the graph in BDIR.
        :BDIR: string, directory to save the plot.
        """
        pmf = self._parent

        if df is None:
            if constrained:
                if pmf.dfBS_profile_c is None:
                    pmf.read.read_constrained_bootstrap()
                df = pmf.dfBS_profile_c
                if pmf.dfprofiles_c is None:
                    pmf.read.read_constrained_profiles()
            else:
                if pmf.dfBS_profile_b is None:
                    pmf.read.read_base_bootstrap()
                df = pmf.dfBS_profile_b
                if pmf.dfprofiles_b is None:
                    pmf.read.read_base_profiles()

        if profiles is None:
            if pmf.profiles is None:
                pmf.read.read_metadata()
            profiles = pmf.profiles

        if species is None:
            if pmf.species is None:
                pmf.read.read_metadata()
            species = pmf.species

        if BDIR is None:
            BDIR = pmf._BDIR

        new_figure = kwargs.pop("new_figure", False)

        sumsp = pd.DataFrame(columns=species, index=['sum'])
        for sp in species:
            sumsp[sp] = df.loc[(sp, slice(None)),:].mean(axis=1).sum()
        for p in profiles:
            self._plot_totalspeciesum(df, profile=p, species=species,
                                      sumsp=sumsp, new_figure=new_figure,
                                      **kwargs)
            plt.subplots_adjust(left=0.1, right=0.9, bottom=0.3, top=0.9)
            if plot_save:
                self._save_plot(DIR=BDIR, name=p+"_profile")

    def plot_contrib(self, dfBS=None, dfDISP=None, dfcontrib=None, profiles=None, specie=None,
                     constrained=True,
                     plot_save=False, BDIR=None, BS=True, DISP=True, BSDISP=False,
                     new_figure=True, **kwargs):
        """Plot temporal contribution in µg/m3.

        Parameters
        ----------

        df : pd.DataFrame, default self.dfBS_profile_c
            DataFrame with multiindex [species, profile] and an arbitrary number
            of column.
        dfcontrib : pd.DataFrame, default self.dfcontrib_c
            Profile as column and specie as index.
        profiles : list of string, default self.profiles
            profile to plot (one figure per profile)
        specie : string, default totalVar.
            specie to plot (y-axis)
        plot_save : boolean, default False
            Save the graph in BDIR.
        BDIR : string
            directory to save the plot
        """
        pmf = self._parent

        if (dfBS is None) and (BS):
            if constrained:
                if pmf.dfBS_profile_c is None:
                    pmf.read.read_constrained_bootstrap()
                dfBS = pmf.dfBS_profile_c
            else:
                if pmf.dfBS_profile_b is None:
                    pmf.read.read_base_bootstrap()
                dfBS = pmf.dfBS_profile_b

        if (dfDISP is None) and (DISP):
            if constrained:
                if pmf.df_uncertainties_summary_c is None:
                    pmf.read.read_constrained_uncertainties_summary()
                dfDISP = pmf.df_uncertainties_summary_c[["DISP Min", "DISP Max"]]
            else:
                if pmf.df_uncertainties_summary_b is None:
                    pmf.read.read_base_uncertainties_summary()
                dfDISP = pmf.df_uncertainties_summary_b[["DISP Min", "DISP Max"]]

        if dfcontrib is None:
            if constrained:
                if pmf.dfcontrib_c is None:
                    pmf.read.read_constrained_contributions()
                dfcontrib = pmf.dfcontrib_c
            else:
                if pmf.dfcontrib_b is None:
                    pmf.read.read_base_contributions()
                dfcontrib = pmf.dfcontrib_b

        if profiles is None:
            if pmf.profiles is None:
                pmf.read.read_metadata()
            profiles = pmf.profiles

        # if pmf.dfprofiles_c is None:
        #     pmf.read.read_constrained_profiles()

        if specie is None:
            if pmf.totalVar is None:
                pmf.read.read_metadata()
            specie = pmf.totalVar
        elif not isinstance(specie, str):
            raise ValueError(
                "`specie` should be a string, got {}.".format(specie)
            )

        if BDIR is None:
            BDIR = pmf._BDIR

        for p in profiles:
            self._plot_contrib(dfBS=dfBS, dfDISP=dfDISP,
                               dfcontrib=dfcontrib, constrained=constrained,
                               profile=p, specie=specie, BS=BS, DISP=DISP,
                               BSDISP=BSDISP, new_figure=new_figure,
                               **kwargs)
            plt.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.9)
            if plot_save:
                self._save_plot(DIR=BDIR, name=p+"_contribution")

    def plot_all_profiles(self, constrained=True, profiles=None, specie=None,
                          BS=True, DISP=True, BSDISP=False, plot_save=False,
                          savedir=None):
        """TODO: Docstring for plot_all_profiles.

        Parameters
        ----------

        constrained: Boolean, default True
            Either to use the constrained run or the base one
        profiles : list of string
            Profiles to plot
        species : ?
        {BS, DISP, BSDISP} : boolean, default True, True, False
            Use them as error estimation
        plot_save : boolean, default False
            Either or not saving the plot
        savedir : str
            Path to save the plot

        """
        pmf = self._parent

        if profiles is None:
            if pmf.profiles is None:
                pmf.read.read_metadata()
            profiles = pmf.profiles

        if BS:
            if constrained:
                if pmf.dfBS_profile_c is None:
                    pmf.read.read_constrained_bootstrap()
                dfBS = pmf.dfBS_profile_c
            else:
                if pmf.dfBS_profile_b is None:
                    pmf.read.read_base_bootstrap()
                dfBS = pmf.dfBS_profile_b
        else:
            dfBS = None

        if DISP:
            if constrained:
                if pmf.df_uncertainties_summary_c is None:
                    pmf.read.read_constrained_uncertainties_summary()
                dfDISP = pmf.df_uncertainties_summary_c[["DISP Min", "DISP Max"]]
            else:
                if pmf.df_uncertainties_summary_b is None:
                    pmf.read.read_base_uncertainties_summary()
                dfDISP = pmf.df_uncertainties_summary_b[["DISP Min", "DISP Max"]]
        else:
            dfDISP = None

        if constrained:
            if pmf.dfcontrib_c is None:
                pmf.read.read_constrained_contributions()
            dfcontrib = pmf.dfcontrib_c
        else:
            if pmf.dfcontrib_b is None:
                pmf.read.read_base_contributions()
            dfcontrib = pmf.dfcontrib_b

        if constrained:
            if pmf.dfprofiles_c is None:
                pmf.read.read_constrained_profiles()
        else:
            if pmf.dfprofiles_b is None:
                pmf.read.read_base_profiles()

        if specie is None:
            if pmf.totalVar is None:
                pmf.read.read_metadata()
            specie = pmf.totalVar

        if savedir is None:
            savedir = pmf._BDIR

        for p in profiles:
            self._plot_profile(
                constrained=constrained, dfcontrib=dfcontrib, dfBS=dfBS, dfDISP=dfDISP, profile=p,
                specie=specie, BS=BS, DISP=DISP, BSDISP=BSDISP
            )
            if plot_save:
                self._save_plot(
                    DIR=savedir,
                    name=pmf._site+"_"+p+"_contribution_and_profiles"
                )

    def plot_stacked_contribution(self, constrained=True, order=None, plot_kwargs=None):
        """Plot a stacked plot for the contribution

        :constrained: TODO
        :returns: TODO

        """
        pmf = self._parent

        df = pmf.to_cubic_meter(constrained=constrained)
        if order:
            if isinstance(order, list):
                df = df.reindex(order, axis=1)
            else:
                df = df.reindex(sorted(df.columns), axis=1)
        labels = df.columns

        y = pd.np.vstack(df.values).T
        colors = [
            get_sourceColor(c) for c in get_sourcesCategories(labels)
        ]
        
        fig, ax = plt.subplots()
        ax.stackplot(df.index, y, colors=colors, labels=labels)
        ax.set_ylabel(pmf.totalVar + "$\mu g/ m^{-3}$")
        ax.set_ylim(0, ax.get_ylim()[1])
        ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1., 1.))
        plt.subplots_adjust(
            top=0.961,
            bottom=0.081,
            left=0.037,
            right=0.887,
            hspace=0.2,
            wspace=0.2
        )

    def plot_seasonal_contribution(self, constrained=True, dfcontrib=None, dfprofiles=None, profiles=None,
            specie=None, plot_save=False, BDIR=None, annual=True,
            normalize=True, ax=None, barplot_kwarg={}):
        """Plot the relative contribution of the profiles.

        :dfcontrib: DataFrame with contribution as column and date as index.
        :dfprofiles: DataFrame with profile as column and specie as index.
        :profiles: list, profile to plot (one figure per profile)
        :specie: string, default totalVar. specie to plot
        :plot_save: boolean, default False. Save the graph in BDIR.
        :BDIR: string, directory to save the plot.
        :annual: plot annual contribution
        :normalize: plot relative contribution or absolute contribution.
        
        :returns: TODO

        """
        from py4pm.dateutilities import add_season
        pmf = self._parent

        if dfcontrib is None:
            if constrained:
                if pmf.dfcontrib_c is None:
                    pmf.read.read_constrained_contributions()
                dfcontrib = pmf.dfcontrib_c
            else:
                if pmf.dfcontrib_b is None:
                    pmf.read.read_base_contributions()
                dfcontrib = pmf.dfcontrib_b

        if dfprofiles is None:
            if constrained:
                if pmf.dfprofiles_c is None:
                    pmf.read.read_constrained_profiles()
                dfprofiles = pmf.dfprofiles_c
            else:
                if pmf.dfprofiles_b is None:
                    pmf.read.read_base_profiles()
                dfprofiles = pmf.dfprofiles_b

        if profiles is None:
            if pmf.profiles is None:
                pmf.read.read_metadata()
            profiles = pmf.profiles

        if specie is None:
            if pmf.totalVar is None:
                pmf.read.read_metadata()
            specie = pmf.totalVar

        if BDIR is None:
            BDIR = pmf._BDIR

        if ax is None:
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 4.7))

        df = pmf.get_seasonal_contribution(specie=specie, normalize=normalize,
                                            annual=annual,
                                           constrained=constrained)
        c = get_sourceColor()
        colors = c.loc["color", get_sourcesCategories(df.columns)]

        df.index = [l.replace("_", " ") for l in df.index]
        axes = df.plot.bar(
            stacked=True,
            rot=0,
            color=colors,
            ax=ax,
            **barplot_kwarg
        )

        ax.set_ylabel("Normalized contribution" if normalize else "$µg.m^{-3}$")
        if normalize:
            ax.set_ylim([0, 1])
            ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
        ax.legend("", frameon=False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left',
                  bbox_to_anchor=(1, 0.5), frameon=False)
        ax.set_xlabel("")
        ax.set_title(specie)
        plt.subplots_adjust(top=0.90, bottom=0.10, left=0.15, right=0.72)
        
        if plot_save:
            title = "_seasonal_contribution_{}".format(
                    "normalized" if normalize else "absolute"
                    )
            self._save_plot(DIR=BDIR, name=pmf._site+title)
        
        return (df)

    def plot_stacked_profiles(self, constrained=True):
        """plot the repartition of the species among the profiles, normalized to
        100%

        Parameters
        ----------
        constrained : boolean, default True
            use the constrained run or not

        Returns
        -------
        ax : the axe
        """
        pmf = self._parent

        df = pmf.get_total_specie_sum(constrained=constrained)

        df = df.sort_index(axis=1)

        colors = [get_sourceColor(c) for c in df.columns]

        fig, ax = plt.subplots(1, 1, figsize=(12, 4))
        df.plot(kind="bar", stacked=True, color=colors, ax=ax)

        xticklabels = [t.get_text() for t in ax.get_xticklabels()]
        ax.set_xticklabels(format_ions(xticklabels), rotation=90)
        ax.set_xlabel("")

        ax.yaxis.set_major_formatter(mticker.PercentFormatter())
        ax.set_ylabel("Normalized contribution (%)")
        ax.set_ylim(0, 100)

        h, l = ax.get_legend_handles_labels()
        h = reversed(h)
        l = reversed(l)
        ax.legend(h, l, loc="upper left", bbox_to_anchor=(1, 1), frameon=False)
        fig.subplots_adjust(bottom=0.275, top=0.96, left=0.09, right=0.83)

        return ax


class PMF(object):

    """PMF are able to read file from US EPA PMF5.0 software output (in xlsx
    format), then parse them in a more handy format (pandas DataFrame).
    Several plot utilities are also available.
    """

    def __init__(self, site, BDIR, program=None):
        """Create a PMF object from the xlsx files output of EPAPMF5.

        :site: str, the name of the site
        :BDIR: str, the directory where the xlsx files live

        """
        self._site = site
        self._BDIR = BDIR
        self._program = program

        self._basename = BDIR+site
        self.profiles = None
        self.nprofiles = None
        self.species = None
        self.nspecies = None
        self.totalVar = None
        self.dfprofiles_b = None
        self.dfcontrib_b = None
        self.dfprofiles_c = None
        self.dfcontrib_c = None
        self.dfBS_profile_b = None
        self.dfBS_profile_c = None
        self.df_disp_swap_b = None
        self.df_disp_swap_c = None
        self.df_uncertainties_summary_b = None
        self.df_uncertainties_summary_c = None

    def to_cubic_meter(self, constrained=True, specie=None, profiles=None):
        """Convert the contribution in cubic meter for the given specie

        Parameters
        ----------

        specie : str, the specie, default totalVar
        profiles : list of profile, default all profiles

        Return
        ------
        df : dataframe

        """
        if specie is None:
            specie = self.totalVar

        if profiles is None:
            profiles = self.profiles

        if constrained:
            df = self.dfcontrib_c
            dfprofiles = self.dfprofiles_c
        else:
            df = self.dfcontrib_b
            dfprofiles = self.dfprofiles_b

        contrib = pd.DataFrame(index=df.index, columns=profiles)

        for profile in profiles:
            contrib[profile] = df[profile] * dfprofiles.loc[specie, profile]

        return contrib

    def to_relative_mass(self, constrained=True, species=None, profiles=None):
        """Compute the factor profile relative mass (i.e. each species divided
        by the totalVar mass)

        :constrained: TODO
        :species: TODO
        :profiles: TODO
        :returns: TODO

        """
        if constrained:
            df = self.dfprofiles_c
        else:
            df = self.dfprofiles_b

        if profiles is None:
            profiles = self.profiles

        if species is None:
            species = self.species

        d = df[profiles] / df.loc[self.totalVar, profiles]

        return d

    def get_total_specie_sum(self, constrained=True):
        """
        Return the total specie sum profiles in %

        Parameters
        ----------
        constrained : boolean, default True
            use the constrained run or not

        Returns
        -------
        df : pd.DataFrame
            The normalized species sum per profiles
        """
        if constrained:
            df = self.dfprofiles_c.copy()
        else:
            df = self.dfprofiles_b.copy()

        # df = (self.dfprofiles_c.T / self.dfprofiles_c.sum(axis=1)).T * 100
        df = (df.T / df.sum(axis=1)).T * 100
        return df

    def get_seasonal_contribution(self, specie=None, annual=True,
                                  normalize=True, constrained=True):
        """
        Get a dataframe of seasonal contribution
        """
        from py4pm.dateutilities import add_season

        if constrained:
            if self.dfprofiles_c is None:
                self.read.read_constrained_profiles()
            if self.dfcontrib_c is None:
                self.read.read_constrained_contributions()
            dfprofiles = self.dfprofiles_c
            dfcontrib = self.dfcontrib_c
        else:
            if self.dfprofiles_b is None:
                self.read.read_base_profiles()
            if self.dfcontrib_b is None:
                self.read.read_base_contributions()
            dfprofiles = self.dfprofiles_b
            dfcontrib = self.dfcontrib_b

        if specie is None:
            if self.totalVar is None:
                self.read.read_metadata()
            specie = self.totalVar


        dfcontribSeason = (dfprofiles.loc[specie] * dfcontrib).sort_index(axis=1)
        ordered_season = ["Winter", "Spring", "Summer", "Fall"]
        if annual:
            ordered_season.append("Annual")

        dfcontribSeason = add_season(dfcontribSeason, month=False)\
                .infer_objects()
        dfcontribSeason = dfcontribSeason.groupby("season")

        if normalize:
            df = (dfcontribSeason.sum().T / dfcontribSeason.sum().sum(axis=1))
            df = df.T
        else:
            df = dfcontribSeason.mean()

        if annual:
            df.loc["Annual", :] = df.mean()

        df = df.reindex(ordered_season)

        return df

    def replace_totalVar(self, newTotalVar):
        """replace the total var to all dataframe

        :newTotalVar: TODO
        :returns: TODO

        """
        DF = [
            self.dfprofiles_b,
            self.dfprofiles_c,
            self.dfBS_profile_b,
            self.dfBS_profile_c,
            self.df_uncertainties_summary_b,
            self.df_uncertainties_summary_c,
        ]
        for df in DF:
            if df is None:
                continue
            df.rename({self.totalVar: newTotalVar}, inplace=True, axis=0)

        self.species = [newTotalVar if x == self.totalVar else x for x in self.species]
        self.totalVar = newTotalVar

    def rename_profile_to_profile_category(self):
        """Rename the factor profile name to match the category
        :returns: TODO

        """
        DF = [
            self.dfprofiles_b,
            self.dfprofiles_c,
            self.dfcontrib_b,
            self.dfcontrib_c,
            self.dfBS_profile_b,
            self.dfBS_profile_c,
            self.df_uncertainties_summary_b,
            self.df_uncertainties_summary_c,
        ]
        for df in DF:
            if df is None:
                continue
            possible_sources = {
                p: get_sourcesCategories([p])[0]
                for p in self.profiles
            }
            df.rename(possible_sources, inplace=True, axis=1)
            df.rename(possible_sources, inplace=True, axis=0)

        self.profiles = [possible_sources[p] for p in self.profiles]
        

    def recompute_new_species(self, specie=None):
        """Recompute a specie given the other species. For instance, recompute OC
        from OC* and a list of organic species.

        It modify inplace both dfprofile_b and dfprofile_c, and update
        self.species.

        :specie: str in ["OC",]

        """
        knownSpecies = ["OC"]
        if specie not in knownSpecies:
            return

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

        if specie == "OC":
            if specie not in self.species:
                self.species.append(specie)
            OCb = self.dfprofiles_b.loc["OC*"].copy()
            OCc = self.dfprofiles_c.loc["OC*"].copy()
            for sp in equivC.keys():
                if sp in self.species:
                    OCb += (self.dfprofiles_b.loc[sp] * equivC[sp]).infer_objects()
                    OCc += (self.dfprofiles_c.loc[sp] * equivC[sp]).infer_objects()
            self.dfprofiles_b.loc[specie] = OCb.infer_objects()
            self.dfprofiles_c.loc[specie] = OCc.infer_objects()


    def _split_df_by_nan(self, df):
        """TODO: Docstring for split_df_by_nan.

        :df: TODO
        :returns: TODO

        """
        d = {}
        dftmp = df.dropna()
        for i, sp in enumerate(self.species):
            d[sp] = dftmp.iloc[self.nprofiles*i:self.nprofiles*(i+1), :]
            d[sp].index = self.profiles
            d[sp].index.name = "profile"
            d[sp].columns = ["Boot{}".format(i) for i in range(len(d[sp].columns))]
        return d


    def print_base_uncertainties_summary(self, profiles=None,
            species=None, BDIR=None):
        """Print the uncertainties given by BS, BS-DISP and DISP
        
        :profiles: list, list of profiles, default all profiles
        :species: list, list of species, default all species
        
        :returns: pd.DataFrame, BS, DISP and BS-DISP ranges
        """
        import seaborn as sns

        if self.df_uncertainties_summary_b is None:
            self.read.read_base_uncertainties_summary()

        if profiles is None:
            if self.profiles is None:
                self.read.read_metadata()
            profiles = self.profiles

        if species is None:
            if self.species is None:
                self.read.read_metadata()
            species = self.species

        if BDIR is None:
            BDIR = self._BDIR

        df = self.df_uncertainties_summary_b

        return df.T.loc[:, (profiles, species)]

    def print_constrained_uncertainties_summary(self, profiles=None,
            species=None, BDIR=None):
        """Plot the uncertainties given by BS, BS-DISP and DISP
        
        :profiles: list, list of profiles, default all profiles
        :species: list, list of species, default all species
        
        :returns: pd.DataFrame, BS, DISP and BS-DISP ranges
        """
        import seaborn as sns

        if self.df_uncertainties_summary_c is None:
            self.read.read_constrained_uncertainties_summary()

        if profiles is None:
            if self.profiles is None:
                self.read.read_metadata()
            profiles = self.profiles

        if species is None:
            if self.species is None:
                self.read.read_metadata()
            species = self.species

        if BDIR is None:
            BDIR = self._BDIR

        df = self.df_uncertainties_summary_c

        return df.T.loc[:, (profiles, species)]

    read = CachedAccessor("read", ReaderAccessor)
    plot = CachedAccessor("plot", PlotterAccessor)
