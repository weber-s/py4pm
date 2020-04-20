Examples
========

Load the PMF
------------

Given an output folder in `/home/myname/Documents/PMF/GRE-cb/MobilAir_woOrga` that looked like:

```
MobilAir_woOrga
├── GRE-cb_BaseErrorEstimationSummary.xlsx
├── GRE-cb_base.xlsx
├── GRE-cb_boot.xlsx
├── GRE-cb_ConstrainedDISPest.dat
├── GRE-cb_ConstrainedDISPres1.txt
├── GRE-cb_ConstrainedDISPres2.txt
├── GRE-cb_ConstrainedDISPres3.txt
├── GRE-cb_ConstrainedDISPres4.txt
├── GRE-cb_ConstrainedErrorEstimationSummary.xlsx
├── GRE-cb_Constrained.xlsx
├── GRE-cb_diagnostics.xlsx
├── GRE-cb_DISPest.dat
├── GRE-cb_DISPres1.txt
├── GRE-cb_DISPres2.txt
├── GRE-cb_DISPres3.txt
├── GRE-cb_DISPres4.txt
├── GRE-cb_Gcon_profile_boot.xlsx
├── GRE-cb_rotational_comments.txt
└── GRE-cb_sourcecontributions.xls

```

in order to convert them to a PMF object, run the following command :

```python
from py4pm.pmfutilities import PMF

grecb = PMF(site="GRE-cb", BDIR="/home/myname/Documents/PMF/GRE-fr/MobilAir_woOrga")

```

Now, `grecb` is an instance of a PMF object, and has a lot of `reader` and `ploter`.

Read the data
-------------

### Organization

The `read` class of the PMF object give access to different reader to retreive data from
the different xlsx files outputed by the EPA PMF5 software.

They all start by `read_base*` or `read_constrained*` name, for the base and constrained
run, respectively.

The special method `read_metadata` is used to retrieve the factors names and species names
from the `_base.xlsx` files, and use them everywhere else. It also try to set the total
variable name if any (one of PM10, PM2.5, PMrecons, PM10rec, PM10recons, otherwise try to
guess), used to convert unit and to be the default variable to plot.

For now, the following readers are implemented :

 - [read_metadata](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_metadata)
 - [read_base_contributions](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_base_contributions)
 - [read_base_profiles](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_base_profiles)
 - [read_base_bootstrap](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_base_bootstrap)
 - [read_base_uncertainties_summary](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_base_uncertainties_summary)
 - [read_constrained_contributions](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_constrained_contributions)
 - [read_constrained_profiles](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_constrained_profiles)
 - [read_constrained_bootstrap](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_constrained_bootstrap)
 - [read_constrained_uncertainties_summary](api_py4pm.html#py4pm.pmfutilities.ReaderAccessor.read_constrained_uncertainties_summary)

### Contribution

The contributions of the factors (`G` matrix) are read from the `_base.xlsx` and
`_Constrained.xlsx` files, sheet `contributions`.
You can read them using the reader `read_base_contributions` and
`read_constrained_contributions`:

```python
grecb.read.read_base_contributions()
grecb.read.read_constrained_contributions()

```

And now, the `grecb` object has a `dfcontrib_b` and `dfcontrib_c` attributes (`_b` for the
base run, `_c` for the constrained run):

```python
>>> grecb.dfcontrib_c

            Sulfate-rich  Nitrate-rich  ...  Biomass burning  Sea/road salt  Mineral dust
Date                                    ...                                              
2017-02-28      0.321580     -0.105980  ...          0.19419       0.606290      0.182880
2017-03-03      0.429480     -0.038802  ...          0.61595       0.050129      0.382890
2017-03-06     -0.098123     -0.151530  ...          0.53346       4.636400      0.272410
2017-03-09      0.643500     -0.002527  ...          1.09060       0.153200      1.083600
2017-03-12      0.664090      0.308390  ...          1.70740      -0.200000      0.846930

```

which is the `G` matrix, in normalized unit. 

### Chemical profiles

The chemical profiles (or simply profiles) is the `F` matrix of the PMF (in `µg/m³`) and
are read from the `_base.xslx` and `_Constrained.xlsx` files, sheet `Profiles`.
You can read them using the reader `read_base_profiles` and `read_constrained_profiles`:

```python
grecb.read.read_base_profile()
grecb.read.read_constrained_profile()

```

and `grecb` has now a not null `dfprofiles_b` and `dfprofiles_c` dataframe :

```python
>>> grecb.dfprofiles_c
              Sulfate-rich  Nitrate-rich ... Biomass burning  Sea/road salt  Mineral dust
specie                                   ...                                             
PMrecons          4.402500      2.421300 ...        3.027900       0.364280      2.009600
OC*               1.225300      0.000000 ...        1.308900       0.041038      0.428110
EC                0.162970      0.000000 ...        0.347050       0.019199      0.030703
Cl-               0.000000      0.002425 ...        0.026819       0.109070      0.000000
NO3-              0.300660      1.702200 ...        0.093396       0.000000      0.000000
SO42-             0.977680      0.010441 ...        0.092800       0.032969      0.189890
...

```

### Uncertainties

You can also read the bootstrap and DISP result :

```python
grecb.read.read_constrained_summary()

```

and now, you have access to :

  * `grecb.df_uncertainties_summary` : the summary of the BS and DISP uncertainties for
      each profiles and species.

And if you want to retreive the whole BS file :

```python
grecb.read.read_constrained_bootstrap()

```

and now, you have access to :

  * `grecb.dfBS_profile_c` : all the bootstrap chemical profiles
  * `grecb.dfbootstrape_mapping_c` : the table of the mapping between base and BS factors

### Utilities

In order to have the contributions in `µg/m³`, which is given by `G⋅F`, we need to know
both the chemical profile `F` and the contribution `G`.
And we can easily reconstruct the time serie in `µg/m³` of each specie for every profile
by simple multiplication of the timeserie by the concentration in the chemical profile.
Since this is a very often computation, the method `to_cubic_metter` does just that :

```python
>>> grecb.to_cubic_metter()
            Sulfate-rich  Nitrate-rich  ... Biomass burning  Sea/road salt  Mineral dust
Date                                    ...                                             
2017-02-28      1.415756     -0.256609  ...        0.587988       0.220859      0.367516
2017-03-03      1.890786     -0.093951  ...        1.865035       0.018261      0.769456
2017-03-06     -0.431987     -0.366900  ...        1.615264       1.688948      0.547435
2017-03-09      2.833009     -0.006120  ...        3.302228       0.055808      2.177603
2017-03-12      2.923656      0.746705  ...        5.169836      -0.072856      1.701991
...                  ...           ...  ...             ...            ...           ...


```

Note that `to_cubic_metter` use by default the contrained run, all the profile and the total variable, but you
can specify other conditions (see [the doc of this method](api_py4pm.html#py4pm.pmfutilities.PMF.to_cubic_meter)).


Plot utilities
--------------

### Chemical profile (per microgram of total variable)

### Chemical profile (in percentage of the sum of each species)

### Contribution time series and uncertainties

```python
grecb.plot.plot_contrib(profiles=["Primary biogenic"])

```

will produce the following graph 

```eval_rst
.. figure:: images/timeseries_POA.png
   :scale: 50 %
   :alt: Time series of POA
   :align: center

   Primary biogenic factor contribution to the total variable.

```

Since the EPA PMF5 does not output the chemical profile (F) matrix of the boostrap, the uncertainties is estimated by computing the species concentration given the F matrix of the reference run and the G matrix of the bootstrap run. As a result, the output is "hacky" since in the bootstrap method, bith the F and G matrix are changing. If you want to remove them, just pass `BS=False` to the method.


