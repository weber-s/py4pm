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

In order to retrieve the contributions of the profiles from the contrained run
(in `GRE-cb_Contrained.xlsx`, sheet "contributions", seconda table for total variable) :

```python
grecb.read.read_constrained_contributions()

```

will go through `GRE-cb_Base.xlsx` file to retreive the list of PMF factor, try to infer
the total variable (if any), and read the contribution from the constrained run.
And now, the `grecb` object has a `dfcontrib_c` attribute (`_c` for the constrained run):

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

which is the `G` matrix, in normalized unit. In order to have it in `µg/m³`, we need to
know the chemical profile from the constrained run, so run :

```python
grecb.read.read_constrained_profile()

```

and `grecb` has now a not null `dfprofile_c` dataframe :

```python
>>> grecb.dfprofile_c
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

And we can easily reconstruct the time serie in `µg/m³` of each specie for every profile
by simple multiplication of the timeserie by the concentration in the chemical profile.
`to_cubic_metter` does just that :

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
can specify other condition (see [the doc of this method](api_py4pm.html#py4pm.pmfutilities.PMF.to_cubic_meter)).

### Retreive uncertainties

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


