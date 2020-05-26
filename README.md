[![Documentation Status](https://img.shields.io/badge/Documentation-API-green)](https://webersa.gricad-pages.univ-grenoble-alpes.fr/py4pm/)

Utilities to handle PMF data
============================

Contents
--------

* [Usage](usage.md)
* [API](api_py4pm.md)

Description
-----------

`py4pm` is a set of utilities to handle the xlsx output of EPA PMF5 software in python.

This project started because I needed to run several PMF for my PhD  and also needed to run some computation on these results.
The raw output of the EPA PMF5 software is a bit messy and hard to understand at a first
glance, and copy/pasting xlsx file is not my taste... So I ended developping this tools
for handling the tasks of maping the xlsx output to nice python objects, on which I can easily run
some computation.

Since I needed to plot the results afterward, I also added some plot utilities in this
package. It then has built in support for ploting :

 * chemical profile (both absolute and normalized)
 * species repartition among factor
 * timeserie contribution (*for all species* and profiles)
 * uncertainties plots (Bootstrap and DISP)
 * seasonal contribution

Moreover, this package include some of the DeltaTool utilities developed by the JRC, notably for profile chemical similarities between several PMF results.

