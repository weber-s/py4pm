import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def add_season(df, month=True, month_to_season=None):
    """
    Add a season column to the DataFrame df.

    parameters
    ----------

    df: Pandas DataFrame.
        The DataFrame to work with.
    month : add month number, default True

    return
    ------

    dfnew: a new pandas DataFrame with a 'season' columns.

    """

    if month_to_season is None:
        month_to_season = {1:'Winter', 2:'Winter', 3:'Spring', 4:'Spring', 5:'Spring', 6:'Summer',
                           7:'Summer', 8:'Summer', 9:'Fall', 10:'Fall', 11:'Fall',
                           12:'Winter'}

    dfnew = df.copy()

    # ensure we have date in index
    if isinstance(dfnew.index, pd.DatetimeIndex):
        dfnew["Date"] = dfnew.index
        dropDate = True
    elif 'Date' in dfnew.columns:
        dfnew["Date"] = pd.to_datetime(dfnew["Date"])
        dropDate = False
    else:
        print("No date given")
        return

    # add a new column with the number of the month (Jan=1, etc)
    dfnew["month"] = dfnew["Date"].apply(lambda x: x.month)
    # sort it. This is not mandatory.
    dfnew.sort_values(by="month", inplace=True)

    # add the season based on the month number
    dfnew["season"] = dfnew["month"].replace(month_to_season)

    if not month:
        dfnew.drop(columns=["month"], inplace=True)
    if dropDate:
        dfnew.drop(columns=["Date"], inplace=True)

    # and return the new dataframe
    return dfnew


def format_xaxis_timeseries(ax):
    """Format the x-axis timeseries with minortick = month and majortick=year

    :ax: the ax to format

    """
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("\n%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))

    plt.setp(ax.xaxis.get_majorticklabels(), rotation=00, ha="center")
