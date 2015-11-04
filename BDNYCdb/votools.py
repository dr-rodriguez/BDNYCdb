#!/usr/bin/python

# David Rodriguez
# Created: November 4, 2015
# Read in a table from an SQL query of the BDNYC Database and convert to a VOTable file for use elsewhere

from BDNYCdb import BDdb
from astropy.table import Table, Column
from astropy.io.votable import from_table

# Function to parse the dictionary and add to the table
def table_add(tab, dat, col):
    x = []
    for i in range(len(dat)):
        temp = dat[i][col]

        # Fix up None elements
        if temp is None: temp = ''

        x.append(temp)

    print 'Adding column', col
    tab.add_column(Column(x, name=col)

# Convert dictionary list to VOTable
def dict_tovot(tabdata, tabname):
    """
    Converts dictionary table **tabdata** to a VOTable with name **tabname**

    Parameters
    ----------
    tabdata: list
      SQL query dictionary list from running query_dict.execute()
    tabname: str
      The name of the VOTable to be created

    Returns
    -------
    None

    """

    # Create an empty table to store the data
    t = Table()

    # Run through all the columns and create them
    colnames = tabdata[0].keys()
    for elem in colnames: table_add(t, tabdata, elem)

    # Output to a file
    votable = from_table(t)
    votable.set_all_tables_format('binary')
    votable.to_xml(tabname)

    print 'Table created:', tabname