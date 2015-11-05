#!/usr/bin/python

# David Rodriguez
# Created: November 4, 2015
# Read in a table from an SQL query of the BDNYC Database and convert to a VOTable file for use elsewhere

from astropy.table import Table, Column
from astropy.io.votable import from_table


def table_add(tab, data, col):
    """
    Function to parse dictionary list **data** and add the data to table **tab** for column **col**

    Parameters
    ----------
    tab: Table class
      Table to store values
    data: list
      Dictionary list from the SQL query
    col: str
      Column name (ie, dictionary key) for the column to add

    Returns
    -------
    None

    """

    x = []
    for i in range(len(data)):
        temp = data[i][col]

        # Fix up None elements
        if temp is None: temp = ''

        x.append(temp)

    print 'Adding column', col
    tab.add_column(Column(x, name=col))


def dict_tovot(tabdata, tabname='votable.xml'):
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

    # Check if input is a dictionary
    if not isinstance(tabdata[0], dict):
        raise TypeError('Table must be a dictionary. Call the SQL query with query_dict.execute()')

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
