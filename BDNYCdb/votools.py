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

        # If the particular key is not present, create an empty value (used for photometry tables)
        if col not in data[i]:
            temp = ''
        else:
            temp = data[i][col]

        # Fix up None elements
        if temp is None:
            temp = ''

        x.append(temp)

    print 'Adding column', col
    tab.add_column(Column(x, name=col))


def dict_tovot(tabdata, tabname='votable.xml', phot=False):
    """
    Converts dictionary table **tabdata** to a VOTable with name **tabname**

    Parameters
    ----------
    tabdata: list
      SQL query dictionary list from running query_dict.execute()
    tabname: str
      The name of the VOTable to be created
    phot: bool
      Parameter specifying if the table contains photometry to be merged

    Returns
    -------
    None

    """

    # Check if input is a dictionary
    if not isinstance(tabdata[0], dict):
        raise TypeError('Table must be a dictionary. Call the SQL query with query_dict.execute()')

    # Create an empty table to store the data
    t = Table()

    colnames = tabdata[0].keys()

    # If this is a photometry table, parse it and make sure to have the full list of columns
    if phot:
        tabdata = photparse(tabdata)

        for i in range(len(tabdata)):
            tmpcol = tabdata[i].keys()
            for elem in tmpcol:
                if elem not in colnames:
                    colnames.append(elem)

    # Run through all the columns and create them
    for elem in colnames:
        table_add(t, tabdata, elem)

    # Output to a file
    votable = from_table(t)
    votable.set_all_tables_format('binary')
    votable.to_xml(tabname)

    print 'Table created:', tabname


def photaddline(tab, id):
    """
    Loop through the dictionary list **tab** creating a line for the source specified in **id**

    :param tab:
      Dictionary list of all the photometry data
    :param id:
      ID of source in the photometry table (source_id)
    :return:
      Dictionary with all the data for the specified source
    """

    colnames = tab[0].keys()
    tmpdict = dict()
    for i in range(len(tab)):

        # If not working on the same source, continue
        if tab[i]['source_id'] != id:
            continue

        for elem in colnames:
            if elem not in ['comments','epoch','instrument_id','magnitude','magnitude_unc','publication_id','system','telescope_id']:
                tmpdict[elem] = tab[i][elem]
            elif elem == 'band':
                continue
            else:
                tmpstr = tab[i]['band']+'.'+elem
                tmpdict[tmpstr] = tab[i][elem]

    return tmpdict


def photparse(tab):
    """
    Parse through a photometry table to group by source_id

    Parameters
    ----------
    tab: list
      SQL query dictionary list from running query_dict.execute()

    Returns
    -------
    newtab: list
      Dictionary list after parsing to group together sources

    """

    # Loop through the table and grab unique band names and source IDs
    bandnames = []
    uniqueid = []
    for i in range(len(tab)):
        tmp = tab[i]['band']
        tmpid = tab[i]['source_id']

        if tmp not in bandnames:
            bandnames.append(tmp)

        if tmpid not in uniqueid:
            uniqueid.append(tmpid)

    # Loop over unique id and create a new table for each element in it
    newtab = []
    for id in uniqueid:
        tmpdict = photaddline(tab, id)
        newtab.append(tmpdict)

    return newtab



