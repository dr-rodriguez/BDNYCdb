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
    # If this is a photometry table, make sure to have the full list of columns
    if phot:
        for i in range(len(tabdata)):
            tmpcol = tabdata[i].keys()
            for elem in tmpcol: if elem not in colnames: colnames.append(elem)


    # Run through all the columns and create them
    for elem in colnames:
        table_add(t, tabdata, elem)

    # Output to a file
    votable = from_table(t)
    votable.set_all_tables_format('binary')
    votable.to_xml(tabname)

    print 'Table created:', tabname


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

    # Create new table, copying over existing columns but adding new ones with the band names
    newtab = []
    colnames = tab[0].keys()
    id0 = tab[0]['source_id']
    for i in range(len(tab)):
        tmpdict = dict()

        # If not working on the same source, save the line and clear out the dictionary
        if tab[i]['source_id'] != id0:
            print tab[i]['source_id']
            id0 = tab[i]['source_id']
            newtab.append(tmpdict)
            tmpdict.clear()

        for elem in colnames:
            if elem not in ['comments','epoch','instrument_id','magnitude','magnitude_unc','publication_id','system','telescope_id']:
                tmpdict[elem] = tab[i][elem]
            elif elem == 'band':
                continue
            else:
                tmpstr = tab[i]['band']+'.'+elem
                print tmpstr
                tmpdict[tmpstr] = tab[i][elem]



