#!/usr/bin/python
# BDNYC database
import sqlite3 as sql, numpy as np, matplotlib.pyplot as plt

def get_db(dbpath, rows=False, modify=False):
  '''
  Creates a conection and cursor object with SQL database at *dbpath*
  '''
  con = sql.connect(dbpath, isolation_level=None)
  con.text_factory = str
  if rows: con.row_factory = sql.Row
  return [con, con.cursor()] if modify else con.cursor()

def spec_arrays(data, rows=False, SNR=False, plot=False):
  '''
  Converts database *data* from list of strings to numpy arrays. If *plot* it plots the spectrum.
  '''
  w, f, e = [np.array(map(float,i.split())) for i in [data['wavelength'], data['flux'], data['unc']]+([data['SNR']] if SNR else [])] if rows else [np.array(map(float,i.split())) for i in list(data)]
  if plot:
    plt.figure(), plt.loglog(w, f, c='b'), plt.fill_between(w, f-e, f+e, color='k', alpha=0.2), plt.xlim(0.4,3.0), plt.grid(True), plt.yscale('log',nonposy='clip')
    if rows: plt.figtext(0.15,0.88,'{}\n{}\n{}\n{}'.format(data['filename'],data['telescope'],data['instrument'],data['obs_date']),verticalalignment='top')
  return [w,f,e]

def inventory(dbpath, ID='', plot=True):
  '''
  Prints a summary of all objects in the database at *dbpath*. If *ID* prints only that object's summary and plots if *plot*
  '''
  db = get_db(dbpath, rows=True)
  if ID: D = db.execute("SELECT id, unum, shortname, ra, dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='optical'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR') FROM sources WHERE id=?", [ID]).fetchall()
  else: D = db.execute("SELECT id, unum, shortname, ra, dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='optical'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR') FROM sources").fetchall()
  if D:
    printer(['id','unum','shortname','ra','dec','Spec Count','Phot Count','Parallax','SpT (opt)','SpT (IR)'], D)
    if ID and plot: [spec_arrays(i, rows=True, plot=True) for i in db.execute("SELECT * FROM spectra WHERE source_id=?", [ID]).fetchall()]
  else: print "No sources found."
  
def printer(labels, values, format='', to_txt=None):
  '''
  Prints a nice table of *values* with *labels* with auto widths. Save *to_txt* filepath, e.g. to_txt='/Users/Joe/Desktop/printout.txt' 
  '''
  print '\r'
  values = [["None" if not i else "{:.10g}".format(i) if isinstance(i,(float,int)) else i if isinstance(i,(str,unicode)) else "{:.10g} {}".format( float(i.magnitude if hasattr(i,'magnitude') else i), str(i.units if hasattr(i,'units') else '').split()[1]) for i in j] for j in values]
  auto, txtFile = [max([len(i) for i in j])+2 for j in zip(labels,*values)], open(to_txt, 'a') if to_txt else None
  lengths = format if isinstance(format,list) else auto
  col_len = [max(auto) for i in lengths] if format=='max' else [150/len(labels) for i in lengths] if format=='fill' else lengths
  for l,m in zip(labels,col_len):
    print str(l).ljust(m),
    if to_txt: txtFile.write(str(l).replace(' ','').ljust(m))
  for v in values:
    print '\n',
    if to_txt: txtFile.write('\n') 
    for k,j in zip(v,col_len):
      print str(k).ljust(j),
      if to_txt: txtFile.write(str(k).replace(' ','').ljust(j))
  print '\n'