# The BDNYC Data Archive

## Installation

To install, just do:

```pip install BDNYCdb```

Then download the `bdnyc198.db` database file [here](https://s3.amazonaws.com/bdnyc/bdnyc198.db). This initial release contains the astrometry, photometry and spectra for the 198 objects in the [Filippazzo et al. (2015)](http://adslabs.org/adsabs/abs/2015ApJ...810..158F/) sample.

## Using the database

To start using the database, launch iPython, import the module, then initialize the database.

```
from BDNYCdb import BDdb
db = BDdb.get_db('/path/to/bdnyc198.db')
```

## Enjoy!
