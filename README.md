# Brown Dwarf NYC Database



## Installation

First you'll need [anaconda](https://store.continuum.io/cshop/anaconda/), which is a high-performance python implementation that has a lot of libraries for scientific computing already setup for you.  You can install it several ways:

- Using [pyenv](https://github.com/yyuu/pyenv): `pyenv install anaconda-2.2.0`
- Using the Continuum installer from: http://continuum.io/downloads

Note that you'll need to install `anaconda`, and **not** `anaconda3`.

Next you'll need to download and setup [Astrolib PySynphot](http://ssb.stsci.edu/pysynphot/docs/).  I recommend creating a `pysyn_cdbs` directory in the same folder as this repo for ease of use and isolation from other programs that require those files:

```sh
cd ~/code/BDNYCdb
mkdir -p pysyn_cdbs/extinction # creating extinction directory supresses warning
curl --compressed ftp://ftp.stsci.edu/cdbs/tarfiles/synphot1.tar.gz | tar xfz - -C pysyn_cdbs
curl --compressed ftp://ftp.stsci.edu/cdbs/tarfiles/synphot2.tar.gz | tar xfz - -C pysyn_cdbs
curl --compressed ftp://ftp.stsci.edu/cdbs/tarfiles/synphot3.tar.gz | tar xfz - -C pysyn_cdbs
curl --compressed ftp://ftp.stsci.edu/cdbs/tarfiles/synphot4.tar.gz | tar xfz - -C pysyn_cdbs
curl --compressed ftp://ftp.stsci.edu/cdbs/tarfiles/synphot5.tar.gz | tar xfz - -C pysyn_cdbs
curl --compressed ftp://ftp.stsci.edu/cdbs/tarfiles/synphot6.tar.gz | tar xfz - -C pysyn_cdbs
```

Now download the `BDNYC.db` database file from the Dropbox link.  I recommend keeping it in the same location as this repo.  Please note that you'll need to get someone from BDNYC to give you access to that file.

## Running the [ipython](http://ipython.org/) HOWTO

You'll need to setup some environment variables so the HOWTO knows where to look for your BDNYC and PySynphot files.

- You can pass in `PYSYN_CDBS` and `BDDB_FILE`: `PYSYN_CDBS=pysyn_cdbs BDDB_FILE=BDNYC.db ipython HOWTO.ipynb`
- You can store them in a `.env` file and use something like [dotenv](https://github.com/bkeepers/dotenv)

```sh
gem install dotenv # rbenv rehash if needed
echo "export PYSYN_CDBS=$(pwd)/pysyn_cdbs
export BDDB_FILE=$(pwd)/BDNYC.db" > .env
dotenv ipython notebook HOWTO.ipynb
```

- You can store them in a `.env` file and just `source` it (to load them into your shell's `ENV`) after you `cd` into the directory:

```sh
echo "export PYSYN_CDBS=$(pwd)/pysyn_cdbs
export BDDB_FILE=$(pwd)/BDNYC.db" > .env
source .env
ipython notebook HOWTO.ipynb
```

- You can add them to your `~/.bash_profile`, `~/.bashrc`, or `~/.zshrc` and open a new shell:

```sh
echo "export PYSYN_CDBS=$(pwd)/pysyn_cdbs
export BDDB_FILE=$(pwd)/BDNYC.db" >> ~/.bash_profile
# open a new shell or source ~/.bash_profile
ipython notebook HOWTO.ipynb
```


## Done!
