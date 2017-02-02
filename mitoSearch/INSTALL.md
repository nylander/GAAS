# Installation

## Dependencies
Python3

## Instructions
```
python3.5 -m virtualenv sbt_mash -p python3.5 --system-site-packages
. sbt_mash/bin/activate
pip install screed pytest PyYAML
pip install git+https://github.com/dib-lab/khmer.git
cd sbt_mash/
git clone https://github.com/dib-lab/sourmash.git -b sbt_search
cd sourmash && make install
```
## Notes
If Python is installed in a different location and not in your path use the following
```
export PATH="</path/to/python>/bin/:$PATH"
export LD_LIBRARY_PATH=</path/to/python>/lib/
```
