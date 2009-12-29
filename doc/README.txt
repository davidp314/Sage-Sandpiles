I. Automated testing
1. cp sandpile.sage sandpile.py
2. add "from sage.all import *" to sandpile.py as the first line
3. sage: -t sandpile.py

II. Documentation
1. Edit sage_sandpiles.rst.  
2. Use utilities.sage to help generate the documentation (by listing the
methods, creating links, etc.
3. **Read utilities.sage file.***
NOTE: cp latest version of sandpile.sage to doc/ directory.
Sphinx stuff:

* create new doc directory
* cd into new doc directory and run sphinx-quickstart
-- instead of "index", use the name "sandpile"
* diff conf.py with the old conf.py
-- default_role = 'math'
-- latex_preamble = '\usepackage{amssymb}'
-- extensions = ['sphinx.ext.pngmath']

Then:

make html
