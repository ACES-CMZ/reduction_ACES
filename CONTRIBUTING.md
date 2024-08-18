Contributing to the ACES software repository
============================================

All code is run through `flake8` code style checks before being merged.  This is primarily to ensure that code is syntactically valid, but it also helps keep the data uniform.  

The check that gets run is `flake8 aces --count --max-line-length=250 --ignore=W504,E265,F401,W292,E261,E262,E124,W503,E266` (https://github.com/ACES-CMZ/reduction_ACES/blob/main/tox.ini#L100).  That may be updated in the future, so be sure to check tox.ini if you're getting odd results.
