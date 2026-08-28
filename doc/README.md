# Dempy Documentation
Documentation of the Dempy toolbox. 
https://dempy.readthedocs.io

*S. Filhol, H. Merzisen, L. Barrois, August 2026

## Documentation TODO:
- [ ] install documentation
- [ ] setup readthedocs
- [ ] add refs.bib file


## Contributing to the Documentation

The documentation is written in Markdown and compiled with the tool https://www.mkdocs.org. Edit the markdown files located in the folder `dempy/doc/docs/`. If you add a page, make sure to add its path to the file [mkdocs.yml](https://github.com/ArcticSnow/dempy/blob/main/doc/mkdocs.yml). The online *readthedocs* website will build automatically when changes are being commited. 

## Build the doc on your local machine

1. Create a local Python VE and install the following packages to build the doc locally

   ```bash
   pip install mkdocs sphinx_rtd_theme mkdocs-material pygments
   
   # cloning using SSH key
   git clone git@github.com:ArcticSnow/dempy.git
   ```

2. Modify a given file or create a new one. If you create a new one add it to `mkdocs.yml`

3. run `mkdocs serve` to visualize the modification on your local machine

4. open the URL: http://127.0.0.1:8000/ in your browser

For customization of the theme, please use the documentation https://squidfunk.github.io/mkdocs-material/. All custom parameters are in the file `mkdocs.yml`

