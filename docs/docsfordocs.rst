1) index.rst : main file / Home page of readthedocs.
    Must cointain all the other .rst files in the toctrees!!!
    eg. : .. toctree::
             :hidden:
             :maxdepth: 6
             :caption: tutorials based on gui
   
             gui (file)
             rvs (file)
             transit (file)
             rvtran (file)
             otbfpu (file)
             stability (file)

 Each toctree is a new section on the Sphinx Template.

2) All the .rst files contained in a toctree must start with
   eg. for gui.rst file 
       .. _gui:
       (.. _filename:)
   in order to be linked with the toctree!!!

3) The headline of a document must be underlined with ... !!!
   eg. in gui.rst file
       GUI Layout
       ..........

    and each sub-header must be underlined with --- or === !!
    eg. in gui.rst file
        Data area
        ---------

4) Insert a picture/gif file
   eg. in gui.rst file
       .. figure:: /images/homepage.png

          (comments)

5) Bold words with the symbol **
   eg. **Exostriker** 

6) Emphasis words with the symbol *
   eg. *Exostriker*

7) Superscripts
   eg. χ\ :sup:`2`: chi-squared

8) Subscripts
   eg. a\ :sub:`pl`\ /(divided) R\ :sub:`*`\ : planet semimajor axis in units of stellar radius.

9) Add extra lines to seperate paragraphs with the symbol --- 
   eg. in gui.rst file
       ----------------------------------------------------------------------------------------------------------

10) Add hyperlinks 
    eg. in gui.rst file
        For more information visit `pyqtgraph documentation`_.

        .. _pyqtgraph documentation : https://pyqtgraph.readthedocs.io/en/latest/index.html

11) Bullet point with the symbol *
    eg. in gui.rst file
       * Stdout/Stderr


For more info check links : 1) https://docutils.sourceforge.io/docs/ref/rst/roles.html#customization
2) https://docutils.sourceforge.io/docs/ref/doctree.html#element-hierarchy
3) https://docutils.sourceforge.io/docs/ref/rst/directives.html#admonitions
4) https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#indentation