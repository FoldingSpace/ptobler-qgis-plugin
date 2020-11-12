### Tasks for getting this plugin operational

1. Plugin offers a Hello, World GUI window
   1. See, as just one example: https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/communicating.html#communicating-with-the-user 
1. Test plugin server workflow by releasing v0.01 and a minorly modified v.0.02
1. Plugin loads a vector layer from QGIS into plugin code using `pyqgis`
   1. See https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/vector.html and https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/vector.html#creating-vector-layers
1. Modifies the geometry of that layer somehow
1. Makes basic connections to our library (via a submodule), importing it, tests out functionality (can equally be done before any of the more substantial pyqgis coding above)
1. Hardcoding all of the above to make the plugin, QGIS vector layer modification, and basic reprojections to work more or less together. Projections do not need to be user-inputted yet.
1. UI design charrette
1. UI implementation
1. Write projstring and WKT parsers (limited implementations of just our subset of WKT projections, not all possible WKT projections, of course! likewise, just handle our projstrings.)
   1. Overview of parsing in Python: https://tomassetti.me/parsing-in-python/
   1. Also perhaps of use: https://www.dabeaz.com/ply/ply.html#ply_nn0 


### QGIS Plugin tutorials
* https://www.e-education.psu.edu/geog489/l4_p9.html 
* https://www.qgistutorials.com/en/docs/3/building_a_python_plugin.html 
* https://digital-geography.com/build-qgis-plugin/
* https://subscription.packtpub.com/book/application_development/9781783984985/1/ch01lvl1sec16/creating-a-qgis-plugin 
* https://gis-ops.com/qgis-3-plugin-tutorial-plugin-development-explained-part-1/ 

### Managing dependencies via packaging wheels
* https://gis.stackexchange.com/a/349502 
* https://realpython.com/python-wheels/#python-packaging-made-better-an-intro-to-python-wheels 

### Adding submodule
* https://github.blog/2016-02-01-working-with-submodules/ 
