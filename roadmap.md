### Tasks for getting this plugin operational

1. Plugin offers a Hello, World GUI window
   1. See, as just one example: https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/communicating.html#communicating-with-the-user
1. Test plugin server workflow by releasing v0.01 and a minorly modified v.0.02
1. Plugin loads a vector layer from QGIS into plugin code using `pyqgis`
   1. See https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/vector.html and https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/vector.html#creating-vector-layers
1. Modifies the geometry of that layer somehow
1. Convert a qgis layer that the user specifies (via a GUI input widget) to a geopandas geodataframe object... _and back again._  These might help:
   1. Use this workflow plus a version edited to run 'in reverse'? https://gis.stackexchange.com/questions/362979/loading-geodataframe-as-qgis-vector-layer-without-exporting-to-shapefile
   1. Perhaps also of use: https://anitagraser.com/2018/11/18/movement-data-in-gis-16-towards-pure-python-trajectories-using-geopandas/
1. Makes basic connections to our library (via a submodule), importing it, tests out functionality (can equally be done before any of the more substantial pyqgis coding above)
1. Hardcoding all of the above to make the plugin, QGIS vector layer modification, and basic reprojections to work more or less together. Projections do not need to be user-inputted yet.
1. UI design charrette
1. UI implementation
1. Write projstring and WKT parsers (limited implementations of just our subset of WKT projections, not all possible WKT projections, of course! likewise, just handle our projstrings.)
   1. Overview of parsing in Python: https://tomassetti.me/parsing-in-python/
   1. Also perhaps of use: https://www.dabeaz.com/ply/ply.html#ply_nn0


### QGIS Plugin tutorials
A bunch of these. This one seems the most useful:
* https://www.qgistutorials.com/en/docs/3/building_a_python_plugin.html

These may also be helpful:
* https://www.e-education.psu.edu/geog489/l4_p9.html
* https://digital-geography.com/build-qgis-plugin/
* https://subscription.packtpub.com/book/application_development/9781783984985/1/ch01lvl1sec16/creating-a-qgis-plugin
* https://gis-ops.com/qgis-3-plugin-tutorial-plugin-development-explained-part-1/

### Managing dependencies:
* Packaging them by including wheels:
   * https://gis.stackexchange.com/a/349502
   * https://realpython.com/python-wheels/#python-packaging-made-better-an-intro-to-python-wheels
* One might want to set the minimum QGIS version for the plugin to be that point at which linux, macos, and windows versions were all bundling `scipy` and `numpy`, which we will need to determine.
  * `3.4.5 MacOS` does not allow `import scipy` from the Python console.
  * For people whose installations still don't have those dependencies and/or for people who want to try to get the plugin to work for earlier versions (for which they'd also have to manually modify a configuration file in the plugin, I assume?), we may want to eventually develop some pointers as to how to get `scipy` etc. installed on various platforms in a way that is visible to QGISes on those platforms, as some past instructions exist. See: https://gis.stackexchange.com/questions/366848/scp-plugin-for-qgis-3-10-on-mac

### Adding submodule
* https://github.blog/2016-02-01-working-with-submodules/
