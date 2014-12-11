__author__ = 'svfilhol'


# Example:
# os.system("qgis-colorramp.py --vmin 0 --vmax 0.18 --extend -1 5 /home/svfilhol/Downloads/nrwc_Sim.cpt")

def qgis_colorramp(Zmin, Zmax, Extmin, Extmax, cptFile):
    """
    Function to produce colorscale based on .cpt file  available at (http://soliton.vm.bytemark.co.uk/pub/cpt-city/)
     Return product can be used in QGIS
    :param Zmin: Minimum of the colorscale
    :param Zmax: Maximum of the colorscale
    :param Extmin: Far extent of colorscale (to include features far beyond the visible range of the colorscale
    :param Extmax: Far extent of colorscale (to include features far beyond the visible range of the colorscale
    :param cptFile: .cpt file directory
    :return: Save .txt, .png file for import in QGIS
    """
    import os
    mycommand = "qgis-colorramp.py --vmin " + str(Zmin) + " --vmax " + str(Zmax) + " --extend " + str(Extmin) + " " + str(Extmax) + " " + cptFile
    os.system(mycommand)

