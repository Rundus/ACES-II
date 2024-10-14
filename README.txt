--- README ---

    UIOWA_CDF_operator is a set of python scripts that process .tad telemetry files exported from an ALTAIR
    tm device to a useful scientific data product in the command data format (CDF) file structure.


--- CODE REQUIREMENTS ---

    (1) Anaconda python interpreter (preferred, not tested if required)
    (2) the .../Data folder must be in the same directory as the python scripts and follow the naming structure
    (3) "model" data files must be placed in the right directory under the correct /Data/TRICEII/... folder in order to produce the ACESII cdf files. The code uses these TRICE files to skip tedious steps in writing out the ACESII files.
    (4) the anaconda interpreter has most of the needed python libraries except for:
        (i) cdflib
        (ii) pycdf
        (iii) tqdm

    NOTE: pycdf requires not only a pip install, but also the NASA CDF library for the specific operating system. These
    can be found here:

        pycdf: https://spacepy.github.io/

        NASA CDF: https://cdf.gsfc.nasa.gov/html/sw_and_docs.html


--- CONTACT ---

    Author: Connor Feltman
    Author email: cfeltman@uiowa.edu


--- GIT ---

    the git repository (without most of the rocket data) for this project is found here:

    https://github.com/Rundus/UIOWA_CDF_operator.git


--- DATA FILE LOCATION ---

    the data files for the ACESII mission are stored on the UIOWA phi/tau server, located here:

    https://phi.physics.uiowa.edu/science/tau/data0/rocket/

