import numpy as np
import io
def import_airfoil(filename):
    # Load the data from the text file
    flnm = str(filename) + ".dat"
    s = io.BytesIO(open(flnm, 'rb').read().replace(b'-', b' -'))
    dataBuffer = np.genfromtxt(s, delimiter='  ', skip_header=1)

    # Extract data from the loaded dataBuffer array
    dataX = dataBuffer[:, 0]
    dataY = dataBuffer[:, 1]

    return dataX, dataY