import glob

from flagships.post_processing.EquilibriumPlotting import *
from flagships.post_processing.ParseFlagshipsFile import FlagshipsParser

equil_dir = 'out'
output_dir = equil_dir + '_plots'
for equil_file in glob.glob(os.path.join(equil_dir, '*.hdf5')): 
    f = FlagshipsParser('', equil_file)
    PlotBoundary(f)
    PlotPsi10PctContours(f)
    basename = equil_file.split(os.sep)[-1]
    plt.title(basename)
    plt.gca().set_aspect('equal')
    plt.savefig(os.path.join(output_dir, basename[:-5]+'.png'))
    plt.clf()
