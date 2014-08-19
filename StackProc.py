# Created 2014, Zack Gainsforth
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import io
import Image
import tifffile as tf
from QuickPlot import QuickPlot

def ReadTifStack(StackFileName, EnergyFileName=None):

    # Extract the file name sans extention.
    BaseName, _ = os.path.splitext(StackFileName)

    # If the user didn't supply a file name for the energies, then we assume the default name.
    if EnergyFileName is None:
        EnergyFileName = BaseName + '_energies.txt'

    # Read the energy axis.
    E = np.genfromtxt(EnergyFileName, dtype=np.float32)

    # Read the stack.
    s = tf.imread(StackFileName)

    return s, E

def SaveStackEnergies(Stack, Energies, EnergiesToSave, PreEdgeEnergy=None, PostEdgeEnergy=None, Path=''):

    # Save all the energies in EnergiesToSave.
    for e in EnergiesToSave:
        tf.imsave(Path + str(e) + ' eV.tif', Stack[Energies.searchsorted(e)])

        # If the user also supplied a pre-edge energy, then we will do each energy minus pre-edge.
        if PreEdgeEnergy is not None:
            tf.imsave(Path + str(e) + ' - ' + str(PreEdgeEnergy) + '.tif', Stack[Energies.searchsorted(e)] - Stack[Energies.searchsorted(PreEdgeEnergy)])

        # If the user also supplied a post-edge energy (edge jump), then produce images that are normalized by the
        # edge jump.  Because this can be noisy, we add one to the edge jump and to the energy itself.
        if PreEdgeEnergy is not None and PostEdgeEnergy is not None:
            EdgeJump = Stack[Energies.searchsorted(PostEdgeEnergy)] - Stack[Energies.searchsorted(PreEdgeEnergy)] + 1
            tf.imsave(Path + str(e) + ' div (' + str(PostEdgeEnergy) + '-' + str(PreEdgeEnergy) + ').tif',
                      (Stack[Energies.searchsorted(e)]+1) / EdgeJump)

def GetMeanSpectrum(Stack):
    return np.mean(np.mean(Stack, axis=2), axis=1)

def GetSumSpectrum(Stack):
    return np.sum(np.sum(Stack, axis=2), axis=1)

def ClickOnPlot(fig):

    #fig, ax = QuickPlot(E, Spec, xlabel='eV', title='Click with mouse to save values.')
    Coords = []

    def onclick(event):
        #print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
        print 'x=%d, y=%d, xdata=%f, ydata=%f'%(event.x, event.y, event.xdata, event.ydata)
        Coords.append((event.xdata, event.ydata))

    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

    return Coords

def ChooseEnergiesSaveAndPlot(S, E):

    print('Drawing sum spectrum for stack.')
    print('First click on the pre-edge. Then click on any other energies.  Click last on the post-edge for edge jump.')

    # Make a sum spectrum.
    Spec = GetSumSpectrum(S)

    #Make a plot
    fig, ax = QuickPlot(E, Spec, xlabel='eV', ylabel='O.D.', title='Click pre-edge, energies, post-edge', boldlevel=1)

    # Get the energies we want to plot.
    Coords = ClickOnPlot(fig)
    EnergiesToSave, _ = zip(*Coords)

    # Save the requested images
    SaveStackEnergies(S, E, EnergiesToSave, PreEdgeEnergy=EnergiesToSave[0], PostEdgeEnergy=EnergiesToSave[-1])

    # The user chose coords.  Now draw the plot with vertical lines on his energies.
    fig, ax = QuickPlot(E, Spec, xlabel='eV', ylabel='O.D.', title='Sum Spectrum', boldlevel=3)
    plt.savefig('Energies Saved No Lines.png', format='png')
    for c in Coords:
        ax.plot([c[0], c[0]], [np.min(Spec), np.max(Spec)], linewidth=3)

    plt.savefig('Energies Saved.png', format='png')

    plt.show()

if __name__ == '__main__':
    StackName = '/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6, S3/11_131121013_Cecil_a6s3_O/Stack aligned OD.tif'
    EnergiesName = '/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6, S3/11_131121013_Cecil_a6s3_O/Stack aligned_energies.txt'
    BaseDir = '/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6, S3/11_131121013_Cecil_a6s3_O'

    cwd = os.getcwd()
    os.chdir(BaseDir)

    S, E = ReadTifStack(StackName, EnergiesName)
    ChooseEnergiesSaveAndPlot(S,E)

    os.chdir(cwd)

    # # Make a sum spectrum.
    # Spec = GetSumSpectrum(S)
    # # Get the energies we want to plot.
    # ClickOnPlot(E, Spec)
