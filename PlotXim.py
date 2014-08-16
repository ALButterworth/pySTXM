# Created 2014, Zack Gainsforth
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import io
import Image
import pandas as pd
import glob2

def SaveTifStack(FileName, Stack):
    try:
        # We can only write multiframe tiffs if the user has downloaded the tifffile module.
        from tifffile import imsave as tiffsave
    except:
        return

    # # tifffile will save it using the order: frames, height, width.  So we have to redimension the array.
    # # Currently the array is height, width, frames.
    # if len(Stack.shape) == 3:
    #     Stack = Stack.swapaxes(1,2)
    #     Stack = Stack.swapaxes(0,1)
    # In the case shape is 2D then it will just write it correctly.

    # Now save it.
    tiffsave(FileName, Stack)


def ReadXim(XimName):
    # ReadXim reads a Xim header, and then reads in the xim files to produce an image, map or stack.
    #
    # XimName is the name of the header.
    # Returns a dictionary: Xim with the following entries:
    # NumRegions: How many regions in the image.  Only supported for stacks.
    # Energies: How many energies.  1 means image.  2 means map, 3 or more is a stack.
    # Image: A single frame representation.  For an image, (1 energy) it is just the image. For the map,
    # it is the PostEdge - PreEdge.  For a stack, it is a sum along the E axis.
    # Plot: Image plotted with the info about it in the title.  This can be saved to a Tif or jpg.
    # RegionN: For each region this contains the stack as a numpy cube.

    # Extract the file name for the header.
    BaseName, _ = os.path.splitext(XimName)
    AxisName = os.path.basename(BaseName)
    hdrName = BaseName + '.hdr'

    # Open the header file for the xim.
    with open(hdrName, 'r') as f:
        hdr = f.read()

    # Find out how many regions there are.
    RegionsMatch = re.search('Regions = \(([0-9]*),', hdr)
    try:
        NumRegions = int(RegionsMatch.group(1))
    except:
        print 'Could not get number of regions from hdr file.'
        return None

    # Find out the energ(ies).  Look in the StackAxis section for the entry Points = (numbers);
    # Pull out numbers.
    EnergyMatch = re.search('StackAxis.*?Points = \((.*?)\);', hdr, re.S)

    # Convert the text numbers into energies.
    try:
        Energies = np.fromstring(EnergyMatch.group(1), dtype = float, sep = ', ')
    except:
        print 'Could not get energy axis for image/stack'
        return None

    # Now the header is too nice.  It tells us how many energies there are in the first entry.  Test it and then
    # discard it.
    if Energies[0] != len(Energies)-1:
        print "Hdr file corrupted.  Energy axis length doesn't match the first entry (which gives the length of the axis)."
        return None

    Energies = np.delete(Energies, 0)

    # Put the data into the XimHeader
    Xim = {}
    Xim['Energies'] = Energies
    Xim['NumRegions'] = NumRegions

    # While Axis probably does it, we don't support multiple regions for anything other than a stack.
    # We can add this in the future.
    if len(Energies) <= 2 and NumRegions > 1:
        print "Multi-region images and maps are not currently supported."
        return None

    # Three options based on how many energies we have:
    # 1: It's just an image.
    # 2: It's a map.
    # >2: It's a stack.
    if len(Energies) == 1:
        # For images, we just load the image.  The stack is rather dull, just one frame.
        Xim['Type'] = 'Image'
        Xim['Image'] = np.loadtxt(BaseName + '_a.xim', dtype='uint16')
        Xim['Plot'] = PlotImage(Xim['Image'], AxisName + ': %0.2f eV'%Xim['Energies'][0])
        Xim['Region1'] = Xim['Image']
    elif len(Energies) == 2:
        Xim['Type'] = 'Map'
        PreEdge = np.loadtxt(BaseName + '_a000.xim', dtype='uint16')
        PostEdge = np.loadtxt(BaseName + '_a001.xim', dtype='uint16')
        Xim['Region1'] = np.array([PreEdge, PostEdge])
        Xim['Image'], Xim['Plot'] = PlotMap(PreEdge, PostEdge, Xim['Energies'], AxisName)
    else:
        Xim['Type'] = 'Stack'

        # For single region stacks, the numbers go 000, 001, 002, ...
        if NumRegions == 1:
            NumberIncrement = 1
        # For multi region stacks, region 1 goes 0000, 0010, 0020, ... and region 2 is 0001, 0011, 0021, ...
        else:
            NumberIncrement = 10

        # Load each region into Xim.
        for n in range(1, NumRegions+1):
            # Name this one.
            RegionName = 'Region%d'%n
            # Load just the first frame so we can allocate the numpy array with the right dimensions.
            if NumRegions == 1:
                ExtensionStr = '_a000.xim'
            else:
                ExtensionStr = '_a%04d.xim'%(n-1)
            FirstFrame = np.loadtxt(BaseName + ExtensionStr, dtype='uint16')

            # Allocate the numpy array.
            ThisRegion = np.zeros((len(Xim['Energies']), FirstFrame.shape[0], FirstFrame.shape[1]), dtype='uint16')
            # Put in the one frame we already loaded.
            ThisRegion[0] = FirstFrame
            # Loop for each of the remaining frames.
            for i in range(1, len(Xim['Energies'])):
                # Read in one frame and store it.
                try:
                    # (i*NumberIncrement+(n-1)) because n-1 gives the 0 based region (region 1 goes 000, then 010, etc...)
                    # and i*NumberIncrement gives us 010, 020, etc.
                    if NumRegions == 1:
                        ExtensionStr = '_a%03d.xim'%(i*NumberIncrement+(n-1))
                    else:
                        ExtensionStr = '_a%04d.xim'%(i*NumberIncrement+(n-1))
                    #ThisRegion[i] = np.loadtxt(BaseName + '_a%03d.xim'%(i*NumberIncrement+(n-1)), dtype='uint16')
                    # We're doing this with pandas read_csv because it is so fast.
                    t = pd.read_csv(BaseName + ExtensionStr, sep='\t', header=None)
                    # But Tolek has a \t at the end of the line and Pandas reads in a last column of NaNs.  So ditch
                    # the last col.
                    ThisRegion[i] = t.values[:,:-1]
                except IOError:
                    # It is common that stacks will be truncated.  In this case, we will fail at reading some file.
                    # Lets trunctate the stack here.
                    Xim['Energies'] = Xim['Energies'][:i]
                    ThisRegion = ThisRegion[:i]
                    # Done.  Go to next region.
                    break

            # Store it.
            Xim[RegionName] = ThisRegion

            # And generate "pretty images" i.e. thumbnails.
            ThisAxisName = AxisName
            if NumRegions > 1:
                ThisAxisName = AxisName + '_Region%d'%n
            Xim['Image'+RegionName], Xim['Plot'+RegionName] = PlotStack(Xim[RegionName], Xim['Energies'], AxisName)

    return Xim

def PlotImage(ImRaw, title):
    # Get the min and max intensity.
    vmin = np.min(ImRaw)
    vmax = np.max(ImRaw)
    # Plot it.
    ax = plt.imshow(ImRaw, interpolation='none', cmap='gray', origin='upper', vmin=vmin, vmax=vmax)
    plt.title(title)
    plt.colorbar()
    # And save that plot as an image in memory.
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    im = Image.open(buf)

    return im

def PlotMap(PreEdge, PostEdge, Energies, AxisName):
    # Compute the OD image.
    Map = -np.log(PostEdge.astype(float)/PreEdge.astype(float))

    # Make a plot in memory with a title.  Png format will be default.
    im = PlotImage(Map, AxisName + ': ln(%0.2f/%0.2f) eV'%(Energies[1], Energies[0]))

    # And pass it back out.
    return Map, im

def PlotStack(Region, Energies, AxisName):
    # Compute the sum image
    Sum = np.sum(Region.astype(float), axis=0)

    # Make a plot in memory with a title.  Png format will be default.
    im = PlotImage(Sum, AxisName + ': sum between %0.2f to %0.2f eV'%(Energies[0], Energies[-1]))

    # And pass it back out.
    return Sum, im



def PlotXim(XimName):
    # Extract the file name.
    BaseName, _ = os.path.splitext(XimName)

    Xim = ReadXim(XimName)


    for n in range(1, Xim['NumRegions']+1):

        RegionName = 'Region%d'%n
        # If there is more than one region, then we need to append a region name to the base name.
        if Xim['NumRegions'] > 1:
            RegionSepName = '_Region%d'%n  # In the fileName we want a separator
        else:
            RegionSepName = ''

        if Xim['Type'] == 'Stack':
            PlotName = 'Plot'+RegionName
        else:
            PlotName = 'Plot'

        # Save a user friendly plot image.
        Xim[PlotName].save(BaseName + RegionSepName + '.png')
        # Save the raw data as tif.
        SaveTifStack(BaseName + RegionSepName + '.tif', Xim[RegionName])

        # Save the energy axis (but not for single energy images).
        if len(Xim['Energies']) > 1:
            np.savetxt(BaseName + RegionSepName + '_energies.txt', Xim['Energies'])

def ProcessAllXimUnderDir(DirName=None):
    import gc

    if DirName is None:
        DirName =  os.getcwd()

    RawFiles = glob2.glob(DirName + '/**/*.hdr')
    for f in RawFiles:
        print 'Processing %s\n' % f
        try:
            PlotXim(f)
        except:
            # If that one Xim failed, go on and process the next.
            pass
        # We are tearing through RAM with these stacks.  Sometimes, the garbage collection can't keep up and we run out of memory.
        # Or we force it to clean up after each stack.
        gc.collect()
        print '\n'


if __name__ == '__main__':
    #PlotXim('/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6,
    # S3/131121/11_131121000.hdr')
    #PlotXim('/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6,
    # S3/131121/11_131121006/11_131121006.hdr')
    # PlotXim('/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6, '
    # 'S3/131121/11_131121013/11_131121013.hdr')
    #PlotXim('/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6,
    # S3/131121/11_131121020/11_131121020.hdr')
    # PlotXim('/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6, '
    # 'S3/131121/11_131121044/11_131121044.hdr')

    ProcessAllXimUnderDir('/Users/Stardust/Zack/MoonBaseBeta/ParticleBase/Track 162a - Cecil/Grid A6/20131121 - STXM - Cecil A6, S3/131121')
    print 'Done'

'''

			QAxis = { Name = "Sample Y"; Unit = "um"; Min = 496.7525; Max = 497.1525; Dir = 1;
				Points = (10, 496.7525, 496.7925, 496.8325, 496.8725, 496.9125, 496.9525, 496.9925, 497.0325, 497.0725, 497.1125);
};
});

	StackAxis = { Name = "Energy"; Unit = "eV"; Min = 445; Max = 485; Dir = 1;
		Points = (122, 445.0000, 445.5000, 446.0000, 446.5000, 447.0000, 447.5000, 448.0000, 448.5000, 449.0000, 449.5000, 450.0000, 450.5000, 451.0000, 451.5000, 452.0000, 452.5000, 453.0000, 453.5000, 454.0000, 454.5000, 455.0000, 455.1010, 455.2020, 455.3030, 455.4040, 455.5050, 455.6060, 455.7070, 455.8080, 455.9090, 456.0100, 456.1110, 456.2120, 456.3130, 456.4140, 456.5150, 456.6160, 456.7170, 456.8180, 456.9190, 457.0200, 457.1210, 457.2220, 457.3230, 457.4240, 457.5250, 457.6260, 457.7270, 457.8280, 457.9290, 458.0300, 458.1310, 458.2320, 458.3330, 458.4340, 458.5350, 458.6360, 458.7370, 458.8380, 458.9390, 459.0400, 459.1410, 459.2420, 459.3430, 459.4440, 459.5450, 459.6460, 459.7470, 459.8480, 459.9490, 460.0500, 460.1510, 460.2520, 460.3530, 460.4540, 460.5550, 460.6560, 460.7570, 460.8580, 460.9590, 461.0600, 461.1610, 461.2620, 461.3630, 461.4640, 461.5650, 461.6660, 461.7670, 461.8680, 461.9690, 462.0700, 462.5000, 463.0000, 463.5000, 464.0000, 464.5000, 465.0000, 465.5000, 466.0000, 466.5000, 467.0000, 467.5000, 468.0000, 468.5000, 469.0000, 469.5000, 470.0000, 471.0000, 472.0000, 473.0000, 474.0000, 475.0000, 476.0000, 477.0000, 478.0000, 479.0000, 480.0000, 481.0000, 482.0000, 483.0000, 484.0000, 485.0000);
};

	Channels =


	'''