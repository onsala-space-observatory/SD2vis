import pylab as pl
import numpy as np
import random as rn
from clean_cli import clean_cli as clean
from clearcal_cli import clearcal_cli as clearcal
from importfits_cli import importfits_cli as importfits
from taskinit import gentools
from simutil import *
import os
import sys

ia = gentools(['ia'])[0]
sm = gentools(['sm'])[0]
me = gentools(['me'])[0]

global __version__
__version__ = '1.4'


if __name__ == '__main__':
    SDimage = 'SDTEST.image'
    SDchannels = -1
    NuReGrid = [93117929541.699982, 30548.095700003487, 3501]
    SDbaseline = 7.0
    nSDvis = 500
    inputvis = 'ACA/X40db/uid___A002_Xb903d6_X40db.ms.split.cal'
    field = "G357"
    inputspw = 0
    inputchan = [0, 0]
    wgtfac = 1.0
    over_resolve = 1.0
    scale = 1.0
    outputvis = "G357.spw17.I.sd.im.fits.SD2VIS"
    Python_DFT = False

    #  NuReGrid = []
    #  SDimage             = 'TX_Psc.spw23.I.sd.im.fits' #  Total-Power (i.e., Single
    #   Dish) CASA image.
    #  SDchannels          = [200, 3300]        #  List of two integers, which give the
    #   range of spectral channels (in the
    #   image) that will be converted into
    #   visibilities. If not a 2-element
    #   list, the whole image is taken.
    #  SDbaseline          =        7.0        #  Maximum baseline for TP visibilities
    #   (ideally, the antenna diameter, but
    #   lower values may work better
    #   sometimes).
    #  nSDvis              =       1000        #  Number of visibilities in synthetic
    #   dataset.
    #  inputvis            = 'TE_7M_spw3.ms'   #  Measurement set with the visibilities
    #   that are going to be concatenated
    #   with this TP data. This is used to
    #   set the appropriate weights to the
    #   TP visibilities.
    #  inputspw            =          0        #  Spectral window that will be used in
    #   the TP+interferometry imaging. This
    #   is used to set the appropriate
    #   weights (and/or scaling factor) to
    #   the TP visibilities.
    #  inputchan           =         -1        #  Range of spectral channels that will
    #   be used in the TP+interferometry
    #   imaging. This is used to set the
    #   appropriate weights (and/or scaling
    #   factor) to the TP visibilities. If
    #   an integer is given, only that
    #   channel is selected, unless the
    #   integer is negative (in that case,
    #   the whole spw is selected).
    #  wgtfac              =        1.0        #  Ratio between the total weight of the
    #   interferometric visibilities and the
    #   total weight of the TP-based
    #   synthetic visibilities. High values
    #   will over-weight the TP data.
    #  over_resolve        =        1.0        #  Over-resolution factor when
    #   unconvolving the TP image from the
    #   TP beam. Default value will use the
    #   exact TP beam, as found from the
    #   header of the TP image. Higher
    #   values will unconvolve narrower
    #   beams.
    #  scale               =        1.0        #  Scaling factor for the amplitudes of
    #   the TP data. If you set over_resolve
    #   != 1.0, you should set scale =
    #   1./(over_resolve)**2, to conserve
    #   the total flux density in the image.
    #   If lower than 0, the scaling factor
    #   will be derived automatically using
    #   baselines (from the interferometry
    #   data) shorter than abs(scale) (in
    #   meters).
    #  outputvis           = 'TP_spw3.ms'  #  Name of output measurement set
    #  Python_DFT          =       True        #  If a FITS file is used, it is very
    #   possible to encounter a (quite)
    #   silly CASA error related to the sm
    #   tool. In such a case, set this to
    #   True, so that TP2vis will compute
    #   the DFT by itself (much slower
    #   approach, since it has to use some
    #   Python loops).

def SD2vis(SDimage='', SDchannels=-1,  # NuReGrid=[],
           SDbaseline=7.0, nSDvis=1000, inputvis='', field='',
           inputspw=0, inputchan=[0, 0], wgtfac=1.0, over_resolve=1.0,
           scale=1.0, outputvis='SD2vis.ms', Python_DFT=False):

    greetings = "\n\n\n  SD2VIS - VERSION %s - NORDIC ARC NODE\n\n\n" % __version__
    units = {'deg': 3600., 'rad': 180./np.pi*3600.,
             'arcsec': 1., 'arcmin': 60., 'mas': 0.001}
    sig2FWHM = 2.35482
    util = simutil('')

    NuReGrid = []  # Feature turned off by now

    if SDimage[-1] == "/":
        SDimage = SDimage[:-1]

    if len(inputvis) > 0 and inputvis[-1] == "/":
        inputvis = inputvis[:-1]

    if outputvis[-1] == "/":
        outputvis = outputvis[:-1]

    def printMsg(msg):
        """Print message on terminal and log"""
        print msg
        casalog.post('SD2vis: '+msg)

    def printErr(msg):
        """Same as printMsg, but raise an exception afterwards"""
        print msg
        casalog.post('SD2vis ERROR: '+msg)
        raise Exception(msg)

    def truncGauss2D(mu, sigma, Rmax):
        """Draw numbers from a 2D Gaussian distribution, but within a range of allowed radii"""
        a = rn.gauss(mu, sigma)
        b = rn.gauss(mu, sigma)
        q = a*a + b*b
        r = Rmax**2.
        while (q > r):
            a = rn.gauss(mu, sigma)
            b = rn.gauss(mu, sigma)
            q = a*a + b*b
        return a, b

    printMsg(greetings)

    # SANITY CHECKS:
    if type(inputchan) is list:
        if len(inputchan) == 2:
            chans = [abs(int(inputchan[0])), abs(int(inputchan[1])+1)]
        else:
            printErr(
                '\n\nERROR! Bad inputchan! Should be a list of two elements!\n\n')
    else:
        ich = int(inputchan)
        if ich < 0:
            chans = [0, -1]  # THE WHOLE SPW
        else:
            chans = [abs(int(inputchan)), abs(int(inputchan)+1)]

    chans = [min(chans), max(chans)]

    if not os.path.exists(SDimage):
        printErr('SD image does NOT exist!')

    if not os.path.exists(inputvis):
        printMsg(
            'WARNING! No inputvis is given! Only one pointing will be simulated!')

    # PREPARE IMAGE (REGRID AND RE-ARRANGE AXES):
    os.system('rm -rf %s.COPY %s.UNCONVOLVED' % (SDimage, SDimage))
    if len(NuReGrid) == 0:
        os.system('cp -r %s %s.COPY' % (SDimage, SDimage))
    else:
        try:
            printMsg('Will regrid image to Nu0: %.8e GHz, dNu: %.8e MHz and %i channels' % (
                NuReGrid[0]/1.e9, NuReGrid[1]/1.e6, NuReGrid[2]))
        except:
            printErr(
                'BAD NuReGrid. It should ba a list of two floats and an integer')

        ia.open(SDimage)
        mycs = ia.coordsys()
        myrec = mycs.torecord()
        Spec = [i for i in myrec.keys() if 'spectral' in i][0]
        myrec[Spec]['wcs']['cdelt'] = float(NuReGrid[1])
        myrec[Spec]['wcs']['crpix'] = 0.0
        myrec[Spec]['wcs']['crval'] = float(NuReGrid[0])
        Shape = ia.shape()
        NuSh = list(ia.summary()['axisnames']).index('Frequency')
        Shape[NuSh] = int(NuReGrid[2])
        os.system('rm -rf %s.COPY' % SDimage)
        ia.regrid(outfile='%s.COPY' %
                  SDimage, csys=myrec, shape=Shape, overwrite=True)
        ia.close()

    # CONVERT INTO A CASA IMAGE:
    os.system('rm -rf SD2VIS.TEMP')
    try:  # A FITS
        import pyfits as pf
        testf = pf.open('%s.COPY' % SDimage)
        testf.close()
        importfits(fitsimage='%s.COPY' %
                   SDimage, imagename='SD2VIS.TEMP', overwrite=True)
        printMsg('Imported FITS')
    except:  # A CASA IMAGE?
        os.system('cp -r %s.COPY SD2VIS.TEMP' % (SDimage))

    os.system('rm -rf %s.UNCONVOLVED' % SDimage)
    success = ia.open('SD2VIS.TEMP')
    if not success:  # NEITHER A CASA IMAGE
        printErr('NOT A VALID SD image!')

    # RE-ARRANGE:
    NuSh = list(ia.summary()['axisnames']).index('Frequency')
    if NuSh == 2:
        ia.transpose(order='0132', outfile='%s.UNCONVOLVED' % SDimage)
        ia.close()
    elif NuSh == 3:
        ia.close()
        os.system('cp -r SD2VIS.TEMP %s.UNCONVOLVED' % SDimage)
    else:
        ia.close()
        printErr('\n\nERROR! Odd Frequency Axis!\n\n')

    #  ia.close()
    ia.open('%s.UNCONVOLVED' % SDimage)
    data = ia.getchunk()
    data[np.isnan(data)] = 0.0  # Unset masked pixels
    summ = ia.summary()
    outframe = (ia.coordsys().referencecode())[-1]

    source = str(field).replace(' ', '_')

    #  ia.close()
    ndim = len(summ['shape'])
    NPIX = summ['shape'][:2]
    INCR = [ic*units[summ['axisunits'][0]]/units['rad']
            for ic in summ['incr'][:2]]

    nudim = 3  # list(summ['axisnames']).index('Frequency')
    #  if nudim not in [2,3]:
    #    printErr('THIS IS NOT AN IMAGE CUBE!')
    Freqs = (np.arange(summ['shape'][nudim]) + summ['refpix']
             [nudim])*summ['incr'][nudim]+summ['refval'][nudim]

    if type(SDchannels) is list and len(SDchannels) == 2:
        try:
            Freqs = Freqs[np.min(SDchannels):np.max(SDchannels)+1]
            SDchans = [SDchannels[0], SDchannels[1]+1]
        except:
            printErr('BAD CHANNEL RANGE FOR IMAGE!')
    else:
        SDchans = [0, len(Freqs)]

    nchan = len(Freqs)  # summ['shape'][nudim]

    try:
        stdim = list(summ['axisnames']).index('Stokes')
        nst = summ['shape'][stdim]
    except:
        stdim = 0
        nst = 0

    # We assume a square image:
    imsize = np.abs(summ['incr'][0])*summ['shape'][0] * \
        units[summ['axisunits'][0]]/units['rad']  # in radians
    UVsize = summ['shape'][0]/imsize/2.  # in lambda

    psize = np.abs(summ['incr'][0]*units[summ['axisunits'][0]]/units['rad'])
    printMsg('PIXEL SIZE: %.3f arcsec' % (psize*180/np.pi*3600.))

    if summ['unit'].lower() == 'jy/beam':
        printMsg('SD image is convolved with beam!\nWill unconvolve in UV space')
        Maj = summ['restoringbeam']['major']
        Min = summ['restoringbeam']['minor']
        PA = summ['restoringbeam']['positionangle']
        BM = Maj['value']*units[Maj['unit']]/sig2FWHM/units['rad']
        Bm = Maj['value']*units[Min['unit']]/sig2FWHM/units['rad']
        Pang = PA['value']*units[PA['unit']]/units['rad']
        printMsg('\nBEAM: %.3f x %.3f arcsec (PA: %.0f deg.)\n' %
                 (BM*180/np.pi*3600., Bm*180/np.pi*3600., Pang*180./np.pi))
        Fsigma = 2.*(np.pi**2.)*BM*Bm
    else:
        Fsigma = 0.0

    # image pixel coordinates:
    xx = (np.arange(summ['shape'][0]) - summ['refpix'][0]) * \
        summ['incr'][0]*units[summ['axisunits'][0]]/units['rad']
    npix = summ['shape'][0]

    # FFT of model:
    datfft = np.zeros(np.shape(data), dtype=np.complex128)
    for sti in range(nst):
        #    if nudim==3:
        for i in range(SDchans[0], SDchans[1]):
            datfft[:, :, sti, i] = np.fft.fftshift(
                np.fft.fft2(data[:, :, sti, i]))

    lam = 3.e8/np.min(Freqs)
    MaxQ2 = (SDbaseline/lam)**2.
    MaxQ = np.sqrt(MaxQ2)

    # UNCONVOLVE TP BEAM:
    if Fsigma > 0.0:
        printMsg('Deconvolving beam')

        UU = np.linspace(-UVsize, UVsize,
                         int(summ['shape'][0]), endpoint=False)
        Uones = np.ones(summ['shape'][1])

        VV = np.linspace(-UVsize, UVsize,
                         int(summ['shape'][1]), endpoint=False)
        Vones = np.ones(summ['shape'][0])

        Q2 = (np.power(np.outer(UU, Uones), 2.) +
              np.power(np.outer(Vones, VV), 2.))

        Q2[Q2 > MaxQ2] = 0.0

        UVTaper = np.exp(Fsigma*Q2/(over_resolve**2.))/(Fsigma/np.pi)*psize**2.
        UVTaper[Q2 > MaxQ2] = 0.0

        for sti in range(nst):
            #      if nudim==3:
            for i in range(SDchans[0], SDchans[1]):
                data[:, :, sti, i] = np.fft.ifft2(
                    np.fft.ifftshift(datfft[:, :, sti, i]*UVTaper)).real

    # Save unconvolved & PBuncorrected image:
    #  ia.open('%s.UNCONVOLVED'%SDimage)
    data2 = ia.getchunk()
    data2[:] = data  # *Pbeam
    ia.setbrightnessunit('Jy/pixel')
    #  ia.setrestoringbeam(major='%.10farcsec'%(psize*units['rad']),minor='%.10farcsec'%(psize*units['rad']),pa='0deg')
    ia.putchunk(data2)
    ia.close()

    # Generate ms:
    TELNAME = 'ACA'
    if os.path.exists(outputvis):
        os.system('rm -rf %s' % outputvis)

    sm.open(outputvis)

    #  arr = os.getenv('CASAPATH').split(' ')[0]+'/data/alma/simmos/aca.all.cfg'
    #  stx, sty, stz, std, nam, nant, antn = util.readantenna(arr)
    NANT = 10
    stx = np.linspace(-100.*NANT, 100.*NANT, NANT)
    sty = np.copy(stx)
    stz = np.ones(len(stx))
    mount = 'alt-az'

    antn = [TELNAME for i in range(len(stx))]
    nam = ['DP%02i' % i for i in range(len(sty))]
    diams = [7.0]*NANT  # for xi in stx]
    refloc = me.observatory(TELNAME)
    rsystem = 'local'

    sm.setconfig(telescopename=TELNAME,
                 x=stx, y=sty, z=stz,
                 dishdiameter=diams,
                 mount=mount,
                 antname=TELNAME, padname=nam,
                 coordsystem=rsystem,
                 referencelocation=refloc)

    sm.setfeed(mode='perfect X Y', pol=[''])
    sm.setspwindow(spwname='TP0', freq='%.8fGHz' % (Freqs[0]/1.e9),
                   deltafreq='%.9fGHz' % (summ['incr'][nudim]/1.e9),
                   freqresolution='%.9fGHz' % (summ['incr'][nudim]/1.e9),
                   nchannels=len(Freqs), refcode=outframe, stokes='XX YY')

    RA = summ['refval'][0]
    Dec = summ['refval'][1]
    arc = summ['axisunits'][0]

    #####################
    # Set pointings:
    if len(inputvis) > 0:
        tb.open(os.path.join(inputvis, 'FIELD'))
        ff = np.where(tb.getcol('NAME') == field)[0]
        if len(ff) == 0:
            printErr('\n\nERROR! Field %s is NOT found in inputvis!\n\n' % field)
        RAs = tb.getcol('PHASE_DIR')[0, 0, :][ff]
        Decs = tb.getcol('PHASE_DIR')[1, 0, :][ff]
        NPOINT = len(RAs)
        printMsg('GOT %i POINTINGS FROM MS' % NPOINT)
        Pointings = ['%.16frad %.16frad' %
                     (RAs[i], Decs[i]) for i in range(NPOINT)]
        tb.close()
    else:
        printMsg('WARNING! No inputvis given! Only one pointing will be used!')
        NPOINT = 1
        Pointings = ['%.16f%s %.16f%s' % (RA, arc, Dec, arc)]
    #####################

    for i in range(NPOINT):
        sm.setfield(sourcename='%s_%i' % (source, i), sourcedirection=Pointings[i],
                    calcode='TARGET', distance='0m')

    sm.setlimits(shadowlimit=0.001, elevationlimit='1deg')
    sm.setauto(autocorrwt=0.0)

    refdate = '2017/01/01/12:00:00'
    mereftime = me.epoch('TAI', refdate)

    tint = 2  # seconds (integer)

    sm.settimes(integrationtime='%is' %
                tint, usehourangle=True, referencetime=mereftime)

    nbas = len(stx)*(len(stx)-1)/2
    ntime = nSDvis/nbas + 1
    times0 = []
    times1 = []
    snames = []

    for j in range(NPOINT):
        times0.append('%is' % (j*ntime*tint))
        times1.append('%is' % ((j+1)*ntime*tint))
        snames.append('%s_%i' % (source, j))

    sm.observemany(sourcenames=snames, spwname='TP0',
                   starttimes=times0, stoptimes=times1)

    # Change baseline coordinates to a random
    # Gaussian distribution:
    U = np.zeros(ntime*NPOINT*nbas)
    V = np.zeros(ntime*NPOINT*nbas)
    if Fsigma > 0.0:
        Gsigma = np.sqrt(1./(2.*Fsigma))*lam
    else:
        AuxF = 2.*(np.pi**2.)*psize**2.
        Gsigma = np.sqrt(1./(2.*AuxF))*lam
    for i in range(ntime*NPOINT*nbas):
        U[i], V[i] = truncGauss2D(0., Gsigma, SDbaseline)

    # Write UV coordinates:
    tb.open(outputvis, nomodify=False)
    uvw = tb.getcol('UVW')
    nvis = np.shape(uvw)[-1]
    tpout = np.logical_or(uvw[0, :] != 0, uvw[1, :] != 0.0)
    uvw[0, tpout] = U
    uvw[1, tpout] = V
    uvw[2, tpout] = 0.0
    tb.putcol('UVW', uvw)

    nivis = sum(tpout)

    # Observe the TP image:
    printMsg('Going to compute DFT of image')
    tb.close()

    # Derive/apply TP-to-interferometry amplitude scaling:
    if scale > 0.:
        printMsg('\n\nApplying a scaling factor of %.2f\n\n' % scale)
    else:
        printMsg(
            '\n\nDeriving scaling factor from baselines shorter than %.2f m.\n\n' % (-scale))
        if os.path.exists('%s.TPOBS' % inputvis):
            os.system('rm -rf %s.TPOBS' % inputvis)
        os.system('cp -r %s %s.TPOBS' % (inputvis, inputvis))
        sm.close()
        sm.openfromms('%s.TPOBS' % inputvis)
        sm.setdata(spwid=int(inputspw), fieldid=range(NPOINTS))
        sm.setvp(dovp=True, usedefaultvp=False)
        sm.predict(imagename='%s.UNCONVOLVED' % SDimage)
        sm.close()
        sm.openfromms(outputvis)
        tb.open('%s.TPOBS' % inputvis)
        intbas = tb.getcol('UVW')
        tpdats = tb.getcol('DATA')
        tb.close()
        tb.open(inputvis)
        visdats = tb.getcol('DATA')
        tb.close()
        spchan = np.shape(visdats)[1]

        # Select the whole spw if chans=[0,-1]
        if chans[1] == -1:
            chans[1] = spchan

        if chans[1] > spchan:
            printErr(
                'ERROR: There are only %i channels! Change the channel range!' % spchan)
        else:
            printMsg('Will use channels from %i to %i (included)' %
                     (chans[0], chans[1]-1))

        Qs = intbas[0, :]**2. + intbas[1, :]**2.
        Q = np.logical_and(Qs > 0., Qs < scale**2.)
        Mask = np.logical_and(
            visdats[:, chans[0]:chans[1], Q] != 0., tpdats[:, chans[0]:chans[1], Q] != 0.)
        if np.sum(Mask) > 0:
            scale = np.average(np.abs(
                visdats[:, chans[0]:chans[1], Q][Mask])/np.abs(tpdats[:, chans[0]:chans[1], Q][Mask]))
        else:
            printMsg(
                'WARNING: Did not find baselines as short as %.1f m.!' % (-scale))
            scale = 1.0

        printMsg('\n\nWill apply a scaling factor of %.1e\n\n' % scale)

    ia.open('%s.UNCONVOLVED' % SDimage)
    data2 = ia.getchunk()
    data2[:] *= scale
    ia.putchunk(data2)
    ia.close()

    if Python_DFT:
        data2sq = np.squeeze(data2)

        ms.open(outputvis, nomodify=False)
        ms.selectinit(datadescid=0)
        freqs = np.squeeze(ms.range('CHAN_FREQ')['chan_freq'])
        tpdata = ms.getdata(['u', 'v', 'data', 'flag'])
        Ulamb = 2.*np.pi*freqs[:, np.newaxis] / \
            299792458.*tpdata['u'][np.newaxis, :]
        Vlamb = 2.*np.pi*freqs[:, np.newaxis] / \
            299792458.*tpdata['v'][np.newaxis, :]
        tpdata['flag'] = False
        TotVis = np.zeros(len(freqs), dtype=np.complex128)
        RASum = np.zeros(len(freqs))
        for k in range(nvis):

            sys.stdout.write('\r Doing Vis. %i of %i' % (k+1, nvis))
            sys.stdout.flush()
            TotVis[:] = 0.0

            for i in range(NPIX[0]):
                dRA = (i-NPIX[0]/2.)*INCR[0]
                RASum[:] = Ulamb[:, k]*dRA

                for j in range(NPIX[1]):
                    dDec = (j-NPIX[0]/2.)*INCR[1]
                    TotVis[:] += data2sq[i, j, SDchans[0]:SDchans[1]] * \
                        np.exp(1.j*(RASum + Vlamb[:, k]*dDec))

            tpdata['data'][0, :, k] = TotVis
            tpdata['data'][1, :, k] = TotVis

        printMsg('\n\n DONE!\n')
        ms.putdata(tpdata)
        ms.close()

    else:
        sm.setdata(fieldid=range(NPOINT))
        sm.setvp(dovp=True, usedefaultvp=False)
        sm.predict(imagename='%s.UNCONVOLVED' % SDimage)

    sm.close()
    sm.done()

    # Set weight:
    printMsg('COMPUTING WEIGHT RATIO')

    if len(inputvis) > 0:
        ms.open(inputvis)

        ms.selectinit(datadescid=int(inputspw))
        msdata = ms.getdata(['weight', 'flag', 'u', 'v'])

        ms.close()

        NTIME = np.shape(msdata['u'])[0]
        NPOL = np.shape(msdata['weight'])[0]
        tpout = np.logical_or(msdata['u'] != 0, msdata['v'] != 0.0)
        if len(np.shape(msdata['flag'])) == 1:
            NCHAN = np.shape(msdata['flag'])[0]/NTIME/NPOL
            fgs = msdata['flag'].reshape(NPOL, NCHAN, NTIME)
        else:
            fgs = msdata['flag']

        fgoutput = np.sum(np.logical_not(fgs), axis=1) > 0
        goods = np.logical_and(tpout, fgoutput)

        wgtsum = np.sum(msdata['weight'][goods])
        nvis = max(1, np.sum(goods))

        printMsg('Average weight/vis in ms: %.3e' % (wgtsum/nvis))

        SDwgt = wgtfac*wgtsum/nivis

    else:
        SDwgt = 1.0

    printMsg('Will set SD visibility weights to: %.3e' % (SDwgt))

    clearcal(outputvis, addmodel=True)
    ms.open(outputvis, nomodify=False)
    ms.selectinit(datadescid=0)
    msdata = ms.getdata(['weight'])
    msdata['weight'][:] = SDwgt
    ms.putdata(msdata)
    ms.close()
    ia.close()
