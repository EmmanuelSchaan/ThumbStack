

for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:

   ###################################################################################
   ###################################################################################
   # kSZ null tests

   rAp = rKsz[catalogKey+'_pactf150daynight20200228maskgal60']
   #
   # fiducial uncertainty
   ksz150 = ksz[catalogKey+'_pactf150daynight20200228maskgal60']
   sKsz150 = sKsz[catalogKey+'_pactf150daynight20200228maskgal60']
   # Mariana's velocities
   if catalogKey=='cmass_kendrick':
      catalogMariana = 'cmass_mariana'
   elif catalogKey=='lowz_kendrick':
      catalogMariana = 'lowz_mariana'
   ksz150Mariana = ksz[catalogMariana+'_pactf150daynight20200228maskgal60']
   sKsz150Mariana = sKsz[catalogMariana+'_pactf150daynight20200228maskgal60']
   #
   # 150 - tilec cmb, as a consistency check
   ksz150MinusTilecCmb = ksz[catalogKey+'_pactf150daynight20200228maskgal60_minus_tilecpactcmbksz']
   sKsz150MinusTilecCmb = sKsz[catalogKey+'_pactf150daynight20200228maskgal60_minus_tilecpactcmbksz']
   #
   # 150 reconv to 90 minus 90, to check consistency
   ksz150Reconv90Minus90 = ksz[catalogKey+'_pactf150reconvto90minus90daynight20200228maskgal60']
   sKsz150Reconv90Minus90 = sKsz[catalogKey+'_pactf150reconvto90minus90daynight20200228maskgal60']
   #
   # tilec y no cmb, to check for tSZ contamination
   kszYNoCmb = ksz[catalogKey+'_tilecpactynocmb'] * yTomuK150
   sKszYNoCmb = sKsz[catalogKey+'_tilecpactynocmb'] * yTomuK150
   #
   # 150 - tilec cmb no cib, to check for dust contamination
   ksz150MinusCmbNoCib = ksz[catalogKey+'_pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib']
   sKsz150MinusCmbNoCib = sKsz[catalogKey+'_pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib']




   # kSZ pipeline null tests
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   #
   ax.axhline(0., c='k', lw=1)
   #
   # Uncertainty band
   ax.fill_between(rAp, - sKsz150, sKsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
   #
   # V-shuffle mean
   ax.errorbar(rAp, ksz150VShuffleMean, yerr=sKsz150VShuffleMean, fmt='--', label='mean of 100 v-shuffles')
   #
   # Average of many mocks
   ax.errorbar(rAp + 0.025, meanStackedKszGRF, yerr=sStackedKszGRF, fmt='--', label=r'mean of '+str(nMocks)+' mocks')
   #
   # Mariana - Kendrick
   #ax.errorbar(rAp, ksz150MKDiff, yerr=sKsz150MKDiff, fmt='-', c='r', label=r'$v_\text{Mariana} - v_\text{Kendrick}$')
   #ax.errorbar(rAp, ksz150, yerr=sKsz150, fmt='-', label='K')
   #ax.errorbar(rAp, ksz150Mariana, yerr=sKsz150Mariana, fmt='-', label='M')
   ax.errorbar(rAp + 0.05, (ksz150-ksz150Mariana), yerr=sKsz150, fmt='-', label=r'$v_\text{Kendrick} - v_\text{Mariana}$')
   #
   # 150 - tilec cmb
   ax.errorbar(rAp + 0.075, ksz150MinusTilecCmb, yerr=sKsz150MinusTilecCmb, fmt='-', label='150 - TileC CMB/kSZ')
   #
   # 150 reconv to 90 minus 90
   ax.errorbar(rAp + 0.1, ksz150Reconv90Minus90, yerr=sKsz150Reconv90Minus90, fmt='-', label='150\' - 90')
   #
   ax.set_ylim((-10., 15.))
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'kSZ pipeline null tests')
   ax.set_ylim((-6., 6.))
   #
   path = pathFig + "pipenulltests_ksz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()




   # kSZ foreground null tests
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   ax.axhline(0., c='k', lw=1)
   #
   # Uncertainty band
   ax.fill_between(rAp, - sKsz150, sKsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
   #
   # kSZ on TileC y no CMB map
   ax.errorbar(rAp, kszYNoCmb, yerr=sKszYNoCmb, label=r'TileC y no CMB')
   #
   # cmbksz no cib, to check for dust
   ax.errorbar(rAp + 0.025, ksz150MinusCmbNoCib, yerr=sKsz150MinusCmbNoCib, fmt='-', label=r'150 - TileC CMB/kSZ no CIB')
   #
   ax.set_ylim((-10., 15.))
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'kSZ foreground null tests')
   #ax.set_ylim((-2., 2.))
   #
   path = pathFig + "fgnulltests_ksz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   ###################################################################################
   # tSZ null tests

   rAp = rTsz[catalogKey+'_pactf150daynight20200228maskgal60']
   #
   # fiducial uncertainty
   sTsz150 = sTsz[catalogKey+'_pactf150daynight20200228maskgal60'] 
   #
   # 150 - y, to check for map consistency
   tsz150MinusY = tsz[catalogKey+'_pactf150daynight20200228maskgal60_minus_tilecpactymuk']
   sTsz150MinusY = sTsz[catalogKey+'_pactf150daynight20200228maskgal60_minus_tilecpactymuk']
   #
   # y - y no CIB
   tszYMinusYNoCib = tsz[catalogKey+'_tilecpactyminusynocib'] * yTomuK150
   sTszYMinusYNoCib = sTsz[catalogKey+'_tilecpactyminusynocib'] * yTomuK150
   #
   # 150' - 90, after rescaling 90 to null tSZ
   # in order to check for the dust contamination
   tsz150Reconv90Minus90NoY = tsz[catalogKey+'_pactf150reconvto90minus90noydaynight20200228maskgal60']
   sTsz150Reconv90Minus90NoY = sTsz[catalogKey+'_pactf150reconvto90minus90noydaynight20200228maskgal60']



   # tSZ pipeline test
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   ax.axhline(0., c='k', lw=1)
   #
   # Uncertainty band
   ax.fill_between(rAp, - sTsz150, sTsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
   #
   # mean of GRF mocks
   ax.errorbar(rAp, meanStackedTszGRF, yerr=sStackedTszGRF, fmt='--', label=r'mean of '+str(nMocks)+' mocks')
   #
   # 150 - tilec y
   ax.errorbar(rAp, tsz150MinusY, yerr=sTsz150MinusY, fmt='--', label=r'150 - TileC y ')
   #
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'tSZ pipeline null tests')
   ax.set_ylim((-6., 6.))
   #
   path = pathFig + "pipenulltests_tsz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   # dust contamination to tSZ
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   ax.axhline(0., c='k', lw=1)
   #
   # Uncertainty band
   ax.fill_between(rAp, - sTsz150, sTsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
   #
   # 150 - tilec y
   ax.errorbar(rAp, tsz150MinusY, yerr=sTsz150MinusY, fmt='--', label=r'150 - TileC y')
   #
   # y - y no CIB
   ax.errorbar(rAp, tszYMinusYNoCib, yerr=sTszYMinusYNoCib, fmt='--', label=r'TileC y - y no CIB')
   #
   # 150' - 90 rescaled to null y
   ax.errorbar(rAp, tsz150Reconv90Minus90NoY, yerr=sTsz150Reconv90Minus90NoY, fmt='--', label=r"150\' - 90 no y")
   #
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{dust}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'Dust emission')
   ax.set_ylim((-6., 6.))
   #
   path = pathFig + "fgnulltests_tsz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   ###################################################################################
   # summary tSZ plot


   rAp = rTsz[catalogKey+'_pactf150daynight20200228maskgal60']

   # PACT 150
   tsz150 = tsz[catalogKey+'_pactf150daynight20200228maskgal60']
   sTsz150 = sTsz[catalogKey+'_pactf150daynight20200228maskgal60']
   #
   # PACT 90
   tsz90 = tsz[catalogKey+'_pactf90daynight20200228maskgal60']
   sTsz90 = sTsz[catalogKey+'_pactf90daynight20200228maskgal60']
   #
   # TileC y no Cib
   tszYNoCib = tsz[catalogKey+'_tilecpactynocib'] * yTomuK150
   sTszYNoCib = sTsz[catalogKey+'_tilecpactynocib'] * yTomuK150


   # tSZ-only from the TileC y no CIB map
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   ax.axhline(0., c='k', lw=1)
   #
   # PACT 150
   #ax.errorbar(rAp, tsz150, yerr=sTsz150, fmt='-', c='royalblue', label='150GHz')
   # PACT 90
   #ax.errorbar(rAp, tsz90, yerr=sTsz90, fmt='-', c='darkviolet', label='90GHz')
   # Tilec y no CIB
   ax.errorbar(rAp, tszYNoCib, yerr=sTszYNoCib, fmt='-', c='r', label='TileC y no CIB')
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_{\text{tSZ}}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'tSZ profile')
   #ax.set_ylim((0., 2.))
   #
   path = pathFig+"summary_tsz_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()


   # tSZ + dust plot at 150 and 90
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   ax.axhline(0., c='k', lw=1)
   #
   # PACT 150
   ax.errorbar(rAp, tsz150, yerr=sTsz150, fmt='-', c='royalblue', label='150GHz')
   # PACT 90
   ax.errorbar(rAp, tsz90, yerr=sTsz90, fmt='-', c='darkviolet', label='90GHz')
   # Tilec y no CIB
   ax.errorbar(rAp, tszYNoCib, yerr=sTszYNoCib, fmt='-', label='TileC y no CIB')
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_{\text{tSZ} + \text{dust}}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'tSZ + dust profile')
   #ax.set_ylim((0., 2.))
   #
   path = pathFig+"comparison_tsz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()






   ###################################################################################
   ###################################################################################





