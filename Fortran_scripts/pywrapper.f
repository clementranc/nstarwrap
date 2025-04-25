C
C======================================================================
C
C This is a wrapper to enable the usage of NSTAR-MCMC from python.
C
C                 Version 1.0 - 6 April 2022
C           Original implementation: Clément Ranc
C
C This code runs NSTAR-MCMC_PYWRAPPER, which is a modification of 
C NSTAR-MCMC. NSTAR-MCMC was initially developped by Sean Terry:
C
C    - Code and documentation: https://github.com/skterry/daophot_mcmc
C    - Article: Terry et al. 2021, AJ, 161, 54.
C 
C NSTAR-MCMC is a modified version of the subroutine NSTAR of DAOPHOT-II
C developped by Peter Stetson:
C    - Code: http://www.star.bris.ac.uk/~mbt/daophot/
C    - Documentation: http://www.astro.wisc.edu/sirtf/daophot2.pdf
C
C=======================================================================
C
C
C=======================================================================
C
C   This subroutine runs the ATTACH command to get the picture size.
C
C=======================================================================
C
      subroutine image_size(file, ncol, nrow, filel) 
     . bind(c, name="imagesize")

      use iso_c_binding
      implicit none

      integer(c_int), intent(inout)    :: ncol
      integer(c_int), intent(inout)    :: nrow
      integer(c_int), intent(in), value :: filel
      character(c_char), intent(inout) :: file(filel)
      character*40 file2
      integer i

      file2=''
      do i = 1, filel
        file2(i:i) = file(i)
      end do

C     Image dimensions
      CALL ATTACH (file2, NCOL, NROW)

      end subroutine
C
C=======================================================================
C
C   This subroutine load the picture and send it to python.
C
C=======================================================================
C
      subroutine load_image(nrow, data, maxcol, ier) 
     . bind(c,name="loadimage")

      use iso_c_binding
      implicit none

      integer(c_int), intent(in), value :: nrow, maxcol
      integer(c_int), intent(inout)     :: ier
      real(c_float), intent(inout)      :: data(maxcol,maxcol)
      integer ly

      ly=1
      call rdsect (1, ly, nrow, data, maxcol, ier)

      end subroutine
C
C=======================================================================
C
C   This subroutine run nstar-mcmc_pywrapper and communicate with python
C
C=======================================================================
C
      subroutine pynstar(nrow, ncol, data, maxcol,
     .  xmcmc1,ymcmc1,xmcmc2,ymcmc2, xmcmc3,ymcmc3,xmcmc4,ymcmc4, 
     .  fratio12, fratio13, fratio14, z0, chi2,
     .  watch, fitrad, e1, e2, psffil, psfl, grpfil, grpl,
     .  xmin, xmax, ymin, ymax, rfac, fit_stars, zpmag, residuals)
     .  bind(c,name="pynstar")

      use iso_c_binding
      implicit none

      integer(c_int), intent(in), value :: nrow, ncol, maxcol
      integer(c_int), intent(in), value :: psfl, grpl
      real(c_float), intent(inout) :: zpmag
      real(c_float), intent(inout) :: data(maxcol,maxcol)
      real(c_double), intent(inout) :: residuals(maxcol*maxcol)
      real(c_float), intent(inout) :: xmcmc1,ymcmc1,xmcmc2,ymcmc2
      real(c_float), intent(inout) :: xmcmc3,ymcmc3,xmcmc4,ymcmc4
      real(c_double), intent(inout) :: fratio12,fratio13,fratio14, z0
      real(c_float), intent(in), value :: watch, fitrad, e1, e2
      real(c_double), intent(in), value :: rfac
      real(c_double), intent(inout) :: chi2
      integer(c_int), intent(inout) :: xmin, xmax, ymin, ymax
      integer(c_int), intent(in), value :: fit_stars
      character(c_char), intent(inout) :: psffil(psfl), grpfil(grpl)

      integer maxpsf, maxexp, maxpar
      parameter (maxpsf=407, maxexp=10, maxpar=6)
      real par(maxpar), psf(maxpsf,maxpsf,maxexp)

C This commented lines are related to a future update of the code
C    call nstar  (par, maxpar, psf, maxpsf, maxexp, data,
C    .     ncol, nrow, maxcol, watch, fitrad, e1, e2,
C    .     xmcmc1,ymcmc1,xmcmc2,ymcmc2,xmcmc3,ymcmc3,xmcmc4,ymcmc4,
C    .     fratio12,fratio13,fratio14,z0,chi2,
C    .     psffil, psfl, grpfil, grpl, xmin, xmax, ymin, ymax, rfac,
C    .     residuals, fit_stars, zpmag)
C
C    end subroutine

      call nstar  (par, maxpar, psf, maxpsf, maxexp, data, 
     .     ncol, nrow, maxcol, watch, fitrad, e1, e2,
     .     xmcmc1,ymcmc1,xmcmc2,ymcmc2,xmcmc3,ymcmc3,xmcmc4,ymcmc4,
     .     fratio12,fratio13,fratio14,z0,chi2,
     .     psffil, psfl, grpfil, grpl, xmin, xmax, ymin, ymax, rfac,
     .     fit_stars, zpmag, residuals)

      end subroutine
