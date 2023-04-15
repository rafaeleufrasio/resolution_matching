pro match_res,imagelist,uncert_list,RESMODE=resmode,GRIDSIZE=gridsize,REPROJECT=reproject,$
              IMAGEREF=imageref,BCKGSUB=bckgsub,TRUNC=trunc,FWHMLW=fwhmlw,TROUBLESHOOT=troubleshoot
;+
; NAME:
;   MATCH_RES
;
; PURPOSE:
;   TO CONVOLVE A LIST OF IMAGES TO THE SAME RESOLUTION AND ALIGN TO THE SAME ASTROMETRY
;
; EXPLANATION:
;   This routine convolves the input images (Spitzer's IRAC and MIPS
;   mosaics, Palomar, VLA, KP, 2MASS, GALEX) to the same resolution
;   and align convolved images to match the astrometry and pixel
;   scale of reference image
;
;   HEADERS WILL BE UPDATED FOR ALL IMAGE PROCESSING
;
; CALLING SEQUENCE:
;   MATCH_RES, imagelist, [uncert_list, RESMODE=, GRIDSIZE=, TRUNC=, IMAGEREF=, /BCKGSUB]
;
; INPUTS:
;   IMAGELIST => Name of the image list file
;   Headers of all (non-Spitzer) images must have PSF keyword
;     'GAUSS' - requires FWHM keyword
;     'GAUSS_ELLIPTICAL' - requires BMAJ, BMIN, and BPA keywords
;     'IDL_VAR' - require PSF_FILE keyword with the filename for the psf
;     'FITS' - require PSF_FILE keyword with the filename for the psf
;
; KEYWORDS
;   RESMODE - The resolution to be matched. Lowest resolution. Options:
;     => 'I1'    - IRAC channel 1  3.6 microns
;     => 'I2'    - IRAC channel 2  4.5 microns
;     => 'I3'    - IRAC channel 3  5.8 microns
;     => 'I4'    - IRAC channel 4  8.0 microns
;     => 'M1'    - MIPS channel 1  24 microns
;     => 'M2w'   - MIPS channel 2  Wide Field 70 microns
;     => 'M2n'   - MIPS channel 2  Narrow Field 70 microns 
;     => 'M3'    - MIPS channel 3  160 microns
;     => 'GAUSS' - Gaussian
;          This keyword requires specifying a fwhmlw keyword as the FWHM resolution
;     => filename for any fits file with 'PSF' header keyword (or a Spitzer image). PSF can be set up to
;          'GAUSS' - requires header keyword FWHM
;          'GAUSS_ELLIPTICAL' - requires additional BMAJ, BMIN and BPA header keywords
;
; OPTIONAL INPUT KEYWORDS : 
;   UNCERT_LIST - List of uncertanty files
;   GRIDSIZE    - Gridsize for PSF model (~280 pixels works ok for all irac/mips PSFs for 1.2"/pix)
;                 (default = 3^6 = 729)
;                 (in case using KP average PSFs (restored as IDL variables) then better to use odd number for gridsize so that psf is well centered)
;                 (in case using GALEX, gridsize should not be smaller than 120 pixels (size of GALEX average PSFs) and better to be even number)
;   TRUNC       - Program will truncate values gt trunc  (this is not good for some images like Palomar, so check the image first)
;                 (default = 500)
;   IMAGEREF    - The image to be used as reference for Headmaster
;                 if IMAGEREF='FIRST' then the first on the input imagelist
;                 You can also specify the image you want, IMAGEREF='imageIwant.fits'
;   /BCKGSUB     - If set, then the backgrounds will be subtracted before performing the convolution
;
; OUTPUTS:
;   For each input image.fits there will be a
;     image_conv+'resmode'+.fits    => convolved to chosen resolution mode.
;
; OPTIONAL OUTPUTS:
;   For each input image.fits there will be a
;     image_conv+'resmode'+ok.fits  => convolved and aligned to headmaster FITS header astrometry.
;
; HEADERS WILL BE UPDATED FOR ALL IMAGE PROCESSING
;
; OPTIONAL OUTPUT KEYWORDS:
;
; EXAMPLES:
;   Convolve images in file list to match resolution of MIPS channel 1, doing first
;   a background subtraction and gridsize=200 and aligning astrometry with 'myimage.fits'  :
;     IDL> MATCH_RES,'list.txt',resmode='M1',gridsize=3^5*5,imageref='myimage.fits',/BCKGSUB
;
;   Convolve images in file list to match resolution of MIPS channel 1, doing first
;     IDL> MATCH_RES,'list.txt',resmode='GAUSS',fwhmlw=15.,gridsize=3^5*5,/BCKGSUB
;
; RESTRICTIONS:
;   
; PROCEDURES CALLED:
;  When Spitzer images are on the list (or the lowest resolution is one of Spitzer's)
;    this routine will call FUNCTION: stinytim. Stinytim requires template files:
;      I1.par,I2.par,I3.par,I4.par,M1.par,M2n.par,M2w.par,M3.par
;    It will also call stiny2: Program that can be downloaded from Spitzer website 
;
; REVISION HISTORY:
;   Written     G. Kober & R. Eufrasio     September, 2010


if (N_elements(troubleshoot) eq 0) then troubleshoot = 0

; reads filenames and uncertainties
readcol,imagelist,filename,format='A'
if n_elements(uncert_list) eq 1 then readcol,uncert_list,unc_filename,format='A'
nel = n_elements(filename)
if n_elements(uncert_list) eq 1 then nelunc = n_elements(unc_filename) else nelunc=0
if (nelunc ne 0) and (nel ne nelunc) then begin
  print, 'ERROR: List of uncertainty files do not have the same size as files list'
  STOP
endif

; sets the default model grid size if no value provided (default = 3^6 = 729)
if (n_elements(gridsize) ne 1) then gridsize = 3^6
gridsize = long(gridsize)

; sets the default reference image if no value is provided 
if (n_elements(imageref) ne 1) then imageref='FIRST'

; reads headmaster for reprojection
if keyword_set(REPROJECT) then begin
  if (n_elements(imageref) eq 1) then imageref = string(imageref)
  if strupcase(imageref) eq 'FIRST' then imageref=filename[0]  ; First image is the default for imageref
  headmaster = headfits(imageref)
  print,'HEADMASTER ='+imageref
endif
  
; Resolution to convolve to
; sets the default resolution mode if no value is provided 
if n_elements(resmode) ne 1 then resmode='FIRST'
resmode = string(resmode)
if strupcase(resmode) eq 'FIRST' then begin
  hlw = headfits(filename[0])  ; First image is the default for imageref
endif

;mode = string(resmode)
mode = strupcase(resmode)

if strtrim(mode) eq 'GAUSS' then begin
  if (n_elements(fwhmlw) eq 1) then begin
    fwhmlw = double(fwhmlw)         ; Gaussian FWHM in arcseconds
    print, 'mode = '+mode+' and FWHM = '+strtrim(fwhmlw,2)+' arcseconds'
  endif
  if (n_elements(fwhmlw) eq 0) then begin
    print, 'ERROR: Resmode "GAUSS" requires specifying a "fwhmlw" keyword' & STOP
  endif
  if (n_elements(fwhmlw) gt 1) then begin
    print, 'ERROR: Resmode "GAUSS" only allows specifying a single "fwhmlw" keyword' & STOP
  endif
endif

if strpos(strupcase(resmode),'FITS') ne -1 then begin          ; If the input for resmode is a FITS file
  hlw = headfits(resmode)
  if strpos(strupcase(sxpar(hlw,'TELESCOP')),'SPITZER') ne -1 then begin
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'IRAC') ne -1 and sxpar(hlw,'CHNLNUM') eq 1 then mode='I1'
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'IRAC') ne -1 and sxpar(hlw,'CHNLNUM') eq 2 then mode='I2'
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'IRAC') ne -1 and sxpar(hlw,'CHNLNUM') eq 3 then mode='I3'
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'IRAC') ne -1 and sxpar(hlw,'CHNLNUM') eq 4 then mode='I4'
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'MIPS') ne -1 and sxpar(hlw,'CHNLNUM') eq 1 then mode='M1'
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'MIPS') ne -1 and sxpar(hlw,'CHNLNUM') eq 2 then begin
      if strpos(strupcase(sxpar(hlw,'EXPTYPE')),'PHT') ne -1 and sxpar(hlw,'FOVID') gt 117 then mode='M2n' else mode='M2w'
    endif
    if strpos(strupcase(sxpar(hlw,'INSTRUME')),'MIPS') ne -1 and sxpar(hlw,'CHNLNUM') eq 3 then mode='M3'
    print, 'mode = '+mode
  endif
  print, strupcase(sxpar(hlw,'PSF'))
  if strtrim(strupcase(sxpar(hlw,'PSF')),2) eq 'GAUSS' then begin
    mode='GAUSS'
    fwhmlw = fxpar(hlw,'FWHM')         ; Gaussian FWHM in arcseconds
    print, 'mode = '+mode+' and FWHM = '+strtrim(fwhmlw,2)+' arcseconds'
  endif
  if strupcase(sxpar(hlw,'PSF')) eq 'GAUSS_ELLIPTICAL' then begin
    mode='GAUSS_ELLIPTICAL'
    bmajlw = 3600.*fxpar(hlw,'BMAJ')   ; Major Axis in arcseconds 
    bminlw = 3600.*fxpar(hlw,'BMIN')   ; Minor Axis in arcseconds
    bpalw = fxpar(hlw,'BPA')           ; Beam Position Angle in degrees
    bpalw_rad = bpalw*!dtor            ; Beam Position Angle in radians
    print, 'mode = GAUSS_ELLIPTICAL, BMAJ = '+strtrim(bmajlw,2)+', BMIN = '+strtrim(bminlw,2)+', BPA = '+strtrim(bpalw,2)
  endif
  if strupcase(strtrim(sxpar(hlw,'PSF'),2)) eq 'IDL_VAR' then begin
    mode='IDL_VAR'
    psflw_file = fxpar(hlw,'PSF_FILE')
    ;-- restore from header the idl variable that contains the normalized psf 
    restore,psflw_file
    psflw_header = psf
    splw = size(psflw_header)
  endif
  if strupcase(strtrim(sxpar(hlw,'PSF'),2)) eq 'FITS' then begin
    mode='FITS'
    psflw_header = readfits(fxpar(hlw,'PSF_FILE'),/SILENT)
    splw = size(psflw_header)
  endif
endif

print, 'RESMODE: '+mode


;---- Reads Lower Resolution Image Pixel Scale ----------
pixsclw = 0.000

if (mode eq 'GAUSS') then pixsclw = fwhmlw/5. else begin
  if (abs(fxpar(hlw,'CD1_1')) gt 0) or (abs(fxpar(hlw,'CD1_2')) gt 0) then pixsclw = 3600.*sqrt(fxpar(hlw,'CD1_1')^2.+fxpar(hlw,'CD1_2')^2.) else begin
    if abs(fxpar(hlw,'CDELT1')) gt 0 then pixsclw = abs(3600.*fxpar(hlw,'CDELT1')) else begin
      ;if (abs(fxpar(hlw,'PXSCAL1')) gt 0) or (fxpar(hlw,'PXSCAL2') gt 0) then pixsclw = fxpar(hlw,'PXSCAL1') else begin
        if fxpar(hlw,'PLTSCALE') gt 0 then pixsclw = fxpar(hlw,'PLTSCALE')
      ;endelse
    endelse
  endelse
endelse

;#########################################################################################################################

window,1 & window,2 & window,3

for i=0L,(nel-1) do begin
  print,'#########################################'
  print,'Image # '+strtrim(i+1,2)+': '+filename[i]

  ;fits_read,filename[i],mosa,hmo
  mosa = mrdfits(filename[i],0,hmo)
  ;if nelunc ne 0 then fits_read,unc_filename[i],unc,hunc
  if nelunc ne 0 then unc = mrdfits(unc_filename[i],0,hunc)
 
  ;---- SUBTRACT BACKGROUND IF /BCKGSUB AND UPDATE HEADER -----------------
  if keyword_set(BCKGSUB) then begin
    sky,mosa,skymode,skysig
    mosa = mosa-skymode
    if nelunc ne 0 then unc = sqrt((unc)^2 + (skysig)^2)
    ; -- update header:
    skym=strtrim(string(skymode),2)
    sxaddhist,'Before Convolution Subtracted Background: '+skym,hmo
    if nelunc ne 0 then sxaddhist,'Before Convolution Subtracted Background: '+skym,hunc
  endif
 
  ;------- CHANGE NaN VALUES TO 0.0 ----------
  problem=where(finite(mosa,/NAN)) 
  if problem[0] gt -1 then mosa[problem]=0.0 

  ;------- CHANGE NaN VALUES TO 2.0 FOR UNCERTAINTIES -----
  if nelunc ne 0 then begin
    problem=where(finite(unc,/NAN)) 
    if problem[0] gt -1 then unc[problem]=2.0    
  endif
     
  sm = size(mosa)
  nx = sm[1]
  ny = sm[2]

  ;-------- CHECKING WHICH INSTRUMENT FOR THE i-th IMAGE ---------------------------------------------------
  if strpos(strupcase(sxpar(hmo,'TELESCOP')),'SPITZER') ne -1 then begin
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'IRAC') ne -1 and sxpar(hmo,'CHNLNUM') eq 1 then instrument='I1'
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'IRAC') ne -1 and sxpar(hmo,'CHNLNUM') eq 2 then instrument='I2'
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'IRAC') ne -1 and sxpar(hmo,'CHNLNUM') eq 3 then instrument='I3'
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'IRAC') ne -1 and sxpar(hmo,'CHNLNUM') eq 4 then instrument='I4'
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'MIPS') ne -1 and sxpar(hmo,'CHNLNUM') eq 1 then instrument='M1'
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'MIPS') ne -1 and sxpar(hmo,'CHNLNUM') eq 2 then begin
      if strpos(strupcase(sxpar(hmo,'EXPTYPE')),'PHT') ne -1 and sxpar(hmo,'FOVID') gt 117 then instrument='M2n' else instrument='M2w'
    endif
    if strpos(strupcase(sxpar(hmo,'INSTRUME')),'MIPS') ne -1 and sxpar(hmo,'CHNLNUM') eq 3 then instrument='M3'
  endif
  if strupcase(strtrim(sxpar(hmo,'PSF'),2)) eq 'GAUSS' then begin
    instrument='GAUSS'
    fwhm = fxpar(hmo,'FWHM')         ; Gaussian FWHM in arcseconds
  endif
  if strupcase(strtrim(sxpar(hmo,'PSF'),2)) eq 'GAUSS_ELLIPTICAL' then begin
    instrument='GAUSS_ELLIPTICAL'
    bmaj = 3.6D3*fxpar(hmo,'BMAJ')   ; Major Axis in arcseconds 
    bmin = 3.6D3*fxpar(hmo,'BMIN')   ; Minor Axis in arcseconds
    bpa = fxpar(hmo,'BPA')           ; Beam Position Angle in degrees
    bpa_rad = bpa*!dtor              ; Beam Position Angle in radians
  endif
  if strupcase(strtrim(sxpar(hmo,'PSF'),2)) eq 'IDL_VAR' then instrument='IDL_VAR'
  if strupcase(strtrim(sxpar(hmo,'PSF'),2)) eq 'FITS' then instrument='FITS'

  
  ;---- Reads i-th Image Pixel Scale ----------
  pixsc = 0.000
  
  if (abs(fxpar(hmo,'CD1_1')) gt 0) or (abs(fxpar(hmo,'CD1_2')) gt 0) then pixsc = 3600.*sqrt(fxpar(hmo,'CD1_1')^2.+fxpar(hmo,'CD1_2')^2.) else begin
    if abs(fxpar(hmo,'CDELT1')) gt 0 then pixsc = abs(3600.*fxpar(hmo,'CDELT1')) else begin
      if fxpar(hmo,'PLTSCALE') gt 0 then pixsc = fxpar(hmo,'PLTSCALE') else begin
        if (abs(fxpar(hmo,'PXSCAL1')) gt 0) or fxpar(hmo,'PXSCAL2') gt 0 then pixsc = abs(fxpar(hmo,'PXSCAL1'))
      endelse
    endelse
  endelse

  if troubleshoot eq 1 then stop
  
  ;-------- CREATING PSF FOR LOWEST RESOLUTION (WITH SAME PIXSCALE OF i-th IMAGE TO CONVOLVE) ---
  CASE mode OF
    'I1' : psflw = stinytim('I1',pixsc,pixsc,gridsize)
    'I2' : psflw = stinytim('I2',pixsc,pixsc,gridsize)
    'I3' : psflw = stinytim('I3',pixsc,pixsc,gridsize)
    'I4' : psflw = stinytim('I4',pixsc,pixsc,gridsize)
    'M1' : psflw = stinytim('M1',pixsc,pixsc,gridsize)
    'M2w' : psflw = stinytim('M2w',pixsc,pixsc,gridsize)
    'M2n' : psflw = stinytim('M2n',pixsc,pixsc,gridsize)
    'M3' : psflw = stinytim('M3',pixsc,pixsc,gridsize)
    'GAUSS' : BEGIN
      sigmalw = fwhmlw/(pixsc*sqrt(8.D0*alog(2.D0)))     ; Gaussian sigma in pixels
      x = dindgen(gridsize) #replicate(1,gridsize)
      y = transpose(x)
      psflw = exp(-0.5D0*((x-gridsize/2)^2 + (y-gridsize/2)^2)/sigmalw^2.)
      psflw = psflw/total(psflw)
    END
    'GAUSS_ELLIPTICAL' : BEGIN
      sigmaxlw = bmajlw/(pixsc*sqrt(8.D0*alog(2.D0)))    ; Major Axis sigma in pixels
      sigmaylw = bminlw/(pixsc*sqrt(8.D0*alog(2.D0)))    ; Minor Axis sigma in pixels
      x = (findgen(gridsize,gridsize) mod gridsize) - gridsize/2  ; makes a grid of the x coord
      y = transpose(x)                                            ; makes a grid of the y coord
      x_r = x*cos(bpalw_rad) + y*sin(bpalw_rad)           ; calculates grid of rotated x
      y_r = -x*sin(bpalw_rad) + y*cos(bpalw_rad)          ; calculates grid of rotated y
      u = (x_r/sigmaxlw)^2. + (y_r/sigmaylw)^2.
      psflw = exp(-0.5D0*u)
      psflw = psflw/total(psflw)
    END
    'IDL_VAR' : BEGIN
      cen = fix(gridsize/2)
      pi = cen - fix(splw[1]/2)
      pf = cen + fix(splw[1]/2)
      psflw=dblarr(gridsize,gridsize)
      psflw[pi:pf,pi:pf] = psflw_header
      psflw=psflw/total(psflw)
    END                  
    'FITS' : BEGIN
      cen = fix(gridsize/2)
      pi = cen - fix(splw[1]/2)
      pf = cen + fix(splw[1]/2)
      psflw=dblarr(gridsize,gridsize)
      psflw[pi:pf,pi:pf] = psflw_header
      psflw=psflw/total(psflw)
    END
  ENDCASE
  
  if troubleshoot eq 1 then stop
  
  ;-------- CHECKING IF GRIDSIZE IS GOOD FOR PSFLW ----------------------------------
  wset,1
  shade_surf,psflw,title='LOWER RESOLUTION PSF',charsize=2.
     
  ;-------- CREATING PSF FOR i-TH IMAGE ---------------------- 
  CASE instrument OF
    'I1'  : psf = stinytim('I1', pixsc,pixsc,gridsize)
    'I2'  : psf = stinytim('I2', pixsc,pixsc,gridsize)
    'I3'  : psf = stinytim('I3', pixsc,pixsc,gridsize)   
    'I4'  : psf = stinytim('I4', pixsc,pixsc,gridsize)
    'M1'  : psf = stinytim('M1', pixsc,pixsc,gridsize)
    'M2n' : psf = stinytim('M2n',pixsc,pixsc,gridsize)
    'M2w' : psf = stinytim('M2w',pixsc,pixsc,gridsize)
    'M3'  : psf = stinytim('M3', pixsc,pixsc,gridsize)
    'GAUSS' : BEGIN
      fwhm = double(fxpar(hmo,'FWHM'))
      fwhmpix = fwhm/pixsc
      sigma = fwhmpix /sqrt(8.D0*alog(2.D0))
      x = dindgen(gridsize) #replicate(1,gridsize)
      y = transpose(x)
      psf = exp(-0.5D0*((x-gridsize/2)^2 + (y-gridsize/2)^2)/sigma^2.)
      psf = psf/total(psf)
    END
    'GAUSS_ELLIPTICAL' : BEGIN
      ;---- PSF = ELLIPTICAL GAUSSIAN ---------------
      fwhmy = fxpar(hmo,'BMAJ')*3.6D3
      fwhmx = fxpar(hmo,'BMIN')*3.6D3
      theta = fxpar(hmo,'BPA')  ; theta in degrees from North to East (counterclockwise), with 0 degrees pointing up
      theta_rad= theta*!dtor    ; theta in radians
      fwhmx_pix = fwhmx/pixsc
      fwhmy_pix = fwhmy/pixsc
      sigmax = double(fwhmx_pix)/sqrt(8.D0*alog(2.D0))
      sigmay = double(fwhmy_pix)/sqrt(8.D0*alog(2.D0))
      x = (findgen(gridsize,gridsize) mod gridsize) - gridsize/2  ; make a grid of the x coord
      y = transpose(x)                                            ; make a grid of the y coord
      xp = x*cos(theta_rad) + y*sin(theta_rad)           ; calculate grid of rotated x
      yp = -x*sin(theta_rad) + y*cos(theta_rad)          ; calculate grid of rotated y
      u = (xp/sigmax)^2. + (yp/sigmay)^2.
      psf = exp(-0.5D0*u)
      psf = psf/total(psf)                            
    END      
    'IDL_VAR' : BEGIN
      psf_file = fxpar(hmo,'PSF_FILE')
      ;-- restore from header the idl variable that contains the normalized psf 
      restore,psf_file
      psf_header = psf
      sp = size(psf_header)
      cen = fix(gridsize/2)
      pi = cen- fix(sp[1]/2)
      pf = cen+ fix(sp[1]/2)
      psf=dblarr(gridsize,gridsize)
      psf[pi:pf,pi:pf] = psf_header
      psf=psf/total(psf)    
    END                  
    'FITS' : BEGIN
      psfgal = readfits(fxpar(hmo,'PSF_FILE'),/SILENT)
      sp = size(psfgal)
      cen = fix(gridsize/2)
      pi = cen- fix(sp[1]/2)
      pf = cen+ fix(sp[1]/2)
      psf=dblarr(gridsize,gridsize)
      psf[pi:pf,pi:pf] = psfgal
      psf=psf/total(psf)    
    END
  ENDCASE

  ;----- DISPLAY THE PSF IMAGE TO CHECK IF GRIDSIZE IS OK FOR PSF IMAGE ----------------------
  wset,2
  shade_surf,psf,title=string('PSF # '+strtrim(i+1,2)),charsize=2.

  ;----- CREATING CONVOLUTION KERNEL ---------------------
  xx = fft(psflw)
  probx = where(xx lt 1.D-5*max(xx))
  ;probx = where(xx lt 1.D-4*max(xx))
  xx[probx] = 0.0
  yy = fft(psf)
  h0 = where(yy eq 0,n0,complement=h1)
  if (n0 ge 1) then yy[h0] = min(abs(yy[h1]))
  kern1 = shift(double(fft(xx/yy,/inverse)),gridsize/2,gridsize/2)

  kk=dblarr(nx,ny)
  if ((gridsize le nx) and (gridsize le ny)) then begin
    kk[0:gridsize-1,0:gridsize-1] = kern1
    kern = shift(kk,-gridsize/2,-gridsize/2)
  endif
  if ((gridsize gt nx) and (gridsize le ny)) then begin
    kk[*,0:(gridsize-1)] = kern1[(gridsize/2-Nx/2):(gridsize/2+Nx-Nx/2-1),*]
    kern = shift(kk,-Nx/2,-gridsize/2)
  endif
  if ((gridsize le nx) and (gridsize gt ny)) then begin
    kk[0:(gridsize-1),*] = kern1[*,(gridsize/2-Ny/2):(gridsize/2+Ny-Ny/2-1)]
    kern = shift(kk,-gridsize/2,-Ny/2)
  endif
  if ((gridsize gt nx) and (gridsize gt ny)) then begin
    kk = kern1[gridsize/2-Nx/2:gridsize/2+Nx-Nx/2-1,gridsize/2-Ny/2:gridsize/2+Ny-Ny/2-1]
    kern = shift(kk,-Nx/2,-Ny/2)
  endif

  ;----- NORMALIZE CONVOLUTION KERNEL -------
  kern = kern / total(kern)

  ;---- DISPLAY THE CONVULUTION KERNEL TO CHECK IF IT LOOKS OK ----------------------
  wset,3
  shade_surf,kern,title='KERNEL',charsize=2.

  ;----- SOME CONSISTENCY CHECKS -------
  eff_FWHM_psf = sqrt(4.*alog(2.d)/!dpi*total(psf)/max(psf))       ;effective FWHM of the psf
  eff_FWHM_psflw = sqrt(4.*alog(2.d)/!dpi*total(psflw)/max(psflw)) ;effective FWHM of the psflw
  eff_FWHM_kern = sqrt(4.*alog(2.d)/!dpi*total(kern)/max(kern))    ;effective FWHM of the convolution kernel
  print,'Image # '+strtrim(i+1,2)+' pixel scale: '+strtrim(pixsc,2)+' arcsec.'
  print,'Image # '+strtrim(i+1,2)+' effective FWHM: '+strtrim(eff_FWHM_psf,2)+' pixels ('+strtrim(eff_FWHM_psf*pixsc,2)+' arcsec)'
  print,'Final resolution effective FWHM: '+strtrim(eff_FWHM_psflw,2)+' pixels ('+strtrim(eff_FWHM_psflw*pixsc,2)+' arcsec)'
  print,'Convolution grid: '+strtrim(gridsize,2)+' x '+strtrim(gridsize,2)+' pixels ('+strtrim(gridsize*float(pixsc),2)+' x '+strtrim(gridsize*float(pixsc),2)+' arcsec)'
  if (eff_FWHM_psf gt eff_FWHM_psflw) then print,'!!!WARNING!!! Effective FWHM of final PSF smaller than original PSF. DECONVOLUTION TAKING PLACE???
  if (eff_FWHM_psf lt 2.)             then print,'!!!WARNING!!! Image # '+strtrim(i+1,2)+' originally undersampled, FWHM < 2 pixels
  if (gridsize lt 10.*eff_FWHM_psflw) then print,'!!!WARNING!!! Consider making gridsize larger, gridsize < 8. x FWHM(psflw)

;  ; sets the default for trunc if no value provided (default = 500)
;  if (n_elements(trunc) ne 1) then trunc = max(mosa)
;  ;---- TRUNCATE WHERE IMAGE VALUE GT TRUNC or LESS THAN -10 (one image had a high neg. value) also unc (DEFAULT=500) --------
;  mosa_before=mosa
;   tr = where(mosa gt trunc)
;   if tr[0] ge 0 then mosa[tr] = trunc
;   tr2 = where(mosa lt -10.)
;   if tr2[0] ge 0 then mosa[tr2] = -1.0
;   if nelunc ne 0 then begin
;     tru = where(unc gt trunc)
;     if tru[0] ge 0 then unc[tru] = trunc
;     tru2 = where(unc lt -10.)
;     if tru2[0] ge 0 then unc[tru2] = 1.0
;   endif

  ;------- CONVOLVE IMAGE WITH KERNEL ---------
  mosaconv = fft(fft(mosa)*fft(kern),/inverse)*nx*ny
  mosaconvR = double(mosaconv)

  ;------ SAVE CONVOLVED IMAGE ----------
  sp = strpos(filename[i],'.fits')
  nam = strmid(filename[i],0,sp)
  iname = nam+'_conv'+mode+'.fits
  ; info for updating header
  sxaddhist,'Convolved by a convolution kernel to match resolution of: ',hmo
  htext=strarr(1)
  if mode eq 'I1' then htext='IRAC channel 1 (3.6 microns) resolution. PSF produced by stinytim.'
  if mode eq 'I2' then htext='IRAC channel 2 (4.5 microns) resolution. PSF produced by stinytim.'
  if mode eq 'I3' then htext='IRAC channel 3 (5.8 microns) resolution. PSF produced by stinytim.'
  if mode eq 'I4' then htext='IRAC channel 4 (8.0 microns) resolution. PSF produced by stinytim.'
  if mode eq 'M1' then htext='MIPS channel 1 (24 microns) resolution. PSF produced by stinytim.'
  if mode eq 'M2n' then htext='MIPS channel 2 narrow field (70 microns) resolution. PSF produced by stinytim.'
  if mode eq 'M2w' then htext='MIPS channel 2 wide field (70 microns) resolution. PSF produced by stinytim.'
  if mode eq 'M3' then htext='MIPS channel 3 (160 microns) resolution. PSF produced by stinytim.'
  if mode eq 'GAUSS' then htext='GAUSSIAN PSF with FWHM of '+strtrim(FWHMlw,2)+' arcseconds
  if mode eq 'GAUSS_ELLIPTICAL' then htext='GAUSSIAN PSF with BMAJ of '+strtrim(BMAJlw,2)+' deg, BMIN of '+strtrim(BMINlw,2)+', and BPA of '+strtrim(BPAlw,2)
  sxaddhist,htext,hmo
  writefits,iname,float(mosaconvR),hmo

  ;------- CONVOLVE UNCERTAINTY WITH KERNEL ----
  if nelunc ne 0 then begin
    unc_conv = fft(fft(unc^2)*fft(kern),/inverse)*nx*ny
    ;unc_conv = fft(fft(unc^2)*fft(kern^2),/inverse)*nx*ny
    uncconvR= double(unc_conv)
    unc_convR = sqrt(abs(uncconvR))

    ;----- SAVE CONVOLVED UNCERTAINTY -----
    spu = strpos(unc_filename[i],'.fits')
    namu = strmid(unc_filename[i],0,spu)
    inameu = namu+'_conv'+mode+'.fits'

    ; update header:
    sxaddhist,htext,hunc
    writefits,inameu,float(unc_convR),hunc
  endif

  if keyword_set(REPROJECT) then begin
    ;----- ALIGN ASTROMETRY AND PIXEL SIZE ACCORDING TO HEADMASTER ------
    hastrom,mosaconvR,hmo,mosaok,hok,headmaster,miss=0.0
      ; update header:
      hist2=strarr(2)
      hist2[0]='Aligned with hastrom to match astrometry of reference image: '
      hist2[1]= imageref
      sxaddhist,hist2,hok
      writefits,nam+'_conv'+mode+'ok.fits',float(mosaok),hok  

    ;----- ALIGN UNCERTAINTIES ACCORDING TO HEADMASTER -----------
    if (nelunc ne 0) then begin
      hastrom,unc_convR^2,hunc,uncok,huncok,headmaster,miss=10.0
        ; update header:
        sxaddhist,hist2,huncok
        writefits,namu+'_conv'+mode+'ok.fits',sqrt(float(uncok)),huncok
    endif   
  endif
;stop ;;; FOR DIAGNOSTIC PURPOSES
endfor 

end
