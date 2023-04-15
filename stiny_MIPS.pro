function stiny_MIPS, channel, pixel1, pixel2, grid

; This calls stiny2 to create a PSF for MIPS 24, 70 and 160 micron
; It requires pre-existing template files: 
; mips_1.par for 24 micron
; mips_2_narrow.par for 70 micron narrow field
; mips_2_wide.par  for 70 micron wide field
; mips_3.par  for 160 micron

; INPUT
; ch = 'M1' or 'M2n' or 'M2w'  or 'M3' ; MIPS Channel
 
; pixel1 and pixel2 = pixel scale in arcsec 
; grid = model grid size in pixels (default = )
;'/media/Master/NGC6946/' = directory where parameter files are. 

;---------------------------------------------------------------

; Check which channel is requested and read the appropriate parameter file
ch = string(channel)
  if ch eq 'M1' then parfile = strarr(249)
  if ch eq 'M2w' then parfile = strarr(275)
  if ch eq 'M2n' then parfile = strarr(275)
  if ch eq 'M3' then parfile = strarr(146)

openr,1,ch+'.par'
readf,1,parfile
close,1

; modify the parameter file for: pixel scale, grid size
parfile[1] = 'MIPS_psf'+strmid(parfile[1],strpos(parfile[1],' ='))
parfile[4] = string(pixel1)+strmid(parfile[4],strpos(parfile[4],' ='))
parfile[5] = string(pixel2)+strmid(parfile[5],strpos(parfile[5],' ='))
parfile[7] = string(grid)+strmid(parfile[7],strpos(parfile[7],' ='))

; write a new parameter file
openw,1,'list_MIPS.par'
for i=0, n_elements(parfile)-1 do printf,1,parfile[i]
close,1

; run stiny2 on this modified parameter file
spawn,'/media/Master/stinytim/stiny2 list_MIPS.par'; spawn,'stiny2 list_MIPS.par'

;read and return the PSF
psf = readfits('MIPS_psf.fits')

return,psf
end
