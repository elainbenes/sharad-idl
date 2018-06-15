;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;	SHARAD RDR Browse Image Generator and Converter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;	filename: rdr2san.pro
; 	ver 1.0 -- 2007 Jul 18 -- Initial Release (Anthony Egan, MATLAB)
; 	ver 1.1 -- 2007 Aug 20 -- Conversion to IDL (Than Putzig)
; 	ver 1.2 -- 2007 Sep 13 -- Fixed bug when RANGE_SHIFT does not vary; 
;				  added options for binary data output
;				  and location (STN-LAT-LON) output
; 	ver 1.3 -- 2008 Feb 15 -- Added summing of frames when decimating
; 	ver 1.4 -- 2008 Apr 08 -- Added DIR, fixed fullpath parsing (AE,TP)
; 	           2008 Apr 09 -- Fixed basename problem (for windows) (TP)
; 	ver 1.5 -- 2008 Jun 13 -- Parameterized NBS (# binary samples) and
;                                 set nulls to file minimum value
; 	ver 1.6 -- 2008 Jun 20 -- Removed default zeroing of RANGE_SHIFT when
;                                 values are all the same; parameterized
;				  null value and set default to -99.
; 	ver 1.7 -- 2008 Jun 27 -- Reduce image size if exceeds 65500 (JPG max);
; 	                          Added RANGE_SHIFT and SHIFT parameters;
; 	                          Replaced TNO with FNO ('trace' vs. 'frame').
; 	           2008 Jul 28 -- Fix to test for nonzero max shift.
; 	           2008 Jul 29 -- Fix to handle case of all values below noise
; 	ver 1.8 -- 2010 May 10 -- Parameterized pid to allow override 
; 	ver 1.9 -- 2011 Jul 21 -- Set pid to long64 (fix for obs>2147401000)
; 	ver 2.0 -- 2013 Feb 14 -- Modified location file format:
;                                 (I12,I8,F16.6,F16.5) -> (I12,I12,F12.6,F12.6)
;	contact: Anthony Egan
;		 Southwest Research Institute
;		 email: anthony@boulder.swri.edu
;	contact: Than Putzig
;		 Southwest Research Institute
;		 email: nathaniel@putzig.com
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;	SUMMARY
;	IDL program that creates a browse jpeg for SHARAD RDR data. The 
;       following parameters are parsed from each RDR frame.
;
;	RDR Field			Variable(s)
;	-----------------		------------------
;	RANGE_SHIFT			sft 
;	ECHO_SAMPLES_REAL		re 
;	ECHO_SAMPLES_IMAGINARY		im
;	SUB_SC_EAST_LONGITUDE		slon (orbit start)
;					elon (orbit end)
;	SUB_SC_PLANETOCENTRIC_LATITUDE	slat (orbit start)
;					elat (orbit end)
;	
;	The real and imaginary parts of the signal are then joined
;	via a traditional log-scale transform:
;
;	sig = 10 * log10 (re^2 + im^2)
;
;	The sig is sorted numerically and the mean value for the lowest
;	100 sig values is added to a running average noise value.
;
;	The mean noise is subtracted from each pixel value, transforming 
;	the "sig" variable into "san", or "signal above noise."
;
;	san = 10*log10(10.^(0.1*sig)-10^(0.1*noise));
;
;	The san values are scaled so the pixel brightness values
;	correspond linearly to an absolute range of 0 to 40 decibels.
;	Finally, the radargram is shifted to adjust for the RANGE_SHIFT
;	value.  The shift is adjusted depending as follows:
;	Default: 
;		Adjust the values by a constant on a per-file basis
;		such that the maximum shift is zero (otherwise,
;		positive values will cause near-range data to wrap to 
;		the bottom of the radargram or large negative values
;		will leave much padding of NULLs at the top).  
;	SHIFT=value
; 		Add the specified value (down-range shift is negative) to the 
; 		RANGE_SHIFT before applying.  Set this value to 0 to
; 	        apply actual RANGE_SHIFT values.
;	/NOSHIFT 
; 		Specify this keyword to turn off all shifts.
;
;	If the RANGE_SHIFT parameter is given, the RANGE_SHIFT values 
;	will be passed back to the calling routine -- useful for determining 
;	a SHIFT value to use for multiple-mode observations.  Use of 
;       this parameter will bypass actual conversion and other output.
;
;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;       IDL version only:
; 	Parameterized the maximum db level (defaults to 40 dB).
;       Optionally, the output JPEG may be compressed along track by
;       specifying an output frame interval larger than 1. Additionally,
;       the adjustment for the RANGE_SHIFT may be bypassed.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
; define procedure and input variables
PRO RDR2SAN,RGRAM,RANGE_SHIFT,MAXDB=MAXDB,INTERVAL=INTERVAL,NOSHIFT=NOSHIFT,$
    DIAG=DIAG,BINARY=BINARY,NBS=NBS,LOCATION=LOCATION,FNO=FNO,DIR=DIR,$
    PID=PID,SHIFT=SHIFT,NULL=NULL,HELP=HELP
prog="RDR2SAN"
prom='['+prog+']: '
vers='2.0'	; set version number

; give syntax if no parameters
IF N_PARAMS() LT 1 OR KEYWORD_SET(HELP) THEN BEGIN     
   PRINT, prog + ' version ',vers
   PRINT,''
   PRINT, prog +   ',RGRAM,RANGE_SHIFT,MAXDB=MAXDB,INTERVAL=INTERVAL,FNO=FNO,'
   PRINT,'         PID=PID,DIR=DIR,NBS=NBS,NULL=NULL,SHIFT=SHIFT,/NOSHIFT,'
   PRINT,'         /BINARY,/LOCATION,/DIAG,/HELP'
   PRINT,' where:  RGRAM       is filename of RDR *.dat file (paths okay)'
   PRINT,'         RANGE_SHIFT is returned array of RANGE_SHIFT values {}'
   PRINT,'         MAXDB       is desired maximum level in dB {40}'
   PRINT,'         INTERVAL    is output frame interval {1}'
   PRINT,'         FNO         is frame number offset for location file {0}'
   PRINT,'         PID         is product id {taken from RGRAM}'
   PRINT,'         DIR         is output directory path {.}'
   PRINT,'         NBS         number of samples/frame in binary output {1000}'
   PRINT,'         NULL        value to use for null data, in dB {-99.}'
   PRINT,'         SHIFT       additional range shift to apply, in samples {}'
   PRINT,'         /NOSHIFT    switch to turn off range shifts'
   PRINT,'         /BINARY     switch to save a binary data file'
   PRINT,'         /LOCATION   switch to save a location file'
   PRINT,'         /DIAG       switch for diagnostic output'
   PRINT,'         /HELP       prints this message'
   PRINT,''
   RETURN
ENDIF    

if n_elements(maxdb) eq 0 then maxdb=40			; maximum level in dB
if n_elements(interval) eq 0 then interval=1		; output frame interval
if n_elements(fno) eq 0 then fno=0			; frame number offset
if n_elements(dir) eq 0 then dir='.'			; set default path
dir=file_search(dir,/mark_dir)
binary=(keyword_set(binary)) ? binary : 0		; binary data switch
if n_elements(nbs) eq 0 then nbs=1000 else binary=1	; default # samples
null=(n_elements(null) eq 0) ? -99. : float(null)	; default null value
location=(keyword_set(location)) ? location : 0		; location file switch
noshift=(keyword_set(noshift)) ? noshift : 0		; RANGE_SHIFT switch
diag=(keyword_set(diag)) ? diag : 0			; diagnostics switch

; set constants
bpf=5822		; bytes per frame
otl=5637		; byte offset to lon/lat info within a frame
nfs=667			; number of frame samples
hdl=194			; header length in bytes
ors=5560		; byte offset to range shift value

; Open the RDR DAT file.
; example for explicitly calling for one RDR dat file:
; openr,ilun,'RDR0169001/R_0169001_001_SS19_700_A.DAT',/get_lun
openr,ilun,rgram,/get_lun
info=file_info(rgram)	; get file info
ate=info.size/bpf 	; number of frames in RDR

; Read start SUB_SC_EAST_LONGITUDE and SUB_SC_PLANETOCENTRIC_LATITUDE from RDR
sloc=assoc(ilun,dblarr(2),otl)
sll=swap_endian(sloc[0],/swap_if_big)
slon=string(round(sll[0]*10),format='(I04)')
slat=string(round(sll[1]*10),format='(I+04)')

; Read end SUB_SC_EAST_LONGITUDE and SUB_SC_PLANETOCENTRIC_LATITUDE from RDR
eloc=assoc(ilun,dblarr(2),info.size-bpf+otl)
ell=swap_endian(eloc[0],/swap_if_big)
elon=string(round(ell[0]*10),format='(I04)')
elat=string(round(ell[1]*10),format='(I+04)')

if diag then begin
	print,'Number of frames is   ',ate
	print,'Start lon and lat are ',sll
	print,'End   lon and lat are ',ell
endif

; Loop through RDR DAT and stream signal (sig) and RANGE_SHIFT (sft) 
; to variables.  Also calculate mean running noise from 1st 100 samples.
; Bypass data read if RANGE_SHIFT is requested by calling program.

prs=arg_present(RANGE_SHIFT)	; check whether RANGE_SHIFT given by caller

sft=intarr(ate)
if prs eq 0 then begin
   noise=0
   sig=fltarr(ate,nfs)
   lat=dblarr(ate)
   lon=lat
endif
for at=long64(0),ate-1 do begin
    isft=assoc(ilun,intarr(1),at*bpf+ors)
    sft[at]=swap_endian(isft[0],/swap_if_big)
    if prs eq 0 then begin
       re=assoc(ilun,fltarr(nfs),at*bpf+hdl)
       im=assoc(ilun,fltarr(nfs),at*bpf+hdl+nfs*4)
       sig[at,*]=10.*alog10(swap_endian(re[0],/swap_if_big)^2+$
	                    swap_endian(im[0],/swap_if_big)^2)
       bnoise=sig(at,[sort(sig[at,*])])
       bnoise=mean(bnoise(0:99))
       if diag and (at-1) mod 1000 eq 0 then $
	   print,'At input frame ',at-1,' bnoise is ', bnoise
       noise=noise+bnoise/ate
       if location then begin
          loc=assoc(ilun,dblarr(2),at*bpf+otl)
          ll=swap_endian(loc[0],/swap_if_big)
          lon[at]=ll[0]
          lat[at]=ll[1]
       endif
    endif
endfor

if diag then begin
	if prs eq 0 then begin
       		print,'Successfully read in the array'
		print,'Average background noise is ',noise
		print,'Minimum raw signal value is ',min(sig)
		print,'Maximum raw signal value is ',max(sig)
	endif
        print,'RANGE_SHIFT values in this file:'
	print, string(sft[uniq(sft)],format='(I6)')
endif

free_lun,ilun
; bail out if return of RANGE_SHIFT values is called for
if prs then begin
	RANGE_SHIFT=sft
	RETURN
endif

; Subtract the noise from each sample. 
; The "sig" is now considered "san", as in "signal above noise".
;sig=sig-noise					; too 'binary'
;sig=10.*alog10(10.^(0.1*sig)-10.^(0.1*noise))	; produces NaNs
sig=10.^(0.1*sig)-10.^(0.1*noise)
wsn=where(sig le 0)
wsd=where(sig gt 0)
if (wsd[0] eq -1) then print, 'WARNING: all data values below noise level'
if wsn[0] ne -1 then sig[wsn]=10.^(null/10.)
sig=10.*alog10(sig)
smn=(wsd[0] ne -1)?min(sig[wsd]):min(sig)	; find min san

if diag then begin
	print, 'Minimum signal-above-noise is ', smn
	print, 'Maximum signal-above-noise is ', max(sig)
endif

; build output file name 
noise=strcompress(round(noise*10),/rem)
bname=file_basename(rgram)
outn=strmid(bname,0,strlen(bname)-6)
if n_elements(pid) eq 0 then pid=long64(strmid(bname,2,7))*1000

; Uncomment one line to include additional information with the filename
;; outname=dir+outn+'_'+noise+'_'+slon+'_'+slat+'_'+elon+'_'+elat
;; outname=dir+outn+'_'+noise+'_'+slon+'_'+slat
outname=dir+outn

; This section was below conversion to byte values; moved before it to allow
; range shift correction prior to saving binary data file
; check RANGE_SHIFT switch
nsft=n_elements(SHIFT)		; check whether SHIFT is specified
if nsft eq 0 then SHIFT=0	; set SHIFT value to zero if not given
if min(sft) eq max(sft) and nsft eq 0 then begin
   noshift=1
   if diag then print,'RANGE_SHIFT values are all the same; no shift applied.'
endif
if noshift eq 0 then begin
   if diag then print,'Beginning RANGE_SHIFT scaling'
   ; transform RANGE_SHIFT into pixels
   ; also tare (zero) against the maximum shift in file or SHIFT value if given
   sft=round(0.5*(sft-max(sft)*(1-nsft)+SHIFT))
   szs=size(sig)
   if max(abs(sft)) gt 0 then begin	; apply sft when max is nonzero
      sig=[[sig],[fltarr(szs[1],max(abs(sft)))+null]]
      for i=0l,n_elements(sft)-1 do begin
          sig[i,*]=shift(sig[i,*],0,0-sft[i]);
      endfor
   endif
   if diag then begin
	   print,'Pixel shift values are:'
	   print,string(sft[uniq(sft)],format='(I6)')
   endif
endif

szs=size(sig)		; get present size of array
nf=szs[1]		; number of frames
ns=szs[2]		; number of samples (after range shifts)

; average over intervals and decimate
if interval gt 1 then begin
   msg=nf mod interval
   ;;sig=reform(sig[0:nf-msg-1,*],nf/interval,interval,ns,/overwrite)
   ;;sig=total(sig,2)/interval
   sig=rebin(sig[0:nf-msg-1,*],nf/interval,ns)
   szs=size(sig)	; get new size of array
   nf=szs[1]		; number of frames
   ns=szs[2]		; number of samples (after range shifts)
endif

if diag then begin
	print, 'number of samples is ', ns
	print, 'number of frames is ', nf
endif

; write binary data file if requested
if binary then begin
   while nbs lt ns do nbs=nbs+nbs
   if nbs gt 1000 then print,'NUMBER OF DATA SAMPLES INCREASED TO '+$
      strcompress(nbs,/rem) +' TO ACCOMMODATE RANGE SHIFTS'
   sig=[[sig],[fltarr(nf,nbs-ns)+null]]
   openw,olun,outname+'.dat',/get_lun
   a=assoc(olun,sig)
   a[0]=sig
   free_lun,olun
   if diag then begin
      print,'Binary array size is ', $
		strcompress(nf,/rem), ' X ', strcompress(nbs,/rem)
      print,'Binary file written as ',outname+'.dat'
   endif
endif

; Scale linearly so pixel value 255 corresponds to a MAXDB dB signal 
; above the mean noise
;sig=1*sig./40;				; original MATLAB code
sig=sig/MAXDB*255
wsl=where(sig lt 0)
if wsl[0] ne -1 then sig[wsl]=0.	; zero out values below noise floor
wsh=where(sig gt 255)
if wsh[0] ne -1 then sig[wsh]=255.	; clip values greater than MAXDB
sig=byte(sig)

if diag then begin
	print,'Minimum byte-scaled signal-above-noise is ',min(sig)
	print,'Maximum byte-scaled signal-above-noise is ',max(sig)
endif

; write output image
mxpx=65500				; max # of pixels for JPEGs
if nf gt mxpx then begin
   sig=congrid(sig,mxpx,nbs)
   print, "WARNING: reduced image frame dimension from ",strcompress(nf,/rem), $
          " to ",strcompress(mxpx,/rem)," (max dimension for JPEGs)."
endif
outjpg=outname+'.jpg'
write_jpeg,outjpg,rotate(sig,7)

if diag then begin
	szs=size(sig)	; get new size of array
	print,'Image array size is ', $
		strcompress(szs[1],/rem), ' X ', strcompress(szs[2],/rem)
	print,'Image file written as ',outjpg
endif

if location then begin
   openw,olun,outname+'.latlon',/get_lun
   for at=long64(0),ate/interval-1 do begin
       atl=at*interval+(interval-1)/2
       printf,ilun,pid,at+1+fno,lat[atl],lon[atl],format='(I12,I12,F12.6,F12.6)'
   endfor
   free_lun,olun
   if diag then print,'Location file written as '+outname+'.latlon'
endif

end
