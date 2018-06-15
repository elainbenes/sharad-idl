; makebrowse.pro

; Wrapper for rdr2san.pro to create SHARAD PDS browse images and ensure that
; RANGE_SHIFT adjustments for multi-mode observations are properly aligned by
; searching for the maximum RANGE_SHIFT value in all modes.
 
; Normally, run in IDL as "makebrowse,168901".  Optional parameters IDIR and 
; ODIR can be used to redirect input and output, INTERVAL defaults to standard
; value of 16 for PDS, and the DIAG flag will be passed to rdr2san.pro.
; 
; Written 2008 Jun 30 by Than Putzig.
; Updated 2008 Jul 18 by Than Putzig - remove "_ss19_700_a" from modes search
;	  string (these items may vary with different observations).
; Updated 2008 Dec 11 by Anthony Egan - change default INTERVAL from 32 to 16
;
; Calls rdr2san.pro
;
pro makebrowse,observation,idir=idir,odir=odir,interval=interval,diag=diag

; Process input parameters
obs=string(observation,format='(I07)'); 	; convert obs to a string
idef='c:\sharad\rdr'				; default input directory
odef='c:\sharad\browse\rdr'+obs 		; default output directory
if n_elements(idir) eq 0  then idir=idef	; set input directory
if n_elements(odir) eq 0  then odir=odef	; set output directory
if file_test(idir)  eq 0 then idir='.'		; use CWD if idir not found
if file_test(odir)  eq 0 then odir='.'		; use CWD if odir not found
idir=file_search(idir,/mark_dir)		; ensure trailing separator
odir=file_search(odir,/mark_dir)		; ensure trailing separator
if n_elements(interval) eq 0 then interval=16	; set default frame interval
diag=keyword_set(diag)				; set diagnostic flag

; Find input modes
modes=idir+'r_'+obs+'_*.dat' 			; modes search string 
if diag then print, 'Searching for files '+modes
file_name=file_search(modes,/fold_case) 	; get list of mode file names
nf=n_elements(file_name)			; count number of files
if file_name[0] eq '' then begin
	print, 'Files '+modes+' not found; exiting.'
	return
	stop
endif

; Determine maximum RANGE_SHIFT values for each mode
maxshift=intarr(nf) 					; initialize maxshift
for n=0,nf-1 do begin
	if diag then print, 'Extracting RANGE_SHIFTs from '+file_name[n]
	rdr2san,file_name[n],range_shift,diag=diag 	; extract RANGE_SHIFTs
	maxshift[n]=max(range_shift) 			; determine maximum
endfor
; set adjustment shift to negative of maximum shift in all modes
shift=-1*max(maxshift) 

; Create the browse images
for n=0,nf-1 do begin
	if diag then print, 'Running RDR2SAN on '+file_name[n]+$
		' with SHIFT of '+strcompress(shift,/rem)
	rdr2san,file_name[n],shift=shift,interval=interval,dir=odir,diag=diag
endfor
end
