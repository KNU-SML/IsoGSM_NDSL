      subroutine msgwrt(lunit,mbay,mbyt)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
      common /bufrmg/ msglen,msgtxt(5000)

      character*4 bufr,sevn
      dimension   mbay(*)

      data bufr/'BUFR'/
      data sevn/'7777'/
      data nmcd/   99 /

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  make sure all sections have even number of bytes
c  ------------------------------------------------

      iad1 = 8
      len1 = iupb(mbay,iad1+1,24)
      len2 = iupb(mbay,iad1+8, 1)
      mtyp = iupb(mbay,iad1+9, 8)
      iad2 = iad1+len1
      len2 = iupb(mbay,iad2+1,24)*len2
      iad3 = iad2+len2
      len3 = iupb(mbay,iad3+1,24)
      iad4 = iad3+len3
      len4 = iupb(mbay,iad4+1,24)

      if(mod(len1,2).ne.0) goto 901
      if(mod(len2,2).ne.0) goto 902
      if(mod(len3,2).ne.0) goto 903
      if(mod(len4,2).ne.0) then
         ibit = iad4*8
         len4 = len4+1
         call pkb(len4,24,mbay,ibit)
         iad5 = iad4+len4
         ibit = iad5*8
         call pkb(0,8,mbay,ibit)
         mbyx = mbyt+1
      else
         mbyx = mbyt
      endif

c  write section 0, nmc code, and section 5
c  ----------------------------------------

      ibit = 0
      jbit = 12*8
      kbit = (mbyx-4)*8

      call pkc(bufr, 4,mbay,ibit)
      call pkb(mbyx,24,mbay,ibit)
      call pkb(   2, 8,mbay,ibit)
      call pkb(nmcd,16,mbay,jbit)
      call pkc(sevn, 4,mbay,kbit)

c  zero out the extra bytes which will be written
c  ----------------------------------------------

      call pkb(0,nbitw,mbay,kbit)

c  write the message plus padding to a word boundary
c  -------------------------------------------------

      imsg = 8/nbytw
      mwrd = (mbyx/8+1)*imsg
      write(lunit) (mbay(i),i=1,mwrd)
c     print*,'msgwrt - lunit=',lunit,' bytes=',mbyx

c  save a memory copy of this message - no bufr tables though
c  ----------------------------------------------------------

      if(mtyp.ne.11) then
         msglen = mwrd
         do i=1,msglen
         msgtxt(i) = mbay(i)
         enddo
      endif

      return
901   call abt('msgwrt - uneven section 1')
902   call abt('msgwrt - uneven section 2')
903   call abt('msgwrt - uneven section 3')
      end
