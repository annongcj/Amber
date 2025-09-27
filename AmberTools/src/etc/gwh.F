      program     GuesWatH
c
c**********************************************************************c
c                                                                      c
c    Program : GuesWatH - create hydrogen positions of water           c
c    Author : Shuichi Miyamoto                                         c
c    Modified: D.A. Case                                               c
c    Date of last update : February 8, 1994                            c
c                                                                      c
c**********************************************************************c
c
      parameter   (maxsph=15000,maxcir=1000)
      parameter   (maxatm=100000,maxwat=3000,maxwa3=9000)
      character*80 title,prmtop,watpdb,arg
      character*4  atnam(maxatm),atnamw(maxatm),rsnamt(maxatm) 
      character*3  rsnam(maxatm),rsnamw(maxatm)
      character*1  chid(maxatm) 
      integer*4   resno(maxatm) 
      dimension   msph(3,maxsph), ncirc(maxsph),icirc(maxcir,maxsph)
      dimension   lp(3,2), lpp(2,maxcir,maxsph)
      dimension   cord(3,maxatm), kres(maxatm)
     &          , chrg(maxatm), elecp(maxsph,maxwat)
     &                        , elecw(maxsph,maxwat)
      dimension   occ(maxwat),temp(maxwat)
      dimension   cordo(3,maxwat),cordh(3,2,maxwat),cordw(3,maxwa3)
      dimension   dmin(maxwat),emin(maxwat),isort(maxwat)
      dimension   iterm(maxwa3)
      dimension   rh(2), hnew(3,2)
      dimension   vech1(3),vech2(3),vecx(3), unset(3,3)
c
      data  oh /0.9572/ , ohcos /0.586077/, hh2 /0.7568/
     &     ,chrg_ow /-15.1974/ , chrg_hw /7.59870/
c
c
      prmtop = 'prmtop'
      watpdb = 'watpdb'
      iddd = 0
      iord = 0
      nargs = iargc()
      iarg = 0
      if (nargs.eq.iarg) go to 20
   10 iarg = iarg + 1
          call getarg(iarg,arg)
          if (arg.eq.'-p') then
            iarg = iarg + 1
            call getarg(iarg,prmtop)
          else if (arg.eq.'-w') then
            iarg = iarg + 1
            call getarg(iarg,watpdb)
          else if (arg.eq.'-consd') then
            iddd = 1
          else if (arg.eq.'-elstat') then
            iord = 1
          else if (arg.eq.'-h' .or. arg.eq.'--help') then
            write(0,*) 'Usage: gwh [-p prmtop ] [-w watpdb ]',
     .       ' [-consd] [-elstat] < input-pdb-file > output-pdb-file'
      write(0,*) '  watpdb contains new waters to have hydrogens fixed'
      write(0,*) '  -consd: use a constant dielectric'
      write(0,*) '  -elstat: sort on electrostatics; otherwise sort on'
      write(0,*) '           minimum H-bond distance'
      write(0,*) 'Note: waters in input-pdb-file will not be modified'
            call mexit(0,1)
          else
            write(0,*) 'Bad argument: ',iarg,arg
            write(0,*) 'Usage: gwh [-p prmtop ] [-w watpdb ]',
     .       ' [-consd] [-elstat] < input-pdb-file > output-pdb-file'
      write(0,*) '  watpdb contains new waters to have hydrogens fixed'
      write(0,*) '  -consd: use a constant dielectric'
      write(0,*) '  -elstat: sort on electrostatics; otherwise sort on'
      write(0,*) '           minimum H-bond distance'
      write(0,*) 'Note: waters in input-pdb-file will not be modified'
            call mexit(0,1)
          end if
      if (iarg.lt.nargs) go to 10
   20 continue
      call amopen(1,prmtop,'O', 'F', 'R')
      call top_3 (1, title,natom,atnam,chrg,nres,rsnamt,kres,ifbox)
      close(1)
      call pdbatm (maxatm,5,0,0, natomp,atnam,
     .      rsnam,chid,resno,cord,iterm)
      call amopen(3,watpdb,'O','F','R')
      call pdbhet2 (3,0,0, no,atnamw,rsnamw,cordo,occ,temp)
      close(3)
c
      if (natom.ne.natomp) then
         write  (0,'(a,2i6)') 'diff in natom: ',natom,natomp
         call mexit(0,1)
      endif
c
      mcycle =5
c
c
c                                             ----- Set up of Sphere -----
      weig = 1.0
      nsph = 0
      do i=1,21
         do j=1,21
            do k=1,21
               ijk = (i-11)*(i-11) + (j-11)*(j-11) + (k-11)*(k-11) 
               if (ijk.gt.81 .and. ijk.le.100) then
                  nsph = nsph + 1
                  msph(1,nsph) = i - 11
                  msph(2,nsph) = j - 11
                  msph(3,nsph) = k - 11
               endif
            end do
         end do
      end do
c
      do i=1,nsph-1
         nhit = 0
         do j=i+1,nsph
            kx = (msph(1,i) - msph(1,j)) ** 2
            ky = (msph(2,i) - msph(2,j)) ** 2
            kz = (msph(3,i) - msph(3,j)) ** 2
            kxyz = kx + ky + kz
            if (kxyz.gt.213 .and. kxyz.le.243) then
               nhit = nhit + 1
               icirc(nhit,i) = j
               vech1(1) = float(msph(1,i)) / 10.
               vech1(2) = float(msph(2,i)) / 10.
               vech1(3) = float(msph(3,i)) / 10.
               vech2(1) = float(msph(1,j)) / 10.
               vech2(2) = float(msph(2,j)) / 10.
               vech2(3) = float(msph(3,j)) / 10.
               vecx(1) = vech1(1) + vech2(1)
               vecx(2) = vech1(2) + vech2(2)
               vecx(3) = vech1(3) + vech2(3)
               call uniset (1,vecx,vech2, unset)
               lp(1,1) = - unset(1,1) * 5.80 + unset(1,3) * 8.15
               lp(2,1) = - unset(2,1) * 5.80 + unset(2,3) * 8.15
               lp(3,1) = - unset(3,1) * 5.80 + unset(3,3) * 8.15
               lp(1,2) = - unset(1,1) * 5.80 - unset(1,3) * 8.15
               lp(2,2) = - unset(2,1) * 5.80 - unset(2,3) * 8.15
               lp(3,2) = - unset(3,1) * 5.80 - unset(3,3) * 8.15 
               do l=1,2
                  kmin = 99999
                  do k=1,nsph
                     kx = (lp(1,l) - msph(1,k)) ** 2
                     ky = (lp(2,l) - msph(2,k)) ** 2
                     kz = (lp(3,l) - msph(3,k)) ** 2  
                     if (kmin .gt. kx + ky + kz) then
                         kmin  =   kx + ky + kz
                        lpp(l,nhit,i) = k
                     endif
                  end do
               end do
            endif
         end do
         ncirc(i) = nhit
      end do
c
c                                               --- Calc Dist from Prot ---
      do io=1,no
         dismin = 999.e10
         do k=1,natom
            dist2 = 0.
            do j=1,3
               dist2 = dist2 + (cordo(j,io) - cord(j,k)) ** 2
            end do
            if (dismin.gt.dist2) dismin = dist2
         end do
         dmin(io) = dismin
      end do
c
c                                               --- Calc Elec from Prot ---
c
      do io=1,no
         enemin = 999.e10
         do i=1,nsph
            xs = float(msph(1,i)) / 10. + cordo(1,io)
            ys = float(msph(2,i)) / 10. + cordo(2,io)
            zs = float(msph(3,i)) / 10. + cordo(3,io)
            esum = 0.
            do k=1,natom
               r2 = ( ( xs - cord(1,k) ) * ( xs - cord(1,k) )
     &              + ( ys - cord(2,k) ) * ( ys - cord(2,k) )
     &              + ( zs - cord(3,k) ) * ( zs - cord(3,k) ) )
               if (iddd.eq.0) then
                  esum = esum + chrg(k) / r2
               else
                  esum = esum + chrg(k) / sqrt(r2)
               endif
            end do
            elecp(i,io) = esum
            if (enemin .gt. elecp(i,io) )  enemin = elecp(i,io)
         end do
         emin(io) = enemin
      end do
c
      if (iord.eq.0) then
         call sortf_as (no,dmin, isort)
      else
         call sortf_as (no,emin, isort)
      endif
c
c                                                 ----- Water Loop -----
      do l=1,mcycle
c
      etot = 0.
      rhh2 = 0.
c
      do io=1,no
         iso = isort(io)
         if (io.eq.1) go to 170
c                                                 --- Calc ElecW ---
         do i=1,nsph
            xs = float(msph(1,i)) / 10. + cordo(1,iso)
            ys = float(msph(2,i)) / 10. + cordo(2,iso)
            zs = float(msph(3,i)) / 10. + cordo(3,iso)
            esum = 0.
            do iio=1,no
               if (l.eq.1 .and. iio.eq.io) exit
               if (iio.eq.io) cycle
               iiso = isort(iio)
               if (k.eq.io) cycle
               r2 = ( xs - cordo(1,iiso) ) ** 2
     &            + ( ys - cordo(2,iiso) ) ** 2
     &            + ( zs - cordo(3,iiso) ) ** 2
               if (iddd.eq.0) then
                  esum = esum + chrg_ow / r2
               else
                  esum = esum + chrg_ow / sqrt(r2)
               endif
               do m=1,2
                  rh(m) = ( xs - cordh(1,m,iiso) ) ** 2
     &               + ( ys - cordh(2,m,iiso) ) ** 2
     &               + ( zs - cordh(3,m,iiso) ) ** 2
                  if (iddd.eq.0) then
                     esum = esum + chrg_hw / rh(m)
                  else
                     esum = esum + chrg_hw / sqrt(rh(m))
                  endif
               end do
            end do
c
            elecw(i,iso) = esum
         end do
c
c
  170    continue
c                                                 --- Calc Min ---
         emint = 999.e10
         do i=1,nsph
            if (ncirc(i).eq.0) cycle
            eminj = 999.e10
            do j=1,ncirc(i)
               elp = elecp(lpp(1,j,i),iso) + elecw(lpp(1,j,i),iso)
     &             + elecp(lpp(2,j,i),iso) + elecw(lpp(2,j,i),iso)
               if (eminj.gt.elecp(icirc(j,i),iso)
     &                    + elecw(icirc(j,i),iso) - elp*weig) then
                   eminj =  elecp(icirc(j,i),iso)
     &                    + elecw(icirc(j,i),iso) - elp*weig
                   nminj = icirc(j,i)
               endif
            end do
            if (emint .gt. elecp(i,iso) + elecw(i,iso) + eminj) then
                emint  =   elecp(i,iso) + elecw(i,iso) + eminj
                nmini = i
                nmink = nminj
            endif
         end do
c
         etot = etot + emint
c
         vech1(1) = float(msph(1,nmini)) / 10.
         vech1(2) = float(msph(2,nmini)) / 10.
         vech1(3) = float(msph(3,nmini)) / 10.
         vech2(1) = float(msph(1,nmink)) / 10.
         vech2(2) = float(msph(2,nmink)) / 10.
         vech2(3) = float(msph(3,nmink)) / 10.
         vecx(1) = vech1(1) + vech2(1)
         vecx(2) = vech1(2) + vech2(2)
         vecx(3) = vech1(3) + vech2(3)
         call uniset (1,vecx,vech2, unset)
         do m=1,3
            hnew(m,1) = cordo(m,iso)
     &                + ohcos * unset(m,1) + hh2 * unset(m,2)
            hnew(m,2) = cordo(m,iso)
     &                + ohcos * unset(m,1) - hh2 * unset(m,2)
            rhh2 = rhh2 + ( hnew(m,1) - cordh(m,1,iso) ) ** 2
     &                  + ( hnew(m,2) - cordh(m,2,iso) ) ** 2
            cordh(m,1,iso) = hnew(m,1)
            cordh(m,2,iso) = hnew(m,2)
         end do
      end do
c
c
      avet = etot / float(no)
      rms = sqrt ( 0.5 * rhh2 / float(no) )
c
      write  (0,220) avet,rms
  220 format ('REMARK   Av.E =',f8.3,' rms =',f8.3)
c
      end do
c
c                                                 ----- End of Water Loop -----
c
      iok = 0
      do 260 io=1,no
         iok = iok + 1
         iio  = iok * 3 - 2
         iih1 = iio + 1
         iih2 = iio + 2 
         do 240 m=1,3
            cordw(m,iio) = cordo(m,io)
            cordw(m,iih1) = cordh(m,1,io)
            cordw(m,iih2) = cordh(m,2,io)
  240    continue
         atnamw(iio)  = ' O  '
         atnamw(iih1) = ' H1 '
         atnamw(iih2) = ' H2 '
         do 250 m=1,3
            rsnamw(iio-1+m) = 'HOH'
            chid  (iio-1+m) = ' '
            resno (iio-1+m) = io + nres
  250    continue
         iterm(iio)  = 0
         iterm(iih1) = 0
         iterm(iih2) = 1
  260 continue
c
      natw = 3 * iok
      call pdbatm(maxatm,0,0,6,
     .      natw,atnamw,rsnamw,chid,resno,cordw,iterm)
c
      call mexit(6,0)
      end
c
c======================================================================c
c
c Function that reads and writes PDB files.
c JPH added Feb 97 - also write PDBQ files - used by Autodock.
c
c
      subroutine pdbatm( maxatm, ird, ipr, iwr,
     .                   natom, atmnam, resnam, chid, resno, cord,
     .                   iterm )
c
      integer maxatm, ird, ipr, iwr, natom
      character*4  atmnam(maxatm)
      character*3  resnam(maxatm)
      character*1  chid(maxatm)
      integer      resno(maxatm), iterm(maxatm)
      real         cord(3,maxatm), chrg(1)
c
c  ---  Read  ---
c
      if ( ird .gt. 0 ) 
     .   call pdbin( ird, maxatm, natom, atmnam, resnam, chid,
     .               resno, cord, iterm )
c
c  ---  Print  ---
c
      if ( ipr .gt. 0 )
     .   call pdbout( ipr, 2, maxatm, natom, atmnam, resnam, chid,
     .                resno, cord, iterm, 0, 0, chrg )
c
c  ---  Print  ---
c
      if ( iwr .gt. 0 )
     .   call pdbout( iwr, 1, maxatm, natom, atmnam, resnam, chid,
     .                resno, cord, iterm, 0, 0, chrg )
c
      return
      end
c
c======================================================================cc
      subroutine pdbin( iun, maxatm, natom, atmnam, resnam, chid,
     .                   resno, cord, iterm )
c
c Read a PDB file
c
c
c
      integer iun, maxatm, natom, resno(maxatm), iterm(maxatm)
      character*4  atmnam(maxatm)
      character*3  resnam(maxatm), tag
      character*1  chid(maxatm)
      real         cord(3,maxatm)
      character*80 line
      equivalence  (line,tag)
c
      natom = 0
  100 continue
      read( iun, 200, end=500 ) line
  200 format( a80 )
      if ( tag .eq. 'END' ) go to 500
      if ( tag .eq. 'TER' ) iterm(natom) = 1
      if ( tag .eq. 'HET' .or. tag .eq. 'ATO' ) then
         natom = natom + 1
         if ( natom .gt. maxatm ) then
            write(0,*) 'pdbatm: maxatm exceeded'
            call mexit(0,1)
         end if
         iterm( natom ) = 0
         read( line, 300 ) atmnam(natom), resnam(natom), chid(natom),
     .                     resno(natom), cord(1,natom),
     .                     cord(2,natom), cord(3,natom)
  300    format( 12x, a4, 1x, a3, 1x, a1, i4, 4x, 3f8.0, 16x )
      endif
      go to 100
c
c
  500 continue
      return
      end
c
c======================================================================c
c
      subroutine pdbout( iun, ipos, maxatm, natom, atmnam, resnam, chid,
     .                   resno, cord, iterm, iq, iqr, chrg )
c
c Write a PDB or PDBQ file
c
c
c
      integer iun, ipos, maxatm, natom, iq, iqr, iend,
     .        resno(maxatm), iterm(maxatm), i
      character*4  atmnam(maxatm)
      character*3  resnam(maxatm)
      character*1  chid(maxatm)
      real         cord(3,maxatm), chrg(maxatm), q1
      character*80 line
c
      if (iun .gt. 0 ) then
         iend = ipos + 77
         if ( iq .ne. 0 ) iend = ipos + 75
         do 100 i = 1, natom
            line = ' '
            write( line(ipos:), 200 ) i, atmnam(i), resnam(i), chid(i),
     &         resno(i),cord(1,i), cord(2,i), cord(3,i), atmnam(i)(1:2)
c
c Add the charges for a PDBQ file
c
            if ( iq .ne. 0 ) then
c
c NB Remove conversion factor from PARM
c
               q1 = chrg(i) / 18.2223d0
c
c iqr describes the format for the charges
c     = 0  write charges in F6.3 format - truncating
c       1  write charges in F6.3 format - rounding
c      -1  write charges in F8.5 format
c
               if ( iqr .lt. 0 ) then
                  write( line(ipos+70:), 400 ) q1
               else
                  if ( iqr .gt. 0 ) q1 = q1 + .0005
                  write( line(ipos+70:), 500 ) q1
               end if
            end if
            write( iun, 1000 ) line(1:iend)
            if ( iterm(i) .eq. 1 ) then
               line(ipos:) = 'TER'
               write( iun, 1000 ) line(1:ipos+2)
            endif
  100    continue
         line(ipos:) = 'END'
         write( iun, 1000 ) line(1:ipos+2)
      end if
c
  200 format( 'ATOM  ', i5, 1x, a4, 1x, a3, 1x, a1, i4, 4x, 3f8.3, 
     &        '  1.00  0.00          ', a2 )
  400 format( f8.5 )
  500 format( f6.3 )
 1000 format( a )
c
      return
      end
c======================================================================c
c
      subroutine pdbhet2 (ird,ipr,iwr, natom,atmnam,resnam,cord
     &                                ,occ,temp)
      character*80 line
      character*4  atmnam(1)
      character*3  resnam, tag
      dimension    cord(3,1),occ(1),temp(1)
      equivalence  (line,tag)
c
c                                                ---  Read  ---
      if (ird.le.0) go to 6
      natom = 0
   51 read  (ird,501,end=6) line
  501 format (a80)
      if (tag.eq.'END') go to 6
      if (tag.eq.'HET' .or. tag.eq.'ATO') then
         natom = natom + 1
         read   (line,511) atmnam(natom),resnam,(cord(j,natom),j=1,3)
     &                    ,occ(natom),temp
  511    format (12x,a4,1x,a3,10x,3f8.0,2f6.0)
      endif
      go to 51
c
c                                                ---  Print  ---
    6 if (ipr.le.0) go to 8
      do 610 i=1,natom
         write  (ipr,601) i,atmnam(i),resnam,(cord(j,i),j=1,3),
     &                    occ(i),temp(i)
  601    format (' HETATM',i5,1x,a4,1x,a3,2x,'   1',4x,3f8.3,2f6.2)
  610 continue
c
c                                                ---  Write  ---
    8 if (iwr.le.0) return
      do 810 i=1,natom
         write  (iwr,801) i,atmnam(i),resnam,(cord(j,i),j=1,3),
     &                    occ(i),temp(i)
  801    format ('HETATM',i5,1x,a4,1x,a3,2x,'   1',4x,3f8.3,2f6.2)
  810 continue
      write  (iwr,809)
  809 format ('TER   '/'END   ')
      return
      end
c======================================================================c
c
      subroutine top_3 (intop, title,natom,atnam,chrg,nres,resnam,kres
     &                        ,ifbox)
      character*80 title
      character*20 line
      character*4  atnam(*),resnam(*)
      dimension    kres(*),chrg(*)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(a80)') title
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(2i8)') natom,ntypes
      read   (intop,'(8x,i8)') nres
      read   (intop,'(56x,i8)') ifbox
      read   (intop,*)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(20a4)') (atnam(i),i=1,natom)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(5e16.8)') (chrg(i),i=1,natom)
c
      do i=1,99999
         read   (intop,'(a20)') line
         if( line(1:19) .eq. '%FLAG RESIDUE_LABEL' ) exit
      end do
      read   (intop,'(a20)') line
      read   (intop,'(20a4)') (resnam(i),i=1,nres)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(10i8)') (kres(i),i=1,nres)
c
c
c      read (NF,9108) ITITL
c      read (NF,9118) NTOTAT,NTYPES,NBONH,NBONA,NTHETH,NTHETS,NPHIH,
c     +          NPHIA,JHPARM,JPARM,NEXT,NTOTRS,MBONA,MTHETS,MPHIA,
c     +          MUMBND,MUMANG,MPTRA,NATYP,NHB,IFPERT,NBPER,NGPER,
c     +          NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP
c      NTYPE = NTYPES*NTYPES
c         WRITE(NF,9108) (IGRAPH(I),I = 1,NTOTAT)
c         WRITE(NF,9128) (CHRG(I),  I = 1,NTOTAT)
c         WRITE(NF,9128) (AMASS(I), I = 1,NTOTAT)
c         WRITE(NF,9118) (IAC(I),   I = 1,NTOTAT)
c         WRITE(NF,9118) (NUMEX(I), I = 1,NTOTAT)
c         WRITE(NF,9118) (NNO(I),   I = 1,NTYPE)
c         WRITE(NF,9108) (LABRES(I),I = 1,NTOTRS)
c         WRITE(NF,9118) (IPRES(I), I = 1,NTOTRS)
c         WRITE(NF,9128) (RK(I),    I = 1,MUMBND)
c         WRITE(NF,9128) (REQ(I),   I = 1,MUMBND)
c         WRITE(NF,9128) (TK(I),    I = 1,MUMANG)
c         WRITE(NF,9128) (TEQ(I),   I = 1,MUMANG)
c         WRITE(NF,9128) (PK(I),    I = 1,MPTRA)
c         WRITE(NF,9128) (PN(I),    I = 1,MPTRA)
c         WRITE(NF,9128) (PHASE(I), I = 1,MPTRA)
c         WRITE(NF,9128) (SOLTY(I), I = 1,NATYP)
c         NTTYP = (NTYPES*(NTYPES+1))/2
c         WRITE(NF,9128) (CN1(I),   I = 1,NTTYP)
c         WRITE(NF,9128) (CN2(I),   I = 1,NTTYP)
c 9108 FORMAT(20A4)
c 9118 FORMAT(12I6)
c 9128 FORMAT(5E16.8)
      return
      end
      SUBROUTINE CORTOR (CORD1,CORD2,CORD3,BL,ALPHA,THETA,CORDO)
      DIMENSION CORD1(3),CORD2(3),CORD3(3),CORDO(3),VEC12(3),VEC23(3),
     &          FCORD(3),ECORD(3),UNSETF(3,3)
C----------------------------------------------------------------------C
C                                                                      C
C                                                 Z                    C
C                           1                     !      O   Y         C
C                           C :                   !     /  /           C
C                               :         ALPHA   !   / /              C
C                                 :               ! //                 C
C                                   2 ----------> 3    -------> X      C
C                                   B             A                    C
C                                                                      C
C----------------------------------------------------------------------C
      DATA DTOR /0.01745329/
C
      DO 10 J=1,3
      VEC23(J)=CORD3(J)-CORD2(J)
      VEC12(J)=CORD2(J)-CORD1(J)
   10 CONTINUE
      FCORD(1)=-BL*COS(ALPHA*DTOR)
      FCORD(2)= BL*SIN(ALPHA*DTOR)*COS(THETA*DTOR)
      FCORD(3)= BL*SIN(ALPHA*DTOR)*SIN(THETA*DTOR)
      CALL UNISET (2,VEC23,VEC12,UNSETF)
      CALL FTOE1 (UNSETF,FCORD,ECORD)
      DO 20 J=1,3
      CORDO(J)=CORD3(J)+ECORD(J)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE FTOE1 (UNSETF,FCORD,ECORD)
      DIMENSION UNSETF(3,3),FCORD(3),ECORD(3)
C----------------------------------------------------------------------C
C          THIS SUB. CHANGES BASES  ( E --> F )                        C
C          THAT IS TO SAY, A COORDINATE ON F IS CHANGED TO             C
C                          A COORDINATE ON E .                         C
C          E IS A SET OF THREE FUNDAMENTAL VECTORS (1,0,0),(0,1,0),    C
C                                          (0,0,1)                     C
C----------------------------------------------------------------------C
      DO 10 J=1,3
      ECORD(J)=UNSETF(J,1)*FCORD(1)+UNSETF(J,2)*FCORD(2)
     &        +UNSETF(J,3)*FCORD(3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE UNISET (IND,VECA,VECB,UNSET)
C
      DIMENSION VECA(3),VECB(3),UNSET(3,3),
     &          VECX(3),VECY(3),VECZ(3),UECX(3),UECY(3),UECZ(3)
C
      IPRINT=0
      JND = IND
      IF (IND.LT.0) THEN
         IPRINT=1
         JND=-IND
      ENDIF
      IF (JND.EQ.0 .OR. JND.GE.8) GO TO 999
C
      GO TO (1,2,3,4,5,6,7),JND
C
    1 CALL VECPRD (VECA,VECB,VECZ)
      CALL VECPRD (VECZ,VECA,VECY)
      GO TO 91
    2 CALL VECPRD (VECB,VECA,VECZ)
      CALL VECPRD (VECZ,VECA,VECY)
      GO TO 91
    3 CALL VECPRD (VECB,VECA,VECY)
      CALL VECPRD (VECA,VECY,VECZ)
      GO TO 91
C
    4 CALL VECPRD (VECB,VECA,VECZ)
      CALL VECPRD (VECA,VECZ,VECX)
      GO TO 94
    5 CALL VECPRD (VECA,VECB,VECX)
      CALL VECPRD (VECX,VECA,VECZ)
      GO TO 94
C
    6 CALL VECPRD (VECA,VECB,VECY)
      CALL VECPRD (VECY,VECA,VECX)
      GO TO 96
    7 CALL VECPRD (VECB,VECA,VECX)
      CALL VECPRD (VECA,VECX,VECY)
      GO TO 96
C
   91 CALL UNIVEC (VECA,UECX)
      CALL UNIVEC (VECY,UECY)
      CALL UNIVEC (VECZ,UECZ)
      GO TO 100
   94 CALL UNIVEC (VECX,UECX)
      CALL UNIVEC (VECA,UECY)
      CALL UNIVEC (VECZ,UECZ)
      GO TO 100
   96 CALL UNIVEC (VECX,UECX)
      CALL UNIVEC (VECY,UECY)
      CALL UNIVEC (VECA,UECZ)
C
  100 DO 10 I=1,3
      UNSET(I,1)=UECX(I)
      UNSET(I,2)=UECY(I)
      UNSET(I,3)=UECZ(I)
   10 CONTINUE
      IF (IPRINT.EQ.0) RETURN
C
      WRITE (6,601) IND,((UNSET(J,I),I=1,3),J=1,3)
  601 FORMAT (/1H ,' *** UNSET WILL BE PRINTED BY IND=',I3,' AT UNISET'
     &            ,' ***'
     &       //1H ,4X,'1(X)',4X,'2(Y)',4X,'3(Z)'
     &        /1H ,'X',3F8.3/1H ,'Y',3F8.3/1H ,'Z',3F8.3/)
      RETURN
C
  999 WRITE(6,699) IND
  699 FORMAT(/1H ,' *** ERROR IS DETECTED AT UNISET : IND=',I4,' ***'/)
      RETURN
      END
      SUBROUTINE UNIVEC (VEC,UEC)
      DIMENSION VEC(3),UEC(3)
      ALENGH=SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2)
      DO 10 J=1,3
      UEC(J)=VEC(J)/ALENGH
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VECPRD (VECA,VECB,VECP)
      DIMENSION VECA(3),VECB(3),VECP(3)
      VECP(1)=VECA(2)*VECB(3)-VECA(3)*VECB(2)
      VECP(2)=VECA(3)*VECB(1)-VECA(1)*VECB(3)
      VECP(3)=VECA(1)*VECB(2)-VECA(2)*VECB(1)
      RETURN
      END
c=======================================================================
c
      subroutine sortf_as (num,x, isort)
c----------------------------------------------------------------------c
c     ascending order                                                  c
c----------------------------------------------------------------------c
      parameter (maxarr=100000)
      dimension  x(*), isort(*)
      dimension  xs(maxarr)
c
      if (num.ge.maxarr) then
        write(6,*) 
     +      'sortf_as.f: sort array too small: redimension/recompile'
        write(6,*) '  size needed is ', num
        stop
      endif
      isort(1) = 1
      if (num.eq.1) return
c
c
      do 10 i=2,num
         isort(i) = i
   10 continue
c
      xs(1) = x(1)
      do 100 i=2,num
         j = i - 1
   55    if (x(i).ge.xs(j)) then
            xs(j+1)    = x(i)
            isort(j+1) = i
            go to 100
         else
            xs(j+1)    = xs(j)
            isort(j+1) = isort(j)
            if (j.eq.1) then
               xs(1)    = x(i)
               isort(1) = i
               go to 100
            else
               j = j - 1
               go to 55
            endif
         endif
  100 continue
c
      return
      end
