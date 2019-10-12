!=======================================================================
! All  subroutines for pre-process
!
! Original author : Aizhong Ye, 05/01/2012 Princeton University
! revise:  Aizhong Ye  
!=======================================================================

! 1.subroutine creleff(cre)
! 2.subroutine creleffa(cre)
! 3.function julday(mm,id,iyyy)
! 4.function gasdev(idum,gset,iset)
! 5.subroutine caldat(caldat1)
! 6.subroutine Wrecord(asRecord)
! 7.subroutine ReadASC(MyASC,asFile)
! 8.subroutine WriteASC(MyASC,asFile)
! 9.subroutine ReadGFS(GFSdata1)
! 10.subroutine ReadObdata(Obsdata1)
! 11.subroutine WriteObdata06(Obsdata1)
! 12.subroutine BinToAdata06(outputfile)
! 13.subroutine Writemapstats(Obsdata1)
! 14.subroutine mapstats06 (stat06a) 
! 15.subroutine mapstats24 (stat24a) 
! 16.subroutine slashchange(anistring) 
! 17.subroutine readDatafcst(rData1)
! 18.1.subroutine readDataob(rData1)
! 18.2.subroutine readDataEp(rData1)
! 19.subroutine WritestatOF(CanE)
! 20.subroutine WritestatOFE(CanE)

! 21.subroutine generateCaE(gCanE) 
! 22.subroutine sort(sort1)
! 23.function GAMMLN(XX)
! 24.subroutine gcf(gammcf1) 
! 25.subroutine gser(gammcf1)
! 26.function gammp(a,x)
! 27.function weibcv(b)
! 28.function weibb(cv)
! 29.function weibp(a,b,x)
! 30.subroutine threshold (threshold1)
! 31.function gausspdf(z)
! 32.function erfcc(x)
! 33.function gaussp(x)
! 34.function gausspi(f)
! 35.subroutine cgauss_old(cgauss1)
! 36.subroutine cgauss(cgauss1)
! 37.function betacf(beta1)
! 38.function betai(beta1)
! 39.function betapi2(beta1)
! 40.function vtrans(vtrans1)
! 41.function binormp2 (bivar1)
! 42.function frho(bivar1)
! 43.function findrho(bivar1)
! 44.subroutine estrho3 (bivar1)
! 45.subroutine epp_fcstparms(fp1) 
! 46.subroutine fill (nparint,x,xts)
! 47.subroutine extractp(extractp1)
! 48.function gammpi(a,f)
! 49.function weibpi(a,b,f)
! 50.function vtransi(vtrans1)
! 51.subroutine fcst2ensem(fe1) 
! 52.function crps(crps1)
! 53.subroutine crps_eval(crpseval1)
! 54.function cumprob(crps1)
! 55.subroutine roc(roc1)
! 56.subroutine ensverify(ensv1)
! 57.subroutine window(win1)
! 58.subroutine EppPara(CanE)
! 59.subroutine NewspaceCE(CanE)
! 60.subroutine deletspaceCE(CanE)
! 61.subroutine WriteCEpara(CanE)
! 62.subroutine ReadCEpara(CanE)
! 63.subroutine fillall(CanE)
! 64.subroutine indexx(sort1)
! 65.subroutine ReadEpara(SrEPara1)
! 66.subroutine NewspaceEP(epp)
! 67.subroutine deletspaceEP(epp)
! 68.subroutine ReadHpara(epp)
! 69.subroutine Hindcast(epp)
! 70.function zero(ix,gset,iset)
! 71.subroutine SchaakeShuffle(epp)
! 72.subroutine Writehincast(epp)
! 73.subroutine ReadCFScoord(CFSData1)
! 74.subroutine ReadCFS(CFSdata1)
! 75.subroutine delblank(mystr)
! 76.subroutine ReadCEparaCG(CanE)
! 77.subroutine ReadCHpara(epp)
! 78.subroutine WriteObdataT(Obsdata1)
! 79.subroutine WritemapstatsT(Obsdata1)
! 80.subroutine CFSinterpolate(CFSData1)
! 81.subroutine readDataobT(rData1)
! 82.subroutine readDataEpT(rData1)
! 83.subroutine ReadCFS4(CFSdata1)
! 84.subroutine Writehincastb(epp)
! 85.subroutine BinToAdata(outputfile)
! 86.subroutine ReadDCFSb(CFSdata1)
! 87.subroutine ReadDCFSb1(CFSdata1)
! 88.subroutine ReadCFSop(CFSdata1)
! 89.subroutine calweight(CFSdata1)
! 90.subroutine ReadCOPpara(epp,CFSData1)
! 91.subroutine Operation(epp)
! 92.subroutine WriteOperation(epp)
! 93.subroutine WriteOperationb(epp)
! 94.subroutine WriteCEparaN(CanE)
! 95.subroutine calMask(CFSdata1)
! 96.subroutine ReadCEparaCC(CanE)
! 97.subroutine readDataobChina(rData1)
! 98.subroutine ReadCEparaCC1(CanE)
! 99.subroutine ReadCEparaCC2(CanE)
! 100.subroutine readmodelbasin(rData1)
! 101.subroutine readstreamflow(rData1)  
! 102.subroutine WriteCEparaF(CanE)
! 103.subroutine ReadCOPparaG(epp)
! 104.subroutine ReadGFSforecast(CFSData1)
! 105.subroutine readDataobChinaEP(rData1)
! 106.subroutine ReadGFSnew(CFSdata1)
! 107.subroutine readDataobZone(rData1)
! 108.subroutine ReadGFSintpara(GFSdata1,GFSdata2,flagC)
! 109.subroutine Read2para(Obsdata1,MyASC,flagC)
! 110.subroutine Read4para(CanE,flagC)
! 111.subroutine Read7para(CFSData1,flagC)
! 112.subroutine Read8para(CFSData2,flagC)
! 113.subroutine Read9para(CanE,flagC)
! 114.subroutine Read12para(CFSData1,flagC)
! 115.

!!********************************!!******************************************************* 
! 114 read  paramater information of system to CFS members correlation
!!********************************!!******************************************************* 
subroutine Read12para(CFSData1,flagC)
use PrecStru
implicit none 
type(CFSData) ::  CFSData1
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=CFSData1.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then    
          
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.ObFile')) read(tvalue,*) CFSData1.ObFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.DFile')) read(tvalue,*) CFSData1.DFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.OutFile')) read(tvalue,*) CFSData1.OutFile  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.leadT')) read(tvalue,*) CFSData1.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.BeginYear')) read(tvalue,*) CFSData1.BeginYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.EndYear')) read(tvalue,*) CFSData1.EndYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.ifile')) read(tvalue,*) CFSData1.ifile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.MName')) read(tvalue,*) CFSData1.MName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.EName')) read(tvalue,*) CFSData1.EName 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.x1')) read(tvalue,*) CFSData1.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.x2')) read(tvalue,*) CFSData1.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.y1')) read(tvalue,*) CFSData1.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.y2')) read(tvalue,*) CFSData1.y2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData5.CFile')) read(tvalue,*) CFSData1.CFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData5.MName')) read(tvalue,*) CFSData1.MName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData5.WFile')) read(tvalue,*) CFSData1.WFile
      endif
    end do 
  10 close(10) 
   
end subroutine
!!********************************!!******************************************************* 
! 113 read  paramater information of system to Calculate the Nash-Sutcliffe efficiency  and correlation coefficient-CFS
!!********************************!!******************************************************* 
subroutine Read9para(CanE,flagC)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then          
         
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint')) read(tvalue,*) CanE.nparint
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x1')) read(tvalue,*) CanE.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x2')) read(tvalue,*) CanE.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y1')) read(tvalue,*) CanE.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y2')) read(tvalue,*) CanE.y2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
          read(tvalue,*) CanE.Events 
          allocate(CanE.Estart(CanE.Events))
          allocate(CanE.Estop(CanE.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)k,(CanE.Estart(i),i=1,CanE.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)k,(CanE.Estop(i),i=1,CanE.Events)  
       
      endif
    end do 
  10 close(10) 
   
end subroutine
!!********************************!!******************************************************* 
! 112 read  paramater information of system to CFS data manage from date form to point form
!!********************************!!******************************************************* 
subroutine Read8para(CFSData2,flagC)
use PrecStru
implicit none 
type(CFSData) ::  CFSData2
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=CFSData2.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then          
                    
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.DFile')) read(tvalue,*) CFSData2.DFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.OutFile')) read(tvalue,*) CFSData2.OutFile 
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.leadT')) read(tvalue,*) CFSData2.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.BeginYear')) read(tvalue,*) CFSData2.BeginYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.EndYear')) read(tvalue,*) CFSData2.EndYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.ifile')) read(tvalue,*) CFSData2.ifile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.MName')) read(tvalue,*) CFSData2.MName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.EName')) read(tvalue,*) CFSData2.EName 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.nparint')) read(tvalue,*) CFSData2.nparint 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.x1')) read(tvalue,*) CFSData2.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.x2')) read(tvalue,*) CFSData2.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.y1')) read(tvalue,*) CFSData2.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.y2')) read(tvalue,*) CFSData2.y2

      endif
    end do 
  10 close(10) 
   
end subroutine
!!********************************!!******************************************************* 
! 111 read  paramater information of system to CFS data merge from Global data
!!********************************!!******************************************************* 
subroutine Read7para(CFSData1,flagC)
use PrecStru
implicit none 
type(CFSData) ::  CFSData1
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=CFSData1.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then          
             
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.CFile')) read(tvalue,*) CFSData1.CFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.DFile')) read(tvalue,*) CFSData1.DFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.OutFile')) read(tvalue,*) CFSData1.OutFile 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.EName')) read(tvalue,*) CFSData1.EName 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.leadT')) read(tvalue,*) CFSData1.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.BeginYear')) read(tvalue,*) CFSData1.BeginYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.EndYear')) read(tvalue,*) CFSData1.EndYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.ifile')) read(tvalue,*) CFSData1.ifile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.MName')) read(tvalue,*) CFSData1.MName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.ihour'))   read(tvalue,*) CFSData1.ihour
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.imonth'))   read(tvalue,*) CFSData1.imonth
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.iday'))   read(tvalue,*) CFSData1.iday
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.WFile')) read(tvalue,*) CFSData1.WFile
      
      endif
    end do 
  10 close(10) 
   
end subroutine
!!********************************!!******************************************************* 
! 110 read  paramater information of system to Calculate the Nash-Sutcliffe efficiency  and correlation coefficient
!!********************************!!******************************************************* 
subroutine Read4para(CanE,flagC)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE 
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then          
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
          read(tvalue,*) CanE.Events 
          allocate(CanE.Estart(CanE.Events))
          allocate(CanE.Estop(CanE.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)k,(CanE.Estart(i),i=1,CanE.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)k,(CanE.Estop(i),i=1,CanE.Events)  
                     
      
      endif
    end do 
  10 close(10) 
   
end subroutine
!!********************************!!******************************************************* 
! 109 read  paramater information of system to daily observed data to 6 hours data
!!********************************!!******************************************************* 
subroutine Read2para(Obsdata1,MyASC,flagC)
use PrecStru
implicit none 
type(Obsdata) ::  Obsdata1
type(ASCFile) ::  MyASC
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=Obsdata1.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then          
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('BeginYear')) read(tvalue,*) Obsdata1.BeginYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('EndYear')) read(tvalue,*) Obsdata1.EndYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Ngrids')) read(tvalue,*) Obsdata1.SumN
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('originalFile')) read(tvalue,*) Obsdata1.BName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('6hourFile')) read(tvalue,*) Obsdata1.BName06
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('obw')) read(tvalue,*) Obsdata1.obw
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.Tim')) read(tvalue,*) Obsdata1.Tim
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.CuVa')) read(tvalue,*) Obsdata1.CuVa
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.ncols')) read(tvalue,*) Obsdata1.ncols
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.nrows')) read(tvalue,*) Obsdata1.nrows
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.NVara')) read(tvalue,*) Obsdata1.NVara
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.undef')) read(tvalue,*) Obsdata1.undef

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('xllcorner')) read(tvalue,*) MyASC.xllcorner
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('yllcorner')) read(tvalue,*) MyASC.yllcorner
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('cellsize')) read(tvalue,*) MyASC.cellsize
            
      endif
    end do 
  10 close(10) 
   
end subroutine

!!********************************!!******************************************************* 
! 108 read  paramater information of system to 2.5*2.5 GFS data interpolate
!!********************************!!******************************************************* 
subroutine ReadGFSintpara(GFSdata1,GFSdata2,flagC)
use PrecStru
implicit none 
type(GFSdata) ::  GFSdata1
type(GFSdata) ::  GFSdata2
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord) 
  open(unit=10,file=GFSdata1.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then          
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('BeginYear')) read(tvalue,*) GFSdata1.BeginYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('EndYear')) read(tvalue,*) GFSdata1.EndYear

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.ncols')) read(tvalue,*) GFSdata1.ncols
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.nrows')) read(tvalue,*) GFSdata1.nrows
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.leadT')) read(tvalue,*) GFSdata1.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.cellsize')) read(tvalue,*) GFSdata1.cellsize
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.xllcorner')) read(tvalue,*) GFSdata1.xllcorner
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.yllcorner')) read(tvalue,*) GFSdata1.yllcorner
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.ncols')) read(tvalue,*) GFSdata2.ncols
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.nrows')) read(tvalue,*) GFSdata2.nrows
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.leadT')) read(tvalue,*) GFSdata2.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.cellsize')) read(tvalue,*) GFSdata2.cellsize
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.xllcorner')) read(tvalue,*) GFSdata2.xllcorner
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.yllcorner')) read(tvalue,*) GFSdata2.yllcorner

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Weightfile')) read(tvalue,*) GFSdata1.Weightfile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('GFSdata1.BName')) read(tvalue,*) GFSdata1.BName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('GFSdata2.BName')) read(tvalue,*) GFSdata2.BName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('GFSdata1.MName')) read(tvalue,*) GFSdata1.MName  
          
            
      endif
    end do 
  10 close(10) 
   
end subroutine


!!********************************!!******************************************************* 
!107 read observed  data  in a zone
!!********************************!!******************************************************* 
subroutine readDataobZone(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
 
integer i,j,year1,year2,ndays 
integer imonth1,iday1 ,imonth2,iday2  
integer julday,jday1,jday2
character*40  tempc  ,ay,am,ad
real xobs(4) ,revise
 
 
  jday1=julday(1,1,rData1.BeginYear)
  rData1.ndays=julday(12,31,rData1.EndYear)-jday1+1
  open (10,file=rData1.ObFile)
    read (10,*)tempc
    read (10,'(a8,a4,a1,a2,a1,a2)') tempc,ay,tempc,am,tempc,ad 
    read(ay,*)year1
    read(am,*)imonth1
    read(ad,*)iday1
    read (10,'(a6,a4,a1,a2,a1,a2)') tempc,ay,tempc,am,tempc,ad 
    read(ay,*)year2
    read(am,*)imonth2
    read(ad,*)iday2
    
   ! read (10,'(a,a4,a1,a2,a1,a2)') tempc,year1,tempc,imonth1,tempc,iday1 
    !read (10,*) tempc,year2,imonth2,iday2  
    read (10,*)
    read (10,*) tempc,revise
    read (10,*)
    jday1=julday(imonth1,iday1,year1)
    jday2=julday(rData1.Bmonth,rData1.Bday,rData1.BeginYear)
    rData1.ob=0
   
 
    if (rData1.BeginYear.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYear,'<',year1
	  return
    end if
    if (rData1.EndYear.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYear,'>',year2
	  return
    end if
    j=jday2-jday1+1
    do i=1,j
      read (10,*)        
    end do 
    !if (rData1.nparint.eq.1)  then
    !  ndays=rData1.ndays      
    !else
    ndays=rData1.ndays+rData1.leadT
    !end if
    do i=1,ndays
	  !write(*,*) i
      read (10,*) (rData1.ob(i,1),j=1,rData1.y)
      if (rData1.flagT .eq.0) then
        rData1.ob(i,1)=max(0.0,rData1.ob(i,1)/revise)
      else
        rData1.ob(i,1)=rData1.ob(i,1)/revise+273.15  
      end if
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ob(i,j)=rData1.ob(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ob(i,j)=rData1.ob(i,j-1)
      end do
    end do
  close(10)

end subroutine


!!********************************!!******************************************************* 
!106 Read GFS File  
!!********************************!!******************************************************* 
subroutine ReadGFSnew(CFSdata1)
use PrecStru
implicit none
integer i,j,k ,ii,jj
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C1000)  asFile,tempc 
integer JULDAY,jday1
logical alive
   

  jday1=JULDAY(1,1,CFSdata1.iYear) 
  caldat1.julian=jday1+CFSdata1.cday-1
  call caldat(caldat1)  
  CFSdata1.iYear=caldat1.iyyy
  CFSdata1.imonth=caldat1.mm
  CFSdata1.iday=caldat1.id
  alive=.false.
  k=CFSdata1.ihour
   
  write(asFile,'(2a,i4,2i2.2,a)') trim(CFSdata1.DFile),'\apcp_sfc_',&
      CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	   '00_mean.grib2' 
  if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
  inquire(file=asFile,exist=alive)
  CFSdata1.alive=alive 
  if (.not.alive) then     
    print *,trim(asFile),' is not exist.'
    !stop
    return
  end if
  !apcp_sfc_2003110100_mean.grib2 trim(CFSdata1.EName),CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',
  if (CFSData1.flagOs.eq.2) then   
    write(asFile,'(3a,i4,2i2.2,3a,i2)')'./wgrib2  ',trim(CFSdata1.DFile),'\apcp_sfc_',&
      CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	   '00_mean.grib2 -bin ',trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile   
  else
    read(CFSdata1.DFile,'(a2,a)') tempc,asFile
    write(asFile,'(a,a1,2a,i4,2i2.2,3a,i2)')'wgrib2  \cygdrive\',tempc,trim(asFile),'\apcp_sfc_',&
      CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
      '00_mean.grib2 -bin ',trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile   
  end if	
  call slashchange(asFile)
  call system(asFile)
  if (CFSData1.flagOs.eq.2) then
    write(asFile,'(2a,i2)')trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile        
    call slashchange(asFile)
  else
    write(asFile,'(2a,i2)') trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile     
  end if
  open(unit=20,file=asFile,form='unformatted',access='sequential')
  j=CFSdata1.leadT 
  CFSdata1.MValue=>CFSdata1.MValue00
  do i=1,45     
   ! write(*,*)"i=",i
    read(20,end=20) CFSdata1.MValue(i,:,:) 
     !write(*,*)CFSdata1.MValue(i,1,1) 
  end do
20  close(20)
  if (CFSdata1.leadT >8) then !open 8-16 days file  
    write(asFile,'(2a,i4,2i2.2,a)') trim(CFSdata1.DFile),'\apcp_sfc_',&
      CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	   '00_mean_t190.grib2' 
    if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
    inquire(file=asFile,exist=alive)
    CFSdata1.alive=alive 
    if (.not.alive) then     
      print *,trim(asFile),' is not exist.'
      !stop
      return
    end if
    !apcp_sfc_2003110100_mean.grib2 trim(CFSdata1.EName),CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',
    if (CFSData1.flagOs.eq.2) then   
      write(asFile,'(3a,i4,2i2.2,3a,i2)')'./wgrib2  ',trim(CFSdata1.DFile),'\apcp_sfc_',&
        CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	     '00_mean_t190.grib2 -bin ',trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile+50   
    else
      read(CFSdata1.DFile,'(a2,a)') tempc,asFile
      write(asFile,'(a,a1,2a,i4,2i2.2,3a,i2)')'wgrib2  \cygdrive\',tempc,trim(asFile),'\apcp_sfc_',&
        CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
        '00_mean_t190.grib2 -bin ',trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile+50  
    end if	 
    call slashchange(asFile)
    call system(asFile)  
    if (CFSData1.flagOs.eq.2) then
      write(asFile,'(2a,i2)')trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile+50      
      call slashchange(asFile)
    else
      write(asFile,'(2a,i2)') trim(CFSdata1.OutFile),'\fort.',CFSData1.ifile+50    
    end if
    open(unit=20,file=asFile,form='unformatted',access='sequential')
    j=CFSdata1.leadT 
    CFSdata1.MValue=>CFSdata1.MValue00
    do i=44,75     
     ! write(*,*)"i=",i
      read(20,end=30) CFSdata1.MValue(i,:,:) 
      !write(*,*)CFSdata1.MValue(i,1,1) 
    end do
30    close(20)
  endif  

 
   
  ! 6 hours data to daily data
  if (trim(CFSdata1.EName).eq.'prate')  then
    do i=1,3
      CFSdata1.MValue(i,:,:) =0;
      do k=(i-1)*8+3,(i-1)*8+9,2    
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(k,:,:) 
      end do
    end do
    do i=4,CFSdata1.leadT  
      !write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =0
      do k=1,4
        j=(2+i)*4+k+1
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(j,:,:) 
      end do       
    end do 
  end if

  if (trim(CFSdata1.EName).eq.'tmp2m') then
  !GEFS format 
    do i=1,3
      CFSdata1.MValue(i,:,:) =0;
      do k=(i-1)*8+3,(i-1)*8+9,2    
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(k,:,:) 
      end do
    end do
    do i=4,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =0
      do k=1,4
        j=(2+i)*4+k+1
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(j,:,:) 
      end do       
    end do     
    do i=1,CFSdata1.leadT  
      CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:)/4
    end do
  end if

  if (trim(CFSdata1.EName).eq.'tmax2m') then
!GEFS format   
    do i=1,3
      CFSdata1.MValue(i,:,:) =-9999
      do k=(i-1)*8+3,(i-1)*8+9,2    
        do ii=1,CFSData1.ncols
          do jj=1,CFSData1.nrows  
            if (CFSdata1.MValue(i,ii,jj) .lt. CFSdata1.MValue(k,ii,jj))  CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(k,ii,jj)
          end do
        end do 
      end do
    end do
    do i=4,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =-9999
      do k=1,4
        j=(2+i)*4+k+1
        do ii=1,CFSData1.ncols
          do jj=1,CFSData1.nrows  
            if (CFSdata1.MValue(i,ii,jj) .lt. CFSdata1.MValue(j,ii,jj))  CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(j,ii,jj)
          end do
        end do       
         
      end do       
    end do      
  end if  
 if (trim(CFSdata1.EName).eq.'tmin2m') then
!GEFS format   
    do i=1,3
      CFSdata1.MValue(i,:,:) =9999
      do k=(i-1)*8+3,(i-1)*8+9,2    
        do ii=1,CFSData1.ncols
          do jj=1,CFSData1.nrows  
            if (CFSdata1.MValue(i,ii,jj) .gt. CFSdata1.MValue(k,ii,jj))  CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(k,ii,jj)
          end do
        end do 
      end do
    end do
    do i=4,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =9999
      do k=1,4
        j=(2+i)*4+k+1
        do ii=1,CFSData1.ncols
          do jj=1,CFSData1.nrows  
            if (CFSdata1.MValue(i,ii,jj) .gt. CFSdata1.MValue(j,ii,jj))  CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(j,ii,jj)
          end do
        end do       
         
      end do       
    end do      
  end if    
  
end subroutine ReadGFSnew





!!********************************!!******************************************************* 
!105 read observed  data  in China for schaake shuffle
!!********************************!!******************************************************* 
subroutine readDataobChinaEP(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays 
integer imonth1,iday1 ,imonth2,iday2  
integer julday,jday1,jday2
character*40  tempc  
real xobs(4) 
  if (rData1.flagOs .eq. 2)   then
    write(inputfile,'(a,i3.3,a)')trim(rData1.EpFile),rData1.y,'.txt'  
  else
    write(inputfile,'(a,i3.3,a)')trim(rData1.EpFile),rData1.y,'.txt'  
  end if
  jday1=julday(1,1,rData1.BeginYearEns)
  rData1.ndays=julday(12,31,rData1.EndYearEns)-jday1+1
  open (10,file=inputfile)
    read (10,*) tempc,year1,imonth1,iday1 
    read (10,*) tempc,year2,imonth2,iday2  
    read (10,*)
    jday1=julday(imonth1,iday1,year1)
    jday2=julday(1,1,rData1.BeginYearEns)
    rData1.ep=0
   
 
    if (rData1.BeginYearEns.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYear,'<',year1
	  return
    end if
    if (rData1.EndYearEns.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYear,'>',year2
	  return
    end if
    j=jday2-jday1+1
    do i=1,j
      read (10,*)        
    end do 
    !if (rData1.nparint.eq.1)  then
    !  ndays=rData1.ndays      
    !else
    ndays=rData1.ndays+rData1.leadT
    !end if
    do i=1,ndays
	  !write(*,*) i
      read (10,*) rData1.ep(i,1)
      rData1.ep(i,1)=max(0.0,rData1.ep(i,1)/100)
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ep(i,j)=rData1.ep(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ep(i,j)=rData1.ep(i,j-1)
      end do
    end do
  close(10)

end subroutine

!!********************************!!******************************************************* 
!104 Read GFS Forecast
!!********************************!!******************************************************* 
subroutine ReadGFSforecast(CFSData1)
use PrecStru
implicit none
type(CFSData) ::  CFSData1
 
character(c100) asRecord
character(20) tempa 
integer i,j,k
logical alive
  
  inquire(file=CFSData1.DFile,exist=alive)   
  if(alive) then
    open(10,File=CFSData1.DFile,status='old')    
      do i=1,CFSData1.x2
        do j=1,CFSData1.y2
          read(10,*) (CFSData1.MValueN(k,i,j),k=1,CFSData1.leadT)	 
          do k=1,CFSData1.leadT
		    CFSData1.MValueN(k,i,j)=CFSData1.MValueN(k,i,j)/100
          end do
        end do
      end do      
    close(10)    
  else
    Write(*,*) "No ",CFSData1.DFile
  return
  endif
  write(asRecord,"('Finish reading file of ',a)")trim(CFSData1.DFile)
  call Wrecord(asRecord)
end subroutine ReadGFSforecast

!!********************************!!******************************************************* 
! 103 read  paramater information of system to operational forecast for GFS
!!********************************!!******************************************************* 
subroutine ReadCOPparaG(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
 
character(c100) ::   tempc,tvalue
integer  i
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
 
  open(unit=10,file=epp.CoFile)
    do 
      read(10,*,end=10)tempc,tvalue     

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.iYear'))   read(tvalue,*) epp.iYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.imonth'))  read(tvalue,*) epp.imonth
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.iday'))   read(tvalue,*) epp.iday

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.EName')) read(tvalue,*) epp.EName 

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.EpFile')) read(tvalue,*) epp.EpFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.ObFile')) read(tvalue,*) epp.ObFile
 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.SiFile')) read(tvalue,*) epp.SiFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.SiFile1')) read(tvalue,*) epp.SiFile1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.PaFile')) read(tvalue,*) epp.PaFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.flagCG')) read(tvalue,*) epp.flagCG
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.flagT')) read(tvalue,*) epp.flagT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.leadT')) read(tvalue,*) epp.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.leadT1')) read(tvalue,*) epp.leadT1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.nparint')) read(tvalue,*) epp.nparint
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.x1')) read(tvalue,*) epp.x1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.x2')) read(tvalue,*) epp.x2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.y1')) read(tvalue,*) epp.y1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.y2')) read(tvalue,*) epp.y2
 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.DoFile')) read(tvalue,*) epp.DoFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.ReFile')) read(tvalue,*) epp.ReFile

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.BeginYear ')) read(tvalue,*) epp.BeginYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.EndYear')) read(tvalue,*) epp.EndYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.BeginYearEns ')) read(tvalue,*) epp.BeginYearEns 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.EndYearEns')) read(tvalue,*) epp.EndYearEns 

 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.iobstran ')) read(tvalue,*) epp.iobstran 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.ifcsttran ')) read(tvalue,*) epp.ifcsttran 
      
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.Events')) then
        read(tvalue,*) epp.Events 
        allocate(epp.Estart(epp.Events))
        allocate(epp.Estop(epp.Events))
      end if 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.Estart'))read(10,*)(epp.Estart(i),i=1,epp.Events)  
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp3.Estop'))read(10,*)(epp.Estop(i),i=1,epp.Events)  
  
    end do 
  10 close(10) 
  epp.nems=epp.EndYearEns-epp.BeginYearEns+1
end subroutine



!!********************************!!******************************************************* 
! 102 Write stream flow paramater information of EPP
!!********************************!!******************************************************* 
subroutine WriteCEparaF(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   outputfile
character(c1000) ::tempc
integer  i,j,k  
character(c100) asRecord 
 
  write(asRecord,*)'Begin to write parameter'
  call Wrecord(asRecord)
  if  (CanE.flagSE==1) then
    write(outputfile,'(4a)')trim(CanE.ReFile),trim(CanE.model(CanE.x)), &
      trim(CanE.basin(CanE.y)), '_re.txt'   
  else
    write(outputfile,'(4a)')trim(CanE.ReFile),trim(CanE.model(CanE.x)), &
      trim(CanE.basin(CanE.y)), '_re1.txt'   
  end if
 
  if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) 
 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('rho_est'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rho_est(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

 	  write (10,'(a,TR1,i3,TR1,i3)') trim('pthresh_obs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.pthresh_obs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('pthresh_fcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.pthresh_fcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
 
 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('popobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.popobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cavgobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cavgobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccvobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccvobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('popfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.popfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cavgfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cavgfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccvfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccvfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

 

  close(10)

end subroutine 
!!********************************!!******************************************************* 
!101 read streamflow  data   
!!********************************!!******************************************************* 
subroutine readstreamflow(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays 
integer julday,jday1,jday2
character*40  tempc  

 
  write(inputfile,'(4a)')trim(rData1.ObFile),trim(rData1.model(rData1.x)), &
      trim(rData1.basin(rData1.y)), '_result1.txt'   
  open (10,file=inputfile) 
    read (10,*)   
    ndays=rData1.ndays
    do i=1,ndays
	  do j=1,rData1.leadT
        read (10,*)tempc,tempc,rData1.ob(i,j),rData1.si(i,j),rData1.ep(i,j)
      end do 
    end do
    
  close(10)

end subroutine

!!********************************!!******************************************************* 
!100 read models name and basins number
!!********************************!!******************************************************* 
subroutine readmodelbasin(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
 
integer i,j 
  
  open (10,file=rData1.SiFile)
 
    do i=1,rData1.x
      read (10,*) rData1.model(i)       
    end do 
    do i=1,rData1.y
      read (10,*) rData1.basin(i)       
    end do 
    
  close(10)

end subroutine


!!********************************!!******************************************************* 
! 99 read  paramater information of system to prepare parameters of EPP for stream flow
!!********************************!!******************************************************* 
subroutine ReadCEparaCC2(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   tempc,tvalue
integer  i,k
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
   
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue  
      if (k.eq.15) then
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo
        !if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.PTmin')) read(tvalue,*) CanE.PTmin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagCG')) read(tvalue,*) CanE.flagCG
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagT')) read(tvalue,*) CanE.flagT

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile1')) read(tvalue,*) CanE.SiFile1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT1')) read(tvalue,*) CanE.leadT1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint')) read(tvalue,*) CanE.nparint
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x1')) read(tvalue,*) CanE.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x2')) read(tvalue,*) CanE.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y1')) read(tvalue,*) CanE.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y2')) read(tvalue,*) CanE.y2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minall')) read(tvalue,*) CanE.minall 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minpos ')) read(tvalue,*) CanE.minpos 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.pop_fraction ')) read(tvalue,*) CanE.pop_fraction 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iobstran ')) read(tvalue,*) CanE.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ifcsttran ')) read(tvalue,*) CanE.ifcsttran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iopt_rho ')) read(tvalue,*) CanE.iopt_rho 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cor_weight ')) read(tvalue,*) CanE.cor_weight 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgobsx ')) read(tvalue,*) CanE.cavgobsx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgfcstx ')) read(tvalue,*) CanE.cavgfcstx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.verstats ')) read(tvalue,*) CanE.verstats 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint ')) read(tvalue,*) CanE.nparint 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nems ')) read(tvalue,*) CanE.nems 
    
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
          read(tvalue,*) CanE.Events 
          allocate(CanE.Estart(CanE.Events))
          allocate(CanE.Estop(CanE.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)k,(CanE.Estart(i),i=1,CanE.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)k,(CanE.Estop(i),i=1,CanE.Events)  
      endif
    end do 
  10 close(10) 
end subroutine



!!********************************!!******************************************************* 
! 98 read  paramater information of system to prepare parameters of EPP for CFS China
!!********************************!!******************************************************* 
subroutine ReadCEparaCC1(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   tempc,tvalue
integer  i
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
   
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)tempc,tvalue        
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo
      !if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.PTmin')) read(tvalue,*) CanE.PTmin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagCG')) read(tvalue,*) CanE.flagCG
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagT')) read(tvalue,*) CanE.flagT

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile1')) read(tvalue,*) CanE.SiFile1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT1')) read(tvalue,*) CanE.leadT1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint')) read(tvalue,*) CanE.nparint
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x1')) read(tvalue,*) CanE.x1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x2')) read(tvalue,*) CanE.x2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y1')) read(tvalue,*) CanE.y1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y2')) read(tvalue,*) CanE.y2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minall')) read(tvalue,*) CanE.minall 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minpos ')) read(tvalue,*) CanE.minpos 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.pop_fraction ')) read(tvalue,*) CanE.pop_fraction 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iobstran ')) read(tvalue,*) CanE.iobstran 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ifcsttran ')) read(tvalue,*) CanE.ifcsttran 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iopt_rho ')) read(tvalue,*) CanE.iopt_rho 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cor_weight ')) read(tvalue,*) CanE.cor_weight 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgobsx ')) read(tvalue,*) CanE.cavgobsx 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgfcstx ')) read(tvalue,*) CanE.cavgfcstx 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.verstats ')) read(tvalue,*) CanE.verstats 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint ')) read(tvalue,*) CanE.nparint 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nems ')) read(tvalue,*) CanE.nems 
    
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
        read(tvalue,*) CanE.Events 
        allocate(CanE.Estart(CanE.Events))
        allocate(CanE.Estop(CanE.Events))
      end if 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)(CanE.Estart(i),i=1,CanE.Events)  
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)(CanE.Estop(i),i=1,CanE.Events)  
  
    end do 
  10 close(10) 
end subroutine


!!********************************!!******************************************************* 
!97 read observed  data  in China
!!********************************!!******************************************************* 
subroutine readDataobChina(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays 
integer imonth1,iday1 ,imonth2,iday2  
integer julday,jday1,jday2
character*40  tempc  
real xobs(4) 
  if (rData1.flagOs .eq. 2)   then
    write(inputfile,'(a,i3.3,a)')trim(rData1.ObFile),rData1.y,'.txt'  
  else
    write(inputfile,'(a,i3.3,a)')trim(rData1.ObFile),rData1.y,'.txt'  
  end if
  jday1=julday(1,1,rData1.BeginYear)
  rData1.ndays=julday(12,31,rData1.EndYear)-jday1+1
  open (10,file=inputfile)
    read (10,*) tempc,year1,imonth1,iday1 
    read (10,*) tempc,year2,imonth2,iday2  
    read (10,*)
    jday1=julday(imonth1,iday1,year1)
    jday2=julday(rData1.Bmonth,rData1.Bday,rData1.BeginYear)
    rData1.ob=0
   
 
    if (rData1.BeginYear.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYear,'<',year1
	  return
    end if
    if (rData1.EndYear.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYear,'>',year2
	  return
    end if
    j=jday2-jday1+1
    do i=1,j
      read (10,*)        
    end do 
    !if (rData1.nparint.eq.1)  then
    !  ndays=rData1.ndays      
    !else
    ndays=rData1.ndays+rData1.leadT
    !end if
    do i=1,ndays
	  !write(*,*) i
      read (10,*) rData1.ob(i,1)
      rData1.ob(i,1)=max(0.0,rData1.ob(i,1)/100)
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ob(i,j)=rData1.ob(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ob(i,j)=rData1.ob(i,j-1)
      end do
    end do
  close(10)

end subroutine



!!********************************!!******************************************************* 
! 96 read  paramater information of system to prepare parameters of EPP for CFS China
!!********************************!!******************************************************* 
subroutine ReadCEparaCC(CanE,flagC)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
   
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue        
      if (k.eq.flagC) then
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo
        !if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.PTmin')) read(tvalue,*) CanE.PTmin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagCG')) read(tvalue,*) CanE.flagCG
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagT')) read(tvalue,*) CanE.flagT

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile1')) read(tvalue,*) CanE.SiFile1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT1')) read(tvalue,*) CanE.leadT1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint')) read(tvalue,*) CanE.nparint
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x1')) read(tvalue,*) CanE.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x2')) read(tvalue,*) CanE.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y1')) read(tvalue,*) CanE.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y2')) read(tvalue,*) CanE.y2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minall')) read(tvalue,*) CanE.minall 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minpos ')) read(tvalue,*) CanE.minpos 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.pop_fraction ')) read(tvalue,*) CanE.pop_fraction 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iobstran ')) read(tvalue,*) CanE.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ifcsttran ')) read(tvalue,*) CanE.ifcsttran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iopt_rho ')) read(tvalue,*) CanE.iopt_rho 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cor_weight ')) read(tvalue,*) CanE.cor_weight 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgobsx ')) read(tvalue,*) CanE.cavgobsx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgfcstx ')) read(tvalue,*) CanE.cavgfcstx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.verstats ')) read(tvalue,*) CanE.verstats 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint ')) read(tvalue,*) CanE.nparint 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nems ')) read(tvalue,*) CanE.nems 
    
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
          read(tvalue,*) CanE.Events 
          allocate(CanE.Estart(CanE.Events))
          allocate(CanE.Estop(CanE.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)k,(CanE.Estart(i),i=1,CanE.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)k,(CanE.Estop(i),i=1,CanE.Events)  
      endif
    end do 
  10 close(10) 
end subroutine


!!********************************!!******************************************************* 
!95 Calculate the MASK of sub-domain in CFS domain 
!!********************************!!******************************************************* 
subroutine calMask(CFSdata1)
use PrecStru
implicit none
integer i,j,k,flagOs
type(CFSdata) ::  CFSdata1 
type(ASCFile) ::  obw
integer,pointer :: mask(:,:)
real sump
character(C1000)  asFile,tempc 
real(r8)  x1,y1,xx,yy,wx1,wy1,wx2,wy2 

    flagOs=CFSData1.flagOs
    if (flagOs .eq. 2)   call slashchange(CFSData1.MName)
    call  ReadASC(obw,CFSData1.MName)
    CFSData1.obw=obw

    allocate(CFSData1.mask(CFSData1.ncols,CFSData1.nrows))
	mask=>CFSData1.mask
    mask=0
    allocate(CFSData1.position(obw.ncols,obw.nrows))
    do i=1,obw.ncols
      j=1
      xx=obw.xllcorner + (i-1)*obw.cellsize
      if (xx.lt.0) xx=xx+360
      do while (CFSData1.x(j).le.xx .and. j.le.CFSData1.ncols)
      j=j+1
      end do
      CFSData1.position(i,:).x=j
    end do
    do i=1,obw.nrows
      j=1
      yy=obw.yllcorner + (i-1)*obw.cellsize 
      do while (CFSData1.y(j).le.yy .and. j.le.CFSData1.nrows)
      j=j+1
      end do
      CFSData1.position(:,i).y=j               
    end do

    ! USA 0.125: X1=X0+i*0.125, Y1=Y0+j*0.125    
    ! GFS 2.5:   X2=X0+i*2.5,   Y2=Y0+j*2.5
    do i=1,obw.ncols
      do j=1,obw.nrows    
        if (abs(obw.MValue(i,j)-obw.NODATA_value)>0.01) &
          mask(CFSData1.position(i,j).x,CFSData1.position(i,j).y)=obw.MValue(i,j)  
      end do
    end do
    if (flagOs .eq. 2)   call slashchange(CFSData1.WFile)
    open(10,file=CFSData1.WFile)  
      do j=CFSData1.nrows,1,-1
        write(10,'(999(i5,a))')(mask(i,j),char(9),i=1,CFSData1.ncols)
	   
      end do
    close(10)

    
  
end subroutine 


!!********************************!!******************************************************* 
! 94 Write necessary paramater information of EPP
!!********************************!!******************************************************* 
subroutine WriteCEparaN(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   outputfile
character(c1000) ::tempc
integer  i,j,k  
character(c100) asRecord 
 
  write(asRecord,*)'Begin to write parameter'
  call Wrecord(asRecord)
   
 
  write(outputfile,'(a,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),CanE.x,'\',CanE.y,'.txt' 
  if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) 
 

 	  write (10,'(a,TR1,i3,TR1,i3)') trim('pthresh_obs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.pthresh_obs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('pthresh_fcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.pthresh_fcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
 
 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('popobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.popobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cavgobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cavgobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccvobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccvobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('popfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.popfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cavgfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cavgfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccvfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccvfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rho_est'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rho_est(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
 

  close(10)

end subroutine 

!!********************************!!******************************************************* 
! 93 Write  operational result -binary
!!********************************!!******************************************************* 
subroutine WriteOperationb(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(Scaldat) ::  caldat1
character(c100) ::   outputfile
integer  i,j,k,jday1,imem,JULDAY ,iyear 
character(c100) asRecord 
character(10) tempc  
character(8) :: Tdate   

  write(asRecord,*)'Begin to write Operational'
  !jday1=JULDAY(1,1,epp.BeginYear)  
  write(Tdate,'(i4,2i2.2)') epp.iYear,epp.imonth,epp.iday
  write(outputfile,'(5a)')trim(epp.ReFile),trim(epp.EName),'\',Tdate,'r.bin' 
  call Wrecord(outputfile)
  if (epp.flagOs .eq. 2)   call slashchange(outputfile)
  !real(r8),pointer ::   PensArea(:,:,:,:) !Ensemble [Lon,Lat,lead time,members]
  
  open (10,file=outputfile,status='REPLACE',form='unformatted',access='direct',recl=epp.ncols*epp.nrows)
  k=0  
    do i=1,epp.leadT        
        do j=1,epp.nems+1     
		  k=k+1    
          write (10,rec=k)epp.PensArea(:,:,i,j)
          !(int(epp.Pens(i,j,imem)*100),imem=1,epp.nems)        
        end do !j      
    end do    
  close(10)
  
  
end subroutine 



!!********************************!!******************************************************* 
! 92 Write  Operational result
!!********************************!!******************************************************* 
subroutine WriteOperation(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp 
character(c100) ::   outputfile
character(8) :: Tdate
integer  i,j,imem 
character(c100) asRecord 
character(c1000):: tempc 
logical alive
real tk
  
  write(asRecord,*)'Begin to write Operational'
  !jday1=JULDAY(1,1,epp.BeginYear)  
  write(Tdate,'(i4,2i2.2)') epp.iYear,epp.imonth,epp.iday
  write(outputfile,'(4a,2i3.3,a)')trim(epp.ReFile),trim(epp.EName),'\',Tdate,epp.x,epp.y,'.txt' 
  call Wrecord(outputfile)
  if (epp.flagOs .eq. 2)   call slashchange(outputfile)

  inquire(file=outputfile,exist=alive)
  
  !if(alive) then 
   ! open(10,File=outputfile,status='old',position='append')
  !else  
    open(10,File=outputfile)
    write (10,'(11a,\)')'Lon',char(9),'Lat',char(9),&
	   'LeadT',char(9),'Obs',char(9),'Fcst',char(9),'EnsM'
    do i=1,epp.nems
	  write (10,'(2a,i3.3,\)')char(9),'Ens',i
    end do
    write (10,*) 
  !endif  

  i=1
  tempc=''
  if (epp.flagT.eq.0)  then
    tk=0
  else
    tk=273.16
  end if
  if (epp.Pens(i,1,1).ge.0.) then	   
    do j=1,epp.leadT          
      !write(tempc,'(i4,100(a,i7))')epp.x,char(9),epp.y,char(9),j,char(9), &
      !int((epp.si(i,j)-tk)*100),char(9),int((epp.Pens(i,j,epp.nems+1)-tk)*100),& 
      !  ((char(9),int((epp.Pens(i,j,imem)-tk)*100)),imem=1,epp.nems) 
      write(tempc,'(i4,2(a,i4),100(a,f7.2))')epp.x,char(9),epp.y,char(9),j,char(9), &
      epp.ep(i,j)-tk,char(9),epp.si(i,j)-tk,char(9),epp.Pens(i,j,epp.nems+1)-tk,& 
        (char(9),epp.Pens(i,j,imem)-tk,imem=1,epp.nems) 
      call delblank(tempc)
      write (10,'(a)') trim(tempc)
	  !write(*,*)epp.Pens(i,j,1)
    end do !j
  endif
       
  close(10)

end subroutine 


!!********************************!!******************************************************* 
! 91  Operational CFS
!!********************************!!******************************************************* 
subroutine Operation(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(SgenerateCaE) ::  gCanE
type(Sfcst2ensem) :: fe1 
type(Scaldat) ::  caldat1 
integer i,j,iyear
integer JULDAY
character(c100) asRecord 

  write(asRecord,*)'Begin to Hindcast CFS'
  call Wrecord(asRecord)

  epp.nems=epp.EndYearEns-epp.BeginYearEns+1
  
  fe1.ifcst=epp.ifcsttran
  fe1.iobs=epp.iobstran
  fe1.npp=epp.nems
  epp.Eens=-99.0
  !1.Generate Canonical Events 
  gCanE.leadT=epp.leadT
  gCanE.ndays=epp.ndays
  gCanE.Events=epp.Events
  gCanE.Estart=>epp.Estart
  gCanE.Estop=>epp.Estop
  !gCanE.ob=>epp.ob
  !gCanE.Eob=>epp.Eob  ![date, Events]
  !call generateCaE(gCanE)
  gCanE.ob=>epp.si
  gCanE.Eob=>epp.Esi    
  call generateCaE(gCanE)

  gCanE.ndays=epp.ndays1
  gCanE.ob=>epp.ep
  gCanE.Eob=>epp.Eep    
  call generateCaE(gCanE)

  iyear=epp.iYear 
  j=JULDAY(epp.imonth,epp.iday,iyear)-JULDAY(1,1,iyear)
  j=min(365,j)
        
  !write(*,'(a,i5,a1,i5)')trim('year,day='),iyear,'/',j

  do i=1,epp.Events
    !2.Generate Canonical Events ensembles
    fe1.pthresh_fcst=epp.pthresh_fcst(j,i)  ![365,Events] 
    fe1.fpop=epp.popfcst(j,i)
    fe1.fcavg= epp.cavgfcst(j,i)
    fe1.fccv= epp.ccvfcst(j,i)
    fe1.pthresh_obs=epp.pthresh_obs(j,i)
    fe1.obspop=epp.popobs(j,i)
    fe1.obscavg=epp.cavgobs(j,i)
    fe1.obsccv=epp.ccvobs(j,i)
    fe1.rho=epp.rho_est(j,i)
    fe1.fcst=epp.Esi(1,i) ![date, Events]
    fe1.pp=>epp.Eens(1,i,:)    ![date,Events,members]
    call fcst2ensem(fe1)
  end do !i Events    
 
  write(asRecord,*)'Begin to SchaakeShuffle CFS'
  call Wrecord(asRecord)

  call SchaakeShuffle(epp)


  write(asRecord,*)'Finished to SchaakeShuffle CFS'
  call Wrecord(asRecord)
end subroutine

!!********************************!!******************************************************* 
! 90 read  paramater information of system to operational forecast for CFS
!!********************************!!******************************************************* 
subroutine ReadCOPpara(epp,CFSData1)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(CFSData) ::  CFSData1
character(c100) ::   tempc,tvalue
integer  i,k
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
 
  open(unit=10,file=epp.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.13) then

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.CFile')) read(tvalue,*) CFSData1.CFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.DFile')) read(tvalue,*) CFSData1.DFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.OutFile')) read(tvalue,*) CFSData1.OutFile 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.EName')) read(tvalue,*) CFSData1.EName 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.leadT')) read(tvalue,*) CFSData1.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.BeginYear')) read(tvalue,*) CFSData1.BeginYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.EndYear')) read(tvalue,*) CFSData1.EndYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.ifile')) read(tvalue,*) CFSData1.ifile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.MName')) read(tvalue,*) CFSData1.MName
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.iYear'))   read(tvalue,*) CFSData1.iYear
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.imonth'))   read(tvalue,*) CFSData1.imonth
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.iday'))   read(tvalue,*) CFSData1.iday
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.ihour'))   read(tvalue,*) CFSData1.ihour
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.ensm'))   read(tvalue,*) CFSData1.ensm
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData4.WFile')) read(tvalue,*) CFSData1.WFile



        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.EpFile')) read(tvalue,*) epp.EpFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.ObFile')) read(tvalue,*) epp.ObFile
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.SiFile')) read(tvalue,*) epp.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.SiFile1')) read(tvalue,*) epp.SiFile1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.PaFile')) read(tvalue,*) epp.PaFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.flagCG')) read(tvalue,*) epp.flagCG
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.flagT')) read(tvalue,*) epp.flagT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.leadT')) read(tvalue,*) epp.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.leadT1')) read(tvalue,*) epp.leadT1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.nparint')) read(tvalue,*) epp.nparint
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.x1')) read(tvalue,*) epp.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.x2')) read(tvalue,*) epp.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.y1')) read(tvalue,*) epp.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.y2')) read(tvalue,*) epp.y2
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.DoFile')) read(tvalue,*) epp.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.ReFile')) read(tvalue,*) epp.ReFile

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.BeginYear ')) read(tvalue,*) epp.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.EndYear')) read(tvalue,*) epp.EndYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.BeginYearEns ')) read(tvalue,*) epp.BeginYearEns 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.EndYearEns')) read(tvalue,*) epp.EndYearEns 

 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.iobstran ')) read(tvalue,*) epp.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.ifcsttran ')) read(tvalue,*) epp.ifcsttran 
      
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.Events')) then
          read(tvalue,*) epp.Events 
          allocate(epp.Estart(epp.Events))
          allocate(epp.Estop(epp.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.Estart'))read(10,*)k,(epp.Estart(i),i=1,epp.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp2.Estop'))read(10,*)k,(epp.Estop(i),i=1,epp.Events)  
      endif
    end do 
  10 close(10) 
  epp.nems=epp.EndYearEns-epp.BeginYearEns+1
end subroutine

!!********************************!!******************************************************* 
!89 Calculate CFS four neighbours weight
!!********************************!!******************************************************* 
subroutine calweight1(CFSdata1)
use PrecStru
implicit none
integer i,j,k,flagOs
type(CFSdata) ::  CFSdata1,CFSdata2 
type(ASCFile) ::  obw
real,pointer :: w(:,:,:)
real sump
character(C1000)  asFile,tempc 
real(r8)  x1,y1,xx,yy,wx1,wy1,wx2,wy2 
real(r8),pointer ::   x(:)  
real(r8),pointer ::   y(:)  
    flagOs=CFSData1.flagOs
    if (flagOs .eq. 2)   call slashchange(CFSData1.MName)
    CFSData2.CFile=CFSData1.MName
    call ReadCFScoord(CFSData2)    !CFSData2.ncols,CFSData2.nrows  
    CFSData1.obw.nrows=CFSData2.nrows
    CFSData1.obw.ncols=1
    allocate(CFSData1.w(1,CFSData2.ncols,4))
	w=>CFSData1.w
 
    allocate(CFSData1.position(1,CFSData2.ncols))
    
    do i=1,CFSData2.ncols
      j=1
      xx=CFSData2.x(i)   
      if (xx.lt.0) xx=xx+360
      do while (CFSData1.x(j).le.xx .and. j.le.CFSData1.ncols)
      j=j+1
      end do
      CFSData1.position(1,i).x=j
    end do
    do i=1,CFSData2.nrows
      j=1
      yy=CFSData2.y(i) 
      do while (CFSData1.y(j).le.yy .and. j.le.CFSData1.nrows)
      j=j+1
      end do
      CFSData1.position(1,i).y=j               
    end do

    ! USA 0.125: X1=X0+i*0.125, Y1=Y0+j*0.125    
    ! GFS 2.5:   X2=X0+i*2.5,   Y2=Y0+j*2.5
    do i=1,1
      do j=1,CFSData2.nrows    
        sump=0     
        wx1= CFSData2.x(j)  
        if (wx1.lt.0.d0) wx1=wx1+360.0 
        wy1= CFSData2.y(j)

        wx2=CFSData1.x(CFSData1.position(i,j).x-1)
        wy2=CFSData1.y(CFSData1.position(i,j).y-1)        
        w(i,j,1)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,1)>0.0001) then 
          w(i,j,1)=1/w(i,j,1)
        else
          w(i,j,1)=9999.0 
        end if
        sump=sump+w(i,j,1)
        if (CFSData1.position(i,j).x.le.CFSData1.ncols) then
          wx2=CFSData1.x(CFSData1.position(i,j).x)
        else 
          wx2=CFSData1.x(1)          
        end if 
        w(i,j,2)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,2)>0.0001) then 
          w(i,j,2)=1/w(i,j,2)
        else
          w(i,j,2)=9999.0 
        end if
        sump=sump+w(i,j,2)
        if (CFSData1.position(i,j).y.le.CFSData1.nrows) then
          wy2=CFSData1.y(CFSData1.position(i,j).y)
        else 
          wy2=CFSData1.y(1)          
        end if                
        w(i,j,3)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,3)>0.0001) then 
          w(i,j,3)=1/w(i,j,3)
        else
          w(i,j,3)=9999.0 
        end if
        sump=sump+w(i,j,3)
        wx2=CFSData1.x(CFSData1.position(i,j).x-1)       
        w(i,j,4)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,4)>0.0001) then 
          w(i,j,4)=1/w(i,j,4)
        else
          w(i,j,4)=9999.0 
        end if
        sump=sump+w(i,j,4)
        
        do k=1,4
          w(i,j,k)=w(i,j,k)/sump
        end do  
      end do
    end do
    if (flagOs .eq. 2)   call slashchange(CFSData1.WFile)
    open(10,file=CFSData1.WFile)  
    do j=1,CFSData2.nrows
      write(10,'(i3,a,4(f6.4,a))')j,char(9),((w(1,j,k),char(9)), k=1,4)
	  !write(CFSData1.ifile,rec=k)i,j,(int(CFSdata1.PValue(l)*100),l=1,CFSdata1.leadT)
    end do
    close(10)

    
  
end subroutine 

!!********************************!!******************************************************* 
!89 Calculate CFS four neighbours weight
!!********************************!!******************************************************* 
subroutine calweight(CFSdata1)
use PrecStru
implicit none
integer i,j,k,flagOs
type(CFSdata) ::  CFSdata1 
type(ASCFile) ::  obw
real,pointer :: w(:,:,:)
real sump
character(C1000)  asFile,tempc 
real(r8)  x1,y1,xx,yy,wx1,wy1,wx2,wy2 
real(r8),pointer ::   x(:)  
real(r8),pointer ::   y(:)  
    flagOs=CFSData1.flagOs
    if (flagOs .eq. 2)   call slashchange(CFSData1.MName)
    call  ReadASC(obw,CFSData1.MName)
    CFSData1.obw=obw

    allocate(CFSData1.w(obw.ncols,obw.nrows,4))
	w=>CFSData1.w
 
    allocate(CFSData1.position(obw.ncols,obw.nrows))
    do i=1,obw.ncols
      j=1
      xx=obw.xllcorner + (i-1)*obw.cellsize
      if (xx.lt.0) xx=xx+360
      do while (CFSData1.x(j).le.xx .and. j.le.CFSData1.ncols)
      j=j+1
      end do
      CFSData1.position(i,:).x=j
    end do
    do i=1,obw.nrows
      j=1
      yy=obw.yllcorner + (i-1)*obw.cellsize 
      do while (CFSData1.y(j).le.yy .and. j.le.CFSData1.nrows)
      j=j+1
      end do
      CFSData1.position(:,i).y=j               
    end do

    ! USA 0.125: X1=X0+i*0.125, Y1=Y0+j*0.125    
    ! GFS 2.5:   X2=X0+i*2.5,   Y2=Y0+j*2.5
    do i=1,obw.ncols
      do j=1,obw.nrows    

        sump=0     
        wx1= obw.xllcorner +(i-1)*obw.cellsize
        if (wx1.lt.0.d0) wx1=wx1+360.0 
        wy1= obw.yllcorner +(j-1)*obw.cellsize

        wx2=CFSData1.x(CFSData1.position(i,j).x-1)
        wy2=CFSData1.y(CFSData1.position(i,j).y-1)        
        w(i,j,1)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,1)>0.0001) then 
          w(i,j,1)=1/w(i,j,1)
        else
          w(i,j,1)=9999.0 
        end if
        sump=sump+w(i,j,1)
        if (CFSData1.position(i,j).x.le.CFSData1.ncols) then
          wx2=CFSData1.x(CFSData1.position(i,j).x)
        else 
          wx2=CFSData1.x(1)          
        end if 
        w(i,j,2)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,2)>0.0001) then 
          w(i,j,2)=1/w(i,j,2)
        else
          w(i,j,2)=9999.0 
        end if
        sump=sump+w(i,j,2)
        if (CFSData1.position(i,j).y.le.CFSData1.nrows) then
          wy2=CFSData1.y(CFSData1.position(i,j).y)
        else 
          wy2=CFSData1.y(1)          
        end if                
        w(i,j,3)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,3)>0.0001) then 
          w(i,j,3)=1/w(i,j,3)
        else
          w(i,j,3)=9999.0 
        end if
        sump=sump+w(i,j,3)
        wx2=CFSData1.x(CFSData1.position(i,j).x-1)       
        w(i,j,4)= (wx1-wx2)*(wx1-wx2)+(wy1-wy2)*(wy1-wy2)
        if (w(i,j,4)>0.0001) then 
          w(i,j,4)=1/w(i,j,4)
        else
          w(i,j,4)=9999.0 
        end if
        sump=sump+w(i,j,4)
        
        do k=1,4
          w(i,j,k)=w(i,j,k)/sump
        end do  
      end do
    end do
    if (flagOs .eq. 2)   call slashchange(CFSData1.WFile)
    open(10,file=CFSData1.WFile,form='unformatted',access='direct',recl=obw.ncols*obw.nrows*4)  
      write(10,rec=1)(((w(i,j,k),i=1,obw.ncols),j=1,obw.nrows), k=1,4)
	  !write(CFSData1.ifile,rec=k)i,j,(int(CFSdata1.PValue(l)*100),l=1,CFSdata1.leadT)

    close(10)

    
  
end subroutine 

!!********************************!!******************************************************* 
!88 Read CFS operational File  
!!********************************!!******************************************************* 
subroutine ReadCFSop(CFSdata1)
use PrecStru
implicit none
integer i,j,k ,ii,jj
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C1000)  asFile,tempc 
character(8) Tdate
integer JULDAY,jday1
logical alive
   
  !jday1=JULDAY(1,1,CFSdata1.iYear) 
  !caldat1.julian=jday1+CFSdata1.cday-1
  !call caldat(caldat1)  
  !CFSdata1.imonth=caldat1.mm
  !CFSdata1.iday=caldat1.id
  alive=.false.
  k=CFSdata1.ihour
  write(Tdate,'(i4,2i2.2)') CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday
  do while (.not.alive.and.k.lt.20)
    write(asFile,'(4a,i2.2,2a,i2.2,3a,i2.2,2a,i2.2,a)') trim(CFSdata1.DFile),'.',Tdate,'\',k,'\',&
		'time_grib_',CFSdata1.ensm,'\',trim(CFSdata1.EName),'.',CFSdata1.ensm,'.',Tdate,k,'.daily.grb2'	  
	if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
    inquire(file=asFile,exist=alive)
    if (.not.alive )k=k+6
  end do
  CFSdata1.alive=alive
 
  if (.not.alive.or.k.gt.18) then     
    print *,trim(asFile),' is not exist.'
    return
  end if
  if (CFSData1.flagOs.eq.2) then 
    if (CFSdata1.EName.eq."prate") then
      write(asFile,'(7a,i2.2)') './wgrib2  ',trim(asFile),' -bin ',&
        trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile    
    else
      write(asFile,'(6a,i2.2)') './wgrib2  ',trim(asFile),' -bin ',&
        trim(CFSdata1.OutFile), 'Tmax\','fort.',CFSData1.ifile    

    endif 		     
  else
    read(asFile,'(a2,a)') tempc,asFile
    write(asFile,'(a,a1,2a,i2.2)')'wgrib2  \cygdrive\',tempc,trim(asFile),' -bin fort.',CFSData1.ifile   
  end if	       
   
  call slashchange(asFile)
  call system(asFile)
  if (CFSData1.flagOs.eq.2) then
    if (CFSdata1.EName.eq."prate") then
      write(asFile,'(4a,i2)')trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile        
    else
      write(asFile,'(3a,i2)')trim(CFSdata1.OutFile),'Tmax\','fort.',CFSData1.ifile        
    end if
    call slashchange(asFile)
  else
    write(asFile,'(a,i2)') 'fort.',CFSData1.ifile       
  end if
  open(unit=20,file=asFile,form='unformatted',access='sequential')

  j=CFSdata1.leadT*4  
  CFSdata1.MValue=>CFSdata1.MValue00
  do i=1,j     
    write(*,*)"i=",i
    read(20,end=20) CFSdata1.MValue(i,:,:) 
  end do
  20 close(20)
  ! 6 hours data to daily data
  if (trim(CFSdata1.EName).eq.'prate')  then
    do k=2,4-CFSdata1.ihour/6    
      CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) +CFSdata1.MValue(k,:,:) 
    end do
    CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) *3600*6
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =0
      do k=1,4
        j=(i-1)*4+k-CFSdata1.ihour/6
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(j,:,:) 
      end do
      CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) *3600*6
    end do 
  end if


  if (trim(CFSdata1.EName).eq.'tmp2m') then

    do k=2,4    
      CFSdata1.MValue1(1,:,:) = CFSdata1.MValue(1,:,:)  
      do ii=1,CFSData1.ncols   
        do jj=1,CFSData1.nrows
          if (CFSdata1.MValue(1,ii,jj).lt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
          if (CFSdata1.MValue1(1,ii,jj).gt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue1(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
        end do
      end do
    end do
    
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =-999
      CFSdata1.MValue1(i,:,:) =999
      do k=1,4
        j=(i-1)*4+k
        do ii=1,CFSData1.ncols   
          do jj=1,CFSData1.nrows
            if (CFSdata1.MValue(i,ii,jj).lt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            if (CFSdata1.MValue1(i,ii,jj).gt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue1(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            
          end do
        end do
      end do       
    end do 

    
  end if
end subroutine ReadCFSop

!!********************************!!******************************************************* 
!87 Read  Date binary CFS File  on one point four members
!!********************************!!******************************************************* 
subroutine ReadDCFSb1(CFSdata1)
use PrecStru
implicit none
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C100)  asFile,tempc 
integer JULDAY,jday1
logical alive
integer i,j,k 
integer ix,iy
   
  do k=0,18,6
    alive=.false.
    write(asFile,'(a,i2.2,a,i4,i2.2,i2.2,a)') trim(CFSdata1.DFile),k,'\prate\',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,'.bin'    
    if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
    inquire(file=asFile,exist=alive) 
    CFSdata1.alive=alive
    if (.not.alive ) then     
      print *,trim(asFile),' is not exist.'
      stop
    end if  
    open(unit=30,file=asFile,form='unformatted',access='direct',recl=CFSdata1.leadT+2 )
      ix=0 !CFSdata1.obw.ncols+1
      iy=0
      do while (ix.lt.CFSdata1.ix.or.iy.lt.CFSdata1.iy)
        read(30,rec=CFSdata1.Crow)ix,iy
        CFSdata1.Crow=CFSdata1.Crow+1
      end do
      CFSdata1.Crow=CFSdata1.Crow-1
	  CFSdata1.Crow1=CFSdata1.Crow
      if (ix.eq.CFSdata1.ix.and.iy.eq.CFSdata1.iy) then
        read(30,rec=CFSdata1.Crow)ix,iy,(CFSdata1.PLValue(CFSdata1.cday,k/6+1,i),i=1,CFSdata1.leadT)	   
      else
        print *,ix,iy,' is not exist.'
        stop
      end if	 
    close (30)
  end do
end subroutine ReadDCFSb1
!!********************************!!******************************************************* 
!86 Read  Date binary CFS File  on one forecast daily
!!********************************!!******************************************************* 
subroutine ReadDCFSb(CFSdata1)
use PrecStru
implicit none
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C100)  asFile,tempc 
integer JULDAY,jday1
logical alive
integer i,j,k 
integer ix,iy
   

  alive=.false.
  write(asFile,'(a,i4,i2.2,i2.2,a)') trim(CFSdata1.DFile),CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,'.bin'    
  if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
  inquire(file=asFile,exist=alive)
    
  !write(*,*)trim(asFile)
  CFSdata1.alive=alive
  if (.not.alive ) then     
    print *,trim(asFile),' is not exist.'
    stop
  end if
  
  open(unit=30,file=asFile,form='unformatted',access='direct',recl=CFSdata1.leadT+2 )
    ix=0 !CFSdata1.obw.ncols+1
    do while (ix.lt.CFSdata1.ix)
      read(30,rec=CFSdata1.Crow)ix
      CFSdata1.Crow=CFSdata1.Crow+1
    end do
    CFSdata1.Crow=CFSdata1.Crow-1
	CFSdata1.Crow1=CFSdata1.Crow
    do while (ix.eq.CFSdata1.ix) !.and.CFSdata1.Crow.le.CFSdata1.obw.SumN       
      read(30,rec=CFSdata1.Crow)ix,iy,(CFSdata1.PLValue(CFSdata1.cday,iy,i),i=1,CFSdata1.leadT)
      CFSdata1.Crow=CFSdata1.Crow+1
      if (CFSdata1.Crow.le.CFSdata1.obw.SumN) then 
        read(30,rec=CFSdata1.Crow)ix
      else
	    ix=0
      end if
    end do 
    !CFSdata1.Crow=CFSdata1.Crow-1
    !read(30,rec=k)ix,iy,(CFSdata1.PLValue(CFSdata1.cday,iy,i),i=1,CFSdata1.leadT)

  close (30)
 
end subroutine ReadDCFSb

!!********************************!!******************************************************* 
!85 binary File to AscII File  
!!********************************!!******************************************************* 

subroutine BinToAdata(outputfile)
use PrecStru
implicit none
type(Scaldat) ::  caldat1
character(c100) ::   outputfile
integer   i,j,pv(282)
 
  open (10,file=outputfile,form='unformatted',access='direct',recl=282) !
  write(outputfile,'(a,a)')trim(outputfile),'a.txt'  
  open (20,file=outputfile)    
   i=1
    do  i=1,73*280*27    
      read (10 ,rec=i) pv    
      write (20,*) pv         
    end do
  100 close(10)
  close(20)

end subroutine
!!********************************!!******************************************************* 
! 84 Write  hindcast result -binary
!!********************************!!******************************************************* 
subroutine Writehincastb(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(Scaldat) ::  caldat1
character(c100) ::   outputfile
integer  i,j,k,jday1,imem,JULDAY ,iyear 
character(c100) asRecord 
character(10) tempc  
  write(asRecord,*)'Begin to write hindcast'
  jday1=JULDAY(1,1,epp.BeginYear)  
  write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(epp.ReFile),epp.x,'\',epp.y,'.bin' 
  call Wrecord(outputfile)
  if (epp.flagOs .eq. 2)   call slashchange(outputfile)  
  open (10,file=outputfile,status='REPLACE',form='unformatted',access='direct',recl=epp.nems+7)!
  k=0  
    do i=1,epp.ndays
      if (epp.Pens(i,1,1).ge.0.) then	 
        caldat1.julian=jday1+i-1
        call caldat(caldat1)   
        do j=1,epp.leadT     
		  k=k+1    
          write (10,rec=k)caldat1.iyyy,caldat1.mm,caldat1.id,j,int(epp.ob(i,j)*100),&  
            int(epp.si(i,j)*100),int(epp.Pens(i,j,epp.nems+1)*100),&        
            (int(epp.Pens(i,j,imem)*100),imem=1,epp.nems)        
        end do !j
      endif
    end do    
  close(10)
  !output average 
  write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(epp.ReFile),epp.x,'\',epp.y,'m.txt' 
  if (epp.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) 
    write (10,'(7a,\)')'Year',char(9),'Month',char(9),'Day',char(9),'Observed'	
    do i=1,epp.leadT,epp.nparint1
	  write (10,'(2a,i3.3,\)')char(9),'Fcst',i
    end do  
    do i=1,epp.leadT,epp.nparint1
	  write (10,'(2a,i3.3,\)')char(9),'Ensmean',i
    end do
    write (10,*)
    do i=365,epp.ndays 
	  if (epp.Pens(i,1,epp.nems+1).gt.0) then
        caldat1.julian=jday1+i-1
        call caldat(caldat1)  
        iyear= caldat1.iyyy
        write(tempc,'(i7)') int(epp.ob(i,1)*100)   
        tempc=adjustl(tempc)
        write (10,'(i4,a,i2,a,i2,2a,\)')caldat1.iyyy,char(9),caldat1.mm,char(9),caldat1.id,char(9),trim(tempc)	  	   
        do j=1,epp.leadT,epp.nparint1
          k=i-j+1
          if (k>0) then
            if (epp.nparint1.ne.5.and.epp.Pens(k,j,epp.nems+1).lt.0) k=k-5
	        if(epp.Pens(k,j,epp.nems+1).lt.0 ) k=k-1
            if (epp.nparint1.ne.5.and.epp.Pens(k,j,epp.nems+1).lt.0) k=k+5
            write(tempc,'(i7)') int(epp.si(k,j) *100)   
          else
            write(tempc,'(i7)') int(epp.si(j,j) *100)   
          endif
          tempc=adjustl(tempc) 
          write (10,'(2a,\)') char(9),trim(tempc) 
        end do !j
        do j=1,epp.leadT,epp.nparint1
          k=i-j+1
          if (k>0) then
            if (epp.nparint1.ne.5.and.epp.Pens(k,j,epp.nems+1).lt.0) k=k-5
	        if(epp.Pens(k,j,epp.nems+1).lt.0 ) k=k-1
            if (epp.nparint1.ne.5.and.epp.Pens(k,j,epp.nems+1).lt.0) k=k+5

            write(tempc,'(i7)') int(epp.Pens(k,j,epp.nems+1)*100)   
          else
            write(tempc,'(i7)') int(epp.Pens(j,j,epp.nems+1)*100)   
          endif
          tempc=adjustl(tempc) 
          write (10,'(2a,\)') char(9),trim(tempc) 
        end do !j
        write(10,*)
      end if
    end do    
  close(10)
end subroutine 


!!********************************!!******************************************************* 
!83 Read CFS File  4 times
!!********************************!!******************************************************* 
subroutine ReadCFS4(CFSdata1)
use PrecStru
implicit none
integer i,j,k,kk ,ii,jj
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C1000)  asFile,tempc 
integer JULDAY,jday1
logical alive
   

  jday1=JULDAY(1,1,CFSdata1.iYear) 
  caldat1.julian=jday1+CFSdata1.cday-1
  call caldat(caldat1)  
  CFSdata1.imonth=caldat1.mm
  CFSdata1.iday=caldat1.id
  alive=.false.
  k=0
  do k=0,18,6
    kk=k
    do while (.not.alive.and.kk.lt.20)
      write(asFile,'(3a,i4,i2.2,3a,i4,3i2.2,a)') trim(CFSdata1.DFile),trim(CFSdata1.EName),'\',&
        CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,&
        CFSdata1.iday,kk,'.time.grb2'    
      if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
      inquire(file=asFile,exist=alive)
      if (.not.alive )kk=kk+6
    end do
    CFSdata1.alive=alive 
    if (.not.alive) then     
      print *,trim(asFile),' is not exist.'      
      return
    end if
    if (CFSData1.flagOs.eq.2) then 
      write(asFile,'(4a,i4,i2.2,3a,i4,3i2.2,5a,i2)')'./wgrib2  ',trim(CFSdata1.DFile),trim(CFSdata1.EName),'\',&
        CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	    kk,'.time.grb2 -bin ',trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile   
    else
      read(CFSdata1.DFile,'(a2,a)') tempc,asFile
      write(asFile,'(a,a1,3a,i4,i2.2,3a,i4,3i2.2,a,i2)')'wgrib2  \cygdrive\',tempc,trim(asFile),trim(CFSdata1.EName),'\',&
        CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
        kk,'.time.grb2 -bin fort.',CFSData1.ifile   
    end if	       

    call slashchange(asFile)
    call system(asFile)
    if (CFSData1.flagOs.eq.2) then
      write(asFile,'(4a,i2)')trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile        
      call slashchange(asFile)
    else
      write(asFile,'(a,i2)') 'fort.',CFSData1.ifile        

    end if
    open(unit=20,file=asFile,form='unformatted',access='sequential')

    j=CFSdata1.leadT*4  
    if (k.eq.0) CFSdata1.MValue=>CFSdata1.MValue00
    if (k.eq.6) CFSdata1.MValue=>CFSdata1.MValue06
    if (k.eq.12) CFSdata1.MValue=>CFSdata1.MValue12
    if (k.eq.18) CFSdata1.MValue=>CFSdata1.MValue18
    do i=1,j     
      write(*,*)"i=",i
      read(20,end=20) CFSdata1.MValue(i,:,:) 
    end do
    20 close(20)
  end do 

  ! 6 hours data to daily data
  if (trim(CFSdata1.EName).eq.'prate')  then
    do i=1,CFSdata1.leadT*4
      CFSdata1.MValue(i,:,:)=CFSdata1.MValue00(i,:,:) +CFSdata1.MValue06(i,:,:) &
        +CFSdata1.MValue12(i,:,:)+CFSdata1.MValue18(i,:,:) 
      CFSdata1.MValue(i,:,:)=CFSdata1.MValue(i,:,:)*3600*6
    end do
    !CFSdata1.MValue(1,:,:)=CFSdata1.MValue00(1,:,:) +CFSdata1.MValue06(1,:,:)+ &
    !  CFSdata1.MValue12(1,:,:)+CFSdata1.MValue18(1,:,:)
    !do k=2,4-CFSdata1.ihour/6    
    !  CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) +CFSdata1.MValue(k,:,:) 
    !end do
    !CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) *3600*6
    !CFSdata1.MValue=CFSdata1.MValue*3600
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      !CFSdata1.MValue(i,:,:) =0
      !do k=1,4
      j=(i-1)*4  !+k-CFSdata1.ihour/6
      CFSdata1.MValue(i,:,:) = CFSdata1.MValue(j,:,:) ! only use one day data each 4 time interval
      !end do
      !CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) *3600*6
    end do 
  end if


  if (trim(CFSdata1.EName).eq.'tmp2m') then

    do k=2,4    
      CFSdata1.MValue1(1,:,:) = CFSdata1.MValue(1,:,:)  
      do ii=1,CFSData1.ncols   
        do jj=1,CFSData1.nrows
          if (CFSdata1.MValue(1,ii,jj).lt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
          if (CFSdata1.MValue1(1,ii,jj).gt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue1(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
        end do
      end do
    end do
    
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =-999
      CFSdata1.MValue1(i,:,:) =999
      do k=1,4
        j=(i-1)*4+k
        do ii=1,CFSData1.ncols   
          do jj=1,CFSData1.nrows
            if (CFSdata1.MValue(i,ii,jj).lt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            if (CFSdata1.MValue1(i,ii,jj).gt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue1(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            
          end do
        end do
      end do       
    end do 

    
  end if
end subroutine ReadCFS4


!!********************************!!******************************************************* 
!82 read observed  data  for shaake schuffle
!!********************************!!******************************************************* 
subroutine readDataEpT(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays
integer month1,month2,day1,day2 
integer julday,jday1,jday2
character*4 units_bmap06
real xobs(4) 
  if (rData1.flagOs.eq.2) then
    if (rData1.flagT.eq.1) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.EpFile),rData1.x,'/',rData1.y,'Tmax.txt' 
    if (rData1.flagT.eq.2) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.EpFile),rData1.x,'/',rData1.y,'Tmin.txt' 
  else
    if (rData1.flagT.eq.1) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.EpFile),rData1.x,'\',rData1.y,'Tmax.txt'  
    if (rData1.flagT.eq.2) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.EpFile),rData1.x,'\',rData1.y,'Tmin.txt'  
  end if
  open (10,file=inputfile)
    read (10,*) year1,month1,day1,year2,month2,day2
    jday1=julday(month1,day1,year1)
    jday2=julday(1,1,rData1.BeginYearEns)
     
    if (rData1.BeginYearEns.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYearEns,'<',year1
	  stop
    end if
    if (rData1.EndYearEns.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYearEns,'>',year2
	  stop
    end if
    j=jday2-jday1
    do i=1,j
      read (10,*)        
    end do 
    !if (rData1.nparint.eq.1)  then
    !  ndays=rData1.ndays      
    !else
    ndays=julday(12,31,rData1.EndYearEns)-jday2 +rData1.leadT
    !end if
    do i=1,ndays
      read (10,*) rData1.ep(i,1)
      !rData1.ep(i,1)=rData1.ep(i,1)  !-200
      !rData1.ob(i,1)=max(0.0,xobs(1)+xobs(2)+xobs(3)+xobs(4))
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ep(i,j)=rData1.ep(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ep(i,j)=rData1.ep(i,j-1)
      end do
    end do
  close(10)

end subroutine

!!********************************!!******************************************************* 
!81 read observed  data  
!!********************************!!******************************************************* 
subroutine readDataobT(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays
integer month1,month2,day1,day2 
integer julday,jday1,jday2
character*4 units_bmap06
real xobs(4) 
  if (rData1.flagOs.eq.2) then
    if (rData1.flagT.eq.1) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.ObFile),rData1.x,'/',rData1.y,'Tmax.txt' 
    if (rData1.flagT.eq.2) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.ObFile),rData1.x,'/',rData1.y,'Tmin.txt' 
  else
    if (rData1.flagT.eq.1) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.ObFile),rData1.x,'\',rData1.y,'Tmax.txt'  
    if (rData1.flagT.eq.2) write(inputfile,'(a,i3.3,a1,i3.3,a)')trim(rData1.ObFile),rData1.x,'\',rData1.y,'Tmin.txt'  
  end if
  open (10,file=inputfile)
    read (10,*) year1,month1,day1,year2,month2,day2
    jday1=julday(month1,day1,year1)
    jday2=julday(1,1,rData1.BeginYear)
 
    if (rData1.BeginYear.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYear,'<',year1
	  stop
    end if
    if (rData1.EndYear.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYear,'>',year2
	  stop
    end if
    j=jday2-jday1
    do i=1,j
      read (10,*)        
    end do 
    !if (rData1.nparint.eq.1)  then
    !  ndays=rData1.ndays      
    !else
    ndays=rData1.ndays+rData1.leadT
    !end if
    do i=1,ndays
      read (10,*) rData1.ob(i,1)
      !rData1.ob(i,1)=rData1.ob(i,1)  !-200
      !rData1.ob(i,1)=max(0.0,xobs(1)+xobs(2)+xobs(3)+xobs(4))
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ob(i,j)=rData1.ob(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ob(i,j)=rData1.ob(i,j-1)
      end do
    end do
  close(10)

    end subroutine

!!********************************!!******************************************************* 
!80 IDW interpolate CFS to 1/8 degree
!!********************************!!******************************************************* 
subroutine CFSinterpolate1(CFSData1)
use PrecStru
implicit none 
type(CFSData) ::  CFSData1
type(ASCFile) ::  obw
integer i,j,k,ii,jj,l  
character(c100) ::  outputfile 
!real   w(464,224,4)
 real  ,pointer :: w(:,:,:) 

  w=>CFSData1.w
  obw=CFSData1.obw
  print*,'year,day= ',CFSData1.iyear,CFSData1.cday ,CFSData1.imonth ,CFSData1.iday  
  write(outputfile,'(3a,i4,2i2.2,a)')trim(CFSData1.OutFile),trim(CFSdata1.EName),'\',&
    CFSData1.iyear,CFSData1.imonth,CFSData1.iday,'.bin'        
  if (CFSData1.flagOs .eq. 2)   call slashchange(outputfile)
  open (CFSData1.ifile,file=trim(outputfile),status='REPLACE',form='unformatted',access='direct',recl=CFSdata1.leadT+2)!
  !open(CFSData1.ifile,file=trim(outputfile))
    k=0
	CFSdata1.MValueN =0
    do i=1,obw.ncols
      do j=1,obw.nrows
        
          do l=1,CFSdata1.leadT
            ! MValue(l,:,:) CFS file raw value (lead time, lon, lat) CFS file point value (lead time)
            ii=CFSData1.position(i,j).x-1
            jj=CFSData1.position(i,j).y-1
            CFSdata1.PValue(l)=CFSdata1.MValue(l,ii,jj)*w(i,j,1)    
            if (CFSData1.position(i,j).x.le.CFSData1.ncols) then
              ii=CFSData1.position(i,j).x
            else
              ii=1
            endif
            CFSdata1.PValue(l)=CFSdata1.PValue(l)+CFSdata1.MValue(l,ii,jj)*w(i,j,2)    
            if (CFSData1.position(i,j).y.le.CFSData1.nrows) then
              jj=CFSData1.position(i,j).y
            else
              jj=1
            endif
            CFSdata1.PValue(l)=CFSdata1.PValue(l)+CFSdata1.MValue(l,ii,jj)*w(i,j,3)    
            ii=CFSData1.position(i,j).x-1
            CFSdata1.PValue(l)=CFSdata1.PValue(l)+CFSdata1.MValue(l,ii,jj)*w(i,j,4)    
            !write(CFSdata1.PIValue(l),'(i6)')int(CFSdata1.PValue(l)*100)
            !CFSdata1.PIValue(l)=adjustl(CFSdata1.PIValue(l))
            CFSdata1.MValueN(l,i,j)=CFSdata1.PValue(l)
          end do 
          !MyASC.MValue(i,j)=MyASC.MValue(i,j)+CFSdata1.PValue(1)
          !2.4. output 0.125 degree CFS data
          !write(CFSData1.ifile,'(i3,a,i3,a,i4,a,i2.2,a,i2.2,800a)')i,char(9),j,char(9),&
          !CFSData1.iyear,char(9),CFSData1.imonth,char(9),CFSData1.iday,(char(9),&
          !trim(CFSdata1.PIValue(l)),l=1,CFSdata1.leadT)
          k=k+1
          write(CFSData1.ifile,rec=k)i,j,(int(CFSdata1.PValue(l)*100),l=1,CFSdata1.leadT)
         
      end do
    end do
  close(CFSData1.ifile)

end subroutine    
    
!!********************************!!******************************************************* 
!80 IDW interpolate CFS to 1/8 degree
!!********************************!!******************************************************* 
subroutine CFSinterpolate(CFSData1)
use PrecStru
implicit none 
type(CFSData) ::  CFSData1
type(ASCFile) ::  obw
integer i,j,k,ii,jj,l  
character(c100) ::  outputfile 
real   w(464,224,4)

  w=CFSData1.w
  obw=CFSData1.obw
  print*,'year,day= ',CFSData1.iyear,CFSData1.cday ,CFSData1.imonth ,CFSData1.iday  
  write(outputfile,'(3a,i4,2i2.2,a)')trim(CFSData1.OutFile),trim(CFSdata1.EName),'\',&
    CFSData1.iyear,CFSData1.imonth,CFSData1.iday,'.bin'        
  if (CFSData1.flagOs .eq. 2)   call slashchange(outputfile)
  open (CFSData1.ifile,file=trim(outputfile),status='REPLACE',form='unformatted',access='direct',recl=CFSdata1.leadT+2)!
  !open(CFSData1.ifile,file=trim(outputfile))
    k=0
	CFSdata1.MValueN =0
    do i=1,obw.ncols
      do j=1,obw.nrows
        if (abs(obw.MValue(i,j)-obw.NODATA_value)>0.01) then
          do l=1,CFSdata1.leadT
            ! MValue(l,:,:) CFS file raw value (lead time, lon, lat) CFS file point value (lead time)
            ii=CFSData1.position(i,j).x-1
            jj=CFSData1.position(i,j).y-1
            CFSdata1.PValue(l)=CFSdata1.MValue(l,ii,jj)*w(i,j,1)    
            if (CFSData1.position(i,j).x.le.CFSData1.ncols) then
              ii=CFSData1.position(i,j).x
            else
              ii=1
            endif
            CFSdata1.PValue(l)=CFSdata1.PValue(l)+CFSdata1.MValue(l,ii,jj)*w(i,j,2)    
            if (CFSData1.position(i,j).y.le.CFSData1.nrows) then
              jj=CFSData1.position(i,j).y
            else
              jj=1
            endif
            CFSdata1.PValue(l)=CFSdata1.PValue(l)+CFSdata1.MValue(l,ii,jj)*w(i,j,3)    
            ii=CFSData1.position(i,j).x-1
            CFSdata1.PValue(l)=CFSdata1.PValue(l)+CFSdata1.MValue(l,ii,jj)*w(i,j,4)    
            !write(CFSdata1.PIValue(l),'(i6)')int(CFSdata1.PValue(l)*100)
            !CFSdata1.PIValue(l)=adjustl(CFSdata1.PIValue(l))
            CFSdata1.MValueN(l,i,j)=CFSdata1.PValue(l)
          end do 
          !MyASC.MValue(i,j)=MyASC.MValue(i,j)+CFSdata1.PValue(1)
          !2.4. output 0.125 degree CFS data
          !write(CFSData1.ifile,'(i3,a,i3,a,i4,a,i2.2,a,i2.2,800a)')i,char(9),j,char(9),&
          !CFSData1.iyear,char(9),CFSData1.imonth,char(9),CFSData1.iday,(char(9),&
          !trim(CFSdata1.PIValue(l)),l=1,CFSdata1.leadT)
          k=k+1
          write(CFSData1.ifile,rec=k)i,j,(int(CFSdata1.PValue(l)*100),l=1,CFSdata1.leadT)
        end if 
      end do
    end do
  close(CFSData1.ifile)

end subroutine

!!********************************!!******************************************************* 
!79 Write temperature map statistical File!!**********************************
!!********************************!!******************************************************* 
subroutine WritemapstatsT(Obsdata1)
use PrecStru
implicit none 
type(Obsdata) ::  Obsdata1
character(c100) ::   outputfile
character(c1000) ::  tempc
integer  jmo,ipd  

  write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(Obsdata1.BName06),Obsdata1.x,'\',Obsdata1.y,'s.txt'  
  if (Obsdata1.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile,access='append') !  
    !write (10,'(2i5)')  Obsdata1.x,Obsdata1.y
    if (Obsdata1.CuVa.eq.2) &
      write (10,'(11a)')"TmaxMonth",char(9),"Nobs",char(9),"avg",char(9),"pop",char(9),"cavg",char(9),"ccv"
    if (Obsdata1.CuVa.eq.3) &
      write (10,'(11a)')"TminMonth",char(9),"Nobs",char(9),"avg",char(9),"pop",char(9),"cavg",char(9),"ccv"
    do jmo=1,12
      write (tempc,'(2(i5,a),4(f10.4,a))') jmo,char(9),Obsdata1.nobs(5,jmo),char(9),&
                            Obsdata1.avg(5,jmo),char(9),Obsdata1.pop(5,jmo),char(9), &
                           Obsdata1.cavg(5,jmo),char(9),Obsdata1.ccv(5,jmo) 
      call delblank(tempc)
      write (10,'(a)')trim(tempc)
    enddo
    

  close(10)    

end subroutine 

!!********************************!!******************************************************* 
!78 Write Observed daily temperature  data to  File
!!********************************!!******************************************************* 
subroutine WriteObdataT(Obsdata1)
use PrecStru
implicit none
integer JULDAY,ndays
type(Obsdata) ::  Obsdata1
character(c100) ::   outputfile
integer i,j
integer   months(12),iyear,imonth,iday     
character*4 units
character(c1000) ::  tempc
 
  units='K'
  ndays=JULDAY(12,31,Obsdata1.EndYear)-JULDAY(1,1,Obsdata1.BeginYear)+1
  if (Obsdata1.CuVa.eq.2) &
    write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(Obsdata1.BName06),Obsdata1.x,'\',Obsdata1.y,'Tmax.txt'  
  if (Obsdata1.CuVa.eq.3) &
    write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(Obsdata1.BName06),Obsdata1.x,'\',Obsdata1.y,'Tmin.txt'  
  
  if (Obsdata1.flagOs .eq. 2)   call slashchange(outputfile)
  if (Obsdata1.y.gt.Obsdata1.nrows/Obsdata1.Tim ) Obsdata1.y=Obsdata1.y-(Obsdata1.k-1)*Obsdata1.nrows/Obsdata1.Tim


  open (10,file=outputfile) !,'(i4,tr1,i4,tr1,i7,tr1,a2)'
  write(10,'(i4,5a,i4,6a)') Obsdata1.BeginYear,char(9),'1',char(9),'1',char(9),&
     Obsdata1.EndYear,char(9),'1',char(9),'1',char(9),trim(units)

  do iyear=Obsdata1.BeginYear,Obsdata1.EndYear   
    j= iyear-Obsdata1.BeginYear+1
    if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
      months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
      months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    endif 
    do imonth=1,12 
      do iday=1,months(imonth) 
        write (10,'(f6.2)')Obsdata1.Obdata(Obsdata1.x,Obsdata1.y,j,imonth,iday)      
        !call delblank(tempc)
        !write (10,'(a)')trim(tempc)          
      end do
    end do
  enddo 
  close(10)   

end subroutine 


!!********************************!!******************************************************* 
! 77 read  paramater information of system to hindcast for CFS
!!********************************!!******************************************************* 
subroutine ReadCHpara(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
character(c100) ::   tempc,tvalue
integer  i,k
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
 
  open(unit=10,file=epp.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.11) then
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.EpFile')) read(tvalue,*) epp.EpFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.ObFile')) read(tvalue,*) epp.ObFile
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.SiFile')) read(tvalue,*) epp.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.SiFile1')) read(tvalue,*) epp.SiFile1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.PaFile')) read(tvalue,*) epp.PaFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.flagCG')) read(tvalue,*) epp.flagCG
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.flagT')) read(tvalue,*) epp.flagT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.leadT')) read(tvalue,*) epp.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.leadT1')) read(tvalue,*) epp.leadT1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.nparint')) read(tvalue,*) epp.nparint
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.nparint1')) read(tvalue,*) epp.nparint1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.x1')) read(tvalue,*) epp.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.x2')) read(tvalue,*) epp.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.y1')) read(tvalue,*) epp.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.y2')) read(tvalue,*) epp.y2
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.DoFile')) read(tvalue,*) epp.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.ReFile')) read(tvalue,*) epp.ReFile

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.BeginYear ')) read(tvalue,*) epp.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.EndYear')) read(tvalue,*) epp.EndYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.BeginYearEns ')) read(tvalue,*) epp.BeginYearEns 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.EndYearEns')) read(tvalue,*) epp.EndYearEns 

 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.iobstran ')) read(tvalue,*) epp.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.ifcsttran ')) read(tvalue,*) epp.ifcsttran 
      
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.Events')) then
          read(tvalue,*) epp.Events 
          allocate(epp.Estart(epp.Events))
          allocate(epp.Estop(epp.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.Estart'))read(10,*)k,(epp.Estart(i),i=1,epp.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp1.Estop'))read(10,*)k,(epp.Estop(i),i=1,epp.Events)  
      endif
    end do 
  10 close(10) 
  epp.nems=epp.EndYearEns-epp.BeginYearEns+1
end subroutine



!!********************************!!******************************************************* 
! 76 read  paramater information of system to prepare parameters of EPP for CFS+GFS
!!********************************!!******************************************************* 
subroutine ReadCEparaCG(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   tempc,tvalue
integer  i,k
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
   
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue  
      if (k.eq.10) then
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo
        !if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.PTmin')) read(tvalue,*) CanE.PTmin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagCG')) read(tvalue,*) CanE.flagCG
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagT')) read(tvalue,*) CanE.flagT

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile1')) read(tvalue,*) CanE.SiFile1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT1')) read(tvalue,*) CanE.leadT1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint')) read(tvalue,*) CanE.nparint
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x1')) read(tvalue,*) CanE.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x2')) read(tvalue,*) CanE.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y1')) read(tvalue,*) CanE.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y2')) read(tvalue,*) CanE.y2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minall')) read(tvalue,*) CanE.minall 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minpos ')) read(tvalue,*) CanE.minpos 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.pop_fraction ')) read(tvalue,*) CanE.pop_fraction 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iobstran ')) read(tvalue,*) CanE.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ifcsttran ')) read(tvalue,*) CanE.ifcsttran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iopt_rho ')) read(tvalue,*) CanE.iopt_rho 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cor_weight ')) read(tvalue,*) CanE.cor_weight 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgobsx ')) read(tvalue,*) CanE.cavgobsx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgfcstx ')) read(tvalue,*) CanE.cavgfcstx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.verstats ')) read(tvalue,*) CanE.verstats 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint ')) read(tvalue,*) CanE.nparint 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nems ')) read(tvalue,*) CanE.nems 
    
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
          read(tvalue,*) CanE.Events 
          allocate(CanE.Estart(CanE.Events))
          allocate(CanE.Estop(CanE.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)k,(CanE.Estart(i),i=1,CanE.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)k,(CanE.Estop(i),i=1,CanE.Events)  
      endif
    end do 
  10 close(10) 
end subroutine




!!********************************!!******************************************************* 
!75 Read  Date CFS File  on one forecast daily
!!********************************!!******************************************************* 
subroutine ReadYCFS(CFSdata1)
use PrecStru
implicit none
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C100)  asFile,tempc 
integer JULDAY,jday1
logical alive
integer i,j,k 
integer ix,iy
integer Maxtrix(285,52455)   

  alive=.false.
  write(asFile,'(a,i4,i2.2,i2.2,a)') trim(CFSdata1.DFile),CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,'.txt'    
  if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
  inquire(file=asFile,exist=alive)
    
  CFSdata1.alive=alive
  if (.not.alive ) then     
    print *,trim(asFile),' is not exist.'
    stop
  end if
 
  open(unit=30,file=asFile )
  do  
    read(30,*,end=30)ix,iy,k,k,k,(CFSdata1.MValue(i,ix,iy),i=1,CFSdata1.leadT) 
  end do     
  30 close(30) 
end subroutine ReadYCFS

!!********************************!!******************************************************* 
!74 Read  Date CFS File  on one forecast daily
!!********************************!!******************************************************* 
subroutine ReadDCFS(CFSdata1)
use PrecStru
implicit none
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C100)  asFile,tempc 
integer JULDAY,jday1
logical alive
integer i,j,k 
integer ix,iy
 

  alive=.false.
  write(asFile,'(a,i4,i2.2,i2.2,a)') trim(CFSdata1.DFile),CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,'.txt'    
  if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
  inquire(file=asFile,exist=alive)
    
  !write(*,*)trim(asFile)
  CFSdata1.alive=alive
  if (.not.alive ) then     
    print *,trim(asFile),' is not exist.'
    stop
  end if

  if (CFSdata1.ix.lt.CFSdata1.ncols*0.05) then !too slowly
    open(unit=30,file=asFile )	 
      read(30,*)ix  
      do while (ix.lt.CFSdata1.ix)  
	  
        read(30,*)ix
      end do
      if (ix.eq.CFSdata1.ix) then
        backspace(30)
        !CFSdata1.nstep=1
	    do while (ix.eq.CFSdata1.ix)          
          read(30,*,end=40)ix,iy,k,k,k,(CFSdata1.PLValue(CFSdata1.cday,iy,i),i=1,CFSdata1.leadT)
          read(30,*,end=40)ix
          backspace(30)
          !CFSdata1.nstep=CFSdata1.nstep+1
        end do 
      end if
    40 close(30) 
  else
    open(unit=30,file=asFile,position='append' )
	  backspace(30)	 
      read(30,*)ix  
      do while (ix.gt.CFSdata1.ix)  
	    backspace(30)	
	    backspace(30)	
        read(30,*)ix
      end do
      if (ix.eq.CFSdata1.ix) then
        backspace(30)
        !CFSdata1.nstep=1
	    do while (ix.eq.CFSdata1.ix)
          read(30,*)ix,iy,k,k,k,(CFSdata1.PLValue(CFSdata1.cday,iy,i),i=1,CFSdata1.leadT)
	      backspace(30)	
	      backspace(30)	
          !CFSdata1.nstep=CFSdata1.nstep+1
        end do 
        read(30,*)
		ENDFILE(30)
      end if
      close(30) 
  end if
end subroutine ReadDCFS

!!********************************!!******************************************************* 
!75 delete all blank space in a string
!!********************************!!******************************************************* 
subroutine delblank(mystr)
use PrecStru
implicit none
integer i,j
character(c1000)::mystr,string1
  string1=""
  string1=adjustl(mystr)
  mystr=""
  j=1
  do i=1,len(trim(string1))
    if(string1(i:i).ne." " ) then
      mystr(j:j)=string1(i:i)
      j=j+1
    end if
  end do
end subroutine

!!********************************!!******************************************************* 
!74 Read CFS File  
!!********************************!!******************************************************* 
subroutine ReadCFS(CFSdata1)
use PrecStru
implicit none
integer i,j,k ,ii,jj
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C1000)  asFile,tempc 
integer JULDAY,jday1
logical alive
   

  jday1=JULDAY(1,1,CFSdata1.iYear) 
  caldat1.julian=jday1+CFSdata1.cday-1
  call caldat(caldat1)  
  CFSdata1.imonth=caldat1.mm
  CFSdata1.iday=caldat1.id
  alive=.false.
  k=CFSdata1.ihour
  do while (.not.alive.and.k.lt.20)
    write(asFile,'(3a,i4,i2.2,3a,i4,3i2.2,a)') trim(CFSdata1.DFile),trim(CFSdata1.EName),'\',&
      CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,&
      CFSdata1.iday,k,'.time.grb2'    
	if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
    inquire(file=asFile,exist=alive)
    if (.not.alive )k=k+6
  end do
  CFSdata1.alive=alive
 
  if (.not.alive.or.k.gt.18) then
     
    print *,trim(asFile),' is not exist.'
    return
  end if
  if (CFSData1.flagOs.eq.2) then 
    write(asFile,'(4a,i4,i2.2,3a,i4,3i2.2,5a,i2)')'./wgrib2  ',trim(CFSdata1.DFile),trim(CFSdata1.EName),'\',&
      CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	   k,'.time.grb2 -bin ',trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile   
  else
    read(CFSdata1.DFile,'(a2,a)') tempc,asFile
    write(asFile,'(a,a1,3a,i4,i2.2,3a,i4,3i2.2,a,i2)')'wgrib2  \cygdrive\',tempc,trim(asFile),trim(CFSdata1.EName),'\',&
      CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
      k,'.time.grb2 -bin fort.',CFSData1.ifile   
  end if	       
   
  call slashchange(asFile)
  call system(asFile)
  if (CFSData1.flagOs.eq.2) then
    write(asFile,'(4a,i2)')trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile        
    call slashchange(asFile)
  else
    write(asFile,'(a,i2)') 'fort.',CFSData1.ifile        

  end if
  open(unit=20,file=asFile,form='unformatted',access='sequential')

  j=CFSdata1.leadT*4  
  CFSdata1.MValue=>CFSdata1.MValue00
  do i=1,j     
    write(*,*)"i=",i
    read(20,end=20) CFSdata1.MValue(i,:,:) 
  end do
  20 close(20)
  ! 6 hours data to daily data
  if (trim(CFSdata1.EName).eq.'prate')  then
    do k=2,4-CFSdata1.ihour/6    
      CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) +CFSdata1.MValue(k,:,:) 
    end do
    CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) *3600*6
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =0
      do k=1,4
        j=(i-1)*4+k-CFSdata1.ihour/6
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(j,:,:) 
      end do
      CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) *3600*6
    end do 
  end if


  if (trim(CFSdata1.EName).eq.'tmp2m') then

    do k=2,4    
      CFSdata1.MValue1(1,:,:) = CFSdata1.MValue(1,:,:)  
      do ii=1,CFSData1.ncols   
        do jj=1,CFSData1.nrows
          if (CFSdata1.MValue(1,ii,jj).lt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
          if (CFSdata1.MValue1(1,ii,jj).gt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue1(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
        end do
      end do
    end do
    
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =-999
      CFSdata1.MValue1(i,:,:) =999
      do k=1,4
        j=(i-1)*4+k
        do ii=1,CFSData1.ncols   
          do jj=1,CFSData1.nrows
            if (CFSdata1.MValue(i,ii,jj).lt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            if (CFSdata1.MValue1(i,ii,jj).gt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue1(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            
          end do
        end do
      end do       
    end do 

    
  end if
end subroutine ReadCFS


!!********************************!!******************************************************* 
!73 Read CFS coordinate
!!********************************!!******************************************************* 
subroutine ReadCFScoord(CFSData1)
use PrecStru
implicit none
type(CFSData) ::  CFSData1
 
character(c100) asRecord
character(20) tempa 
integer i,j
logical alive
  
  inquire(file=CFSData1.CFile,exist=alive)
  if (CFSData1.flagAF==1) then  
  deallocate(CFSData1.MValue)  
  endif
  CFSData1.flagAF=0
   
  if(alive) then
    open(10,File=CFSData1.CFile,status='old')
      read(10,*)CFSData1.ncols,CFSData1.nrows   
      allocate(CFSData1.x(CFSData1.ncols))
      allocate(CFSData1.y(CFSData1.nrows))
      do i=1,CFSData1.ncols
        read(10,*) CFSData1.x(i)
		!if (CFSData1.x(i).gt.180) CFSData1.x(i)=CFSData1.x(i)-360
      end do
      do i=1,CFSData1.nrows
        read(10,*) CFSData1.y(i)
      end do 
     close(10)   
     CFSData1.flagAF=1 
  else
    Write(*,*) "No ",CFSData1.CFile
  return
  endif
  write(asRecord,"('Finish reading file of ',a)")trim(CFSData1.CFile)
  call Wrecord(asRecord)
end subroutine ReadCFScoord

!!********************************!!******************************************************* 
! 72 Write  hindcast result
!!********************************!!******************************************************* 
subroutine Writehincast(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(Scaldat) ::  caldat1
character(c100) ::   outputfile
integer  i,j,k,jday1,imem,JULDAY ,iyear 
character(c100) asRecord 
character(c1000) tempc 
 
  write(asRecord,*)'Begin to write hindcast'
  jday1=JULDAY(1,2,epp.BeginYear)  
  write(outputfile,'(a,i3.3,a1,i3.3,a4)')trim(epp.ReFile),epp.x,'\',epp.y,'.txt' 
  call Wrecord(outputfile)
  if (epp.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) 
    write (10,'(13(a),\)')'Year',char(9),'Month',char(9),'Day',char(9),&
	   'LeadT',char(9),'Obs',char(9),'Fcst',char(9),'EnsM'
    do i=1,epp.nems
	  write (10,'(2a,i3.3,\)')char(9),'Ens',i
    end do
    write (10,*) 

    do i=1,epp.ndays
      if (epp.Pens(i,1,1).ge.0.) then	 
        caldat1.julian=jday1+i-1
        call caldat(caldat1)   
        do j=1,epp.leadT
          write(tempc,'(i7)')	int(epp.ob(i,j)*100) 
          tempc=adjustl(tempc)             
          write (10,'(i4,a,i2,a,i2,a,i3,2a,\)')caldat1.iyyy,char(9),&
		   caldat1.mm,char(9),caldat1.id,char(9),j,char(9),trim(tempc) 
          
          write (tempc,'(200(a,i7))') char(9),int(epp.si(i,j)*100),&
            char(9),int(epp.Pens(i,j,epp.nems+1)*100),&
            (char(9),int(epp.Pens(i,j,imem)*100),imem=1,epp.nems)            
          call delblank(tempc)
          write (10,'(a)')trim(tempc)
      
          !write(tempc,'(i7)')	int(epp.si(i,j)*100) 
          !tempc=adjustl(tempc)             
          !write (10,'(2a,\)') char(9),trim(tempc) 
          !write(tempc,'(i7)')	int(epp.Pens(i,j,epp.nems+1)*100) 
          !tempc=adjustl(tempc)             
          !write (10,'(2a,\)') char(9),trim(tempc)  
          !do imem=1,epp.nems     
          !  write(tempc,'(i7)')	int(epp.Pens(i,j,imem)*100) 
          !  tempc=adjustl(tempc)  
          !  write (10,'(2a,\)') char(9),trim(tempc)  !epp.Pens(i,j,imem)            
          !end do
          !write (10,*) 
        end do !j
      endif
    end do    
  close(10)

  !output average 
   write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(epp.ReFile),epp.x,'\',epp.y,'m.txt' 
  if (epp.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) 
    write (10,'(7a,\)')'Year',char(9),'Month',char(9),'Day',char(9),'Observed'
	
    do i=1,epp.leadT,epp.nparint
	  write (10,'(2a,i3.3,\)')char(9),'Fcst',i
    end do  

    do i=1,epp.leadT,epp.nparint
	  write (10,'(2a,i3.3,\)')char(9),'Ensmean',i
    end do
    write (10,*)

    do i=365,epp.ndays 
	  if (epp.Pens(i,1,epp.nems+1).gt.0) then
        caldat1.julian=jday1+i-1
        call caldat(caldat1)  
        iyear= caldat1.iyyy
        write(tempc,'(i7)') int(epp.ob(i,1)*100)   
        tempc=adjustl(tempc)
        write (10,'(i4,a,i2,a,i2,2a,\)')caldat1.iyyy,char(9),caldat1.mm,char(9),caldat1.id,char(9),trim(tempc)
	  	   
        do j=1,epp.leadT,epp.nparint
          k=i-j+1
          if (k>0) then
	        if(epp.si(k,j).lt.0 ) k=k-1
            write(tempc,'(i7)') int(epp.si(k,j) *100)   
          else
            write(tempc,'(i7)') int(epp.si(j,j) *100)   
          endif
          tempc=adjustl(tempc) 
          write (10,'(2a,\)') char(9),trim(tempc) 
        end do !j
        do j=1,epp.leadT,epp.nparint
          k=i-j+1
          if (k>0) then
	        if(epp.Pens(k,j,epp.nems+1).lt.0 ) k=k-1
            write(tempc,'(i7)') int(epp.Pens(k,j,epp.nems+1)*100)   
          else
            write(tempc,'(i7)') int(epp.Pens(j,j,epp.nems+1)*100)   
          endif
          tempc=adjustl(tempc) 
          write (10,'(2a,\)') char(9),trim(tempc) 
        end do !j
        write(10,*)
      end if
    end do    
  close(10)


end subroutine 

!!********************************!!******************************************************* 
! 71 Schaake Shuffle
!!********************************!!******************************************************* 
subroutine SchaakeShuffle(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(SgenerateCaE) ::  gCanE
type(Sfcst2ensem) :: fe1 
type(Scaldat) ::  caldat1
type(Ssort) ::  sort1
integer i,j,k,iday,cday,iyear,jday,n
integer jday1,jday2,JULDAY,imem,m
real(r8) cmult, GSET,zero
real(r8),pointer ::   x(:)
integer ipp,ix,ISET 
    !epp2.iyear=iYear   
    !epp2.imonth=imonth  
    !epp2.iday=iday   

  jday1=JULDAY(epp.imonth,epp.iday,epp.iyear)
  jday2=JULDAY(1,1,epp.BeginYearEns)
  ISET=0
  allocate(x(epp.nems))
  !1.get the order of observed ensemble members
  allocate(sort1.a(epp.nems))
  sort1.n=epp.nems
  do i=1,365  !,epp.nparint   
    do j=1,epp.Events
	  do k=epp.BeginYearEns,epp.EndYearEns
        iday= JULDAY(1,1,k)+i-jday2	 	!闰年向后推一天，跳过2.29日			  
	    !if(((mod(k,4)==0.and.mod(k,100)/=0).or. mod(k,400)==0).and.&
	    !  i.gt.59 ) then
		!    iday=iday+1  !leap year 闰年向后推一天 跳过2月29日 后处理中没有2月29日          
        !end if 
        sort1.a(k-epp.BeginYearEns+1)=epp.Eep(iday,j)          
      end do
      sort1.indx=>epp.indEep(i,j,:)  ![365,Events,members]
      call indexx(sort1)	   
    end do
  end do
  deallocate(sort1.a)
  !2.get the order of coefficient of correlation between z(fcst) and z(obs) 
  sort1.n=epp.Events
  do i=1,365   
    sort1.a=>epp.rho_est(i,:)
    sort1.indx=>epp.indrho(i,:)
    call indexx(sort1)    
    
  end do
  !3.historical observed data to construct   ensemble members initial value
 
  do i=1,epp.ndays   
	caldat1.julian=jday1+i-1
    call caldat(caldat1)    
    !cday=min(365,caldat1.julian-JULDAY(1,1,caldat1.iyyy)+1) !(1-365)
    cday=caldat1.julian-JULDAY(1,1,caldat1.iyyy)+1 !(1-366)    
    do j=1,epp.leadT
	  do k=epp.BeginYearEns,epp.EndYearEns  
        iday= JULDAY(1,1,k)+cday-jday2	
	     
        epp.Pens(i,j,k-epp.BeginYearEns+1)=epp.ep(iday,j)  
        !epp.Pens(i,j,k-epp.BeginYearEns+1)=epp.si(i,j)  !forecast for divided compsite events.
        
		!if (i.eq.4 .and. j.eq.5.and.(k-epp.BeginYearEns+1).eq.9)   then
		!  n=n
		!end if       
      end do
    end do
  end do
  !4.Canonical Events to daily data  Eens=>Pens
  !real(r8),pointer ::   Eens(:,:,:) !Events Ensemble [date,Events,members]
  !real(r8),pointer ::   Pens(:,:,:) !Ensemble [date,lead time,members]
  sort1.n=epp.nems
  do i=1,epp.ndays 
    if (epp.Eens(i,1,1).lt.0) then
      epp.Pens(i,:,:)=-99.0
	else 
      caldat1.julian=jday1+i-1
      call caldat(caldat1)    
      cday=min(365,caldat1.julian-JULDAY(1,1,caldat1.iyyy)+1) !(1-365)  
     

      do j=1,epp.Events  
          
        k=epp.indrho(cday,j)  ! 从小到大排序，优先计算相关系数小的
        
        k=epp.Events -j+1     !从最后一个事件开始算起 此处修改
        do imem=1,epp.nems  
          x(imem)=0.d0
          do m=epp.Estart(k),epp.Estop(k)  
            if (m.lt.1) then
              n=i+m
              if (n.lt.1) n=1
              if (m.ne.0) x(imem)=x(imem)+epp.Pens(n,1,imem)  
            else
	          x(imem)=x(imem)+epp.Pens(i,m,imem)
            end if
          end do !m
          if (epp.Estart(k).gt.0) then
            x(imem)=x(imem)/(epp.Estop(k)-epp.Estart(k)+1)
          else
            x(imem)=x(imem)/(epp.Estop(k)-epp.Estart(k))   
          end if
          
        end do !imem
        sort1.a=>x
        sort1.indx=>epp.indEep(cday,k,:)    ![365,Events,members]
        call indexx(sort1)  

        do imem=1,epp.nems  
          ipp=epp.indEep(cday,k,imem)	![365,Events,members]此处采用了观测rank
          iday= JULDAY(1,1,ipp+epp.BeginYearEns-1)+cday-jday2		
						         
          !if (epp.Eep(iday,k)>0.01) then 	   
          if (x(ipp) >0.01) then 	   
            cmult=epp.Eens(i,k,imem)/x(ipp) 
          else
            cmult=(epp.Estop(k)-epp.Estart(k)+1)
          endif
          do m=epp.Estart(k),epp.Estop(k)     ! k,k    !此处确定是否修改历史观测           
            if (m.gt.0) then 
              if (cmult<(epp.Estop(k)-epp.Estart(k)+1))	then	 
                epp.Pens(i,m,ipp)=epp.Pens(i,m,ipp)*cmult +zero(ix,gset,iset)
              else
                epp.Pens(i,m,ipp)=epp.Eens(i,k,imem)  +zero(ix,gset,iset)
              endif
            end if            
            
          end do !m
        end do !j      
      end do !imem
    endif
  end do
  !5.compute the average of ensemble members
  do i=1,epp.ndays 
    do j=1,epp.leadT
      epp.Pens(i,j,epp.nems+1)=0.d0
      do imem=1,epp.nems     
        epp.Pens(i,j,epp.nems+1)=epp.Pens(i,j,epp.nems+1)+epp.Pens(i,j,imem)    
      end do
      epp.Pens(i,j,epp.nems+1)=epp.Pens(i,j,epp.nems+1)/epp.nems
    end do !j
  end do
  deallocate(x)
end subroutine

!!********************************!!******************************************************* 
! 70 get a random data
!!********************************!!******************************************************* 
function zero(ix,gset,iset)
use PrecStru
implicit none 
real(r8) zero,zavg,zstd,gasdev, GSET
integer ix,ISET 
  if (ix.eq.0) ix = 129874633       
  zavg = .0001
  zstd = .00001
  zero = zavg +zstd*gasdev(ix,gset,iset)  
  return
end

!!********************************!!******************************************************* 
! 69  hindcast
!!********************************!!******************************************************* 
subroutine Hindcast(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
type(SgenerateCaE) ::  gCanE
type(Sfcst2ensem) :: fe1 
type(Scaldat) ::  caldat1 
integer i,j,k,iday,iyear,ii
integer jday1,JULDAY
character(c100) asRecord 
  write(asRecord,*)'Begin to Hindcast CFS'
  call Wrecord(asRecord)

  epp.nems=epp.EndYearEns-epp.BeginYearEns+1
  jday1=JULDAY(1,1,epp.BeginYear)
  fe1.ifcst=epp.ifcsttran
  fe1.iobs=epp.iobstran
  fe1.npp=epp.nems
  epp.Eens=-99.0
  !1.Generate Canonical Events 
  gCanE.leadT=epp.leadT
  gCanE.ndays=epp.ndays
  gCanE.Events=epp.Events
  gCanE.Estart=>epp.Estart
  gCanE.Estop=>epp.Estop
  gCanE.ob=>epp.ob
  gCanE.hob=>epp.ob
  gCanE.Eob=>epp.Eob  ![date, Events]
  call generateCaE(gCanE)
  gCanE.ob=>epp.si
  gCanE.Eob=>epp.Esi    
  call generateCaE(gCanE)
  gCanE.ndays=epp.ndays1
  gCanE.ob=>epp.ep
  gCanE.Eob=>epp.Eep    
  call generateCaE(gCanE)

  do iyear=epp.BeginYear,epp.EndYear  
    ii=365  !leap year 闰年计算366
    if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0).and.&
      epp.nparint.le.1 ) ii=366		  
    do j=1,ii,epp.nparint1 
      !do iday=1,epp.ndays   
      write(*,'(a,i5,a1,i5)')trim('year,day='),iyear,'/',j
      iday=JULDAY(1,1,iyear)+j-jday1
      k=min(j,365)
	  if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0).and.&
	    j.gt.59.and.epp.nparint.gt.1 ) then
	  	  iday=iday+1  !leap year 闰年向后推一天 跳过2月29日 后处理中没有2月29日
          k=k+1
      end if  
      k=min(j,365)  
      !caldat1.julian=jday1+iday-1
      !call caldat(caldat1)    
      !j=min(365,caldat1.julian-JULDAY(1,1,caldat1.iyyy)+1) !(1-365)
      do i=1,epp.Events
        !2.Generate Canonical Events ensembles
        fe1.pthresh_fcst=epp.pthresh_fcst(k,i)  ![365,Events] 
        fe1.fpop=epp.popfcst(k,i)
        fe1.fcavg= epp.cavgfcst(k,i)
        fe1.fccv= epp.ccvfcst(k,i)
        fe1.pthresh_obs=epp.pthresh_obs(k,i)
        fe1.obspop=epp.popobs(k,i)
        fe1.obscavg=epp.cavgobs(k,i)
        fe1.obsccv=epp.ccvobs(k,i)
        fe1.rho=epp.rho_est(k,i)
        fe1.fcst=epp.Esi(iday,i) ![date, Events]
        fe1.pp=>epp.Eens(iday,i,:)    ![date,Events,members]
		!if (iday.eq.4 .and. i.eq.13) then
		!  j=j
		!end if
        call fcst2ensem(fe1)
    end do !i Events
  end do   !j days 
  end do
  write(asRecord,*)'Begin to SchaakeShuffle CFS'
  call Wrecord(asRecord)

  call SchaakeShuffle(epp)
  write(asRecord,*)'Finished to SchaakeShuffle CFS'
  call Wrecord(asRecord)
end subroutine

!!********************************!!******************************************************* 
! 68 read  paramater information of system to hindcast
!!********************************!!******************************************************* 
subroutine ReadHpara(epp,flagC)
use PrecStru
implicit none 
type(SEPP) ::  epp
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
 
  open(unit=10,file=epp.CoFile)

    do 
      read(10,*,end=10)k,tempc,tvalue       
      if (k.eq.flagC) then
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.EpFile')) read(tvalue,*) epp.EpFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.ObFile')) read(tvalue,*) epp.ObFile
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.SiFile')) read(tvalue,*) epp.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.PaFile')) read(tvalue,*) epp.PaFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.leadT')) read(tvalue,*) epp.leadT
 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.DoFile')) read(tvalue,*) epp.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.ReFile')) read(tvalue,*) epp.ReFile

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.x1')) read(tvalue,*) epp.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.x2')) read(tvalue,*) epp.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.y1')) read(tvalue,*) epp.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.y2')) read(tvalue,*) epp.y2
      
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.BeginYear ')) read(tvalue,*) epp.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.EndYear')) read(tvalue,*) epp.EndYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.BeginYearEns ')) read(tvalue,*) epp.BeginYearEns 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.EndYearEns')) read(tvalue,*) epp.EndYearEns 

 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.iobstran ')) read(tvalue,*) epp.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.ifcsttran ')) read(tvalue,*) epp.ifcsttran 
      
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.Events')) then
          read(tvalue,*) epp.Events 
          allocate(epp.Estart(epp.Events))
          allocate(epp.Estop(epp.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.Estart'))read(10,*)k,(epp.Estart(i),i=1,epp.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('epp.Estop'))read(10,*)k,(epp.Estop(i),i=1,epp.Events)  
      endif
    end do 
  10 close(10) 
  epp.nems=epp.EndYearEns-epp.BeginYearEns+1
end subroutine

!!********************************!!******************************************************* 
!67  delete   space for EPP model to hindcast
!
!!********************************!!******************************************************* 
subroutine deletspaceEP(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
  !1.input data
  deallocate(epp.ep)
  deallocate(epp.ob)
  deallocate(epp.si)
  deallocate(epp.pthresh_obs) 
  deallocate(epp.pthresh_fcst)
  deallocate(epp.popobs)
  deallocate(epp.cavgobs)
  deallocate(epp.ccvobs)
  deallocate(epp.popfcst)
  deallocate(epp.cavgfcst)
  deallocate(epp.ccvfcst)
  deallocate(epp.rho_est) 
 
  !2.output data   
  deallocate(epp.Eep)
  deallocate(epp.Eob)
  deallocate(epp.Esi)
  deallocate(epp.Eens)
  deallocate(epp.Pens)
  deallocate(epp.indEep)
  deallocate(epp.indrho)
  deallocate(epp.PensArea)  

end subroutine
!!********************************!!******************************************************* 
!66  allocate space for EPP  model to hindcast
!
!!********************************!!******************************************************* 
subroutine NewspaceEP(epp)
use PrecStru
implicit none 
type(SEPP) ::  epp
integer julday,i,j
  !1.input data
     
  !i=JULDAY(12,31,epp.EndYearEns)-JULDAY(1,1,epp.BeginYearEns)+1
  !j=max(i,epp.ndays)
  allocate(epp.ep(epp.ndays+epp.leadT,epp.leadT))
  allocate(epp.ob(epp.ndays+epp.leadT,epp.leadT))
  allocate(epp.si(epp.ndays,epp.leadT))
  allocate(epp.pthresh_obs(365,epp.Events)) 
  allocate(epp.pthresh_fcst(365,epp.Events))
  allocate(epp.popobs(365,epp.Events))
  allocate(epp.cavgobs(365,epp.Events))
  allocate(epp.ccvobs(365,epp.Events))
  allocate(epp.popfcst(365,epp.Events))
  allocate(epp.cavgfcst(365,epp.Events))
  allocate(epp.ccvfcst(365,epp.Events))
  allocate(epp.rho_est(365,epp.Events)) 

  !2.output data
   
  allocate(epp.Eep(epp.ndays1,epp.Events))
  allocate(epp.Eob(epp.ndays,epp.Events))
  allocate(epp.Esi(epp.ndays,epp.Events))
  allocate(epp.Eens(epp.ndays,epp.Events,epp.nems+1))  !+1 to save average of ensembles
  allocate(epp.Pens(epp.ndays,epp.leadT,epp.nems+1))
  allocate(epp.indEep(365,epp.Events,epp.nems))
  allocate(epp.indrho(365,epp.Events))
  allocate(epp.PensArea(epp.ncols,epp.nrows,epp.leadT,epp.nems+1))
  !  real(r8),pointer ::   PensArea(:,:,:,:) !Ensemble [Lon,Lat,lead time,members]

end subroutine
!!********************************!!******************************************************* 
!65 Read EPP parameter
!!********************************!!******************************************************* 
subroutine ReadEpara(SrEPara1)
use PrecStru
implicit none 
type(SreadEppPara) ::  SrEPara1
character(c100) ::   tempc
integer  i,j  
real(r8),pointer :: a(:,:)
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read EPP parameter'
  call Wrecord(asRecord)
 
  allocate(a(365,SrEPara1.Events))
  open(unit=10,file=SrEPara1.PaFile)
    do 
      read(10,*,end=10)tempc    
	  do i=1,365   
        read(10,*)(a(i,j),j=1,SrEPara1.Events) !(,i=1,365)  
      end do
      if (trim(tempc).eq.trim('pthresh_obs'))   SrEPara1.pthresh_obs=a 
      if (trim(tempc).eq.trim('pthresh_fcst'))  SrEPara1.pthresh_fcst=a 
      if (trim(tempc).eq.trim('popobs'))   SrEPara1.popobs=a    
      if (trim(tempc).eq.trim('cavgobs'))  SrEPara1.cavgobs=a   
      if (trim(tempc).eq.trim('ccvobs'))   SrEPara1.ccvobs=a    
      if (trim(tempc).eq.trim('popfcst'))  SrEPara1.popfcst=a   
      if (trim(tempc).eq.trim('cavgfcst')) SrEPara1.cavgfcst=a  
      if (trim(tempc).eq.trim('ccvfcst'))  SrEPara1.ccvfcst=a   
      if (trim(tempc).eq.trim('rho_est'))  SrEPara1.rho_est=a   
  
    end do 
  10 close(10) 
  deallocate(a)
end subroutine

!!********************************!!******************************************************* 
!64 Indexes arrin(1:n) so that arrin(indx(i)) is in increasing order
!   for i=1,n.  Values of arrin are not changed.
!!********************************!!******************************************************* 
subroutine indexx(sort1)
use PrecStru
integer i,j,k,n,tn
real(r8) temp
real(r8),pointer :: a(:)
type(Ssort) ::  sort1
  n=sort1.n
  allocate(a(n))
  a=sort1.a
  sort1.indx=(/(i,i=1,n)/)  
  do i=1,n-1  
    k=i
    do j=i+1,n    
      if (a(k)>a(j))  k=j
    end do
    if (k.ne.i) then     
      temp=a(k)
      a(k)=a(i)
      a(i)=temp
      tn=sort1.indx(k)
      sort1.indx(k)=sort1.indx(i)
      sort1.indx(i)=tn
    end if
  end do
  deallocate(a)
end subroutine 


!!********************************!!******************************************************* 
! 63 read  paramater information of system
!!********************************!!******************************************************* 
subroutine fillall(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
integer i,n,j
real(r8),pointer ::   x(:) 
  allocate(x(365))
  n=CanE.nparint
  do i=1,CanE.Events

    x=1.0*CanE.nall(:,i)
    call fill (n,x,x)
    CanE.nall(:,i)=x

    call fill (n,CanE.pthresh_obs(:,i) ,CanE.pthresh_obs(:,i))
    call fill (n,CanE.pthresh_fcst(:,i) ,CanE.pthresh_fcst(:,i))
    call fill (n,CanE.avgobs(:,i) ,CanE.avgobs(:,i))
    call fill (n,CanE.popobs(:,i) ,CanE.popobs(:,i))
    call fill (n,CanE.cavgobs(:,i) ,CanE.cavgobs(:,i))
    call fill (n,CanE.ccvobs(:,i) ,CanE.ccvobs(:,i))
    call fill (n,CanE.avgfcst(:,i) ,CanE.avgfcst(:,i))
    call fill (n,CanE.popfcst(:,i) ,CanE.popfcst(:,i))
    call fill (n,CanE.cavgfcst(:,i) ,CanE.cavgfcst(:,i))
    call fill (n,CanE.ccvfcst(:,i) ,CanE.ccvfcst(:,i))
    call fill (n,CanE.rho_est(:,i) ,CanE.rho_est(:,i))
    call fill (n,CanE.rmsfcst(:,i) ,CanE.rmsfcst(:,i))
    call fill (n,CanE.effnsfcst(:,i) ,CanE.effnsfcst(:,i))
    call fill (n,CanE.ratio(:,i) ,CanE.ratio(:,i))
    call fill (n,CanE.eqts_fcst(:,i) ,CanE.eqts_fcst(:,i))
    call fill (n,CanE.cor(:,i) ,CanE.cor(:,i))
    call fill (n,CanE.ccor(:,i) ,CanE.ccor(:,i))
    call fill (n,CanE.trccor(:,i) ,CanE.trccor(:,i))
    call fill (n,CanE.rhoopt(:,i) ,CanE.rhoopt(:,i))

    x=1.0*CanE.nvalcal(:,i,1)
    call fill (n,x,x)
    CanE.nvalcal(:,i,1)=x
    x=1.0*CanE.nvalcal(:,i,2)
    call fill (n,x,x)
    CanE.nvalcal(:,i,2)=x
    x=1.0*CanE.nvalcal(:,i,3)
    call fill (n,x,x)
    CanE.nvalcal(:,i,3)=x
    x=1.0*CanE.nvalcal(:,i,4)
    call fill (n,x,x)
    CanE.nvalcal(:,i,4)=x

    if (CanE.verstats.eq.2) then
      do j=1,CanE.Events
        call fill (n,CanE.coro(:,i,j),CanE.coro(:,i,j))
        call fill (n,CanE.corf(:,i,j),CanE.corf(:,i,j))
        call fill (n,CanE.corX(:,i,j),CanE.corX(:,i,j)) 
        call fill (n,CanE.xm(:,i,j),CanE.xm(:,i,j)) 
      end do
    end if

    call fill (n,CanE.avgobscal(:,i,1) ,CanE.avgobscal(:,i,1))
    call fill (n,CanE.avgobscal(:,i,2) ,CanE.avgobscal(:,i,2))
    call fill (n,CanE.avgobscal(:,i,3) ,CanE.avgobscal(:,i,3))
    call fill (n,CanE.avgobscal(:,i,4) ,CanE.avgobscal(:,i,4))

    call fill (n,CanE.stdobscal(:,i,1) ,CanE.stdobscal(:,i,1))
    call fill (n,CanE.stdobscal(:,i,2) ,CanE.stdobscal(:,i,2))
    call fill (n,CanE.stdobscal(:,i,3) ,CanE.stdobscal(:,i,3))
    call fill (n,CanE.stdobscal(:,i,4) ,CanE.stdobscal(:,i,4))

    call fill (n,CanE.avgfcstcal(:,i,1) ,CanE.avgfcstcal(:,i,1))
    call fill (n,CanE.avgfcstcal(:,i,2) ,CanE.avgfcstcal(:,i,2))
    call fill (n,CanE.avgfcstcal(:,i,3) ,CanE.avgfcstcal(:,i,3))
    call fill (n,CanE.avgfcstcal(:,i,4) ,CanE.avgfcstcal(:,i,4))

    call fill (n,CanE.stdfcstcal(:,i,1) ,CanE.stdfcstcal(:,i,1))
    call fill (n,CanE.stdfcstcal(:,i,2) ,CanE.stdfcstcal(:,i,2))
    call fill (n,CanE.stdfcstcal(:,i,3) ,CanE.stdfcstcal(:,i,3))
    call fill (n,CanE.stdfcstcal(:,i,4) ,CanE.stdfcstcal(:,i,4))
    if (CanE.verstats.eq.1) then
      call fill (n,CanE.crps_avgobs(:,i) ,CanE.crps_avgobs(:,i))
      call fill (n,CanE.crps_clim(:,i) ,CanE.crps_clim(:,i))
      call fill (n,CanE.crps_fcst(:,i) ,CanE.crps_fcst(:,i))
      call fill (n,CanE.crps_ensavg(:,i) ,CanE.crps_ensavg(:,i))
      call fill (n,CanE.crps_ens(:,i) ,CanE.crps_ens(:,i))
      call fill (n,CanE.crps_ensmed(:,i) ,CanE.crps_ensmed(:,i))
      call fill (n,CanE.crps_cavgobs(:,i) ,CanE.crps_cavgobs(:,i))
      call fill (n,CanE.crps_cclim(:,i) ,CanE.crps_cclim(:,i))
      call fill (n,CanE.crps_cfcst(:,i) ,CanE.crps_cfcst(:,i))
      call fill (n,CanE.crps_censavg(:,i) ,CanE.crps_censavg(:,i))
      call fill (n,CanE.crps_cens(:,i) ,CanE.crps_cens(:,i))
      call fill (n,CanE.crps_censmed(:,i) ,CanE.crps_censmed(:,i))

      call fill (n,CanE.roc_area(:,i) ,CanE.roc_area(:,i))


      call fill (n,CanE.biasens(:,i) ,CanE.biasens(:,i))
      call fill (n,CanE.biasfcst(:,i) ,CanE.biasfcst(:,i))
      call fill (n,CanE.corfcst(:,i) ,CanE.corfcst(:,i))
      call fill (n,CanE.corens(:,i) ,CanE.corens(:,i))
      call fill (n,CanE.corensfcst(:,i) ,CanE.corensfcst(:,i))
      call fill (n,CanE.rmsens(:,i) ,CanE.rmsens(:,i))
      call fill (n,CanE.effnsens(:,i) ,CanE.effnsens(:,i))
      call fill (n,CanE.bss(:,i) ,CanE.bss(:,i))
      call fill (n,CanE.rmspop(:,i) ,CanE.rmspop(:,i))
      call fill (n,CanE.bsscpc1(:,i) ,CanE.bsscpc1(:,i))
      call fill (n,CanE.bsscpc2(:,i) ,CanE.bsscpc2(:,i))
      call fill (n,CanE.bsscpc3(:,i) ,CanE.bsscpc3(:,i))
      call fill (n,CanE.cdfdiffmax(:,i) ,CanE.cdfdiffmax(:,i))
      call fill (n,CanE.rmscdf(:,i) ,CanE.rmscdf(:,i))
    end if
  end do
  deallocate(x)
end subroutine
!!********************************!!******************************************************* 
! 62 read  paramater information of system to prepare parameters of EPP
!!********************************!!******************************************************* 
subroutine ReadCEpara(CanE,flagC)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   tempc,tvalue
integer  i,k,flagC
character(c100) asRecord 
 
  write(asRecord,*)'Begin to read parameter'
  call Wrecord(asRecord)
   
  open(unit=10,file=CanE.CoFile)
    do 
      read(10,*,end=10)k,tempc,tvalue     
      if (k.eq.flagC) then
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ObFile')) read(tvalue,*) CanE.ObFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.flagPo')) read(tvalue,*) CanE.flagPo

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.SiFile')) read(tvalue,*) CanE.SiFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ReFile')) read(tvalue,*) CanE.ReFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.leadT')) read(tvalue,*) CanE.leadT
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x1')) read(tvalue,*) CanE.x1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.x2')) read(tvalue,*) CanE.x2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y1')) read(tvalue,*) CanE.y1
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.y2')) read(tvalue,*) CanE.y2
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minwin')) read(tvalue,*) CanE.minwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.maxwin')) read(tvalue,*) CanE.maxwin
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.DoFile')) read(tvalue,*) CanE.DoFile
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.BeginYear ')) read(tvalue,*) CanE.BeginYear 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.EndYear')) read(tvalue,*) CanE.EndYear 

        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minall')) read(tvalue,*) CanE.minall 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.minpos ')) read(tvalue,*) CanE.minpos 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.pop_fraction ')) read(tvalue,*) CanE.pop_fraction 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iobstran ')) read(tvalue,*) CanE.iobstran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.ifcsttran ')) read(tvalue,*) CanE.ifcsttran 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.iopt_rho ')) read(tvalue,*) CanE.iopt_rho 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cor_weight ')) read(tvalue,*) CanE.cor_weight 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgobsx ')) read(tvalue,*) CanE.cavgobsx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.cavgfcstx ')) read(tvalue,*) CanE.cavgfcstx 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.verstats ')) read(tvalue,*) CanE.verstats 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nparint ')) read(tvalue,*) CanE.nparint 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.nems ')) read(tvalue,*) CanE.nems 
    
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Events')) then
          read(tvalue,*) CanE.Events 
          allocate(CanE.Estart(CanE.Events))
          allocate(CanE.Estop(CanE.Events))
        end if 
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estart'))read(10,*)k,(CanE.Estart(i),i=1,CanE.Events)  
        if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE.Estop'))read(10,*)k,(CanE.Estop(i),i=1,CanE.Events)  
      endif
    end do 
  10 close(10) 
end subroutine

!!********************************!!******************************************************* 
! 61 Write  paramater information of EPP
!!********************************!!******************************************************* 
subroutine WriteCEpara(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   outputfile
character(c1000) ::tempc
integer  i,j,k  
character(c100) asRecord 
 
  write(asRecord,*)'Begin to write parameter'
  call Wrecord(asRecord)
   
 
  write(outputfile,'(a,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),CanE.x,'\',CanE.y,'.txt' 
  if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) 
	  write (10,'(a,TR1,i3,TR1,i4)') trim('nall'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,i4))') (char(9),CanE.nall(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,i3))') (CanE.nall(i,j),j=1,CanE.Events)   
      end do 

 	  write (10,'(a,TR1,i3,TR1,i3)') trim('pthresh_obs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.pthresh_obs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.pthresh_obs(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('pthresh_fcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.pthresh_fcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.pthresh_fcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('avgobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.avgobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.avgobs(i,j),j=1,CanE.Events)   
      end do 
 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('popobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.popobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.popobs(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cavgobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cavgobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.avgobs(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccvobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccvobs(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.ccvobs(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('avgfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.avgfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.avgfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('popfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.popfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.popfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cavgfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cavgfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
       !write (10,'(20(TR1,f5.2))') (CanE.cavgfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccvfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccvfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.ccvfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rho_est'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rho_est(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.rho_est(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rmsfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rmsfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.rmsfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('effnsfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.effnsfcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.effnsfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ratio'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ratio(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.ratio(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('eqts_fcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.eqts_fcst(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.eqts_fcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cor'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cor(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.cor(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('ccor'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.ccor(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.ccor(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('trccor'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.trccor(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.trccor(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rhoopt'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rhoopt(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.rhoopt(i,j),j=1,CanE.Events)   
      end do 

      do k=1,4
	    write (10,'(a,i1,TR1,i3,TR1,i3)') trim('nvalcal'),k,365,CanE.Events   
        do i=1,365   
          write (tempc,'(365(a,i5))') (char(9),CanE.nvalcal(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)
          !write (10,'(20(TR1,i4))') (CanE.nvalcal(i,j,k),j=1,CanE.Events)   
        end do 
      end do

      do k=1,4
	    write (10,'(a,i1,TR1,i3,TR1,i3)') trim('avgobscal'),k,365,CanE.Events   
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.avgobscal(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)
          !write (10,'(20(TR1,f5.2))') (CanE.avgobscal(i,j,k),j=1,CanE.Events)   
        end do 
      end do 

      do k=1,4
	    write (10,'(a,i1,TR1,i3,TR1,i3)') trim('avgfcstcal'),k,365,CanE.Events   
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.avgfcstcal(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)
          !write (10,'(20(TR1,f5.2))') (CanE.avgfcstcal(i,j,k),j=1,CanE.Events)   
        end do 
      end do 

      do k=1,4
	    write (10,'(a,i1,TR1,i3,TR1,i3)') trim('avgfcstcal'),k,365,CanE.Events   
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.avgfcstcal(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)
          !write (10,'(20(TR1,f5.2))') (CanE.avgfcstcal(i,j,k),j=1,CanE.Events)   
        end do 
      end do 

      do k=1,4
	    write (10,'(a,i1,TR1,i3,TR1,i3)') trim('stdfcstcal'),k,365,CanE.Events   
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.stdfcstcal(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)
          !write (10,'(20(TR1,f5.2))') (CanE.stdfcstcal(i,j,k),j=1,CanE.Events)   
        end do 
      end do 

    if (CanE.verstats.eq.1) then

 	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_avgobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_avgobs(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_avgobs(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_clim'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_clim(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_clim(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_fcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_fcst(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_fcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_ensavg'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_ensavg(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_ensavg(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_ens'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_ens(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_ens(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_ensmed'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_ensmed(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_ensmed(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_cavgobs'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_cavgobs(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_cavgobs(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_cclim'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_cclim(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_cclim(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_cfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_cfcst(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_cfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_censavg'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_censavg(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_censavg(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_cens'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_cens(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_cens(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('crps_censmed'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.crps_censmed(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.crps_censmed(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('biasens'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.biasens(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.biasens(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('biasfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.biasfcst(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.biasfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('corens'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.corens(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.corens(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('corensfcst'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.corensfcst(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.corensfcst(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rmsens'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rmsens(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.rmsens(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('effnsens'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.effnsens(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
       !write (10,'(20(TR1,f5.2))') (CanE.effnsens(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('bss'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.bss(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.bss(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rmspop'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rmspop(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.rmspop(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('bsscpc1'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.bsscpc1(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.bsscpc1(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('bsscpc2'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.bsscpc2(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.bsscpc2(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('bsscpc3'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.bsscpc3(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.bsscpc3(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('cdfdiffmax'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.cdfdiffmax(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.cdfdiffmax(i,j),j=1,CanE.Events)   
      end do 

	  write (10,'(a,TR1,i3,TR1,i3)') trim('rmscdf'),365,CanE.Events   
      do i=1,365   
        write (tempc,'(365(a,f8.2))') (char(9),CanE.rmscdf(i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        !write (10,'(20(TR1,f5.2))') (CanE.rmscdf(i,j),j=1,CanE.Events)   
      end do 

    end if

    if (CanE.verstats.eq.2) then

      do k=1,CanE.Events
	    write (tempc,'(a,i3,a,i3,a,i3)') trim('corf'),k,char(9),365,char(9),CanE.Events   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.corf(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)            
        end do 
      end do
      do k=1,CanE.Events
	    write (tempc,'(a,i3,a,i3,a,i3)') trim('coro'),k,char(9),365,char(9),CanE.Events   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.coro(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)            
        end do 
      end do


      do k=1,CanE.Events
	    write (tempc,'(a,i3,a,i3,a,i3)') trim('corX'),k,char(9),365,char(9),CanE.Events   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
	    
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.corX(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)             
        end do 
      end do
      do k=1,CanE.Events
	    write (tempc,'(a,i3,a,i3,a,i3)') trim('xm'),k,char(9),365,char(9),CanE.Events   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
	   
        do i=1,365   
          write (tempc,'(365(a,f8.2))') (char(9),CanE.xm(i,j,k),j=1,CanE.Events)   
          call delblank(tempc)
          write (10,'(a)')trim(tempc)             
        end do 
      end do 


	  write (10,'(a,TR1,i3,TR1,i3)') trim('Meancoro'),CanE.Events,CanE.Events   
      do i=1,CanE.Events
        write (tempc,'(365(a,f8.2))') (char(9),CanE.coro(366,i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)           
      end do 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('Meancorf'),CanE.Events,CanE.Events   
      do i=1,CanE.Events
        write (tempc,'(365(a,f8.2))') (char(9),CanE.corf(366,i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)           
      end do 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('MeancorX'),CanE.Events,CanE.Events   
      do i=1,CanE.Events
        write (tempc,'(365(a,f8.2))') (char(9),CanE.corX(366,i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)           
      end do 
	  write (10,'(a,TR1,i3,TR1,i3)') trim('Meanxm'),CanE.Events,CanE.Events   
      do i=1,CanE.Events
        write (tempc,'(365(a,f8.2))') (char(9),CanE.xm(366,i,j),j=1,CanE.Events)   
        call delblank(tempc)
        write (10,'(a)')trim(tempc)           
      end do 


    end if


  close(10)

end subroutine 
 

!!********************************!!******************************************************* 
!60  delete   space for EPP model to prepare parameter
!
!!********************************!!******************************************************* 
subroutine deletspaceCE(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE

  deallocate(CanE.ob)
  deallocate(CanE.si)
  deallocate(CanE.Eob)
  deallocate(CanE.Esi)
  deallocate(CanE.model)
  deallocate(CanE.basin)
 
  deallocate(CanE.nall)
  deallocate(CanE.pthresh_obs)
  deallocate(CanE.pthresh_fcst)
  deallocate(CanE.avgobs)
  deallocate(CanE.popobs)
  deallocate(CanE.cavgobs)
  deallocate(CanE.ccvobs)
  deallocate(CanE.avgfcst)
  deallocate(CanE.popfcst)
  deallocate(CanE.cavgfcst)
  deallocate(CanE.ccvfcst)
  deallocate(CanE.rho_est)
  deallocate(CanE.rmsfcst)
  deallocate(CanE.effnsfcst)
  deallocate(CanE.ratio)
  deallocate(CanE.eqts_fcst)
  deallocate(CanE.cor)
  deallocate(CanE.ccor)
  deallocate(CanE.trccor)
  deallocate(CanE.rhoopt)
  deallocate(CanE.nvalcal)
  deallocate(CanE.avgobscal)
  deallocate(CanE.stdobscal)
  deallocate(CanE.avgfcstcal)
  deallocate(CanE.stdfcstcal)

  deallocate(CanE.crps_avgobs)
  deallocate(CanE.crps_clim)
  deallocate(CanE.crps_fcst)
  deallocate(CanE.crps_ensavg)
  deallocate(CanE.crps_ens)
  deallocate(CanE.crps_ensmed)
  deallocate(CanE.crps_cavgobs)
  deallocate(CanE.crps_cclim)
  deallocate(CanE.crps_cfcst)
  deallocate(CanE.crps_censavg)
  deallocate(CanE.crps_cens)
  deallocate(CanE.crps_censmed)
  deallocate(CanE.roc_area)

  deallocate(CanE.biasfcst)
  deallocate(CanE.biasens)
  deallocate(CanE.corfcst)
  deallocate(CanE.corens)
  deallocate(CanE.corensfcst)
  deallocate(CanE.rmsens)
  deallocate(CanE.effnsens)
  deallocate(CanE.bss)
  deallocate(CanE.rmspop)
  deallocate(CanE.bsscpc1)
  deallocate(CanE.bsscpc2)
  deallocate(CanE.bsscpc3)
  deallocate(CanE.cdfdiffmax)
  deallocate(CanE.rmscdf)
  if (CanE.verstats.eq.2) then
    deallocate(CanE.coro)
    deallocate(CanE.corf)
    deallocate(CanE.corX)
    deallocate(CanE.xm)
  endif
end subroutine
!!********************************!!******************************************************* 
!59  allocate space for EPP  model to prepare parameter
!
!!********************************!!******************************************************* 
subroutine NewspaceCE(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
  !1.input data
  allocate(CanE.ob(CanE.ndays+CanE.leadT,CanE.leadT))
  allocate(CanE.si(CanE.ndays,CanE.leadT))
  allocate(CanE.model(CanE.x2))
  allocate(CanE.basin(CanE.y2))
 
  allocate(CanE.Eob(CanE.ndays,CanE.Events))
  allocate(CanE.Esi(CanE.ndays,CanE.Events))
  CanE.Eob=0
  CanE.Esi=0
  !2.output data
   
  allocate(CanE.nall(365,CanE.Events))
  allocate(CanE.pthresh_obs(365,CanE.Events))
  allocate(CanE.pthresh_fcst(365,CanE.Events))
  allocate(CanE.avgobs(365,CanE.Events))
  allocate(CanE.popobs(365,CanE.Events))
  allocate(CanE.cavgobs(365,CanE.Events))
  allocate(CanE.ccvobs(365,CanE.Events))
  allocate(CanE.avgfcst(365,CanE.Events))
  allocate(CanE.popfcst(365,CanE.Events))
  allocate(CanE.cavgfcst(365,CanE.Events))
  allocate(CanE.ccvfcst(365,CanE.Events))
  allocate(CanE.rho_est(365,CanE.Events))
  allocate(CanE.rmsfcst(365,CanE.Events))
  allocate(CanE.effnsfcst(365,CanE.Events))
  allocate(CanE.ratio(365,CanE.Events))
  allocate(CanE.eqts_fcst(365,CanE.Events))
  allocate(CanE.cor(365,CanE.Events))
  allocate(CanE.ccor(365,CanE.Events))
  allocate(CanE.trccor(365,CanE.Events))
  allocate(CanE.rhoopt(365,CanE.Events))
  if (CanE.verstats.eq.2) then
    allocate(CanE.coro(366,CanE.Events,CanE.Events))
    allocate(CanE.corf(366,CanE.Events,CanE.Events))
    allocate(CanE.corX(366,CanE.Events,CanE.Events))
    allocate(CanE.xm(366,CanE.Events,CanE.Events))
  end if
  allocate(CanE.nvalcal(365,CanE.Events,4))
  allocate(CanE.avgobscal(365,CanE.Events,4))
  allocate(CanE.stdobscal(365,CanE.Events,4))
  allocate(CanE.avgfcstcal(365,CanE.Events,4))
  allocate(CanE.stdfcstcal(365,CanE.Events,4))

  allocate(CanE.crps_avgobs(365,CanE.Events))
  allocate(CanE.crps_clim(365,CanE.Events))
  allocate(CanE.crps_fcst(365,CanE.Events))
  allocate(CanE.crps_ensavg(365,CanE.Events))
  allocate(CanE.crps_ens(365,CanE.Events))
  allocate(CanE.crps_ensmed(365,CanE.Events))
  allocate(CanE.crps_cavgobs(365,CanE.Events))
  allocate(CanE.crps_cclim(365,CanE.Events))
  allocate(CanE.crps_cfcst(365,CanE.Events))
  allocate(CanE.crps_censavg(365,CanE.Events))
  allocate(CanE.crps_cens(365,CanE.Events))
  allocate(CanE.crps_censmed(365,CanE.Events))
  allocate(CanE.roc_area(365,CanE.Events))

  allocate(CanE.biasfcst(365,CanE.Events))
  allocate(CanE.biasens(365,CanE.Events))
  allocate(CanE.corfcst(365,CanE.Events))
  allocate(CanE.corens(365,CanE.Events))
  allocate(CanE.corensfcst(365,CanE.Events))
  allocate(CanE.rmsens(365,CanE.Events))
  allocate(CanE.effnsens(365,CanE.Events))
  allocate(CanE.bss(365,CanE.Events))
  allocate(CanE.rmspop(365,CanE.Events))
  allocate(CanE.bsscpc1(365,CanE.Events))
  allocate(CanE.bsscpc2(365,CanE.Events))
  allocate(CanE.bsscpc3(365,CanE.Events))
  allocate(CanE.cdfdiffmax(365,CanE.Events))
  allocate(CanE.rmscdf(365,CanE.Events))



end subroutine

!!********************************!!******************************************************* 
!58   from historical observed data and forecast to get parameter of EPP model
!
! Input arguments:
! 
! Output arguments:      
!
!!********************************!!******************************************************* 
subroutine EppPara(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
type(Swindow) ::  win1
type(Sensverify) :: ensv1
type(Sfcst2ensem) :: fe1  
type(Sfcstparm) :: fp1 
type(Scrpseval) ::  crpseval1
type(Sroc) ::  roc1 
type(SgenerateCaE) ::  gCanE 
type(CRelEffP) ::  CRE
character(c100) asRecord 
integer i,j,k
  write(asRecord,*)'Begin to compute parameter'
  call Wrecord(asRecord)
  win1.BeginYear=CanE.BeginYear
  win1.EndYear=CanE.EndYear
  win1.maxwin=CanE.maxwin
  win1.minwin=CanE.minwin
  win1.minall=CanE.minall
  win1.minpos=CanE.minpos
  win1.pop_fraction=CanE.pop_fraction  
  win1.ndays=CanE.ndays  
  k=win1.maxwin*(win1.EndYear-win1.BeginYear+1) 
  allocate(win1.obs(k))
  allocate(win1.fcst(k))
 

  if (CanE.verstats.eq.2) then
    allocate(win1.Enall(365,CanE.Events))
    allocate(win1.Eobs(365,CanE.Events,k))
    allocate(win1.Efcst(365,CanE.Events,k))
    allocate(CRE.yy(k))
    allocate(CRE.yc(k))
  end if

  if (CanE.verstats.eq.1) then
    fe1.ifcst=CanE.ifcsttran
    fe1.iobs=CanE.iobstran
    fe1.npp=CanE.nems
    allocate(crpseval1.ensem(fe1.npp+1,k))
    crpseval1.nmem=fe1.npp

    roc1.nmem=fe1.npp
    roc1.nthresh=9
    roc1.xx=>crpseval1.ensem
    roc1.obs=>win1.obs
    allocate(roc1.hr(roc1.nthresh))
    allocate(roc1.far(roc1.nthresh))
    allocate(roc1.obsthresh(roc1.nthresh))

    ensv1.iobstran=CanE.iobstran
    ensv1.ifcsttran=CanE.ifcsttran
    ensv1.nmem=fe1.npp
    ensv1.obs=>win1.obs
    ensv1.fcst=>win1.fcst
    ensv1.xppx=>crpseval1.ensem
    allocate(ensv1.icatobs(k))
    allocate(ensv1.avgpx(k))
    allocate(ensv1.stdpx(k))
    allocate(ensv1.cavgpx(k))
    allocate(ensv1.cstdpx(k))
    allocate(ensv1.ccvpx(k))
    allocate(ensv1.poppx(k))
    allocate(ensv1.cumestprob(k))
    allocate(ensv1.xobs(k))
    allocate(ensv1.cdfobs(k))
    allocate(ensv1.prhit33(3,k))
  endif 

  fp1.obsx=>win1.obs
  fp1.fcstx=>win1.fcst    
  

  fp1.iobstran=CanE.iobstran  
  fp1.ifcsttran=CanE.ifcsttran
  fp1.iopt_rho=CanE.iopt_rho
  fp1.cor_weight=CanE.cor_weight


  !1.Generate Canonical Events 

  gCanE.leadT=CanE.leadT
  gCanE.ndays=CanE.ndays
  gCanE.Events=CanE.Events
  gCanE.ob=>CanE.ob
  gCanE.hob=>CanE.ob
  gCanE.Estart=>CanE.Estart
  gCanE.Estop=>CanE.Estop
  gCanE.Eob=>CanE.Eob

  call generateCaE(gCanE)
  gCanE.ob=>CanE.si
  gCanE.Eob=>CanE.Esi    
  gCanE.ndays=CanE.ndays
 
  call generateCaE(gCanE)

  win1.Eob=>CanE.Eob
  win1.Esi=>CanE.Esi
  do i=1,CanE.Events
    do j=1,365,CanE.nparint
	  !write(*,'(a,i4,tr1,a,i4)')trim('Events='),i,trim('day='),j
      !if (j.eq.226) &
	  ! win1.ne=win1.ne
      !2.Compute the window of data
      win1.ne=i
      win1.day=j
      call  window(win1) 
	   
      if (CanE.verstats.eq.2) then
        win1.Enall(j,i)=win1.nall
        win1.Eobs(j,i,:)=win1.obs
        win1.Efcst(j,i,:)=win1.fcst
      end if    
	      
      !3.Estimate the parameter of EPP                           
      fp1.nobs=win1.nall
      fp1.pthresh1=win1.pthresh_obs
      fp1.pthresh2=win1.pthresh_fcst

      fp1.cavgobsx=CanE.cavgobsx
      fp1.cavgfcstx=CanE.cavgobsx

      call epp_fcstparms(fp1) 

      CanE.cavgobsx=fp1.cavgobsx
      CanE.cavgfcstx=fp1.cavgfcstx
      CanE.nall(j,i)=win1.nall
      CanE.pthresh_obs(j,i)=win1.pthresh_obs
      CanE.pthresh_fcst(j,i)=win1.pthresh_fcst
      CanE.avgobs(j,i)=fp1.avgobs
      CanE.popobs(j,i)=fp1.popobs
      CanE.cavgobs(j,i)=fp1.cavgobs
      CanE.ccvobs(j,i)=fp1.ccvobs
      CanE.avgfcst(j,i)=fp1.avgfcst
      CanE.popfcst(j,i)=fp1.popfcst
      CanE.cavgfcst(j,i)=fp1.cavgfcst
      CanE.ccvfcst(j,i)=fp1.ccvfcst
      CanE.rho_est(j,i)=fp1.rho_est
      CanE.rmsfcst(j,i)=fp1.rmsfcst
      CanE.effnsfcst(j,i)=max(0.,fp1.effnsfcst)
      CanE.ratio(j,i)=fp1.ratio
      CanE.eqts_fcst(j,i)=fp1.eqts_fcst
      CanE.cor(j,i)=fp1.cor
      CanE.ccor(j,i)=fp1.ccor
      CanE.trccor(j,i)=fp1.trccor
      CanE.rhoopt(j,i)=fp1.rhoopt

      CanE.nvalcal(j,i,1)=fp1.nvalcal(1,1)
      CanE.nvalcal(j,i,2)=fp1.nvalcal(1,2)
      CanE.nvalcal(j,i,3)=fp1.nvalcal(2,1)
      CanE.nvalcal(j,i,4)=fp1.nvalcal(2,2)

      CanE.avgobscal(j,i,1)=fp1.avgobscal(1,1)
      CanE.avgobscal(j,i,2)=fp1.avgobscal(1,2)
      CanE.avgobscal(j,i,3)=fp1.avgobscal(2,1)
      CanE.avgobscal(j,i,4)=fp1.avgobscal(2,2)

      CanE.stdobscal(j,i,1)=fp1.stdobscal(1,1)
      CanE.stdobscal(j,i,2)=fp1.stdobscal(1,2)
      CanE.stdobscal(j,i,3)=fp1.stdobscal(2,1)
      CanE.stdobscal(j,i,4)=fp1.stdobscal(2,2)

      CanE.avgfcstcal(j,i,1)=fp1.avgfcstcal(1,1)
      CanE.avgfcstcal(j,i,2)=fp1.avgfcstcal(1,2)
      CanE.avgfcstcal(j,i,3)=fp1.avgfcstcal(2,1)
      CanE.avgfcstcal(j,i,4)=fp1.avgfcstcal(2,2)

      CanE.stdfcstcal(j,i,1)=fp1.stdfcstcal(1,1)
      CanE.stdfcstcal(j,i,2)=fp1.stdfcstcal(1,2)
      CanE.stdfcstcal(j,i,3)=fp1.stdfcstcal(2,1)
      CanE.stdfcstcal(j,i,4)=fp1.stdfcstcal(2,2)
   
      if (CanE.verstats.eq.1) then
       
        fe1.pthresh_fcst=CanE.pthresh_fcst(j,i)

        fe1.fpop=CanE.popfcst(j,i)
        fe1.fcavg= CanE.cavgfcst(j,i)
        fe1.fccv= CanE.ccvfcst(j,i)

        fe1.pthresh_obs=CanE.pthresh_obs(j,i)
        fe1.obspop=CanE.popobs(j,i)
        fe1.obscavg=CanE.cavgobs(j,i)
        fe1.obsccv=CanE.ccvobs(j,i)
        fe1.rho=CanE.rho_est(j,i)
		
        do k=1,CanE.nall(j,i)
          fe1.fcst=fp1.fcstx(k)
          fe1.pp=>crpseval1.ensem(:,k) 
          call fcst2ensem(fe1)
		end do
      
        crpseval1.fcst=>win1.fcst 
        crpseval1.obs=>win1.obs 
        crpseval1.nobs=CanE.nall(j,i) 
        crpseval1.pthresh=CanE.pthresh_obs(j,i) 

        call  crps_eval(crpseval1)
        CanE.crps_avgobs(j,i)=crpseval1.crps_avgobs         
        CanE.crps_clim(j,i)=crpseval1.crps_clim
        CanE.crps_fcst(j,i)=crpseval1.crps_fcst
        CanE.crps_ensavg(j,i)=crpseval1.crps_ensavg
        CanE.crps_ens(j,i)=crpseval1.crps_ens
        CanE.crps_ensmed(j,i)=crpseval1.crps_ensmed
        CanE.crps_cavgobs(j,i)=crpseval1.crps_cavgobs
        CanE.crps_cclim(j,i)=crpseval1.crps_cclim
        CanE.crps_cfcst(j,i)=crpseval1.crps_cfcst
        CanE.crps_censavg(j,i)=crpseval1.crps_censavg
        CanE.crps_cens(j,i)=crpseval1.crps_cens
        CanE.crps_censmed(j,i)=crpseval1.crps_censmed

   
        roc1.nfcst=CanE.nall(j,i) 
        call roc(roc1) 
        CanE.roc_area(j,i)=roc1.roc_score

        ensv1.nobs=CanE.nall(j,i) 
        ensv1.pthresh=CanE.pthresh_obs(j,i) 
        call ensverify(ensv1)

 
        CanE.biasfcst(j,i)=CanE.avgfcst(j,i)/CanE.avgobs(j,i)
        CanE.biasens(j,i)=ensv1.avgens/CanE.avgobs(j,i)
        !CanE.corfcst(j,i)=

        CanE.corens(j,i)=ensv1.corens
        CanE.corensfcst(j,i)=ensv1.corensfcst
        CanE.rmsens(j,i)=ensv1.rmsens
        CanE.effnsens(j,i)=max(0.,1.-(ensv1.rmsens*ensv1.rmsens/ensv1.stdobs/ensv1.stdobs))
        CanE.bss(j,i)=ensv1.bss
        CanE.rmspop(j,i)=ensv1.rmspop
        CanE.bsscpc1(j,i)=ensv1.bsscpc(1)
        CanE.bsscpc2(j,i)=ensv1.bsscpc(2)
        CanE.bsscpc3(j,i)=ensv1.bsscpc(3)
        CanE.cdfdiffmax(j,i)=ensv1.cdfdiffmax
        CanE.rmscdf(j,i)=ensv1.rmscdf
    
 
      end if


    end do 
  end do
  !calculate the cross correlation among events

  if (CanE.verstats.eq.2) then
    do j=1,365,CanE.nparint
      do i=1,CanE.Events
	    write(*,'(a,i4,tr1,a,i4)')trim('Cross Coro: Events='),i,trim('day='),j
        CRE.yy=win1.Eobs(j,i,:)
        CRE.TimeSum=win1.Enall(j,i)
        do k=i,CanE.Events
          CRE.yc=win1.Eobs(j,k,:)
          call CRelEff(CRE)
          CanE.coro(j,i,k)=CRE.r 
        end do 
      end do 
    end do
   
    do j=1,365,CanE.nparint
      do i=1,CanE.Events
	    write(*,'(a,i4,tr1,a,i4)')trim('Cross Corf: Events='),i,trim('day='),j
        CRE.yy=win1.Efcst(j,i,:)
        CRE.TimeSum=win1.Enall(j,i)
        do k=i,CanE.Events
          CRE.yc=win1.Efcst(j,k,:)
          call CRelEff(CRE)
          CanE.corf(j,i,k)=CRE.r 
        end do 
      end do 
    end do
    do j=1,365,CanE.nparint
      do i=1,CanE.Events 
        win1.Efcst(j,i,:)=win1.Efcst(j,i,:)-win1.Eobs(j,i,:) 
      end do 
    end do
    do j=1,365,CanE.nparint
      do i=1,CanE.Events
        write(*,'(a,i4,tr1,a,i4)')trim('Cross CorX: Events='),i,trim('day='),j
        CRE.yy=win1.Efcst(j,i,:) 
        CRE.TimeSum=win1.Enall(j,i)
        do k=i,CanE.Events
          CRE.yc=win1.Efcst(j,k,:) 
          call CRelEff(CRE)
          CanE.corX(j,i,k)=CRE.r 
          CanE.xm(j,i,k)=CRE.xm 
        end do 
      end do 
    end do

    do j=1,365,CanE.nparint
      do i=1,CanE.Events
        do k=i+1,CanE.Events
          CanE.coro(j,k,i)= CanE.coro(j,i,k) 
          CanE.corf(j,k,i)= CanE.corf(j,i,k) 
          CanE.corX(j,k,i)= CanE.corX(j,i,k)
          CanE.xm(j,k,i)= CanE.xm(j,i,k)
        end do 
      end do 
    end do
  end if

  if (CanE.nparint>1) call fillall(CanE)      

  if (CanE.verstats.eq.2) then
    ! compute the average of Cross correlation
    do i=1,CanE.Events
      do k=1,CanE.Events
        CanE.coro(366,k,i)= 0.
        CanE.corf(366,k,i)= 0.
        CanE.corX(366,k,i)= 0.
        CanE.xm(366,k,i)= 0.
        do j=1,365
          CanE.coro(366,k,i)= CanE.coro(366,k,i)+abs(CanE.coro(j,k,i)) 
          CanE.corf(366,k,i)= CanE.corf(366,k,i)+abs(CanE.corf(j,k,i)) 
          CanE.corX(366,k,i)= CanE.corX(366,k,i)+abs(CanE.corX(j,k,i))
          CanE.xm(366,k,i)= CanE.xm(366,k,i)+abs(CanE.xm(j,k,i))
        end do 
        CanE.coro(366,k,i)= CanE.coro(366,k,i)/365
        CanE.corf(366,k,i)= CanE.corf(366,k,i)/365
        CanE.corX(366,k,i)= CanE.corX(366,k,i)/365
        CanE.xm(366,k,i)= CanE.xm(366,k,i)/365
      end do 
    end do
  end if

  deallocate(win1.obs)   
  deallocate(win1.fcst)   
  if (CanE.verstats.eq.2) then
    deallocate(win1.Enall)
    deallocate(win1.Eobs)
    deallocate(win1.Efcst)
    deallocate(CRE.yy)
    deallocate(CRE.yc)
  end if
  if (CanE.verstats.eq.1) then

    deallocate(crpseval1.ensem)   
    deallocate(roc1.hr)   
    deallocate(roc1.far) 
	deallocate(ensv1.icatobs)
    deallocate(ensv1.avgpx)
    deallocate(ensv1.stdpx)
    deallocate(ensv1.cavgpx)
    deallocate(ensv1.cstdpx)
    deallocate(ensv1.ccvpx)
    deallocate(ensv1.poppx)
    deallocate(ensv1.cumestprob)
    deallocate(ensv1.xobs)
    deallocate(ensv1.cdfobs)
    deallocate(ensv1.prhit33 )
    deallocate(roc1.obsthresh)
  end if
  write(asRecord,*)'Finish to compute parameter'
  call Wrecord(asRecord)
 
end subroutine
!!********************************!!******************************************************* 
!57 Compute the window of data

!!********************************!!******************************************************* 

subroutine window(win1)
use PrecStru
implicit none 
type(Swindow) ::  win1
type(Sthreshold) ::  threshold1
integer i,j,k,kwin,iyear,jday1,cday
integer nall,npos,julday
  win1.win=win1.minwin
  win1.nall=0
  npos=0
  kwin=0
  jday1=JULDAY(1,1,win1.BeginYear)
  
  do while (win1.win.lt.win1.maxwin.and.&
           (win1.nall.lt.win1.minall.or.npos.lt.win1.minpos))            
    npos = 0  
    win1.win = win1.minwin + 5*kwin
    i=win1.ne  
	j=win1.day
    win1.nall=0
    do iyear=win1.BeginYear,win1.EndYear
      cday=julday(1,1,iyear)+j-jday1
      do k=max(1,cday-win1.win/2),min(win1.ndays,cday+win1.win/2)
        if(win1.Eob(k,i).ge.0.0.and.win1.Esi(k,i).ge.0.0) then
          win1.nall=win1.nall+1
          win1.obs(win1.nall)=win1.Eob(k,i)
          win1.fcst(win1.nall)=win1.Esi(k,i)
        end if
      end do
    end do
    threshold1.obs=>win1.obs
    threshold1.nobs=win1.nall
    threshold1.pop_fraction=win1.pop_fraction
    call threshold (threshold1)
    win1.pthresh_obs=threshold1.pthresh

    threshold1.obs=>win1.fcst  
    call threshold (threshold1)
    win1.pthresh_fcst=threshold1.pthresh

    do i=1,win1.nall
      !if (win1.obs(i).gt.win1.pthresh_obs) np = np + 1
      !if (win1.fcst(i).gt.win1.pthresh_fcst) nposfcst = nposfcst + 1
      if (win1.obs(i).gt.win1.pthresh_obs.and.win1.fcst(i).gt.win1.pthresh_fcst) npos = npos + 1               
    enddo
    kwin = kwin + 1    
  enddo

end subroutine


!!********************************!!******************************************************* 
!56  ensemble verify
!
! Input arguments:
!   
!   maxobs = maximum obs dimension
!   maxmem = maximum nmem dimension
!   nobs = number of observations/forecasts
!   obs(i) = vector of observed values, obs(i), i=1,nobs
!   nmem(i) = vector specifying number of ensemble members in each forecast,
!             nmem(i), i=1,nobs
!   xppx(i,j) = array of forecast member values, xppx(i,j), j=1,nmem(i), i=1,nobs
!   pthresh = threshold value to define conditional event 
!               (obs>pthresh or xppx>pthresh)
!
! Output arguments:      
!   avgpx(i) = ensemble average for forecast i (unconditional)
!   stdpx(i) = ensemble standard deviation of forecast i (unconditional)
!   cavgpx(i) = conditional average of forecast i (> pthresh)
!   cstdpx(i) = conditional stand deviation of forecast i (> pthresh)
!   ccvpx(i) = conditional coefficient of variation of forecast i  (> pthresh)
!                
!   poppx(i) = propability of precipitation of forecast i
 
!   avgobs = average value of obs(i)
!   stdobs = standard deviation of obs(i)
!   cvobs = coefficient of variation of obs(i)
!   cavgobs = conditional average value of obs(i) (obs>pthresh)
!   cstdobs = conditional stanedard deviation of obs(i) (obs>pthresh)
!   ccvobs = conditional coefficient of variation of obs(i) (obs>pthresh)
!   popobs = observed probability of precipitation

!   avgens = average of the ensemble forecast averages
!   stdens = standard deviation of the ensemble forecast averages
!   cvens = coefficient of variation of the ensemble forecast averages
!   popens = fraction of avgpx(i)>pthresh
!   cavgens = average of the conditional ensemble forecast averages
!   cstdens = standard deviation of the conditional ensemble forecast averages 
!            
!   ccvens = coefficient of variation of the conditional ensemble forecast averages 
!           
!   avgpopens = average of the ensemble forecast pop
!   stdpopens = standard deviation of ensemble forecast pop
!   avgcavgens = average of the ensemble conditional average
!   stdcavgens = standard deviation of the ensemble conditional average
!   avgccvens = average of the ensemble forecast ccv
!   stdccvens = standard deviation of the ensemble forecast ccv
!   biasens = (avgens - avgobs)/avgobs + 1
!   cbiasens = (cavgens - cavgobs)/cavgobs + 1
!   corens = coeffient of correlation between ensemble mean and 
!            observed value
!   ccorens = coefficient of correlation between conditional ensemble 
!             mean and observed value, given obs>pthresh
!   tccorens = coefficient of correlation between transformed 
!              conditional ensemble mean and transformed observed value,
!              given obs>pthresh
!   rhoopt = contingency estimate of bivariate normal correlation 
!            parameter for avgpx(i)
!   effens = Nash-Sutcliffe efficiency of ensemble mean forecast
!   ceffens = Nash-Sutcliffe efficiency of ensemble conditional mean 
!             forecast, given obs>pthresh
!   rmsens = rms error of ensemble mean
!   crmsens =  rms error of conditional ensemble mean, given obs>pthresh

!   avgensmem = average value of all ensemble members
!   stdensmem = standard deviation of all ensemble members
!   cvensmem = coefficient of variation of all ensemble members
!   cavgensmem = conditional average value of all ensemble members
!   cstdensmem = conditional standard deviation of all ensemble members
!   ccvensmem = conditional coefficient of variation of all ensemble 
!               members
!   popensmem = pop of all ensemble memmbers
!   biasmem = (avgensmem - avgobs)/avgobs + 1
!   cbiasmem = (cavgensmem - cavgobs)/cavgobs + 1

!   bs = brier score of ensemble pop
!   bss = brier skill score of ensemble pop
!   rmspop = rms difference between avgpclass and cobspop
!   popclass(k) = limits of pop categories (increasing from 0 to 1) k=1,4
!   npclass(k) = number of forecasts in category k
!   avgpclass(k) = average forecast pop in class k
!   nposobs(k) = number of obs>pthresh given pop fcst in category k
!   cobspop(k) = obs avg pop in category k, (nposobs(k)/npclass(k))
!   expclass(k) = class limit of conditional cdf category k, k=1,4
!   ctrpclass(k) = centroid value of category k
!   nexpclass(k) = number of forecasts when obs fell in conditional 
!                  forecast proability category k
!   avgexpclass(k) = average value of forecast probability in category k
!   fracpclass(k) = fraction of forecasts falling in category k
!   cdfpclass(k) = cumulative fraction of forecasts falling in category k
!   cdfdiffmax = maximum difference between cdfpclass and uniform 
!                distribution (= Kolmogov statistic)
!   rmscdf = rms difference between cdfdiffmax and uniform distribution
!   eq_ts = equitable threat score for avgpx(i)>pthresh
!   cumestprob(i) = est probability of event being less than or equal 
!                   to obs(i)
!   npopobs = number of obs > pthresh
!   xobs(i) = list of obs > pthresh (i=1,npopobs)
!   cdfobs(i) = cum dist fn of xobs(i) (i=1,npopobs)
!   bsscpc(k) = brier skill score for cpc tercile conditional 
!               prob fcsts (k=1,3)
!   icatobs(i) = tercile category (1,2,3) of the i-th observation
!   nhits33 = number of hits of tercile categorical forecasts
!   prhit33(k,i) = forecast prob of i-th positive obs being in 
!                  k-th tercile
!   heidke = heidke skill score
!   sskuipers = Kuiper's skill score for categorical preditiction of 
!               positive obs
!   prob_err(k) = transformed rank histogram quartile value
!                 k=1 => bias error
!                 k=2 => spread error
!                 k=3 => random error

!   pval(i,j) = fraction of events, class = i,j
!   avgobsval(i,j) = average observation, class = i,j
!   stdobsval(i,j) = std dev observation, class = i,j
!   avgensval(i,j) = average ens mean, class = i,j
!   stdensval(i,j) = std dev ens mean, class = i,j
!!********************************!!******************************************************* 
subroutine ensverify(ensv1) ! (maxobs,maxmem,&
                          ! nobs,obs,nmem,xppx,pthresh,&
                          ! avgpx,stdpx,cavgpx,cstdpx,ccvpx,poppx,&
                          ! avgobs,stdobs,cvobs,&
                          ! cavgobs,cstdobs,ccvobs,popobs,&
                          ! avgens,stdens,cvens,popens,&
                          ! cavgens,cstdens,ccvens,&
                          ! avgpopens,stdpopens,&
                          ! avgcavgens,stdcavgens,&
                          ! avgccvens,stdccvens,&
                          ! biasens,cbiasens,&
                          ! corens,ccorens,rhoopt,&  !,tccorens
                          ! effens,ceffens,rmsens,crmsens,&
                          ! avgensmem,stdensmem,cvensmem,&
                          ! cavgensmem,cstdensmem,ccvensmem,popensmem,&
                          ! biasmem,cbiasmem,&
                          ! bs,bss,rmspop,&
                          ! popclass,npclass,avgpclass,nposobs,&
                          ! cobspop,&
                          ! expclass,ctrpclass,nexpclass,avgexpclass,&
                          ! fracpclass,cdfpclass,&
                          ! cdfdiffmax,rmscdf,eq_ts,&
                          ! cumestprob,npopobs,xobs,cdfobs,&
                          ! bsscpc,icatobs,&
                          ! nhits33,nfcsts33,prhit33,heidke,sskuipers,&
                          ! prob_err,pval,&
                          ! avgobsval,stdobsval,avgensval,stdensval,ifcsttran,iobstran)


use PrecStru
implicit none 
type(Ssort) ::  sort1
type(Scrps) ::  crps1
type(Sroc) ::  roc1 
type(Sbivar) :: bivar1  
type(Svtrans) :: vtrans1
type(Sensverify) :: ensv1
 
real(r8),pointer :: cpp(:)
real(r8),pointer :: pp(:)
character(c100) asRecord

integer npopobs,nensmem,minobs,nobs,nccor
integer i,j,k,ncpc,nposcpc,n33,n66,n,npop,kpos,ic
integer icx,nzero,nposensmem,nposens
integer ncs,ndp,np
integer nfcsts33,nhits33,kmax,i1,ntot
integer nens

real(r8) cstdobs,cavgobs,stdobs,avgobs,bs,stdfcst,avgfcst
real(r8) stdccvens,avgccvens,stdpopens,avgpopens,crmsens 
real(r8) rmsens,ccorens,corens,corensfcst,avg1,avg2,std1,std2,stdcavgens
real(r8) avgcavgens,finc,obs66,obs33,f33,f66
real(r8) cumpr,cumprcpc33,cumprcpc66,cumprcpc100
real(r8) pr33,pr66,cumprob,cstdens,stdensmem,popobs, ccvobs,cvobs
real(r8) cstdensmem,popensmem,ccvens,popens,rmspop,bss,bsc,cvens
real(r8) stdens,ccvensmem,ceffens,effens,cbiasmem,a,b,c,d
real(r8) sskuipers,eq_ts,exp_h,biasmem,cbiasens,biasens,travg
real(r8) trcstdobs,trcavgobs,trcstd,trcavg,pzfcst,pzobs
real(r8) trfcst,trobs,travgobs,trccor,bsscpcx,rhoopt,heidke
real(r8) rmscdf,prmax,rho22,cdfdiff,cdf,cdfdiffmax,vtrans 
real(r8) cvensmem,pr100,cavgens,avgens,cavgensmem,avgensmem
!real*4 avgobsval(2,2),stdobsval(2,2),pval(2,2)
!real*4 avgensval(2,2),stdensval(2,2),prob_err(3)
logical iclass
integer*4 nval(2,2) !nposobs(4),npclass(4),nexpclass(4),
!real*4 cobspop(4),popclass(4) !avgpclass(4),
!real*4 expclass(4),avgexpclass(4),ctrpclass(4) 
!real*4 fracpclass(4),cdfpclass(4) 
integer*4 nper(4)
real(r8) bscpc(3) !,bsscpc(3)

!parameter (maxnobs=5000)
!real*4  poppx(maxobs),avgpx(maxobs),pp(maxmem)
!real*4  cavgpx(maxobs),cstdpx(maxobs), ccvpx(maxobs),stdpx(maxobs)
      
!integer*4 nmem(maxobs),npos(maxnobs)
!real*4 obs(maxobs),xppx(maxmem,maxobs),
!real*4 cdfobs(maxobs),cumestprob(maxobs)
!real*4 prhit33(3,maxobs)
!integer*4 icatobs(maxobs)

 
 
  nobs=ensv1.nobs
  if (ensv1.nobs.le.0) then
    asRecord="number of observations/forecasts is less than 1 in  ensemble verify "
    call Wrecord(asRecord)
    return
  endif
  allocate(cpp(nobs))
  allocate(pp(ensv1.nmem))
    !  initialize statistics
  minobs=20
  ensv1.popclass(1) = 0.25
  ensv1.popclass(2) = 0.50
  ensv1.popclass(3) = 0.75
  ensv1.popclass(4) = 1.00
  do i=1,4
    ensv1.npclass(i) = 0
    ensv1.avgpclass(i) = 0.
    ensv1.nposobs(i) = 0
  enddo
  ensv1.expclass(1) = 0.25
  ensv1.expclass(2) = 0.50
  ensv1.expclass(3) = 0.75
  ensv1.expclass(4) = 1.00
  ensv1.ctrpclass(1) = 0.125
  ensv1.ctrpclass(2) = 0.375
  ensv1.ctrpclass(3) = 0.635
  ensv1.ctrpclass(4) = 0.875
  do i=1,4
    ensv1.nexpclass(i) = 0
    ensv1.avgexpclass(i) = 0.
  enddo

  bs = 0.

  avgfcst = 0.
  stdfcst = 0.
  
  avgobs = 0.
  stdobs = 0.
  npopobs = 0
  cavgobs = 0.
  cstdobs = 0.
   
  nensmem = 0
  avgensmem = 0.
  stdensmem = 0.
  cavgensmem = 0.
  cstdensmem = 0.
  nposensmem = 0

  avgens = 0.
  stdens = 0.
  nens = 0
  nposens = 0
  cavgens = 0.
  cstdens = 0.
  avgcavgens = 0.
  stdcavgens = 0.
  avgccvens = 0.
  stdccvens = 0.
  avgpopens = 0.
  stdpopens = 0.

  rmsens = 0.
  crmsens = 0.
  corens = 0.
  corensfcst=0.
  ccorens = 0.
  nccor = 0
  avg1 = 0.
  std1 = 0.
  avg2 = 0.
  std2 = 0.

  bscpc(1) = 0.
  bscpc(2) = 0.
  bscpc(3) = 0.

  ncpc = 0
  k = 0
  do i=1,nobs
    ensv1.icatobs(i) = 0
    if (ensv1.obs(i).gt.ensv1.pthresh) then
      k = k + 1
      ensv1.xobs(k) = ensv1.obs(i)
      nposcpc = k
    endif
  enddo
  sort1.n=nposcpc
  sort1.a=>ensv1.xobs
  call sort (sort1)
  !call sort (nposcpc,xobs)
  n33 = nposcpc/3
  n66 = (2*nposcpc)/3
  finc = 1./float(nposcpc + 1)
  f33 = n33*finc
  f66 = n66*finc
  if (n33.gt.0) then
    obs33 = ensv1.xobs(n33) + (ensv1.xobs(n33+1)-ensv1.xobs(n33))*((.3333-f33)/finc)
  endif
  if (n66.gt.0) then
    obs66 = ensv1.xobs(n66) + (ensv1.xobs(n66+1)-ensv1.xobs(n66))*((.6667-f66)/finc)
  endif
  do k=1,3
    bscpc(k) = 0.
  enddo

  !c  ...accumulate statistics

  kpos = 0
  do i=1,nobs
    avgobs = avgobs + ensv1.obs(i)
    stdobs = stdobs + ensv1.obs(i)**2
    avgfcst = avgfcst + ensv1.fcst(i)
    stdfcst = stdfcst + ensv1.fcst(i)**2

    if (ensv1.obs(i).gt.ensv1.pthresh) then
      npopobs = npopobs + 1
      cavgobs = cavgobs +ensv1.obs(i)
      cstdobs = cstdobs +ensv1.obs(i)**2
    endif
    do j=1,ensv1.nmem
      pp(j) =ensv1.xppx(j,i)
      nensmem = nensmem + 1
      avgensmem = avgensmem +pp(j)
      stdensmem = stdensmem +pp(j)**2
      if (pp(j).gt.ensv1.pthresh) then
        cavgensmem = cavgensmem +pp(j)
        cstdensmem = cstdensmem +pp(j)**2
        nposensmem = nposensmem + 1
      endif
    enddo
    ensv1.avgpx(i) = 0.
    ensv1.stdpx(i) = 0.
    npop = 0
    ensv1.cavgpx(i) = 0.
    ensv1.cstdpx(i) = 0.
    do k=1,ensv1.nmem
      ensv1.avgpx(i) =ensv1.avgpx(i) +pp(k)
      ensv1.stdpx(i) = ensv1.stdpx(i) +pp(k)**2
      if (pp(k).gt.ensv1.pthresh) then
        npop = npop + 1
        cpp(npop) =pp(k)
        ensv1.cavgpx(i) = ensv1.cavgpx(i) +pp(k)
        ensv1.cstdpx(i) = ensv1.cstdpx(i) +pp(k)**2
 !=====================================================================================
 !     else
 !       npop =npop 
 !=====================================================================================
      endif
    enddo
    if (npop.gt.1) then
    sort1.n=npop
    sort1.a=>cpp
    call sort (sort1)
    !call sort (npop,cpp)
    endif
    if (ensv1.nmem.gt.0) then
      ensv1.avgpx(i) =ensv1.avgpx(i)/ensv1.nmem
      ensv1.stdpx(i) = sqrt(abs(ensv1.stdpx(i)/ensv1.nmem -ensv1.avgpx(i)**2))
      ensv1.poppx(i) = float(npop)/float(ensv1.nmem)
    else
      ensv1.avgpx(i) = 0.
      ensv1.stdpx(i) = 0.
      ensv1.poppx(i) = 0.
    endif
    if (npop.gt.0) then
      ensv1.cavgpx(i) = ensv1.cavgpx(i)/npop
      ensv1.cstdpx(i) = sqrt(abs(ensv1.cstdpx(i)/npop-ensv1.cavgpx(i)**2))
      ensv1.ccvpx(i) = ensv1.cstdpx(i)/ensv1.cavgpx(i)
    else
      ensv1.cavgpx(i) = 0.
      ensv1.cstdpx(i) = 0.
      ensv1.ccvpx(i) = 0.
    endif
     
    avgens = avgens +ensv1.avgpx(i)
    stdens = stdens +ensv1.avgpx(i)**2
    corens = corens +ensv1.avgpx(i)*ensv1.obs(i)
    corensfcst = corensfcst +ensv1.avgpx(i)*ensv1.fcst(i)
    rmsens = rmsens + (ensv1.avgpx(i)-ensv1.obs(i))**2
    avgpopens = avgpopens +ensv1.poppx(i)
    stdpopens = stdpopens +ensv1.poppx(i)**2
    if (ensv1.avgpx(i).gt.ensv1.pthresh) then
      nposens = nposens + 1
      cavgens = cavgens +ensv1.avgpx(i)
      cstdens = cstdens +ensv1.avgpx(i)**2
      avgcavgens = avgcavgens + ensv1.cavgpx(i)  !  cond avg fcst when obs>pthresh
      stdcavgens = stdcavgens + ensv1.cavgpx(i)**2
      avgccvens = avgccvens + ensv1.ccvpx(i)
      stdccvens = stdccvens + ensv1.ccvpx(i)**2
      if (ensv1.obs(i).gt.ensv1.pthresh) then
        ccorens = ccorens + ensv1.cavgpx(i)*ensv1.obs(i)
        avg1 = avg1 + ensv1.cavgpx(i)
        std1 = std1 + ensv1.cavgpx(i)**2
        avg2 = avg2 +ensv1.obs(i)
        std2 = std2 +ensv1.obs(i)**2
        nccor = nccor + 1
      endif
    endif

    if (ensv1.obs(i).gt.ensv1.pthresh) then
      bs = bs + (ensv1.poppx(i) - 1)**2
    else
      bs = bs +ensv1.poppx(i)**2
    endif
    nzero = ensv1.nmem - npop
    iclass = .true.
    do ic=1,4
      if (ensv1.poppx(i).le.ensv1.popclass(ic).and.iclass)then
        iclass = .false.
        ensv1.npclass(ic) = ensv1.npclass(ic) + 1
        ensv1.avgpclass(ic) = ensv1.avgpclass(ic) +ensv1.poppx(i)
        if (ensv1.obs(i).gt.ensv1.pthresh) ensv1.nposobs(ic) = ensv1.nposobs(ic) + 1
      endif
    enddo
    ! write (iout2,'(3f8.2)') fcst2(i),avgpp,obs(i)
    if (ensv1.obs(i).gt.ensv1.pthresh) then
      kpos = kpos + 1
      crmsens = crmsens + (ensv1.cavgpx(i)-ensv1.obs(i))**2
      if (npop.ge.1) then
	    crps1.n=npop
        crps1.x0=ensv1.obs(i)
        crps1.xx=>cpp
        cumpr = cumprob(crps1) !(obs(i),npop,cpp)
      else
        cumpr = 1.0
      endif
      ensv1.cumestprob(kpos) = cumpr
      ensv1.cdfobs(kpos) = cumpr
      icx = 0
      do ic=1,3
        if (cumpr.gt.ensv1.expclass(ic)) icx = ic
      enddo  
      icx = icx + 1
      ensv1.nexpclass(icx) = ensv1.nexpclass(icx) + 1
      ensv1.avgexpclass(icx) = ensv1.avgexpclass(icx) + cumpr
      ncpc = ncpc + 1
      if (npop.gt.0) then
  	    crps1.n=npop
        crps1.x0=obs33
        crps1.xx=>cpp
        cumprcpc33 = cumprob(crps1) !(obs33,npop,cpp)
        crps1.x0=obs66
        cumprcpc66 = cumprob(crps1) !(obs66,npop,cpp)
        cumprcpc100 = 1.
      else
        cumprcpc33 = 0.
        cumprcpc66 = 0.
        cumprcpc100 = 0.
      endif
      pr33 = cumprcpc33
      pr66 = cumprcpc66 - cumprcpc33
      pr100 = cumprcpc100 - cumprcpc66
      ensv1.prhit33(1,kpos) = pr33
      ensv1.prhit33(2,kpos) = pr66
      ensv1.prhit33(3,kpos) = pr100
      if (ensv1.obs(i).le.obs33) then
        bscpc(1) = bscpc(1) + (pr33-1.)**2
        bscpc(2) = bscpc(2) + (pr66)**2
        bscpc(3) = bscpc(3) + (pr100)**2
        ensv1.icatobs(i) = 1    
      elseif (ensv1.obs(i).le.obs66) then
        bscpc(1) = bscpc(1) + (pr33)**2
        bscpc(2) = bscpc(2) + (pr66-1.)**2
        bscpc(3) = bscpc(3) + (pr100)**2
        ensv1.icatobs(i) = 2    
      else
        bscpc(1) = bscpc(1) + (pr33)**2
        bscpc(2) = bscpc(2) + (pr66)**2
        bscpc(3) = bscpc(3) + (pr100-1.)**2
        ensv1.icatobs(i) = 3    
      endif
    endif
  enddo
  do icx=1,4
    if (ensv1.nexpclass(icx).gt.0) then 
      ensv1.avgexpclass(icx) = ensv1.avgexpclass(icx)/ensv1.nexpclass(icx)
    else
      ensv1.avgexpclass(icx) = 0.
    endif
  enddo
  sort1.n=nobs
  sort1.a=>ensv1.cumestprob
  call sort (sort1)
  !call sort (nobs,cumestprob)
  sort1.n=kpos
  sort1.a=>ensv1.cdfobs
  call sort (sort1)
  !call sort (kpos,cdfobs)
  !c   ...compute statistics
  if (nobs.gt.0) then
    avgfcst = avgfcst/nobs
    stdfcst = sqrt(abs(stdfcst/nobs - avgfcst**2))
    avgobs = avgobs/nobs
    stdobs = sqrt(abs(stdobs/nobs - avgobs**2))
    cvobs = stdobs/avgobs
  else
    avgfcst = 0.
    stdfcst = 0.
    avgobs = 0.
    stdobs = 0.
    cvobs = 0.
  endif
  ccvobs = 0.
  popobs = 0.
  if (npopobs.gt.0) then
    cavgobs = cavgobs/npopobs
    cstdobs = sqrt(abs(cstdobs/npopobs - cavgobs**2))
    ccvobs = cavgobs/cstdobs
    popobs = float(npopobs)/float(nobs)
  else
    cavgobs = 0.
    cstdobs = 0.
    ccvobs = 0.
    popobs = 0.
  endif

  if (nensmem.gt.0) then
    avgensmem = avgensmem/nensmem
    stdensmem = sqrt(abs(stdensmem/nensmem - avgensmem**2))
    cvensmem = stdensmem/avgensmem
    popensmem = float(nposensmem)/float(nensmem)
  else
    avgensmem = 0.
    stdensmem = 0.
    cvensmem = 0.
    popensmem = 0.
  endif
  if (nposensmem.gt.0) then
    cavgensmem = cavgensmem/nposensmem
    cstdensmem = sqrt(abs(cstdensmem/nposensmem - cavgensmem**2))       
    ccvensmem = cstdensmem/cavgensmem
  else
    cavgensmem = 0.
    cstdensmem = 0.
    ccvensmem = 0.
  endif

  if (nobs.gt.0) then
    avgens = avgens/nobs
    stdens = sqrt(abs(stdens/nobs - avgens**2))
    cvens = stdens/avgens
    popens = float(nposens)/float(nobs)
    avgpopens = avgpopens/nobs
    stdpopens = sqrt(abs(stdpopens/nobs - avgpopens**2))
  else
    avgens = 0.
    stdens = 0.
    cvens = 0.
    popens = 0.
    avgpopens = 0.
    stdpopens = 0.
  endif

  if (nposens.gt.0) then
    cavgens = cavgens/nposens
    cstdens = sqrt(abs(cstdens/nposens-cavgens**2))
    ccvens = cstdens/cavgens
    avgccvens = avgccvens/nposens
    stdccvens = sqrt(abs(stdccvens/nposens -avgccvens**2)) 
         
  else
    cavgens = 0.
    cstdens = 0.
    ccvens = 0.
    avgccvens = 0.
    stdccvens = 0.
  endif

  if (nobs.gt.0.and.stdens*stdobs*stdfcst.gt.0) then
    corens = (corens/nobs - avgens*avgobs)/(stdens*stdobs)   
	if (corens.lt.0.1)   &
	  corens=corens
    corensfcst = (corensfcst/nobs - avgens*avgfcst)/(stdens*stdfcst)     	 
  else
    corens = 0.
    corensfcst = 0.
  endif

  if (nccor.gt.0) then
    avg1 = avg1/nccor
    std1 = sqrt(abs(std1/nccor - avg1**2))
    avg2 = avg2/nccor
    std2 = sqrt(abs(std2/nccor - avg2**2))    
    ccorens = (ccorens/nccor - avg1*avg2)/(std1*std2)
  else
    avg1 = 0.
    std1 = 0.
    avg2 = 0.
    std2 = 0.
    ccorens = 0.
  endif

  bsc = popobs*(1.-popobs)
  if (nobs.gt.0) then
    bs = bs/nobs
  else
    bs = 0.
  endif
  if (popobs.lt.1.0) then
    bss = 1 - (bs/bsc)
  else
    bss = 1.0
  endif

  rmspop = 0.
  ncs = 0
  do ic=1,4
    if (ensv1.npclass(ic).gt.0) then
      ensv1.avgpclass(ic) = ensv1.avgpclass(ic)/ensv1.npclass(ic)
      ensv1.cobspop(ic) = float(ensv1.nposobs(ic))/float(ensv1.npclass(ic))
    endif
    if (ensv1.npclass(ic).gt.10) then
      rmspop = rmspop + (ensv1.avgpclass(ic) - ensv1.cobspop(ic))**2
      ncs = ncs + 1
    endif
  enddo
  if (ncs.gt.0) then
    rmspop = sqrt(rmspop/ncs)
  else
    rmspop = 0.
  endif

  ntot = 0
  do ic=1,4
    ntot = ntot + ensv1.nexpclass(ic)
  enddo
  do ic=1,4
    ensv1.fracpclass(ic) = float(ensv1.nexpclass(ic))/float(ntot)
    if (ic.eq.1) then
      ensv1.cdfpclass(ic) = ensv1.fracpclass(ic)
    else
      ensv1.cdfpclass(ic) = ensv1.cdfpclass(ic-1) + ensv1.fracpclass(ic)
    endif
  enddo
  ensv1.prob_err(1) = (ensv1.fracpclass(1) + ensv1.fracpclass(2) - ensv1.fracpclass(3) - ensv1.fracpclass(4))/2.     
  ensv1.prob_err(2) = (ensv1.fracpclass(1) - ensv1.fracpclass(2) -ensv1.fracpclass(3) + ensv1.fracpclass(4))/2.     
  ensv1.prob_err(3) = (ensv1.fracpclass(1) - ensv1.fracpclass(2) +ensv1.fracpclass(3) - ensv1.fracpclass(4))/2.     
 
  biasens = avgens/avgobs  !(avgens - avgobs)/avgobs + 1
  cbiasens = cavgens/cavgobs !(cavgens - cavgobs)/cavgobs + 1
  biasmem = avgensmem/avgobs !(avgensmem - avgobs)/avgobs + 1
  cbiasmem = cavgensmem/cavgobs ! (cavgensmem - cavgobs)/cavgobs + 1
  
  if (nobs.gt.0) then
    rmsens = sqrt(rmsens/nobs)
!*********************************************************************************
!	if (rmsens.gt.10) &   
!      rmsens =rmsens 
!*********************************************************************************
  else
    rmsens = 0.
  endif
  if (npopobs.gt.0) then
    crmsens = sqrt(crmsens/npopobs)
  else
    crmsens = 0.
  endif
 
  if (stdobs.gt.0.) then
    effens = 1. - rmsens**2/stdobs**2
  else
    effens = 0.
  endif
  if (cstdobs.gt.0.) then
    ceffens = 1. - crmsens**2/cstdobs**2
  else
    ceffens = 0.
  endif
 
  !   ...compute contingency table...
 
  do i=1,2
    do j=1,2
    nval(i,j) = 0
    ensv1.avgobsval(i,j) = 0.
    ensv1.avgensval(i,j) = 0.
    ensv1.stdobsval(i,j) = 0.
    ensv1.stdensval(i,j) = 0.
    enddo
  enddo
  do k=1,nobs
    if (ensv1.avgpx(k).gt.ensv1.pthresh) then
      i = 2
    else
      i = 1
    endif
    if (ensv1.obs(k).gt.ensv1.pthresh) then
      j = 2
    else
      j = 1
    endif
    nval(i,j) = nval(i,j) + 1
    ensv1.avgensval(i,j) = ensv1.avgensval(i,j) +ensv1.avgpx(i)
    ensv1.avgobsval(i,j) = ensv1.avgobsval(i,j) +ensv1.obs(i)
    ensv1.stdensval(i,j) = ensv1.stdensval(i,j) +ensv1.avgpx(i)**2
    ensv1.stdobsval(i,j) = ensv1.stdobsval(i,j) +ensv1.obs(i)**2
  enddo
  do i=1,2
    do j=1,2
      if (nval(i,j).gt.0) then
        ensv1.avgensval(i,j) = ensv1.avgensval(i,j)/nval(i,j)
        ensv1.avgobsval(i,j) = ensv1.avgobsval(i,j)/nval(i,j)
        ensv1.stdensval(i,j) = sqrt(abs(ensv1.stdensval(i,j)/nval(i,j) -ensv1.avgensval(i,j)**2))           
        ensv1.stdobsval(i,j) = sqrt(abs(ensv1.stdobsval(i,j)/nval(i,j) -ensv1.avgobsval(i,j)**2))            
      endif
    enddo
  enddo
  ndp = 0
  do i=1,2
    do j=1,2
     ndp = ndp + nval(i,j)
    enddo
  enddo
  do i=1,2
    do j=1,2
      if (ndp.gt.0) then
        ensv1.pval(i,j) = float(nval(i,j))/float(ndp)
      else
        ensv1.pval(i,j) = 0.
      endif
    enddo
  enddo
  a = nval(2,2)
  b = nval(2,1)
  c = nval(1,2)
  d = nval(1,1)
  exp_h = ((a+b)*(a+c))/(a+b+c+d)
  eq_ts = (a - exp_h)/(a + b + c - exp_h)
  sskuipers = (ensv1.pval(1,1)*ensv1.pval(2,2) - ensv1.pval(1,2)*ensv1.pval(2,1))/&
	       ((ensv1.pval(1,1)+ensv1.pval(2,1))*(ensv1.pval(1,2)+ensv1.pval(2,2)))
   
  !   ...transform obs (Weibull) and ens mean (gamma distr) to snd

  trccor = 0.
  trcavg = 0.
  trcstd = 0.
  trcavgobs = 0.
  trcstdobs = 0.
  travg = 0.
  travgobs = 0.
  np = 0  
  pzobs = 1. - popobs
  pzfcst = 1. - popens
  do i=1,nobs
    if (ensv1.obs(i).gt.ensv1.pthresh.and.ensv1.avgpx(i).gt.ensv1.pthresh) then
 
      vtrans1.avg=cavgobs
      vtrans1.cv=ccvobs
      vtrans1.pz=pzobs
      vtrans1.obs=ensv1.obs(i)
      vtrans1.ivopt=ensv1.iobstran    
      trobs = vtrans(vtrans1)!(cavgobs,ccvobs,pzobs,obs(i),iobstran)

      vtrans1.avg=cavgens !cavgfcst  here change the Schaake code 5-15-2012 by Aizhong Ye
      vtrans1.cv=ccvens !ccvfcst
      vtrans1.pz=pzfcst
      vtrans1.obs=ensv1.avgpx(i)
      vtrans1.ivopt=ensv1.ifcsttran
      trfcst = vtrans(vtrans1) !(cavgfcst,ccvfcst,pzfcst,avgpx(i),ifcsttran)
      
      trcavg = trcavg + trfcst
      trcstd = trcstd + trfcst**2
      trcavgobs = trcavgobs + trobs
      trcstdobs = trcstdobs + trobs**2
      trccor = trccor + trfcst*trobs
      np = np + 1
    endif
  enddo
  if (np.gt.minobs) then
    trcavg = trcavg/np
    trcstd = sqrt(abs(trcstd/np - trcavg**2))
    trcavgobs = trcavgobs/np
    trcstdobs = sqrt(abs(trcstdobs/np - trcavgobs**2))
    trccor = (trccor/np - trcavg*trcavgobs)/(trcstd*trcstdobs)    
  else
    trccor = corens
  endif
  rho22=trccor  
  if (ensv1.pval(1,1).lt.0.05.or. &
    ensv1.pval(1,1).gt.0.95.or. &
    ensv1.pval(2,2).lt.0.05.or. &
    ensv1.pval(2,2).gt.0.95) then
    rhoopt = corens
  else
    bivar1.p11=ensv1.pval(1,1)
    bivar1.p12=ensv1.pval(1,2)
    bivar1.p21=ensv1.pval(2,1)
    bivar1.p22=ensv1.pval(2,2)
    bivar1.rho22=rho22
    bivar1.rho=rhoopt
    call estrho3(bivar1)  
    rhoopt=bivar1.rho  
    !call estrho3(iout,pval(1,1),pval(1,2),pval(2,1),pval(2,2),rho22,rhoopt)      
  endif
 
  bsscpcx = .3333*.6667
  if (ncpc.gt.0) then
    do k=1,3
      bscpc(k) = bscpc(k)/ncpc
      ensv1.bsscpc(k) = 1. - (bscpc(k)/bsscpcx)
    enddo
  else
    do k=1,3
      bscpc(k) = 0.
      ensv1.bsscpc(k) = 0.
    enddo
  endif
 
  nhits33 = 0
  nfcsts33 = 0
  do i=1,npopobs
    prmax = 0.
    kmax = 0
    do k=1,3
      if (ensv1.prhit33(k,i).gt.prmax) then
        prmax = ensv1.prhit33(k,i)
        kmax = k
      endif
    enddo
    if (kmax.eq.1.or.kmax.eq.3) then
      nfcsts33 = nfcsts33 + 1
      if (kmax.eq.1.and.ensv1.cdfobs(i).lt.0.3333) then
        nhits33 = nhits33 + 1
      endif
      if (kmax.eq.3.and.ensv1.cdfobs(i).gt.0.6666) then
        nhits33 = nhits33 + 1
      endif
    endif
  enddo
  heidke = 100*(float(nhits33)-.3333)/(float(nfcsts33)-.3333)
 
  rmscdf = 0.
  cdfdiffmax = 0.
  cdf = 0.
  i1 = nobs - npopobs
  do i=1,npopobs
    cdf = float(i)/(float(npopobs) + 1)
    cdfdiff = abs(ensv1.cdfobs(i) - cdf)
    rmscdf = rmscdf + cdfdiff**2
    if (cdfdiff.gt.cdfdiffmax) cdfdiffmax = cdfdiff
  enddo
  if (npopobs.gt.0) then
    rmscdf = sqrt(rmscdf/npopobs)
  else
    rmscdf = 0.
  endif
  !write (iout,'(i5,11f7.3)') nobs,avgobs,avgens,stdobs,stdens,corens,rmsens, &
  !       effens,rmscdf,(bsscpc(i),i=1,3) 
 
   
  deallocate(cpp)
  deallocate(pp)
  ensv1.avgobs=avgobs
  ensv1.stdobs=stdobs
  ensv1.cvobs=cvobs
  ensv1.cavgobs=cavgobs
  ensv1.cstdobs=cstdobs
  ensv1.ccvobs=ccvobs
  ensv1.popobs=popobs
  ensv1.avgens=avgens
  ensv1.stdens=stdens
  ensv1.cvens=cvens
  ensv1.popens=popens
  ensv1.cavgens=cavgens
  ensv1.cstdens=cstdens
  ensv1.ccvens=ccvens
  ensv1.avgpopens=avgpopens
  ensv1.stdpopens=stdpopens
  ensv1.avgcavgens=avgcavgens
  ensv1.stdcavgens=stdcavgens
  ensv1.avgccvens=avgccvens
  ensv1.stdccvens=stdccvens
  ensv1.biasens=biasens
  ensv1.cbiasens=cbiasens
  ensv1.corens=corens
  ensv1.ccorens=ccorens
  ensv1.corensfcst=corensfcst

  ensv1.rhoopt=rhoopt
  ensv1.effens=effens
  ensv1.ceffens=ceffens
  ensv1.rmsens=rmsens
  ensv1.crmsens=crmsens
  ensv1.avgensmem=avgensmem
  ensv1.stdensmem=stdensmem
  ensv1.cvensmem=cvensmem
  ensv1.cavgensmem=cavgensmem
  ensv1.cstdensmem=cstdensmem
  ensv1.ccvensmem=ccvensmem
  ensv1.popensmem=popensmem
  ensv1.biasmem=biasmem
  ensv1.cbiasmem=cbiasmem
  ensv1.bs=bs
  ensv1.bss=bss
  ensv1.rmspop=rmspop
  ensv1.cdfdiffmax=cdfdiffmax
  ensv1.rmscdf=rmscdf
  ensv1.eq_ts=eq_ts
  ensv1.npopobs=npopobs
  ensv1.nhits33=nhits33
  ensv1.nfcsts33=nfcsts33
  ensv1.heidke=heidke
  ensv1.sskuipers=sskuipers
  ensv1.heidke=heidke
  ensv1.sskuipers=sskuipers


  return
end


!!********************************!!******************************************************* 
!55 Description: Subroutine to compute realive operating characteristic statistics
!  
! Original Author:  John Schaake
! File Creation Date:  March 14, 2007
! Development Group:  HEP
! 
! Calling Arguments:
!   Name     Input/Output   Type        Description
!   nfcst       input       integer*4   number of ensemble forecasts w/ obs
!   obs         input       real*4      obs(i) = value of observation for i-th forecast
!   maxmem      input       integer*4   maximum number of members for array xx
!   nmem        input       integer*4   actual number of members in array xx
!   xx          input       real*4      xx(j,i) = j-th member value for ensemble forecast i 
!   nthresh     input       integer*4   number of observations thresholds  to be used

!   obsthresh   output      real*4      obsthresh(i) = observation threshold (i=1,nthresh)
!   hr          output      real*4      hr(i) = hit rate for i-th observation threshold 
!   far         output      real*4      far(i) = false alarm rate for i-th observation threshold 
!   roc_area    output      real*4      area of skill for relative operating characteristics curve 
!   roc_score   output      real*4      relative ro! skill = roc_area/0.5
!   istatus     output      integer*4   error status (=0, no errors)
! Required Files/Databases:  
!       cumprob.for     sort.for
! Error Codes/Exceptions:  Listed at end of module
! Modification History:
!  Date      Developer      Action
!    
!!********************************!!******************************************************* 
subroutine roc(roc1) !(nfcst,obs,maxmem,nmem,xx, &
                     !nthresh,obsthresh, &
                     !hr,far,roc_area,roc_score,istatus)
use PrecStru
implicit none 
type(Ssort) ::  sort1
type(Scrps) ::  crps1
type(Sroc) ::  roc1 
!compute probabilisti! verification statistics
integer maxobs,nmem,maxlevels,nfcst,nthresh,nobs
parameter (maxobs=5000)
parameter (maxlevels=30)
real(r8) cumprob 
!real(r8) xfcst(maxobs)
real(r8),pointer :: obsthresh(:)
real(r8),pointer :: sortobs(:)
real(r8),pointer :: sortfcst(:)

real(r8),pointer :: xx(:,:)! xx(j,i) = j-th member value for ensemble forecast i xx(maxmem,nfcst),
real(r8),pointer :: obs(:) ! obs(i) = value of observation for i-th forecast
real(r8),pointer :: hr(:)  ! hr(i) = hit rate for i-th observation threshold hr(nthresh),far(nthresh)
real(r8),pointer :: far(:) ! far(i) = false alarm rate for i-th observation threshold 

!real  sortobs(maxobs)
real*4  sortenspr(maxobs),enspr(maxobs)
real*4 a(maxlevels),b(maxlevels),c(maxlevels),d(maxlevels)
real*4 da,roc_area ,roc_score
integer i,j,k,istatus,nobsbx

  nfcst=roc1.nfcst
  nthresh=roc1.nthresh
  obsthresh=>roc1.obsthresh
  nobsbx=0
  nobs=0
  nmem=roc1.nmem

  xx=>roc1.xx
  obs=>roc1.obs
  hr=>roc1.hr
  far=>roc1.far 
  istatus = 0
  allocate(sortobs(nfcst))
  allocate(sortfcst(nfcst))

  !allocate(obsthresh(nthresh))
  !allocate(hr(nthresh))
  !allocate(far(nthresh))

  ! establish values of the observation thresholds
  do i=1,nfcst
    sortobs(i) = obs(i)
  enddo
  sort1.n=nfcst
  sort1.a=>sortobs
  call sort (sort1)
  !call sort (nfcst,sortobs)
  do i=1,nthresh
    j = nint(float(nfcst)*((float(i))/float(nthresh+1)))
    if (j.lt.1) j = 1
    if (j.gt.nfcst) j = nfcst
    obsthresh(i) = sortobs(j)
  enddo
  !For each obs threshold, compute hit rate and false alarm rate
  do i=1,nthresh
    ! initialize contingency table
    a(i) = 0.
    b(i) = 0.
    c(i) = 0.
    d(i) = 0.
    !  process each forecast
    do j=1,nfcst
      ! get cumulative distribution of ensemble members
      do k=1,nmem
        sortfcst(k) = xx(k,j)
      enddo
      sort1.n=nmem
      sort1.a=>sortfcst
      call sort (sort1)
      !call sort (nmem,sortfcst)
      ! get non-exceedance probability of observation threshold
      crps1.n=nmem
      crps1.x0=obsthresh(i)
      crps1.xx=>sortfcst
      enspr(j) = cumprob(crps1)!(obsthresh(i),nmem,sortfcst)
      if (obs(j).le.obsthresh(i)) then
        nobs = nobs + 1
        a(i) = a(i) + enspr(j)
        b(i) = b(i) + (1. - enspr(j))
      else
        nobsbx = nobsbx + 1
        c(i) = c(i) + enspr(j)
        d(i) = d(i) + (1. - enspr(j))
      endif
    enddo
    hr(i) = a(i)/(a(i) + b(i))
    far(i) = c(i)/(c(i) + d(i))
  enddo
  !c compute roc statistics
  da = 0.5*hr(1)*far(1)
  roc_area = da
  do i=2,nthresh
    da = 0.5*(far(i)-far(i-1))*(hr(i)+hr(i-1))
    roc_area = roc_area + da
  enddo
  da = 0.5*(1.-far(nthresh))*(1.+hr(nthresh))
  roc_area = roc_area + da
  roc_area = roc_area - 0.5
  roc_area = max(roc_area,0.)
  roc_score = roc_area/0.5

  roc1.roc_area=roc_area
  roc1.roc_score=roc_score
  roc1.istatus=istatus 

  deallocate(sortobs)   
  deallocate(sortfcst)   
   
  
  return
end

!!********************************!!******************************************************* 
!54 find cumulative probability position of x in xx
!   n = number of values of xx (sorted in increasing order)
!!********************************!!******************************************************* 
function cumprob(crps1) ! (x,n,xx)
use PrecStru
implicit none 
type(Ssort) ::  sort1
type(Scrps) ::  crps1
real(r8),pointer :: xx(:)
integer n,m
real(r8) x,pr1,pr2,x1,x2,dx,cumprob
!real*4 xx(n)
logical ifound

  n=crps1.n
  x=crps1.x0
  xx=>crps1.xx

  if (x.le.xx(1)) then
    pr1 = 0.
    pr2 = 1./float(2*n) 
    if (x.ge.0.and.xx(1).ge.0.) then
      x1 = 0.
    else
      x1 = xx(1) - (xx(2) - xx(1))
    endif
    x2 = xx(1)
  elseif (x.ge.xx(n)) then
    pr1 = 1. - 1./float(2*n)
    pr2 = 1.
    x1 = xx(n)
    x2 = 2*xx(n)
  else
    ifound = .true.
    m = 2
    do while (ifound)
      if (x.gt.xx(m)) then
        m = m+1
        if (m.ge.n) ifound = .false.
      else
        ifound = .false.
      endif
    enddo
    pr1 = (m - 1 - 0.5)/float(n)
    pr2 = (m - 0.5)/float(n)
    x1 = xx(m-1)
    x2 = xx(m)
  endif
  dx = x2 - x1
  if (x.lt.x1) then
    cumprob = 0.
  elseif (x.gt.x2) then
      cumprob = 1.
  else
    if (dx.gt.0.) then
      cumprob = pr1 + (x - x1)*((pr2 - pr1)/dx)
    else
      cumprob = pr1
    endif
  endif
  if (cumprob.lt.0.) cumprob = 0.
  if (cumprob.gt.1.) cumprob = 1.
  return
end
!!********************************!!******************************************************* 
!53 compute continuous rank probability skill scores
!!********************************!!******************************************************* 
subroutine crps_eval(crpseval1) ! (nobs,fcst,obs,ensem,nmemdim,nmem,pthresh, &
                            !crps_avgobs,crps_clim,crps_fcst,&
                            !crps_ensavg,crps_ens,crps_ensmed,&
                            !crps_cavgobs,crps_cclim,crps_cfcst,&
                            !crps_censavg,crps_cens,crps_censmed)
use PrecStru
implicit none 
type(Ssort) ::  sort1
type(Scrps) ::  crps1
type(Scrpseval) ::  crpseval1
integer nobs,n,i,j,npos,nn 
real(r8),pointer :: fcst(:)
real(r8),pointer :: obs(:)
real(r8),pointer :: ensem(:,:)

!real  fcst(nobs),obs(nobs),ensem(nmemdim,nobs) 
!real  xobs(5000),obspos(5000),ensavg(5000),censavg(5000)
real(r8),pointer :: xobs(:)
real(r8),pointer :: obspos(:)
real(r8),pointer :: ensavg(:)
real(r8),pointer :: censavg(:)
integer  nmem

!real*4 x(5000),cx(5000)
integer,pointer ::  nposens(:)
real(r8),pointer :: xsort(:)
real(r8),pointer :: xensem(:)
real(r8),pointer :: ensmedian(:)
real(r8),pointer :: censmedian(:)
!real  xsort(5000),xensem(5000),ensmedian(5000),censmedian(5000)
 
real(r8),pointer :: x(:)
real(r8),pointer :: cx(:)
real(r8) crpsx,crpsx1,crpsx2,crpsx3,crpsx4,crpsx5,crpsx6
real(r8) avgobs,cavgobs,crps_fcst,crps_avgobs,crps_clim,crps_ensavg,crps_ens
real(r8) crps,crps_censmed,crps_cens,crps_censavg,crps_cclim,crps_cavgobs
real(r8) crps_cfcst,crps_ensmed,pthresh

  nobs=crpseval1.nobs
   
  nmem=crpseval1.nmem
  pthresh=crpseval1.pthresh
  fcst=>crpseval1.fcst
  obs=>crpseval1.obs
  ensem=>crpseval1.ensem


  allocate(x(nobs))
  allocate(cx(nobs))

  allocate(xsort(nobs))
  allocate(xensem(nobs))
  allocate(ensmedian(nobs))
  allocate(censmedian(nobs))

  allocate(xobs(nobs))
  allocate(obspos(nobs))
  allocate(ensavg(nobs))
  allocate(censavg(nobs))
  allocate(nposens(nobs))
 
 
  avgobs = 0.
  cavgobs = 0.
  npos = 0
 
  do i=1,nobs
    avgobs = avgobs + obs(i)
    nposens(i) = 0
    ensavg(i) = 0.
    censavg(i) = 0.
    censmedian(i) = 0.
    ensmedian(i) = 0.
    n = 0
    do j=1,nmem
      xensem(j) = ensem(j,i)
      ensavg(i) = ensavg(i) + ensem(j,i)
      x(j) = ensem(j,i)
      if (ensem(j,i).gt.pthresh) then
        n = n + 1
        nposens(i) = nposens(i) + 1
        censavg(i) = censavg(i) + ensem(j,i)
        cx(n) = ensem(j,i)
      endif
    enddo
    if (n.gt.1) then
	  sort1.n=n
      sort1.a=>cx
      call sort (sort1)
      j = n/2 + 1
      censmedian(i) = cx(j)
      censavg(i) = censavg(i)/n
    endif
	sort1.n=nmem
    sort1.a=>x
    call sort (sort1)
    j = nmem/2 + 1
    ensmedian(i) = x(j)
    ensavg(i) = ensavg(i)/nmem
    if (obs(i).gt.pthresh) then
      cavgobs = cavgobs + obs(i)
      npos = npos + 1
      xobs(npos) = obs(i)
    endif
  enddo
  avgobs = avgobs/nobs
  if (npos.gt.0) cavgobs = cavgobs/npos
 
  crps_fcst = 0.
  crps_avgobs = 0.
  crps_clim = 0.
  crps_ensavg = 0.
  crps_ens = 0.
  crps_ensmed = 0.
  crps_cfcst = 0.
  crps_cavgobs = 0.
  crps_cclim = 0.
  crps_censavg = 0.
  crps_cens = 0.
  crps_censmed = 0.
  
  nn = 0     
  do i=1,nobs
    !crps1=Scrps(obs(i),1,fcst(i))
    !crpsx1 = crps(crps1)
	crpsx1 = abs(obs(i) - fcst(i))
    crps_fcst = crps_fcst + crpsx1
    !crps1=Scrps(obs(i),1,avgobs)
	crpsx2 = abs(obs(i) - avgobs)
    !crpsx2 = crps(crps1)
    crps_avgobs = crps_avgobs + crpsx2
    crps1=Scrps(obs(i),nobs,obs)
    crpsx3 = crps(crps1)
    crps_clim = crps_clim + crpsx3
    !crps1=Scrps(obs(i),1,ensavg(i))
    !crpsx4 = crps(crps1)
	crpsx4 = abs(obs(i) - ensavg(i))
    crps_ensavg = crps_ensavg + crpsx4
    crps1=Scrps(obs(i),nmem,ensem(:,i))
    crpsx5 = crps(crps1)
    crps_ens = crps_ens + crpsx5
    !crps1=Scrps(obs(i),1,ensmedian(i))
    !crpsx6 = crps(crps1)
	crpsx6 = abs(obs(i) - ensmedian(i))
    crps_ensmed = crps_ensmed + crpsx6
    n = 0
    do j=1,nmem
      if (ensem(j,i).gt.pthresh) then
        n = n + 1
        x(n) = ensem(j,i)
      endif
    enddo
    if (obs(i).gt.pthresh) then
      nn = nn + 1
      n = 0
      do j=1,nmem 
        if (ensem(j,i).gt.pthresh) then
          n = n + 1
          x(n) = ensem(j,i)
        endif
      enddo
      if (n.eq.0) then
        n = 1
        x(n) = 0
      endif
	  crpsx = abs(obs(i) - fcst(i))
      !crpsx = crps(obs(i),1,fcst(i))
      crps_cfcst = crps_cfcst + crpsx
      !crps1=Scrps(obs(i),1,cavgobs)
      !crpsx = crps(crps1)
	  crpsx = abs(obs(i) - cavgobs)
      crps_cavgobs = crps_cavgobs + crpsx
      crps1=Scrps(obs(i),npos,xobs)
      crpsx = crps(crps1)
      crps_cclim = crps_cclim + crpsx
      !crps1=Scrps(obs(i),1,censavg(i))
      !crpsx = crps(crps1)
	  crpsx = abs(obs(i) - censavg(i))
      crps_censavg = crps_censavg + crpsx
      crps1=Scrps(obs(i),n,x)
      crpsx = crps(crps1)
      crps_cens = crps_cens + crpsx
      !crps1=Scrps(obs(i),1,censmedian(i))
      !crpsx = crps(crps1)
	  crpsx = abs(obs(i) - censmedian(i))
      crps_censmed = crps_censmed + crpsx
    endif
  enddo
 
  crps_fcst = crps_fcst/nobs
  crps_avgobs = crps_avgobs/nobs
  crps_clim = crps_clim/nobs
  crps_ensavg = crps_ensavg/nobs
  crps_ens = crps_ens/nobs
  crps_ensmed = crps_ensmed/nobs
   
  if (nn.gt.0) then
    crps_cfcst = crps_cfcst/nn
    crps_cavgobs = crps_cavgobs/nn
    crps_cclim = crps_cclim/nn
    crps_censavg = crps_censavg/nn
    crps_cens = crps_cens/nn
    crps_censmed = crps_censmed/nn
  endif

  crpseval1.crps_avgobs=crps_avgobs
  crpseval1.crps_clim=crps_clim
  crpseval1.crps_fcst=crps_fcst
  crpseval1.crps_ensavg=crps_ensavg
  crpseval1.crps_ens=crps_ens
  crpseval1.crps_ensmed=crps_ensmed
  crpseval1.crps_cavgobs=crps_cavgobs
  crpseval1.crps_cclim=crps_cclim
  crpseval1.crps_cfcst=crps_cfcst
  crpseval1.crps_censavg=crps_censavg
  crpseval1.crps_cens=crps_cens
  crpseval1.crps_censmed=crps_censmed
 

    
  crps_clim = 1. - crps_clim/crps_avgobs
  crps_fcst = 1. - crps_fcst/crps_avgobs
  crps_ensavg = 1. - crps_ensavg/crps_avgobs
  crps_ens = 1. - crps_ens/crps_avgobs
  crps_ensmed = 1. - crps_ensmed/crps_avgobs
 
  crps_cavgobs= 1. - crps_cavgobs/crps_cavgobs
  crps_cclim = 1. - crps_cclim/crps_cavgobs
  crps_cfcst= 1. - crps_cfcst/crps_cavgobs
  crps_censavg = 1. - crps_censavg/crps_cavgobs
  crps_cens = 1. - crps_cens/crps_cavgobs
  crps_censmed = 1. - crps_censmed/crps_cavgobs
           

  deallocate(x)
  deallocate(cx)
  deallocate(xsort)
  deallocate(xensem)
  deallocate(ensmedian)
  deallocate(censmedian)
  deallocate(xobs)
  deallocate(obspos)
  deallocate(ensavg)
  deallocate(censavg)
  deallocate(nposens)
  
  return
end
      

!!********************************!!******************************************************* 
!52 compute continuous rank probability score for a single forecast
!!********************************!!******************************************************* 
function crps(crps1) !(x0,n,xx)
use PrecStru
implicit none 
type(Ssort) ::  sort1
type(Scrps) ::  crps1
 
real(r8)  sum
real x0
integer n,i,ifound
real(r8) crps 
real(r8),pointer :: x(:)
real(r8),pointer :: fx(:)
real(r8),pointer :: xx(:)
real(r8),pointer :: h(:)

  n=crps1.n
  x0=crps1.x0
  xx=>crps1.xx
  allocate(x(n))
  allocate(fx(n))
  allocate(h(n))

  crps = -1.
  if (n.eq.1) then  !  single-value forecast
    crps = abs(x0 - xx(1))
    deallocate(x)
    deallocate(fx)
    deallocate(h)
    return
  elseif (n.gt.1) then
    do i=1,n
      x(i) = xx(i)
      fx(i) = float(i)/(float(n+1))
    enddo
	sort1.n=n
    sort1.a=>x
    call sort (sort1)

    !call sort (n,x)
    ifound = 0
    do i=1,n
      if (x(i).gt.x0.and.ifound.eq.0.) then
        ifound = i
      endif
    enddo
    if (ifound.eq.1) then  !  x0 < x(1)
      do i=1,n
        h(i) = 1.
      enddo
      sum = x(1) - x0
      do i=2,n
        sum = sum + (x(i) - x(i-1))*(((h(i) - fx(i))**2 + (h(i-1) - fx(i-1))**2)/2.)
         
      enddo
    elseif (ifound.eq.0) then  !  x0 > x(n)
      do i=1,n
        h(i) = 0.
      enddo
      sum = x0 - x(n)
      do i=2,n
        sum = sum + (x(i) - x(i-1))*(((h(i) - fx(i))**2 + (h(i-1) - fx(i-1))**2)/2.)
         
      enddo
    else
      do i=1,n
        if (i.lt.ifound) then
          h(i) = 0.
        else
          h(i) = 1.
        endif
      enddo
      sum = 0.
      do i=2,n
        sum = sum + (x(i) - x(i-1))*(((h(i) - fx(i))**2 + (h(i-1) - fx(i-1))**2)/2.)
         
      enddo
    endif
  endif
  crps = sum
  deallocate(x)
  deallocate(fx)
  deallocate(h)
  return
end
      


!!********************************!!******************************************************* 
!51 produce a distribution of precipitation values, pp, given the deterministic forecast, fcst.
!  Input Arguments:
!    
!     fcst    = deterministic forecast
!     fcstopt ifcst= climatological forecast distribution
!     fpop    = climatological probability of fcst > pthreshold
!     fcavg   = climatological mean of fcst when fcst > pthreshold
!     fccv    = climatological cv of fcst when fcst > pthreshold
!     obsopt  iobs= climatological observations distribution
!     obspop  = climatological probability of obs > pthreshold
!     obscavg = climatological mean of obs when obs > pthreshold
!     obsccv  = climatological cv of obs when obs > pthreshold
!     rho     = coefficient of correlation between z(fcst) and z(obs)
!     npp     = number of values of pp (number of ensemble members)
!
!  Output Arguments:
!     genavg  = average value of forecast ensemble distribution (internal)
!     genpop  = pop of forecast ensemble distribution (internal)
!     gencavg = conditional mean of fcst ensemble distribution (internal)
!     genccv  = conditional cv of forecast ensemble distribution (internal)
!     vmax    = maximum value of forecast ensemble distribution (internal)
!     exprobmin = excedence probability associated with vmax
!         note:  The output ensemble distribution is pp and has npp members.  
!                This is computed from a more detailed internal distribution 
!                with nvmax members that has a maximum value of vmax
!     pp      = vector of precipitation values sorted with pp(1) = smallest value
!     pmin    = value climatological cdf probability of observing a value 
!               corresponding to the first (lowest) point on the internal 
!               ensemble forecast distribution
!     pmax    = value climatological cdf probability of observing a value 
!               corresponding to the last (highest) point on the internal 
!               ensemble forecast distribution
!!!********************************!!******************************************************* 
subroutine fcst2ensem(fe1) ! (fcst, &
                               ! fcstopt,pthresh_fcst, &
                               ! fpop,fcavg,fccv, &
                               ! obsopt,pthresh_obs, &
                               ! obspop,obscavg,obsccv, &
                               ! rho,npp, &
                               ! genavg,genpop,gencavg,genccv, &
                               ! vmax,exprobmin,pp,pmin,pmax)
use PrecStru
implicit none 
type(Scgauss) ::  cgauss1
type(Svtrans) :: vtrans1  
type(Sfcst2ensem) :: fe1  
type(Sextractp) :: extractp1  
integer nvmax,npp,nv
parameter (nvmax=1000)
real(r8),pointer :: vval(:)  != values of standard normal deviate v
real(r8),pointer :: cpdfv(:) != probability density function of v
real(r8),pointer :: ccdfv(:) != cumulative distribution function of v
real(r8),pointer :: p(:)     ! input conditional mean variate in interval <j-1,j>
real(r8),pointer :: pp(:)    ! output conditional mean variate in interval i:i!1,npp

!character*8 fcstopt,obsopt
real(r8) f1,ufcst,exprobmin,obscdfmax,xrho,obsccvx,xfcst
real(r8) fcst,pthresh_fcst, fpop,fcavg,fccv, pthresh_obs
real(r8)   obspop,obscavg,obsccv,rho
integer i,ifcst,iobs,npos,ifound
real(r8)  genavg,genpop,gencavg,genccv, vmax,  pmin,pmax
real(r8) popthresh,fpz,obspz,u0 ,gausspi,vtrans,vtransi
real(r8) gaussp,avg,cavg,ccv,pz,estpop,pop

  fcst=fe1.fcst
  ifcst=fe1.ifcst
  pthresh_fcst=fe1.pthresh_fcst
  fpop=fe1.fpop
  fcavg=fe1.fcavg
  fccv=fe1.fccv

  iobs=fe1.iobs
  pthresh_obs=fe1.pthresh_obs
  obspop=fe1.obspop
  obscavg=fe1.obscavg
  obsccv=fe1.obsccv
  rho=fe1.rho
  npp=fe1.npp

  pp=>fe1.pp

  cgauss1.nv=nvmax
  allocate(vval(cgauss1.nv))
  allocate(cpdfv(cgauss1.nv))
  allocate(ccdfv(cgauss1.nv))
  allocate(p(cgauss1.nv))
  !allocate(pp(npp))

  obsccvx = obsccv
  xrho = rho
  if (rho.lt.0.) xrho = 0.
  xfcst = fcst
  if (fcst.lt.pthresh_fcst) then
    xfcst = 0.
  else
    xfcst = fcst - pthresh_fcst
  endif
  genavg = 0.
  genpop = 0.
  gencavg = 0.
  genccv = 0.
  popthresh = .01
  if (fpop.ge.popthresh) then
    fpz = 1. - fpop
  else
    fpz = 1. - popthresh
  endif
  obspz = 1. - obspop
  u0 = gausspi(fpz)
 

  vtrans1.avg=fcavg
  vtrans1.cv=fccv
  vtrans1.pz=fpz
  vtrans1.obs=xfcst
  vtrans1.ivopt=ifcst

  ufcst = vtrans(vtrans1)
!  write (ilog,1001) ufcst,u0
 1001 format (/'subroutine fcst2ensem'/'ufcst = ',d20.10/'u0 =',f10.4)
      
  nv = nvmax
!  vmin = -6.
!  vmax = 6.

  cgauss1.u=ufcst
  cgauss1.u0=u0
  cgauss1.nv=nv
  cgauss1.rho=xrho
 
  cgauss1.vval=>vval
  cgauss1.cpdfv=>cpdfv
  cgauss1.ccdfv=>ccdfv

  call cgauss (cgauss1)
  vmax = vval(nv)
  obscdfmax = gaussp(vmax)
  exprobmin = 1. - obscdfmax
  avg = 0.
  npos = 0
  cavg = 0.
  ccv = 0.
  f1 = 0.
 
  ifound = 0
  pz = 0.
  vtrans1.avg=obscavg
  vtrans1.cv=obsccvx
  vtrans1.pz=obspz
  vtrans1.ivopt=iobs
  do i=1,nv

    vtrans1.obs=vval(i)
    p(i) =max(0.,vtransi(vtrans1))
    if (p(i).gt.0.) p(i) = p(i) + pthresh_obs
    if (p(i).eq.0. .and. pthresh_obs.gt.10) p(i) = p(i) + pthresh_obs+0.1

    avg = avg + p(i)*cpdfv(i)
    estpop = 0.
    if (p(i).eq.0..and.ifound.eq.0) pz = ccdfv(i)
    if (p(i).gt.pthresh_obs) then
      cavg = cavg + p(i)*cpdfv(i)
      ccv = ccv + (p(i)**2)*cpdfv(i)
      ifound = 1
    endif
    f1 = ccdfv(i)
  enddo
  estpop = 1. - pz
  pop = 1. - fpz
  if (estpop.gt.0.) then
    cavg = cavg/estpop
    ccv = ccv/estpop
    ccv = sqrt(abs(ccv - cavg**2))
  else
    cavg = 0.
    ccv = 0.
  endif
  if (cavg.gt.0.) ccv = ccv/cavg
  genavg = avg
  genpop = estpop
  gencavg = cavg
  genccv = ccv
  pmin = p(1)
  pmax = p(nv)

  extractp1.nv=nv
  extractp1.obspz=pz
  extractp1.npp=npp
  extractp1.p=>p
  extractp1.pp=>pp
  extractp1.ccdfv=>ccdfv

  call extractp (extractp1)

  fe1.genavg=genavg
  pp(npp+1)=genavg
  fe1.genpop=genpop
  fe1.gencavg=gencavg
  fe1.genccv=genccv

  fe1.vmax=vmax
  fe1.exprobmin=exprobmin
  fe1.pmin=pmin
  fe1.pmax=pmax


!  write (ilog,*) 'leaving epp3_fcst2ensem_v2'
  deallocate(vval)
  deallocate(cpdfv)
  deallocate(ccdfv)
  deallocate(p)

  return
end



!!********************************!!******************************************************* 
!50 convert std normal deviate, zz, to random variable, xvtransi
!             pz = prob(obs=0.)
!             cavg = mean of the conditional pdf of obs>0.
!             ccv = coefficient of variation of conditional pdf of obs>0.
!
!             ivopt defines the pdf of obs>0. as follows:
!
!                   0 = no transformation:       v = obs
!                   1 = gamma distribution
!                   2 = lognormal distribution
!                   3 = exponential distribution
!                   4 = normal distribution
!                   5 = weibull distribution
!                   6 = beta distribution
!                   7 = uniform distribution
!!!********************************!!******************************************************* 
function vtransi(vtrans1)
! (cavg,ccv,pz,zz,iopt)=(AVG,CV,PZ,OBS,IVOPT) ,pz8,a8,b8,x8
use PrecStru
implicit none 
type(Sbeta) :: beta1  
type(Svtrans) :: vtrans1  
real(r8) betapi2,vtransi
real(r8) cdf,ccdf,gammpi
real(r8) cavg,ccv,pz,zz,z,ccvmax,cavglogp,x
integer iopt
real(r8) zlimit,a,b,xx,f,stdlogp,std
real(r8) gammln,gausspi,gaussp,betai,weibp,weibpi,weibb
  
  cavg=vtrans1.avg
  ccv=vtrans1.cv
  pz=vtrans1.pz
  zz=vtrans1.obs
  iopt=vtrans1.ivopt
  
!     ...check limits...
!
  z = zz
  if (z.lt.-10.) then
    vtransi = -9999.
    return
  endif
  zlimit = 10.
  if (z.lt.-zlimit) z = -zlimit
  if (z.gt.zlimit) z = zlimit
  vtransi = 0.
  if (cavg.le.0.) return
  if (ccv.le.0.02) then
    vtransi = cavg
    return
  endif
  cdf = gaussp(z)
!
!     ...check if variate is below threshold
!
  if (cdf.lt.pz.or.cdf.gt.1.0) then
    vtransi = 0.
    return
  else
!    f = (f-pz)/(1.-pz)
    ccdf = (cdf-pz)/(1.-pz)
  endif
!
!     ...do transformation...
!
  select case (iopt)
  case(0) !    no transformation
    vtransi = z
  case(1) !    gamma distribution
    a = 1./ccv**2
    b = cavg/a
    x = gammpi(a,ccdf)
    vtransi = b*x
  case(5)
    if (ccv.le.0.75) then
      a = 1./ccv**2
      b = cavg/a
      x = gammpi(a,ccdf)
      vtransi = b*x
    else   !    weibull distribution
      b = weibb(ccv)
      a = cavg/dexp(gammln(1. + 1./b))
      vtransi = weibpi (a,b,ccdf)
	end if
  case(2)  !    lognormal distribution
    zz = gausspi(ccdf)
    stdlogp = sqrt(dlog(ccv**2+1.))
    cavglogp = dlog(cavg) - 0.5*stdlogp**2
    x = cavglogp + zz*stdlogp
    vtransi = dexp(x)
  case(3)  !    exponential distribution
    b = 1./cavg
    vtransi = -dlog(1.-ccdf)/b
  case(4)  !   normal distribution
    z = gausspi(ccdf)
    std = ccv*cavg
    vtransi = cavg + z*std
    if (vtransi.lt.0.) vtransi = 0.
  case(6)  !    beta distribution
    if (cavg.gt.1.0) then
      vtransi = -9999.
    else
      a = (1.-cavg*(1.+ccv**2))/ccv**2
      b = ((1.-cavg)/cavg)*a
      f = ccdf
	  beta1.a=a
	  beta1.b=b
	  beta1.f=f
      vtransi = betapi2(beta1)
    endif
  case(7)  !  uniform distribution
    ccvmax = 2./sqrt(12.)
    if (ccv.gt.ccvmax) ccv = ccvmax  !  re-scale to avoid negative values
    a = cavg*(1. - ccv/ccvmax)
    b = cavg*(1. + ccv/ccvmax)
    vtransi = a + ccdf*(b-a)
  case default 
    vtransi = -9999.
  end select
  return
end
 

!!********************************!!******************************************************* 
!49  compute weibull variate corresponding to cumulative probability = f.
!!********************************!!******************************************************* 
function weibpi(a,b,f)
use PrecStru
implicit none 
real(r8)  weibpi,a,b,f
  if (a.lt.0..or.b.lt.0..or.f.le.0..or.f.ge.1.) then
    weibpi = 0.
  else
    weibpi = dexp((dlog(-dlog(1-f))+b*dlog(a))/b)
  endif
  return
end


!!********************************!!******************************************************* 
!48  compute standardized gamma variate corresponding to cumulative probability = f.
!!********************************!!******************************************************* 
function gammpi(a,f)
use PrecStru
implicit none 
real(r8) gammpi,gammp
real(r8) x1,x2,f1,f2,avgf,a,f,df,dfdx,tol,exprob
integer maxit,it
  gammpi=0.d0
  it=0
  maxit=100
  tol=.00000001
  x1=.9*a
  x2=1.1*a
  f1=gammp(a,x1)
  do while (it.le.maxit)
    !10  
    f2=gammp(a,x2)
    df=f2-f1
    avgf=(f1+f2)/2.
    exprob=1.  -  avgf
    if  (abs(avgf).gt.0.)  then
      if  (abs(df/exprob).lt.tol)  then
        gammpi=x2
        return
      endif
    else
      return
    endif
    !if  (it.gt.maxit)  return
    it=it+1
    dfdx=df/(x2-x1)
    x1=x2
    f1=f2
    x2=x2  -  (f2-f)/dfdx
    if  (x2.lt.0.)  x2=0.
    !go  to  10
  end do
  if (it.ge.maxit) gammpi=x2
end



!!********************************!!******************************************************* 
!47 extract the distribution of conditional mean variate pp from the
!     distribution of the conditional mean variate p
!
!     nv       = number of values of input distribution of p(.)
!     p(j)     = input conditional mean variate in interval <j-1,j>
!     ccdfv(j) = cdf of p (cdf at end of interval j)
!     obspz    = probability that p = 0.
!     npp      = number of values of pp(.)
!     pp(i)    = output conditional mean variate in interval i:i=1,npp
!
!     note:  pp and p are sorted with smallest value first
!
!!********************************!!******************************************************* 
subroutine extractp(extractp1) !(nv,p,ccdfv,obspz,npp,pp)     
use PrecStru
implicit none 
type(Sextractp) :: extractp1  

integer i,npp,nv,jlast,jnext
real(r8),pointer :: p(:),pp(:)
real(r8),pointer :: ccdfv(:)
real*8 obspz
logical inside
real  f1,f2,pj,dfv,dff,p1,plast,flast,pnext,fnext,df,avg,sumdf 
  nv=extractp1.nv
  obspz=extractp1.obspz
  npp=extractp1.npp
  p=>extractp1.p
  pp=>extractp1.pp
  ccdfv=>extractp1.ccdfv


  f1 = 0. 
  p1 = 0.
  jlast = 0
  plast = 0.
  flast = 0.
  jnext = 1
  pnext = p(1)
  fnext = ccdfv(1)
  df = 1./float(npp)
  avg = 0.
  sumdf = 0.
  
  do i=1,npp  !  compute value of pp(i)


    f2 = i*df  !  end of probability interval
    if (f2.le.obspz) then
      pp(i) = 0.
    else  !  f2 >= obspz
      inside = .true.
      pj = 0.
      do while (inside)
        jnext = jlast + 1
        if (flast.gt.f2.or.jnext.gt.nv) then
          inside = .false.  !  this should not happen
        else
          if (jnext.le.nv) then  !  jnext is in bounds
            pnext = p(jnext) !  
            fnext = ccdfv(jnext)  !  end of p(jnext) interval
            if (fnext.gt.f2) then  !  p(next) ends after f2
              if (flast.ge.f1) then
                dfv = f2 - flast
              else
                dfv = f2 - f1
              endif
              dff = max(0.,min(dfv,f2-obspz))
              sumdf = sumdf + dfv
              avg = avg + dff*pnext
              if (sumdf.gt.0.) avg = avg/sumdf
              pp(i) = avg
              dfv = fnext - f2
              sumdf = dfv
              dff = max(0.,min(dfv,f2-obspz))
              avg = dff*pnext
              inside = .false.
            else
              if (fnext.le.f1) then  !  p(jnext) is before <f1,f2> interval
                flast = fnext
                plast = pnext
                jlast = jnext
               else       !  p(jnext) is in <f1,f2> interval
                 if (flast.lt.f1) then  !  p(jnext) interval starts before f1
                   if (f1.le.obspz) then
                     if (fnext.le.f2) then  !  p(jnext) ends before f2, update jlast
                      dfv = fnext - f1
                      dff = max(0.,fnext-obspz)
                      sumdf = sumdf + dfv
                      avg = avg + dff*pnext
                      jlast = jnext
                      plast = pnext
                      flast = fnext
                    else       !  p(jnext) ends after f2, don't update jlast
                      dfv = f2 - f1
                      dff = max(0.,f2-obspz)
                      sumdf = sumdf + dfv
                      avg = avg + dff*pnext
                    endif
                  else      !  p(jnext) starts after f1
                    if (fnext.le.f2) then  !  p(jnext) ends before f2, update jlast
                      dfv = fnext - f1
                      sumdf = sumdf + dfv
                      if (fnext.lt.obspz) then
                        dff = 0.
                      else
                        if (f1.lt.obspz) then
                          dff = fnext - obspz
                        else
                          dff = fnext - f1
                        endif
                      endif
                      avg = avg + dff*pnext
                      jlast = jnext
                      plast = pnext
                      flast = fnext
                    else       !  p(jnext) ends after f2, don't update jlast
                      dff = f2 - f1
                      sumdf = sumdf + df
                      avg = avg + dff*pnext
                    endif
                  endif
                else  !  flast is inside of the <f1.f2> interval
                  if (fnext.lt.f2) then  !  all of p(j) is in the interval, update jlast
                    dfv = fnext - flast
                    if (flast.le.obspz) then
                      dff = max(0.,fnext-obspz)
                    else
                      dff = fnext - flast
                    endif
                    sumdf = sumdf + dfv
                    avg = avg + dff*pnext
                    jlast = jnext
                    plast = pnext
                    flast = fnext
                  else  !  fnext is outside of the interval ending at f2, don't update jlast
                    inside = .false.
                    if (flast.ge.f1) then
                      if (flast.lt.obspz) then
                        dff = max(0.,f2-obspz)
                      else
                        dff = f2 - flast
                      endif
                      dfv = f2 - flast
                    else  !  flast < f1
                      if (f1.lt.obspz) then
                        dff = max(0.,f2-obspz)
                      else
                        dff = f2 - f1
                      endif
                      dfv = f2 - f1
                    endif
                    sumdf = sumdf + dfv
                    avg = avg + pnext*dff
                    if (sumdf.gt.0.) then
                      avg = avg/sumdf
                    endif
                    pp(i) = avg
                    sumdf = 0.
                    avg = 0.
                  endif
                endif
              endif
            endif
          endif
        endif  
      enddo
    endif
    f1 = f2  !  set start of next interval
	!if (pp(i)>1000) then
	!	    pp(i)=pp(i)
	!	  end if
  enddo
  return
end

!!********************************!!******************************************************* 
!46 fill daily values from nparint-day time series
!!********************************!!******************************************************* 
subroutine fill (nparint,x,xts)
use PrecStru
implicit none 
integer k,iday,jday,nparint,nint,lday 
real(r8) x(365),xts(365) 
  lday=(365-1)/nparint*nparint+1
  do jday=1,lday-nparint,nparint
    do k=1,nparint
	  iday = jday+k-1       
      xts(iday) = x(jday)+(x(jday+nparint)-x(jday))*(k-1)/nparint
    end do
  end do
  nint=366-lday   
  do iday=lday,365      
    xts(iday) = x(lday)+(x(1)-x(lday))*(iday-lday)/nint     
  end do   
  return
end

!!********************************!!******************************************************* 
!45 analyze joint distribution to estimate parameters
! INPUT ARGUMENTS:  
!   nobs     = number of observations = number of forecasts
!   pthresh1 = precipitation threshold
!   pthresh2 = forecast threshold
!   obstran  = observations probability distribution
!                GAMM = Gamma
!                EXPO = Exponetial
!                LOGN = Log Normal
!                WEIB = Weibull
!   fcsttran = forecast probability distribution
!                GAMM = Gamma
!                EXPO = Exponetial
!                LOGN = Log Normal
!                WEIB = Weibull
!   obsx     = observed values
!   fcstx    = deterministic forecasts
!   verstats = verification statistics switch (compute if 'yes')
!
! OUTPUT ARGUMENTS:
!   avgobs     = average obervation
!   popobs     = climatological probability of obs > pthreshold
!   cavgobs    = climatological mean of obs when obs > pthreshold
!   ccvobs     = climatological cv of obs when obs > pthreshold
!   avgfcst    = average forecast
!   popfcst    = climatological probability of fcst > pthreshold
!   cavgfcst   = climatological mean of fcst when fcst > pthreshold
!   ccvfcst    = climatological cv of fcst when fcst > pthreshold
!   rho        = recommended coefficient of correlation between z(fcst) and z(obs)
!   rmsfcst    = rms difference between fcst and obs
!   effnfcst   = Nash-Sutcliffe efficiency of forecast
!   ratio      = avgobs/avgfcst
!   eqts_fcst  = equitable threat score for wet day forecast
!   cor        = correlation coefficient between obs and fcst (includes zeros)
!   ccor       = correlation coefficient between wet obs and fcsts
!   trccor     = correlation coefficient between transformed wet obs and fcsts
!   rhoopt     = correlation coefficient from transformed forecasts and obs
!!!********************************!!******************************************************* 
subroutine epp_fcstparms(fp1) 

use PrecStru
implicit none 
type(Sbivar) :: bivar1  
type(Svtrans) :: vtrans1  
type(Sfcstparm) :: fp1  
integer i,minobs,j,k,nccor,nobspx,nfcstpx
integer ndp,np
real(r8) a,b,c,d,stdobs,stdfcst,cstdobs,exp_h
real(r8) cstdfcst,ccavgobs,ccstdobs,ccavgfcst,ccstdfcst,qscore
real(r8) rho,eff,ratio,cvobs,cvfcst,rho22,bias,eqts_fcst
real(r8) trcstd,trcavg ,trcavgobs,trfcst,trobs,vtrans
real(r8) pzfcst,trcstdobs,travg,travgobs,pzobs
real(r8),pointer :: obs(:)
real(r8),pointer :: fcst(:)
real(r8),pointer :: obsxx(:)
real(r8),pointer :: fcstxx(:)
real(r8),pointer :: ensfcstx(:)
real(r8) pval(2,2)
integer*4 nvalval(2,2)
real(r8) avgobsval(2,2),stdobsval(2,2)
real(r8) avgensval(2,2),stdensval(2,2)

  allocate(obs(fp1.nobs))
  allocate(fcst(fp1.nobs))
  allocate(obsxx(fp1.nobs))
  allocate(fcstxx(fp1.nobs))
  allocate(ensfcstx(fp1.nobs))

  minobs = 20  
  fp1.avgobs = 0.
  stdobs = 0.
  fp1.avgfcst = 0.
  stdfcst = 0.
  fp1.cor = 0.
  nobspx = 0
  fp1.cavgobs = 0.
  cstdobs = 0.
  nfcstpx = 0
  fp1.cavgfcst = 0.
  cstdfcst = 0.
  fp1.ccor = 0.
  ccavgobs = 0.
  ccstdobs = 0.
  ccavgfcst = 0.
  ccstdfcst = 0.
  nccor = 0
  fp1.rmsfcst = 0.
  qscore = 0.
  rho = 0.
  eff = 0.
  ratio = 0.

  do i=1,fp1.nobs
    if (fp1.obsx(i).gt.fp1.pthresh1) then
      obs(i) = fp1.obsx(i) - fp1.pthresh1
    else
      obs(i) = 0.
    endif
    if (fp1.fcstx(i).gt.fp1.pthresh2) then
      fcst(i) = fp1.fcstx(i) - fp1.pthresh2
    else
      fcst(i) = 0.
    endif
  enddo

  if (fp1.nobs.ge.minobs) then
    do i=1,fp1.nobs
      fp1.avgobs = fp1.avgobs + fp1.obsx(i)
      stdobs = stdobs + fp1.obsx(i)**2
      fp1.avgfcst = fp1.avgfcst + fp1.fcstx(i)
      stdfcst = stdfcst + fp1.fcstx(i)**2
      fp1.rmsfcst = fp1.rmsfcst + (fp1.obsx(i)-fp1.fcstx(i))**2
      fp1.cor = fp1.cor + fp1.obsx(i)*fp1.fcstx(i)
      if (obs(i).gt.0.) then
        nobspx = nobspx + 1
        fp1.cavgobs = fp1.cavgobs + obs(i)
        cstdobs = cstdobs + obs(i)**2
      endif
      if (fcst(i).gt.0.) then
        nfcstpx = nfcstpx + 1
        fp1.cavgfcst = fp1.cavgfcst + fcst(i)
        cstdfcst = cstdfcst + fcst(i)**2
      endif
      if (obs(i).gt.0..and.fcst(i).gt.0.) then
        nccor = nccor + 1
        ccavgobs = ccavgobs + obs(i)
        ccstdobs = ccstdobs + obs(i)**2
        ccavgfcst = ccavgfcst + fcst(i)
        ccstdfcst = ccstdfcst + fcst(i)**2
        fp1.ccor = fp1.ccor + obs(i)*fcst(i)
      endif
    enddo
    fp1.avgobs = fp1.avgobs/fp1.nobs
    stdobs = sqrt(abs(stdobs/fp1.nobs - fp1.avgobs**2))
    cvobs = stdobs/fp1.avgobs
    fp1.avgfcst = fp1.avgfcst/fp1.nobs
    stdfcst = sqrt(abs(stdfcst/fp1.nobs-fp1.avgfcst**2))
    cvfcst = stdfcst/fp1.avgfcst
    fp1.rmsfcst = sqrt(fp1.rmsfcst/fp1.nobs)
    bias = (fp1.avgfcst - fp1.avgobs)/((fp1.avgfcst+fp1.avgobs)/2.)
    fp1.effnsfcst = 1. - fp1.rmsfcst**2/stdobs**2
    if (stdfcst.gt.0..and.stdobs.gt.0.) then
      fp1.cor = (fp1.cor/fp1.nobs-fp1.avgobs*fp1.avgfcst)/(stdfcst*stdobs)
      fp1.ratio = fp1.avgobs/fp1.avgfcst
    else
      fp1.cor = 0.    
      fp1.ratio = 0.
    endif

    if (nobspx.gt.0) then
      fp1.popobs = float(nobspx)/float(fp1.nobs)
      if (nobspx.gt.20) then
        fp1.cavgobs = fp1.cavgobs/nobspx
        cstdobs = sqrt(abs(cstdobs/nobspx-fp1.cavgobs**2))
        fp1.ccvobs = cstdobs/fp1.cavgobs
        fp1.cavgobsx = fp1.cavgobs
      else
        fp1.cavgobs = fp1.cavgobsx
        cstdobs = fp1.cavgobsx
        fp1.ccvobs = 1.
      endif
!    if (fp1.cavgobs.gt.0..and.fp1.cavgobs.ge.avgobs) then
!      popobs = avgobs/fp1.cavgobs
!    endif
    else
      fp1.popobs = 0.
      fp1.cavgobs = fp1.cavgobsx
      cstdobs = fp1.cavgobsx
      fp1.ccvobs = 1.
    endif

    if (nfcstpx.gt.0) then
      fp1.popfcst = float(nfcstpx)/float(fp1.nobs)
      if (nfcstpx.gt.20) then
        fp1.cavgfcst = fp1.cavgfcst/nfcstpx
        cstdfcst = sqrt(abs(cstdfcst/nfcstpx-fp1.cavgfcst**2))
        fp1.ccvfcst = cstdfcst/fp1.cavgfcst
        fp1.cavgfcstx = fp1.cavgfcst
      else
        fp1.cavgfcst = fp1.cavgfcstx
        cstdfcst = fp1.cavgfcstx
        fp1.ccvfcst = 1.
      endif
    else
      fp1.cavgfcst = fp1.cavgfcstx
      cstdfcst = fp1.cavgfcstx
      fp1.ccvfcst = 1.
      fp1.popfcst = 0.
    endif  

    if (nccor.gt.0) then
      ccavgfcst = ccavgfcst/nccor
      ccstdfcst = sqrt(abs(ccstdfcst/nccor-ccavgfcst**2))
      ccavgobs = ccavgobs/nccor
      ccstdobs = sqrt(abs(ccstdobs/nccor-ccavgobs**2))
	  if (abs(ccstdobs*ccstdfcst).lt.0.0001) then
        fp1.ccor =1
      else
        fp1.ccor = (fp1.ccor/nccor - ccavgobs*ccavgfcst)/(ccstdobs*ccstdfcst) 
      end if   
    endif
!
!   ...compute contingency table...
!
    do i=1,2
      do j=1,2
        fp1.nvalcal(i,j) = 0
        fp1.avgobscal(i,j) = 0.
        fp1.stdobscal(i,j) = 0.
        fp1.avgfcstcal(i,j) = 0.
        fp1.stdfcstcal(i,j) = 0.
      enddo
    enddo
    do k=1,fp1.nobs
      if (fp1.fcstx(k).gt.fp1.pthresh2) then
        i = 2
      else
        i = 1
      endif
      if (fp1.obsx(k).gt.fp1.pthresh1) then
        j = 2
      else
        j = 1
      endif
      fp1.nvalcal(i,j) = fp1.nvalcal(i,j) + 1
      fp1.avgobscal(i,j) = fp1.avgobscal(i,j) + fp1.obsx(k)
      fp1.stdobscal(i,j) = fp1.stdobscal(i,j) + fp1.obsx(k)**2
      fp1.avgfcstcal(i,j) = fp1.avgfcstcal(i,j) + fp1.fcstx(k)
      fp1.stdfcstcal(i,j) = fp1.stdfcstcal(i,j) + fp1.fcstx(k)**2
    enddo
    do i=1,2
      do j=1,2
        if (fp1.nvalcal(i,j).gt.0.) then
          fp1.avgobscal(i,j) = fp1.avgobscal(i,j)/fp1.nvalcal(i,j)
          fp1.stdobscal(i,j) = sqrt(abs(fp1.stdobscal(i,j)/fp1.nvalcal(i,j) -fp1.avgobscal(i,j)**2))
          fp1.avgfcstcal(i,j) = fp1.avgfcstcal(i,j)/fp1.nvalcal(i,j)
          fp1.stdfcstcal(i,j) = sqrt(abs(fp1.stdfcstcal(i,j)/fp1.nvalcal(i,j) -fp1.avgfcstcal(i,j)**2))
                  
        endif
      enddo
    enddo
    ndp = 0
    do i=1,2
      do j=1,2
        ndp = ndp + fp1.nvalcal(i,j)
      enddo
    enddo
    do i=1,2
      do j=1,2
        pval(i,j) = float(fp1.nvalcal(i,j))/float(ndp)
      enddo
    enddo
    a = fp1.nvalcal(2,2)
    b = fp1.nvalcal(2,1)
    c = fp1.nvalcal(1,2)
    d = fp1.nvalcal(1,1)
    exp_h = ((a+b)*(a+c))/(a+b+c+d)
    eqts_fcst = (a - exp_h)/(a + b + c - exp_h)
        
    !   ...transform obs (Weibull) and ens mean (gamma distr) to snd
 
    fp1.trccor = 0.
    trcavg = 0.
    trcstd = 0.
    trcavgobs = 0.
    trcstdobs = 0.
    travg = 0.
    travgobs = 0.
    np = 0  

    if (fp1.cor.lt.0.9995) then
      pzobs = 1. - fp1.popobs
      pzfcst = 1. - fp1.popfcst
      do i=1,fp1.nobs
        if (obs(i).gt.0.) then
          if (fcst(i).gt.0.) then
		    vtrans1=Svtrans(fp1.cavgobs,fp1.ccvobs,pzobs,obs(i),fp1.iobstran)
            trobs = vtrans(vtrans1)  
            vtrans1=Svtrans(fp1.cavgfcst,fp1.ccvfcst,pzfcst,fcst(i),fp1.ifcsttran)
            trfcst = vtrans(vtrans1)
            
            trcavg = trcavg + trfcst
            trcstd = trcstd + trfcst**2
            trcavgobs = trcavgobs + trobs
            trcstdobs = trcstdobs + trobs**2
            fp1.trccor = fp1.trccor + trfcst*trobs
            np = np + 1
          endif
        endif
      enddo
      if (np.gt.minobs) then
        trcavg = trcavg/np
        trcstd = sqrt(abs(trcstd/np - trcavg**2))
        trcavgobs = trcavgobs/np
        trcstdobs = sqrt(abs(trcstdobs/np - trcavgobs**2))
        fp1.trccor = (fp1.trccor/np - trcavg*trcavgobs)/(trcstd*trcstdobs)        
	  else
        fp1.trccor = fp1.cor
      endif
      rho22=fp1.trccor  
      if (pzobs.lt.0.10.or. &
        pzobs.gt.0.90.or. &
        pzfcst.lt.0.10.or. &
        pzfcst.gt.0.90) then
        fp1.rhoopt = fp1.cor
      else
        if (pval(1,1).lt.0.05.or. &
           pval(1,1).gt.0.95.or. &
           pval(2,2).lt.0.05.or. &
           pval(2,2).gt.0.95) then
           fp1.rhoopt = fp1.cor
        else
			!real(r8) ax,bx,cx,tol,rho22,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,rho
		  bivar1.p11=pval(1,1)
		  bivar1.p12=pval(1,2)
		  bivar1.p21=pval(2,1)
		  bivar1.p22=pval(2,2)
		  bivar1.rho22=rho22
		  bivar1.rho=fp1.rhoopt
          call estrho3(bivar1)  
          fp1.rhoopt=bivar1.rho        
        endif
      endif
    else
      fp1.ccor = fp1.cor
      fp1.trccor = fp1.cor
      fp1.rhoopt = fp1.cor
    endif
  else
    fp1.cavgobs = fp1.cavgobsx
    fp1.cavgfcst = fp1.cavgfcstx
  endif

  !1 cor        ! correlation coefficient between obs and fcst (includes zeros)
  !2 rhoopt     ! correlation coefficient from transformed forecasts and obs
  !3 w*cor + (1. - w)*rhoopt
  !4 ccor       ! correlation coefficient between wet obs and fcsts (>0)
  !5 trccor     ! correlation coefficient between transformed wet obs and fcsts (>0)
  select case (fp1.iopt_rho)
  case(1) 
    fp1.rho_est = fp1.cor
  case(2) 
    fp1.rho_est = fp1.rhoopt
  case(3) 
    fp1.rho_est = fp1.cor_weight*fp1.cor + (1. - fp1.cor_weight)*fp1.rhoopt
  case(4) 
    fp1.rho_est = fp1.ccor
  case(5) 
    fp1.rho_est = fp1.trccor
  case default  
    fp1.rho_est = fp1.cor
  end select
   
  deallocate(obs)  
  deallocate(fcst)  
  deallocate(obsxx)  
  deallocate(fcstxx)  
  deallocate(ensfcstx)  

  return
end


!!********************************!!******************************************************* 
!44 estimate correlation coefficient of bivariate normal distribution 
!     consistent with the following parameters:
!     u0 = probability that u <= 0 (input std normal deviate)
!     v0 =      "       "   v <= 0 (output "    "       "   )
!     p11 = prob v<=0 if u<=0
!     p12 = prob v>0 if u<=0
!     p21 = prob v<=0 if u>0
!     p22 = prob v>0 if u>0
!!********************************!!******************************************************* 
subroutine estrho3(bivar1)
use PrecStru
implicit none 
type(Sbivar) :: bivar1  
real(r8) pu0,pv0,findrho,f1,GAUSSPI
!real(r8) minfunc, minrho
!real(r8) ax,bx,cx,tol,rho22,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,rho
!p11,p12,p21,p22,rho22,rho 
!  write (iout,'(2f6.3)') u0,v0
!  write (iout,*)
!  write (iout,'(2f6.3)') p21,p22
!  write (iout,'(2f6.3)') p11,p12

  pu0 = bivar1.p11+bivar1.p12
  pv0 = bivar1.p11+bivar1.p21
  bivar1.u0 = GAUSSPI(pu0)
  bivar1.v0 = GAUSSPI(pv0)
!  
! write (iout,*)
! write (iout,'(2f6.3)') pu0,pv0
! write (iout,'(2f6.3)') u0,v0
! write (iout,*)
!
  bivar1.umin = -5.
  bivar1.umax = 5.
  bivar1.vmin = -5.
  bivar1.vmax = 5.

  bivar1.ax = 0.02
  bivar1.bx = bivar1.rho22
  bivar1.cx = 0.98
  bivar1.tol = 0.001
  !bivar1=Sbivar(ax,bx,cx,tol,rho22,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,rho)
  f1 = findrho(bivar1)
          
!      AX,BX,CX,TOL,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,rho
!     write (iout,*) 'rho ='
!     write (iout,'(5x,f7.3)') rho
!
  return
end
!!********************************!!******************************************************* 
!43 find  correlation coefficient of bivariate normal distribution
!!********************************!!******************************************************* 
function findrho(bivar1)
use PrecStru
implicit none 
type(Sbivar) :: bivar1  
type(Sbivar) :: bivar2 !ax,bx,cx,tol,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,xmin
real(r8) findrho,x0,x1,x2,x3,f0,f1,f2,f3
real(r8) r,c,frho
integer it,maxit
  bivar2=bivar1   
           
  r=.61803399
  c=.38196602
 				          
  maxit = 5
  it = 0
  x0=bivar1.ax
  x3=bivar1.cx
  if(abs(bivar1.cx-bivar1.bx).gt.abs(bivar1.bx-bivar1.ax))then
    x1=bivar1.bx
    x2=bivar1.bx+c*(bivar1.cx-bivar1.bx)
  else
    x2=bivar1.bx
    x1=bivar1.bx-c*(bivar1.bx-bivar1.ax)
  endif
  bivar2.rho=x1
  f1=frho(bivar2)!(p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,x1)
  bivar2.rho=x2
  f2=frho(bivar2)!(p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,x2)
  do while (it.lt.maxit .and. abs(x3-x0).gt.bivar1.tol*(abs(x1)+abs(x2)))
    it = it + 1
    !if (it.lt.maxit) then
    !if(abs(x3-x0).gt.bivar1.tol*(abs(x1)+abs(x2)))then
    if(f2.lt.f1)then
      x0=x1
      x1=x2
      x2=r*x1+c*x3
      if (x2.lt.bivar1.ax) then
        x2 = bivar1.ax
      endif
      if (x2.gt.bivar1.cx) then
        x2 = bivar1.cx
      endif
      f0=f1
      f1=f2
      bivar2.rho=x2
      f2=frho(bivar2)!(p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,x2)
    else
      x3=x2
      x2=x1
      x1=r*x2+c*x0
      f3=f2
      f2=f1
      if (x1.lt.bivar1.ax) then
        x1 = bivar1.ax
      endif
      if (x1.gt.bivar1.cx) then
        x1 = bivar1.cx
      endif
      bivar2.rho=x1
      f1=frho(bivar2)!(p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,x1)
    endif
  end do
  if(f1.lt.f2)then
    findrho=f1
    bivar1.rho=x1
  else
    findrho=f2
    bivar1.rho=x2
  endif
  return
end

!!********************************!!******************************************************* 
!42 evaluate binormal fit to a contingency table...
!!********************************!!******************************************************* 
function frho(bivar1)
use PrecStru
implicit none 
type(Sbivar) :: bivar1  
type(Sbivar) :: bivar2 
real(r8) frho !,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,rho
real(r8) p11hat,p12hat,p21hat,p22hat,f1,binormp2
  !bivar1=Sbivar(p11,p12,p21,p22,umin,u0,u0,vmin,v0,v0,rho) ! (umin,u0,vmin,v0,rho)
  bivar2=bivar1
  bivar2.umax=bivar1.u0  
  bivar2.vmax=bivar1.v0  
  p11hat = binormp2 (bivar2)
  !bivar1=Sbivar(p11,p12,p21,p22,umin,u0,u0,v0,v0,vmax,rho) !(umin,u0,v0,vmax,rho)
  bivar2=bivar1
  bivar2.umax=bivar1.u0  
  bivar2.vmin=bivar1.v0  
  p12hat = binormp2 (bivar2)
  !bivar1=Sbivar(p11,p12,p21,p22,u0,u0,umax,vmin,v0,v0,rho) !(u0,umax,vmin,v0,rho)
  bivar2=bivar1
  bivar2.umin=bivar1.u0  
  bivar2.vmax=bivar1.v0  
  p21hat = binormp2 (bivar2)
  !bivar1=Sbivar(p11,p12,p21,p22,u0,u0,umax,v0,v0,vmax,rho) !(u0,umax,v0,vmax,rho)  
  bivar2=bivar1
  bivar2.umin=bivar1.u0  
  bivar2.vmin=bivar1.v0  
  p22hat = binormp2 (bivar2)  
  f1=abs(bivar1.p11-p11hat)+abs(bivar1.p12-p12hat)+abs(bivar1.p21-p21hat)+abs(bivar1.p22-p22hat)
  frho = f1
  return
end

!!********************************!!******************************************************* 
!41 function to integrate bivariate normal density function
!!********************************!!******************************************************* 
function binormp2 (bivar1)
use PrecStru
implicit none 
type(Sbivar) :: bivar1  
real(r8) sum,pi
real(r8) binormp2,umin,umax,vmin,vmax,rho
real(r8) dv,evar,std,v,uavg,z1,z2,gaussp,probu,probuv,gausspdf
integer nv,i
real(r8) uminx,umaxx,vminx,vmaxx

  umin=bivar1.umin
  umax=bivar1.umax
  vmin=bivar1.vmin
  vmax=bivar1.vmax
  rho=bivar1.rho

  uminx = max(umin,-5.)
  umaxx = min(umax,5.)
  vminx = max(vmin,-5.)
  vmaxx = min(vmax,5.)
  nv = 1000
  dv = (vmaxx - vminx)/nv
  if (abs(rho).lt.1.0) then
    evar = 1. - rho**2
    if (evar.lt.0.) then
      evar = 1.0
    endif
    if (evar.gt.1.0) then
      evar = 1.0 
    endif  
    std = sqrt(abs(evar))
    sum = 0.0
  else
    std = 0.001
  endif
  do i=1,nv
    v = vminx + (i-.5)*dv
    uavg = rho*v
    z1 = (uminx-uavg)/std
    z2 = (umaxx-uavg)/std
    probu = gaussp(z2) - gaussp(z1)
    probuv = probu*gausspdf(v)*dv
    sum = sum + probuv
  enddo
  binormp2 = sum
  return
end

!!********************************!!******************************************************* 
!40 compute standard normal deviate corresponding to OBS 
!   when OBS is from a distribution having the following properties:
!             PZ=Prob(OBS=0.)
!             AVG=mean of the conditional PDF of OBS>0.
!             CV=coefficient of variation of conditional PDF of OBS>0.
!
!             IVOPT defines the PDF of OBS>0. as follows:
!
!                   0=no transformation:       v=obs
!                   1=gamma distribution
!                   2=lognormal distribution
!                   3=exponential distribution
!                   4=normal distribution
!                   5=weibull distribution
!                   6=beta distribution
!                   7=uniform distribution
!
!       NOTES:
!             If OBS<0, the transformation is not defined and
!             VTRANS is set equal to the missing data value -9999.
!
!             Values of VTRANS satisfy:  -3.0 < VTRANS < 3.0
!
!********************************!!******************************************************* 
function vtrans (vtrans1 )
use PrecStru
implicit none 
type(Sbeta) :: beta1  
type(Svtrans) :: vtrans1  

real(r8) vtrans,avg,cv,pz,obs
real(r8) zlimit,theta,cvt,avgt,a,b,xx,f,stdlogp,avglogp
real(r8) cvmax,ccv,z0,std,z
real(r8) gammp,gammln,weibb,gausspi,gaussp,betai,weibp
integer ivoptmax,ivopt
  avg=vtrans1.avg
  cv=vtrans1.cv
  pz=vtrans1.pz
  obs=vtrans1.obs
  ivopt=vtrans1.ivopt  

  zlimit=5.
  theta=0.25
  ivoptmax=6
!
!  ...check  input  arguments...
!
  if  (ivopt.lt.0.or.ivopt.gt.ivoptmax)  then
    vtrans=-9999.
  endif
  if  (obs.lt.0.)  then
    vtrans=-9999.
    return
  endif
  if  (cv.lt.0.02)  then
    cvt=.02
  else
    cvt=cv
  endif
  if  (avg.lt.0.01)  then
    avgt=0.01
  else
    avgt=avg
  endif
  !
  !  ...do  transformation...
  !
  if  (obs.eq.0.)  then
    if  (pz.gt.0)  then    !    set  vtrans=expected  value  for  obs=0.
      z0=gausspi(pz)
      vtrans=theta*z0  +(1.-theta)*(-0.3984228*exp(-z0*z0/2.)/pz)        
    else
      if  (avg.eq.0.)  then
        vtrans=0.
      else
        vtrans=-zlimit
      endif
    endif
    return
  endif

  !if  (obs.gt.0.)  then
  select case (ivopt)

  case(0) !    no  transformation
  !if  (ivopt.eq.0)  then      
    vtrans=obs
    return
  !endif
  case(1) ! Gamma distribution
  !if  (ivopt.eq.1)  then
  ! a=1./sqrt(cvt)
    a=1./cvt**2
    b=avgt/a
    xx=obs/b
    f=gammp(a,xx)
    f=pz  +  (1.-pz)*f 
  case(2) !Log Normal distribution
  !elseif  (ivopt.eq.2)  then
    stdlogp=sqrt(dlog(cvt**2+1.))
    avglogp=dlog(avgt)  -  0.5*stdlogp**2
    xx=dlog(obs)
    z=(xx-avglogp)/stdlogp
    f=gaussp(z)
    f=pz  +  (1.-pz)*f 
  case(3) !Exponetial distribution
  !elseif  (ivopt.eq.3)  then
    b=1./avgt
    f=1.  -  dexp(-b*obs)
    f=pz  +  (1.-pz)*f
  case(4) 
  !elseif  (ivopt.eq.4)  then    !    normal  distribution
    std=avg*cvt
    z=(obs-avg)/std
    f=gaussp(z)
    f=pz  +  (1.-pz)*f
  case(5) 
    if (cvt.le.0.75) then  !  gamma distribution
      a = 1./cvt**2
      b=avgt/a
      xx=obs/b
      f=gammp(a,xx)
	  f=pz  +  (1.-pz)*f  
	else
      !elseif  (ivopt.eq.5)  then      !      weibull  distribution
      b=weibb(cvt)
      a=avg/dexp(gammln(1.  +  1./b))
      f=weibp(a,b,obs)
      f=pz  +  (1.-pz)*f
	end if
  case(6) 
  !elseif  (ivopt.eq.6)  then      !      beta  distribution
    if  (avg.gt.1.0.or.obs.gt.1.0)  then
      vtrans=-9999.
      return
    else
      a=(1.-avg*(1.+cv**2))/cv**2
      b=((1.-avg)/avg)*a
      beta1.a=a
      beta1.b=b
      beta1.f=obs
      f=betai(beta1)
    endif
    f=pz  +  (1.-pz)*f
  case(7) 
  !elseif  (ivopt.eq.7)  then    !    uniform  distribution
    cvmax=2./sqrt(12.)
    ccv=cvt
    if  (ccv.gt.cvmax)  ccv=cvmax    !    rescale  to  avoid  negative  values
    a=avg*(1.  -  ccv/cvmax)
    b=avg*(1.  +  ccv/cvmax)
    if  (obs.lt.a)  then
      f=0.
    elseif  (obs.gt.b)  then
      f=1.0
    else
      f=pz  +  (1.-pz)*(obs-a)/(b+a)  
    endif

  case default  
    f=0.5   
  end select
  vtrans=gausspi(f)      !  compute  standard  normal  deviate
!
!  ...check  limits  on  vtrans...
!
  if  (vtrans.gt.zlimit)  vtrans=zlimit
  if  (vtrans.lt.-zlimit)  vtrans=-zlimit
  return
end

!!********************************!!******************************************************* 
!39 compute beta variate corresponding to cumulative probability=f.
!!********************************!!******************************************************* 
function betapi2(beta1)
use PrecStru
implicit none 
type(Sbeta) :: beta1  
type(Sbeta) :: beta2  
real(r8)  betapi2,a,b,f,betai
real(r8) r,c,ax,bx,cx,fa,fb,fc,ea,eb,ec,tol
real(r8) betap12,xmin,x0,x1,x2,x3,golden,f0,f1,f2,f3
integer it,maxit
  a=beta1.a
  b=beta1.b
  f=beta1.f

  r=.61803399
  c=.38196602
 
  maxit=100
  tol=0.0001
  it=0
  ax=0.
  bx=0.50
  cx=1.0
  beta2.a=a
  beta2.b=b
  beta2.f=ax
  fa=betai(beta2)
  beta2.f=bx
  fb=betai(beta2)
  beta2.f=cx
  fc=betai(beta2)
  ea=abs(fa-f)
  eb=abs(fb-f)
  ec=abs(fc-f)
  do while  (eb.gt.ea.or.eb.gt.ec)  
    it=it + 1
    if (it.gt.maxit) then
      betapi2=bx
      return
    endif
    if (ea.lt.eb) then
      bx=(ax + bx)/2.
    else
      bx=(bx + cx)/2.
    endif
    beta2.a=a
    beta2.b=b
    beta2.f=bx
    fb=betai(beta2)
    eb=abs(fb-f)
  end do
  
  if (eb.lt.ea.and.eb.lt.ec) then
    it=0
    x0=ax
    x3=cx
    if(abs(cx-bx).gt.abs(bx-ax))then
      x1=bx
      x2=bx+c*(cx-bx)
    else
      x2=bx
      x1=bx-c*(bx-ax)
    endif
    beta2.a=a
    beta2.b=b
    beta2.f=x1
    f1=abs(betai(beta2)-f)
    beta2.f=x2
    f2=abs(betai(beta2)-f)
    do while (it.lt.maxit.and.abs(x3-x0).gt.0.00001)
      it=it + 1
      if(f2.lt.f1)then
        x0=x1
        x1=x2
        x2=r*x1+c*x3
        f0=f1
        f1=f2
        beta2.a=a
        beta2.b=b
        beta2.f=x2
        f2=abs(betai(beta2)-f)
      else
        x3=x2
        x2=x1
        x1=r*x2+c*x0
        f3=f2
        f2=f1
        beta2.a=a
        beta2.b=b
        beta2.f=x1
        f1=(betai(beta2)-f)
        f1=abs(f1)
      endif
    end do 
    if(f1.lt.f2)then
      golden=f1
      xmin=x1
    else
      golden=f2
      xmin=x2
    endif
    betapi2=xmin
  else
    betap12=bx
    !  stop 'error in betapi2'
  endif
  return
end


!!********************************!!******************************************************* 
!38    compute beta distribution
!!********************************!!******************************************************* 

function betai(beta1)
use PrecStru
implicit none 
type(Sbeta) :: beta1
type(Sbeta) :: beta2
 
real(r8) betai,betacf,a,b,x,bt
real(r8) gammln
  a=beta1.a
  b=beta1.b
  x=beta1.f
  if (a.le.0..or.b.le.0.) then
    betai=999
    return
  end if
!     if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
      if (x.lt.0..or.x.gt.1.) then
        betai=0.
        return
      endif
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=dexp(gammln(a+b)-gammln(a)-gammln(b)+a*dlog(x)+b*dlog(1.-x))
          
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
	    beta2.a=a
	    beta2.b=b
	    beta2.f=x
        betai=bt*betacf(beta2)/a
        return
      else
	    beta2.a=b
	    beta2.b=a
	    beta2.f=1.-x
        betai=1.-bt*betacf(beta2)/b
        return
      endif
end


!!********************************!!******************************************************* 
!37 compute beta distribution
!!********************************!!******************************************************* 

function betacf(beta1)
use PrecStru
implicit none 
type(Sbeta) :: beta1
real(r8) betacf,a,b,x
integer itmax,m,em,tem
real(r8) eps,am,bm,az,qab,qap,qam,bz
real(r8) d,app,bpp,ap,bp,aold
  a=beta1.a
  b=beta1.b
  x=beta1.f

  itmax=100
  eps=3.e-7
  am=1.
  bm=1.
  az=1.
  aold=0.
  qab=a+b
  qap=a+1.
  qam=a-1.
  bz=1.-qab*x/qap
  m=1  
  do while (m.le.itmax.and.abs(az-aold).ge.eps*abs(az))
    em=m
    tem=em+em
    d=em*(b-m)*x/((qam+tem)*(a+tem))
    ap=az+d*am
    bp=bz+d*bm
    d=-(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
    app=ap+d*az
    bpp=bp+d*bz
    aold=az
    am=ap/bpp
    bm=bp/bpp
    az=app/bpp
    bz=1.
    !if(abs(az-aold).lt.eps*abs(az))go to 1
    m=m+1
  end do 
  betacf=az
  beta1.beta=betacf
  return
end
!!********************************!!******************************************************* 
!36   compute conditional distribution of bivariate std normal
!     variable, v, given the other std normal value is equal to u.
!     (u,u0,rho,nv,vval,cpdfv,ccdfv)
!     u     = given value of correlated standard normal deviate
!     u0    = value of u corresponding to g(u) = 0, u<=u0
!     rho   = coefficient of correlation between u and v
!     nv    = number of values of vval
!     vmin  = minimum value of vval range
!     vmax  = maximum value of vval range
!     vval  = values of standard normal deviate v
!     cpdfv = probability density function of v
!     ccdfv = cumulative distribution function of v
!!********************************!!******************************************************* 
subroutine cgauss(cgauss1)
use PrecStru
implicit none 
type(Scgauss) ::  cgauss1
real(r8),pointer :: vval(:),cpdfv(:),ccdfv(:)
integer nv,i
real(r8) u,u0,rho  !,vmin,vmax
real(r8) dv,std,sum,vavg,avgzv,avgval
real(r8) z0,vv,uavg,fpz,fpop,zv,ustd
real(r8) gausspdf,gaussp
real(r8) dz,z1,z2,nv2,vstd,vstdthresh,theta,rhothresh 

  u=cgauss1.u
  u0=cgauss1.u0
  nv=cgauss1.nv
  rho=cgauss1.rho
 
  vval=>cgauss1.vval
  cpdfv=>cgauss1.cpdfv
  ccdfv=>cgauss1.ccdfv

  rhothresh = 0.9999d0
  vstdthresh = sqrt(1. - rhothresh**2)
  theta = 0.25d0
  z1 = -4.d0
  z2 = 4.d0
  dz = (z2 - z1)/float(nv)
  if (rho.lt.0.) then
    vavg = 0.d0
    vstd = 1.d0
  elseif (rho.lt.rhothresh) then
    vavg = rho*u  !  conditional mean of v given u
    vstd = sqrt(1. - rho**2)  !  conditional std dev of v given u
  else
    vavg = rhothresh*u
    vstd = vstdthresh
  endif
  fpz = gaussp(u0)
  fpop = 1. - fpz
  nv2 = nv/2
  if (rho.le.rhothresh) then
    dv = vstd*dz
  else
    dv = 0.
  endif
  do i=1,nv
    vval(i) = vavg + (i-1-nv2)*dv
  enddo
  sum = 0.
  if (u.gt.u0) then  !  v is normal conditional dististribution, given u
    if (rho.lt.rhothresh) then
      do i=1,nv
        zv = (vval(i) - vavg)/vstd  !  standardized conditional v
        cpdfv(i) = gausspdf(zv)*dv
        sum = sum + cpdfv(i)
        ccdfv(i) = sum
      enddo
    else
      do i=1,nv
        cpdfv(i) = 1./float(nv)
        sum = sum + cpdfv(i)
        ccdfv(i) = sum
      enddo
    endif
  else  !  conditional distribution of v is integrated over u<=u0
    if (rho.gt.rhothresh) then
      do i=1,nv
        cpdfv(i) = 1./float(nv)
        sum = sum + cpdfv(i)
        ccdfv(i) = sum
      enddo
    else
      ustd = sqrt(1. - rho**2)
      do i=1,nv
        uavg = rho*vval(i)  !  conditional mean of u, given v
        z0 = (u0 - uavg)/ustd  !  standardized value of u0 conditioned on v
        cpdfv(i) = gausspdf(vval(i))*gaussp(z0)*dv
        sum = sum + cpdfv(i)
        ccdfv(i) = sum
      enddo
    endif
  endif
  sum = ccdfv(nv)
  do i=1,nv
    ccdfv(i) = ccdfv(i)/sum
    cpdfv(i) = cpdfv(i)/sum
  enddo
  return
end

!!********************************!!******************************************************* 
!35   compute conditional distribution of bivariate std normal
!     variable, v, given the other std normal value is equal to u.
!!********************************!!******************************************************* 

subroutine cgauss_old(cgauss1)
use PrecStru
implicit none 
type(Scgauss) ::  cgauss1
!     u    =given value of correlated standard normal deviate
!     u0   =value of u corresponding to g(u)=0, u<=u0
!     rho  =coefficient of correlation between u and v
!     nv   =number of values of vval
!     vmin =minimum value of vval range
!     vmax =maximum value of vval range
!     vval =values of standard normal deviate v
!     cpdfv=probability density function of v
!     ccdfv=cumulative distribution function of v
real(r8),pointer :: vval(:),cpdfv(:),ccdfv(:)
integer nv,i
real(r8) u,u0,rho,vmin,vmax
real(r8) dv,std,sum,vavg,avgzv,avgval
real(r8) z0,vv,uavg,fpz,fpop,zv
real(r8) gausspdf,gaussp
 
  u=cgauss1.u
  u0=cgauss1.u0
  nv=cgauss1.nv
  rho=cgauss1.rho
  vmin=cgauss1.vmin
  vmax=cgauss1.vmax

  vval=>cgauss1.vval
  cpdfv=>cgauss1.cpdfv
  ccdfv=>cgauss1.ccdfv

  dv=(vmax - vmin)/float(nv)
  std=sqrt(abs(1. - rho**2))
 
  sum=0.
  if (u.gt.u0) then  !  v is normal conditional dististribution, given u
    vavg=rho*u
    avgzv=0.
    avgval=0.
    do i=1,nv
      zv=vmin + (i-.5)*dv
      vval(i)=vavg + zv*std
      cpdfv(i)=gausspdf(zv)*dv
      sum=sum + cpdfv(i)
      ccdfv(i)=sum
      avgzv=avgzv + zv*cpdfv(i)
      avgval=avgval + vval(i)*cpdfv(i)
    enddo
    avgzv=avgzv
    avgval=avgval
    !write (*,"('vavg=',f10.4/'avgzv=',f10.4/'avgval=',f10.4/'ccdfv(nv)=',f10.4)") vavg,avgzv,avgval,ccdfv(nv)
               
  else  !  conditional distribution of v is integrated over u<=u0
    fpz=gaussp(u0)
    fpop=1. - fpz
    avgval=0.
    do i=1,nv
          vval(i)=vmin + (i-.5)*dv
          vv=vval(i)
          avgval=avgval + vv
    enddo
    avgval=avgval/nv
    !write (*,"(/'subroutine cgauss'//'std=',f10.4/'avgval=',f10.4/'nv=',i5/1x)") std,avgval,nv
  
    do i=1,nv
      uavg=rho*vval(i)  !  conditional mean of u, given v
      z0=(u0 - uavg)/std  !  standardized value of u0 conditioned on v
      cpdfv(i)=gausspdf(vval(i))*gaussp(z0)*dv/fpz
      sum=sum + cpdfv(i)
      ccdfv(i)=sum
    enddo
    !write (*,"('avgval=',f10.4/'ccdfv(nv)=',f10.4)") avgval,ccdfv(nv)
 
  endif
  do i=1,nv
    ccdfv(i)=ccdfv(i)/ccdfv(nv)
    cpdfv(i)=cpdfv(i)/ccdfv(nv)
  enddo
  return
end

!!********************************!!******************************************************* 
!34 compute standardized normal variate corresponding to cumulative probability=f.
!     
!!********************************!!******************************************************* 
function gausspi(f)
use PrecStru
implicit none 
real(r8)  gausspi,f,f1,f2,df,avgf,x1,x2,dfdx,tol
real(r8) gaussp
integer it,maxit

  gausspi=0.
  it=0
  maxit=100
  tol=.00001
  x1=-1.
  x2=1.
  f1=gaussp(x1)
  do while (it.le.maxit)
    f2=gaussp(x2)
    df=f2-f1
    avgf=(f1 + f2)/2.
    if (abs(avgf).gt.0.) then
      if (abs(df/avgf).lt.tol) then
        gausspi=x2
        return
      endif
    else
      return
    endif
    !if (it.gt.maxit) return
    it=it+1
    dfdx=df/(x2-x1)
    x1=x2
    f1=f2
    x2=x2 - (f2-f)/dfdx
  end do
end
!!********************************!!******************************************************* 
!33  compute cumulative probability that a standard normal deviate is less than x.
!     
!!********************************!!******************************************************* 

function gaussp(x)
use PrecStru
implicit none 
real(r8)  x,y,erfcc,gaussp
  y=x/sqrt(2.)
  gaussp=(2.-erfcc(y))/2.
  return
end



!!********************************!!******************************************************* 
!32
!!********************************!!******************************************************* 

function erfcc(x)
use PrecStru
implicit none 
real(r8)   z,x,erfcc,t
      z=abs(x)      
      t=1./(1.+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
         t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
         t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.-erfcc
      return
end

!!********************************!!******************************************************* 
!31 compute gauss probability density of standard normal deviate z
!!********************************!!******************************************************* 
function gausspdf(z)
use PrecStru
implicit none 
real(r8) gausspdf,z
  gausspdf=0.3989423*exp(-(z**2)/2.)
  return
end
!!********************************!!******************************************************* 
!30 Compute "zero precipitation" threshold
! Original Author:  John Schaake 20080815
! Revise : Aizhong Ye 2012-5-7
!!********************************!!******************************************************* 
subroutine threshold (threshold1)
use PrecStru
implicit none 
type(Sthreshold) ::  threshold1
type(Ssort) ::  sort1
real(r8),pointer :: obs(:)  !observations array      
real(r8) pthresh,sum,amt,fraction
integer n,i    
  allocate(obs(threshold1.nobs))  
  obs=threshold1.obs
  pthresh=0.
  n=threshold1.nobs
  fraction=0.d0
  !if (n.gt.maxobs) n=maxobs
  !do i=1,n
  !  sortobs(i)=obs(i)
  !enddo
  sort1.a=> obs
  sort1.n=n
  call sort (sort1)
  sum=0.
  do i=1,threshold1.nobs
    sum=sum + obs(i)
  enddo
  amt=0.
  i=threshold1.nobs
  do while( i>=1.and.fraction.lt.threshold1.pop_fraction)
    amt=amt + obs(i)
    fraction=amt/sum    
    i=i-1     
  enddo
  if (i.eq.0) then 
    pthresh=obs(i+1)-1
  else
    pthresh=obs(i+1)
  end if
  threshold1.pthresh=pthresh
  deallocate (obs)
end
      
!!********************************!!******************************************************* 
!29 compute cumulative probability that a weib deviate is less than x.
!!********************************!!******************************************************* 

function weibp(a,b,x)
use PrecStru
implicit none 
real(r8) weibp,a,b,x
  if (x.lt.0..or.a.lt.0..or.b.lt.0.) then
    weibp=0.
  else
    weibp=1. - exp(-(x/a)**b)
  endif
  return
end

!!********************************!!******************************************************* 
!28 compute b-parameter of weibull distribution for given cv
!!********************************!!******************************************************* 

function weibb(cv)
use PrecStru
implicit none 
real(r8) weibb,cv,b1,b2,cv1,cv2
real(r8) weibcv,dcv,dcvdb
integer i,maxit
  maxit=100
  if (cv.le.0.) then
    weibb=0.
  else
    b2=cv**(-0.875)
    b1=0.95*b2
    cv1=weibcv(b1)
    cv2=weibcv(b2) 
	i=1
	dcv=cv2 - cv1
	do while (i.ne.maxit.and.abs(dcv).ge.0.0001 )
      i=i+1  	  
	  dcvdb=dcv/(b2-b1)     
      b1=b2
      b2=b2 - (cv2-cv)/dcvdb 
      cv1=cv2
      cv2=weibcv(b2)   
	  dcv=cv2 - cv1       
    end do 
    weibb=b2
  endif
  return
end


!!********************************!!******************************************************* 
!27  compute cv for weibull distribution with parameter=b
!!********************************!!******************************************************* 

function weibcv(b)
use PrecStru
implicit none 
real(r8) weibcv,b,g1,g2,gammln
  if (b.le.0.) then
        weibcv=0.01
  else
        g1=dexp(gammln(1.+2./b))
        g2=dexp(gammln(1.+1./b))
        weibcv=sqrt(abs((g1/(g2**2)) - 1.))
  endif
  return
end


!!********************************!!******************************************************* 
!*****************************26 Gamma function*****************************gammp,a,x,gln
!!********************************!!******************************************************* 

function gammp(a,x)
use PrecStru
implicit none 
type(Sgammcf) ::  gammcf1
real(r8) a,x,gammp
  gammcf1.a=a
  gammcf1.x=x
  if(x.lt.0..or.a.le.0.) stop 'error in gammp'
  if(x.lt.a+1.)then
    call gser(gammcf1)
    gammp=gammcf1.gammcf
  else
    call gcf(gammcf1)
    gammp=1.-gammcf1.gammcf
  endif 
end
!!********************************!!******************************************************* 
!*****************************25 Gamma function*****************************
!!********************************!!******************************************************* 

subroutine gser(gammcf1) !gamser,a,x,gln
use PrecStru
implicit none 
type(Sgammcf) ::  gammcf1
integer n,itmax
real(r8) eps
parameter (itmax=100,eps=3.e-7)
real(r8)   a,x,gln,gammln
real(r8)  ap,sum,del
  a=gammcf1.a
  x=gammcf1.x   
 
  gln=gammln(a)
  gammcf1.gln=gln
  if(x.le.0.)then
    !if(x.lt.0.)pause
    gammcf1.gammcf=0.
    return
  endif
  ap=a
  sum=1./a
  del=sum
  n=1 
  do while (n.ne.itmax .and. abs(del).ge.abs(sum)*eps)
    ap=ap+1.
    del=del*x/ap
    sum=sum+del
	n=n+1
        !if(abs(del).lt.abs(sum)*eps)go to 1
  end do 
  !if (x.le.0) then
  !  print *,'gser<0 x',x
  !  pause 4
  !end if
  gammcf1.gammcf=sum*dexp(-x+a*dlog(x)-gln)
end

!!********************************!!******************************************************* 
!*****************************24 Gamma function*****************************
!!********************************!!******************************************************* 
subroutine gcf(gammcf1)  !
use PrecStru
implicit none 
type(Sgammcf) ::  gammcf1
integer n,itmax
real(r8) eps
parameter (itmax=100)
real(r8)   a,x,gln  ,GAMMLN
real(r8) a0,a1,b0,b1,fac,anf,an,ana,g,gold
  eps=3.e-7
  a=gammcf1.a
  x=gammcf1.x     
  gln=gammln(a)
  gammcf1.gln=gln
  gold=0.
  a0=1.
  a1=x
  b0=0.
  b1=1.
  fac=1.
  n=1
  g=1.
  do while (n.ne.itmax .and. abs((g-gold)/g).ge.eps)
    n=n+1
    an=float(n)
    ana=an-a
    a0=(a1+a0*ana)*fac
    b0=(b1+b0*ana)*fac
    anf=an*fac
    a1=x*a0+anf*a1
    b1=x*b0+anf*b1
    if(a1.ne.0.)then
      fac=1./a1
      gold=g
      g=b1*fac
      !if(abs((g-gold)/g).lt.eps)go to 1
    endif
  end do 
  !if (x.le.0) then
  !  print *,'gcf<0 x',x
  !  pause 3
  !end if
  gammcf1.gammcf=dexp(-x+a*dlog(x)-gln)*g
  gammcf1.gln=gln 
end
!!********************************!!******************************************************* 
!*****************************23 Gamma function*****************************
!!********************************!!******************************************************* 
function gammln(xx)
use PrecStru
implicit none 
real(r8) cof(6),stp,half,one,fpf,x,tmp,ser
real(r8) xx,gammln
integer j
  data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0, &
         -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
  data half,one,fpf/0.5d0,1.0d0,5.5d0/
    x=xx-one
    tmp=x+fpf
	if (tmp.le.0) then
	  print *,'gammaln<0', tmp
	  pause 1
	endif
    tmp=(x+half)*dlog(tmp)-tmp
    ser=one
    do  j=1,6
      x=x+one
      ser=ser+cof(j)/x
    end do
    !if (stp*ser.le.0) then
    !  print *,'gammaln<0 stp*ser', stp*ser
	!  pause 2
    !end if
    gammln=tmp+dlog(stp*ser)
    return
end function


!!********************************!!******************************************************* 
!*****************************22 sort a *****************************
!!********************************!!******************************************************* 
subroutine sort(sort1)
use PrecStru
integer i,j,k,n
real(r8) temp
real(r8),pointer :: a(:)
type(Ssort) ::  sort1
  n=sort1.n
  a=>sort1.a
  do i=1,n-1  
    k=i
    do j=i+1,n    
      if (a(k)>a(j))  k=j
    end do
    if (k.ne.i) then     
      temp=a(k)
      a(k)=a(i)
      a(i)=temp
    end if
  end do
end subroutine 
!!********************************!!******************************************************* 
!***************************21 generate Canonical Events !!**********************************
!!********************************!!******************************************************* 
subroutine generateCaE(gCanE)
use PrecStru
implicit none 
type(SgenerateCaE) ::  gCanE 
integer  i,j,k,n  

  do i=1,gCanE.ndays
    do j=1,gCanE.Events
      gCanE.Eob(i,j)=0
      do k=gCanE.Estart(j),gCanE.Estop(j)
        if (k.lt.1) then
          n=i+k
          if (n.lt.1) n=1
          if (k.ne.0) gCanE.Eob(i,j)=gCanE.Eob(i,j)+gCanE.hob(n,1)  
        else
          gCanE.Eob(i,j)=gCanE.Eob(i,j)+gCanE.ob(i,k)
        end if
      end do
      if (gCanE.Estart(j).gt.0) then
          gCanE.Eob(i,j)=gCanE.Eob(i,j)/(gCanE.Estop(j)-gCanE.Estart(j)+1)
      else
          gCanE.Eob(i,j)=gCanE.Eob(i,j)/(gCanE.Estop(j)-gCanE.Estart(j))
      end if
    end do
  end do !anonical Events data
 
end subroutine
!!********************************!!******************************************************* 
! 20 Write  statistical information between observed and forecasting events 
!!********************************!!******************************************************* 
!Nash-Sutcliffe efficiency  and correlation coefficient 
subroutine WritestatOFE(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   outputfile
character(c1000) ::   tempc
integer  i,j,k
  k=365 
   
  write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(CanE.ReFile),CanE.x,'\',CanE.y,'E.txt' 
  if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
    open (10,file=outputfile) 
	  write (10,'(2a,i3,a,i3)') trim('correlation-coefficient'),char(9),365,char(9),CanE.Events   
      do i=1,k  	   
        write (tempc,'(365(a,f5.2))') (char(9),CanE.r(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
    !close(10)

    !write(outputfile,'(a,a3,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),'EE\',CanE.x,'.',CanE.y,'.txt'  
    !if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
    !open (10,file=outputfile)     
	  write (10,'(2a,i3,a,i3)') trim('Nash-Sutcliffe-efficiency'),char(9),365,char(9),CanE.Events   
      do i=1,k   
        write (tempc,'(365(a,f5.2))') (char(9),CanE.DC(i,j),j=1,CanE.Events)            
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 
    !close(10)

    !write(outputfile,'(a,a3,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),'BE\',CanE.x,'.',CanE.y,'.txt'  
    !if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
    !open (10,file=outputfile)     
	  write (10,'(a,a,i3,a,i3)') trim('Bias'),char(9),365,char(9),CanE.Events   
      do i=1,k   
        write (tempc,'(365(a,f5.2))') (char(9),CanE.Bias(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

	  write (10,'(a,a,i3,a,i3)') trim('RMSE'),char(9),365,char(9),CanE.Events   
      do i=1,k   
        write (tempc,'(365(a,f7.2))') (char(9),CanE.RMSE(i,j),j=1,CanE.Events)   
		call delblank(tempc)
        write (10,'(a)')trim(tempc)
      end do 

    close(10)

end subroutine 

!!********************************!!******************************************************* 
! 19 Write  statistical information between observed and forecasting 
!!********************************!!******************************************************* 
subroutine WritestatOF(CanE)
use PrecStru
implicit none 
type(CanonicalE) ::  CanE
character(c100) ::   outputfile
character(c1000) ::   tempc
integer  i,j  
   
  write(outputfile,'(a,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),CanE.x,'\',CanE.y,'.txt' 
  if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile)    
	write (10,'(a,a,i3,a,i3)') trim('correlation-coefficient'),char(9),365,char(9),CanE.leadT   
    do i=1,365   
      write (tempc,'(365(a,f5.2))') (char(9),CanE.r(i,j),j=1,CanE.leadT)   
      call delblank(tempc)
      write (10,'(a)')trim(tempc)
    end do 
  !lose(10)

  !write(outputfile,'(a,a2,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),'E\',CanE.x,'.',CanE.y,'.txt'  
  !if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
  !open (10,file=outputfile)     
    write (10,'(a,a,i3,a,i3)') trim('Nash-Sutcliffe-efficiency'),char(9),365,char(9),CanE.leadT   
    do i=1,365   
      write (tempc,'(365(a,f5.2))') (char(9),CanE.DC(i,j),j=1,CanE.leadT)   
      call delblank(tempc)
      write (10,'(a)')trim(tempc)
    end do 
  !lose(10)

  !write(outputfile,'(a,a2,i3.3,a1,i3.3,a4)')trim(CanE.ReFile),'B\',CanE.x,'.',CanE.y,'.txt'  
  !if (CanE.flagOs .eq. 2)   call slashchange(outputfile)
  !open (10,file=outputfile)     
    write (10,'(a,a,i3,a,i3)') trim('Bias'),char(9),365,char(9),CanE.leadT   
    do i=1,365   
    write (tempc,'(365(a,f5.2))') (char(9),CanE.Bias(i,j),j=1,CanE.leadT)   
      call delblank(tempc)
      write (10,'(a)')trim(tempc)
    end do 

	write (10,'(a,a,i3,a,i3)') trim('RMSE'),char(9),365,char(9),CanE.leadT   
    do i=1,365   
      write (tempc,'(365(a,f7.2))') (char(9),CanE.RMSE(i,j),j=1,CanE.leadT)   
      call delblank(tempc)
      write (10,'(a)')trim(tempc)
    end do 

  close(10)

end subroutine 


!!********************************!!******************************************************* 
!18 read observed  data  for shaake shuffle
!!********************************!!******************************************************* 
subroutine readDataEp(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays 
integer julday,jday1,jday2
character*4 units_bmap06
real xobs(4) 

  if (rData1.flagOs .eq. 2)   then
    write(inputfile,'(a,i3.3,a1,i3.3)')trim(rData1.EpFile),rData1.x,'/',rData1.y 
  else
    write(inputfile,'(a,i3.3,a1,i3.3)')trim(rData1.EpFile),rData1.x,'\',rData1.y 
  end if
  !write(inputfile,'(a,i3.3,a1,i3.3)')trim(rData1.EpFile),rData1.x,'\',rData1.y 
  open (10,file=inputfile,form='unformatted')
    read (10) year1,year2,ndays, units_bmap06
    jday1=julday(1,1,year1)
    jday2=julday(1,1,rData1.BeginYearEns)
    if (rData1.BeginYearEns.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYearEns,'<',year1
	  stop
    end if
    if (rData1.EndYearEns.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYearEns,'>',year2
	  stop
    end if
    j=jday2-jday1
    do i=1,j
      read (10)        
    end do 
    jday1=julday(12,31,rData1.EndYearEns)
    ndays=jday1-jday2+1+rData1.leadT
  
    do i=1,ndays
      read (10) (xobs(j),j=1,4)
      rData1.ep(i,1)=max(0.0,xobs(1)+xobs(2)+xobs(3)+xobs(4))
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ep(i,j)=rData1.ep(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ep(i,j)=rData1.ep(i,j-1)
      end do
    end do
  close(10)

end subroutine

!!********************************!!******************************************************* 
!18 read observed  data  
!!********************************!!******************************************************* 
subroutine readDataob(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
integer i,j,year1,year2,ndays 
integer julday,jday1,jday2
character*4 units_bmap06
real xobs(4) 
  if (rData1.flagOs .eq. 2)   then
    write(inputfile,'(a,i3.3,a1,i3.3)')trim(rData1.ObFile),rData1.x,'/',rData1.y 
  else
    write(inputfile,'(a,i3.3,a1,i3.3)')trim(rData1.ObFile),rData1.x,'\',rData1.y 
  end if
  open (10,file=inputfile,form='unformatted')
    read (10) year1,year2,ndays, units_bmap06
    jday1=julday(1,1,year1)
    jday2=julday(1,1,rData1.BeginYear)
 
    if (rData1.BeginYear.lt.year1) then
      write(*,*)"the start year is too early, stop",rData1.BeginYear,'<',year1
	  stop
    end if
    if (rData1.EndYear.ge.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYear,'>',year2
	  stop
    end if
    j=jday2-jday1
    do i=1,j
      read (10)        
    end do 
    !if (rData1.nparint.eq.1)  then
    !  ndays=rData1.ndays      
    !else
    ndays=rData1.ndays+rData1.leadT
    !end if
    do i=1,ndays
      read (10) (xobs(j),j=1,4)
      rData1.ob(i,1)=max(0.0,xobs(1)+xobs(2)+xobs(3)+xobs(4))
    end do
    do j=2,rData1.leadT
      do i=1,ndays-j+1
        rData1.ob(i,j)=rData1.ob(i+j-1,1)
      end do
    end do
    do j=2,rData1.leadT
      do i=ndays-j+2,ndays
        rData1.ob(i,j)=rData1.ob(i,j-1)
      end do
    end do
  close(10)

end subroutine
!!********************************!!******************************************************* 
!17 read forecast  data  
!!********************************!!******************************************************* 

subroutine readDatafcst(rData1)
use PrecStru
implicit none
type(SreadData) :: rData1 
character(c100) ::   inputfile
character(len=50) ::   tempc
integer i,j,k,year1,imonth1,iday1 ,year2,imonth2,iday2 
integer julday,jday1,jday2,cday
 
  if (rData1.flagOs .eq. 2)   then
    write(inputfile,'(2(a,i3.3),a)')trim(rData1.SiFile),rData1.x,'/',rData1.y,'.txt'          
  else
    write(inputfile,'(2(a,i3.3),a)')trim(rData1.SiFile),rData1.x,'\',rData1.y,'.txt'          
  end if

  !write(inputfile,'(a,2(i3.3,a1),i3.3,a4)')trim(rData1.SiFile),rData1.x,'\',rData1.x,'.',rData1.y,'.txt'  
  open (10,file=inputfile )
    read (10,*) tempc,year1,imonth1,iday1 
    read (10,*) tempc,year2,imonth2,iday2  
    read (10,*)
    jday1=julday(imonth1,iday1,year1)
    jday2=julday(1,1,rData1.BeginYear)
 
    if (jday2.lt.jday1) then
      write(*,*)"the start date is too early, stop",rData1.BeginYear,'<',year1
	  stop
    end if 
    if (rData1.EndYear.Gt.year2) then
      write(*,*)"the end year is too late, stop",rData1.EndYear,'>',year2
	  stop
    end if 

    !rData1.si=0.0
    if (rData1.nparint.eq.0)  then  !GFS data
      j=jday2-jday1
      do i=1,j
        read (10,*)        
      end do 
      do i=1,rData1.ndays 
        read (10,*)(rData1.si(i,j),j=1,rData1.leadT)  
        do j=1,rData1.leadT
          rData1.si(i,j)=rData1.si(i,j)/100.0 !(:,:)
		  if (rData1.si(i,j)>1000) rData1.si(i,j)=0
        end do
      end do  
    else   !CFS data    
      if (rData1.nparint.eq.1) then
        j=jday2-jday1
      else
        j=(rData1.BeginYear-year1)*(365/rData1.nparint)
      end if
      do i=1,j
        read (10,*)        
      end do 

      if (rData1.nparint.eq.1) then
        do i=1,rData1.ndays                   
          read (10,*)(rData1.si(i,k),k=1,rData1.leadT)  
          do k=1,rData1.leadT
            if (rData1.flagT.eq.0) then
              rData1.si(i,k)=rData1.si(i,k)/100.0  
            else
              rData1.si(i,k)=rData1.si(i,k)/100.0+273.16 !-200 
            end if
          end do     
        end do          
      else
        do i=rData1.BeginYear,rData1.EndYear
          do j=1,365,rData1.nparint
            cday=JULDAY(1,1,i)+j-jday2  
		    if(((mod(i,4)==0.and.mod(i,100)/=0).or. mod(i,400)==0).and.j.gt.59) &  !leap year
                 cday=cday+1       
            read (10,*)(rData1.si(cday,k),k=1,rData1.leadT)  
            do k=1,rData1.leadT
              if (rData1.flagT.eq.0) then
                rData1.si(cday,k)=rData1.si(cday,k)/100.0  
              else
                rData1.si(cday,k)=rData1.si(cday,k)/100.0+273.16 
              end if
            end do     
          end do
        end do
      end if
    
    end if
  close(10) 
end subroutine


!!********************************!!******************************************************* 
!***************************16 change '\' to '/' in a character string!*********************************
!!********************************!!******************************************************* 
subroutine slashchange(anistring) 
use PrecStru
implicit none
character(c100) ::   anistring
integer i

      do i=1,c100   !for linux
        if (anistring(i:i)=='\') anistring(i:i)='/'
      enddo  

end subroutine 
!!********************************!!******************************************************* 
!***************************15 Statistical information of daily data !!**********************************
!!********************************!!******************************************************* 
subroutine mapstats24 (stat24a) 
use PrecStru
implicit none
integer  ndays
!type(Obsdata) ::  Obsdata1
type(stat24) ::  stat24a

real(r8)  avg(12),pop(12),cavg(12),ccv(12)
integer*4 npos(12),nobs(12)
real(r8)  pthresh,val,tempf
integer j
integer   months(12),iyear,imonth,iday,ipd  
  pthresh=stat24a.pthresh
  do j=1,12
    avg(j)=0.
        npos(j)=0
        nobs(j)=0
        cavg(j)=0.
        ccv(j)=0.
        pop(j)=0.
  enddo
  do iyear=stat24a.BeginYear,stat24a.EndYear   
    j= iyear-stat24a.BeginYear+1
    if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
      months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
      months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    endif 
    do imonth=1,12 
      do iday=1,months(imonth) 
            val=stat24a.Obdata(j,imonth,iday)
            if (val.ge.0.) then
              nobs(imonth)=nobs(imonth) + 1  ! number of observed data
              avg(imonth)=avg(imonth) + val  !average all data
              if (val.gt.pthresh) then
                npos(imonth)=npos(imonth) + 1     ! Number more than pthresh
                cavg(imonth)=cavg(imonth) + val !average more than pthresh
				ccv(imonth)=ccv(imonth) + val**2 
              endif
            endif
        enddo
    enddo
  enddo  
  do imonth=1,12
        if (nobs(imonth).gt.0) then
          avg(imonth)=avg(imonth)/nobs(imonth)
          if (npos(imonth).gt.0) then
            pop(imonth)=float(npos(imonth))/float(nobs(imonth))
            cavg(imonth)=cavg(imonth)/npos(imonth)
			tempf=max(ccv(imonth)/npos(imonth)-cavg(imonth)*cavg(imonth),0.0)
			ccv(imonth)=sqrt(tempf)
			ccv(imonth) =ccv(imonth) /cavg(imonth) 
          endif
        endif
  enddo
    !V=σ/|μ|，where σ=√∑(xi-u)^2/(n-1)，u=(∑xi)/n
  !do iyear=Obsdata1.BeginYear,Obsdata1.EndYear   
  !  j= iyear-Obsdata1.BeginYear+1
  !  if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
  !    months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
  !  else
  !    months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  !  endif 
  !  do imonth=1,12 
  !    do iday=1,months(imonth) 
  !          val=Obsdata1.Obdata(Obsdata1.x,Obsdata1.y,j,imonth,iday)
  !          if (val.ge.pthresh) then                
  !               ccv(imonth)=ccv(imonth) + (val-cavg(imonth))**2  !averrage square data value               
  !          endif
  !     enddo
  !  enddo
  !enddo  
  !    do imonth=1,12
  !                 
  !        if (npos(imonth).gt.1) then 
  !          ccv(imonth)=sqrt(ccv(imonth)/(npos(imonth) -1))
  !          ccv(imonth)=ccv(imonth)/cavg(imonth)
  !        endif
         
  !    enddo
  stat24a.avg=avg
  stat24a.pop=pop
  stat24a.cavg=cavg
  stat24a.ccv=ccv
  stat24a.npos=npos
  stat24a.nobs=nobs 

 end subroutine
!!********************************!!******************************************************* 
!***************************14 Statistical information of 6 hours data !!**********************************
!!********************************!!******************************************************* 

subroutine mapstats06 (stat06a) 
use PrecStru
implicit none
type(stat06) ::  stat06a

integer  ndays
type(Obsdata) ::  Obsdata1
real(r8)  avg(4,12),pop(4,12),cavg(4,12),ccv(4,12)
integer*4 npos(4,12),nobs(4,12)
real(r8)  pthresh,val,tempf
integer j
integer   months(12),iyear,imonth,iday,ipd  
  pthresh=stat06a.pthresh
 
  avg =0.
  npos =0
  nobs=0
  cavg =0.
  ccv =0.
  pop =0.
  do ipd=1,4 
    do iyear=stat06a.BeginYear,stat06a.EndYear   
    j= iyear-stat06a.BeginYear+1
    if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
      months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
      months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    endif 
        do imonth=1,12 
      do iday=1,months(imonth) 
            val=stat06a.Obdata06(j,imonth,iday,ipd)     
            if (val.ge.0.) then
              nobs(ipd,imonth)=nobs(ipd,imonth) + 1  ! number of observed data
              avg(ipd,imonth)=avg(ipd,imonth) + val  !average all data
              if (val.gt.pthresh) then
                npos(ipd,imonth)=npos(ipd,imonth) + 1     ! Number more than pthresh
                cavg(ipd,imonth)=cavg(ipd,imonth) + val !average more than pthresh
                ccv(ipd,imonth)=ccv(ipd,imonth) + val**2  !averrage square data value               

              endif
            endif
          enddo
        enddo
      enddo  
  
      do imonth=1,12
        if (nobs(ipd,imonth).gt.0) then
          avg(ipd,imonth)=avg(ipd,imonth)/nobs(ipd,imonth)
          if (npos(ipd,imonth).gt.0) then
            pop(ipd,imonth)=float(npos(ipd,imonth))/float(nobs(ipd,imonth))
            cavg(ipd,imonth)=cavg(ipd,imonth)/npos(ipd,imonth)
            tempf=max(ccv(ipd,imonth)/npos(ipd,imonth)-cavg(ipd,imonth)*cavg(ipd,imonth),0.0)
            !(ccv(ipd,imonth)-npos(ipd,imonth)*cavg(ipd,imonth)*cavg(ipd,imonth) )
            ccv(ipd,imonth)=sqrt(tempf)
            ccv(ipd,imonth)=ccv(ipd,imonth)/cavg(ipd,imonth)
          endif
        endif
      enddo
    !V=σ/|μ|，where σ=√∑(xi-u)^2/(n-1)，u=(∑xi)/n
  !do iyear=Obsdata1.BeginYear,Obsdata1.EndYear   
  !  j= iyear-Obsdata1.BeginYear+1
  !  if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
  !    months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
  !  else
  !    months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  !  endif 
  !      do imonth=1,12 
  !    do iday=1,months(imonth) 
  !          val=Obsdata1.Obdata06(j,imonth,iday,ipd)  
  !          if (val.ge.pthresh) then                
  !               ccv(ipd,imonth)=ccv(ipd,imonth) + (val-cavg(ipd,imonth))**2  !averrage square data value               
  !          endif
  !        enddo
  !      enddo
  !    enddo  
  !    do imonth=1,12
  !                
  !        if (npos(ipd,imonth).gt.1) then 
  !          ccv(ipd,imonth)=sqrt(ccv(ipd,imonth)/(npos(ipd,imonth) -1))
  !          ccv(ipd,imonth)=ccv(ipd,imonth)/cavg(ipd,imonth)
  !        endif
  !       
  !    enddo
    stat06a.avg(ipd,:)=avg(ipd,:)
    stat06a.pop(ipd,:)=pop(ipd,:)
    stat06a.cavg(ipd,:)=cavg(ipd,:)
    stat06a.ccv(ipd,:)=ccv(ipd,:)
    stat06a.npos(ipd,:)=npos(ipd,:)
    stat06a.nobs(ipd,:)=nobs(ipd,:)
  enddo

 end subroutine

!!********************************!!******************************************************* 
!***************************13 Write map statistical File!!**********************************
!!********************************!!******************************************************* 
subroutine Writemapstats(Obsdata1)
use PrecStru
implicit none 
type(Obsdata) ::  Obsdata1
character(c100) ::   outputfile
character(c1000) ::  tempc
integer  jmo,ipd  

  write(outputfile,'(a,i3.3,a1,i3.3,a)')trim(Obsdata1.BName06),Obsdata1.x,'\',Obsdata1.y,'s.txt'  
  if (Obsdata1.flagOs .eq. 2)   call slashchange(outputfile)
  open (10,file=outputfile) !,access='append'  
    !write (10,'(2i5)')  Obsdata1.x,Obsdata1.y
    write (10,'(11a)')"Month",char(9),"Nobs",char(9),"avg",char(9),"pop",char(9),"cavg",char(9),"ccv"
    do jmo=1,12
      write (tempc,'(2(i5,a),4(f10.4,a))') jmo,char(9),Obsdata1.nobs(5,jmo),char(9),&
                            Obsdata1.avg(5,jmo),char(9),Obsdata1.pop(5,jmo),char(9), &
                           Obsdata1.cavg(5,jmo),char(9),Obsdata1.ccv(5,jmo) 
      call delblank(tempc)
      write (10,'(a)')trim(tempc)
    enddo
    write (10,'(13a)')"Month",char(9),"Hour",char(9),"Nobs",char(9),"avg",char(9),"pop",char(9),"cavg",char(9),"ccv"
    do jmo=1,12
      do ipd=1,4
        write (tempc,'(3(i5,a),4(f10.4,a))') jmo,char(9),ipd,char(9),Obsdata1.nobs(ipd,jmo),char(9), &
                              Obsdata1.avg(ipd,jmo),char(9),Obsdata1.pop(ipd,jmo),char(9),&                             
                              Obsdata1.cavg(ipd,jmo),char(9), Obsdata1.ccv(ipd,jmo)
                                 
        call delblank(tempc)
        write (10,'(a)')trim(tempc)
      enddo
    enddo

  close(10)    

end subroutine 

!!********************************!!******************************************************* 
!***************************12 binary File to AscII File !!**********************************
!!********************************!!******************************************************* 

subroutine BinToAdata06(outputfile)
use PrecStru
implicit none
type(Scaldat) ::  caldat1
character(c100) ::   outputfile
integer i,j,year1,year2,ndays,iyear,imonth,iday
integer julday,jday1
character*4 units_bmap06
real  xobs(4)
  open (10,file=outputfile,form='unformatted')
  write(outputfile,'(a,a4)')trim(outputfile),'.txt'  
  open (20,file=outputfile)
    
    read (10) year1,year2,ndays, units_bmap06
    write(20,'(a4,TR1,a2,TR1,a2,TR1,a3,TR1,a3,a2)') &
        'date','h0','h6','h12','h18',trim(units_bmap06)
    jday1=julday(1,1,year1)
    do i=1,ndays  
      j=jday1 + i - 1
	  caldat1.julian=j
      call caldat(caldat1)   
	  iyear=caldat1.iyyy
	  imonth=caldat1.mm
	  iday=caldat1.id 
      read (10) (xobs(j),j=1,4)
      write (20,'(i4,a1,i2.2,a1,i2.2,3(f6.2,tr1),f6.2)') &
        iyear,'-',imonth,'-',iday,(xobs(j),j=1,4)
    end do

  close(10)
  close(20)

end subroutine


!!********************************!!******************************************************* 
!***************************11 Write Observed 6 hours data binary File!!**********************************
!!********************************!!******************************************************* 
subroutine WriteObdata06(Obsdata1)
use PrecStru
implicit none
integer JULDAY,ndays
type(Obsdata) ::  Obsdata1
character(c100) ::   outputfile
integer i,j
integer   months(12),iyear,imonth,iday,ipd     
character*4 units
real xobs(4)
  units='mm  '
  ndays=JULDAY(12,31,Obsdata1.EndYear)-JULDAY(1,1,Obsdata1.BeginYear)+1
  write(outputfile,'(a,i3.3,a1,i3.3)')trim(Obsdata1.BName06),Obsdata1.x,'\',Obsdata1.y  
  if (Obsdata1.flagOs .eq. 2)   call slashchange(outputfile)

  open (10,file=outputfile,form='unformatted') !,'(i4,tr1,i4,tr1,i7,tr1,a2)'
  write(10) Obsdata1.BeginYear,Obsdata1.EndYear,ndays,units

  do iyear=Obsdata1.BeginYear,Obsdata1.EndYear   
    j= iyear-Obsdata1.BeginYear+1
    if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
      months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
      months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    endif 
    do imonth=1,12 
      do iday=1,months(imonth) 
        do ipd=1,3
          xobs(ipd)=Obsdata1.Obdata06(j,imonth,iday,ipd+1)
        end do         
        if (iday+1>months(imonth)) then        
          if (imonth+1>12) then
            if(iyear+1>Obsdata1.EndYear ) then
            xobs(4)=Obsdata1.Obdata06(j,12,31,4)
          else
            xobs(4)=Obsdata1.Obdata06(j+1,1,1,1)
          endif

          else
            xobs(4)=Obsdata1.Obdata06(j,imonth+1,1,1)
          endif
        else
          xobs(4)=Obsdata1.Obdata06(j,imonth,iday+1,1)
        endif
        write (10) (xobs(i),i=1,4) !,'(3(f6.2,tr1),f6.2)'
      end do
    end do
  enddo 
  close(10)   

end subroutine 

!!********************************!!******************************************************* 
!***************************10 Read Observed daily data binary File!!**********************************
!!********************************!!******************************************************* 

subroutine ReadObdata(Obsdata1)
use PrecStru
implicit none
type(Obsdata) ::  Obsdata1
       
real :: daily(464,224,31) !(,Obsdata1.nrows,Obsdata1.NVara,31)  
integer   months(12),iyear,imonth
integer i,j,ii,jj
character(C100) ::  inputfile
 
  do iyear=Obsdata1.BeginYear,Obsdata1.EndYear    
    if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
      months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
      months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    endif 
    do imonth=1,12 
      print *,'year=',iyear,'   month=',imonth
      write(inputfile,'(a,i4,i2.2,a2)')trim(Obsdata1.BName),iyear,imonth,"01"  
      j=Obsdata1.ncols*Obsdata1.nrows   !*Obsdata1.NVara*months(imonth)
      if (Obsdata1.flagOs .eq. 2)   call slashchange(inputfile)
      open(10,file=trim(inputfile),status='old',form='unformatted',access='direct',recl=j)      
        do  i=1,months(imonth)  
          j=(i-1)*4+Obsdata1.CuVa		  
          read(10,rec=j)daily(:,:,i)   !(daily(:,:,:,i),i=1,months(imonth))!read month data (464:224:4:nday) months(imonth)*    
        end do		     
      close(10)  
      ii=(Obsdata1.k-1)*Obsdata1.nrows/Obsdata1.Tim+1
      jj=Obsdata1.k*Obsdata1.nrows/Obsdata1.Tim
      Obsdata1.Obdata(:,:,iyear-Obsdata1.BeginYear+1,imonth,:)=daily(:,ii:jj,:)
    end do
  enddo
end subroutine ReadObdata


!!********************************!!******************************************************* 
!***************************9 Read GFS File!!**********************************
!!********************************!!******************************************************* 
subroutine ReadGFS(GFSdata1)
use PrecStru
implicit none
integer i,j,k,x,y
type(GFSdata) ::  GFSdata1
character(c100)  asFile
real(r8) longtitude,latitude,tempf(14)
integer BeginYear,EndYear,cyear
integer months(12),iyear,imonth,iday

  BeginYear=GFSdata1.BeginYear
  EndYear=GFSdata1.EndYear

    do iyear=BeginYear,EndYear   
    cyear=iyear-BeginYear+1     
    if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
      months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
      months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    endif
        do imonth=1,12
      do iday=1,months(imonth)
        !2.2. read GFS precipitation data          G:\USA\usa_apcp_gfs_ensmean
            write(asFile,'(a,i4,i2.2,i2.2,a4)')trim(GFSdata1.BName),&
        iyear,imonth,iday,".txt"  
            print*,iyear,imonth, iday        
      if (GFSdata1.flagOs .eq. 2)   call slashchange(asFile)
      open(10,File=asFile)
      j=GFSdata1.ncols*GFSdata1.nrows
      do i=1,j
        read(10,*) longtitude,longtitude,latitude,tempf(1:14)
        x=int((longtitude-GFSdata1.xllcorner)/GFSdata1.cellsize+1)
        y=int((latitude-GFSdata1.yllcorner)/GFSdata1.cellsize+1)
        do k=1,14
          GFSdata1.Gdatat(x,y,cyear,imonth,iday,k)=tempf(k) 
        end do
      end do
      close(10)
    end do
    end do
  end do
end subroutine ReadGFS

!!********************************!!******************************************************* 
!***************************8 Write ASCII File!!**********************************
!!********************************!!******************************************************* 
subroutine WriteASC(MyASC,asFile)
use PrecStru
implicit none
type(ASCFile) ::  MyASC
character(c100)  asFile
character(c100) asRecord
character(20) tempa 
integer i,j
  if (MyASC.flagAF/=1) then  
    print *,"No data for output! "
    return
  endif
  open(10,File=asFile)
  write(asRecord,"('Begin writing file of ',A30)")asFile
  call Wrecord(asRecord)
  write(10,"('ncols',TR1,I6)")MyASC.ncols
  write(10,"('nrows',TR1,I6)")MyASC.nrows
  write(10,"('xllcorner',TR1,F15.4)")MyASC.xllcorner
  write(10,"('yllcorner',TR1,F15.4)")MyASC.yllcorner
  write(10,"('cellsize',TR1,F8.3)")MyASC.cellsize
  write(10,"('NODATA_value',TR1,F10.3)")MyASC.NODATA_value
  do j=MyASC.nrows,1,-1
    do  i=1,MyASC.ncols 
      select case (int(MyASC.MValue(i,j)+0.5))
    case (-999:-1)
      write(10,"(TR1,F6.1,\)")MyASC.MValue(i,j) 
    case (-9998)
      write(10,"(TR1,I5,\)")int(MyASC.MValue(i,j)) 
    case (0:9)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
           write(10,"(TR1,I2,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F5.2,\)")MyASC.MValue(i,j)        
        endif
    case (10:99)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I2,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F5.2,\)")MyASC.MValue(i,j)        
        endif
    case (100:999)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I3,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F6.2,\)")MyASC.MValue(i,j)        
        endif
    case (1000:9999)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I4,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F7.2,\)")MyASC.MValue(i,j)        
        endif
    case (10000:99999)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I5,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F8.2,\)")MyASC.MValue(i,j)        
        endif
    case (100000:999999)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I6,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F9.2,\)")MyASC.MValue(i,j)        
        endif
    case (1000000:9999999)
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I7,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F10.2,\)")MyASC.MValue(i,j)        
        endif
      case default
        if (abs(int(MyASC.MValue(i,j))-MyASC.MValue(i,j))<0.0001) then
        write(10,"(TR1,I9,\)")int(MyASC.MValue(i,j)) 
        else
          write(10,"(TR1,F12.2,\)")MyASC.MValue(i,j)        
        endif
      end select     
    end do   
    write(10,*)
  end do
  close(10) 
  write(asRecord,"('Finish writing file of ',A20)")asFile
  call Wrecord(asRecord)
  MyASC.flagAF=0
end subroutine WriteASC
!!********************************!!******************************************************* 
!****************************7 Read ASCII File!!**********************************
!!********************************!!******************************************************* 
subroutine ReadASC(MyASC,asFile)
use PrecStru
implicit none
type(ASCFile) ::  MyASC
character(c100)  asFile
character(c100) asRecord
character(20) tempa 
integer i,j
logical alive
  inquire(file=asFile,exist=alive)
  if (MyASC.flagAF==1) then  
  deallocate(MyASC.MValue)  
  endif
  MyASC.flagAF=0
  MyASC.SumN=0
  if(alive) then
    open(10,File=asFile,status='old')
    read(10,*)tempa,MyASC.ncols  
    if (tempa/='ncols') then     
      Write(*,*)"Wrong! ",asFile
     return
    endif    
  write(asRecord,"('Begin reading file of ',A20)")asFile
    call Wrecord(asRecord)
    read(10,*)tempa,MyASC.nrows
  read(10,*)tempa,MyASC.xllcorner
  read(10,*)tempa,MyASC.yllcorner
  read(10,*)tempa,MyASC.cellsize
  read(10,*)tempa,MyASC.NODATA_value
    allocate(MyASC.MValue(MyASC.ncols,MyASC.nrows))
    do j=MyASC.nrows,1,-1           
      read(10,*)(MyASC.MValue(i,j),i=1,MyASC.ncols)
      do i=1,MyASC.ncols             
        if (abs(MyASC.MValue(i,j)-MyASC.NODATA_value)>0.01) MyASC.SumN=MyASC.SumN+1           
      end do
    end do
     close(10)   
     MyASC.flagAF=1 
  else
    Write(*,*) "No ",asFile
  return
  endif
  write(asRecord,"('Finish reading file of ',A20)")asFile
  call Wrecord(asRecord)
end subroutine ReadASC
!!********************************!!******************************************************* 
!*************************6 the log:record the work of system*****************************
!!********************************!!******************************************************* 
subroutine Wrecord(asRecord)
use PrecStru
implicit none
character(c100) asRecord
character(c100) asFile
character(10) Temp
INTEGER DATE_TIME (8)
logical alive
  asFile="Record.txt"
  inquire(file=asFile,exist=alive)
  CALL DATE_AND_TIME (Temp,Temp,Temp,DATE_TIME)
  if(alive) then 
    open(999,File=asFile,status='old',position='append')
  else  
    open(999,File=asFile)
    write(999,"('Date',TR7,'Time',TR5,'Program')")  
  endif     
    write(999,"(I4.4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2,1x,A60)")&
    DATE_TIME(1),DATE_TIME(2),DATE_TIME(3),&
    DATE_TIME(5),DATE_TIME(6),DATE_TIME(7),asRecord   
  close(999)
end subroutine Wrecord


!!********************************!!******************************************************* 
!*************************5 the log:CALDAT*****************************
!!********************************!!******************************************************* 
subroutine caldat(caldat1)
use PrecStru
implicit none
type(Scaldat) ::  caldat1
integer julian,mm,id,iyyy
integer jalpha,ja,jb,jc,jd,je
integer igreg

  julian=caldat1.julian
  igreg=2299161
  if(julian.ge.igreg)then
    jalpha=int(((julian-1867216)-0.25)/36524.25)
    ja=julian+1+jalpha-int(0.25*jalpha)
  else
    ja=julian
  endif
  jb=ja+1524
  jc=int(6680.+((jb-2439870)-122.1)/365.25)
  jd=365*jc+int(0.25*jc)
  je=int((jb-jd)/30.6001)
  id=jb-jd-int(30.6001*je)
  mm=je-1
  if(mm.gt.12)mm=mm-12
  iyyy=jc-4715
  if(mm.gt.2)iyyy=iyyy-1
  if(iyyy.le.0)iyyy=iyyy-1
  caldat1.iyyy=iyyy
  caldat1.mm=mm
  caldat1.id=id

end subroutine


!!********************************!!******************************************************* 
! 4 the log:GASDEV 
!!********************************!!******************************************************* 
function gasdev(idum,gset,iset)
use PrecStru
implicit none
integer iset 
real(r8) gasdev, gset,r,fac,v1,v2
integer idum
  r=2
  if (iset.eq.0) then
    do while (r.ge.1.)
      v1=2.*ran(idum)-1.
      v2=2.*ran(idum)-1.
      r=v1**2+v2**2
    end do
    if(r.le.0.) r=0.1
    fac=sqrt(-2.*dlog(r)/r)
    gset=v1*fac
    gasdev=v2*fac
    iset=1
  else
    gasdev=gset
    iset=0
  endif
  return
end function

 

!!********************************!!******************************************************* 
! 3 the log:JULDAY 
!!********************************!!******************************************************* 

function julday(mm,id,iyyy)
implicit none
integer julday,mm,id,iyyy
integer jy,jm,ja
integer igreg
  igreg=15+31*(10+12*1582)
  if (iyyy.eq.0) pause 'there is no year zero.'
  if (iyyy.lt.0) iyyy=iyyy+1
  if (mm.gt.2) then
    jy=iyyy
    jm=mm+1
  else
    jy=iyyy-1
    jm=mm+13
  endif
  julday=int(365.25*jy)+int(30.6001*jm)+id+1720995  
  if (id+31*(mm+12*iyyy).ge.igreg) then
    ja=int(0.01*jy)
    julday=julday+2-ja+int(0.25*ja)
  endif
  return
end

!!********************************!!******************************************************* 
!  2 Calculate the Cumulative Nash-Sutcliffe efficiency  and correlation coefficient 
!!********************************!!******************************************************* 
subroutine CRelEffa(CRE)
use PrecStru
implicit none
type(CRelEffP) ::  CRE
type(CRelEffP) ::  CRE1
integer i,j,k
 
  call CRelEff(CRE)
  allocate(CRE1.yy(CRE.TimeSum))
  allocate(CRE1.yc(CRE.TimeSum))   
  !Sliding
  do i=CRE.TimeSum,1,-1   
    if (i>CRE.aSum )then
    k=CRE.aSum
    else            
    k=i
  end if
    CRE1.yy(i)=CRE.yy(i)
    CRE1.yc(i)=CRE.yc(i)
    do j=i-1,i-k,-1     
      CRE1.yy(i)=CRE1.yy(i)+CRE.yy(j)
      CRE1.yc(i)=CRE1.yy(i)+CRE.yc(j)
    end do
  end do
  CRE1.TimeSum=CRE.TimeSum
  call CRelEff(CRE1)
  CRE.sr=CRE1.r
  CRE.sDC=CRE1.DC
  CRE.sBias=CRE1.Bias
  CRE.sRMSE=CRE1.RMSE

  !cumulative
  do  i=1,CRE.TimeSum/CRE.aSum
    CRE1.yy(i)=0
    CRE1.yc(i)=0
    do  j=1,CRE.aSum     
      CRE1.yy(i)=CRE1.yy(i)+CRE.yy((i-1)*CRE.aSum+j)
      CRE1.yc(i)=CRE1.yc(i)+CRE.yc((i-1)*CRE.aSum+j)
    end do
  end do
  CRE1.TimeSum=CRE.TimeSum/CRE.aSum  
  call CRelEff(CRE1)
  CRE.ar=CRE1.r
  CRE.aDC=CRE1.DC
  CRE.aBias=CRE1.Bias
  CRE.aRMSE=CRE1.RMSE
  deallocate (CRE1.yy) 
  deallocate (CRE1.yc) 
end subroutine

!!********************************!!******************************************************* 
! 1 Calculate the Nash-Sutcliffe efficiency  and correlation coefficient  
!!********************************!!******************************************************* 
subroutine CRelEff(CRE)
use PrecStru
implicit none
type(CRelEffP) ::  CRE
real(r8),pointer ::  yy(:)
real(r8),pointer ::  yc(:) 
integer TimeSum
real(r8) r,DC
real(r8) ay,ac,Sc
real(r8) lcc,lyy,lcy
integer i,j,k
 yy=>CRE.yy
 yc=>CRE.yc
 TimeSum=CRE.TimeSum
!1.mean
 ay=0
 ac=0
 CRE.xm=0.0
 do i=1,TimeSum  
  CRE.xm=CRE.xm+ yc(i)*yy(i)
  ay=ay+yy(i)
  ac=ac+yc(i)
 end do
 CRE.xm=CRE.xm/TimeSum
 ay=ay/TimeSum
 ac=ac/TimeSum
!2.square deviation
 lcc=0
 lyy=0
 lcy=0
 Sc=0
 do i=1,TimeSum  
  lcc=lcc+(yc(i)-ac)*(yc(i)-ac)
  lyy=lyy+(yy(i)-ay)*(yy(i)-ay)
  lcy=lcy+(yc(i)-ac)*(yy(i)-ay)
  Sc=Sc+(yy(i)-yc(i))*(yy(i)-yc(i))
 end do
 if  (lcc*lyy>0.00001) then  
  r=lcy/sqrt(lcc*lyy)  
 else  
  r=0
 end if
 if (lyy==0) lyy=1
 DC=1-Sc/lyy
 CRE.r=r
 CRE.DC=DC
 if (CRE.TimeSum>0 ) CRE.RMSE=sqrt(Sc/CRE.TimeSum)
 if (ay==0) ay=1
 CRE.Bias=ac/ay
end subroutine