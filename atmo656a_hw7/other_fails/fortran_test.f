
c     PROGRAM GETPHASE
  
      COMMON / SIZDIS / NEWSD,         RGV,           SIGMAG,
     2                  RFRS,          RFIS,          RFRC,  
     3                  RFIC,          ALAMB,         IPHASE
 
      COMMON / STUFF /  OMEGA,         ASY,           EXT,   
     2                  SCAT,          QEXT,          QSCAT, 
     3                  QBS
      DIMENSION         COSPHI(181),       SCTPHS(181)
  
      do rr = 0.01, 1 ,0.01
            write(*,*) rr
      end do

      END
