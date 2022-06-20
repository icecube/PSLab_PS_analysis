#include "GeneralFC_Functions.h"

#include <stdio.h>


double SetDecGeneralFC(GeneralFC& gfc, char *mainDir, int nSigLoad,
		       double sigMax, double decDeg)
{
  char buffer[1000];
  sprintf(buffer,"%s/range.txt",mainDir);
  FILE *fp = fopen(buffer,"r");
  if (!fp) { 
    printf("ERROR: %s not found.  Exit.\n",buffer);
    exit(1);
  }


  double decMinDeg;
  double decMaxDeg;
  int nDecBins;
  fscanf(fp,"%lg %lg %d", &decMinDeg, &decMaxDeg, &nDecBins);
  fclose(fp);

  // find dec file which is closest to the requested decDeg
  int nBin;

  if (decDeg<decMinDeg) { 
    nBin = 0;
    printf("WARNING: dec=%lg < decMin=%lg\n",decDeg,decMinDeg);
  }
  else if (decDeg>=decMaxDeg) { 
    nBin = nDecBins-1;
    printf("WARNING: dec=%lg > decMin=%lg\n",decDeg,decMinDeg);
  }
  else {
    double fraction = (decDeg-decMinDeg)/(decMaxDeg-decMinDeg);
    nBin = int( fraction * nDecBins );
  }

  // nBin ranges from 0 to nDecBins-1
  // now calculate dec at center of bin, which is name of FC directory
  double srcDecDeg = decMinDeg + 
    (decMaxDeg-decMinDeg) * (nBin+0.5) / nDecBins;

  sprintf(buffer,"%s/dec_%08.4f/",mainDir,srcDecDeg);

  gfc.Load(buffer,nSigLoad,sigMax);

  return srcDecDeg;
}

void SetNameGeneralFC(GeneralFC& gfc, char *mainDir, int nSigLoad,
               double sigMax, char *nameDir)
{
  char buffer[1000];

  sprintf(buffer,"%s/%s/",mainDir,nameDir);

  gfc.Load(buffer,nSigLoad,sigMax);

}

  
