
bool WriteFluxToFile(FILE *fp, double eMinGeV, double eMaxGeV, 
		     double index, double flux, double nEv,
		     char* fluxString) {
  fprintf(fp, "%lg %lg %lg %lg %lg\n",
	  eMinGeV, eMaxGeV, index, flux, nEv);
  fprintf(fp, "%s\n",fluxString);
}

