#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;

double GetError(double bf,double dely){
	double error=fabs(0.01*(1.0-dely/1.8));
	error+=0.075*fabs(bf);
	return error;
}

int main(int argc, char * const argv[]){
	int nruns=1024,irun,iread;
	double rdummy1,balnorm,error;
	double bf[19],dely[17];
	char dummy[100];
	char command[100],readfilename[100],writefilename[100];
	FILE *wfptr,*fptr;

		sprintf(writefilename,"experimental_results.dat");
		wfptr=fopen(writefilename,"w");

		sprintf(readfilename,"star_pipi.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>0) balnorm+=0.1*bf[iread];
		}
		//fprintf(wfptr,"Bpipi_norm %8.6f %8.6f\n",balnorm,0.15);
		for(iread=1;iread<18;iread++){
			error=GetError(bf[iread],dely[iread]);
			fprintf(wfptr,"Bpipi_ybin%d %8.5f %8.5f\n",iread,bf[iread],error);
		}
		fclose(fptr);

		sprintf(readfilename,"star_KK.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>0) balnorm+=0.1*bf[iread];
		}
		//fprintf(wfptr,"BKK_norm %8.6f %8.6f\n",balnorm,0.15*balnorm);
		for(iread=1;iread<18;iread++){
			error=GetError(bf[iread],dely[iread]);
			fprintf(wfptr,"BKK_ybin%d %8.5f %8.5f\n",iread,bf[iread],error);
		}
		fclose(fptr);

		sprintf(readfilename,"star_ppbar.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>4) balnorm+=0.1*bf[iread];
		}
		//fprintf(wfptr,"Bppbar_norm %8.6f %8.6f\n",balnorm,0.15*balnorm);
		for(iread=1;iread<18;iread++){
			error=GetError(bf[iread],dely[iread]);
			fprintf(wfptr,"Bppbar_ybin%d %8.5f %8.5f\n",iread,bf[iread],error);
		}
		fclose(fptr);
		
		sprintf(readfilename,"star_pK.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>4) balnorm+=0.1*bf[iread];
		}
		//fprintf(wfptr,"BpK_norm %8.6f %8.6f\n",balnorm,0.15*balnorm);
		for(iread=1;iread<18;iread++){
			error=GetError(bf[iread],dely[iread]);
			fprintf(wfptr,"BpK_ybin%d %8.5f %8.5f\n",iread,bf[iread],error);
		}
		fclose(fptr);

		fclose(wfptr);

	return 0;
}
