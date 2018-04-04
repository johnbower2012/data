#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;


int main(int argc, char * const argv[]){
	double ERROR_FACTOR=0.04,y;
	int nruns=1024,irun,iread;
	double rdummy1,balnorm,error;
	double bf[19],dely[17];
	char dummy[100];
	char command[100],readfilename[100],writefilename[100];
	FILE *wfptr,*fptr;

		sprintf(writefilename,"results.dat");
		wfptr=fopen(writefilename,"w");

		sprintf(readfilename,"star_pipi.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>0) balnorm+=bf[iread]*0.1;
		}
		fprintf(wfptr,"double balnormpipi %8.6f %8.6f\n",balnorm*0.1,0.15*balnorm*0.1);
		for(iread=1;iread<18;iread++){
			y=(iread+1.5)*0.1;
			error=ERROR_FACTOR*(1.8-y);
			fprintf(wfptr,"double ypipi%d %8.5f %8.5f\n",iread,bf[iread]/balnorm,error);
		}
		fclose(fptr);

		sprintf(readfilename,"star_KK.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>0) balnorm+=bf[iread]*0.1;
		}
		fprintf(wfptr,"double balnormKK %8.6f %8.6f\n",balnorm*0.1,0.15*balnorm*0.1);
		for(iread=1;iread<18;iread++){
			y=(iread+1.5)*0.1;
			error=ERROR_FACTOR*(1.8-y);
			fprintf(wfptr,"double yKK%d %8.5f %8.5f\n",iread,bf[iread]/balnorm,error);
		}
		fclose(fptr);

		sprintf(readfilename,"star_ppbar.dat");
		fptr=fopen(readfilename,"r");
		balnorm=0.0;
		for(iread=0;iread<18;iread++){
			fscanf(fptr,"%lf %lf %lf",&dely[iread],&bf[iread],&rdummy1);
			if(iread>4) balnorm+=bf[iread]*0.1;
		}
		fprintf(wfptr,"double balnormppbar %8.6f %8.6f\n",balnorm*0.1,0.15*balnorm*0.1);
		for(iread=1;iread<18;iread++){
			y=(iread+1.5)*0.1;
			error=ERROR_FACTOR*(1.8-y);
			fprintf(wfptr,"double yppbar%d %8.5f %8.5f\n",iread,bf[iread]/balnorm,error);
		}
		fclose(fptr);

		fclose(wfptr);

	return 0;
}
