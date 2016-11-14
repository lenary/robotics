#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

void
write_laser_data( FILE *fout,
		  int numLaserv1, int *laserval1,
		  int numLaserv2, int *laserval2 )
{
  int i;
  fprintf( fout, "#LASER %d %d:", numLaserv1, numLaserv2 );
  for (i=0; i<numLaserv1; i++) {
    fprintf( fout, " %d", laserval1[i] );
  }
  for (i=0; i<numLaserv2; i++) {
    fprintf( fout, " %d", laserval2[i] );
  }
  fprintf( fout, "\n" );
}

int
min( int val1, int val2 )
{
  if (val1<=val2)
    return(val1);
  else
    return(val2);
}


int
main( int argc, char *argv[] )
{
  
  FILE   * fpin;
  FILE   * fpout;
  
  char     infile[256];
  char     outfile[256];

  char     line[4096];
  char     command[4096];
  char     dummy[4096];
  time_t   tsec;
  
  float    f1, f2, f3;
  
  double   fsec, cur_posx = 0.0, cur_posy = 0.0, cur_poso = 0.0;
  float    arange;

  int      lasernr, numValues, numValues1, numValues2=0;
  int      laser1[1000];
  int      laser2[1000];
  char   * running;

  unsigned long sec, usec;
  struct tm *actual_date;
  
  int i, rline = 0;
  
  if (argc==3) {
    strncpy(infile,argv[1],256);
    strncpy(outfile,argv[2],256);
    if ((fpin=fopen(infile,"r"))==0) {
      fprintf( stderr, "Can't open input file %s !\n", infile );
      exit(0);
    }
    if ((fpout=fopen(outfile,"w"))==0) {
      fprintf( stderr, "Can't open output file %s !\n", outfile );
      exit(0);
    }

    while (fgets(line,4096,fpin) != NULL) {
      rline++;
      if (rline%100==0)
	fprintf( stderr, "\r-> read line: %d", rline );
      if ((strlen(line)>1) && (sscanf(line,"%s",command)!=0)) {
	if (!strcmp( command, "POS")) {
	  sscanf(line, "%s %ld %ld: %f %f %f",
		 dummy, &sec, &usec, &f1, &f2, &f3 );
	  cur_posx = (double) f1;
	  cur_posy = (double) f2;
	  cur_poso = (double) (90.0-f3);
	  tsec = sec;
	  actual_date = localtime( &tsec );
	  fsec =  (double) (actual_date->tm_sec+(usec/1000000.0));
	  fprintf( fpout, "@SENS %s%d-%s%d-%d %s%d:%s%d:%s%f\n",
		   actual_date->tm_mday<10?"0":"",
		   actual_date->tm_mday,
		   actual_date->tm_mon<10?"0":"",
		   actual_date->tm_mon,
		   actual_date->tm_year,
		   actual_date->tm_hour<10?"0":"",
		   actual_date->tm_hour,
		   actual_date->tm_min<10?"0":"",
		   actual_date->tm_min,
		   fsec<10?"0":"",
		   fsec );
	  fprintf( fpout, "#ROBOT %f %f %f\n",
		   (float) cur_posx, (float) cur_posy, (float) cur_poso );
	  
	} else if (!strcmp( command, "LASER-RANGE")) {
	  sscanf(line, "%s %ld %ld %d %d %f:",
		 dummy, &sec, &usec, &lasernr, &numValues, &arange );
	  running = line;
	  strtok( running, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  if (1 || lasernr==0) {
	    for (i=0; i<numValues; i++) {
	      laser1[i] = atoi( strtok( NULL, " ") );
	    }
	    numValues1 = numValues;
	    /* if you are using 0.5 degrees resolution, this
	       will convert it to 1 degree resolution (using
	       minimum of the both values)! */
	    if (0 && numValues1==361) {
	      for (i=0; i<180; i++) {
		laser1[i] = min( laser1[i*2],laser1[i*2+1] );
	      }
	      numValues1 = 180;
	    }
	    tsec = sec;
	    actual_date = localtime( &tsec );
	    fsec =  (double) (actual_date->tm_sec+(usec/1000000.0));
	    fprintf( fpout, "@SENS %s%d-%s%d-%d %s%d:%s%d:%s%f\n",
		     actual_date->tm_mday<10?"0":"",
		     actual_date->tm_mday,
		     actual_date->tm_mon<10?"0":"",
		     actual_date->tm_mon,
		     1900+actual_date->tm_year,
		     actual_date->tm_hour<10?"0":"",
		     actual_date->tm_hour,
		     actual_date->tm_min<10?"0":"",
		     actual_date->tm_min,
		     fsec<10?"0":"",
		     fsec );
	    write_laser_data( fpout,
			      numValues1, laser1,
			      numValues2, laser2 );
	    numValues2 = 0;
	  } else {
	    for (i=0; i<numValues; i++) {
	      laser2[numValues-1-i] = atoi( strtok( NULL, " ") );
	    }
	    numValues2 = numValues;
	    if (0 && numValues2==361) {
	      for (i=0; i<180; i++) {
		laser2[i] = min( laser2[i*2],laser2[i*2+1] );
	      }
	      numValues2 = 180;
	    }
	  }

	  
	} else if (!strcmp( command, "LASER")) {

	  sscanf(line, "%s %ld %ld %d %d:",
		 dummy, &sec, &usec, &lasernr, &numValues );
	  running = line;
	  strtok( running, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  strtok( NULL, " ");
	  if (lasernr==0) {
	    for (i=0; i<numValues; i++) {
	      laser1[i] = atoi( strtok( NULL, " ") );
	    }
	    numValues1 = numValues;
	    /* if you are using 0.5 degrees resolution, this
	       will convert it to 1 degree resolution (using
	       minimum of the both values)! */
	    if (0 && numValues1==361) {
	      for (i=0; i<180; i++) {
		laser1[i] = min( laser1[i*2],laser1[i*2+1] );
	      }
	      numValues1 = 180;
	    }
	    tsec = sec;
	    actual_date = localtime( &tsec );
	    fsec =  (double) (actual_date->tm_sec+(usec/1000000.0));
	    fprintf( fpout, "@SENS %s%d-%s%d-%d %s%d:%s%d:%s%f\n",
		     actual_date->tm_mday<10?"0":"",
		     actual_date->tm_mday,
		     actual_date->tm_mon<10?"0":"",
		     actual_date->tm_mon,
		     1900+actual_date->tm_year,
		     actual_date->tm_hour<10?"0":"",
		     actual_date->tm_hour,
		     actual_date->tm_min<10?"0":"",
		     actual_date->tm_min,
		     fsec<10?"0":"",
		     fsec );
	    write_laser_data( fpout,
			      numValues1, laser1,
			      numValues2, laser2 );
	    numValues2 = 0;
	  }
	}
      }
    }
    fclose(fpin);
    fprintf( stderr, "\n" );
  } else {
    fprintf( stderr, "usage: %s rec-file script-file\n", argv[0] );
  }
  return(0);
}


