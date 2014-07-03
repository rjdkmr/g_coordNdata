/*
 * This file is part of g_coordNdata
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014  Rajendra Kumar
 *
 * g_coordNdata is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_coordNdata is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_coordNdata.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "do_fit.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "index.h"
#include "tpxio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rmpbc.h"

void CopyRightMsg() {

    const char *copyright[] = {
            "                                                                        ",
            "               :-)  g_coordNdata (-:                                     ",
            "                                                                        ",
            "               Author: Rajendra Kumar                                  ",
            "                                                                        ",
            "         Copyright (C) 2014  Rajendra Kumar                             ",
            "                                                                        ",
            "                                                                        ",
            "g_coordNdata is a free software: you can redistribute it and/or modify      ",
            "it under the terms of the GNU General Public License as published by    ",
            "the Free Software Foundation, either version 3 of the License, or       ",
            "(at your option) any later version.                                     ",
            "                                                                        ",
            "g_coordNdata is distributed in the hope that it will be useful,             ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License       ",
            "along with g_coordNdata.  If not, see <http://www.gnu.org/licenses/>.       ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                           :-)  g_coordNdata (-:                       ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<35; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}


int gmx_coordNdata(int argc, char *argv[])
{
	const char *desc[] = {
			"This program can be used to extract coordinates from the GROMACS MD Trajectory."
			" The output format of coordinate data file is as follows:\n",
			"--------------------------------------------------------------------\n",
			"[number of frames]       [number of atoms]                          \n",
			"X11 Y11 Z11 X12 Y12 Z12 X13 Y13 Z13 ... ... ... ... ... X1n Y1n Z1n \n",
			"X21 Y21 Z21 X22 Y22 Z22 X23 Y23 Z23 ... ... ... ... ... X2n Y2n Z2n \n",
			"X31 Y31 Z31 X32 Y32 Z32 X33 Y33 Z33 ... ... ... ... ... X3n Y3n Z3n \n",
			"... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... \n",
			"... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... \n",
			"... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... \n",
			"... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... \n",
			"Xm1 Ym1 Zm1 Xm2 Ym2 Zm2 Xm3 Ym3 Zm3 ... ... .. ... ... Xmn Ymn Zmn  \n",
			"--------------------------------------------------------------------\n",
			"where \"m\" is number of frames and \"n\" is number of atoms. This file",
			"could be used in other program or scripts for MATLAB, R or Python for ",
			"further analysis and coordinates transformation.                 [PAR]",
			"The program could read coordinates data from an external file with the",
			"above described file format to create a GROMACS format MD trajectory. ",
			"This is useful when after analysis and coordinates transformations ",
			"visualization of trajectory is required [PAR]"
	};

	//COMMAND LINE ARGUMENT STUFF
	gmx_bool bFit=TRUE, bM = FALSE, bPBC=TRUE, bTrajWrite=FALSE, bSame= FALSE;
	output_env_t oenv;
	real init_time = 0, dt = 10;

	t_pargs pa[] = {
			  { "-fit", TRUE, etBOOL, {&bFit}, "To fit structure" },
			  { "-mass", FALSE,etBOOL, {&bM}, " Mass weighted fitting" },
			  { "-pbc",  TRUE,  etBOOL, {&bPBC}, "Apply corrections for periodic boundary conditions" },
			  { "-write", FALSE, etBOOL, {&bTrajWrite}, "To write trajectory from coordinates data"},
			  { "-same", FALSE, etBOOL, {&bSame}, "If frame length and time is similar to the input trajectory"},
			  { "-it", FALSE, etREAL, {&init_time}, "initial time, in the output trajectory"},
			  { "-diff_t", FALSE, etREAL, {&dt}, "time difference between each frame for output trajectory"}
	};

	t_filenm   fnm[] = {
			{ efTRX, "-f",   NULL,            ffOPTRD },
			{ efTPS, NULL,   NULL,            ffREAD  },
			{ efNDX, NULL,   NULL,            ffOPTRD },
			{ efDAT, "-in",  "in_coord",      ffOPTRD },
			{ efDAT, "-o",  "out_coord",      ffOPTWR },
			{ efTRX, "-v",  "trajout.xtc",    ffOPTWR }
	};

#define NFILE asize(fnm)
   int npargs;
   int i=0, j=0;
   CopyRightMsg();
   npargs = asize(pa);
   parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE ,
 		             NFILE,fnm,npargs,pa, asize(desc),desc,0,NULL,&oenv);

   //TOPOLOGY RELATED VARIABLE
   char title[STRLEN];
   t_topology	top;
   int		ePBC;
   rvec	*xtop;
   matrix	box;

   //INDEX RELATED VARIABLES
   int isize, nfit;	//Number of index group
   atom_id *index, *ifit;
   char *grpnm,*fitname;

   //TRAJECTORY RELATED VARIABLES
   int natoms, nframe=0; //Number of atoms in the trajectory
   real t, *time, ***x; //Time
   rvec *xread; //Coordinate
   t_trxstatus *status, *out;
   const char *TrajIn, *fnCoordIn;
   gmx_bool bReadTraj=TRUE;

	//Fitting related variable
	real *w_rls;

	TrajIn = opt2fn_null("-f",NFILE,fnm);
	if(TrajIn == NULL)
		bReadTraj = FALSE;

	fnCoordIn = opt2fn_null("-in",NFILE,fnm);
	if(fnCoordIn == NULL){
		if(bTrajWrite)
			gmx_fatal(FARGS,"-write is on but input coordinate data file is missing\n");
	}
	else{
		if(!bTrajWrite)	{
			fprintf(stderr, "WARNING: \"-write\" option was switched off, switching it on\n");
			bTrajWrite = TRUE;
		}
	}


	if(bTrajWrite && bSame && !bReadTraj){
		fprintf(stderr,"WARNING: Initial time = %lf and time step = %lf\n",init_time,dt);
		bSame = FALSE;
	}

	if(!bTrajWrite && !bReadTraj)
		gmx_fatal(FARGS,"Need either trajectory file or coordinate file in the input\n");

   read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,FALSE);

   if (bFit && bReadTraj && !bTrajWrite) {
      printf("\nChoose a group for the least squares fit\n");
      get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nfit,&ifit,&fitname);
      if (nfit < 3)
        gmx_fatal(FARGS,"Need >= 3 points to fit!\n");
    }
   else
	   nfit = 0;

   printf("\nChoose a group for the output\n");
   get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isize,&index,&grpnm);

   if(bFit && bReadTraj && !bTrajWrite)	{
	   snew(w_rls,top.atoms.nr);
	   for(i=0; (i<nfit); i++) {
		   if (bFit && bM)
			   w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;
		   if(bFit && !bM)
			   w_rls[ifit[i]] = 1.0;
	   }
   }

	//Removing PBC from first frame
	gmx_rmpbc_t  gpbc=NULL;
	   if (bPBC) {
	     gpbc = gmx_rmpbc_init(&top.idef,ePBC,top.atoms.nr,box);
	     gmx_rmpbc(gpbc,top.atoms.nr,box,xtop);
	   }

	if (bFit && bReadTraj && !bTrajWrite)
		reset_x(nfit,ifit,top.atoms.nr,NULL,xtop,w_rls);

	if(bReadTraj)	{
		natoms=read_first_x(oenv,&status,TrajIn,&t,&xread,box);
		snew(time,1);
	}

	if(bReadTraj && !bTrajWrite)	{
		snew(x,3);
		for(i=0;i<3;i++)	{
			snew(x[i],isize);
			for(j=0;j<isize;j++)	{
				snew(x[i][j],1);
			}
		}
	}

	if(bReadTraj)	{
		do {
			if(ePBC)
				gmx_rmpbc(gpbc,natoms,box,xread);

			if (bFit && bReadTraj && !bTrajWrite) {
			  reset_x(nfit,ifit,natoms,NULL,xread,w_rls);
			  do_fit(natoms,w_rls,xtop,xread);
			}

			if(bReadTraj && !bTrajWrite){
				for(i=0;i<isize;i++)
				{
					srenew(x[0][i],nframe+1);
					srenew(x[1][i],nframe+1);
					srenew(x[2][i],nframe+1);
					x[0][i][nframe] = xread[index[i]][0];
					x[1][i][nframe] = xread[index[i]][1];
					x[2][i][nframe] = xread[index[i]][2];
					//printf( "%10.3f%10.3f%10.3f",xread[index[i]][0],xread[index[i]][1],xread[index[i]][0]);
					//printf("%10.3f%10.3f%10.3f",x[0][i][nframe],x[1][i][nframe],x[2][i][nframe]);
				}
			}
			srenew(time,nframe+1);
			time[nframe] = t;
			nframe++;
		}
		while(read_next_x(oenv,status,&t,natoms,xread,box));
		sfree(xread);
	}


	if(bReadTraj && !bTrajWrite){
		//Opening coordinate output file
		const char *fnCoordOut;
		FILE *fnmCoordOut;
		fnCoordOut = opt2fn("-o",NFILE,fnm);
		fnmCoordOut = ffopen(fnCoordOut,"w");
		printf("\nWriting coordinate file....\n");
		fprintf(fnmCoordOut,"%10d%10d\n",nframe,isize);
		for(i=0;i<nframe;i++){
			for(j=0;j<isize;j++)
				fprintf(fnmCoordOut, "%10.3f%10.3f%10.3f",x[0][j][i],x[1][j][i],x[2][j][i]);
			fprintf(fnmCoordOut,"\n");
		}
		ffclose(fnmCoordOut);
	}

	if(bTrajWrite)
	{
		const char *fnOutTraj;
		FILE *fnmCoordIn;
		int infr, InAtomnum,d;
		fnmCoordIn = ffopen(fnCoordIn,"r");
		fscanf(fnmCoordIn,"%d%d",&infr,&InAtomnum);
		if(bSame){
			if(((infr != nframe) || (InAtomnum != isize)) || ((infr != nframe) && (InAtomnum != isize)))
				gmx_fatal(FARGS, "Frame number or atom number or both did not match in input Trajectory and input coordinates data\n");
		}

		if(!bReadTraj)	{
			snew(time,infr);
			for(i=0;i<infr;i++)
				time[i] = init_time + (i*dt);
		}
		snew(xread,top.atoms.nr);
		fnOutTraj = opt2fn_null("-v", NFILE, fnm);
		out=open_trx(fnOutTraj,"w");
		printf("\nWriting trajectory file....\n");
	    for(j=0; j<infr; j++) {
	    	for(i=0; i<isize; i++)
	    		for(d=0; d<DIM; d++)
	    			 fscanf(fnmCoordIn,"%g",&xread[index[i]][d]);
	    	write_trx(out,isize,index,&top.atoms,0,time[j],box,xread,NULL,NULL);
	    }
	    close_trx(out);
	}

	fprintf(stdout, "\nThanks for using g_coordNdata!!!\n");
   return 0;
}

int main(int argc, char *argv[])
{
	gmx_coordNdata(argc, argv);
    return 0;
}

