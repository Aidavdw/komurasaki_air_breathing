
/*  
List of functions that are related to export of simulation data. Export format can be set by the user and is .DAT by default. This can easily be changed in the future.
*/


/* Create a subfolder. */
void make_dir(char* name)
{
    int check;
    char *dirname = name;
    DIR* dir = opendir(dirname);
    if (!dir)
    {
        check = mkdir(dirname,0777);
        if (check)
        {
            fprintf(stderr,"\nERROR CREATING DIRECTORY");exit(0);
            // exit(1);  
        } 
    }
    else closedir(dir);
}


/* Export time and a double type number at the current time step in a file of specified extension. 
    data: double type number to append to existing file
    time: Time at which data is appended
    filename: Name of file in which data is appended
    extension: File extension (typically .dat)
*/
void append_data_to_file(double data,double time,char *filename,char *extension)
{
    char fullname[50];
    strcpy(fullname,filename);
    strcat(fullname,extension);
    FILE *fp;

    if (time==0.0)
    {
        fp = fopen(fullname,"w");
        fprintf(fp, "Time ");
        fprintf(fp,"%s\n",filename);
    }
    else
    {
        fp = fopen(fullname,"a");
    }

    // Data is written in the file in the following format: "time value"
    fprintf(fp,"%.7g %.12g\n",time,data);
    fclose(fp);
}


/* Append time and a double type vector at the end of an existing file */
void append_multidata_to_file(int n_data, double *data, double time, char *filename, char *extension)
{
    char fullname[50];
    strcpy(fullname,filename);
    strcat(fullname,extension);
    FILE *fp;

    if (time==0.0)
    {
        fp = fopen(fullname,"w");
        fprintf(fp, "Time ");
        for (int i = 0; i < n_data; ++i)
        {
            fprintf(fp,"valve_%i ",i);
        }
        fprintf(fp, "\n");
    }
    else
    {
        fp = fopen(fullname,"a");
    }

    // Data is written in the file in the following format: "time value1 value2 ... valueN"
    fprintf(fp, "%.7g ",time);
    for (int i = 0; i < n_data; ++i)
    {
        fprintf(fp,"%.12g ",data[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
}


/* Export function for fluid data. This function distinguishes data between time steps, domains, and parameters (all in separate folders and files).
    Ndom: Number of domains
    Nx[],Ny[]: Sizes of domains stored in lists
    x,y: Nodal coordinates of all nodes
    xc,yc: Centre-cell coordinates of cell cells
    rho,u,v,p,E,T,H: Variables to export
    time: Time at which data is exported
    outputFolderName: Name of output folder
    extension: Same as before
    format: Format of time for directory names
*/
void export_fluid_data(int Ndomain,int *nx,int *ny,double ***x,double ***y,double ***xc,double ***yc,double ***rho,double ***u,double ***v,double ***p,double ***E,double ***T,double ***H,double time,char* outputFolderName,char* extension,char* format)
{    
    char folderName[20],timeStr[10];
    sprintf(timeStr,format,time);
    strcpy(folderName,outputFolderName);
    make_dir(folderName);
    strcat(folderName,"/");
    strcat(folderName,timeStr);
    make_dir(folderName);
    strcat(folderName,"/"); 

    char *var[] = { "x","y","xc","yc","rho","u","v","p","E","T","H" };

    #pragma omp parallel for num_threads(Ndomain)
    for (int k = 0; k < Ndomain; ++k)
    {
        FILE *fp[11];
        char curDomain[10],tempName[50],domainName[50];

        for (int i = 0; i < 11; ++i)
        {
            sprintf(curDomain,"%d",k);
            strcpy(domainName,folderName);
            strcat(domainName,curDomain);
            make_dir(domainName);
            strcat(domainName,"/");

            strcpy(tempName,domainName);
            strcat(tempName,var[i]);
            strcat(tempName,extension);
            fp[i] = fopen(tempName,"w");
        }

        #pragma omp parallel sections
        {
            #pragma omp section
            for (int i = 2; i < nx[k]-1; ++i)
            {
                for (int j = 2; j < ny[k]-1; ++j)
                {
                    fprintf(fp[0],"%g ",x[k][i][j]);
                    fprintf(fp[1],"%g ",y[k][i][j]);
                }   
                fprintf(fp[0],"\n"); 
                fprintf(fp[1],"\n");    
            }

            #pragma omp section
            for (int i = 2; i < nx[k]-2; ++i)
            {
                for (int j = 2; j < ny[k]-2; ++j)
                {   
                    fprintf(fp[2],"%g ",xc[k][i][j]);
                    fprintf(fp[3],"%g ",yc[k][i][j]);
                    fprintf(fp[4],"%.12g ",rho[k][i][j]);
                    fprintf(fp[5],"%.12g ",u[k][i][j]);
                    fprintf(fp[6],"%.12g ",v[k][i][j]);
                    fprintf(fp[7],"%.12g ",p[k][i][j]);
                    fprintf(fp[8],"%.12g ",E[k][i][j]);
                    fprintf(fp[9],"%.12g ",T[k][i][j]);
                    fprintf(fp[10],"%.12g ",H[k][i][j]);
                }

                fprintf(fp[2],"\n");
                fprintf(fp[3],"\n");
                fprintf(fp[4],"\n");
                fprintf(fp[5],"\n");
                fprintf(fp[6],"\n");
                fprintf(fp[7],"\n");
                fprintf(fp[8],"\n");
                fprintf(fp[9],"\n");
                fprintf(fp[10],"\n");
            }
        }  

        for (int i = 0; i < 11; ++i)
        {
            fclose(fp[i]);
        }
    }
}


/* Export function for the FEM model (solid) data. Variables at nodal sections, global stiffness and mass matrices, as well as degrees of freedom (DOF) are all exported in the output folder of the current simulation. */
void export_fem_data(int N_elem,int Ndof,double *x_node,double *b,double *h,double *A,double *I,int *DOF,double **K,double **M,double **C,char *filename_sections, char *filename_K, char *filename_M, char *filename_C, char* extension)
{
    char filename[20];
    strcpy(filename,filename_sections);
    strcat(filename,extension);
    FILE *file1 = fopen(filename,"w");

    strcpy(filename,filename_M);
    strcat(filename,extension);
    FILE *file2 = fopen(filename,"w");

    strcpy(filename,filename_K);
    strcat(filename,extension);
    FILE *file3 = fopen(filename,"w");

    strcpy(filename,filename_C);
    strcat(filename,extension);
    FILE *file4 = fopen(filename,"w");

    for (int i = 0; i < N_elem+1; ++i)
    {
        fprintf(file1,"%.10g %.20g %.20g %.20g %.20g\n",x_node[i],b[i],h[i],A[i],I[i]);
    }
    fclose(file1);

    for (int i = 0; i < Ndof; ++i)
    {
        for (int j = 0; j < Ndof; ++j)
        {
            fprintf(file2,"%.20g ",M[i][j]);
            fprintf(file3,"%.20g ",K[i][j]);
            fprintf(file4,"%.20g ",C[i][j]);
        }
        fprintf(file2,"\n");
        fprintf(file3,"\n");
        fprintf(file4,"\n");
    }
    fclose(file2);
    fclose(file3);
    fclose(file4);
}


/* This function exports the parameters that are necessary for dispay and post-processing under Python. Please modify this function freely to take more parameters into account. */
void export_parameters(int Ndom, double Tsim, double dt, double *xstart, double *ystart, double *xlength, double *ylength, int n_valve, int n_fem, int *nx, int *ny, int nghost, int n_export, char *format, char *filename)
{
    FILE *fp=fopen(filename,"w");
    if(fp == NULL)
    {
        sprintf("ERROR: File '",filename,"' could not be opened for writing!");
    }
    else
    {
        fprintf(fp,"%d\n", Ndom);    // Number of domains
        fprintf(fp,"%g\n", Tsim);    // Simulation time
        fprintf(fp,"%g\n", dt);      // Time step
        fprintf(fp,"%d\n", n_valve); // Number of reed valves
        fprintf(fp,"%d\n", n_fem);   // Number of finite elements per valve

        // And for each domain...
        for (int i = 0; i < Ndom; ++i)
        {
            fprintf(fp,"%g ",xstart[i]); // Starting X location of domain
        }
        fprintf(fp,"\n");
        for (int i = 0; i < Ndom; ++i)
        {
            fprintf(fp,"%g ",ystart[i]); // Starting Y location of domain
        }
        fprintf(fp,"\n");
        for (int i = 0; i < Ndom; ++i)
        {
            fprintf(fp,"%g ",xlength[i]); // X-wise length of domain
        }
        fprintf(fp,"\n");
        for (int i = 0; i < Ndom; ++i)
        {
            fprintf(fp,"%g ",ylength[i]); // Y-wise length of domain
        }
        fprintf(fp,"\n");
        for (int i = 0; i < Ndom; ++i)
        {
            fprintf(fp,"%d ",nx[i]-2*nghost); // X-wise number of cells in domain
        }
        fprintf(fp,"\n");
        for (int i = 0; i < Ndom; ++i)
        {
            fprintf(fp,"%d ",ny[i]-2*nghost); // Y-wise number of cells in domain
        }
        fprintf(fp,"\n");

        fprintf(fp,"%d\n", n_export);
        fprintf(fp,"%s",format);
    }
    fclose(fp);
}

/* Export the position of each of the nodes of a given FEM reed petal. Coordinates are referenced with regard to the center axis of the rocket. */
void export_valve_data(int n_valve, double n_node,double **x,double **y, double radius, double time,char* outputFolderName,char* extension,char* format,int nghost)
{    
    char folderName[20],tempName[30],timeStr[10];
    sprintf(timeStr,format,time);
    strcpy(folderName,outputFolderName);
    make_dir(folderName);
    strcat(folderName,"/");
    strcat(folderName,timeStr);
    make_dir(folderName);
    strcat(folderName,"/"); 

    FILE *fp;


    for (int k = 0; k < n_valve; ++k)
    {
        char var[20];
        sprintf(var,"valve_%i%s",k,extension);
        strcpy(tempName,folderName);
        strcat(tempName,var);
        fp = fopen(tempName,"w");

        for (int i = 0; i < n_node; ++i)
        {
            fprintf(fp, "%g %g\n", x[k][i], radius - y[k][i]);
        }
        fclose(fp);
    }
}

/* End of file */