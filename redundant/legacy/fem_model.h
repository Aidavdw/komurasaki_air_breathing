/*  A Finite Element Model (FEM) for the double-tapered reed valve
*/

/* Compute the structural characteristics of the sections of the tapered reed valve. N is the number of nodes, which should coincide with the size of x_node, b, A, and I. */
int * compute_solid(const int n_v, const int n_node, double **x_node, const double *x_start, double *b, double *h, double *A, double *I, const double b0, const double b1, const double h0, const double h1, const double l, const double l_fix, const int n_fix, const int n_clamp, const int n_dof_per_node,int *n_active,int *n_inactive)
{
	// The beam is split into two parts: A fixed part (0->l_fix, n_fix nodes), and a freely deforming part.


	// populates the matrix of the x-locations of each node.
	for (int k = 0; k < n_v; ++k) // iters over the amount of nodes
	{
		for (int i = 0; i < n_node; ++i)
		{
			// Equidistantly places over the interval up until the l_fix. Note that this means that the spacing between nodes in the 'fixed' part and the free-moving part are not the same!
			if (i<n_fix)
			{

				x_node[k][i] = x_start[k] + i*(l_fix)/(n_fix-1);
			}
			else
			{
				x_node[k][i] = x_node[k][n_fix-1] + (i-(n_fix-1))*(l)/(n_node-n_fix);
			}
		}
	}


	for (int i = 0; i < n_node; ++i)
	{
		// Geometry of sections referenced by each node
		b[i] = b0 + i*(b1 - b0)/(n_node-1); // Sets the width by reference. Note that this assumes that the reed valves are tapered perfectly.
		if (x_node[0][i] - x_start[0]<= l_fix)
		{
			h[i] = h0;
		}
		else
		{
			h[i] = h0 + (h1 - h0)*(x_node[0][i]-x_start[0]-l_fix)/l;
		}

		// Cross-section and second moment of intertia
		A[i] = h[i]*b[i];
		I[i] = b[i]*pow(h[i],3)/12.0;
	}

	// Deduce DOFs by assuming that half of the fixation length is fixed
	// IMPORTANT FOR GOOD RESPONSE OF THE VALVE
	int dof_index=0;
	// int n_clamp = n_fix-2;
	// printf("%i\n",n_clamp);
	int *dof_vector = malloc((n_node-n_clamp)*n_dof_per_node*sizeof(int));
	for (int i = n_clamp; i < n_node; ++i)
	{
		// if (x_node[0][i]-x_start[0] > l_fix/2.0)
		// {
			for (int j = 0; j < n_dof_per_node; ++j)
			{
				dof_vector[dof_index] = n_dof_per_node*i + j;
				// printf("DOF AT %d IS: %d\n",dof_index,dof_vector[dof_index]);
				dof_index++;
			}
		// }
	}
	// I think this just creates a vector like
	/*
	* [12, 13, 14, 15, 16]
	* 
	*/

	// Returns n_active and n_inactive by reference. This is done directly with the input params though, doesn't have to do anything with the rest of the calculation.
	*n_inactive = n_clamp*n_dof_per_node;
	*n_active = (n_node - n_clamp)*n_dof_per_node;

	return dof_vector;
}

/* Define the cells that mark the beginning and end of each reed valve.
 * Sets *coef and *p_neighbour to the index of the cell that is 'next to it'.
 * 
 */
void build_fem_interface(const int n_node, const int n_interface, const int n_start, const double **x_fluid, const double **xc_fluid, const double *x_fem, int *p_neighbour,  double *coef, const int nghost)
{
	int cur_index = n_start;
	for (int i = 0; i < n_node; ++i)
	{
		// Locate first cell (x(j),x(j+1)) such that x(j)<x_fem(i) and x(j+1)>x_fem(i)
		while ( xc_fluid[cur_index+1][nghost]< x_fem[i] )
		{
			cur_index ++;
		}
		p_neighbour[i] = cur_index;
		// I have no idea what exactly this does. It's related to Florian (2017) eq 3.22, specifically the x/L_elem term, but what it actually represents is too vague.
		coef[i] = (x_fem[i]-xc_fluid[cur_index][nghost])/(xc_fluid[cur_index+1][nghost]-xc_fluid[cur_index][nghost]);
		//printf("NODE %i: NEIGHBOUR %i COEF %f\n", i,cur_index,coef[i]);
	}
}

/* Calculate the force on a tapered plate due to pressure difference (linear pressure distribution assumption) */
double elem_force(double x1, double x2, double b1, double b2, double p1, double p2)
{
	return (p1*b1*(x2-x1) + (b2-b1)*(p2-p1)*(x2-x1)/3.0 + 0.5*(x2-x1)*(p1*(b2-b1)+b1*(p2-p1)));
	// return (x2-x1)*0.35*101325.0*0.5*(b1+b2);
	// printf("%f %f %f %f %f\n",x1,x2,b1,b2,dp);
	// return (x2-x1)*(p1+p2)*0.25*(b1+b2);
}

/* Valve width at a distance x from the point where the tapered shape begins */
double v_width(double b0, double b1, double l, double x, double x0)
{
	// if (x-x0>l)
	// {
	// 	fprintf(stderr,"\nERROR IN VALVE WIDTH CALCULATION: %f is TOO LARGE\n",x-x0); exit(0);
	// }
	return b0 + (b1-b0)/l*(x-x0);
}

// void fem_load(int n_elem, int n_dof, int n_dof_per_node, int n_clamp, int n_ghost, double x_start, int *p_neighbour, int n_interface, double *dp, double **x_fluid, double *x_fem, double w1, double w2, double l, double *f_dof)
// {
// 	// DOF force vector
// 	double f_elem = 0.0;

// 	// Set to zero
// 	for (int i = 0; i < n_dof; ++i)
// 	{
// 		f_dof[i] = 0.0;
// 	}

// 	// For each element
// 	int n1,n2,n_start=p_neighbour[0];
// 	double b1,b2,x1,x2;
// 	for (int i = n_clamp; i < n_elem; ++i)
// 	{
// 		n1 = p_neighbour[i];
// 		n2 = p_neighbour[i+1];
// 		// First case: n1=n2, i.e. the full element is in one single fluid cell. FEM coordinates can be used to compute the load on the element.
// 		if (n1 == n2)
// 		{
// 			x1=x_fem[i];
// 			x2=x_fem[i+1];
// 			b1=v_width(w1,w2,l,x1,x_start);
// 			b2=v_width(w1,w2,l,x2,x_start);
// 			// printf("B0 AND B1: %f %f\n",b1,b2);
// 			// printf("X0 AND X1: %f %f\n",x1-x_start,x2-x_start);
// 			f_elem = elem_force(x1,x2,b1,b2,dp[n1-n_start]);

// 			f_dof[i*n_dof_per_node] += f_elem/2.0;
// 			f_dof[(i+1)*n_dof_per_node] += f_elem/2.0;
// 		}
// 		// General code: the current FEM element is facing several fluid cells
// 		else if (n1<n2)
// 		{
// 			for (int j = n1; j < n2+1; ++j)
// 			{
// 				if (j==n1)
// 				{
// 					x1=x_fem[i];
// 					x2=x_fluid[n1+1][n_ghost];
// 					b1=v_width(w1,w2,l,x1,x_start);
// 					b2=v_width(w1,w2,l,x2,x_start);
// 					// printf("B0 AND B1: %f %f\n",b1,b2);
// 					// printf("X0 AND X1: %f %f\n",x1-x_start,x2-x_start);
// 					f_elem = elem_force(x1,x2,b1,b2,dp[j-n_start]);

// 					f_dof[i*n_dof_per_node] += f_elem/2.0;
// 					f_dof[(i+1)*n_dof_per_node] += f_elem/2.0;

// 				}
// 				else if (j==n2)
// 				{
// 					x1=x_fluid[j][n_ghost];
// 					x2=x_fem[i+1];
// 					b1=v_width(w1,w2,l,x1,x_start);
// 					b2=v_width(w1,w2,l,x2,x_start);
// 					// printf("B0 AND B1: %f %f\n",b1,b2 );
// 					// printf("X0 AND X1: %f %f\n",x1-x_start,x2-x_start);
// 					f_elem = elem_force(x1,x2,b1,b2,dp[j-n_start]);

// 					f_dof[i*n_dof_per_node] += f_elem/2.0;
// 					f_dof[(i+1)*n_dof_per_node] += f_elem/2.0;
// 				}
// 				else if (j>n1 && j<n2)
// 				{
// 					x1=x_fluid[j][n_ghost];
// 					x2=x_fluid[j+1][n_ghost];
// 					b1=v_width(w1,w2,l,x1,x_start);
// 					b2=v_width(w1,w2,l,x2,x_start);
// 					// printf("B0 AND B1: %f %f\n",b1,b2 );
// 					// printf("X0 AND X1: %f %f\n",x1-x_start,x2-x_start);
// 					f_elem = elem_force(x1,x2,b1,b2,dp[j-n_start]);

// 					f_dof[i*n_dof_per_node] += f_elem/2.0;
// 					f_dof[(i+1)*n_dof_per_node] += f_elem/2.0;
// 				}
// 				else
// 				{
// 					fprintf(stderr,"\nERROR: FEM INDEX NOT BETWEEN N1 AND N2: %i\n",j); exit(0);
// 				}

// 				// f_elem = elem_force(x1,x2,b1,b2,dp[n1]);
// 			}
// 		}
// 		// f_elem = elem_force(x_fem[i],x_fem[i+1],b[i],b[i+1],p_fem[i],p_fem[i+1]);

// 		// For pressure distribution case
// 		// f_dof[i*n_dof_per_node] += f_elem/2.0;
// 		// f_dof[(i+1)*n_dof_per_node] += f_elem/2.0;

// 		// For tip load case
// 		// f_dof[(n_elem)*n_dof_per_node] += 0.5*f_elem;
// 		// printf("%d %f\n",i,f_dof[(n_elem+1)*n_dof_per_node]);

// 		if (isnan(f_dof[i*n_dof_per_node]))
// 		{
// 			fprintf(stderr,"\nERROR: COMPUTED FEM LOAD IS NaN\n"); exit(0);
// 			// printf("FEM LOAD IS NAN AT %d\n",i );
// 		}
// 	}

// 	// for (int i = 0; i < n_dof/2; ++i)
// 	// {
// 	// 	printf("F(%i)=%f\n",i,f_dof[2*i]);
// 	// }
// 	// for (int i = 0; i < n_dof; ++i)
// 	// {
// 	// 	printf("LOAD AT NODE %d is %f \n",i,f_dof[i]);
// 	// }
// }

/* Generate the load vector applied to solve the FEM equation (Newmark scheme). Populates f_dof */
void fem_load(const int n_elem, const int n_dof, const int n_clamp, const int n_dof_per_node, const double *p_fem, const double *u_dof, double *f_dof, const double *x_fem, const double *b)
{
	// For each element, compute total pressure-induced load and distribute between nodes
	double f_elem = 0.0;
	double theta_elem = 0.0;
	for (int i = 0; i < n_dof; ++i)
	{
		f_dof[i] = 0.0;
	}
	for (int i = n_clamp; i < n_elem; ++i)
	{
		// f_elem = elem_force(x_fem[i],x_fem[i+1],b[i],b[i+1],10000,10000);
		
		//
		theta_elem = 0.5*(u_dof[i*n_dof_per_node+1]+u_dof[(i+1)*n_dof_per_node+1]);
		f_elem = elem_force(x_fem[i],x_fem[i+1],b[i],b[i+1],p_fem[i],p_fem[i+1])*cos(theta_elem);

		// For distributed case
		f_dof[i*n_dof_per_node] += f_elem/2.0;
		f_dof[(i+1)*n_dof_per_node] += f_elem/2.0;

		// For tip load case
		// f_dof[(n_elem)*n_dof_per_node] += f_elem;
		// printf("%d %f\n",i,f_dof[(n_elem+1)*n_dof_per_node]);

		if (isnan(f_dof[i*n_dof_per_node]))
		{
			fprintf(stderr,"\nERROR: COMPUTED FEM LOAD IS NaN\n"); exit(0);
			// printf("FEM LOAD IS NAN AT %d\n",i );
		}
	}
	// for (int i = 0; i < n_dof; ++i)
	// {
	// 	printf("LOAD AT NODE %d is %f \n",i,f_dof[i]);
	// }
}

/* Apply the aerodynamic damping correction term (should be improved because might not be physical!) */
void fem_flow_damping(const int nElem, int n_dof, const int nClamp, const int nDofPerNode, const double *u1, const double *u2, const double *b, const double *h, const double rho_v, const double f0, const double dt, const double c1, const double c2, const double c3, double *f_dof)
{
	double f_damping = 0.0;
	double eps_eff = 0.0;
	double m_valve = 0.0;
	double dy = 0.0;
	double y = 0.0;

	for (int i = nClamp; i < nElem; ++i)
	{
		m_valve = 0.25*rho_v*(b[i]+b[i+1])*(h[i]+h[i+1]); // Isn't this forgetting the thicknesss?
		dy = 0.5*(u2[(i+1)*nDofPerNode]-u1[(i+1)*nDofPerNode] + (u2[i*nDofPerNode]-u1[i*nDofPerNode]));
		y = 0.5*(u2[i*nDofPerNode]+u2[(i+1)*nDofPerNode]);

		if (dy>= 0.0)
		{
			eps_eff = c1 + c2*dy/dt;
		}
		else if (dy < 0.0)
		{
			eps_eff = c3*y;
			// eps_eff = c3*pow(u2[i*n_dof_per_node],2);
		}

		f_damping = -2.0*f0/m_valve*dy/dt*eps_eff;
		
		f_dof[i*nDofPerNode] += f_damping/2.0;
		f_dof[(i+1)*nDofPerNode] += f_damping/2.0;

		if (isnan(f_dof[i*nDofPerNode]))
		{
			fprintf(stderr,"\nERROR: COMPUTED FEM LOAD IS NaN\n"); exit(0);
		}
	}
}

// void fem_pressure(int n_interf, int index_min, double *x_fem, int ny_inf, int n_cell_p_fem, double **p_in_bot, double **p_in_top, int nghost, double *dp_output)
// {
// 	// double dp_1 = 0.0;
// 	// double dp_2 = 0.0;
// 	for (int i = 0; i < n_interf; ++i)
// 	{
// 		dp_output[i] = 0.0;
// 		// dp_1 = 0.0;
// 		// dp_2 = 0.0;
// 		for (int m = 0; m < n_cell_p_fem; ++m)
// 		{
// 			dp_output[i] += (p_in_top[index_min+i][nghost+m]-p_in_bot[index_min+i][ny_inf-nghost-1-m])/n_cell_p_fem;
// 			// dp_1 += (p_in_top[p_neighbour[i]][nghost+m]-p_in_bot[p_neighbour[i]][ny_inf-nghost-1-m])/n_cell_p_fem;
// 			// dp_2 += (p_in_top[p_neighbour[i]+1][nghost+m]-p_in_bot[p_neighbour[i]+1][ny_inf-nghost-1-m])/n_cell_p_fem;
// 			// dp_1 = 0.35*101325.0;
// 			// dp_2 = 0.35*101325.0;
// 		}
// 		// printf("PRESSURE RATIO AT %i:%f\n",i,dp_output[i]);
// 		// p_output[i] = dp_1 + p_coef[i]*(dp_2-dp_1);
// 		if (isnan(dp_output[i]))
// 		{
// 			fprintf(stderr,"\nERROR: COMPUTED FEM PRESSURE IS NaN AT %d!\n",i);
// 			exit(0);
// 			//p_output[i] = 0.0;
// 		}
// 	}
// }


/* Generate tbe pressure difference that will be used by the FEM method at each element centre by finite difference with its given pressure coefficient.  */
void fem_pressure(const int n_node, int n_interf, int index_min, double *x_fem, const int *p_neighbour, const double *p_coef, const int ny_input, const double **p_in_bot, const double **p_in_top, const int nghost, double *dp_output)
{
	double dp_1 = 0.0;
	double dp_2 = 0.0;
	for (int i = 0; i < n_node; ++i)
	{
		dp_output[i] = 0.0;
		// Gets the difference in pressure between the ambient and cell value for two different x-coordinates. Interpolates the two then, because the node is not exactly in one x-coordinate!
		// the [][nghost+2] implies that it reads (arbitrarily) 2 cells above the ghost cells. Maybe to combat circular dependency and add some diffusion?
		dp_1 = p_in_top[p_neighbour[i]][nghost+2]-p_in_bot[p_neighbour[i]][ny_input-nghost-1-2];
		dp_2 = p_in_top[p_neighbour[i]+1][nghost+2]-p_in_bot[p_neighbour[i]+1][ny_input-nghost-1-2]; 
		// This is the (dp1 + (dp2 - dp1) * x/L_elem) in Florian (2017) eq 3.22
		dp_output[i] = dp_1 + p_coef[i]*(dp_2-dp_1);
		if (isnan(dp_output[i]))
		{
			fprintf(stderr,"\nERROR: COMPUTED FEM PRESSURE IS NaN AT %d: %f %f %f %f!\n",i,p_in_top[p_neighbour[i]][nghost],p_in_bot[p_neighbour[i]][ny_input-nghost-1],p_in_top[p_neighbour[i]+1][nghost],p_in_bot[p_neighbour[i]+1][ny_input-nghost-1]);
			exit(0);
			//p_output[i] = 0.0;
		}
	}
}


/* Assemble Mass matrix */
void build_mass_mat(int N_elem_tot, double **M, double *x_node, double *b, double *h, double *A, double *I, double b0, double b1, double h0, double h1, double l, double rho, int N_dof_per_node)
{
	double l_elem=0, M_elem[4][4], A0=0, A1=0;
	for (int i = 0; i < N_elem_tot; ++i)
	{
		A0 = A[i];
		A1 = A[i+1];
		l_elem = fabs(x_node[i+1] - x_node[i]);

		M_elem[0][0] = l_elem*(10.0*A1 + 3.0*A0)/35.0;
		M_elem[0][1] = pow(l_elem,2)*(15.0*A1 + 7.0*A0)/420.0;
		M_elem[0][2] = 9.0*l_elem*(A1 + A0)/140.0;
		M_elem[0][3] = -pow(l_elem,2)*(7.0*A1 + 6.0*A0)/420.0;

		M_elem[1][0] = M_elem[0][1];
		M_elem[1][1] = pow(l_elem,3)*(3.0*A1 + 5.0*A0)/840.0;
		M_elem[1][2] = -M_elem[0][3];
		M_elem[1][3] = -pow(l_elem,3)*(A0 + A1)/280.0;

		M_elem[2][0] = M_elem[0][2];
		M_elem[2][1] = M_elem[1][2];
		M_elem[2][2] = l_elem*(10.0*A0 + 3.0*A1)/35.0;
		M_elem[2][3] = - pow(l_elem,2)*(15.0*A0 + 7.0*A1)/420.0;

		M_elem[3][0] = M_elem[0][3];
		M_elem[3][1] = M_elem[1][3];
		M_elem[3][2] = M_elem[2][3];
		M_elem[3][3] = M_elem[1][1];

		// printf("M_ELEM:\n");
		for (int k = 0; k < 4; ++k)
		{
			for (int j = 0; j < 4; ++j)
			{
				M[N_dof_per_node*i + k][N_dof_per_node*i + j] = M[N_dof_per_node*i + k][N_dof_per_node*i + j] + rho*M_elem[k][j];
				// printf("M edited at (%d,%d)\n", N_dof_per_node*i + k,N_dof_per_node*i + j);
			}
		}
	}
}

/* Assemble Stiffness matrix */
void build_stiff_mat(int N_elem_tot, double **K, double *x_node, double *b, double *h, double *A, double *I, double b0, double b1, double h0, double h1, double l, double E, int N_dof_per_node)
{
	double l_elem=0, K_elem[4][4]={}, I0=0, I1=0;
	for (int i = 0; i < N_elem_tot; ++i)
	{
		I0 = I[i];
		I1 = I[i+1];
		l_elem = fabs(x_node[i+1] - x_node[i]);

		K_elem[0][0] = 6.0*(I0 + I1)/pow(l_elem,3);
		K_elem[0][1] = 2.0/pow(l_elem,2)*(I0 + 2.0*I1);
		K_elem[0][2] = - K_elem[0][0];
		K_elem[0][3] = 2.0/pow(l_elem,2)*(2.0*I0 + I1);

		K_elem[1][0] = K_elem[0][1];
		K_elem[1][1] = (I0 + 3.0*I1)/l_elem;
		K_elem[1][2] = - K_elem[0][1];
		K_elem[1][3] = (I0 + I1)/l_elem;

		K_elem[2][0] = K_elem[0][2];
		K_elem[2][1] = K_elem[1][2];
		K_elem[2][2] = K_elem[0][0];
		K_elem[2][3] = - K_elem[0][3];

		K_elem[3][0] = K_elem[0][3];
		K_elem[3][1] = K_elem[1][3];
		K_elem[3][2] = K_elem[2][3];
		K_elem[3][3] = (3*I0 + I1)/l_elem;

		for (int k = 0; k < 4; ++k)
		{
			for (int j = 0; j < 4; ++j)
			{
				K[N_dof_per_node*i + k][N_dof_per_node*i + j] = K[N_dof_per_node*i + k][N_dof_per_node*i + j] + E*K_elem[k][j];
			}
		}
	}
}

/* Build Damping matrix based on Rayleigh's damping model. ALPHA and BETA factors (respectively of M and K) should be specified by the user. These coefficients can be known from experiment. */
void build_damp_mat(double N_dof,double **C,double **K,double **M,double alpha,double beta)
{
	for (int i = 0; i < N_dof; ++i)
	{
		for (int j = 0; j < N_dof; ++j)
		{
			C[i][j] = alpha*M[i][j] + beta*K[i][j];
		}
	}
}

/* Computes matrices needed for resolution according to Newmark's scheme */
void build_newmark_mat(int N_dof, double **C, double **K, double **M, double dt, double **R1, double **R2, double **R3)
{
	for (int i = 0; i < N_dof; ++i)
	{
		for (int j = 0; j < N_dof; ++j)
		{
			R1[i][j] = M[i][j]/dt/dt + C[i][j]/2.0/dt + K[i][j]/3.0;
			R2[i][j] = 2.0*M[i][j]/dt/dt - K[i][j]/3.0;
			R3[i][j] = -M[i][j]/dt/dt + C[i][j]/2.0/dt - K[i][j]/3.0;
		}
	}
}

/* Given a symmetric positive definite matrix M (size NxN) (should be checked by user), this function computes the Cholesky decomposition of this matrix and fills a matrix L (that should be allocated and initialized by the user) so that A = L*LT (LT is the transpose matrix of L). The outputed L matrix is a superior triangular matrix. */
void cholesky_decomposition(const int n_act, const double **M,double **L, const int *act_dof)
{

	// A: I think this is a variation on the Cholesky-Crout algorithm?
	
   	for (int i = 0; i < n_act; ++i)
   	{
   		L[i][i] = M[act_dof[i]][act_dof[i]]; // Everything is shifted up-left; the first 'free' node is now the top-left node.
      	for (int k=0; k < i; ++k)
      	{
      		L[i][i] -= L[k][i]*L[k][i];
      	}
      	if (L[i][i] <= 0) 
      	{
	 		fprintf(stderr,"\nERROR: non-positive definite matrix!\n");
	 		printf("\nproblem from %d %.20f\n",i,L[i][i]);
      	}
      	L[i][i] = sqrt(L[i][i]);

      	for (int j = i + 1; j < n_act; ++j) 
      	{
	 		L[i][j] = M[act_dof[i]][act_dof[j]];
	 		for (int k = 0; k < i; ++k)
	 		{
	    		L[i][j] -= L[k][i]*L[k][j];
	 		}
	 		L[i][j] /= L[i][i];
      	}
   	}

   	// printf("Cholesky decomposition: \n");
   	// for (int i = 0; i < n_act; ++i)
   	// {
   	// 	for (int j = 0; j < n_act; ++j)
   	// 	{
   	// 		printf("%g ",L[i][j]);
   	// 	}
   	// 	printf("\n");
   	// }
}

/* The system "Ax = b", with A being symmetric definite positive of size NxN, is solved by using A's Cholesky decomposition. The system is only solved along the dimensions which are activated in the FEM model, i.e. referenced in active_dof. The returned solution is also of size NxN and includes all nodal displacements. The L matrix is lower triangular, and therefore A = L*LT. */
void cholesky_solve(const int n_dof,const int n_active,const double **L, const double *load, const int *active_dof, double *u_dof)
{
	
	
	double sum;
	// double *solution = malloc(n_dof*sizeof(double));

	// Initialize solution
	for (int i = 0; i < n_dof; ++i)
	{
		// solution[i] = 0.0;
		u_dof[i] = 0.0;
	}

	// printf("Check that L is an inferior triangular matrix: \n");
	// for (int i = 0; i < n_active; ++i)
	// {
	// 	for (int j = 0; j < n_active; ++j)
	// 	{
	// 		printf("%f ",L[i][j]);
	// 	}
	// 	printf("\n");
	// }

	// First step : Forward substitution starting from first index
	for (int i = 0; i < n_active; ++i)
	{
		sum = 0.0;
		for (int j = 0; j < i; ++j)
		{
			sum += L[i][j]*u_dof[active_dof[j]];
		}
		u_dof[active_dof[i]] = (load[active_dof[i]] - sum)/L[i][i];
	}

	// Second step : Backward subsitution starting from last index
	for (int i = n_active - 1; i > -1; --i)
	{
		sum = 0.0;
		for (int j = i + 1; j < n_active; ++j)
		{
			sum += L[j][i]*u_dof[active_dof[j]];
		}
		u_dof[active_dof[i]] = (u_dof[active_dof[i]] - sum)/L[i][i];
	}

	// for (int i = 0; i < n_dof; ++i)
	// {
	// 	printf("U_DOF(%d): %f\n",i,u_dof[i] );
	// }

	// Replace the computed solution in the global vector of DOFs

	for (int i = 0; i < n_dof; ++i)
	{
		if (isnan(u_dof[i]))
		{
			fprintf(stderr,"\nERROR: CHOLESKY SOLUTION IS NaN!\n");
			exit(0);
		}
	// 	printf("Solution at DOF %d for load %.10f is : %.10f\n",i,load[i],u_dof[i]);
	}
}

void newmark_solve(const int n_dof, const int n_active, const double **LR1, const double **R2, const double **R3, const double *load, const int *active_dof, const double *u0, const double *u1, double *u2)
{
	double sum;
	double b[n_dof];

	// Initialize solution and right-hand array, memorizing current state
	for (int i = 0; i < n_dof; ++i)
	{
		// u0[i] = u1[i];
		// u1[i] = u2[i];
		u2[i] = 0.0;
		b[i] = load[i];
	}

	// Compute right-hand term of "A1*U(n+1) = F(n) + A2*U(n) + A3*U(n-1)"
	for (int i = 0; i < n_dof; ++i)
	{
		for (int j = 0; j < n_dof; ++j)
		{
			b[i] += R2[i][j]*u1[j] + R3[i][j]*u0[j];
		}
	}

	// First step : Forward substitution starting from first index
	for (int i = 0; i < n_active; ++i)
	{
		sum = 0.0;
		for (int j = 0; j < i; ++j)
		{
			sum += LR1[i][j]*u2[active_dof[j]];
		}
		u2[active_dof[i]] = (b[active_dof[i]] - sum)/LR1[i][i];
	}

	// Second step : Backward subsitution starting from last index
	for (int i = n_active - 1; i > -1; --i)
	{
		sum = 0.0;
		for (int j = i + 1; j < n_active; ++j)
		{
			sum += LR1[j][i]*u2[active_dof[j]];
		}
		u2[active_dof[i]] = (u2[active_dof[i]] - sum)/LR1[i][i];
	}

	// Replace the computed solution in the global vector of DOFs
	for (int i = 0; i < n_dof; ++i)
	{
		if (isnan(u2[i]))
		{
			fprintf(stderr,"\nERROR: CHOLESKY SOLUTION IS NaN!\n");
			exit(0);
		}
	}

	// for (int i = 0; i < n_dof; ++i)
	// {
	// 	printf("U_DOF(%d): %f\n",i,u2[i] );
	// }

}

/* Provided an input matrix of size NxN, this function computes its transpose matrix as a pointer allocated by the user beforehand. */
void transpose(int n,double **M, double **MT)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			MT[i][j] = M[j][i];
		}
	}
}
// sets the deflection of the valve; y_fem
void update_valve(const int n_node, const int n_dof_per_node,const double *u0,const double *u1,double *u2, double *y_fem)
{
	if (u2[(n_node-1)*n_dof_per_node] <= 0.012)
	{
		for (int i = 0; i < n_node; ++i)
		{
			if (u2[i*n_dof_per_node] <= 0.0)
			{
				for (int j = 0; j < (i+1)*n_dof_per_node; ++j)
				{
					u2[j] = 0.0;
				}
			}
			y_fem[i] = u2[i*n_dof_per_node];
			// printf("%i %f\n", i,y_fem[i]);
			// printf("Deflection at %d is %f\n",i,y_fem[i]);
		}	
	}
	else
	{
		for (int i = 0; i < n_node; ++i)
		{
			y_fem[i] = u1[i*n_dof_per_node];
		}
	}
}
/* The solution vector of the FEM problem is used along with theoretical quadrature functions in order to derive the exact new Y-wise coordinate of a given fluid node. By convention, this function only returns the deflection with regard to the standstill position. A positive displacement corresponds to a motion of the valve towards the inside of the tube (i.e. plenum pressure is higher than thruster pressure). */
double interpolate_deflection(double x_out,int n_elem,int ndof,int nghost,int n_dof_per_node,double *x_node, double *valve_dof_solution,double xstart,double l_t)
{
	int index=0.0;
	double slratio,l;
	double N[4];
	double y_out=0.0;

	// Error case if outside range of FEM nodes
	if (x_out-xstart < 0.0 || x_out-xstart > l_t)
	{
		fprintf(stderr,"\nERROR: Cell out of FEM range: X=%f should be between %f and %f!\n",x_out-xstart, xstart,l_t);
	}

	// Locate relevant element for interpolation
    while (index < n_elem && !(x_node[index]+xstart <= x_out && x_node[index + 1]+xstart >= x_out))
    {
    	++index;
    }

	// Define x/L ratio in the current element
	l = fabs(x_node[index+1]-x_node[index]);
	slratio = (x_out-x_node[index]-xstart)/l;

    // Compute quadrature functions and deduce interpolated coordinate
	N[0] = 1.0 - 3.0*pow(slratio,2) + 2.0*pow(slratio,3);
	N[1] = l*(slratio - 2.0*pow(slratio,2) + pow(slratio,3));
	N[2] = 3.0*pow(slratio,2) - 2.0*pow(slratio,3);
	N[3] = l*(-pow(slratio,2) + pow(slratio,3));

	for (int k = 0; k < 2*n_dof_per_node; ++k)
	{
		y_out = y_out - N[k]*valve_dof_solution[n_dof_per_node*index + k];
	}

    // printf("Deflection at solid element %d (X=%f and S/L=%f) : %f.\n",index,x_out-xstart,slratio,y_out);

	return y_out;
}


/* End of file */