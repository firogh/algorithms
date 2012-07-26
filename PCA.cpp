/**
 *principal component analysis(PCA)
 *always be used to reduce dimensionality
 *@author songtianyi630@163.com
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
using namespace std;

#define SIG_EXP 30 //if the exponent of double eeigen value less than SIGN_EXP,
                    //corresponding eigen vector will be discarded
#define EPS            0.000001
#define ITERATION      60

void mySwap(int *a,int *b){

    if(a == b)return;
    assert(a != b);
	(*a) = (*a) ^ (*b);
	(*b) = (*a) ^ (*b);
	(*a) = (*a) ^ (*b);

    return ;
}

/**
 *quickSort ,[left,right]
 *sort the index(id) instead of original array
 *so we needn't exchange the eigen vector
 *@author songtianyi630@163.com
 */
void quickSort(int left,int right,double array[],int id[]){
    int i = left,j = right;double x = array[id[(left+right)/2]];
    do{
        while(array[id[i]] < x)i++;
        while(array[id[j]] > x)j--;
        if(i <= j) {mySwap(&id[i++],&id[j--]);}
    }while(i < j);//i >= j
    if(i < right)quickSort(i,right,array,id);
    if(j > left)quickSort(left,j,array,id);
}


/**
 *Householder_Tri_Symetry_Diagonal
 *@author (copy from internet)
 */
void Householder_Tri_Symetry_Diagonal(double a[], int n, double q[], double b[], double c[])
{
	int i, j, k, u;
	double h, f, g, h2;

	for ( i = 0; i <= n-1; i++ )
	{
		for ( j = 0; j <= n-1; j++ )
		{
			u = i * n + j;
			q[u] = a[u];
		}
	}
	for ( i = n-1; i >= 1; i-- )
	{
		h = 0.0;
		if ( i > 1 )
		{
			for ( k = 0; k <= i-1; k++ )
			{
				u = i * n + k;
				h = h + q[u] * q[u];
			}
		}
		if ( h + 1.0 == 1.0 )
		{
			c[i] = 0.0;
			if ( i == 1 ) c[i] = q[i*n+i-1];
			b[i] = 0.0;
		}
		else
		{
			c[i] = sqrt( h );
			u = i * n + i - 1;
			if ( q[u] > 0.0 ) c[i] = -c[i];
			h = h - q[u] * c[i];
			q[u] = q[u] - c[i];
			f = 0.0;
			for ( j = 0; j <= i - 1; j++ )
			{
				q[j*n+i] = q[i*n+j] / h;
				g = 0.0;
				for ( k = 0; k <= j; k++ )
					g = g + q[j*n+k] * q[i*n+k];
				if ( j + 1 <= i-1 )
					for ( k = j+1; k <= i-1; k++ )
						g = g + q[k*n+j] * q[i*n+k];
				c[j] = g / h;
				f = f + g * q[j*n+i];
			}
			h2 = f / ( h + h );
			for ( j = 0; j <= i-1; j++ )
			{
				f = q[i*n+j];
				g = c[j] - h2 * f;
				c[j] = g;
				for ( k = 0; k <= j; k++ )
				{
					u = j * n + k;
					q[u] = q[u] - f * c[k] - g * q[i*n+k];
				}
			}
			b[i] = h;
		}
	}
	for ( i = 0; i <= n-2; i++ )
	{
		c[i] = c[i+1];
	}
	c[n-1] = 0.0;
	b[0] = 0.0;
	for ( i = 0; i <= n-1; i++ )
	{
		if ( ( b[i] != 0.0 ) && ( i-1 >= 0 ) )
		{
			for ( j = 0; j <= i-1; j++ )
			{
				g = 0.0;
				for ( k = 0; k <= i-1; k++ )
					g = g + q[i*n+k] * q[k*n+j];
				for ( k = 0; k <= i-1; k++ )
				{
					u = k * n + j;
					q[u] = q[u] - g * q[k*n+i];
				}
			}
		}
		u = i * n + i;
		b[i] = q[u]; q[u] = 1.0;
		if ( i - 1 >= 0 )
		{
			for ( j = 0; j <= i - 1; j++ )
			{
				q[i*n+j] = 0.0; q[j*n+i] = 0.0;
			}
		}
	}

	return;
}


/**
 *Tri_Symmetry_Diagonal_Eigenvector
 *@return if the algorithm fail to find the eigen vector, it will return -1
 *@author (copy from internet)
 */

int Tri_Symmetry_Diagonal_Eigenvector(int n, double b[], double c[], double q[], double eps, int l)
{
	int i, j, k, m, it, u, v;
	double d, f, h, g, p, r, e, s;

	c[n-1] = 0.0; d = 0.0; f = 0.0;
	for ( j = 0; j <= n-1; j++ )
	{
		it = 0;
		h = eps * ( fabs( b[j] ) + fabs( c[j] ) );
		if ( h > d )
		{
			d = h;
		}
		m = j;
		while ( ( m <= n-1 ) && ( fabs( c[m] ) > d ) )
		{
			m = m+1;
		}
		if ( m != j )
		{
			do
			{
				if ( it == l )
				{
#ifdef ALGO_DEBUG
					printf( "fail\n" );
#endif
					return( -1 );
				}
				it = it + 1;
				g = b[j];
				p = ( b[j+1] - g ) / ( 2.0 * c[j] );
				r = sqrt( p * p + 1.0 );
				if ( p >= 0.0 )
					b[j] = c[j] / ( p + r );
				else
					b[j] = c[j] / ( p - r );
				h = g - b[j];
				for ( i = j+1; i <= n-1; i++ )
					b[i] = b[i] - h;
				f = f + h; p = b[m]; e = 1.0; s = 0.0;
				for ( i = m-1; i >= j; i-- )
				{
					g = e * c[i]; h = e * p;
					if ( fabs( p ) >= fabs( c[i] ) )
					{
						e = c[i] / p; r = sqrt( e * e + 1.0 );
						c[i+1] = s * p * r; s = e / r; e = 1.0 / r;
					}
					else
					{
						e = p / c[i]; r = sqrt( e * e + 1.0 );
						c[i+1] = s * c[i] * r;
						s = 1.0 / r; e = e / r;
					}
					p = e * b[i] - s * g;
					b[i+1] = h + s * ( e * g + s * b[i] );
					for ( k = 0; k <= n-1; k++ )
					{
						u = k * n + i + 1; v = u - 1;
						h = q[u]; q[u] = s * q[v] + e * h;
						q[v] = e * q[v] - s * h;
					}
				}
				c[j] = s * p; b[j] = e * p;
			}
			while ( fabs( c[j] ) > d );
		}
		b[j] = b[j] + f;
	}
	for ( i = 0; i <= n-1; i++ )
	{
		k = i; p = b[i];
		if ( i+1 <= n-1 )
		{
			j = i+1;
			while ( ( j <= n-1 ) && ( b[j] <= p ) )
			{
				k = j; p = b[j]; j = j+1;
			}
		}
		if ( k != i )
		{
			b[k] = b[i]; b[i] = p;
			for ( j = 0; j <= n-1; j++ )
			{
				u = j * n + i; v = j * n + k;
				p = q[u]; q[u] = q[v]; q[v] = p;
			}
		}
	}

	return( 1 );
}



/**
 *SymmetricRealMatrixEigen
 *@author (copy from internet)
 */
int calEigenVector(double CovMatrix[], int n, double Eigen[], double EigenVector[])
{
	int k;
	double * subDiagonal;

	subDiagonal = ( double * )malloc( sizeof( double ) * n );
	Householder_Tri_Symetry_Diagonal( CovMatrix, n, EigenVector, Eigen, subDiagonal );
	k = Tri_Symmetry_Diagonal_Eigenvector( n, Eigen, subDiagonal, EigenVector, EPS, ITERATION );
	free( subDiagonal );

	return( k );
}

inline int getExponent(double v)
{
    assert(sizeof(short) == 2);
    short *t = ((short *)&v);
    t += 3;
    short tt = *t;
    tt = tt & (32767);
    tt >>= 4;
    tt -= 1023;
    return tt;
}

void copyEigenVector(int col,int dim,double *retained_eig,const double *eigen_vector,const int *eigen_id)
{
    //col*col -> col*dim
    int it = 0;
    for(int i = 0;i < col;i++)
    {
        for(int j = col-1,c = 0;c < dim;c++,j--)
        {
            retained_eig[it++] = eigen_vector[i*col + eigen_id[j]];
        }
    }
    assert(it == col*dim);
}

/**
 *[ii][kk] X [kk][jj] = [ii][jj]
 */
void matrixMulti(double *rs,const double *a,const double *b,int ii,int kk,int jj)
{
    double *tmp = new double[ii*jj];
    memset(tmp,0,sizeof(double)*ii*jj);

    for(int i = 0;i < ii;i++)
    {
        for(int j = 0;j < jj;j++)
        {
            for(int k = 0;k < kk;k++)
            {
                tmp[i*jj + j] += a[i*kk+k]*b[k*jj+j];
            }
        }
    }
    memcpy(rs,tmp,sizeof(double)*ii*jj);
    delete [] tmp;
}

/**
 *matrix transposition
 */
void matrixTransposition(double *m,int row,int col)
{
    //row*col -> col*row;
    double *tmp = new double[row*col];
    int it = 0;
    for(int i = 0;i < col;i++)
    {
        for(int j = 0;j < row;j++)
        {
            tmp[it++] = m[j*col+i];
        }
    }
    assert(it == row*col);
    memcpy(m,tmp,sizeof(double)*it);
    delete [] tmp;
}

/**
 *principal component analysis and dimensionality reducing
 */
void PCA(double *mat,const int row,const int col)
{

    double *exp_value    = new double[col]; // expectation value of each column
	double *eigen_vector = new double[col*col];
	double *retained_eig = new double[col*col];
    double *eigen_value  = new double[col];
	int    *eigen_id     = new int[col];

    //calcualte expection value
    memset(exp_value,0,sizeof(double)*col);
    for(int j = 0;j < col;j++)
    {
        for(int i = 0;i < row;i++)
        {
            exp_value[j] += mat[i*col + j];
        }
        exp_value[j] /= row;
    }

    //calculate covariance matrix
    //symetric matrix cov(x,y) = cov(y,x)
    assert(row > 0 && col > 0);
    double *cova_mat = new double[ col * col];
	for(int j=0; j < col; j++)
	{
		for(int  k= j; k < col; k++)
		{
			double lMjk = 0;
			for(int i = 0; i <row; i++)
			{
				lMjk += (mat[i*col + j ] - exp_value[j])*(mat[ i*col + k ] - exp_value[k]);
			}
			cova_mat[j * col + k] = lMjk / (row - 1);
			cova_mat[k * col + j] = cova_mat[j*col+k];
		}
	}

    //calculate eigen vector
	int indi = calEigenVector(cova_mat,col,eigen_value,eigen_vector);
	if(indi == -1) return;

	delete [] cova_mat;

	//index sorting
    for(int i = 0;i < col;i++) eigen_id[i] = i;
	quickSort(0,col-1,eigen_value,eigen_id);

	int re_dim = col;
	for(int i = 0;i < col;i++)
	{
	    if(getExponent(eigen_value[eigen_id[i]]) < SIG_EXP)
	    {
	        //discard the eigen vector
	        re_dim = col - i - 1;
	    }
	    else break;
	}
	assert(re_dim > 0);

	copyEigenVector(col,re_dim,retained_eig,eigen_vector,eigen_id);

    //free memory
    delete [] eigen_vector; eigen_vector    = 0;
	delete [] eigen_value;  eigen_value     = 0;
	delete [] eigen_id;     eigen_id        = 0;
	delete [] exp_value;    exp_value       = 0;

    //p = retained_eig
    //N*M X M*P = N*P
    //row * col X col*dim = row * dim

    matrixMulti(mat,mat,retained_eig,row,col,re_dim);

    /**
    restore mat
    //N*P X P*M = N*M
    //row*dim X dim*col = row*col
    //rotate matrix ,col*dim -> dim*col;
    matrixTransposition(retained_eig,col,re_dim);

    matrixMulti(mat,mat,retained_eig,row,dim,col);
    */

    delete [] retained_eig; retained_eig = 0;
}


int main(){return 0;}
