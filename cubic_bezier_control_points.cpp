/**
 *@brief given a set of points, first and last point would be
         endpoints(p0,p3), we need find two control points(p1,p2)
         to do cubic bezier fitting.

         the curve that derived from p0,p1,p2p3 with cubic bezier interpolation,
         would traverse the set of points
 *@author songtianyi630@163.com
 */
#include <cassert>
#include <cmath>
using namespace std;
struct _3DPoint
{
    double x,y,z;
    _3DPoint()
    {
        x = y = z = 0;
    }
    _3DPoint(double xx,double yy,double zz)
    {
        x = xx; y = yy; z = zz;
    }
};

//double*_3Dpoint
_3DPoint operator * (const double v,const _3DPoint &p)
{
    _3DPoint r;
    r.x = p.x*v; r.y = p.y*v; r.z = p.z*v;
    return r;
}
_3DPoint operator - (const _3DPoint &a,const _3DPoint &b)
{
    return _3DPoint(a.x - b.x,a.y - b.y,a.z - b.z);
}
void operator += (_3DPoint &a, const _3DPoint &b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
//_3Dpoint / double
_3DPoint operator / (const _3DPoint &a,const double v)
{
    _3DPoint r;
    r.x = a.x / v;
    r.y = a.y / v;
    r.z = a.z / v;
    return r;
}

double euclideanDistance(const _3DPoint &a,const _3DPoint &b)
{
    return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
}
//p0 p1 p2 p3
void cubicBezierControlPoint(_3DPoint *q,int n,_3DPoint &p1,_3DPoint &p2)
{
    assert(n > 4);//it is nonsense when n less than or equal to 4
    _3DPoint p0 = q[0],p3 = q[n-1];

    double *t = new double[n];

    //calculate chord length
    t[0] = 0;
    for(int i = 1;i < n;i++)
    {
        //accumulate
        t[i] = t[i-1] + euclideanDistance(q[i],q[i-1]);
    }
    //calculate time
    for(int i = 1;i < n;i++)
    {
        t[i] /= t[n-1];
    }

    double A1 = 0,A2 = 0, A12 = 0;
    _3DPoint C1,C2;//0 0 0

    for(int i = 0;i < n;i++)
    {
        double ti2   = t[i]*t[i] , ti3   = ti2*t[i] , ti4   = ti3*t[i];

        double _ti   = 1-t[i] , _ti2  = _ti*_ti , _ti3  = _ti2*_ti, _ti4  = _ti3*_ti;

        A1 += ti2 * _ti4;
        A2 += ti4 * _ti2;
        A12+= ti3 * _ti3;

        _3DPoint tmp;
        tmp = (q[i] - _ti3*p0 - ti3*p3);
        C1 += ( 3*t[i]*_ti2*tmp );
        C2 += ( 3*ti2*_ti*tmp );
    }
    delete [] t; t = 0;
    //P1 = (A2C1?A12C2)/(A1A2?A12A12)
    //P2 = (A1C2?A12C1)/(A1A2?A12A12)
    double base = A1*A2 - A12*A12;
    p1 = (A2*C1 - A12*C2)/base;
    p2 = (A1*C2 - A12*C1)/base;
    assert(t == 0);
}

int main(){return 0;}
