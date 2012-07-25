/**
   Piecewise cubic bezier curve as defined by Adobe in Postscript
   The two end points are p0 and p3
   Their associated control points are p1 and p2
   p(t) = a t^3 + b t^2 + c t + p0
   or
   p(t) = p0 * (1-t)^3 + p1 * 3*t(t-1)^2 + p2 * 3*(1-t)*t^2 + p3 * t^3
 *@author  http://local.wasp.uwa.edu.au/~pbourke/geometry/bezier/cubicbezier.html
 *	 edited by songtianyi630@163.com
 */
struct _3DPoint
{
    double x,y,z;
};
_3DPoint CubicBezier(_3DPoint p0,_3DPoint p1,_3DPoint p2,_3DPoint p3,double mu)
{
   _3DPoint a,b,c,p;

   c.x = 3 * (p1.x - p0.x);
   c.y = 3 * (p1.y - p0.y);
   c.z = 3 * (p1.z - p0.z);
   b.x = 3 * (p2.x - p1.x) - c.x;
   b.y = 3 * (p2.y - p1.y) - c.y;
   b.z = 3 * (p2.z - p1.z) - c.z;
   a.x = p3.x - p0.x - c.x - b.x;
   a.y = p3.y - p0.y - c.y - b.y;
   a.z = p3.z - p0.z - c.z - b.z;

   p.x = a.x * mu * mu * mu + b.x * mu * mu + c.x * mu + p0.x;
   p.y = a.y * mu * mu * mu + b.y * mu * mu + c.y * mu + p0.y;
   p.z = a.z * mu * mu * mu + b.z * mu * mu + c.z * mu + p0.z;

   return(p);
}