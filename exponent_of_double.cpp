
/**
 *@brief get the exponent of double value
        double  Sign bit: 1 bit
                Exponent width: 11 bits
                Significand precision: 53 bits (52 explicitly stored)
 *@author songtianyi630@163.com
 */
inline int exponentOf(const double v)
{
    assert(sizeof(short) == 2);
    short *t = ((short *)&v);
    t += 3;
    short t_v = *t;
    t_v = t_v & (32767);t_v >>= 4;t_v -= 1023;
    return t_v;
}
int main(){return 0;}
