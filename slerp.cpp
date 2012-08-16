/**
 *@author http://www.oschina.net/code/snippet_65674_781
 */
void CSlerp::slerp( float result[4], float starting[4], float ending[4], float t)
{
	float cosa = starting[0]*ending[0] + starting[1]*ending[1] + starting[2]*ending[2] + starting[3]*ending[3];
	if ( cosa < 0.0f ) {
		ending[0] = -ending[0];
		ending[1] = -ending[1];
		ending[2] = -ending[2];
		ending[3] = -ending[3];
		cosa = -cosa;
	}
	float k0, k1;
	if ( cosa > 0.9999f ) {
		k0 = 1.0f - t;
		k1 = t;
	}
	else {
		float sina = sqrt( 1.0f - cosa*cosa );
		float a = atan2( sina, cosa );
		float invSina = 1.0f / sina;
		k0 = sin((1.0f - t)*a) * invSina;
		k1 = sin(t*a) * invSina;
	}
	result[0] = starting[0]*k0 + ending[0]*k1;
	result[1] = starting[1]*k0 + ending[1]*k1;
	result[2] = starting[2]*k0 + ending[2]*k1;
	result[3] = starting[3]*k0 + ending[3]*k1;
}