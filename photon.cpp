#include <omp.h>

#include "extlength.h"
#include "indicatrix.h"
#include "partition.h"
#include "optics.h"
#include "coords.h"
#include "photon.h"


Float Photon::probs[kThetaIterations*kPhiIterations] = {0};
ExtLength* Photon::s_length       = NULL;
Partition* Photon::s_partition    = NULL;


using namespace std::tr1;

mt19937 Photon::rng_core = mt19937();
uniform_real<Float> Photon::dist = uniform_real<Float> (0., 1.); 

variate_generator<mt19937, uniform_real<Float> > Photon::rng =
		variate_generator<mt19937, uniform_real<Float> > (rng_core, dist);




void Photon::init(	ExtLength* length_,
					Partition* partition_,
					unsigned long seed_)
{
	s_length    = length_;
	s_partition = partition_;

	Photon::rng_core.seed(seed_);
}


Photon::Photon() :
	pos(0.,0.,0.),
	s_i(0., 0., 1.),
	a_i(Angle(s_i, Optics::n)),
	scatterings(0),
	weight(1.),
	fullIntegral(0.),
	length(*s_length),
	partition(*s_partition)
{
}


void Photon::move()
{
	Float rnd;
	#pragma omp critical
	{
		rnd = 1. - rng();
	}

	Float d   = -log(rnd) * length(Angle(s_i, Optics::n));
	//fprintf(stderr, "d_pos: %.17f\t%.17f\t%.17f\n", d*s_i.x(), d*s_i.y(), d*s_i.z());

	pos += d*s_i;
}

void Photon::scatter()
{
	//coordinate system
	Vector3 v2;

	if (fabs(Angle(s_i, Optics::n).sintheta) > kMachineEpsilon) {

		v2 = crossProduct(s_i, Optics::n).normalize();
	}
	else {

		v2 = createSomePerpendicular(s_i).normalize();
	}

	Vector3 v3 = crossProduct(s_i, v2).normalize();
	Vector3 v1 = Vector3(1., 0., 0.);

	Matrix3 mtx = createTransformMatrix(s_i, v2, v3);

	Vector3 nn = mtx*Optics::n;

	Angle a_i = Angle(v1, nn);
	Vector3 k_i = Optics::ke(v1, a_i);

	Indicatrix ind = Indicatrix(v1, nn);

	//compute integrals
	{
		KnotsVector& knots = partition.m_knots;
		KnotsVector::iterator i;	
		for (i = knots.begin(); i != knots.end(); ++i) {

			Float   sintheta = sin(i->x);
			Vector3 s_s      = Vector3(  cos(i->x),
									     sintheta*sin(i->y),
									     sintheta*cos(i->y));
	
			i->val = sintheta*ind(s_s);
		}
	}

	RectsVector& rects = partition.m_rects;

	{
		RectsVector::iterator i;

		for (i = rects.begin(); i != rects.end(); ++i) {

			fullIntegral += i->integral();
			i->val = fullIntegral;
		}
	}

	//random value to choose rect

	Float randRect;
	Float randX;
	Float randY;
	Float randPhi;

	#pragma omp critical
	{
		randRect = rng()*fullIntegral;
		randX = rng();
		randY = rng();
		randPhi = rng();
	}


	//binary search
	int rectIdx = 0;
	//int first = 0;
	//int last = rects.size() - 1;


/*	while (first <= last) {

    	int mid = (first + last) / 2;
       	if (randRect > rects[mid].val) {

        	first = mid + 1;
	   	}
       	else if (randRect < rects[mid].val) {

        	last = mid - 1; 
	   	}
       	else {

        	rectIdx = mid;     
           	break;
	   	}
    }
*/
	while (randRect > rects[rectIdx].val) {

		rectIdx++;
	}


	Float phi_s, theta_s;
	//fprintf(stderr, "x=%f\ty=%f\tw=%f\th=%f\n", (*Rect::s_knots)[rects[rectIdx].tl].x, (*Rect::s_knots)[rects[rectIdx].tl].y, rects[rectIdx].width, rects[rectIdx].height);
	rects[rectIdx].choosePointInRect(theta_s, phi_s, randX, randY);
	//fprintf(stderr, "%f\t%f\n", phi_s, theta_s);

	if (randPhi > 0.5)   //choose one of symmetrical cases
		phi_s = 2*M_PI - phi_s;

	Float sintheta = sin(theta_s);
	Vector3 s_s =  Vector3(  cos(theta_s),
						     sintheta*sin(phi_s),
						     sintheta*cos(phi_s));
	//fprintf(stderr, "%f\t%f\n", sintheta, phi_s, cos(phi_s));

	//mtx = invert(mtx);

	s_i = invert(mtx)*s_s;

	/*Matrix3 testMtx = invert(mtx);
	testMtx = mtx*testMtx;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			fprintf(stderr, "%f\t", testMtx(i,j));

	fprintf(stderr, "\n");

	fprintf(stderr, "%f\t%f\t%f\n", s_s.x(), s_s.y(), s_s.z());
	fprintf(stderr, "%f\t%f\t%f\n", s_i.x(), s_i.y(), s_i.z());
*/

	scatterings++;
}
