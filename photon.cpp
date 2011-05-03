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
	fullEscIntegral(0.),
	length(*s_length),
	partition(*s_partition),
	m_knotValues(),
	m_rectValues(),
	m_knotEscValues()//,
//	m_rectEscValues()
{
	m_knotValues.reserve(s_partition->m_knots.size());
	m_rectValues.reserve(s_partition->m_rects.size());

	m_knotEscValues.reserve(s_partition->m_knots.size());
//    m_rectEscValues.reserve(s_partition->m_knots.size());
}


void Photon::move()
{
	Float rnd;
	Float extLength = length(Angle(s_i, Optics::n));

    Float c1 = (s_i.z() >= 0) ? 1. : -expm1(pos.z()/s_i.z()/extLength);

	#pragma omp critical
	{
		rnd = rng();
	}

	Float d = -log1p(-c1*rnd)*extLength;

	pos += d*s_i;
}

void Photon::scatter()
{
 /*   //roughly estimate escape function
    Float esc = 0.;

    {
        const Float MIN_PHI = 0.;
        const Float MAX_PHI = 2*M_PI;
        const Float MIN_THETA = 0.5*M_PI + 1e-3;
        const Float MAX_THETA = M_PI - 1e-3;

        const int PHI_ITERATIONS = 10;
        const int THETA_ITERATIONS = 10;

        const Float PHI_STEP = (MAX_PHI - MIN_PHI) / PHI_ITERATIONS;
        const Float THETA_STEP = (MAX_THETA - MIN_THETA) / THETA_ITERATIONS;

        Float theta, phi;

        Indicatrix ind = Indicatrix(s_i, Optics::n);

        for (int i = 0; i < PHI_ITERATIONS; ++i) {
            
            phi = MIN_PHI + i*PHI_STEP;
           
            for (int j = 0; j < THETA_ITERATIONS; ++j) {

                theta = MIN_THETA + j*THETA_STEP;

                Vector3 direction = Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

                Angle angle = Angle(direction, Optics::n);

                esc += ind(direction)*sin(theta)*exp(pos.z()/length(angle)/direction.z());

                theta += THETA_STEP;        
            }
        }

        esc *= PHI_STEP*THETA_STEP;
    }
*/

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

	Vector3 nn = mtx*Optics::n;            //director in k_i based coordinate system
	Vector3 oz = mtx*Vector3(0., 0., 1.);  //z axis in k_i based coordinate system

	Angle a_i = Angle(v1, nn);
	Vector3 k_i = Optics::ke(v1, a_i);

	Indicatrix ind = Indicatrix(v1, nn);

	//compute integrals
	{
		KnotsVector& knots = partition.m_knots;
		KnotsVector::iterator k;

		m_knotValues.clear();
		m_knotEscValues.clear();

		Float res;
		Float dist;

		for (k = knots.begin(); k != knots.end(); ++k) {

			Float   sintheta = sin(k->x);
			Vector3 s_s      = Vector3(  cos(k->x),
									     sintheta*sin(k->y),
									     sintheta*cos(k->y));
	
	        res = sintheta*ind(s_s);
			m_knotValues.push_back(res);

			Angle sz_angle  = Angle(s_s, oz);

            if (sz_angle.costheta >= -kMachineEpsilon) {

                m_knotEscValues.push_back(0.);
            }
            else {

                dist =  pos.z() / sz_angle.costheta;   
                m_knotEscValues.push_back(res*exp(dist/length(Angle(s_s, nn))));
            }
		}
	}

	
	RectsVector& rects = partition.m_rects;

	{
		RectsVector::iterator i;

		m_rectValues.clear();

		for (i = rects.begin(); i != rects.end(); ++i) {

			Float rectIntegral = 0.25* (m_knotValues[(*i).tl] +
										m_knotValues[(*i).tr] +
										m_knotValues[(*i).bl] +
										m_knotValues[(*i).br]) * 
								(*i).square;

			fullIntegral += rectIntegral;

            Float rectEscIntegral = 0.25 * (m_knotEscValues[(*i).tl] +
                                            m_knotEscValues[(*i).tr] +
                                            m_knotEscValues[(*i).bl] + 
                                            m_knotEscValues[(*i).br]) *
                                (*i).square;

            fullEscIntegral += rectEscIntegral;

			m_rectValues.push_back(fullIntegral);
		//	m_rectEscValues.push_back(fullEscIntegral);
		}
	}


    //calc weight
    weight *= 1. - fullEscIntegral/fullIntegral;


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
	while (randRect > m_rectValues[rectIdx]) {

		rectIdx++;
	}


	Float p, t;
	choosePointInRect(t, p, rectIdx, randX, randY);

	if (randPhi > 0.5)   //choose one of symmetrical cases
		p = 2*M_PI - p;

	Float sintheta = sin(t);
	Vector3 s_s =  Vector3(  cos(t),
						     sintheta*sin(p),
						     sintheta*cos(p));


	s_i = invert(mtx)*s_s;

	scatterings++;
}

void Photon::choosePointInRect(Float& x, Float& y, const int rectNum, const Float randX, const Float randY)
{

	Rect& rect = partition.m_rects[rectNum];


	Float b1 = m_knotValues[rect.tl];
	Float b2 = m_knotValues[rect.tr] - b1;
	Float b3 = m_knotValues[rect.bl] - b1;
	Float b4 = b1 - m_knotValues[rect.tr] - m_knotValues[rect.bl] + m_knotValues[rect.br];

	int roots;

	{
		Float x1, x2;
		roots = solveQuadric(b2 + 0.5*b4, 0.5*(b1 + 0.5*b3), -randX*(b2+0.5*b4+0.5*b1+0.25*b3), x1, x2);

		if (roots == 1)
			x = x1;
		else if (roots == 2) {

			if (x1 >= 0. && x1 < 1.)
				x = x1;
			else if (x2 >= 0. && x2 < 1.)
				x = x2;
			else
				fprintf(stderr, "x out of range, %f\t%f\n", x1, x2);
		}
	}

	{
		Float y1, y2;
		roots = solveQuadric(0.5*(b3 + b4*x), b1 + b2*x, -randY*(b1+b2*x + 0.5*(b3+b4*x)), y1, y2);

		if (roots == 1)
			y = y1;
		else if (roots == 2) {

			if (y1 >= 0. && y1 < 1.)
				y = y1;
			else if (y2 >= 0. && y2 < 1.)
				y = y2;
			else
				fprintf(stderr, "y out of range, %f\t%f\n", y1, y2);
		}
	}

	x *= rect.width;
	y *= rect.height;

	x += partition.m_knots[rect.tl].x;
	y += partition.m_knots[rect.tl].y;
}

