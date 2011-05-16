#include <omp.h>

#include "freepath.h"
#include "indicatrix.h"
#include "partition.h"
#include "partchunk.h"
#include "optics.h"
#include "coords.h"
#include "photon.h"


FreePath* Photon::s_length        = NULL;
Partition* Photon::s_partition   = NULL;


using namespace std::tr1;

mt19937 Photon::rng_core = mt19937();
uniform_real<Float> Photon::dist = uniform_real<Float> (0., 1.); 

variate_generator<mt19937, uniform_real<Float> > Photon::rng =
		variate_generator<mt19937, uniform_real<Float> > (rng_core, dist);




void Photon::init(	FreePath* length_,
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
	m_chunk(NULL),
	m_knotValues(),
	m_rectValues(),
	m_knotEscValues()
{
	m_knotValues.reserve(s_partition->getMaxKnotsCount());
	m_rectValues.reserve(s_partition->getMaxRectsCount());
	m_knotEscValues.reserve(s_partition->getMaxKnotsCount());
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
	//coordinate system
	Vector3 v2;
	Vector3 nn = Optics::n;

	if (Angle(s_i, nn).costheta < 0)
	    nn = -nn;

	if (fabs(Angle(s_i, nn).sintheta) > kMachineEpsilon) {

		v2 = crossProduct(s_i, nn).normalize();
	}
	else {

		v2 = createSomePerpendicular(s_i).normalize();
	}

	Vector3 v3 = crossProduct(s_i, v2).normalize();
	Vector3 v1 = Vector3(1., 0., 0.);

	Matrix3 mtx = createTransformMatrix(s_i, v2, v3);

	nn = mtx*nn;                            //director in k_i based coordinate system
	Vector3 oz = mtx*Vector3(0., 0., 1.);   //z axis in k_i based coordinate system

	Angle a_i = Angle(v1, nn);
	Vector3 k_i = Optics::ke(v1, a_i);

	Indicatrix ind = Indicatrix(v1, nn);
	m_chunk = partition.getChunk(a_i.theta);

	//compute integrals
	{
		KnotsVector& knots = m_chunk->m_knots;
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

#if 1

    Float esc = 0;

    {
        Indicatrix ind1 = Indicatrix(s_i, Optics::n);
        Float phi, theta, dist;
       
        for (int i = 0; i < 1000; ++i) {

            phi = 2.*M_PI/1000 * i;

            for (int j = 0; j < 1000; ++j) {

                theta = M_PI/1000 * j;
                dist = pos.z() /  cos(theta);
                Vector3 s1 = Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
                esc += ind1(s1)*sin(theta);
            }
        }

        esc *= 2.*M_PI * M_PI / 1000 / 1000;
    }
#endif

	RectsVector& rects = m_chunk->m_rects;

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
		}
	}

    fprintf(stderr, "esc=%.17e\tapprox=%.17e\t%.17e\n", esc, fullIntegral, esc/fullIntegral);
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

	Rect& rect = m_chunk->m_rects[rectNum];


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

	x += m_chunk->m_knots[rect.tl].x;
	y += m_chunk->m_knots[rect.tl].y;
}
