#include <omp.h>


#include "partition.h"
#include "partchunk.h"
#include "escfunction.h"
#include "optics.h"
#include "indicatrix.h"
#include "coords.h"
#include "photon.h"



LinearInterpolation*	Photon::s_oLength       = NULL;
LinearInterpolation*	Photon::s_eLength       = NULL;

Partition*		        Photon::s_oePartition   = NULL;
Partition*		        Photon::s_eoPartition   = NULL;
Partition*		        Photon::s_eePartition   = NULL;

EscFunction*	        Photon::s_oEscFunction	= NULL;
EscFunction*	        Photon::s_eEscFunction	= NULL;

LinearInterpolation*    Photon::s_eChannelProb  = NULL;



void Photon::init(	LinearInterpolation*    oLength_,
                    LinearInterpolation*    eLength_,
                    Partition*              oePartition_,
                    Partition*              eoPartition_,
                    Partition*              eePartition_,
					EscFunction*            oEscFunction_,
                    EscFunction*            eEscFunction_,
                    LinearInterpolation*    eChannelProb_)
{
	s_oLength     = oLength_;
	s_eLength     = eLength_;

	s_oePartition = oePartition_;
	s_eoPartition = eoPartition_;
	s_eePartition = eePartition_;

	s_oEscFunction = oEscFunction_;
	s_eEscFunction = eEscFunction_;

	s_eChannelProb = eChannelProb_;
}


Photon::Photon(RngEngine& rng_engine_, const Vector3& s, const int channel_) :
	pos(0.,0.,0.),
	s_i(s),
	scatterings(0),
	weight(1.),
	fullIntegral(0.),
	channel(channel_),
	rng_engine(rng_engine_),
	oLength(*s_oLength),
	eLength(*s_eLength),
	oePartition(*s_oePartition),
	eoPartition(*s_eoPartition),
	eePartition(*s_eePartition),
	oEscFunction(*s_oEscFunction),
	eEscFunction(*s_eEscFunction),
	eChannelProb(*s_eChannelProb),
	m_chunk(NULL),
	m_knotValues(),
	m_rectValues()
{
	m_knotValues.reserve(std::max(s_oePartition->getMaxKnotsCount(),
	                    std::max(s_eoPartition->getMaxKnotsCount(), s_eePartition->getMaxKnotsCount())));
	m_rectValues.reserve(std::max(s_oePartition->getMaxRectsCount(),
	                    std::max(s_eoPartition->getMaxRectsCount(), s_eePartition->getMaxRectsCount())));
}


void Photon::move()
{
	Float rnd;
    Float theta = symmetrizeTheta(Angle(s_i, Optics::director).theta);

	Float meanFreePath = (channel == Optics::OCHANNEL) ? oLength(theta) : eLength(theta);

    Float c1 = (s_i.z() >= 0) ? 1. : -expm1((pos.z()/s_i.z())/meanFreePath);

	rnd = random();

	Float d = -log1p(-c1*rnd)*meanFreePath;

	pos += d*s_i;
}

void Photon::scatter()
{
    //reduce weight
    
    if (Optics::OCHANNEL == channel)
        weight *= (1. - oEscFunction(acos(s_i.z()), atan2(s_i.y(),s_i.x()), pos.z()));
    else
        weight *= (1. - eEscFunction(acos(s_i.z()), atan2(s_i.y(),s_i.x()), pos.z()));


    //prepare random numbers
    Float randChannel = random();
	Float randRect    = random();
	Float randX       = random();
	Float randY       = random();
	Float randPhi     = random();

    Vector3 nn;
    Matrix3 mtx;
    Angle   a_i;

    createTransformToPartitionCoords(mtx, nn, a_i);


    //select channel
    int newChannel;

    if (Optics::OCHANNEL == channel) {

        newChannel = Optics::ECHANNEL;
        m_chunk = eoPartition.getChunk(a_i.theta);
        calcPartitionValues<IndicatrixOE>(nn);
    }
    else {

        newChannel = randChannel > eChannelProb(symmetrizeTheta(Angle(s_i, Optics::director).theta)) ? 
            Optics::OCHANNEL : Optics::ECHANNEL;

        if (Optics::OCHANNEL == newChannel) {

            m_chunk = eoPartition.getChunk(a_i.theta);
            calcPartitionValues<IndicatrixEO>(nn);
        }
        else {

            m_chunk = eePartition.getChunk(a_i.theta);
            calcPartitionValues<IndicatrixEE>(nn);
        }
    }


    //renormalization  
    randRect    *= fullIntegral;


	//binary search of partition rect
	int rectIdx = 0;
	int first = 0;
	int last = m_chunk->m_rects.size() - 1;

	while (first < last) {

    	rectIdx = (first + last) / 2;
       	if (randRect > m_rectValues[rectIdx]) {

        	first = rectIdx + 1;
	   	}
       	else if (randRect < m_rectValues[rectIdx]) {

        	last = rectIdx - 1;
	   	}
       	else {

           	break;
	   	}
    }

    //adjust point
	Float p, t;
	choosePointInRect(t, p, rectIdx, randX, randY);

	if (randPhi > 0.5)   //choose one of symmetrical cases
		p = 2*M_PI - p;

	Float sintheta = sin(t);
	Vector3 s_s =  Vector3(sintheta*cos(p), sintheta*sin(p), cos(t));

	s_i = invert(mtx)*s_s;

	s_i.normalize();  //to be sure

	channel = newChannel;

	scatterings++;
}

void Photon::createTransformToPartitionCoords(Matrix3& mtx, Vector3& nn, Angle& a_i)
{
	Vector3 v2;
	nn = Optics::director;

	if (Angle(s_i, nn).costheta < 0) {

	    nn = -nn;
    }

    a_i = Angle(s_i, nn);
   
	if (fabs(a_i.sintheta) > kMachineEpsilon) {

		v2 = crossProduct(s_i, nn).normalize();
	}
	else {

		v2 = createSomePerpendicular(s_i).normalize();
	}

	Vector3 v3 = crossProduct(s_i, v2).normalize();

	mtx = createTransformMatrix(v2, v3, s_i);

	nn = mtx*nn;                            //director in s_i -based coordinate system
}

template <class T>
void Photon::calcPartitionValues(const Vector3& nn)
{
    T ind(Vector3(0., 0., 1.), nn);
	
	{
		KnotsVector& knots = m_chunk->m_knots;
		KnotsVector::iterator k;

		m_knotValues.clear();

		Float res;

		for (k = knots.begin(); k != knots.end(); ++k) {

			Float   sintheta = sin(k->x);
			Vector3 s_s      = Vector3(  sintheta*cos(k->y),
									     sintheta*sin(k->y),
									     cos(k->x));

	        res = sintheta*ind(s_s);
			m_knotValues.push_back(res);
		}
	}

	{
    	RectsVector& rects = m_chunk->m_rects;

		RectsVector::iterator i;
		fullIntegral = 0.;

		m_rectValues.clear();

		for (i = rects.begin(); i != rects.end(); ++i) {

			Float rectIntegral = 0.25* (m_knotValues[(*i).tl] +
										m_knotValues[(*i).tr] +
										m_knotValues[(*i).bl] +
										m_knotValues[(*i).br]) *
								(*i).square;

			fullIntegral += rectIntegral;

			m_rectValues.push_back(fullIntegral);
		}
	}
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

			if (x1 >= 0. && x1 < 1.) {

				x = x1;
			}
			else if (x2 >= 0. && x2 < 1.) {

				x = x2;
			}
			else {

				x = 0.5;
				fprintf(stderr, "x out of range, %f\t%f\n", x1, x2);
			}
		}
	}

	{
		Float y1, y2;
		roots = solveQuadric(0.5*(b3 + b4*x), b1 + b2*x, -randY*(b1+b2*x + 0.5*(b3+b4*x)), y1, y2);

		if (roots == 1)
			y = y1;
		else if (roots == 2) {

			if (y1 >= 0. && y1 < 1.) {

				y = y1;
			}
			else if (y2 >= 0. && y2 < 1.) {

				y = y2;
			}
			else {

				y = 0.5;
				fprintf(stderr, "y out of range, %f\t%f\n", y1, y2);
			}
		}
	}

	x *= rect.width;
	y *= rect.height;

	x += m_chunk->m_knots[rect.tl].x;
	y += m_chunk->m_knots[rect.tl].y;

}
