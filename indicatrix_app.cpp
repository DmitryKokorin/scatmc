#include <cstdio>
#include <cmath>
#include <memory.h>

#include "common.h"

#include "direction.h"
#include "extlength.h"
#include "photon.h"
#include "indicatrix_app.h"

#include "indicatrix.h"


//test
#include "partition.h"





int main(int /*argc*/, char ** /*argv*/)
{
    IndicatrixApp app;
    app.run();
    
	return 0;
}


void IndicatrixApp::run()
{
	/*//calculating extinction lengths...
		ExtLength length;

	//scattering...
	Photon::init(&length); 

	int cnt = 0;
	
	for (int i = 0; i < 100; ++i) {
		
		Photon ph;
		
		while (ph.pos.z() >= 0 && ph.scatterings < 1000) {

			ph.move();
			ph.scatter();
		}
			
		if (ph.pos.z() < 0) {
			cnt++;
		}
		
		printf("Catched: %d (%d scatterings)\n", cnt, ph.scatterings);
		
	}*/


	//	Float data[Partition::size][Partition::size];
	Float ** data = allocate_2d_array<Float>(Partition::size, Partition::size);

	//memset(&data[0][0], 0, sizeof(data[0][0])*Partition::size*Partition::size);
	Float phiStep   = 2*M_PI/Partition::size;
	Float thetaStep =   M_PI/Partition::size;

	printf("Precalculating greed...\n");
	
	Indicatrix ind(Direction(0.5*M_PI, 0.));
	for (int i = 0; i < Partition::size; ++i)
		for (int j = 0; j < Partition::size; ++j) {

			Direction d = Direction(i*thetaStep, j*phiStep);
			data[j][i] = ind(d)*d.sintheta;
		}


	printf("creating partition...\n");

	Partition p;
	p.setData(data, (M_PI/Partition::size) * (2*M_PI/Partition::size));

	p.refine();

	printf("%d\n", p.rectCount);

	free_2d_array(data);
}
