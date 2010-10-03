#include <cstdio>
#include <cmath>

#include "common.h"

#include "direction.h"
#include "extlength.h"
#include "photon.h"
#include "indicatrix_app.h"

#include "indicatrix.h"


int main(int /*argc*/, char ** /*argv*/)
{
    IndicatrixApp app;
    app.run();
    
	return 0;
}


void IndicatrixApp::run()
{
	//precalculate q-tree for indicatrix integration

	//calculating extinction lengths...
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
		
	}
}
