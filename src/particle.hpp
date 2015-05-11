#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_


#include <ofMain.h>
#include "ofxColorMap.h"

class Controler;

//TODO: add window parameters:
class Particle
{
public:

    Particle();
    Particle(ofVec2f pos, ofVec2f vel);


    void draw(ofxColorMap &cm);
    void update(Controler &controler);

    void onResize(int w, int h);

    ofVec2f position;
    ofVec2f velocity;
    ofPolyline trail;
};

#endif //_PARTICLE_HPP_
