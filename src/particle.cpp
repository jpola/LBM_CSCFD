#include "particle.hpp"
#include "interpolation.hpp"

Particle::Particle()
{

}

Particle::Particle(ofVec2f pos, ofVec2f vel) :
    position(pos), velocity(vel)
{

}

void Particle::draw(ofxColorMap& cm)
{
    static float pmax = -1000;
    static float pmag_scale = 1.0;

    float mag = velocity.length();

    if (mag > pmax)
    {
        pmax = mag;
        pmag_scale = float(NCOLORS) / pmax;
    }

    ofSetColor(cm.use( mag * pmag_scale));
    ofDrawCircle(position, 6);
    trail.draw();
}

void Particle::update(Controler &controler)
{
    float maxx = ofGetWindowWidth() - 5;
    float minx = 5;

    float miny = 5;
    float maxy = ofGetWindowHeight() - 5;

    for(int i = 0; i < 25; i++)
    {

        if (position.x > maxx)
        {
            trail.clear();
            position.x = minx;
        }
        if (position.x < minx)
        {
            trail.clear();
            position.x = maxx;
        }

        if (position.y > maxy)
        {
            trail.clear();
            position.y = miny;
        }

        if (position.y < miny)
        {
            trail.clear();
            position.y = maxy;
        }


//        position.x = ofWrap(position.x, 1, 1280 - 10);
//        position.y = ofWrap(position.y, 1, 720 - 10);
        velocity = -interpolate(controler, position);
        position += velocity*100;

        trail.addVertex(position);

        if (trail.size() > 35) trail.getVertices().erase(trail.getVertices().begin());
    }
    trail.simplify();
}

void Particle::onResize(int w, int h)
{
    position.set(w/2,h/2);
}


