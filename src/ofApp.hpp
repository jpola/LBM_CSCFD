#pragma once
#ifndef _OFAPP_HPP_
#define _OFAPP_HPP_

#include <ofMain.h>
#include <ofxGui.h>

#include "lbmController.hpp"
#include "colormap/ofxColorMap.h"
#include "particle.hpp"

class ofApp : public ofBaseApp
{
  public:
    void setup  ();
    void update ();
    void draw   ();
    void exit   ();

    void keyPressed    (ofKeyEventArgs&);
    void keyReleased   (ofKeyEventArgs&);
    void mouseMoved    (ofMouseEventArgs&);
    void mouseDragged  (ofMouseEventArgs&);
    void mousePressed  (ofMouseEventArgs&);
    void mouseReleased (ofMouseEventArgs&);
    void windowResized (ofResizeEventArgs&);
    void gotMessage    (ofMessage);
    void dragEvent     (ofDragInfo);

    Controler lbm_controler;

    ofxColorMap colormap;

    //GUI
    ofxPanel gui;
    ofParameter<float> viscosity;
    ofParameter<ofVec2f> acceleration;

    ofxToggle model;

    std::vector<Particle> particles;


    void viscosityChanged(float& v);
    void accelerationChanged(ofVec2f &a);
    void modelChanged(bool &b);

};


#endif // _OFAPP_HPP
