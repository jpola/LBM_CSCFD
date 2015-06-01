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
    void saveVTK ();

    Controler lbm_controler;
    int Lx, Ly;

    ofxColorMap colormap;

    //GUI
    ofxPanel gui;
    ofParameter<float> viscosity;
    ofParameter<ofVec2f> acceleration;
    //ofParameter<float> reynolds;
    ofxLabel reynolds;

    ofxToggle model;

    ofxButton saveBtn;

    std::vector<Particle> particles;


    void viscosityChanged(float& v);
    void accelerationChanged(ofVec2f &a);
    void modelChanged(bool &b);
    void reynoldsChanged();

};


#endif // _OFAPP_HPP
