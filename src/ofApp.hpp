#pragma once
#ifndef _OFAPP_HPP_
#define _OFAPP_HPP_


#include <ofMain.h>
#include <ofxGui.h>

#include "lbmController.hpp"
#include "colormap/ofxColorMap.h"

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

    //input parameters
    float viscosity;

    ofVec2f acceleration;

    //GUI
    ofxPanel gui;
    ofParameter<float> visc;

    //void viscosityChanged(float visc);
};


#endif // _OFAPP_HPP
